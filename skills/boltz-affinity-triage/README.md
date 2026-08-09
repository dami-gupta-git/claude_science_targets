# boltz-affinity-triage

Ranks candidate small molecules against a protein target by co-folding them with
Boltz-2 and calibrating the output against internal controls. Boltz-2 returns an
affinity estimate for any protein-ligand pair it is given, including pairs that
do not bind, so a score alone does not indicate whether a run has discriminative
power on the target in question. This skill retrieves known actives and
property-matched decoys from ChEMBL, BindingDB and PubChem BioAssay, folds them
alongside the query compounds under identical settings, and reports the ranking
only when the controls separate. Each run writes a directory of result files
with a README of its own. The mechanics of running the model are covered by the
`boltz` and `chai1` skills; this skill supplies the experimental design around
them.

## Contents

| File | Description |
| --- | --- |
| `SKILL.md` | the workflow: target resolution, control retrieval from each knowledge base, decoy matching, batched co-folding, gating, writing the run directory, and the report format. |
| `kernel.py` | helper functions, loaded into the Python kernel when the skill loads. |
| `tests/test_kernel.py` | unit tests covering every helper. Run with `pytest tests/ -q` from the skill directory; no GPU, network, or Boltz install is required, and the suite passes with or without RDKit. |

## Run outputs

A run writes to `results/boltz_affinity_triage/<target>/`, one directory per
target under a shared topic so that runs stay comparable. `kba_write_run`
assembles it:

| file | contents |
| --- | --- |
| `README.md` | the run report: calibration verdict, ranked queries, files, data sources, limits |
| `calibration.json` | the `kba_gate` dict, AUC, size-only null AUC, EF@10%, leakage split, compounds skipped at the atom cap |
| `controls.csv` | one row per control compound: knowledge-base activity, leakage flag, co-folding scores |
| `ranked_queries.csv` | query compounds ranked by `binder_score`, written only when the gate passes |
| `scripts/` | the target-specific wiring that produced the run |

`ranked_queries.csv` is withheld on a failing gate rather than written with a
caveat, because a file on disk outlives the message that accompanied it and a
ranking from an uncalibrated run orders compounds by noise. `calibration.json`
records `ranking_written: false` in that case, and the README states why the
section is empty. Co-folding output — the YAMLs, the CIFs, the `out/` tree — is
not copied into the run directory; the three files above carry everything the
report rests on.

## Helpers

**Retrieval and parsing**

| Function | Description |
| --- | --- |
| `kba_parse_nm(value)` | parses a BindingDB affinity string to `(nanomolar, qualifier)`, preserving `>` and `<` as censored measurements. |
| `kba_pactivity(row, source='chembl')` | unified pActivity from a ChEMBL or BindingDB row, falling back to `standard_value` when `pchembl_value` is absent. |
| `kba_leakage_flag(row, cutoff_year=2024)` | classifies a ChEMBL row as `likely_train`, `likely_novel` or `unknown` from its publication year. |

**Control construction**

| Function | Description |
| --- | --- |
| `kba_heavy_atoms(smiles)` | heavy-atom count, via RDKit when importable. |
| `kba_props(smiles)` | heavy atoms, molecular weight, cLogP and ring count. |
| `kba_match_decoys(actives, pool, per_active=1, tol=None)` | greedy nearest-neighbour selection of decoys matched on heavy atoms and molecular weight, sampling without replacement. |
| `kba_null_auc(actives, decoys)` | AUC of a size-only scorer on the same controls. |

**Co-folding input and output**

| Function | Description |
| --- | --- |
| `kba_write_yaml(path, target_sequence, smiles, ...)` | one Boltz-2 input YAML with the `properties: affinity:` block; `None` for ligands at or above the 128-atom affinity cap. |
| `kba_read_result(pred_dir, model_index=0, ...)` | the affinity and confidence JSON of one ranked diffusion sample, normalised to `binder_score` where higher indicates a more likely binder. |

**Scoring and gating**

| Function | Description |
| --- | --- |
| `kba_auc(pos, neg)` | Mann-Whitney AUC with tie correction. |
| `kba_enrichment(scored, frac=0.1)` | enrichment factor over the top fraction of a ranked list. |
| `kba_percentile(value, reference)` | position of a query score within the control distribution. |
| `kba_gate(auc, n_pos, n_neg, null_auc=None, ...)` | a verdict of `interpretable` or `not interpretable` with the reasons, from control counts, AUC and the margin over the null. |
| `kba_tanimoto(smiles_a, smiles_b, ...)` | Morgan-fingerprint similarity, or `None` when RDKit is absent or a SMILES does not parse. |
| `kba_nearest_active(smiles, actives)` | the closest known active to a query compound as `(compound_id, tanimoto)`. |
| `kba_pose_flags(row, iptm_floor=None)` | the caveat strings for one scored compound: low or absent pose confidence, or a missing `binder_score`. |

**Run outputs**

| Function | Description |
| --- | --- |
| `results_root(root=None)` | this repo's `results/` directory, from `$SCIENCE_RESULTS_ROOT` or an explicit `root=`; raises when neither is set. |
| `kba_run_dir(target, root=None, ...)` | the run directory for one target, slugged to snake_case under the `boltz_affinity_triage` topic, with `scripts/` beside it. |
| `kba_write_table(path, rows, headers=None)` | CSV writer that keeps keys not named in `headers`; `None` for empty input. |
| `kba_markdown_table(rows, headers=None)` | the same rows rendered for the README. |
| `kba_check_words(text, cap, label)` | raises when a README field exceeds its word cap. |
| `kba_run_readme(target, gate, summary, controls, ...)` | the run report, leading with the calibration verdict and substituting the diagnostic for the ranked table when the gate fails. |
| `kba_write_run(out_dir, target, gate, controls, ...)` | writes the whole run directory and returns `{name: path}`. |

## Calibration of the decoy-matching thresholds

Set against one real retrieval, EGFR (CHEMBL203), via `decoy_confounder.png` /
`kernel.py`. Each item below is a way this workflow's numbers can look
meaningful while carrying no binding information.

**Size confound.** A naive weak-activity pull is smaller than actives (16.0 vs
30.8 heavy atoms), so a scorer reading size alone gets AUC 0.968. Size-matching
against a wider pool brings that null to 0.494. `kba_null_auc` reports it;
`kba_gate` fails a run whose margin over it is below 0.1.

**Retrieval quirks.** 38/50 weak EGFR rows lack `pchembl_value` —
`kba_pactivity` derives it from `standard_value` + unit. All 50 actives predate
the affinity head's training cutoff — `kba_leakage_flag` marks them
`likely_train` and reporting stratifies by it.

**Ranking column.** `optimization_score`'s sign convention is disputed, so
ranking uses `affinity_probability_binary`/`binding_confidence` instead; an
`optimization_score`-only output gets `binder_score=None` rather than a
possibly-inverted order.

**Per-sample confidence.** Boltz writes one confidence file per diffusion
sample, not per input — misattributing one can swing ipTM across the 0.5
threshold (0.81 vs 0.28-0.44 seen on the same run).

Without RDKit, `kba_heavy_atoms` falls back to a regex estimate (MAE 0.36
atoms) and molecular-weight matching is skipped.

## Tests

`tests/test_kernel.py` runs on literal dicts and a synthesised `predictions/`
layout — no GPU, network, or Boltz install needed. RDKit-only assertions skip
when it's absent.

Coverage centers on what conclusions depend on: ranked-sample selection, the
pActivity fallback, qualifier preservation, cross-version score-key discovery,
`binder_score` direction, the size-only null, the gate's three failure modes
(tested at exactly-satisfied floating-point margins), and the output layer's
writes/withholdings against a temp directory.

## Scope

This skill does not run Boltz-2 or Chai-1 itself, does not predict structures
without a co-folder, and does not convert `affinity_pred_value` into a
calibrated potency — the quantity is treated throughout as an ordering. Model
invocation and output-file layout are documented in `boltz` and `chai1`; pose
generation without an affinity head is covered by `diffdock`; receptor
preparation and PDB entry selection by `get-protein-structure`. Protein-protein
affinity is outside the affinity head's scope, which accepts a single
small-molecule `ligand` chain as the binder.

The co-folding step requires a GPU, which is not available for verifying this
skill directly. The retrieval, matching, gating and output-parsing functions
are validated against live ChEMBL and BindingDB responses; the `boltz predict`
invocation between them is specified from the `boltz` skill and the upstream
documentation.
