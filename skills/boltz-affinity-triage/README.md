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

- `SKILL.md` — the workflow: target resolution, control retrieval from each
  knowledge base, decoy matching, batched co-folding, gating, writing the run
  directory, and the report format.
- `kernel.py` — 21 helper functions, loaded into the Python kernel when the
  skill loads.
- `tests/test_kernel.py` — unit tests covering all 21 helpers. Run with
  `pytest tests/ -q` from the skill directory; no GPU, network, or Boltz install
  is required, and the suite passes with or without RDKit.

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

Retrieval and parsing:

- `kba_parse_nm(value)` — parses a BindingDB affinity string to
  `(nanomolar, qualifier)`, preserving `>` and `<` as censored measurements.
- `kba_pactivity(row, source='chembl')` — unified pActivity from a ChEMBL or
  BindingDB row, falling back to `standard_value` when `pchembl_value` is
  absent.
- `kba_leakage_flag(row, cutoff_year=2024)` — classifies a ChEMBL row as
  `likely_train`, `likely_novel` or `unknown` from its publication year.

Control construction:

- `kba_heavy_atoms(smiles)` — heavy-atom count, via RDKit when importable.
- `kba_props(smiles)` — heavy atoms, molecular weight, cLogP and ring count.
- `kba_match_decoys(actives, pool, per_active=1, tol=None)` — greedy
  nearest-neighbour selection of decoys matched on heavy atoms and molecular
  weight, sampling without replacement.
- `kba_null_auc(actives, decoys)` — AUC of a size-only scorer on the same
  controls.

Co-folding input and output:

- `kba_write_yaml(path, target_sequence, smiles, ligand_id='L', protein_id='A',
  msa=None, request_affinity=True)` — writes one Boltz-2 input YAML with the
  `properties: affinity:` block, returning `None` for ligands at or above the
  128-atom affinity cap.
- `kba_read_result(pred_dir, model_index=0, optimization_score_direction=None)`
  — reads the affinity JSON and the
  confidence JSON of one ranked diffusion sample, discovering score keys rather
  than assuming them, and normalises to `binder_score` where higher indicates a
  more likely binder. Confidence is per-sample and affinity is per-input, so
  the sample is selected by rank; `per_sample` returns every sample's
  confidence and `model_index` records which was used.

Scoring and gating:

- `kba_auc(pos, neg)` — Mann-Whitney AUC with tie correction.
- `kba_enrichment(scored, frac=0.1)` — enrichment factor over the top fraction
  of a ranked list.
- `kba_percentile(value, reference)` — position of a query score within the
  control distribution.
- `kba_gate(auc, n_pos, n_neg, null_auc=None, min_auc=0.65, min_per_class=5,
  min_margin=0.1)` — returns a verdict of `interpretable` or
  `not interpretable` with the reasons, failing on too few controls, an AUC at
  chance, or an AUC that does not clear the size-only null by the margin.
- `kba_tanimoto(smiles_a, smiles_b, radius=2, n_bits=2048)` — Morgan-fingerprint
  similarity, or `None` when RDKit is absent or either SMILES does not parse.
- `kba_nearest_active(smiles, actives)` — the closest known active to a query
  compound as `(compound_id, tanimoto)`.
- `kba_pose_flags(row, iptm_floor=None)` — the caveat strings for one scored
  compound: pose confidence below the floor (0.5 when unset), confidence
  absent, or no `binder_score` with the reason it could not be assigned.

Run outputs:

- `results_root(root=None)` — resolves this repo's `results/` directory from
  `$SCIENCE_RESULTS_ROOT` or an explicit `root=`; raises naming the variable
  when neither is set, rather than silently creating a `results/` folder
  wherever the kernel session's cwd happens to be.
- `kba_run_dir(target, root=None, topic=None, make=True)` — the run
  directory for one target, slugged to snake_case under the
  `boltz_affinity_triage` topic, with `scripts/` created beside it. `root`
  resolves via `results_root()`.
- `kba_write_table(path, rows, headers=None)` — CSV writer that keeps keys not
  named in `headers` rather than dropping them; returns `None` for empty input.
- `kba_markdown_table(rows, headers=None)` — the same rows rendered for the
  README.
- `kba_check_words(text, cap, label)` — raises when a README field exceeds its
  word cap, naming the field and the count.
- `kba_run_readme(target, gate, summary, controls, queries=None, files=(),
  data_sources=(), limits=(), enrichment=None, leakage_split=None, skipped=(),
  title=None)` — renders the run report, leading with the calibration verdict
  and replacing the ranked table with the diagnostic when the gate fails.
- `kba_write_run(out_dir, target, gate, controls, queries=None, summary=None,
  ...)` — writes the whole run directory and returns `{name: path}`, filling in
  query percentiles, nearest actives and pose flags from the controls.

## Calibration of the decoy-matching thresholds

The default matching tolerances and the `min_margin` gate were set against
EGFR (ChEMBL target CHEMBL203, UniProt P00533) using the retrieval recorded in
`decoy_confounder.png` and reproducible from `kernel.py`.

Known actives at pChEMBL >= 7 average 30.8 heavy atoms (n = 50). A default
weak-activity pull at >= 10 uM averages 16.0 heavy atoms, a difference of 14.8
atoms. Scoring these two sets with a function reading molecular size alone and
carrying no binding information gives AUC 0.968. Widening the weak pull to
1000 rows raises the pool mean to 26.1 heavy atoms and provides molecules large
enough to pair with the actives; matching against it brings the decoy mean to
30.9, and the same size-only scorer then gives AUC 0.494. `kba_null_auc`
reports this quantity so it can be compared against the real AUC, and
`kba_gate` fails a run whose margin over it is below 0.1.

Two retrieval properties constrain the workflow and are handled in the helpers.
Of 50 weak EGFR rows, 38 carry no `pchembl_value`, so a pipeline reading that
field alone retains 12 of its negatives; `kba_pactivity` derives the value from
`standard_value` and the unit instead. Separately, all 50 actives in this
retrieval trace to publications between 1997 and 2003, which precede the
affinity head's mid-2023 training cutoff, so `kba_leakage_flag` marks them
`likely_train` and the report format stratifies the AUC accordingly.

Ranking uses `affinity_probability_binary`, or `binding_confidence` on
Boltz-2.1. `optimization_score` is read and reported but not used for ranking
unless `optimization_score_direction` is supplied, because its sign convention is
undocumented and the published descriptions of the analogous Boltz-2 quantity
disagree: the model documentation defines `affinity_pred_value` as log10(IC50 in
µM), where lower values indicate tighter binding, while third-party summaries
describe the same field as a pIC50, where higher values do. An output carrying
only `optimization_score` returns `binder_score` as `None` together with an
`unscored_reason`, rather than an ordering that may be inverted.

Boltz writes one confidence file per diffusion sample, ranked so that
`model_0` is its best, while the affinity head runs once per input. Confidence
must therefore be attributed to a named sample: on a five-sample run whose
ranked pose has ipTM 0.81, the remaining samples scored 0.28 to 0.44, so
attributing an unranked sample's value to the prediction can move the reported
ipTM across the 0.5 threshold used to flag unreliable poses.

Without RDKit, `kba_heavy_atoms` falls back to a regular-expression estimate
with a mean absolute error of 0.36 atoms and a maximum of 2 over 100 ChEMBL
structures, and molecular-weight matching is skipped.

## Tests

`tests/test_kernel.py` synthesises co-folding output from the documented
`predictions/` layout and uses literal dicts shaped like ChEMBL and BindingDB
rows, so the suite runs without a GPU, a network connection, or a Boltz install.
RDKit-dependent assertions skip when it is absent, and the regular-expression
fallback for heavy-atom counting is exercised in both configurations.

Coverage concentrates on the properties the workflow's conclusions depend on:
selection of the ranked diffusion sample, the `standard_value` fallback in
pActivity, qualifier preservation on censored BindingDB measurements, score-key
discovery across Boltz versions, the direction of `binder_score`, the size-only
null, and each of the gate's three failure modes. Threshold comparisons are
tested at exactly-satisfied margins, where binary floating point places a
difference such as `0.60 - 0.50` marginally below its nominal value.

The output layer is tested against a temporary directory rather than a real
run: the gate dicts and score rows are literals, so the suite covers what is
written and what is withheld without a co-folding step. Inverting the gate check
so that a ranking is always written, overwriting a caller-supplied percentile,
moving the ipTM floor from exclusive to inclusive, dropping the control-role
validation, removing the summary word cap, and discarding CSV columns not named
in `headers` each fail exactly one test.

Removing the `sorted()` call from the output-directory traversal does not change
any result, because the sample is selected by rank rather than by traversal
position; the sort is retained as redundant ordering rather than as the
mechanism. Suppressing rank selection, the pActivity fallback, qualifier
parsing, the null-margin check, or the 128-atom cap each fails between 1 and 9
tests.

## Scope

This skill does not run Boltz-2 or Chai-1 itself, does not predict structures
without a co-folder, and does not convert `affinity_pred_value` into a
calibrated potency — the quantity is treated throughout as an ordering. Model
invocation and output-file layout are documented in `boltz` and `chai1`; pose
generation without an affinity head is covered by `diffdock`; receptor
preparation and PDB entry selection by `get-protein-structure`. Protein-protein
affinity is outside the affinity head's scope, which accepts a single
small-molecule `ligand` chain as the binder.

The co-folding step has not been executed against a GPU in the environment where
this skill was authored. The retrieval, matching, gating and output-parsing
functions are validated against live ChEMBL and BindingDB responses; the
`boltz predict` invocation between them is specified from the `boltz` skill and
the upstream documentation.
