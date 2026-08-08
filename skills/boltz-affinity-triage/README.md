# boltz-affinity-triage

Ranks candidate small molecules against a protein target by co-folding them with
Boltz-2 and calibrating the output against internal controls. Boltz-2 returns an
affinity estimate for any protein-ligand pair it is given, including pairs that
do not bind, so a score alone does not indicate whether a run has discriminative
power on the target in question. This skill retrieves known actives and
property-matched decoys from ChEMBL, BindingDB and PubChem BioAssay, folds them
alongside the query compounds under identical settings, and reports the ranking
only when the controls separate. The mechanics of running the model are covered
by the `boltz` and `chai1` skills; this skill supplies the experimental design
around them.

## Contents

- `SKILL.md` — the workflow: target resolution, control retrieval from each
  knowledge base, decoy matching, batched co-folding, gating, and the report
  format.
- `kernel.py` — 13 helper functions, loaded into the Python kernel when the
  skill loads.
- `tests/test_kernel.py` — unit tests covering all 13 helpers. Run with
  `pytest tests/ -q` from the skill directory; no GPU, network, or Boltz install
  is required, and the suite passes with or without RDKit.

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
