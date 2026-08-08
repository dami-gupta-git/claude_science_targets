---
name: depmap-local
description: Query a locally downloaded DepMap release - CRISPR gene effect/dependency (Chronos), cell-line metadata, hotspot mutations, and PRISM Repurposing compound sensitivity. Use when asked whether cancer cells depend on a gene, which lineages or genotypes are selectively dependent, whether a gene is common-essential or selective, for synthetic-lethal / biomarker contrasts (e.g. WRN in MSI, KRAS in KRAS-mutant lines), or to find PRISM compounds annotated against a target. Handles the wide 18k-column matrices without loading them whole.
---

# DepMap (local release)

Answers gene-level dependency questions from DepMap files on disk. No network:
depmap.org sits behind a bot-verification wall, so files are downloaded manually.

## Data location

No built-in default: set `$DEPMAP_ROOT` to the directory holding the release,
or pass `root=` to any helper - `depmap_root()` raises, naming the variable,
when neither is set. Call `depmap_inventory()` first - it lists present files
and whether the cache is built.

| File | Shape | Notes |
|---|---|---|
| `CRISPRGeneEffect.csv` | 1208 lines x 18531 genes | Chronos score. **0 = no effect, -1 ~ median common essential** |
| `CRISPRGeneDependency.csv` | same | probability of dependency (0-1) |
| `Model.csv` | 2154 lines x 49 | `OncotreeLineage`, `OncotreePrimaryDisease`; index `ModelID` (`ACH-######`) |
| `OmicsSomaticMutationsMatrixHotspot.csv` | 1968 x 553 | hotspot counts; **filter `IsDefaultEntryForModel=="Yes"`** |
| `OmicsFusionFiltered.csv` | 79030 calls / 1719 models | filter `IsDefaultEntryForModel=="Yes"` |
| `CRISPRInferredCommonEssentials.csv` | 1827 genes | DepMap's own call - see caveat below |
| `Repurposing_Public_<REL>_*` | 24Q2: 6790 x 919; 23Q2: 6658 x 919 | PRISM primary screen, **single-dose** LFC at 2.5 uM |
| `prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv` | 1489 compounds x ~500 lines | **8-point dose-response**: fitted `auc`/`ic50`/`ec50`. Higher AUC = MORE resistant |

Only 1208 of 2154 models have CRISPR data - always report the joined n.

## Build the cache first (one-off, ~15 s)

The wide CSVs cost ~1 s per gene per query. `build_depmap_cache()` writes
column-addressable parquet to `_cache/`, taking a single-gene read to ~0.2 s
(5x) and 8 genes to ~0.17 s. Requires `pyarrow`. Values are float32 - agreement
with the CSV is to ~1e-7, which is far below any biologically meaningful
difference, but do not use the cache for exact bit-reproduction.

## Wide-matrix reads go through one contract

`depmap_read_matrix(kind, columns)` is the single entry point for every wide
matrix (`'effect'`, `'dependency'`, `'hotspot'`, `'damaging'`). It picks the
parquet cache when present and the raw CSV otherwise, and both backends obey the
same contract: columns returned in request order with the Entrez suffix stripped,
a `missing` list for absent symbols, `KeyError` when none exist, and
`FileNotFoundError` when no backend is available. `depmap_matrix_backend(kind)`
reports which one served the read.

This consolidation is load-bearing: three separate cache-or-CSV implementations
previously drifted apart, and every divergence was a bug (silent empty frame on
one path vs `KeyError` on the other; a leaked pyarrow error instead of a domain
error). Add new matrices to `MATRIX_SPECS` rather than writing a fourth reader.
`depmap_read_cached()` remains as a deprecated shim.

## Workflow

1. `depmap_inventory()` - confirm files and cache.
2. `depmap_selectivity(gene)` - is it common-essential, selective, or non-essential?
3. `depmap_lineage_enrichment(gene)` - which lineages are selectively dependent (BH-adjusted).
4. `depmap_mutation_contrast(dep_gene, marker_gene)` - genotype-driven dependency.
4b. `depmap_fusion_contrast(gene, canonical_only=True)` - dependency in
   fusion-positive lines; separates *untestable* (no such line) from a real negative.
5. `depmap_prism_compounds(target)` - compounds annotated to the target.
   Defaults to the **newest** release on disk; `depmap_prism_releases()` lists
   them, `release='23Q2'` pins an older one.

## Interpreting gene effect

- `< -1.0` strong dependency; `< -0.5` the conventional dependent cut (`DEP_THRESHOLD`).
- `frac_dependent >= 0.90` -> common essential: usually a poor therapeutic window.
- **Positive effect means knockout HELPS growth** - typical of tumour suppressors
  (TP53 mean +0.42). Never read a positive score as a weak dependency.
- A flat, low-variance profile near zero *may* mean the panel lacks the relevant
  context - but check with `depmap_fusion_contrast()` before assuming so. For ALK
  the panel does contain 19 canonical-fusion lines with CRISPR data and they are
  **still not dependent** (mean -0.16, p = 0.69); ROS1 likewise (n=6, p = 0.97).
  RET is genuinely untestable (n=1).

## Validated behaviour

Reproduced on this release, and worth re-running after any change:

- Pan-essentials RAN/RPL23/PSMA1 -> 100% of lines dependent (RAN mean -4.055).
- Olfactory controls OR5A1/CSN1S1 -> ~1% dependent.
- **WRN in MSI lines** (Chan 2019): 52% vs 4% dependent, p = 2e-12.
- **KRAS by lineage** (effect size): pancreas -2.14, bowel -1.43, biliary -1.00,
  lung -0.88. Ranked by BH-adjusted p the order differs - lung precedes biliary
  (126 lines vs 34), so report effect size and significance separately.
- **KRAS-hotspot-mutant vs WT KRAS dependency**: -1.92 vs -0.52, d = -3.3.
- Negative control TP53-hotspot vs KRAS dependency: d = -0.14, p = 0.45.

## Pitfalls

- Matrix columns are `SYMBOL (ENTREZ)`; helpers strip the Entrez suffix. The
  first header field is empty - select by position, not by name.
- Omics files key on `SequencingID`/`ModelConditionID`; collapse to `ModelID`
  with `IsDefaultEntryForModel=="Yes"` or rows multiply per model.
- `Model.csv` must be indexed by `ModelID` before joining. A bare `read_csv`
  leaves it as a column and every lineage join silently returns **zero rows** -
  assert the joined n before interpreting.
- **DepMap's `CRISPRInferredCommonEssentials.csv` (1827 genes) is not a frequency
  threshold.** It includes KRAS (52% of lines dependent) and genes at 0%. Our
  `frac_dependent >= 0.90` rule gives 974 genes, a strict SUBSET of theirs.
  Adopting DepMap's list wholesale would flag KRAS as having no therapeutic
  window; `depmap_common_essentials()` exposes it for cross-reference instead.
- Very small p-values underflow to 0.0; `floor_pvalue()` clamps them.
- **The secondary screen is a separate, older release.** `Repurposing_Public_23Q2`/`24Q2` are
  primary single-dose only, whatever the version number suggests; the dose-response data lives in
  the `prism-repurposing-20q2-secondary-*` files. Only primary-screen HITS were carried into the
  secondary screen (1489 of ~6800 compounds), so check membership before promising the analysis -
  read `usecols=["depmap_id","name","auc","ic50","r2"]`, lowercase `name`, and average replicate
  rows per (`depmap_id`, `name`). **AUC polarity is inverted relative to the primary screen's LFC**:
  higher AUC = more resistant. AUC is also less proliferation-confounded than single-dose viability
  (on one worked target, panel-wide correlation with expression r = -0.13 for AUC vs -0.41 for GDSC
  lnIC50), so prefer it whenever the compound is present.
- **PRISM and CRISPR are different releases.** 23Q2 and 24Q2 share the *same*
  919-line panel, of which **732 overlap the CRISPR panel (61%)**. Coverage is
  uneven by lineage (Lymphoid 44%, Bone 43%, Eye 20%), so a drug-side answer
  rests on roughly half the lines in those lineages. Always report the joint n.
- **One drug name can span several compound IDs** - 207 of 6575 names in 24Q2,
  up to 4 (DOXYCYCLINE). AZ-628 has two BRD IDs differing only by batch suffix,
  and they behave differently (median LFC -0.85 vs +0.15). Never aggregate PRISM
  by `Drug.Name`; the compound ID is the unit.
- A compound screened in two screens shares one matrix row (AZ-628 in
  REP.1M+REP.300); `depmap_prism_compounds()` collapses those into `screens`.
  In 24Q2 exactly one matrix row has no compound-list annotation and is
  unreachable by target.
