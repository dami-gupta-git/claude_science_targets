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

`depmap_inventory()` reports the shape of every file present; read current
dimensions from it rather than from this table.

| File | Contents |
|---|---|
| `CRISPRGeneEffect.csv` | Chronos score per line x gene. **0 = no effect, -1 ~ median common essential** |
| `CRISPRGeneDependency.csv` | the same matrix as probability of dependency (0-1) |
| `Model.csv` | cell-line metadata: `OncotreeLineage`, `OncotreePrimaryDisease`; index `ModelID` (`ACH-######`) |
| `OmicsSomaticMutationsMatrixHotspot.csv` | hotspot counts; **filter `IsDefaultEntryForModel=="Yes"`** |
| `OmicsFusionFiltered.csv` | fusion calls; filter `IsDefaultEntryForModel=="Yes"` |
| `CRISPRInferredCommonEssentials.csv` | DepMap's own essential call - see caveat below |
| `Repurposing_Public_<REL>_*` | PRISM primary screen, **single-dose** LFC at 2.5 uM |
| `prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv` | **8-point dose-response**: fitted `auc`/`ic50`/`ec50`. Higher AUC = MORE resistant |

Only about half the models in `Model.csv` have CRISPR data - always report the
joined n.

## Build the cache first

The wide CSVs cost roughly a second per gene per query. `build_depmap_cache()`
writes column-addressable parquet to `_cache/`, which cuts a single-gene read by
several fold and amortises further across multiple genes. Requires `pyarrow`.
Values are float32, so cache and CSV agree far below any biologically meaningful
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
- **Positive effect means knockout HELPS growth** - typical of tumour
  suppressors such as TP53. Never read a positive score as a weak dependency.
- A flat, low-variance profile near zero *may* mean the panel lacks the relevant
  context - but check with `depmap_fusion_contrast()` before assuming so. The
  panel holds enough canonical-fusion lines to test ALK and ROS1, and both are
  **still not dependent**; RET has too few to test at all, which the contrast
  reports as `untestable` rather than as a negative.

## Validated behaviour

These reproduce on this release and are worth re-running after any change. Run
each call for its current values.

- Pan-essentials RAN/RPL23/PSMA1 are dependent in every line; olfactory controls
  OR5A1/CSN1S1 in almost none. Together they bracket the usable range of
  `depmap_selectivity()`.
- **WRN in MSI lines** (Chan 2019): MSI-high lines are dependent, MSS lines are
  not - `depmap_msi_contrast("WRN")`.
- **KRAS by lineage**: pancreas carries the largest effect, then bowel, biliary
  and lung. Ranked by BH-adjusted p the order differs, because lung contributes
  several times as many lines as biliary and a smaller effect clears
  significance more easily - report effect size and significance separately.
- **KRAS-hotspot-mutant vs WT** shows a large dependency difference, while the
  TP53-hotspot negative control against the same gene shows none.

## Pitfalls

- Matrix columns are `SYMBOL (ENTREZ)`; helpers strip the Entrez suffix. The
  first header field is empty - select by position, not by name.
- Omics files key on `SequencingID`/`ModelConditionID`; collapse to `ModelID`
  with `IsDefaultEntryForModel=="Yes"` or rows multiply per model.
- `Model.csv` must be indexed by `ModelID` before joining. A bare `read_csv`
  leaves it as a column and every lineage join silently returns **zero rows** -
  assert the joined n before interpreting.
- **DepMap's `CRISPRInferredCommonEssentials.csv` is not a frequency
  threshold.** It includes KRAS, dependent in only about half of lines, and
  genes dependent in none. Our `frac_dependent >= 0.90` rule yields a strict
  SUBSET of theirs. Adopting DepMap's list wholesale would flag KRAS as having
  no therapeutic window; `depmap_common_essentials()` exposes it for
  cross-reference instead.
- Very small p-values underflow to 0.0; `floor_pvalue()` clamps them.
- **The secondary screen is a separate, older release.** `Repurposing_Public_23Q2`/`24Q2` are
  primary single-dose only, whatever the version number suggests; the dose-response data lives in
  the `prism-repurposing-20q2-secondary-*` files. Only primary-screen HITS were carried into the
  secondary screen - a minority of the primary set - so check membership before promising the
  analysis; read `usecols=["depmap_id","name","auc","ic50","r2"]`, lowercase `name`, and average
  replicate rows per (`depmap_id`, `name`). **AUC polarity is inverted relative to the primary
  screen's LFC**: higher AUC = more resistant. AUC is also less proliferation-confounded than
  single-dose viability, so prefer it whenever the compound is present.
- **PRISM and CRISPR are different releases.** 23Q2 and 24Q2 share the *same*
  cell-line panel, and only about three-fifths of it overlaps the CRISPR panel.
  Coverage is uneven by lineage - Lymphoid, Bone and Eye are among the
  thinnest - so a drug-side answer in those lineages rests on a fraction of the
  lines. Always report the joint n.
- **One drug name can span several compound IDs**, and they can behave
  differently: AZ-628 has two BRD IDs differing only by batch suffix whose
  median LFCs fall on opposite sides of zero. Never aggregate PRISM by
  `Drug.Name`; the compound ID is the unit.
- A compound screened in two screens shares one matrix row (AZ-628 in
  REP.1M+REP.300); `depmap_prism_compounds()` collapses those into `screens`.
  In 24Q2 exactly one matrix row has no compound-list annotation and is
  unreachable by target.
