# depmap-fusion

Joins human genetic evidence from Open Targets with cell-line dependency from
DepMap and returns one verdict per candidate drug target. Open Targets scores a
gene's association with a disease from 0 to 1 but encodes no direction, so
tumour suppressors such as `TP53` rank near the top of most cancer lists and
inhibiting them would be harmful. DepMap measures whether ~1,200 cancer cell
lines still grow after CRISPR knockout of each of ~18,000 genes, but carries no
human disease evidence and cannot on its own separate a tumour-selective
dependency from a gene every cell needs. Fusing the two classifies each target
into one of ten verdicts and sets a `knockout_actionable` flag. In the two
cancers tested here, **about two-thirds** of the top-ranked Open Targets genes
turned out not to be viable knockout targets (36 of 55; lung 19/30, pancreatic
17/25).

## Functions

- `fuse_target_row(symbol, ot_score, effect_stats, lineage_row=None,
  tractable_modalities=None, dep_threshold=-0.5)` — classifies one target and
  returns the verdict plus the numbers behind it. `effect_stats` must be the
  **full** `depmap_selectivity()` dict: the verdict reads `sd_effect` and
  `frac_positive`, not just the mean. Pass `lineage_row` (from
  `depmap_lineage_enrichment()`) to enable the `context-restricted` verdict.
- `fusion_notes()` — the interpretation guide for the verdict vocabulary, as a
  dict keyed by verdict.

Module constants: `VERDICTS`, `EVIDENCE_FLOOR`, `SUPPRESSOR_MIN_FRAC_POS`,
`SUPPRESSOR_MIN_SD`, `INERT_MAX_SD`, `INERT_MAX_ABS_MEAN`.

## Companion skills

Evidence retrieval and dependency statistics live in two separate skills; this
one only adds the join and the verdict. Load all three for a verdict workflow.

| Skill | Answers | Source | Helpers |
|---|---|---|---|
| `depmap-local` | Do cancer cells need this gene? | Local DepMap files | 22 |
| `opentargets-evidence` | Is it linked to disease in humans, and is it druggable? | Open Targets GraphQL | docs only |
| `depmap-fusion` | Do those two agree? | both | 2 |

## Workflow

1. **Evidence** — in a `repl` cell, `open_targets_disease_targets(efo_id,
   size=N)` for a disease, or the dossier query for one target; collect
   `ot_score` and tractability modalities and write them to `handoff/*.json`.
2. **Dependency** — in a `python` cell, `depmap_selectivity(gene)` per target,
   plus the disease-relevant lineage slice.
3. **Fuse** — `fuse_target_row(...)` per gene; report the verdict mix and the
   joined n.

```python
st = depmap_selectivity("TP53")          # from depmap-local
fuse_target_row("TP53", ot_score=0.74, effect_stats=st)
# verdict='growth-suppressive-mismatch', knockout_actionable=False

st = depmap_selectivity("KRAS")
fuse_target_row("KRAS", ot_score=0.82, effect_stats=st)
# verdict='concordant-dependency', knockout_actionable=True
```

## The ten verdicts

| Verdict | Meaning | Actionable |
|---|---|---|
| `concordant-dependency` | Human evidence and cell dependency agree. Strongest KO/degrader case. | yes |
| `context-restricted` | Dependency confined to a lineage or genotype. Needs a biomarker. | yes |
| `dependency-with-thin-evidence` | `ot_score` in [0.10, 0.50) with `frac_dependent >= 0.10` — both halves real, evidence weaker than concordant. | yes |
| `common-essential-no-window` | Real dependency, but every cell has it — normal tissue included. | no |
| `growth-suppressive-mismatch` | Knockout **helps** growth (positive tail, high variance). A suppressor; pivot to synthetic lethality. | no |
| `inert-in-panel` | Flat near zero. Test with `depmap_fusion_contrast()` before blaming the panel. | no |
| `evidence-without-dependency` | Real human evidence (`ot_score >= EVIDENCE_FLOOR`) but no cell-autonomous dependency. | no |
| `dependency-without-evidence` | Cells need it, no human evidence at all (`ot_score` below `EVIDENCE_FLOOR` or absent). Possible novel target. | no |
| `no-evidence-no-dependency` | Neither side supports it — `ot_score` below `EVIDENCE_FLOOR` (0.10) or absent, AND no dependency. | no |
| `indeterminate` | Dispersion inputs missing, so the suppressive/inert tests could not run. | no |

Branch order is load-bearing: common-essential, then suppressive, then context,
then concordant, then inert; first match wins. Below the strong-evidence bar,
the thin-evidence tier is tested before either single-sided verdict, so a row
with both real evidence and a dependency signal is never labelled as lacking
one of them.

## Reading the numbers

- Gene effect `0` means knockout changed nothing, `-1` that cells grow badly
  without the gene (≈ median common essential), `-2` or lower that they largely
  die. A **positive** effect means knockout *helps* growth; never read it as a
  weak dependency. `TP53` sits at +0.42.
- `< -0.5` is the conventional dependent cut (`DEP_THRESHOLD`).
- `frac_dependent >= 0.90` marks a common essential — a poor therapeutic window.
- Separating a suppressor from an inert gene needs dispersion, not the mean:
  TP53 (sd 0.58, 42% of lines positive) is a suppressor; NRG1 (sd 0.08, 8%) is
  inert in this panel. Hence `SUPPRESSOR_MIN_FRAC_POS = 0.15` and
  `SUPPRESSOR_MIN_SD = 0.20`. These are conventions, not decision boundaries —
  report the underlying numbers alongside any verdict.

## Data layout

No built-in default root: `depmap-local` (the dependency-statistics half this
skill joins against) requires `$DEPMAP_ROOT` or an explicit `root=`.
`depmap_inventory()` reports what is present.

| File | Shape |
|---|---|
| `CRISPRGeneEffect.csv` | 1208 lines × 18531 genes (Chronos) |
| `CRISPRGeneDependency.csv` | same, as probabilities |
| `Model.csv` | 2154 lines × 49 metadata columns |
| `OmicsSomaticMutationsMatrixHotspot.csv` | 1968 × 553 |
| `OmicsFusionFiltered.csv` | 79030 fusion calls / 1719 models |
| `CRISPRInferredCommonEssentials.csv` | 1827 genes (DepMap's own call) |
| `Repurposing_Public_<REL>_*` | 24Q2: 6790 × 919 · 23Q2: 6658 × 919 |

Only **1208 of 2154** models have CRISPR data — always report the joined n.

`build_depmap_cache()` writes column-addressable parquet to `_cache/`
(872 MB CSV → 244 MB, already built). A single-gene read goes from ~1 s to
~0.2 s. Values are float32; agreement with the CSV is ~1e-7. Requires
`pyarrow`.

## Data constraints

- `Model.csv` must be indexed by `ModelID` before joining. A bare `read_csv`
  leaves it as a column and every lineage join returns **zero rows silently** —
  assert the joined n before interpreting anything.
- Omics matrices key on `SequencingID`/`ModelConditionID`; collapse to
  `ModelID` with `IsDefaultEntryForModel == "Yes"` or rows multiply per model.
- Open Targets `tractability` exposes `label`/`modality`/`value`; there is no
  `id` field, and `drugAndClinicalCandidates` takes no arguments.
- One PRISM drug name can span several compound IDs — 207 of 6575 in 24Q2, up
  to 4 for DOXYCYCLINE, and AZ-628's two batches give median LFC −0.85 and
  +0.15. Aggregate by compound ID, never by `Drug.Name`.
- Common essentials are excluded from `knockout_actionable`: POLD1, POLE, RRM1
  and PRIM2 are replication genes with no therapeutic window, so check
  `depmap_class` before trusting the flag.

## Validation

Checks that reproduce on this release; re-run after a DepMap version bump.

- Pan-essentials RAN/RPL23/PSMA1: 100% of lines dependent (RAN mean −4.055).
  Olfactory controls OR5A1/CSN1S1: ~1%.
- `WRN`: classification `selective`, mean_effect −0.18, frac_dependent 0.068.
  In MSI-high lines (Chan 2019) 52% are dependent against 3.8% of MSS lines,
  p = 2e-12 (`depmap_msi_contrast`) — the 7% that depend on it are the MSI
  lines. Negative control KRAS/MSI: p = 0.61.
- KRAS by lineage (by effect size): pancreas −2.14 (p_bh 1e-20), bowel −1.43,
  biliary −1.00, lung −0.88. Ranked by BH-adjusted p, lung precedes biliary —
  it has 126 lines against 34, so a smaller effect clears significance more
  easily.
- KRAS-hotspot-mutant vs WT KRAS dependency: −1.92 vs −0.52, Cohen's d = −3.3.
  Negative control (TP53-hotspot vs KRAS dependency): d = −0.14, p = 0.45.

## Coverage and limits

- **PRISM and CRISPR are different releases, and upgrading PRISM does not fix
  the overlap.** 23Q2 and 24Q2 use the *identical* 919-line panel; 24Q2 adds
  132 compounds (the REP.300 screen), not cell lines. Only **732 lines (61% of
  the CRISPR panel)** have both, and coverage is uneven by lineage — Lymphoid
  44%, Bone 43%, Eye 20%. Drug-side conclusions in those lineages rest on about
  half the panel. `depmap_prism_compounds()` defaults to the newest release on
  disk — `BRAF` returns 21 compounds with median_lfc, frac_killed_lfc_lt_1,
  screens and release; pass `release=` to pin one.
- DepMap measures proliferation in 2D culture. Immune, stromal, angiogenic and
  differentiation-dependent targets are invisible by construction.
- `inert-in-panel` is testable and often a real negative.
  `depmap_fusion_contrast()` separates "no fusion-positive line exists"
  (untestable — RET, n=1) from "fusion-positive lines exist and still show no
  dependency" (ALK n=19 p=0.69, ROS1 n=6 p=0.97). The ALK result holds even
  restricted to canonical EML4–ALK/NPM1–ALK drivers.
- Two essentiality definitions disagree. DepMap's inferred list (1827 genes)
  includes KRAS at 52% dependency; the `frac_dependent >= 0.90` rule gives 974,
  a strict subset. The stricter rule drives the no-window verdict — adopting
  DepMap's would discard KRAS — and DepMap's call is surfaced as the
  `depmap_inferred_essential` cross-reference field.
- `evidence-without-dependency` and `dependency-with-thin-evidence` both
  require `ot_score >= EVIDENCE_FLOOR` (0.10). Below the floor, or `None`, a
  row returns `no-evidence-no-dependency` or `dependency-without-evidence`
  according to its dependency signal.
- The thin-evidence tier is invisible in a top-25/30 sweep — every target in
  the two worked examples scores >= 0.509 — but dominates deeper slices: 269 of
  the top 300 lung targets score in [0.10, 0.50), and 57 of them carry
  `frac_dependent >= 0.10`, taking the actionable set from 12 to 69.
- No verdict is validation. `concordant-dependency` is a prioritisation signal.

## Scope

This skill does not retrieve evidence or compute dependency statistics. Open
Targets queries (`open_targets_disease_targets` and the target dossier query,
via the `clinical-genomics` connector in a `repl` cell) belong to
`opentargets-evidence`; gene-effect statistics, lineage and mutation contrasts,
fusion contrasts, PRISM compound summaries, the file inventory and the parquet
cache belong to `depmap-local`. It performs no chemistry, structure or
expression analysis.
