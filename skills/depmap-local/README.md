# depmap-local

[DepMap](https://depmap.org) used CRISPR to switch off each of ~18,000 genes,
one at a time, across ~1,200 cancer cell lines, then measured whether the cells
still grew. This skill reads a locally downloaded DepMap release from disk and
turns those wide matrices into target-level answers: whether cancer cells depend
on a gene, how broadly, which lineages or genotypes carry the dependency, and
which PRISM compounds hit the target. The released files are CSV matrices of
thousands of columns keyed on identifiers that do not join directly to the
metadata table, so each question otherwise costs a parse, a rename, and a join
before any statistics. All reads are local; no network service is queried.

## Functions

Inventory and cache:

- `depmap_inventory()` — which release files are present, and whether the
  parquet cache is built.
- `depmap_root()` — resolved data root.
- `build_depmap_cache()` — one-off parquet build, ~15 s.

Matrix and metadata reads:

- `depmap_models(columns=None)` — `Model.csv` indexed by `ModelID`.
- `depmap_gene_effect(genes, dataset='effect')` — cell lines x genes matrix;
  `dataset` selects `effect` (Chronos) or `dependency` (probability) only.
  `hotspot`/`damaging` are mutation-count matrices, not gene-effect data —
  read them via `depmap_read_matrix` directly, as `depmap_mutation_contrast`
  does.
- `depmap_read_matrix(kind, columns)` — column subset of any named matrix,
  `kind` one of `effect`, `dependency`, `hotspot`, `damaging`.
- `depmap_read_cached(stem, genes)` — read the parquet cache directly.
- `depmap_matrix_backend(kind)` — whether a matrix resolves to cache or CSV.
- `strip_entrez(name)` — drops the trailing Entrez id from a column name.

Gene effect and selectivity:

- `depmap_selectivity(gene)` — `classification`, `mean_effect`,
  `frac_dependent`, `sd_effect`, `frac_positive`, effect percentiles.
- `depmap_common_essentials()` — DepMap's inferred essential set (1827 genes).

Contrasts:

- `depmap_lineage_enrichment(gene, min_lines=15)` — per-lineage mean effect,
  Mann-Whitney vs all other lines, `p_bh`.
- `depmap_mutation_contrast(dep_gene, marker_gene)` — mutant vs wild-type,
  `cohens_d`, p.
- `depmap_msi_contrast(gene)` — MSI-high vs MSS (MSI is a score, not a gene).
- `depmap_fusion_contrast(gene, canonical_only=False)` — `status` plus
  fusion-positive and fusion-negative means, or `'untestable'`.
- `depmap_fusion_models(gene)` — fusion-positive `ModelID`s.
- `depmap_fusion_profiled_models()` — models with fusion calls, the denominator
  for a fusion contrast.
- `depmap_stratified_lineage_contrast(gene, stratum)` — dependency within a
  lineage, split by a marker.

PRISM:

- `depmap_prism_releases()` — PRISM releases on disk (`['23Q2','24Q2']`).
- `depmap_prism_compounds(target, release=None)` — compounds annotated against
  a target, with a potency summary.

Statistics helpers:

- `bh_adjust(pvals)` — Benjamini-Hochberg adjustment.
- `floor_pvalue(p)` — floors p at the representable minimum.

## Gene effect

Every other result is derived from gene effect.

| Score | Meaning |
|---|---|
| `0` | Switching the gene off changed nothing |
| `-1` | Cells grow badly without it |
| `-2` or lower | Cells more or less die without it |
| **positive** | Cells grow **better** without it (a brake was removed) |

`-0.5` is the conventional "this line is dependent" cut-off (`DEP_THRESHOLD`).

A gene whose loss kills cancer cells is a plausible drug target — a drug
blocking that protein should do something similar. Two qualifications: a gene
that kills **every** cell (ribosomes, polymerases) has no therapeutic window,
and a **positive** score means inhibiting it would help the tumour.

## Reference results

These calls reproduce published findings and serve as regression checks.

```python
depmap_selectivity("WRN")
# classification='selective', mean_effect=-0.18, frac_dependent=0.068

depmap_msi_contrast("WRN")            # synthetic lethality, Chan et al. 2019
# MSI-high: 52% of lines dependent | MSS: 3.8% | p = 2e-12

depmap_lineage_enrichment("KRAS")
# Pancreas -2.14 (p_bh 1e-20) | Bowel -1.43 (5e-13) | Lung -0.88 (1e-02)

depmap_mutation_contrast("KRAS", "KRAS")
# mutant -1.92 vs WT -0.52, cohens_d = -3.3

depmap_fusion_contrast("ALK", canonical_only=True)
# status='not-dependent-despite-fusion', n=19, p=0.69   <- a real negative
depmap_fusion_contrast("RET", canonical_only=True)
# status='untestable', n=1                              <- absence of evidence

depmap_prism_compounds("BRAF")        # 21 compounds, newest release on disk
# top: KG-5, median_lfc -1.36
```

Only 6.8% of lines depend on WRN, and that 6.8% is the MSI-high subset — a
selective dependency is read together with the stratum that explains it.
`not-dependent-despite-fusion` and `untestable` are separate statuses because a
null result and an underpowered one differ; check `n` against
`depmap_fusion_profiled_models()` before reading either.

## Setup

No built-in default root: set `$DEPMAP_ROOT` to the directory holding the
release, or pass `root=` to any helper. `depmap_root()` raises, naming the
variable, when neither is set. Run `build_depmap_cache()` once — it takes a
single-gene read from
~1 s to ~0.2 s and requires `pyarrow`. Only **1208 of 2154** models have CRISPR
data; always report the joined n.

## Data properties

- `Model.csv` must be indexed by `ModelID`. Left as a column, lineage joins
  return zero rows with no error.
- Matrix columns are `SYMBOL (ENTREZ)` and the first header field is empty.
- Omics files key on sequencing IDs, not model IDs; collapse them with
  `IsDefaultEntryForModel=="Yes"`.
- One drug name can span several PRISM compound IDs (207 of 6575 in 24Q2), so
  aggregate by compound ID rather than by name.
- `CRISPRInferredCommonEssentials.csv` (1827 genes) is not a frequency
  threshold — it includes KRAS at 52%. The stricter rule applied by
  `depmap_selectivity()` gives 974 genes, a strict subset of it.

## Scope

This skill reads one local DepMap release: no downloading or release
management, no expression, copy-number, or proteomics matrices, no tractability,
genetic association, or literature evidence, and no dose-response fitting beyond
the single-dose PRISM summary. Target-level evidence outside DepMap comes from
Open Targets; fusion-partner structure and breakpoint work belongs to
`depmap-fusion`, whose README carries the cross-skill context.
