# depmap-local

[DepMap](https://depmap.org) CRISPR-knocked-out ~18,000 genes one at a time
across ~1,200 cell lines and measured growth. This skill reads a local
DepMap release and turns those wide, oddly-keyed CSV matrices into
target-level answers — dependency, breadth, lineage/genotype, PRISM hits —
without a parse/rename/join step each time. All reads are local.

## Functions

**Inventory and cache**

| Function | Description |
| --- | --- |
| `depmap_inventory()` | which release files are present, and whether the parquet cache is built. |
| `depmap_root()` | resolved data root. |
| `build_depmap_cache()` | one-off parquet build. |

**Matrix and metadata reads**

| Function | Description |
| --- | --- |
| `depmap_models(columns=None)` | `Model.csv` indexed by `ModelID`. |
| `depmap_gene_effect(genes, dataset='effect')` | cell lines x genes matrix; `dataset` selects `effect` (Chronos) or `dependency` (probability) only. `hotspot`/`damaging` are mutation-count matrices, not gene-effect data — read them via `depmap_read_matrix` directly, as `depmap_mutation_contrast` does. |
| `depmap_read_matrix(kind, columns)` | column subset of any named matrix, `kind` one of `effect`, `dependency`, `hotspot`, `damaging`. |
| `depmap_read_cached(stem, genes)` | read the parquet cache directly. |
| `depmap_matrix_backend(kind)` | whether a matrix resolves to cache or CSV. |
| `strip_entrez(name)` | drops the trailing Entrez id from a column name. |

**Gene effect and selectivity**

| Function | Description |
| --- | --- |
| `depmap_selectivity(gene)` | `classification`, `mean_effect`, `frac_dependent`, `sd_effect`, `frac_positive`, effect percentiles. |
| `depmap_common_essentials()` | DepMap's inferred essential set. |

**Contrasts**

| Function | Description |
| --- | --- |
| `depmap_lineage_enrichment(gene, min_lines=15)` | per-lineage mean effect, Mann-Whitney vs all other lines, `p_bh`. |
| `depmap_mutation_contrast(dep_gene, marker_gene)` | mutant vs wild-type, `cohens_d`, p. |
| `depmap_msi_contrast(gene)` | MSI-high vs MSS (MSI is a score, not a gene). |
| `depmap_fusion_contrast(gene, canonical_only=False)` | `status` plus fusion-positive and fusion-negative means, or `'untestable'`. |
| `depmap_fusion_models(gene)` | fusion-positive `ModelID`s. |
| `depmap_fusion_profiled_models()` | models with fusion calls, the denominator for a fusion contrast. |
| `depmap_stratified_lineage_contrast(gene, stratum)` | dependency within a lineage, split by a marker. |

**PRISM**

| Function | Description |
| --- | --- |
| `depmap_prism_releases()` | PRISM releases on disk (`['23Q2','24Q2']`). |
| `depmap_prism_compounds(target, release=None)` | compounds annotated against a target, with a potency summary. |

**Statistics helpers**

| Function | Description |
| --- | --- |
| `bh_adjust(pvals)` | Benjamini-Hochberg adjustment. |
| `floor_pvalue(p)` | floors p at the representable minimum. |

## Gene effect

Every other result is derived from gene effect.

| Score | Meaning |
|---|---|
| `0` | Switching the gene off changed nothing |
| `-1` | Cells grow badly without it |
| `-2` or lower | Cells more or less die without it |
| **positive** | Cells grow **better** without it (a brake was removed) |

`-0.5` is the conventional dependency cut-off (`DEP_THRESHOLD`). A gene whose
loss kills cancer cells is a plausible target, unless it kills **every** cell
(ribosomes, polymerases — no window) or scores **positive** (inhibiting it
would help the tumour).

## Reference results

These calls reproduce published findings and serve as regression checks. Run
them against the release on disk for the current values.

```python
depmap_selectivity("WRN")
# selective, not common-essential

depmap_msi_contrast("WRN")            # synthetic lethality, Chan et al. 2019
# MSI-high lines are dependent, MSS lines are not

depmap_lineage_enrichment("KRAS")
# pancreas strongest, then bowel, biliary and lung

depmap_mutation_contrast("KRAS", "KRAS")
# hotspot-mutant lines far more dependent than wild-type

depmap_fusion_contrast("ALK", canonical_only=True)
# status='not-dependent-despite-fusion'   <- a real negative
depmap_fusion_contrast("RET", canonical_only=True)
# status='untestable'                     <- absence of evidence

depmap_prism_compounds("BRAF")        # newest release on disk
```

WRN's small dependent fraction is the MSI-high subset — read a selective
dependency together with the stratum that explains it. `not-dependent-despite-fusion`
and `untestable` are separate statuses (a null differs from an underpowered
result); check `n` against `depmap_fusion_profiled_models()` before reading either.

## Setup

No built-in default root: set `$DEPMAP_ROOT` or pass `root=`; `depmap_root()`
raises, naming the variable, if neither is set. Run `build_depmap_cache()`
once to speed up per-gene reads (needs `pyarrow`). Only about half the models
in `Model.csv` have CRISPR data — always report the joined n.

## Data properties

- `Model.csv` must be indexed by `ModelID`. Left as a column, lineage joins
  return zero rows with no error.
- Matrix columns are `SYMBOL (ENTREZ)` and the first header field is empty.
- Omics files key on sequencing IDs, not model IDs; collapse them with
  `IsDefaultEntryForModel=="Yes"`.
- One drug name can span several PRISM compound IDs, and their potencies can
  disagree, so aggregate by compound ID rather than by name.
- `CRISPRInferredCommonEssentials.csv` is not a frequency threshold — it
  includes KRAS, which is dependent in only about half of lines. The stricter
  rule applied by `depmap_selectivity()` yields a strict subset of it.

## Scope

This skill reads one local DepMap release: no downloading or release
management, no expression, copy-number, or proteomics matrices, no tractability,
genetic association, or literature evidence, and no dose-response fitting beyond
the single-dose PRISM summary. Target-level evidence outside DepMap comes from
Open Targets; fusion-partner structure and breakpoint work belongs to
`depmap-fusion`, whose README carries the cross-skill context.
