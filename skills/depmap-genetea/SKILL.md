---
name: depmap-genetea
description: Find a gene's codependencies in DepMap CRISPR data - other genes whose knockout profiles track it across cell lines, which proposes complex members and pathway partners - then name what a gene set shares using GeneTEA overrepresentation analysis over free-text gene descriptions, with each term traced back to the RefSeq/UniProt/Alliance sentence that produced it. Use when asked what a gene's functional partners are, what a hit list or cluster has in common, to annotate an unnamed gene set or CRISPR screen hits, for guilt-by-association or "what pathway is this", or when a gene needs a plain-language functional summary rather than a GO term.
---

# DepMap codependency x GeneTEA

Two questions this answers that the other DepMap skills do not: **which genes
behave like this one across the CRISPR panel**, and **what does a set of genes
have in common, in words**.

Codependency is guilt-by-association on 1,208 cell lines. If knocking out gene A
and gene B hurts the same lines and spares the same lines, they are likely in the
same complex or pathway. That proposes partners but does not name them. GeneTEA
(Boyle et al. 2025, Genome Biology 26:376, doi:10.1186/s13059-025-03844-8) supplies the name by treating gene
descriptions as a corpus and running overrepresentation analysis over the words
themselves, so a "term" is a phrase like `~ orotate` rather than a GO id, and
every hit traces back to the sentence that produced it.

**Scope.** Per-gene dependency, selectivity, lineage and genotype contrasts, and
PRISM compounds are `depmap-local`'s. Open Targets evidence joins are
`depmap-fusion`'s. This skill adds the gene-to-gene axis and the naming step.

## Setup

Requires `environment="genetea"` (Python 3.10, pandas 1.5.3, scikit-learn pinned
to 1.4.0 — the version the model was trained with; 1.4.2 loads but warns on six
estimators). Two inputs, each named in the error when absent:

| Input | Path | Note |
|---|---|---|
| Trained GeneTEA model | `$DEPMAP_ROOT/genetea/GeneTEA.pkl` | 1.11 GB, md5 `31fea7c4d9d56d59263bd982d426508f`. Figshare [10.6084/m9.figshare.28635317](https://doi.org/10.6084/m9.figshare.28635317), file 53289398 |
| DepMap parquet cache | `$DEPMAP_ROOT/_cache/` | Build with `depmap-local`'s `build_depmap_cache()` |

The model loads in 2.6 s and holds ~1.7 GB resident; `genetea_load()` memoises
per path so a session pays that once. `codep_prepare()` needs **all** 18,531
columns, so the column-addressable CSV fallback in `depmap-local` cannot serve
it — the parquet cache is required, not an optimisation.

Downloading the model needs both `ndownloader.figshare.com` and
`s3-eu-west-1.amazonaws.com` allowlisted: Figshare 302-redirects to S3 with a
signed URL valid for 10 seconds, so the redirect must be followed in one request.

`depmap.org/portal/*` is behind a bot-verification wall and returns a
verification page for every path, including the documented API routes; the
hosted `genetea-api` returns 502. Neither is a usable programmatic surface, which
is why the model runs locally here.

## Workflow

1. `codep_query_quality(gene)` — is this query interpretable at all? Refuse to
   read a partner list from a gene that fails.
2. `depmap_codependencies(gene, n=25)` — ranked partners with `r`, `z`, `n_lines`.
3. `codep_flag_proximity(table, gene)` — mark partners within 5 Mb of the query
   locus as possible copy-number artefacts.
4. `codep_enrich(gene)` — steps 2 and 3 plus GeneTEA, background set correctly.
   This is the one-call entry point.
5. `genetea_term_context(term, genes)` — the sentences behind a term.

For a gene set you already have (screen hits, a cluster), skip to
`genetea_enrich(genes, background=...)`. For a whole-genome ranked vector, use
`genetea_enrich_continuous(x)`.

## Read rank and z, never absolute r

**There is no value of r that separates real partners from background, and none
is shipped.** Over a 300-gene random null, the median top-1 r for an arbitrary
screened gene is **0.270**, and 88.7% reach 0.20, 37% reach 0.30. True curated
partners run down to r = 0.164 while non-partners inside the same top-15 slices
reach 0.550. The olfactory-receptor control OR5A1 returned OR5AN1 at r = 0.345 —
**higher than SMARCA4's best true partner** SMARCB1 at 0.303.

`CODEP_STRONG_R = 0.2` exists only as a convenience cutoff for the shared-partner
view and its docstring says so. Within one query, rank and `z` are interpretable;
across queries only `z` is, since per-query background r s.d. runs 0.043–0.091.

What does discriminate is the **query's own profile**: `profile_sd >= 0.20 and
min_effect <= -1.0` admits 9/9 complex-forming queries, 0/1 control, and 18.3% of
random screened genes. A gene that is never a dependency anywhere has a profile
made of screen noise plus residual copy number, and correlating that against
18,531 genes still returns large r. `codep_enrich()` enforces this gate by
default; override with `check_query=False` only deliberately.

## The background must be the screened genes

GeneTEA's corpus is 35,818 genes; this release screens 18,531. Enriching a
codependency list against the whole corpus credits terms for being *screenable* —
`~ DNA`, `~ transcription`, `~ kinase`, `~ cancer` — rather than for being shared
by the partners. The correct background is screened **intersected with** the
corpus (17,905 genes at `min_lines=200`); 28 screened symbols are absent from the
corpus and inflate the denominator if left in.

Across six real queries the whole-corpus background returns 427 terms against
276: **185 drop out**, 34 are gained, and FDR shifts by a median of +1.77 log10
with 236 of 242 surviving pairs becoming less significant. For SOX10 the top term
changes outright. `codep_enrich()` applies this by default.

The correction costs real power and is not purely a cleanup: VPS4A loses
`~ endosome` (corpus FDR 3.3e-3) and SOX10 loses `~ waardenburg`, both relevant.
Report the corrected run as primary; `corpus_background=True` gives the
uncorrected one for context.

## Validated behaviour

Ten queries x top 15 partners, expected sets written before inspecting output
(`codependency_validation.csv`). Correlations agree with `scipy.stats.pearsonr`
to max |Δr| = 5.3e-8.

| Query | Complex | Top partner | r | Curated recovered in top 15 |
|---|---|---|---|---|
| SMARCA4 | SWI/SNF (BAF) | SMARCB1 | 0.303 | 6/12 |
| ARID1B | SWI/SNF (BAF) | SMARCE1 | 0.270 | 6/11 |
| SMC2 | condensin | SMC4 | 0.259 | 5/7 |
| MRPL11 | mitoribosome | MRPL23 | 0.588 | 7 |
| BMS1 | 90S pre-ribosome | RPS2 | 0.277 | 6/31 |
| VPS4A | ESCRT-III | CHMP1A | 0.269 | 3/12 (paralogue VPS4B rank 2) |
| CTNNB1 | WNT | TCF7L2 | 0.698 | 2/18 |
| SOX10 | melanocyte lineage | BRAF | 0.572 | 6/14 |
| RAD21 | cohesin | CDC73 | 0.298 | **2/11 — see below** |
| OR5A1 | negative control | OR5AN1 | 0.345 | rejected by the gate, not by r |

`codep_enrich("SMARCA4")` returns SMARCB1, SMARCD1, SMARCC1, ARID1A, SMARCE1,
BRD9 as top partners and names them `~ SWI`, `SNF`, `~ BAF155`, `~ chromatin`.

**Pan-essential queries give poor lists.** RAD21 is a dependency in 1,202 of
1,208 lines, which flattens profile variance and buried most of cohesin (SMC1A
rank 1,105). It passes the gate but the list is weak — check `frac_dep_lines`
before trusting a pan-essential query. And **mutually buffered paralogues can be
anti-correlated**: STAG1 sits at r = −0.030, rank 11,603 for RAD21. Absence from
the positive tail is not absence of a relationship.

## Continuous mode

`genetea_enrich_continuous(x)` takes a whole-genome ranked vector — a
codependency vector is the natural input. Three things to know:

- The floor is **>10,000 values after ID alignment** to the corpus, not as
  supplied (18,530 → 18,507 for TP53), so a 10,050-gene input can fail a check
  the caller thought it passed. The wrapper's `ValueError` names both counts.
- **GeneTEA's own default `corr_thresh=0.2` returns zero terms on any real
  DepMap vector.** Term correlations are computed against a sparse TF-IDF column
  and are an order of magnitude below gene-gene correlations: the maximum
  attainable over all 24,831 terms was 0.104, 99.9th percentile 0.064.
  `DEFAULT_CORR_THRESH = 0.03` is shipped instead (281 terms → 74 rows).
  `get_enriched_continuous` **re-applies** `corr_thresh` to a precomputed
  `corr_terms` frame, so pass the same value to both calls or it silently
  returns empty.
- `method="permutation"` in upstream GeneTEA hangs — it pickles the 24 MB TF-IDF
  matrix into every `multiprocessing.Pool` task and never completes under the
  spawn start method. Reimplemented serially at 14 ms/permutation. The f-test
  path is 0.26 s and stays the default.

On the TP53 vector the top groups include `~ MDM2|~ MDM4++` and
`~ fanconi|~ repair|FA++`, all with negative correlation — the correct sign,
since MDM2/MDM4 dependency is anti-correlated with TP53's (−0.718, −0.522). Note
that **ribosome, nucleolus, mitoribosome and large no-expression gene families
(olfactory receptors, VCX/NBPF) dominate any Chronos correlation vector** and
appeared in 5 of 10 TP53 term groups; they are structure in the data, not a
query-specific finding.

## Terms carry a "~ " prefix

17,393 of 24,831 terms are synonym groups written `~ orotate`, and
`get_context`/`validate_terms` require that prefix verbatim — the bare form
raises `RuntimeError("no valid terms")`. `genetea_resolve_term()` accepts either
and resolves group members (`orotidine` → `~ orotate`); an unknown term raises
`KeyError` listing the closest vocabulary entries rather than substituting a
neighbour.

`genetea_term_context("orotate", ["CAD","DHODH","UMPS"])` returns six rows over
four sources. CAD is absent from the excerpts even though the term is enriched at
2 matching genes — its descriptions name the carbamoyl-phosphate steps upstream
of orotate. That asymmetry is the function working: the excerpt list shows which
genes actually carry the term, and a term supported only by `Names/Aliases` (a
symbol list, not curated prose) is visibly weaker evidence than the same term
from RefSeq or UniProt.

Excerpts are deduplicated on the way out. A synonym group matches once per
member, so `get_context` repeats the same sentence joined by `"..."` — UMPS's
UniProt paragraph came back four times, 1,357 characters for 337 characters of
distinct text. The repeat count carries no information the enriched table's
`n Matching Genes` columns do not already give.

## Other pitfalls

- **GeneTEA's `n` does not cap returned rows.** `n=25` returned 86 rows;
  `group_subterms=True` caps *term groups*, and every member row of an admitted
  group comes back. Slice top-N yourself or reported counts will be wrong.
- **Enriched frames carry trailing annotation columns** (`Stopword`, `IDF`,
  `Total Info`). Only the leading `DISCRETE_COLUMNS` / `CONTINUOUS_COLUMNS` are a
  stable contract — assert on the prefix, not exact column equality.
- `sort_by` changes order only in discrete mode, but **changes membership in
  continuous mode** (74 rows under `Effect Size` vs 67 under `Significance`),
  because the graph filter runs on the sort order.
- **Neighbouring-locus partners are copy number, not function.** OR5A1/OR5AN1 sit
  ~78 kb apart and share zero Brunello guides, so this is residual co-deletion
  Chronos does not fully remove, not guide crosstalk. `codep_flag_proximity()`
  flags rather than filters, since real tandem paralogues are also neighbours.
- Pairwise-complete NaN handling is required: 600 genes are screened in only
  ~70–80 lines by library version, so a complete-case intersection would collapse
  the usable line count and mean-imputation would shrink r toward zero for
  exactly those genes.
- `codep_prepare()` once per session — first query 0.51 s, warm queries 0.012 s.
  Do not hold the model (1.7 GB) alongside several matrix copies.

## Re-deriving the thresholds

`calibration.png` shows all three. Panel a is the r overlap, panel b the gate,
panel c the continuous threshold. The scripts that produced them ship in
`scripts/`; run them from the skill root with the root on `PYTHONPATH`, since
they import `kernel` and read the CSVs beside it:

```bash
cd <skill dir> && PYTHONPATH=. python scripts/gate_codependency.py
```

- `calibrate_codependency.py` — the 300-gene random null (writes
  `codependency_calibration_null.csv`). Slowest of the five.
- `validate_codependency.py` — the 10-query panel (writes
  `codependency_validation.csv`, `codependency_calibration_panel.csv`).
- `gate_codependency.py` — the gate grid, from the two files above.
- `background_comparison.py` — corpus vs screened background.
- `bench_continuous.py TP53` — continuous-mode timings and term correlations.

`gate_codependency.py` needs `codependency_calibration_panel.csv` and
`codependency_calibration_null.csv` present; both ship with the skill, so it
runs standalone, but the two scripts that generate them must run first if you
change the panel.

Tests: `python -m pytest -q` in the skill directory, environment `genetea`;
those needing the model or the release skip when absent.
