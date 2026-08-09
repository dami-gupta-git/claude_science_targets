# depmap-genetea

Answers two questions about groups of genes rather than single genes. Given one
gene, it ranks every other screened gene by how closely its CRISPR knockout
profile tracks the query across the DepMap cell-line panel, which proposes
protein-complex members and pathway partners without knowing anything about them
in advance. Given a set of genes, it names what they have in common by running
overrepresentation analysis over free-text gene descriptions with GeneTEA (Boyle
et al. 2025, Genome Biology 26:376, doi:10.1186/s13059-025-03844-8), so the
answer is a phrase such as `~ orotate` backed by the RefSeq or UniProt sentence
it came from, rather than an ontology identifier. The two steps compose:
partners in, named biology out.

Three DepMap skills already exist in this project and this one does not overlap
them. `depmap-local` answers whether cells need a given gene; `depmap-fusion`
reconciles that with human genetic evidence; `target-triage-public-data` sets
out the triage workflow. All of those work one gene at a time against a fixed
question. This skill adds the gene-to-gene axis and the naming step.

## Functions

**Codependency**

| Function | Description |
| --- | --- |
| `depmap_codependencies(gene, n=25, min_lines=200, dataset='effect')` | partners ranked by Pearson r of gene-effect profiles, with `z`, `rank`, `p`, `fdr` and the pairwise-complete `n_lines` behind each row. |
| `codep_query_quality(gene)` | whether the query's own profile can support a partner list at all: `profile_sd`, `min_effect`, `n_dep_lines`, `passes`, `reason`. |
| `codep_enrich(gene, n_partners=50)` | both halves in one call, with the background set correctly. Returns `(partners, terms)`. |
| `codep_enrich_background(min_lines=200)` | screened genes intersected with the GeneTEA corpus, the list `codep_enrich` passes as background. |
| `depmap_screened_genes(min_lines=None)` | symbols screened in this release. |
| `depmap_shared_codependencies(genes, min_queries=2)` | partners recurring across several queries, one r column per query. |
| `codep_flag_proximity(table, query, window_bp=5_000_000)` | adds `near_query_locus`, flagging possible copy-number artefacts. |
| `codep_prepare(dataset='effect')` | centred matrix and masks, cached per session; call once, then queries are fast. |
| `codep_gene_loci()` | symbol to chromosome and start, for the proximity check. |

**GeneTEA**

| Function | Description |
| --- | --- |
| `genetea_load(path=None)` | resolve and load the trained model, memoised per path so a session pays the load cost once. |
| `genetea_enrich(genes, n=10, max_fdr=0.05, background=None)` | enriched terms for a discrete gene list as a tidy frame; empty rather than `None` when nothing survives. |
| `genetea_correlated_terms(x, method='f-test')` | per-term correlation against a whole-genome vector. |
| `genetea_enrich_continuous(x, corr_thresh=0.03)` | the filtered, grouped view of the above. |
| `genetea_resolve_term(term)` | bare term or synonym-group member to the exact vocabulary string. |
| `genetea_term_context(term, genes, source=None)` | per-gene source excerpts behind a term, with the source named. |

**Run output**, writing to `results/depmap_genetea/<slug>/`

| Function | Description |
| --- | --- |
| `genetea_run_dir(name, results_root=None, topic=None, make=True)` | resolves the canonical output directory. |
| `genetea_write_run(out_dir, name, summary, gene=None, codependencies=None, ...)` | writes the tables and run README together. |
| `genetea_write_table(path, table, headers=None)` | a single table. |
| `genetea_markdown_table(rows, cols, top_n=10, fmt=None)` | rows as markdown. |
| `genetea_run_readme(name, summary, gene=None, codependencies=None, ...)` | the run README text alone. |
| `genetea_results_root(results_root=None)` | resolves the results root. |

Support: `depmap_root()`, `strip_entrez()`, `bh_adjust()` — duplicated from
`depmap-local` so this module imports standalone.

The module also defines helpers that carry no leading underscore only because
the skill sidecar loader reserves that prefix: `codep_matrix_path`,
`codep_pairwise_r`, `genetea_tidy`, `genetea_empty_frame`,
`genetea_split_excerpt`, `genetea_dedupe_sentences`, `genetea_check_words`,
`genetea_permutation_corr_serial` and `genetea_correlated_terms_permutation`.
They are implementation detail; call the functions above instead.

## Files

`SKILL.md` carries the guidance an agent loads mid-task. `kernel.py` holds both
halves, because `codep_enrich()` spans them. `calibration.png` shows the panels
behind the thresholds, and `scripts/` holds the scripts that produced them,
kept runnable rather than described — run from the skill root with
`PYTHONPATH=.`. The `*.csv` files and `continuous_bench_TP53.json` are what
those scripts read and write.

`pytest.ini` pins `rootdir`; without it pytest's upward discovery walks out of
the skill directory and aborts on a parent it cannot stat. Run the tests with
`python -m pytest -q` from the skill directory in the `genetea` environment.
Tests that need the trained model or the DepMap release skip when either is
absent, so the suite stays runnable on a bare checkout.

## What calibrated the thresholds

The numbers behind each finding below live in the `*.csv` files in this
directory and are plotted in `calibration.png`. Re-derive them by running
`scripts/` in order: `calibrate_codependency.py`, `validate_codependency.py`,
`gate_codependency.py`, `background_comparison.py`. The continuous-mode figures
in `continuous_bench_TP53.json` and
`continuous_term_correlations_TP53_top300.csv` have no script in `scripts/`.

**No absolute correlation threshold is shipped, because none works.** Against a
random-gene null the correlation ranges overlap in both directions: true curated
partners fall below what an olfactory-receptor control returns for its own best
hit. Rank and `z` are the reportable statistics, and `CODEP_STRONG_R` is
retained only for the shared-partner view.

**The query-side gate is what discriminates instead.** `profile_sd >= 0.20` with
`min_effect <= -1.0` admits every complex-forming query tested and rejects the
control; loosening or tightening either bound drops a real query.

**The enrichment background is the screened genes, not the corpus.** GeneTEA's
corpus is roughly twice the screened set, so using it inflates the term count
with generic screenability vocabulary and shifts FDR throughout — VPS4A loses
`~ endosome`.

**Continuous mode needs its own threshold.** GeneTEA's default `corr_thresh`
returns zero terms on any real DepMap vector, because the maximum attainable
correlation over the vocabulary sits well below it; `DEFAULT_CORR_THRESH = 0.03`
yields a usable set. `CODEP_MIN_LINES = 200` is not knife-edge — the per-gene
valid-line count is bimodal with nothing in the middle, so any floor across a
wide range selects the same set.

## Scope

This skill does not score a single gene's dependency, classify it as
common-essential or selective, or contrast dependency by lineage, genotype or
fusion status — that is `depmap-local`. It does not join to human genetic
evidence or return a target verdict, which is `depmap-fusion`, and it does not
cover compound sensitivity or the wider triage sequence, which are
`depmap-local`'s PRISM helpers and `target-triage-public-data`. Codependency is
correlation across a panel of proliferating 2D cultures: it proposes partners
and does not establish physical interaction, and relationships absent from that
panel are invisible here. GeneTEA describes what has been written about genes,
so it inherits the literature's biases and says nothing about genes nobody has
characterised.
