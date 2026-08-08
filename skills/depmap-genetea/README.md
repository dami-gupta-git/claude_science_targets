# depmap-genetea

This skill answers two questions about groups of genes rather than single genes.
Given one gene, it ranks every other screened gene by how closely its CRISPR
knockout profile tracks the query across 1,208 DepMap cell lines, which proposes
protein-complex members and pathway partners without knowing anything about them
in advance. Given a set of genes, it names what they have in common by running
overrepresentation analysis over free-text gene descriptions with GeneTEA (Boyle
et al. 2025, Genome Biology 26:376, doi:10.1186/s13059-025-03844-8), so the answer
is a phrase such as `~ orotate` backed by the RefSeq or UniProt sentence it came
from, rather than an ontology identifier. The two steps compose: partners in,
named biology out.

Three DepMap skills already exist in this project and this one does not overlap
them. `depmap-local` answers whether cells need a given gene; `depmap-fusion`
reconciles that with human genetic evidence; `target-triage-public-data` sets out
the triage workflow. All of those work one gene at a time against a fixed
question. This skill adds the gene-to-gene axis and the naming step.

## Files

- `SKILL.md` — the guidance an agent loads mid-task: setup, workflow, the two
  calibration findings, and the pitfalls.
- `kernel.py` — 26 module-level functions, auto-loaded into the kernel when the
  skill is loaded: the 18 documented below plus 8 internal helpers that carry no
  leading underscore only because the sidecar loader reserves that prefix. Both
  halves live here because `codep_enrich()` spans them.
- `tests/` — three files. `test_codependency.py` checks correlation
  correctness on a synthesised matrix with a known answer;
  `test_genetea_enrich.py` covers the wrappers; `test_bridge.py` covers the join.
- `pytest.ini` — pins `rootdir`. Without it pytest's upward discovery walks out
  of the skill directory and aborts on a parent it cannot stat.
- `calibration.png` — the three panels behind every threshold below.
- `scripts/` — the five scripts that produced every calibration number, kept
  runnable rather than described. Run from the skill root with `PYTHONPATH=.`.
- `*.csv`, `continuous_bench_TP53.json` — the calibration and validation data
  those scripts read and write.

Run the tests with `python -m pytest -q` from the skill directory in the
`genetea` environment. Tests that need the 1.1 GB model or the DepMap release
skip when either is absent, so the suite stays runnable on a bare checkout.

## Functions

Codependency:

- `depmap_codependencies(gene, n=25, min_lines=200, dataset='effect')` — partners
  ranked by Pearson r of gene-effect profiles, with `z`, `rank`, `p`, `fdr` and
  the pairwise-complete `n_lines` behind each row.
- `codep_query_quality(gene)` — whether the query's own profile can support a
  partner list at all: `profile_sd`, `min_effect`, `n_dep_lines`, `passes`,
  `reason`.
- `codep_enrich(gene, n_partners=50)` — the two halves in one call, with the
  background set correctly. Returns `(partners, terms)`.
- `codep_enrich_background(min_lines=200)` — screened genes intersected with the
  GeneTEA corpus, the list `codep_enrich` passes as background.
- `depmap_screened_genes(min_lines=None)` — symbols screened in this release.
- `depmap_shared_codependencies(genes, min_queries=2)` — partners recurring
  across several queries, one r column per query.
- `codep_flag_proximity(table, query, window_bp=5_000_000)` — adds
  `near_query_locus`, flagging possible copy-number artefacts.
- `codep_prepare(dataset='effect')` — centred matrix and masks, cached per
  session; call once, then queries cost ~12 ms.
- `codep_gene_loci()` — symbol to chromosome and start, for the proximity check.

GeneTEA:

- `genetea_load(path=None)` — resolve and load the trained model, memoised per
  path so a session pays the 2.6 s and 1.7 GB once.
- `genetea_enrich(genes, n=10, max_fdr=0.05, background=None)` — enriched terms
  for a discrete gene list as a tidy frame; empty rather than `None` when nothing
  survives.
- `genetea_correlated_terms(x, method='f-test')` — per-term correlation against a
  whole-genome vector.
- `genetea_enrich_continuous(x, corr_thresh=0.03)` — the filtered, grouped view
  of the above.
- `genetea_resolve_term(term)` — bare term or synonym-group member to the exact
  vocabulary string.
- `genetea_term_context(term, genes, source=None)` — per-gene source excerpts
  behind a term, with the source named.

Support: `depmap_root()`, `strip_entrez()`, `bh_adjust()` — duplicated from
`depmap-local` so this module imports standalone.

Internal helpers (`codep_matrix_path`, `codep_pairwise_r`, `genetea_empty_frame`,
`genetea_tidy`, `genetea_split_excerpt`, `genetea_dedupe_sentences`,
`genetea_permutation_corr_serial`, `genetea_correlated_terms_permutation`)
carry no leading underscore because the skill sidecar loader reserves
that prefix. They are implementation detail; call the functions above instead.

## What calibrated the thresholds

Every number here is read from a file in this directory and shown in
`calibration.png`.

**No absolute correlation threshold is shipped, because none works.** Over a
300-gene random null (`codependency_calibration_null.csv`) the median top-1 r for
an arbitrary screened gene is 0.270, with 88.7% reaching 0.20 and 37% reaching
0.30. True curated partners run down to r = 0.164
(`codependency_validation.csv`) while non-partners in the same top-15 slices
reach 0.550. The olfactory-receptor control OR5A1 returned OR5AN1 at r = 0.345,
above SMARCA4's best true partner at 0.303. `CODEP_STRONG_R = 0.2` is retained
only for the shared-partner view. Rank and `z` are the reportable statistics.

**The query-side gate is what discriminates.** `profile_sd >= 0.20` and
`min_effect <= -1.0` admits 9 of 9 complex-forming queries, 0 of 1 control, and
18.3% of random screened genes (`codependency_gate_grid.csv`). The grid also
shows the alternatives are worse: loosening `min_effect` to −1.5 or tightening
`profile_sd` to 0.25 each drop real-query admission below 1.0. Re-derive by
running `calibrate_codependency.py` then `gate_codependency.py`.

**The enrichment background is the screened genes, not the corpus.** GeneTEA
carries 35,818 genes against 18,531 screened, a 1.94x inflation
(`background_gene_counts.csv`). Across six queries the whole-corpus background
returns 427 terms against 276, dropping 185 — mostly generic screenability
vocabulary — and shifting FDR by a median +1.77 log10
(`background_comparison_summary.csv`). It also costs real power: VPS4A loses
`~ endosome` at corpus FDR 3.3e-3. Re-derive with `background_comparison.py`.

**Continuous mode needed its own threshold.** GeneTEA's default
`corr_thresh = 0.2` returns zero terms on any real DepMap vector: the maximum
attainable over all 24,831 terms is 0.104
(`continuous_term_correlations_TP53_top300.csv`, the 300 strongest of 24,831; regenerate the full table with `bench_continuous.py TP53`). `DEFAULT_CORR_THRESH = 0.03` yields
281 correlated terms and 74 enriched rows. Re-derive with
`bench_continuous.py TP53`, which also writes the timings in
`continuous_bench_TP53.json`.

`CODEP_MIN_LINES = 200` is not knife-edge: the per-gene valid-line count is
bimodal, with 17,931 genes screened in at least 100 lines (17,787 of them in at
least 1,000) and the remaining 600 in about 70 to 80, and none in between, so any
floor from 100 to 200 selects the same set.

## Scope

This skill does not score a single gene's dependency, classify it as
common-essential or selective, or contrast dependency by lineage, genotype or
fusion status — that is `depmap-local`. It does not join to human genetic
evidence or return a target verdict, which is `depmap-fusion`, and it does not
cover compound sensitivity or the wider triage sequence, which are
`depmap-local`'s PRISM helpers and `target-triage-public-data`. Codependency is
correlation across a panel of proliferating 2D cultures: it proposes partners and
does not establish physical interaction, and relationships absent from that panel
are invisible here. GeneTEA describes what has been written about genes, so it
inherits the literature's biases and says nothing about genes nobody has
characterised.
