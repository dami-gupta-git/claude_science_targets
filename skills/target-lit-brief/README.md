# target-lit-brief

Builds a recent-literature brief for a single drug target. Given a gene symbol,
it retrieves the most recently published PubMed papers about that gene, removes
the ones that only use the symbol as an ordinary word, summarises each paper
against a target-triage lens, and writes both a per-paper table and a synthesis
grouped by theme. It exists because the two obvious shortcuts for this task —
taking the first N results of a `sort=pub_date` search, and trusting a bare gene
symbol as a query — both return the wrong papers.

The retrieval itself runs against the `pubmed` connector, which is reachable
only from the control-plane kernel. The helpers here are pure: they build the
query, rank and screen the returned records, and render the outputs, so the
whole module is testable without network access.

## Functions

- `pubmed_query(symbol, aliases=None, extra=None, restrict_to_title_abstract=True)`
  — the PubMed query string for one symbol. Aliases are OR-grouped and
  deduplicated against the symbol; `extra` is ANDed as raw PubMed syntax.
  Raises on a blank symbol.
- `plan_fetch(symbol, aliases=None, extra=None, n_papers=100, overfetch=200, date_from=None, date_to=None)`
  — the retrieval parameters as a dict, for handing to the control-plane kernel.
  Raises when `overfetch` exceeds the connector's 200-per-call limit or falls
  below `n_papers`.
- `parse_pubmed_date(date_field)` — PubMed's `publication_date` to a sortable
  `(year, month, day)`, accepting both `"08"` and `"Aug"` as the month. Absent
  components become 0.
- `format_date(date_tuple)` — that tuple to `"2026-08-04"`, with `??` for
  unknown parts.
- `flatten_article(article)` — one connector record to a flat row: PMID, date,
  journal, title, lead author, author count, DOI, PMCID, article types,
  abstract.
- `rank_recent(articles, n_papers=100)` — the most recent papers, newest first,
  deduplicated by PMID and ranked on the parsed date rather than the order the
  connector returned. Undated records sort last but are retained.
- `screen_relevance(rows, symbol, gene_name="", llm=None)` — sets `on_target`
  and `screen_note` per row by judging each abstract against the gene identity.
  Rows the model could not judge are kept and marked `unscreened`.
- `summarize_papers(rows, symbol, gene_name="", llm=None, model=None)` — adds a
  one-line `summary` and a `category` drawn from `TRIAGE_CATEGORIES`.
  Unparseable replies yield an empty summary marked `unsummarized`, never an
  invented one.
- `papers_frame(rows)` / `write_paper_table(rows, path)` — the per-paper table;
  the writer drops abstracts and orders PMID and date first.
- `coverage_stats(all_rows, kept_rows)` — retrieved, screened, excluded and
  covered counts plus the date span of the kept papers.
- `synthesis_prompt(rows, symbol, gene_name="", stats=None)` — the cross-paper
  prompt, grouped by category with every PMID listed so citations come from the
  retrieved set.
- `write_brief(path, symbol, synthesis_text, stats, gene_name="", query="")` —
  the markdown brief with its provenance header and method note.

Module constants: `SEARCH_MAX`, `METADATA_CHUNK`, `DEFAULT_OVERFETCH`,
`DEFAULT_N_PAPERS`, `TRIAGE_CATEGORIES`, `MONTHS`.

## Calibration

Every limit below was measured against the live `pubmed` connector in August
2026, not taken from documentation. The measurements are reproducible by
re-running the queries named.

`METADATA_CHUNK = 20`. `get_article_metadata` returns at most 20 records and
reports success for larger requests: asking for 21, 25 and 50 PMIDs each
returned `count: 20`, with the surplus dropped and no error raised.

`SEARCH_MAX = 200`. `search_articles` accepts `max_results` up to 200 and
returns `INVALID_PARAMETERS` at 300.

`DEFAULT_OVERFETCH = 200`. Ranking cannot rely on the connector's own ordering.
For `USP1[Title/Abstract]` and `WRN[Title/Abstract]`, the API's first 100 hits
under `sort="pub_date"` contained only 96 and 97 of the 100 most recent papers
by the `publication_date` field in the same payload, and the returned sequences
held 56 and 53 descending-order violations across 199 adjacent pairs. Over-
fetching the full 200 and re-ranking recovers the correct set. Re-derive by
fetching 200 for either query, parsing the dates, and comparing the date-ranked
top 100 against the first 100 ids returned.

Query-side disambiguation was tested and rejected. Adding a gene-context
qualifier (`genes`/`proteins` MeSH plus gene, protein, expression, kinase and
mutation in title/abstract) reduced hit counts by 6–36% across 14 unambiguous
symbols and 53–83% across 12 word-like symbols — but the ranges overlap, with
`SRC` at 9% and `ATM` at 36% against `AR` at 53%, so no threshold separates the
two groups. Abstract-level screening replaced it. The measurement is in
`handoff/ambiguity_calibration.json` of the authoring session and reproducible
by running both query forms for each symbol and comparing `total_count`.

Alias expansion is off by default on the same basis: `USP1` plus its MyGene
alias `UBP` goes from 373 to 771 hits, while `WRN` plus three aliases gains 10
and `CDK12` plus three gains 11.

## Tests

38 tests, no network and no LLM — the connector payload shape is synthesised
from live responses and the LLM is injected as a stub.

```
cd skills/target-lit-brief && python -m pytest
```

The suite was checked against deliberate inversions of the behaviours it
constrains, confirming each is detected rather than merely documented. A
mutation may fail more than one test, since several tests share a behaviour:

| Inversion applied to `kernel.py` | Tests that fail |
| --- | --- |
| `rank_recent` no longer sorts by date | 3 (`ranking_ignores_api_order`, `undated_records_sort_last`, `same_day_papers_order_deterministically`) |
| abbreviated months treated as unparseable | 1 (`numeric_and_abbreviated_months_agree`) |
| verdict matched as a substring, unbounded | 2 (`verdict_matching_is_word_bounded`, `unparseable_verdict_is_treated_as_unscreened`) |
| unscreened papers dropped instead of kept | 2 (`llm_failure_keeps_paper_and_flags_it`, `coverage_counts_only_screened_papers_as_off_target`) |
| overfetch cap guard removed | 1 (`overfetch_above_connector_cap_is_rejected`) |

Reproduce by applying one inversion at a time and running the suite. An
inversion that fails nothing would mean the named behaviour is unconstrained;
none of the five is.

## Scope

This skill covers recent publications for one gene. It does not judge whether
the target is worth pursuing — `target-triage-public-data` and `depmap-fusion`
handle dependency and human-genetic evidence, `opentargets-evidence` handles
tractability and known drugs. Topic-level reviews, citation graphs and preprints
belong to `literature-review`. Full text is never fetched: summaries are built
from titles and abstracts, so a paper whose contribution appears only in its
results section will be summarised from what the abstract claims.
