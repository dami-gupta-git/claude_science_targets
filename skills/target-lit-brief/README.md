# target-lit-brief

Builds a recent-literature brief for a single drug target. Given a gene symbol,
it retrieves the most recently published PubMed papers about that gene, removes
the ones that only use the symbol as an ordinary word, summarises each paper
against a target-triage lens, and writes both a per-paper table and a synthesis
grouped by theme. It exists because the two obvious shortcuts for this task —
taking the first N results of a `sort=pub_date` search, and trusting a bare gene
symbol as a query — both return the wrong papers.

The brief is written for a scientist who knows basic genetics and drug
discovery but does not work on the target: a plain-prose overview first, then
themed sections that open with a plain sentence and continue as bullets, and
numbered references at the end rather than PubMed ids inline.

The retrieval itself runs against the `pubmed` connector, which is reachable
only from the control-plane kernel. The helpers here are pure: they build the
query, rank and screen the returned records, and render the outputs, so the
whole module is testable without network access.

## Functions

- `pubmed_query(symbol, aliases=None, extra=None, restrict_to_title_abstract=True)`
  — the PubMed query string for one symbol. Aliases are OR-grouped and
  deduplicated against the symbol; `extra` is ANDed as raw PubMed syntax.
  Raises on a blank symbol.
- `plan_fetch(symbol, aliases=None, extra=None, n_papers=5, overfetch=200, date_from=None, date_to=None)`
  — the retrieval parameters as a dict, for handing to the control-plane kernel.
  Raises when `overfetch` exceeds the connector's 200-per-call limit or falls
  below `n_papers`.
- `parse_pubmed_date(date_field)` — PubMed's `publication_date` to a sortable
  `(year, month, day)`, accepting both `"08"` and `"Aug"` as the month. Absent
  components become 0.
- `format_date(date_tuple)` — that tuple to `"2026-08-04"`, with `??` for
  unknown parts.
- `clean_text(value)` — decodes the HTML entities PubMed leaves in titles and
  abstracts, where Greek letters arrive as numeric character references.
- `flatten_article(article)` — one connector record to a flat row: PMID, date,
  journal, title, lead author, author count, DOI, PMCID, article types,
  abstract.
- `rank_recent(articles, n_papers=None)` — papers newest first, deduplicated by
  PMID and ranked on the parsed date rather than the order the connector
  returned. Undated records sort last but are retained. The default keeps the
  whole pool, since screening runs down it.
- `screen_relevance(rows, symbol, gene_name="", llm=None, target_n=None, batch_factor=2, batch_min=20)`
  — sets `on_target` and `screen_note` per row by judging each abstract against
  the gene identity. With `target_n` it screens in batches and stops once that
  many on-target papers are found, leaving the rest of the pool unjudged. Rows
  the model could not judge are kept and marked `unscreened`.
- `select_for_brief(rows, n_papers=5)` — the on-target, actually-screened rows
  that go into the brief, newest first. Truncation happens here rather than
  before screening, so off-target papers near the top of the ranking do not
  reduce the count.
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
- `renumber_citations(text, rows)` — inline `(PMID 12345678)` citations to `[1]`
  markers, numbered by first appearance. Returns the rewritten text, the
  reference entries, and any cited PMID absent from `rows`.
- `format_references(references)` — the reference list as markdown, one
  numbered line per paper with a PubMed link and the DOI where present.
- `write_brief(path, symbol, synthesis_text, stats, gene_name="", query="", rows=None)`
  — the markdown brief with its provenance header, references and method note.
  Passing `rows` converts the citations; omitting it leaves them inline.

Module constants: `SEARCH_MAX`, `METADATA_CHUNK`, `DEFAULT_OVERFETCH`,
`DEFAULT_N_PAPERS`, `SCREEN_BATCH_FACTOR`, `SCREEN_BATCH_MIN`,
`TRIAGE_CATEGORIES`, `MONTHS`, `CITATION_GROUP`, `CITATION_SINGLE`.

### Saving a run

- `brief_run_dir(symbol, root="results", topic="target_lit_brief")` — the
  path for one brief, `<root>/<topic>/<slug>/`, with `scripts/` created
  beside it.
- `brief_write_run(out_dir, symbol, kept_rows, synthesis_text, stats,
  gene_name="", query="", scripts=())` — wraps `write_paper_table` and
  `write_brief` rather than replacing them: writes
  `<symbol>_recent_papers.csv` and `README.md` into `out_dir`, matching the
  file names the existing runs already use, and copies `scripts` into
  `scripts/`. Returns the paths written.

Runs land in `results/target_lit_brief/<slug>/`, the topic the existing USP1
and WRN runs already use.

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

`DEFAULT_N_PAPERS = 5` is an editorial default, not a measurement: it is sized
for "what is new on this target", and any run can override it with `n_papers=`.
The two earlier runs in `results/target_lit_brief` used 95 and 100 papers and
produced theme sections of 400-plus words carrying twenty citations each, which
is what the current prompt structure and this default are set against.

`SCREEN_BATCH_FACTOR = 2` with `SCREEN_BATCH_MIN = 20` sets how far past the
requested count screening goes before giving up on filling it. The measured
off-target rate was 21% for `USP1` (25 of 118 screened) and 33% for `WRN` (61
of 186), so screening twice the quota fills it in a single batch at both rates,
and the floor keeps a small run to one batch. Re-derive from the
`n_off_target` and `n_screened` counts that `coverage_stats` reports in any
run's brief. A symbol that is also an ordinary word discards more and warrants
a higher factor.

## Tests

63 tests, no network and no LLM — the connector payload shape is synthesised
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
| `rank_recent` truncates before screening | 1 (`ranking_keeps_the_whole_pool_by_default`) |
| screening never stops early | 2 (`screening_stops_once_the_quota_is_filled`, `screening_continues_when_the_first_batch_falls_short`) |
| screening stops after one batch regardless of quota | 1 (`screening_continues_when_the_first_batch_falls_short`) |
| unscreened rows admitted to the brief | 1 (`selection_ignores_an_on_target_flag_this_screen_did_not_set`) |
| a cited PMID outside the retrieved set dropped silently | 1 (`pmid_outside_the_retrieved_set_is_reported_not_renumbered`) |
| citations numbered by count rather than first appearance | 1 (`inline_pmids_become_numbered_markers_in_order_of_appearance`) |
| single-citation second pass removed | 1 (`pmids_separated_by_prose_inside_one_parenthetical_are_converted`) |

Reproduce by applying one inversion at a time and running the suite. An
inversion that fails nothing would mean the named behaviour is unconstrained;
none of these is.

The citation regex is also checked against the two briefs already in
`results/target_lit_brief`, which between them carry 260 inline citations in
four formats — `(PMID n)`, a line-wrapped variant, `(PMID n, PMID m)`, and
`(PMID n, correction at PMID m)`. All convert, resolving to 95 and 100 distinct
references with none cited from outside the retrieved set and none left inline.
The prose phrase "one row per paper: PMID, date, journal" correctly does not
convert.

## Scope

This skill covers recent publications for one gene. It does not judge whether
the target is worth pursuing — `target-triage-public-data` and `depmap-fusion`
handle dependency and human-genetic evidence, `opentargets-evidence` handles
tractability and known drugs. Topic-level reviews, citation graphs and preprints
belong to `literature-review`. Full text is never fetched: summaries are built
from titles and abstracts, so a paper whose contribution appears only in its
results section will be summarised from what the abstract claims.
