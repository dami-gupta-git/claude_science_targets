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
numbered references at the end.

The retrieval itself runs against the `pubmed` connector, which is reachable
only from the control-plane kernel. The helpers here build the query, rank and
screen the returned records, and render the outputs; none of them make a
network call directly, so the whole module is testable without network access.

## Functions

| Function | Description |
| --- | --- |
| `pubmed_query(symbol, aliases=None, extra=None)` | the PubMed query string for one symbol, with aliases OR-grouped. |
| `plan_fetch(symbol, ..., n_papers=None, overfetch=None)` | the retrieval parameters as a dict, for handing to the control-plane kernel. |
| `parse_pubmed_date(date_field)` | PubMed's `publication_date` to a sortable `(year, month, day)`. |
| `format_date(date_tuple)` | that tuple to `"2026-08-04"`, with `??` for unknown parts. |
| `clean_text(value)` | decodes the HTML entities PubMed leaves in titles and abstracts. |
| `flatten_article(article)` | one connector record to a flat row of PMID, date, journal, title, authors, ids and abstract. |
| `rank_recent(articles, n_papers=None)` | papers newest first, deduplicated by PMID and ranked on the parsed date. |
| `screen_relevance(rows, symbol, ..., target_n=None)` | sets `on_target` and `screen_note` per row, screening in batches until `target_n` papers are found. |
| `select_for_brief(rows, n_papers=None)` | the on-target, screened rows that go into the brief, newest first. |
| `summarize_papers(rows, symbol, ...)` | adds a one-line `summary` and a `category` drawn from `TRIAGE_CATEGORIES`. |
| `papers_frame(rows)` / `write_paper_table(rows, path)` | the per-paper table, and the CSV writer for it. |
| `coverage_stats(all_rows, kept_rows)` | retrieved, screened, excluded and covered counts plus the date span. |
| `synthesis_prompt(rows, symbol, ..., stats=None)` | the cross-paper prompt, grouped by category with every PMID listed. |
| `renumber_citations(text, rows)` | inline `(PMID 12345678)` citations to `[1]` markers, numbered by first appearance. |
| `format_references(references)` | the reference list as markdown, with a PubMed link and DOI per entry. |
| `write_brief(path, symbol, synthesis_text, stats, ...)` | the markdown brief with its provenance header, references and method note. |
| `brief_run_dir(symbol)` / `brief_write_run(out_dir, ...)` | the run directory under `results/target_lit_brief/`, and both output files written into it. |

Module constants: `SEARCH_MAX`, `METADATA_CHUNK`, `DEFAULT_OVERFETCH`,
`DEFAULT_N_PAPERS`, `SCREEN_BATCH_FACTOR`, `SCREEN_BATCH_MIN`,
`TRIAGE_CATEGORIES`, `MONTHS`, `CITATION_GROUP`, `CITATION_SINGLE`.

## Tests

No network and no LLM — the connector payload shape is synthesised
from live responses and the LLM is injected as a stub.

```
cd skills/target-lit-brief && python -m pytest
```

## Scope

This skill covers recent publications for one gene. It does not judge whether
the target is worth pursuing — `target-triage-public-data` and `depmap-fusion`
handle dependency and human-genetic evidence, `opentargets-evidence` handles
tractability and known drugs. Topic-level reviews, citation graphs and preprints
belong to `literature-review`. Full text is never fetched: summaries are built
from titles and abstracts, so a paper whose contribution appears only in its
results section will be summarised from what the abstract claims.
