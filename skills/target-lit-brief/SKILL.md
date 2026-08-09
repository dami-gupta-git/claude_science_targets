---
name: target-lit-brief
description: Pull the most recent PubMed papers on one drug target and turn them into a per-paper table plus a triage-oriented synthesis brief. Use when asked what is new on a gene or target, to catch up on recent literature for a specific gene, to build the literature section of a target dossier, or when a target assessment needs current publications rather than a general topic search. Handles gene symbols that collide with ordinary words (SET, MAX, REST, CS).
---

# Recent-literature brief for a target

Collects the most recently published PubMed papers about one gene, screens out
the ones that merely use the symbol as a word, summarises each against a
target-triage lens, and writes a synthesis. Output is a CSV of papers and a
markdown brief.

The default is 5 papers, set by `DEFAULT_N_PAPERS` and overridable per run with
`n_papers=`. It is sized for "what is new on this target", which is the common
question; a target dossier or a first look at an unfamiliar gene warrants more.
The candidate pool stays at 200 whatever the count, because the connector's
ordering cannot be trusted — see the recency section below.

## Workflow

1. **Resolve the symbol.** `query_genes` on the `genes-ontologies` connector
   (`terms=[symbol], scopes="symbol", fields="symbol,name,alias"`) confirms the
   gene exists and gives the full name used in the prompts. A symbol that
   returns no record is a typo or a non-human gene — stop and say so rather
   than searching for it anyway.
2. **Plan the fetch.** `plan_fetch(symbol, n_papers=5)` returns the query and
   paging parameters. Write it to `handoff/plan.json`.
3. **Retrieve** in the `repl` tool, since PubMed is reachable only there. Read
   the plan, call `search_articles`, then `get_article_metadata` **in chunks of
   20** — see the caps below — and write the raw articles to
   `handoff/articles.json`.
4. **Rank and screen.** `rank_recent(articles)` orders the whole pool by parsed
   date. `screen_relevance(rows, symbol, gene_name, target_n=n_papers)` then
   screens down that list and stops once `n_papers` on-target papers are found,
   so a 5-paper brief costs about 20 screening calls rather than 200.
   `select_for_brief(rows, n_papers)` returns the papers to write up.
5. **Summarise.** `summarize_papers(rows, symbol, gene_name)` attaches a
   one-line `summary` and a `category` to each paper. Run it on the selected
   rows only — summarising the whole pool pays for papers no one will read.
6. **Synthesise.** `synthesis_prompt(rows, symbol, gene_name, stats)` builds the
   prompt; pass it through `host.llm(..., model=host.reasoning_model())`.
7. **Write.** `out_dir = brief_run_dir(symbol)`, then `brief_write_run(out_dir,
   symbol, rows, synthesis_text, stats, gene_name, query)` — not
   `write_paper_table`/`write_brief` called by hand with a hand-built path.
   Lands in `results/target_lit_brief/<slug>/` as `<slug>_recent_papers.csv`
   and `README.md`, matching the existing runs — `<slug>` is the symbol
   lowercased (`CDK12` writes to `cdk12/`), not the case it was given in.
   `brief_run_dir` needs
   `$SCIENCE_RESULTS_ROOT` set to this repo's `results/` directory (or an
   explicit `root=`) — it raises naming the variable rather than silently
   creating a `results/` folder wherever the session's cwd happens to be.
   `brief_write_run` passes `rows=` through to `write_brief` internally,
   which is what converts inline PMIDs to numbered references — see below.

Steps 4 and 5 both cost one LLM call per paper, so the count drives the cost of
a run. Retrieval does not: search and metadata are flat in the pool size.

## Connector limits that bite

Both were measured against the live connector, not read from documentation.

`get_article_metadata` **silently truncates to the first 20 PMIDs**. A 100-id
request returns HTTP success and `count: 20`; the other 80 vanish with no error.
Always chunk at `METADATA_CHUNK`:

```python
articles = []
for i in range(0, len(pmids), 20):
    articles.extend(host.mcp("pubmed", "get_article_metadata",
                             pmids=pmids[i:i + 20]).get("articles", []))
assert len(articles) >= 0.9 * len(pmids), f"got {len(articles)} of {len(pmids)}"
```

`search_articles` caps `max_results` at 200 and rejects more with
`INVALID_PARAMETERS`. For more than ~180 papers after screening, page with
`retstart`.

## Recency is not what the API returns

`sort="pub_date"` does **not** order records by the `publication_date` field in
the same payload. Measured on two genes, the API's first 100 hits missed 3–4 of
the genuinely most recent 100, and the returned sequence held roughly 55
descending-order violations out of 199. `rank_recent` therefore over-fetches 200
candidates and re-ranks on the parsed date. Do not shortcut this by taking the
first N ids from the search result.

The `month` field arrives as both `"08"` and `"Aug"` **within one result set**,
so any string comparison over the raw date is wrong — `"12" < "Aug"`.
`parse_pubmed_date` handles both.

## Gene symbols that are also words

`SET`, `MAX`, `REST`, `CS`, `IMPACT` and `SHE` each retrieve hundreds of
thousands of papers that have nothing to do with the gene. Adding a
gene-context qualifier to the query does not separate these cleanly: the
reduction in hit count it produces overlaps between clean and ambiguous symbols
(`SRC` 9%, `ATM` 36%, `AR` 53%), so no threshold on it is safe.

Filtering is done **after** retrieval instead, by `screen_relevance`, which
judges each abstract against the gene identity. This costs one cheap LLM call
per paper and is reliable because the abstract carries the disambiguating
context the query lacks. When a symbol is ambiguous, over-fetch the full 200 —
screening may discard most of them.

A paper the screen could not judge is **kept** and marked `unscreened`, so an
LLM failure degrades to unfiltered output rather than a silently empty brief.
The counts of screened and excluded papers are reported in the brief.

## Aliases are opt-in

`pubmed_query` takes aliases but does not add them by default. Expanding `USP1`
with its single MyGene alias `UBP` takes the result set from 373 to 771, almost
all of the gain unrelated. `WRN` gains 10 papers from three aliases. Pass
aliases only when the gene was recently renamed and the older symbol dominates
the literature.

## The brief is written to be read

The reader is a scientist who knows basic genetics and drug discovery but does
not work on this target. `synthesis_prompt` states that in the prompt and
constrains the shape of the output:

- `## Overview` first — a few sentences of plain prose, no citations, no
  bullets. Someone who reads only this section should still know where the
  target stands.
- `## What changed` and each `## By theme` subsection lead with a plain
  sentence, then bullets. A theme written as one long paragraph is the failure
  mode this replaced; the earlier version produced 400-word blocks carrying
  twenty citations each.
- Sentences under ~25 words, no paragraph over four sentences, at most three
  citations per bullet, specialist terms glossed on first use.

Check the rendered brief against those constraints rather than assuming the
model honoured them. A theme section that came back as a single paragraph
should be sent back, not published.

## Citations become numbered references

Pass `rows=` to `write_brief` and inline `(PMID 12345678)` citations are
converted to `[1]` markers with a `## References` list at the end, each entry
linking to PubMed. Eight-digit ids inline are most of what makes the prose hard
to scan.

The numbering is done in `renumber_citations`, not by the model: a model
numbering its own references has to hold the whole mapping while writing prose,
and a mis-numbered citation points at the wrong paper without looking wrong.

A cited PMID that is not in the retrieved set is left inline unconverted and
named in the method note. The synthesis prompt forbids outside citations, so
one appearing means the model invented it — treat the claim carrying it as
unsupported rather than editing the number out.

## Reporting

The brief must state the number of papers covered, their date span, and how
many were excluded by screening — `coverage_stats` returns all of these and
`write_brief` renders them. PubMed is cited as the source and DOIs are carried
in the table.

Do not describe a target as promising or attractive. Report what the papers
establish and where they disagree.

## Scope

This skill covers recent literature for one gene. It does not assess whether the
target is worth pursuing — `target-triage-public-data` and `depmap-fusion` cover
dependency and human-genetic evidence, and `opentargets-evidence` covers
tractability and known drugs. For a full review of a topic rather than a target,
including citation-graph work and preprints, use `literature-review`. Full text
is not retrieved; summaries are built from titles and abstracts.
