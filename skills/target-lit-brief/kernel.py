"""Recent-literature brief for one drug target.

Takes a gene symbol, collects the most recently published PubMed papers that
are actually about that gene, and turns them into a per-paper table and a
triage-oriented synthesis.

Three properties of the PubMed connector shape this module, all measured
rather than assumed (see README ## Calibration):

- ``search_articles`` accepts at most 200 ids per call and
  ``get_article_metadata`` silently truncates to the first 20 -- it returns
  ``count: 20`` and HTTP success for a 100-id request, so an unchunked call
  drops the rest without raising.
- ``sort="pub_date"`` does not return records ordered by the
  ``publication_date`` field it reports. Over two genes the API's first 100
  hits missed 3-4 of the genuinely most recent 100 and the returned sequence
  held ~55 descending-order violations out of 199. Recency is therefore
  re-derived here from the parsed dates.
- ``publication_date.month`` arrives as both ``"08"`` and ``"Aug"`` in the
  same result set, so a string sort over the raw field is wrong.

A gene symbol is also not a search term. ``SET``, ``MAX``, ``REST`` and ``CS``
each match hundreds of thousands of papers using the word rather than the
gene, and no query-side qualifier separates them from clean symbols reliably
(measured: ``SRC`` and ``ATM`` straddle every threshold tried). On-target
filtering is therefore done against the abstract, after retrieval.
"""
import csv
import html
import json
import re

import pandas as pd

# PubMed caps, both observed from the live connector rather than documented.
SEARCH_MAX = 200
METADATA_CHUNK = 20

# Over-fetch this many candidates to rank down to the requested count. The
# default equals SEARCH_MAX: recency must be re-derived from parsed dates, and
# on-target screening removes some candidates, so the pool has to exceed the
# target count. Raising it needs paging (see plan_fetch).
DEFAULT_OVERFETCH = 200
DEFAULT_N_PAPERS = 100

MONTHS = {"jan": 1, "feb": 2, "mar": 3, "apr": 4, "may": 5, "jun": 6,
          "jul": 7, "aug": 8, "sep": 9, "oct": 10, "nov": 11, "dec": 12}

# Triage lens. Every paper is assigned exactly one, so the categories are
# ordered by what a target assessment acts on first.
TRIAGE_CATEGORIES = (
    "mechanism",       # molecular function, structure, interactors, pathway
    "dependency",      # knockout/knockdown/CRISPR, synthetic lethality
    "disease_link",    # human genetics, expression-outcome association
    "chemical_matter", # inhibitors, degraders, tool compounds, screens
    "biomarker",       # patient selection, stratification, response prediction
    "resistance",      # acquired or intrinsic resistance mechanisms
    "clinical",        # trials, patient cohorts, clinical outcomes
    "other",
)


def parse_pubmed_date(date_field):
    """PubMed ``publication_date`` to a sortable ``(year, month, day)``.

    The month arrives as either a zero-padded number or a three-letter
    abbreviation within a single result set, so both are handled. Absent
    components become 0, which sorts a month-only record before any dated
    record in the same month rather than after it -- deliberate, so an
    imprecise date never outranks a precise one for "most recent".
    """
    if not date_field:
        return (0, 0, 0)
    raw_year = date_field.get("year")
    year = int(raw_year) if str(raw_year or "").strip().isdigit() else 0

    raw_month = str(date_field.get("month") or "").strip().lower()
    if raw_month.isdigit():
        month = int(raw_month)
    else:
        month = MONTHS.get(raw_month[:3], 0)

    raw_day = str(date_field.get("day") or "").strip()
    day = int(raw_day) if raw_day.isdigit() else 0
    return (year, month, day)


def format_date(date_tuple):
    """``(2026, 8, 4)`` to ``"2026-08-04"``; unknown parts render as ``??``."""
    year, month, day = date_tuple
    ystr = f"{year:04d}" if year else "????"
    mstr = f"{month:02d}" if month else "??"
    dstr = f"{day:02d}" if day else "??"
    return f"{ystr}-{mstr}-{dstr}"


def pubmed_query(symbol, aliases=None, extra=None, restrict_to_title_abstract=True):
    """Build the PubMed query string for one gene symbol.

    Aliases are opt-in and off by default. Expanding them is not free: adding
    the single MyGene alias ``UBP`` to ``USP1`` takes the hit count from 373 to
    771, nearly all of the gain being unrelated papers, while ``WRN`` gains
    only 10 papers from three aliases. Pass them only when the symbol was
    recently renamed and the older name dominates the literature.

    ``extra`` is appended with AND and should already be valid PubMed syntax,
    e.g. ``'"Clinical Trial"[Publication Type]'``.
    """
    if not str(symbol or "").strip():
        raise ValueError("symbol is required and cannot be blank")
    symbol = str(symbol).strip()
    tag = "[Title/Abstract]" if restrict_to_title_abstract else ""

    terms = [symbol] + [a for a in (aliases or [])
                        if str(a).strip() and str(a).strip().upper() != symbol.upper()]
    clause = " OR ".join(f'"{t}"{tag}' for t in terms)
    if len(terms) > 1:
        clause = f"({clause})"
    if extra:
        clause = f"{clause} AND ({extra})"
    return clause


def plan_fetch(symbol, aliases=None, extra=None, n_papers=None, overfetch=None,
               date_from=None, date_to=None):
    """The parameters for the retrieval step, as a plain dict.

    Kept separate from the fetch itself because PubMed is reachable only from
    the control-plane kernel, which has no third-party packages: the plan is
    built here, written to a handoff file, and executed there. That keeps the
    query and paging logic in one place instead of duplicated into a cell.

    Raises when the requested count exceeds what one search call can supply,
    rather than silently returning fewer papers than asked for.
    """
    n_papers = DEFAULT_N_PAPERS if n_papers is None else int(n_papers)
    overfetch = DEFAULT_OVERFETCH if overfetch is None else int(overfetch)
    if n_papers < 1:
        raise ValueError(f"n_papers must be >= 1, got {n_papers}")
    if overfetch > SEARCH_MAX:
        raise ValueError(
            f"overfetch={overfetch} exceeds the connector's per-call limit of "
            f"{SEARCH_MAX}; page with retstart to go deeper")
    if overfetch < n_papers:
        raise ValueError(
            f"overfetch={overfetch} is below n_papers={n_papers}; the pool must "
            "exceed the target count because off-target papers get screened out")
    return {
        "symbol": symbol,
        "query": pubmed_query(symbol, aliases=aliases, extra=extra),
        "max_results": overfetch,
        "n_papers": n_papers,
        "metadata_chunk": METADATA_CHUNK,
        "sort": "pub_date",
        "date_from": date_from,
        "date_to": date_to,
    }


def clean_text(value):
    """Decode the HTML entities PubMed leaves in titles and abstracts.

    Greek letters and other symbols arrive as numeric character references
    (``&#x3b2;`` for beta), which would otherwise reach the CSV and the LLM
    prompts verbatim. Applied to titles and abstracts only.
    """
    if not value:
        return ""
    return html.unescape(str(value)).strip()


def flatten_article(article):
    """One connector article record to a flat dict.

    ``authors`` is collapsed to a first-author string plus a count, and the
    date to both a sortable tuple and a display string. The PMID is read from
    ``identifiers``, which is the only field guaranteed to carry it.
    """
    identifiers = article.get("identifiers") or {}
    journal = article.get("journal") or {}
    authors = article.get("authors") or []

    first_author = ""
    if authors:
        lead = authors[0] or {}
        last = str(lead.get("last_name") or "").strip()
        initials = str(lead.get("initials") or "").strip()
        first_author = f"{last} {initials}".strip()

    date_tuple = parse_pubmed_date(article.get("publication_date"))
    return {
        "pmid": str(identifiers.get("pmid") or ""),
        "date": format_date(date_tuple),
        "_date_key": date_tuple,
        "journal": journal.get("iso_abbreviation") or journal.get("title") or "",
        "title": clean_text(article.get("title")),
        "first_author": first_author,
        "n_authors": len(authors),
        "doi": identifiers.get("doi") or article.get("doi") or "",
        "pmcid": identifiers.get("pmc") or "",
        "article_types": "; ".join(article.get("article_types") or []),
        "abstract": clean_text(article.get("abstract")),
    }


def rank_recent(articles, n_papers=None):
    """Most recent ``n_papers`` records, newest first, de-duplicated by PMID.

    Ranking is done on the parsed date rather than the order the connector
    returned, because that order disagrees with the dates in the same payload.
    PMID descending breaks ties, so same-day papers order deterministically.
    Records whose date is entirely unparseable sort last but are not dropped.
    """
    n_papers = DEFAULT_N_PAPERS if n_papers is None else int(n_papers)
    rows = [flatten_article(a) for a in articles]

    seen, unique = set(), []
    for row in rows:
        if row["pmid"] and row["pmid"] in seen:
            continue
        if row["pmid"]:
            seen.add(row["pmid"])
        unique.append(row)

    def sort_key(row):
        pmid = row["pmid"]
        return (row["_date_key"], int(pmid) if pmid.isdigit() else 0)

    unique.sort(key=sort_key, reverse=True)
    return unique[:n_papers]


def _resolve_host_llm():
    """Find ``host.llm`` whether this module was sidecar-loaded or imported.

    Loading the skill mirrors these names into the session globals, where
    ``host`` is already bound; importing ``kernel`` directly (as the tests and
    any standalone script do) gives the module its own namespace, where it is
    not. Checking builtins covers the injected case without making ``host`` a
    module-level import that would fail at load time.
    """
    import builtins
    injected = getattr(builtins, "host", None) or globals().get("host")
    if injected is None or not hasattr(injected, "llm"):
        raise RuntimeError(
            "host.llm is unavailable; pass llm=<callable> explicitly. The "
            "callable takes a list of request dicts and returns a list of "
            "{'text': ...} in the same order.")
    return injected.llm


def _batch_llm(prompts, llm=None, max_tokens=400, model=None):
    """Fan a list of prompts out to the LLM, returning text or None per slot.

    A failed slot yields None rather than raising, so one malformed abstract
    cannot abort a 200-paper run. Callers must treat None as "not judged".
    """
    if not prompts:
        return []
    if llm is None:
        llm = _resolve_host_llm()
    requests = [{"prompt": p, "max_tokens": max_tokens} for p in prompts]
    if model:
        for request in requests:
            request["model"] = model
    replies = llm(requests, max_concurrency=8)
    out = []
    for reply in replies:
        if isinstance(reply, dict) and reply.get("text") is not None:
            out.append(reply["text"])
        else:
            out.append(None)
    return out


def screen_relevance(rows, symbol, gene_name="", llm=None):
    """Mark which rows are actually about the gene, judged from the abstract.

    Needed because a symbol that is also an English word retrieves mostly
    unrelated papers, and no query-side qualifier separates the two cases
    reliably. Sets ``on_target`` (bool) and ``screen_note`` on every row and
    returns the same list.

    A row the model did not judge is kept (``on_target=True``,
    ``screen_note="unscreened"``) so an LLM failure degrades to the unfiltered
    behaviour rather than silently emptying the brief.
    """
    if not rows:
        return rows
    label = f"{symbol} ({gene_name})" if gene_name else symbol
    prompts = [
        f"Gene symbol of interest: {label}\n\n"
        f"Title: {row['title']}\n"
        f"Abstract: {row['abstract'][:1500]}\n\n"
        "Is this paper about that specific gene or its protein product? "
        "Papers that merely use the symbol as an ordinary word, an acronym for "
        "something else, or a different species' unrelated gene are NOT about it.\n"
        "Answer with exactly one word: YES or NO."
        for row in rows
    ]
    for row, reply in zip(rows, _batch_llm(prompts, llm=llm, max_tokens=5)):
        if reply is None:
            row["on_target"] = True
            row["screen_note"] = "unscreened"
            continue
        verdict = re.search(r"\b(yes|no)\b", reply.strip().lower())
        if verdict is None:
            row["on_target"] = True
            row["screen_note"] = "unscreened"
        else:
            row["on_target"] = verdict.group(1) == "yes"
            row["screen_note"] = "screened"
    return rows


def _parse_summary_reply(text):
    """Pull ``{category, summary}`` out of a model reply, tolerating fences."""
    if not text:
        return None
    match = re.search(r"\{.*\}", text, flags=re.DOTALL)
    if not match:
        return None
    try:
        payload = json.loads(match.group(0))
    except (ValueError, TypeError):
        return None
    if not isinstance(payload, dict):
        return None
    category = str(payload.get("category") or "").strip().lower()
    if category not in TRIAGE_CATEGORIES:
        category = "other"
    return {"category": category, "summary": str(payload.get("summary") or "").strip()}


def summarize_papers(rows, symbol, gene_name="", llm=None, model=None):
    """Add a one-line ``summary`` and a ``category`` to each row.

    The summary is asked to state what the paper contributes about this target
    specifically, not what the paper is about in general -- a distinction that
    matters when the target appears in a supporting role. ``category`` is one
    of ``TRIAGE_CATEGORIES``; anything unrecognised becomes ``"other"`` so the
    column stays closed-vocabulary and groupable.
    """
    if not rows:
        return rows
    label = f"{symbol} ({gene_name})" if gene_name else symbol
    allowed = ", ".join(TRIAGE_CATEGORIES)
    prompts = [
        f"Target: {label}\n\n"
        f"Title: {row['title']}\n"
        f"Abstract: {row['abstract'][:2500]}\n\n"
        f"Classify this paper into exactly one category from: {allowed}\n"
        f"Then write ONE sentence (max 30 words) stating what this paper "
        f"contributes about {symbol} specifically -- the finding, not the topic. "
        f"If {symbol} is only mentioned in passing, say so.\n\n"
        'Reply with only JSON: {"category": "...", "summary": "..."}'
        for row in rows
    ]
    for row, reply in zip(rows, _batch_llm(prompts, llm=llm, max_tokens=250, model=model)):
        parsed = _parse_summary_reply(reply)
        if parsed is None:
            row["category"] = "other"
            row["summary"] = ""
            row["summary_note"] = "unsummarized"
        else:
            row["category"] = parsed["category"]
            row["summary"] = parsed["summary"]
            row["summary_note"] = "ok"
    return rows


def papers_frame(rows):
    """Rows to a DataFrame with the internal sort key dropped."""
    frame = pd.DataFrame(rows)
    return frame.drop(columns=[c for c in ["_date_key"] if c in frame.columns])


def write_paper_table(rows, path):
    """Write the per-paper table, abstracts excluded, and return the frame.

    Abstracts are dropped from the CSV because they are the bulk of the bytes
    and are reproducible from the PMID; the summary column carries what the
    table is read for.
    """
    frame = papers_frame(rows)
    frame = frame.drop(columns=[c for c in ["abstract"] if c in frame.columns])
    preferred = ["pmid", "date", "journal", "title", "first_author", "n_authors",
                 "category", "summary", "doi", "pmcid", "article_types",
                 "on_target", "screen_note", "summary_note"]
    ordered = [c for c in preferred if c in frame.columns]
    frame = frame[ordered + [c for c in frame.columns if c not in ordered]]
    frame.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)
    return frame


def coverage_stats(all_rows, kept_rows):
    """Retrieval and screening counts that the brief must report.

    ``date_span`` is taken from the kept rows, so it describes the papers the
    brief actually covers rather than everything retrieved.
    """
    screened = [r for r in all_rows if r.get("screen_note") == "screened"]
    off_target = [r for r in screened if not r.get("on_target", True)]
    dates = [r["date"] for r in kept_rows if not r["date"].startswith("????")]
    return {
        "n_retrieved": len(all_rows),
        "n_screened": len(screened),
        "n_off_target": len(off_target),
        "n_in_brief": len(kept_rows),
        "date_newest": max(dates) if dates else None,
        "date_oldest": min(dates) if dates else None,
        "n_unsummarized": sum(1 for r in kept_rows if r.get("summary_note") == "unsummarized"),
    }


def synthesis_prompt(rows, symbol, gene_name="", stats=None):
    """Build the prompt for the cross-paper synthesis.

    Papers are grouped by category and each is given with its PMID, so the
    model writes citations from the retrieved set instead of recalling them.
    """
    label = f"{symbol} ({gene_name})" if gene_name else symbol
    blocks = []
    for category in TRIAGE_CATEGORIES:
        members = [r for r in rows if r.get("category") == category]
        if not members:
            continue
        lines = "\n".join(
            f"- PMID {r['pmid']} ({r['date']}, {r['journal']}): {r.get('summary') or r['title']}"
            for r in members)
        blocks.append(f"## {category} ({len(members)} papers)\n{lines}")
    body = "\n\n".join(blocks)
    span = ""
    if stats and stats.get("date_oldest"):
        span = f" spanning {stats['date_oldest']} to {stats['date_newest']}"
    return (
        f"You are writing a target-triage literature brief on {label}, based on "
        f"the {len(rows)} most recent PubMed papers about it{span}.\n\n"
        f"{body}\n\n"
        "Write a brief in markdown with these sections:\n"
        "1. `## What changed` -- 3-5 sentences on the most consequential recent "
        "findings for this target. Lead with substance, not volume.\n"
        "2. `## By theme` -- one short subsection per category present, "
        "summarising what that body of work establishes.\n"
        "3. `## Open questions` -- 3-5 bullets naming what the recent "
        "literature does not settle.\n\n"
        "Cite every specific claim with its PMID inline, e.g. (PMID 12345678). "
        "Only cite PMIDs from the list above. Do not invent findings, and where "
        "the papers disagree, say so. Write descriptively; do not editorialise "
        "about how exciting or promising the target is."
    )


def write_brief(path, symbol, synthesis_text, stats, gene_name="", query=""):
    """Assemble the markdown brief: provenance header, synthesis, method note."""
    label = f"{symbol} -- {gene_name}" if gene_name else symbol
    span = "unknown"
    if stats.get("date_oldest"):
        span = f"{stats['date_oldest']} to {stats['date_newest']}"
    lines = [
        f"# {label}: recent literature brief",
        "",
        f"{stats['n_in_brief']} most recent PubMed papers, published {span}. "
        f"Source: PubMed (NCBI). Retrieved {stats['n_retrieved']} candidates; "
        f"{stats['n_off_target']} of {stats['n_screened']} screened were judged "
        "not to be about this gene and excluded.",
        "",
        synthesis_text.strip(),
        "",
        "## Method",
        "",
        f"- Query: `{query}`",
        f"- Ranked by parsed publication date, newest first; the connector's own "
        f"`sort=pub_date` order disagrees with the dates it returns and is not relied on.",
        f"- Off-target papers removed by abstract-level screening against the gene identity.",
    ]
    if stats.get("n_unsummarized"):
        lines.append(f"- {stats['n_unsummarized']} paper(s) could not be summarised "
                     "and carry an empty summary.")
    lines.append("")
    text = "\n".join(lines)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    return text
