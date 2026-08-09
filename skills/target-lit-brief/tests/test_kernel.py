"""Tests for target-lit-brief.

No network and no LLM: the connector payload shape is synthesised from what
the live API actually returned during calibration (mixed month formats, the
20-record metadata truncation, non-monotonic sort order), and the LLM is
injected as a stub.
"""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import kernel as k  # noqa: E402


# --------------------------------------------------------------------------
# date parsing -- the mixed-format field
# --------------------------------------------------------------------------

def test_numeric_and_abbreviated_months_agree():
    """'08' and 'Aug' are both August; PubMed emits both in one result set."""
    assert k.parse_pubmed_date({"year": "2026", "month": "08", "day": "04"}) == (2026, 8, 4)
    assert k.parse_pubmed_date({"year": "2026", "month": "Aug", "day": "04"}) == (2026, 8, 4)


def test_string_sort_over_raw_field_would_be_wrong():
    """The bug this parser exists to prevent: 'Aug' < '12' as strings."""
    numeric = k.parse_pubmed_date({"year": "2026", "month": "12", "day": "01"})
    abbrev = k.parse_pubmed_date({"year": "2026", "month": "Aug", "day": "01"})
    assert numeric > abbrev          # December after August, correctly
    assert "12" < "Aug"              # but the raw strings sort the other way


def test_missing_components_sort_oldest_within_period():
    assert k.parse_pubmed_date({"year": "2026"}) == (2026, 0, 0)
    assert k.parse_pubmed_date({}) == (0, 0, 0)
    assert k.parse_pubmed_date(None) == (0, 0, 0)
    # A month-only record must not outrank a fully dated one in that month.
    assert k.parse_pubmed_date({"year": "2026", "month": "06"}) < \
           k.parse_pubmed_date({"year": "2026", "month": "06", "day": "01"})


def test_garbage_components_do_not_raise():
    assert k.parse_pubmed_date({"year": "n/a", "month": "Smarch", "day": "??"}) == (0, 0, 0)


def test_format_date_marks_unknown_parts():
    assert k.format_date((2026, 8, 4)) == "2026-08-04"
    assert k.format_date((2026, 0, 0)) == "2026-??-??"
    assert k.format_date((0, 0, 0)) == "????-??-??"


# --------------------------------------------------------------------------
# query construction
# --------------------------------------------------------------------------

def test_plain_symbol_query_is_title_abstract_scoped():
    assert k.pubmed_query("USP1") == '"USP1"[Title/Abstract]'


def test_aliases_are_or_grouped_and_deduplicated():
    q = k.pubmed_query("WRN", aliases=["RECQL2", "wrn", "  "])
    assert q == '("WRN"[Title/Abstract] OR "RECQL2"[Title/Abstract])'


def test_extra_clause_is_anded():
    q = k.pubmed_query("USP1", extra='"Clinical Trial"[Publication Type]')
    assert q == '"USP1"[Title/Abstract] AND ("Clinical Trial"[Publication Type])'


def test_blank_symbol_is_rejected():
    for bad in ["", "   ", None]:
        with pytest.raises(ValueError, match="symbol is required"):
            k.pubmed_query(bad)


# --------------------------------------------------------------------------
# fetch planning -- the connector caps
# --------------------------------------------------------------------------

def test_plan_defaults_overfetch_above_target_count():
    plan = k.plan_fetch("USP1")
    assert plan["n_papers"] == k.DEFAULT_N_PAPERS == 5
    assert plan["max_results"] == 200
    assert plan["metadata_chunk"] == 20      # the silent-truncation limit


def test_pool_stays_at_the_cap_for_a_small_brief():
    """A 5-paper brief still over-fetches 200: the API's ordering is unreliable
    across the whole result set, not just its tail, so the top 5 by date are
    not recoverable from a small fetch."""
    assert k.plan_fetch("USP1", n_papers=5)["max_results"] == 200


def test_requested_count_is_carried_through():
    assert k.plan_fetch("USP1", n_papers=25)["n_papers"] == 25


def test_overfetch_above_connector_cap_is_rejected():
    with pytest.raises(ValueError, match="exceeds the connector's per-call limit"):
        k.plan_fetch("USP1", overfetch=300)


def test_overfetch_below_target_is_rejected():
    """Screening removes papers, so the pool must exceed what is asked for."""
    with pytest.raises(ValueError, match="below n_papers"):
        k.plan_fetch("USP1", n_papers=150, overfetch=120)


def test_nonsense_paper_count_is_rejected():
    with pytest.raises(ValueError, match="n_papers must be"):
        k.plan_fetch("USP1", n_papers=0)


# --------------------------------------------------------------------------
# flattening and ranking
# --------------------------------------------------------------------------

def _article(pmid, year, month, day, title="T", abstract="A", authors=2):
    return {
        "identifiers": {"pmid": pmid, "doi": f"10.1/{pmid}", "pmc": f"PMC{pmid}"},
        "title": title,
        "abstract": abstract,
        "journal": {"title": "Journal of Things", "iso_abbreviation": "J Things"},
        "authors": [{"last_name": f"A{i}", "fore_name": "X", "initials": "X"}
                    for i in range(authors)],
        "publication_date": {"year": year, "month": month, "day": day},
        "article_types": ["Journal Article"],
    }


def test_flatten_reads_identifiers_and_lead_author():
    row = k.flatten_article(_article("111", "2026", "08", "04"))
    assert row["pmid"] == "111"
    assert row["date"] == "2026-08-04"
    assert row["journal"] == "J Things"
    assert row["first_author"] == "A0 X"
    assert row["n_authors"] == 2
    assert row["doi"] == "10.1/111"


def test_html_entities_are_decoded_in_title_and_abstract():
    """PubMed returns Greek letters as numeric character references."""
    article = _article("113", "2026", "01", "01",
                       title="USP18-driven macrophage &#x3b2;-oxidation",
                       abstract="TGF-&beta; signalling &amp; USP1")
    row = k.flatten_article(article)
    assert row["title"] == "USP18-driven macrophage \u03b2-oxidation"
    assert row["abstract"] == "TGF-\u03b2 signalling & USP1"


def test_clean_text_handles_empty_input():
    assert k.clean_text(None) == ""
    assert k.clean_text("") == ""


def test_flatten_survives_empty_author_list():
    article = _article("112", "2026", "01", "01", authors=0)
    row = k.flatten_article(article)
    assert row["first_author"] == ""
    assert row["n_authors"] == 0


def test_ranking_ignores_api_order_and_uses_dates():
    """The API returns these out of order; the ranker must not trust that."""
    articles = [
        _article("1", "2024", "01", "15"),
        _article("2", "2026", "08", "04"),
        _article("3", "2025", "Nov", "02"),   # abbreviated month
        _article("4", "2026", "01", "20"),
    ]
    ranked = k.rank_recent(articles, n_papers=4)
    assert [r["pmid"] for r in ranked] == ["2", "4", "3", "1"]


def test_ranking_truncates_to_requested_count():
    articles = [_article(str(i), "2026", "01", f"{i:02d}") for i in range(1, 21)]
    assert len(k.rank_recent(articles, n_papers=5)) == 5


def test_ranking_keeps_the_whole_pool_by_default():
    """Truncating before screening would leave the brief short by however many
    of the top rows turn out to be off-target."""
    articles = [_article(str(i), "2026", "01", f"{i:02d}") for i in range(1, 21)]
    assert len(k.rank_recent(articles)) == 20


def test_duplicate_pmids_collapse():
    articles = [_article("9", "2026", "01", "01"), _article("9", "2026", "01", "01")]
    assert len(k.rank_recent(articles, n_papers=10)) == 1


def test_undated_records_sort_last_but_survive():
    articles = [
        {"identifiers": {"pmid": "5"}, "title": "no date", "publication_date": None},
        _article("6", "2020", "01", "01"),
    ]
    ranked = k.rank_recent(articles, n_papers=10)
    assert [r["pmid"] for r in ranked] == ["6", "5"]


def test_same_day_papers_order_deterministically():
    articles = [_article("100", "2026", "05", "05"), _article("200", "2026", "05", "05")]
    first = [r["pmid"] for r in k.rank_recent(articles, n_papers=2)]
    second = [r["pmid"] for r in k.rank_recent(list(reversed(articles)), n_papers=2)]
    assert first == second == ["200", "100"]


# --------------------------------------------------------------------------
# screening -- the symbol-collision guard
# --------------------------------------------------------------------------

def _stub_llm(replies):
    """Return a callable matching host.llm's list-in/list-out contract."""
    def fake(requests, max_concurrency=8):
        return [{"text": r} if r is not None else {"error": "boom"}
                for r in replies[:len(requests)]]
    return fake


def _rows(n):
    return k.rank_recent([_article(str(i), "2026", "01", f"{i:02d}")
                          for i in range(1, n + 1)], n_papers=n)


def test_screening_marks_off_target_papers():
    rows = _rows(3)
    k.screen_relevance(rows, "SET", llm=_stub_llm(["YES", "NO", "yes"]))
    assert [r["on_target"] for r in rows] == [True, False, True]
    assert all(r["screen_note"] == "screened" for r in rows)


def test_llm_failure_keeps_paper_and_flags_it():
    """A failed screen must not silently delete papers from the brief."""
    rows = _rows(2)
    k.screen_relevance(rows, "SET", llm=_stub_llm([None, "NO"]))
    assert rows[0]["on_target"] is True
    assert rows[0]["screen_note"] == "unscreened"
    assert rows[1]["on_target"] is False


def test_unparseable_verdict_is_treated_as_unscreened():
    rows = _rows(1)
    k.screen_relevance(rows, "SET", llm=_stub_llm(["I am not sure about this"]))
    assert rows[0]["on_target"] is True
    assert rows[0]["screen_note"] == "unscreened"


def test_verdict_matching_is_word_bounded():
    """'no' inside 'Notably' must not be read as a NO verdict.

    The distractor is placed BEFORE the real verdict, so an unbounded
    substring match would return the wrong answer rather than the right one
    by luck.
    """
    rows = _rows(2)
    k.screen_relevance(rows, "SET", llm=_stub_llm(["Notably, YES", "Unknown gene: NO"]))
    assert rows[0]["on_target"] is True
    assert rows[1]["on_target"] is False


def test_screening_empty_input_is_a_noop():
    assert k.screen_relevance([], "SET", llm=_stub_llm([])) == []


def _counting_llm(replies):
    """Stub that records how many prompts it was asked to judge."""
    calls = {"n": 0}

    def fake(requests, max_concurrency=8):
        start = calls["n"]
        calls["n"] += len(requests)
        return [{"text": r} for r in replies[start:start + len(requests)]]
    return fake, calls


def test_screening_stops_once_the_quota_is_filled():
    """The pool is 200 for ranking reasons, but a 5-paper brief must not pay
    200 screening calls to fill it."""
    rows = _rows(60)
    llm, calls = _counting_llm(["YES"] * 60)
    k.screen_relevance(rows, "USP1", llm=llm, target_n=5)
    assert calls["n"] == 20                      # one batch at SCREEN_BATCH_MIN
    assert sum(1 for r in rows if r.get("screen_note")) == 20
    assert all("screen_note" not in r for r in rows[20:])


def test_screening_continues_when_the_first_batch_falls_short():
    """A symbol that is also an ordinary word discards most candidates; the
    screen must keep going rather than return a short brief."""
    rows = _rows(60)
    replies = ["NO"] * 18 + ["YES"] * 2 + ["YES"] * 40
    llm, calls = _counting_llm(replies)
    k.screen_relevance(rows, "SET", llm=llm, target_n=5)
    assert calls["n"] == 40                      # a second batch was needed
    assert sum(1 for r in rows if r.get("on_target") and r.get("screen_note")) >= 5


def test_screening_everything_when_no_target_is_given():
    rows = _rows(30)
    llm, calls = _counting_llm(["YES"] * 30)
    k.screen_relevance(rows, "USP1", llm=llm, target_n=None)
    assert calls["n"] == 30


def test_selection_drops_unscreened_and_off_target_rows():
    rows = _rows(30)
    k.screen_relevance(rows, "SET", llm=_stub_llm(["YES", "NO"] + ["YES"] * 28),
                       target_n=5)
    chosen = k.select_for_brief(rows, n_papers=5)
    assert len(chosen) == 5
    assert all(r["on_target"] and r["screen_note"] for r in chosen)
    assert "2" not in [r["pmid"] for r in chosen]      # the NO verdict


def test_selection_counts_after_screening_not_before():
    """Two off-target papers at the top must not reduce a 5-paper brief to 3."""
    rows = _rows(30)
    k.screen_relevance(rows, "SET", llm=_stub_llm(["NO", "NO"] + ["YES"] * 28),
                       target_n=5)
    assert len(k.select_for_brief(rows, n_papers=5)) == 5


def test_selection_ignores_an_on_target_flag_this_screen_did_not_set():
    """A pool merged across two paged fetches can carry an on_target set
    elsewhere; only a verdict from this screen counts as one."""
    rows = _rows(4)
    k.screen_relevance(rows, "USP1", llm=_stub_llm(["YES", "YES"]), target_n=2,
                       batch_min=2, batch_factor=1)
    assert sum(1 for r in rows if r.get("screen_note")) == 2
    rows[2]["on_target"] = True            # asserted by something else, unjudged here
    chosen = k.select_for_brief(rows, n_papers=5)
    assert [r["pmid"] for r in chosen] == ["4", "3"]


def test_selection_returns_what_exists_when_the_pool_is_short():
    rows = _rows(3)
    k.screen_relevance(rows, "USP1", llm=_stub_llm(["YES"] * 3), target_n=5)
    assert len(k.select_for_brief(rows, n_papers=5)) == 3


# --------------------------------------------------------------------------
# summarisation
# --------------------------------------------------------------------------

def test_summaries_and_categories_are_attached():
    rows = _rows(2)
    replies = [json.dumps({"category": "dependency", "summary": "Knockout kills."}),
               json.dumps({"category": "chemical_matter", "summary": "New inhibitor."})]
    k.summarize_papers(rows, "USP1", llm=_stub_llm(replies))
    assert [r["category"] for r in rows] == ["dependency", "chemical_matter"]
    assert rows[0]["summary"] == "Knockout kills."


def test_unknown_category_falls_back_to_other():
    rows = _rows(1)
    k.summarize_papers(rows, "USP1", llm=_stub_llm(
        [json.dumps({"category": "spectacular", "summary": "x"})]))
    assert rows[0]["category"] == "other"


def test_json_wrapped_in_prose_or_fences_is_recovered():
    rows = _rows(2)
    replies = ['```json\n{"category": "mechanism", "summary": "Binds PP2A."}\n```',
               'Here you go: {"category": "clinical", "summary": "Phase I."} Hope that helps!']
    k.summarize_papers(rows, "USP1", llm=_stub_llm(replies))
    assert [r["category"] for r in rows] == ["mechanism", "clinical"]
    assert rows[0]["summary"] == "Binds PP2A."


def test_malformed_reply_is_flagged_not_invented():
    rows = _rows(2)
    k.summarize_papers(rows, "USP1", llm=_stub_llm(["not json at all", None]))
    assert all(r["summary"] == "" for r in rows)
    assert all(r["summary_note"] == "unsummarized" for r in rows)
    assert all(r["category"] == "other" for r in rows)


# --------------------------------------------------------------------------
# outputs
# --------------------------------------------------------------------------

def test_table_excludes_abstract_and_leads_with_pmid(tmp_path):
    rows = _rows(3)
    k.screen_relevance(rows, "USP1", llm=_stub_llm(["YES"] * 3))
    k.summarize_papers(rows, "USP1", llm=_stub_llm(
        [json.dumps({"category": "mechanism", "summary": "s"})] * 3))
    path = tmp_path / "papers.csv"
    frame = k.write_paper_table(rows, path)
    assert "abstract" not in frame.columns
    assert "_date_key" not in frame.columns
    assert list(frame.columns)[:4] == ["pmid", "date", "journal", "title"]
    assert path.read_text(encoding="utf-8").splitlines()[0].startswith("pmid,date")


def test_table_survives_commas_and_quotes_in_titles(tmp_path):
    import pandas as pd
    articles = [_article("1", "2026", "01", "01",
                         title='A "quoted", comma-laden title')]
    rows = k.rank_recent(articles, n_papers=1)
    path = tmp_path / "papers.csv"
    k.write_paper_table(rows, path)
    assert pd.read_csv(path).loc[0, "title"] == 'A "quoted", comma-laden title'


def test_coverage_counts_only_screened_papers_as_off_target():
    rows = _rows(4)
    k.screen_relevance(rows, "SET", llm=_stub_llm(["YES", "NO", None, "NO"]))
    kept = [r for r in rows if r["on_target"]]
    stats = k.coverage_stats(rows, kept)
    assert stats["n_retrieved"] == 4
    assert stats["n_screened"] == 3        # the None was not judged
    assert stats["n_off_target"] == 2
    assert stats["n_in_brief"] == 2


def test_coverage_date_span_describes_kept_papers_only():
    rows = _rows(3)
    kept = [r for r in rows if r["pmid"] in {"1", "2"}]
    stats = k.coverage_stats(rows, kept)
    assert stats["date_newest"] == "2026-01-02"
    assert stats["date_oldest"] == "2026-01-01"


def test_coverage_handles_all_undated():
    rows = k.rank_recent([{"identifiers": {"pmid": "1"}, "publication_date": None}], n_papers=1)
    stats = k.coverage_stats(rows, rows)
    assert stats["date_newest"] is None


def test_synthesis_prompt_carries_pmids_and_bans_invention():
    rows = _rows(2)
    k.summarize_papers(rows, "USP1", llm=_stub_llm(
        [json.dumps({"category": "dependency", "summary": "Knockout kills."})] * 2))
    prompt = k.synthesis_prompt(rows, "USP1", gene_name="ubiquitin specific peptidase 1")
    assert "PMID 1" in prompt and "PMID 2" in prompt
    assert "Only cite PMIDs from the list above" in prompt
    assert "dependency (2 papers)" in prompt


def test_synthesis_prompt_omits_absent_categories():
    rows = _rows(1)
    k.summarize_papers(rows, "USP1", llm=_stub_llm(
        [json.dumps({"category": "mechanism", "summary": "s"})]))
    prompt = k.synthesis_prompt(rows, "USP1")
    assert "## mechanism" in prompt
    assert "## resistance" not in prompt


def test_synthesis_prompt_states_the_reader_and_the_length_limits():
    rows = _rows(1)
    k.summarize_papers(rows, "USP1", llm=_stub_llm(
        [json.dumps({"category": "mechanism", "summary": "s"})]))
    prompt = k.synthesis_prompt(rows, "USP1")
    assert "## Overview" in prompt
    assert "does not work on this target" in prompt      # the reader
    assert "One idea per sentence" in prompt


# --------------------------------------------------------------------------
# citation renumbering
# --------------------------------------------------------------------------

def test_inline_pmids_become_numbered_markers_in_order_of_appearance():
    rows = _rows(3)
    text = "Second paper (PMID 2). First paper (PMID 1). Again (PMID 2)."
    out, refs, unknown = k.renumber_citations(text, rows)
    assert out == "Second paper [1]. First paper [2]. Again [1]."
    assert [r["pmid"] for r in refs] == ["2", "1"]
    assert unknown == []


def test_multiple_pmids_in_one_parenthetical_collapse_to_one_bracket():
    rows = _rows(3)
    out, refs, _ = k.renumber_citations("Shown (PMID 1, PMID 3).", rows)
    assert out == "Shown [1, 2]."
    assert len(refs) == 2


def test_citation_variants_are_matched():
    rows = _rows(2)
    for variant in ["(PMID 1)", "(PMID: 1)", "(PMIDs 1)", "(pmid 1)"]:
        out, _, _ = k.renumber_citations(f"Claim {variant}.", rows)
        assert out == "Claim [1].", variant


def test_pmid_outside_the_retrieved_set_is_reported_not_renumbered():
    """The prompt forbids outside citations, so one appearing means the model
    invented it; silently dropping it would hide an unsupported claim."""
    rows = _rows(2)
    out, refs, unknown = k.renumber_citations("Real [x] (PMID 1). Fake (PMID 99999).", rows)
    assert "[1]" in out
    assert "(PMID 99999)" in out          # left inline, visibly unconverted
    assert unknown == ["99999"]
    assert len(refs) == 1


def test_mixed_known_and_unknown_in_one_parenthetical_keeps_the_known():
    rows = _rows(2)
    out, refs, unknown = k.renumber_citations("Both (PMID 1, PMID 88888).", rows)
    assert out == "Both [1]."
    assert unknown == ["88888"]


def test_pmids_separated_by_prose_inside_one_parenthetical_are_converted():
    """The WRN brief contains '(PMID 1, correction at PMID 2)', which matches
    no whole-group pattern; both ids must still be numbered."""
    rows = _rows(3)
    out, refs, _ = k.renumber_citations("Found (PMID 1, correction at PMID 2).", rows)
    assert out == "Found ([1], correction at [2])."
    assert [r["pmid"] for r in refs] == ["1", "2"]


def test_the_word_pmid_without_a_number_is_left_alone():
    """Briefs describe their own CSV columns as 'PMID, date, journal'; that is
    prose, not a citation."""
    text = "one row per paper: PMID, date, journal, title."
    out, refs, unknown = k.renumber_citations(text, _rows(2))
    assert out == text
    assert refs == [] and unknown == []


def test_text_without_citations_is_unchanged():
    out, refs, unknown = k.renumber_citations("No citations here.", _rows(2))
    assert out == "No citations here."
    assert refs == [] and unknown == []


def test_reference_list_links_each_pmid_to_pubmed():
    rows = _rows(2)
    _, refs, _ = k.renumber_citations("A (PMID 1).", rows)
    listing = k.format_references(refs)
    assert listing.startswith("1. ")
    assert "https://pubmed.ncbi.nlm.nih.gov/1/" in listing
    assert "PMID 1" in listing


def test_brief_moves_citations_to_a_reference_section(tmp_path):
    rows = _rows(3)
    k.screen_relevance(rows, "USP1", llm=_stub_llm(["YES"] * 3))
    stats = k.coverage_stats(rows, rows)
    text = k.write_brief(tmp_path / "b.md", "USP1",
                         "## Overview\n\nIt does things (PMID 2).", stats, rows=rows)
    assert "(PMID 2)" not in text.split("## References")[0]
    assert "[1]" in text
    assert "## References" in text


def test_brief_without_rows_leaves_citations_inline(tmp_path):
    rows = _rows(2)
    k.screen_relevance(rows, "USP1", llm=_stub_llm(["YES"] * 2))
    stats = k.coverage_stats(rows, rows)
    text = k.write_brief(tmp_path / "b.md", "USP1", "Claim (PMID 1).", stats)
    assert "(PMID 1)" in text
    assert "## References" not in text


def test_brief_flags_an_invented_citation_in_the_method_note(tmp_path):
    rows = _rows(2)
    k.screen_relevance(rows, "USP1", llm=_stub_llm(["YES"] * 2))
    stats = k.coverage_stats(rows, rows)
    text = k.write_brief(tmp_path / "b.md", "USP1", "Claim (PMID 77777).", stats,
                         rows=rows)
    assert "not in the retrieved set" in text
    assert "77777" in text


def test_brief_reports_provenance_and_exclusions(tmp_path):
    rows = _rows(3)
    k.screen_relevance(rows, "USP1", llm=_stub_llm(["YES", "NO", "YES"]))
    kept = [r for r in rows if r["on_target"]]
    stats = k.coverage_stats(rows, kept)
    path = tmp_path / "brief.md"
    text = k.write_brief(path, "USP1", "## What changed\n\nThings.", stats,
                         gene_name="ubiquitin specific peptidase 1",
                         query='"USP1"[Title/Abstract]')
    assert "PubMed (NCBI)" in text
    assert "1 of 3 screened" in text
    assert '`"USP1"[Title/Abstract]`' in text
    assert path.exists()


def test_brief_flags_unsummarized_papers(tmp_path):
    rows = _rows(2)
    k.screen_relevance(rows, "USP1", llm=_stub_llm(["YES", "YES"]))
    k.summarize_papers(rows, "USP1", llm=_stub_llm(["garbage", "garbage"]))
    stats = k.coverage_stats(rows, rows)
    text = k.write_brief(tmp_path / "b.md", "USP1", "x", stats)
    assert "2 paper(s) could not be summarised" in text
