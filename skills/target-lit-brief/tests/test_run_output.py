"""Coverage for the results/<topic>/<run>/ writer added to target-lit-brief.

Uses a tmp_path root and a synthetic paper row rather than a live PubMed
fetch or LLM call.
"""
import os

import pytest

import kernel as k


def paper_row(pmid="12345678"):
    return {"pmid": pmid, "date": "2025-01-01", "journal": "J. Test",
            "title": "A test paper", "first_author": "Smith", "n_authors": 1,
            "category": "mechanism", "summary": f"Did a thing (PMID {pmid}).",
            "doi": "", "pmcid": "", "article_types": "", "on_target": True,
            "screen_note": "screened", "summary_note": "", "abstract": "..."}


def test_run_dir_slugs_symbol_and_makes_scripts(tmp_path):
    out_dir = k.brief_run_dir("USP1", root=str(tmp_path))
    assert out_dir == str(tmp_path / "target_lit_brief" / "usp1")
    assert os.path.isdir(os.path.join(out_dir, "scripts"))


def test_run_dir_rejects_empty_symbol(tmp_path):
    with pytest.raises(ValueError):
        k.brief_run_dir("   ", root=str(tmp_path))


def test_write_run_uses_symbol_prefixed_table_name(tmp_path):
    rows = [paper_row()]
    stats = k.coverage_stats(rows, rows)
    out_dir = k.brief_run_dir("USP1", root=str(tmp_path))
    written = k.brief_write_run(out_dir, "USP1", rows,
                                "USP1 does a thing (PMID 12345678).", stats,
                                gene_name="ubiquitin specific peptidase 1",
                                query="USP1[tiab]")
    assert os.path.basename(written["table"]) == "usp1_recent_papers.csv"
    assert os.path.isfile(written["table"])
    assert os.path.isfile(written["readme"])
    text = open(written["readme"]).read()
    assert "[1]" in text  # citation renumbered via rows=
    assert "## References" in text


def test_write_run_copies_scripts(tmp_path):
    script = tmp_path / "build.py"
    script.write_text("# wiring\n")
    rows = [paper_row()]
    stats = k.coverage_stats(rows, rows)
    out_dir = k.brief_run_dir("USP1", root=str(tmp_path / "results"))
    written = k.brief_write_run(out_dir, "USP1", rows, "Summary.", stats,
                                scripts=[str(script)])
    assert os.path.isfile(written["scripts"][0])
    assert written["scripts"][0].endswith("scripts/build.py")
