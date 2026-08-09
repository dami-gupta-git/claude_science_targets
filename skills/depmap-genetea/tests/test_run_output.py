"""Coverage for the results/<topic>/<run>/ writer added to depmap-genetea.

Uses a tmp_path root and synthetic frames shaped like depmap_codependencies()'s
and genetea_enrich()'s real return columns, so no trained model or DepMap
release is required.
"""
import os

import pandas as pd
import pytest

import kernel as ge


def test_results_root_raises_when_unconfigured(monkeypatch):
    monkeypatch.delenv("SCIENCE_RESULTS_ROOT", raising=False)
    with pytest.raises(FileNotFoundError, match="SCIENCE_RESULTS_ROOT"):
        ge.genetea_results_root()


def test_run_dir_raises_and_creates_nothing_when_root_unconfigured(monkeypatch, tmp_path):
    monkeypatch.delenv("SCIENCE_RESULTS_ROOT", raising=False)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError):
        ge.genetea_run_dir("TP53")
    assert not (tmp_path / "results").exists()


def test_run_dir_honours_env_var(monkeypatch, tmp_path):
    monkeypatch.setenv("SCIENCE_RESULTS_ROOT", str(tmp_path))
    out_dir = ge.genetea_run_dir("TP53")
    assert out_dir == str(tmp_path / "depmap_genetea" / "tp53")


def codependencies():
    return pd.DataFrame([
        {"gene": "MDM2", "r": -0.718, "n_lines": 1180, "z": -6.2, "rank": 1,
         "p": 1e-30, "fdr": 1e-28},
        {"gene": "MDM4", "r": -0.522, "n_lines": 1150, "z": -4.1, "rank": 2,
         "p": 1e-15, "fdr": 1e-13},
    ])


def terms():
    return pd.DataFrame([
        {"Term": "~ p53", "Synonyms": "p53|TP53", "p-val": 1e-20, "FDR": 1e-18,
         "Effect Size": 4.2, "n Matching Genes in List": 12,
         "n Matching Genes Overall": 40, "Matching Genes in List": "MDM2 MDM4",
         "Term Group": "~ p53"},
    ])


def provenance():
    return pd.DataFrame([
        {"Gene": "MDM2", "Term": "~ p53", "Source": "UniProt",
         "Excerpt": "Binds and inhibits p53 transactivation."},
    ])


def test_run_dir_slugs_name_and_makes_scripts(tmp_path):
    out_dir = ge.genetea_run_dir("TP53", results_root=str(tmp_path))
    assert out_dir == str(tmp_path / "depmap_genetea" / "tp53")
    assert os.path.isdir(os.path.join(out_dir, "scripts"))


def test_run_dir_rejects_empty_name(tmp_path):
    with pytest.raises(ValueError):
        ge.genetea_run_dir("   ", results_root=str(tmp_path))


def test_write_run_needs_at_least_one_table(tmp_path):
    out_dir = ge.genetea_run_dir("EMPTY", results_root=str(tmp_path))
    with pytest.raises(ValueError):
        ge.genetea_write_run(out_dir, "EMPTY", summary="nothing to report")


def test_write_run_codependencies_only(tmp_path):
    out_dir = ge.genetea_run_dir("TP53", results_root=str(tmp_path))
    written = ge.genetea_write_run(
        out_dir, "TP53", gene="TP53",
        summary="TP53's strongest codependencies are MDM2 and MDM4, its "
                "known negative regulators.",
        codependencies=codependencies(),
        data_sources=["DepMap 24Q2 CRISPRGeneEffect"])
    assert os.path.isfile(written["codependencies"])
    assert "terms" not in written
    text = open(written["readme"]).read()
    assert "MDM2" in text and "### Codependencies" in text
    assert "## Limits" in text and "## Files" in text and "## Data sources" in text


def test_write_run_terms_and_provenance(tmp_path):
    out_dir = ge.genetea_run_dir("TP53 partners", results_root=str(tmp_path))
    written = ge.genetea_write_run(
        out_dir, "TP53 partners", gene="TP53",
        summary="The MDM2/MDM4 codependency partners are named by GeneTEA as "
                "p53-pathway terms.",
        terms=terms(), provenance=provenance(),
        data_sources=["GeneTEA (Boyle et al. 2025)"])
    assert os.path.isfile(written["terms"])
    assert os.path.isfile(written["provenance"])
    text = open(written["readme"]).read()
    assert "GeneTEA enriched terms" in text
    assert "Provenance checked for 1 term(s)" in text


def test_write_run_provenance_as_list_of_dicts(tmp_path):
    """provenance passed as a plain list (not a DataFrame) must not crash."""
    out_dir = ge.genetea_run_dir("TP53 list provenance", results_root=str(tmp_path))
    written = ge.genetea_write_run(
        out_dir, "TP53 list provenance", gene="TP53",
        summary="Terms named for the MDM2/MDM4 codependency partners.",
        terms=terms(), provenance=provenance().to_dict("records"),
        data_sources=["GeneTEA (Boyle et al. 2025)"])
    text = open(written["readme"]).read()
    assert "Provenance checked for 1 term(s)" in text


def test_write_run_continuous_term_kind_uses_correlation_columns(tmp_path):
    out_dir = ge.genetea_run_dir("TP53 continuous", results_root=str(tmp_path))
    cont = pd.DataFrame([
        {"Term": "~ MDM2", "Synonyms": "", "Correlation": -0.62, "p-val": 1e-10,
         "FDR": 1e-8, "Effect Size": 3.1, "Directed Effect Size": -3.1,
         "n Matching Genes Overall": 5, "Term Group": "~ MDM2"},
    ])
    written = ge.genetea_write_run(
        out_dir, "TP53 continuous", term_kind="continuous", terms=cont,
        summary="Continuous-mode terms correlated against the TP53 codependency vector.")
    text = open(written["readme"]).read()
    assert "correlated terms" in text


def test_run_readme_rejects_bulleted_summary():
    with pytest.raises(ValueError):
        ge.genetea_run_readme("TP53", "- not prose", codependencies=codependencies())


def test_run_readme_enforces_summary_word_cap():
    with pytest.raises(ValueError):
        ge.genetea_run_readme("TP53", " ".join(["word"] * 200),
                              codependencies=codependencies())
