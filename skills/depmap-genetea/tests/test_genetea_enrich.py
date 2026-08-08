"""Tests for the GeneTEA half of kernel.py.

The trained model is 1.11 GB and takes ~2.6 s to load, so it is a session-scoped
fixture and every model-backed test shares one copy. Tests that do not need the
model (path resolution, excerpt splitting) run without it, so a machine with no
model file still exercises those.

Run:  pytest -q        (environment: genetea)
Skip the model-backed tests by unsetting the file / pointing $GENETEA_MODEL
elsewhere; they are marked with the `needs_model` fixture and will skip.
"""

import os

import pandas as pd
import pytest

import kernel as ge

PYRIMIDINE = ["CAD", "DHODH", "UMPS"]


@pytest.fixture(scope="session")
def tea():
    path = os.environ.get("GENETEA_MODEL", ge.DEFAULT_MODEL_PATH)
    if not path or not os.path.isfile(path):
        pytest.skip(f"trained model not present at {path!r}")
    return ge.genetea_load(path)


def test_load_missing_path_names_the_path(tmp_path):
    missing = str(tmp_path / "absent.pkl")
    with pytest.raises(FileNotFoundError, match="absent.pkl"):
        ge.genetea_load(missing)


def test_load_is_memoised(tea):
    assert ge.genetea_load() is tea


def test_enrich_reproduces_pyrimidine_terms(tea):
    out = ge.genetea_enrich(PYRIMIDINE, n=10)
    assert isinstance(out, pd.DataFrame) and len(out) > 0
    # The documented columns lead; GeneTEA's extra annotation columns trail.
    assert list(out.columns[: len(ge.DISCRETE_COLUMNS)]) == ge.DISCRETE_COLUMNS
    assert (out["FDR"] < 0.05).all()
    # The de novo pyrimidine pathway must be named, via the orotate synonym group.
    assert "~ orotate" in set(out["Term"])
    assert "~ pyrimidine" in set(out["Term"])


def test_enrich_returns_empty_frame_not_none(tea):
    # A single gene cannot satisfy min_n_list=2, so nothing is enriched and
    # GeneTEA returns None upstream.
    out = ge.genetea_enrich(["TP53"], min_n_list=2)
    assert isinstance(out, pd.DataFrame)
    assert out.empty
    assert list(out.columns) == ge.DISCRETE_COLUMNS


def test_enrich_rejects_empty_gene_list(tea):
    with pytest.raises(ValueError, match="empty"):
        ge.genetea_enrich([])


def test_sort_by_changes_the_ordering(tea):
    es = ge.genetea_enrich(PYRIMIDINE, n=10, sort_by="Effect Size")
    sig = ge.genetea_enrich(PYRIMIDINE, n=10, sort_by="Significance")
    assert set(es["Term"]) == set(sig["Term"])
    assert es["Term"].tolist() != sig["Term"].tolist()
    assert es["Effect Size"].is_monotonic_decreasing
    assert sig["FDR"].is_monotonic_increasing


def test_group_subterms_caps_groups_not_rows(tea):
    big = ["CAD", "DHODH", "UMPS", "TYMS", "DHFR", "RRM1", "RRM2", "TK1",
           "CTPS1", "NME1", "CMPK1", "UCK2"]
    grouped = ge.genetea_enrich(big, n=3, group_subterms=True)
    flat = ge.genetea_enrich(big, n=3, group_subterms=False)
    assert len(flat) <= 3
    assert grouped["Term Group"].nunique() <= 3
    # Capping groups admits every member row, so it can exceed n.
    assert len(grouped) > len(flat)


def test_resolve_term_adds_synonym_prefix(tea):
    assert ge.genetea_resolve_term("orotate") == "~ orotate"
    assert ge.genetea_resolve_term("~ orotate") == "~ orotate"
    # A member of the group resolves to the group name.
    assert ge.genetea_resolve_term("orotidine") == "~ orotate"


def test_resolve_term_raises_with_suggestions(tea):
    with pytest.raises(KeyError, match="Closest entries"):
        ge.genetea_resolve_term("pyrimidinne")


def test_term_context_traces_to_named_sources(tea):
    ctx = ge.genetea_term_context("orotate", PYRIMIDINE)
    assert list(ctx.columns) == ["Gene", "Term", "Source", "Excerpt"]
    assert len(ctx) > 0
    assert set(ctx["Gene"]) <= set(PYRIMIDINE)
    known = {"Names/Aliases", "Gene Location", "RefSeq", "UniProt", "Alliance", "CIViC"}
    assert set(ctx["Source"]) <= known
    # The excerpt must actually contain the term that was queried.
    assert ctx["Excerpt"].str.lower().str.contains("orotate|orotidine").any()
    assert (ctx["Excerpt"].str.len() > 0).all()


def test_term_context_rejects_unknown_gene(tea):
    # All genes unknown: upstream raises RuntimeError, converted to ValueError.
    with pytest.raises(ValueError, match="GeneTEA corpus"):
        ge.genetea_term_context("orotate", ["NOT_A_GENE_XYZ"])
    # Some genes unknown: upstream returns them in `invalid`, also a ValueError,
    # and the unknown gene must be named rather than silently dropped.
    with pytest.raises(ValueError, match="NOT_A_GENE_XYZ"):
        ge.genetea_term_context("orotate", ["DHODH", "NOT_A_GENE_XYZ"])


def test_term_context_source_filter(tea):
    only = ge.genetea_term_context("orotate", PYRIMIDINE, source="UniProt")
    assert len(only) > 0
    assert set(only["Source"]) == {"UniProt"}
    assert len(only) < len(ge.genetea_term_context("orotate", PYRIMIDINE))
    # An unrecognised source is a caller error, not a no-match result.
    with pytest.raises(ValueError, match="source must be one of"):
        ge.genetea_term_context("orotate", PYRIMIDINE, source="Refseq")


def test_term_context_empty_when_no_gene_has_term(tea):
    out = ge.genetea_term_context("orotate", ["ACTB"])
    assert isinstance(out, pd.DataFrame) and out.empty
    assert list(out.columns) == ["Gene", "Term", "Source", "Excerpt"]


def test_split_excerpt_keeps_pipes_inside_text():
    # The discriminating case: a description sentence that itself contains both
    # " | " and a colon. Splitting on " | " and taking whatever precedes ": " as
    # the source would invent a source named "beta".
    sources = {"RefSeq", "UniProt"}
    raw = "RefSeq: alpha | beta: gamma | UniProt: delta"
    out = ge.genetea_split_excerpt(raw, sources)
    assert out == [("RefSeq", "alpha | beta: gamma"), ("UniProt", "delta")]
    assert [src for src, _ in out] == ["RefSeq", "UniProt"]


def test_resolve_term_prefixes_when_model_has_no_synonym_index(tea, monkeypatch):
    # Every synonym group in the shipped model has its base word in all_syns, so
    # the group-membership branch normally answers first. A model trained without
    # synonyms has all_syns None, and the "~ " prefix must still be added.
    monkeypatch.setattr(tea, "all_syns", None, raising=False)
    assert ge.genetea_resolve_term("orotate", tea=tea) == "~ orotate"


def test_continuous_rejects_short_vector_naming_the_constraint(tea):
    short = pd.Series(
        {g: float(i) for i, g in enumerate(list(tea.entities)[:500])}, dtype=float
    )
    with pytest.raises(ValueError, match="10,000"):
        ge.genetea_correlated_terms(short)
    with pytest.raises(ValueError, match="aligned"):
        ge.genetea_correlated_terms(short)


def test_continuous_rejects_non_float_series(tea):
    ints = pd.Series({g: 1 for g in list(tea.entities)[:20000]})
    with pytest.raises(TypeError, match="float"):
        ge.genetea_correlated_terms(ints)


def test_continuous_rejects_non_series(tea):
    with pytest.raises(TypeError, match="Series"):
        ge.genetea_correlated_terms([0.1] * 20000)


def test_enrich_continuous_requires_exactly_one_input(tea):
    with pytest.raises(ValueError, match="exactly one"):
        ge.genetea_enrich_continuous()
    with pytest.raises(ValueError, match="exactly one"):
        ge.genetea_enrich_continuous(x=pd.Series(dtype=float), corr_terms=pd.DataFrame())


def test_continuous_f_test_on_whole_genome(tea):
    # A whole-genome vector built from the corpus itself, so the test needs no
    # DepMap file: a deterministic ramp over every corpus gene.
    genes = list(tea.entities)
    x = pd.Series({g: float(i % 997) / 997.0 for i, g in enumerate(genes)}, dtype=float)
    corr = ge.genetea_correlated_terms(x, max_fdr=None, corr_thresh=None)
    assert len(corr) == len(tea.terms)
    assert corr["Correlation"].abs().max() <= 1.0
    assert {"Term", "Correlation", "p-val", "FDR", "Effect Size"} <= set(corr.columns)


def test_default_corr_thresh_is_below_genetea_default():
    # GeneTEA's own 0.2 default is unreachable for term-level correlations.
    assert ge.DEFAULT_CORR_THRESH < 0.2
