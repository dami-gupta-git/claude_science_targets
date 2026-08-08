"""Tests for codep_enrich / codep_enrich_background -- the join between halves.

These are the only functions that need BOTH the DepMap cache and the 1.1 GB
model, so every test here skips when either is absent. The background-size and
gate assertions are the ones that matter: they encode the two findings the
skill's thresholds rest on.
"""
import os

import pandas as pd
import pytest

import kernel as K

_depmap_root = os.environ.get("DEPMAP_ROOT", K.DEPMAP_ROOT)
_genetea_model = os.environ.get("GENETEA_MODEL", K.DEFAULT_MODEL_PATH)
HAVE_DEPMAP = bool(_depmap_root) and os.path.isdir(
    os.path.join(_depmap_root, K.CACHE_DIRNAME))
HAVE_MODEL = bool(_genetea_model) and os.path.isfile(_genetea_model)
needs_both = pytest.mark.skipif(
    not (HAVE_DEPMAP and HAVE_MODEL),
    reason="needs both the DepMap parquet cache and the trained GeneTEA model")


@needs_both
def test_background_is_the_intersection_not_either_side():
    """The background must be screened AND in-corpus, never one alone.

    28 screened symbols are absent from the corpus; including them inflates the
    denominator with genes that cannot carry a term.
    """
    bg = K.codep_enrich_background()
    tea = K.genetea_load()
    screened = set(K.depmap_screened_genes(min_lines=K.CODEP_MIN_LINES))
    corpus = set(pd.Index(tea.entities).astype(str))
    assert set(bg) == screened & corpus
    # Strictly smaller than either input, i.e. the intersection really bites.
    assert len(bg) < len(screened)
    assert len(bg) < len(corpus)
    # Sorted, deduped, and far below the full corpus (the artefact we avoid).
    assert bg == sorted(set(bg))
    assert len(bg) < 0.6 * len(corpus)


@needs_both
def test_codep_enrich_recovers_the_complex_it_queries():
    """SMARCA4 partners must be named as SWI/SNF, not generic chromatin words."""
    partners, terms = K.codep_enrich("SMARCA4", n_partners=50, n=10)
    assert len(partners) == 50
    assert not terms.empty
    top = set(partners["gene"].head(6))
    assert {"SMARCB1", "SMARCC1"} <= top, top
    vocab = " ".join(terms["Term"].astype(str)).upper()
    assert "SNF" in vocab or "SWI" in vocab, vocab


@needs_both
def test_gate_refuses_the_negative_control_and_names_both_criteria():
    """OR5A1 fails on profile_sd AND min_effect; the error must say so."""
    with pytest.raises(ValueError) as exc:
        K.codep_enrich("OR5A1", n_partners=50)
    msg = str(exc.value)
    assert "OR5A1" in msg
    assert "profile_sd" in msg and "min_effect" in msg
    assert "check_query=False" in msg, "error must name the override"


@needs_both
def test_gate_is_overridable():
    partners, terms = K.codep_enrich(
        "OR5A1", n_partners=50, n=10, check_query=False)
    assert len(partners) == 50
    assert isinstance(terms, pd.DataFrame)   # empty is a legitimate answer


def test_dedupe_matches_upstream_separator_exactly():
    """A synonym group repeats the same sentence, joined by "..." after a ".".

    The doubled full stop is the whole difficulty: splitting "UMP)....FUNCTION"
    on "..." leaves a leading "." on the second copy, so naive equality treats
    the copies as distinct. Pure-string test, no model needed.
    """
    one = "FUNCTION: Converts orotate to OMP."
    assert K.genetea_dedupe_sentences(one + "..." + one) == one
    # When the source itself has no final period there is none to restore, so
    # the collapsed form is the unpunctuated sentence.
    bare = one.rstrip(".")
    assert K.genetea_dedupe_sentences(bare + "..." + bare) == bare
    # Genuinely different, non-duplicate sentences are both preserved, in
    # order, WITH the trailing period the first one actually has in the
    # source -- the split's boundary artifact must be reattached to the
    # fragment that owns it, not merely stripped off the next fragment.
    # (A naive fix that just discards the stray leading "." silently drops
    # the first sentence's period instead of restoring it.)
    two = "FUNCTION: Decarboxylates OMP to UMP."
    result = K.genetea_dedupe_sentences(one + "..." + two)
    assert result == one + "..." + two, result
    assert result.startswith(one), "first sentence must keep its own period"


def test_dedupe_does_not_drop_a_non_duplicate_sentences_period():
    """Regression: a period landing on the split boundary must be restored to
    the sentence before it, not merely stripped off the fragment after it.

    Two distinct period-ending sentences (no repeat at all) previously came
    back as "...three dots..." with the first sentence's period silently
    gone -- the old code only restored a lost period when merging an actual
    duplicate, never for unique content that happened to sit before more text.
    """
    a = "FUNCTION: Converts orotate to OMP."
    b = "DISEASE: Deficiency causes orotic aciduria."
    assert K.genetea_dedupe_sentences(a + "..." + b) == a + "..." + b

    # Same check with a real duplicate of `a` after `b`, so the boundary fix
    # and the dedup logic have to cooperate rather than just one running.
    assert K.genetea_dedupe_sentences(a + "..." + b + "..." + a) == a + "..." + b


def test_continuous_floor_matches_upstream_at_the_boundary():
    """Upstream raises on < 10_000, so exactly 10_000 must be ACCEPTED.

    The wrapper previously used <=, rejecting a vector the underlying function
    would have run. Asserted on the constant and the comparison rather than by
    building a 10,000-gene vector, which would need the model.
    """
    import inspect
    src = inspect.getsource(K.genetea_correlated_terms)
    line = [l.strip() for l in src.splitlines() if "MIN_CONTINUOUS_GENES" in l
            and "if " in l][0]
    assert "< MIN_CONTINUOUS_GENES" in line and "<=" not in line, line
    assert K.MIN_CONTINUOUS_GENES == 10_000


@needs_both
def test_term_context_excerpts_are_deduplicated():
    ctx = K.genetea_term_context("orotate", ["CAD", "DHODH", "UMPS"])
    assert not ctx.empty
    for excerpt in ctx["Excerpt"]:
        parts = [p.strip().lower().strip(".;, ") for p in str(excerpt).split("...")]
        assert len(parts) == len(set(parts)), excerpt[:120]
        assert not str(excerpt).startswith("."), excerpt[:60]


@needs_both
def test_corpus_background_returns_more_terms_than_screened():
    """The whole-corpus background is the artefact, so it must inflate.

    Pinning the direction rather than a count: the inflation is the documented
    finding (427 vs 276 terms over six queries) and a change of sign here would
    mean the background argument stopped being passed through.
    """
    _, screened_terms = K.codep_enrich("VPS4A", n_partners=50, n=None)
    _, corpus_terms = K.codep_enrich(
        "VPS4A", n_partners=50, n=None, corpus_background=True)
    assert len(corpus_terms) > len(screened_terms)
