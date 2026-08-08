"""Branch coverage for fuse_target_row's verdict logic.

fuse_target_row is pure (no filesystem/network), so every branch and the
precedence between them can be pinned down exactly with synthetic inputs
mirroring depmap_selectivity()'s output shape.
"""
from kernel import VERDICTS, fuse_target_row, fusion_notes


def stats(mean_effect=None, frac_dependent=0.0, classification=None,
          frac_positive=None, sd_effect=None):
    return {"mean_effect": mean_effect, "frac_dependent": frac_dependent,
            "classification": classification, "frac_positive": frac_positive,
            "sd_effect": sd_effect}


def test_indeterminate_when_dispersion_missing():
    row = fuse_target_row("X", ot_score=0.9,
                          effect_stats=stats(mean_effect=-0.8, frac_dependent=0.9,
                                              classification="common essential",
                                              frac_positive=None, sd_effect=None))
    assert row["verdict"] == "indeterminate"
    assert row["knockout_actionable"] is False


def test_common_essential_beats_concordant_evidence():
    # Strong evidence + high frac_dependent would otherwise look concordant,
    # but classification == "common essential" must take precedence: no
    # therapeutic window regardless of how well evidence and dependency agree.
    row = fuse_target_row("KRAS", ot_score=0.9,
                          effect_stats=stats(mean_effect=-0.8, frac_dependent=0.95,
                                              classification="common essential",
                                              frac_positive=0.01, sd_effect=0.15))
    assert row["verdict"] == "common-essential-no-window"
    assert row["knockout_actionable"] is False


def test_growth_suppressive_mismatch_beats_strong_evidence():
    # TP53-like: strong OT evidence AND frac_dependent high enough to pass the
    # concordant threshold, but the positive-tail/high-variance suppressor
    # signature must still win.
    row = fuse_target_row("TP53", ot_score=0.9,
                          effect_stats=stats(mean_effect=0.5, frac_dependent=0.5,
                                              classification="selective",
                                              frac_positive=0.3, sd_effect=0.3))
    assert row["verdict"] == "growth-suppressive-mismatch"
    assert row["knockout_actionable"] is False


def test_context_restricted_requires_lineage_more_dependent_than_overall():
    row = fuse_target_row("X", ot_score=None,
                          effect_stats=stats(mean_effect=-0.1, frac_dependent=0.01,
                                              classification="selective",
                                              frac_positive=0.02, sd_effect=0.15),
                          lineage_row={"mean_effect": -0.7, "frac_dependent": 0.6,
                                       "p_bh": 0.001})
    assert row["verdict"] == "context-restricted"
    assert row["knockout_actionable"] is True


def test_context_restricted_not_triggered_when_lineage_not_more_dependent():
    # lineage passes the dependency/significance bar but isn't meaningfully
    # more dependent than the overall mean (mean_eff - 0.2 margin) -> should
    # fall through to a different verdict, not context-restricted.
    row = fuse_target_row("X", ot_score=None,
                          effect_stats=stats(mean_effect=-0.55, frac_dependent=0.5,
                                              classification="selective",
                                              frac_positive=0.02, sd_effect=0.3),
                          lineage_row={"mean_effect": -0.55, "frac_dependent": 0.5,
                                       "p_bh": 0.001})
    assert row["verdict"] != "context-restricted"


def test_concordant_dependency():
    row = fuse_target_row("X", ot_score=0.7,
                          effect_stats=stats(mean_effect=-0.7, frac_dependent=0.5,
                                              classification="selective",
                                              frac_positive=0.01, sd_effect=0.3))
    assert row["verdict"] == "concordant-dependency"
    assert row["knockout_actionable"] is True


def test_inert_in_panel():
    row = fuse_target_row("NRG1", ot_score=None,
                          effect_stats=stats(mean_effect=0.05, frac_dependent=0.01,
                                              classification="non-essential",
                                              frac_positive=0.02, sd_effect=0.1))
    assert row["verdict"] == "inert-in-panel"
    assert row["knockout_actionable"] is False


def test_evidence_without_dependency_strong_evidence_not_inert():
    # High variance disqualifies "inert" even though frac_dependent is low,
    # so strong evidence with no dependency signal must land here instead.
    row = fuse_target_row("X", ot_score=0.8,
                          effect_stats=stats(mean_effect=0.05, frac_dependent=0.02,
                                              classification="non-essential",
                                              frac_positive=0.05, sd_effect=0.3))
    assert row["verdict"] == "evidence-without-dependency"


def test_dependency_without_evidence():
    row = fuse_target_row("X", ot_score=None,
                          effect_stats=stats(mean_effect=-0.6, frac_dependent=0.15,
                                              classification="selective",
                                              frac_positive=0.02, sd_effect=0.3))
    assert row["verdict"] == "dependency-without-evidence"


def test_evidence_without_dependency_weak_evidence():
    # ot_score above EVIDENCE_FLOOR but below the strong-evidence threshold.
    row = fuse_target_row("X", ot_score=0.2,
                          effect_stats=stats(mean_effect=0.05, frac_dependent=0.02,
                                              classification="non-essential",
                                              frac_positive=0.05, sd_effect=0.3))
    assert row["verdict"] == "evidence-without-dependency"


def test_no_evidence_no_dependency():
    row = fuse_target_row("X", ot_score=0.05,
                          effect_stats=stats(mean_effect=0.05, frac_dependent=0.02,
                                              classification="non-essential",
                                              frac_positive=0.05, sd_effect=0.3))
    assert row["verdict"] == "no-evidence-no-dependency"


def test_ot_score_rounded_and_none_preserved():
    row = fuse_target_row("X", ot_score=0.733333333,
                          effect_stats=stats(mean_effect=0.05, frac_dependent=0.02,
                                              classification="non-essential",
                                              frac_positive=0.05, sd_effect=0.3))
    assert row["ot_score"] == 0.733

    row_none = fuse_target_row("X", ot_score=None,
                               effect_stats=stats(mean_effect=0.05, frac_dependent=0.02,
                                                   classification="non-essential",
                                                   frac_positive=0.05, sd_effect=0.3))
    assert row_none["ot_score"] is None


def test_depmap_inferred_essential_passthrough():
    es = stats(mean_effect=0.05, frac_dependent=0.02, classification="non-essential",
              frac_positive=0.05, sd_effect=0.3)
    es["depmap_inferred_essential"] = True
    row = fuse_target_row("X", ot_score=None, effect_stats=es)
    assert row["depmap_inferred_essential"] is True


def test_thin_evidence_plus_dependency_is_its_own_verdict():
    # ot_score above EVIDENCE_FLOOR but below the strong bar, WITH a dependency
    # signal: both halves are real, so the row must not be labelled as missing
    # either one.
    row = fuse_target_row("X", ot_score=0.3,
                          effect_stats=stats(mean_effect=-0.6, frac_dependent=0.15,
                                              classification="selective",
                                              frac_positive=0.02, sd_effect=0.3))
    assert row["verdict"] == "dependency-with-thin-evidence"
    assert row["knockout_actionable"] is True


def test_thin_evidence_dependency_at_both_boundaries():
    # ot_score exactly at EVIDENCE_FLOOR and frac exactly at the dependency
    # cut: the lowest inputs that still satisfy both halves.
    row = fuse_target_row("X", ot_score=0.10,
                          effect_stats=stats(mean_effect=-0.6, frac_dependent=0.10,
                                              classification="selective",
                                              frac_positive=0.02, sd_effect=0.3))
    assert row["verdict"] == "dependency-with-thin-evidence"

    # Just below the strong-evidence bar stays in the thin tier; at the bar it
    # becomes concordant.
    assert fuse_target_row("X", ot_score=0.49,
                           effect_stats=stats(mean_effect=-0.6, frac_dependent=0.15,
                                               classification="selective",
                                               frac_positive=0.02, sd_effect=0.3)
                           )["verdict"] == "dependency-with-thin-evidence"
    assert fuse_target_row("X", ot_score=0.50,
                           effect_stats=stats(mean_effect=-0.6, frac_dependent=0.15,
                                               classification="selective",
                                               frac_positive=0.02, sd_effect=0.3)
                           )["verdict"] == "concordant-dependency"


def test_dependency_without_evidence_requires_absent_evidence():
    # Below EVIDENCE_FLOOR the same dependency signal keeps the
    # no-evidence label; None behaves identically.
    for score in (0.05, None):
        row = fuse_target_row("X", ot_score=score,
                              effect_stats=stats(mean_effect=-0.6, frac_dependent=0.15,
                                                  classification="selective",
                                                  frac_positive=0.02, sd_effect=0.3))
        assert row["verdict"] == "dependency-without-evidence"
        assert row["knockout_actionable"] is False


def test_thin_evidence_without_dependency_stays_evidence_without_dependency():
    row = fuse_target_row("X", ot_score=0.3,
                          effect_stats=stats(mean_effect=-0.6, frac_dependent=0.05,
                                              classification="selective",
                                              frac_positive=0.02, sd_effect=0.3))
    assert row["verdict"] == "evidence-without-dependency"


def test_higher_precedence_verdicts_still_win_over_thin_tier():
    # Each of common-essential / suppressive / context must outrank the thin
    # tier even when its own conditions (thin evidence + dependency) hold.
    common = fuse_target_row("X", ot_score=0.3,
                             effect_stats=stats(mean_effect=-0.8, frac_dependent=0.95,
                                                 classification="common essential",
                                                 frac_positive=0.01, sd_effect=0.15))
    assert common["verdict"] == "common-essential-no-window"

    suppressive = fuse_target_row("X", ot_score=0.3,
                                  effect_stats=stats(mean_effect=0.5, frac_dependent=0.5,
                                                      classification="selective",
                                                      frac_positive=0.3, sd_effect=0.3))
    assert suppressive["verdict"] == "growth-suppressive-mismatch"

    context = fuse_target_row("X", ot_score=0.3,
                              effect_stats=stats(mean_effect=-0.3, frac_dependent=0.15,
                                                  classification="selective",
                                                  frac_positive=0.02, sd_effect=0.3),
                              lineage_row={"mean_effect": -0.7, "frac_dependent": 0.6,
                                           "p_bh": 0.001})
    assert context["verdict"] == "context-restricted"


def test_verdict_vocabulary_and_notes_agree():
    assert "dependency-with-thin-evidence" in VERDICTS
    assert set(fusion_notes()) == set(VERDICTS)


def test_explicit_none_frac_dependent_does_not_raise():
    # depmap_selectivity() can return frac_dependent=None (key present, value
    # unset) rather than omitting the key -- dict.get(key, 0.0) does not catch
    # that case, so the fallback must coalesce None explicitly or every
    # frac >= / < comparison below raises TypeError.
    row = fuse_target_row("X", ot_score=0.8,
                          effect_stats=stats(mean_effect=0.05, frac_dependent=None,
                                              classification="non-essential",
                                              frac_positive=0.05, sd_effect=0.3))
    assert row["frac_dependent"] == 0.0
    assert row["verdict"] == "evidence-without-dependency"
