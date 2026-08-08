"""Cross-KB fusion: Open Targets human evidence x DepMap cell-line dependency."""

VERDICTS = (
    "concordant-dependency",
    "common-essential-no-window",
    "growth-suppressive-mismatch",
    "context-restricted",
    "inert-in-panel",
    "evidence-without-dependency",
    "dependency-with-thin-evidence",
    "dependency-without-evidence",
    "no-evidence-no-dependency",
    "indeterminate",
)

# Minimum OT association score for a target to count as having ANY human
# evidence. Below this (or ot_score=None) the `evidence-without-dependency`
# label would assert evidence that was never supplied - such rows get
# `no-evidence-no-dependency` instead.
EVIDENCE_FLOOR = 0.10

# A real growth suppressor shows a POSITIVE tail across many lines; an inert
# gene sits flat near zero. Calibrated on TP53/RB1 (suppressors) vs
# ROS1/NRG1/ALK (inert in the DepMap panel) - see SKILL.md ## Calibration.
SUPPRESSOR_MIN_FRAC_POS = 0.15
SUPPRESSOR_MIN_SD = 0.20
INERT_MAX_SD = 0.20
INERT_MAX_ABS_MEAN = 0.25

# DepMap ships CRISPRInferredCommonEssentials.csv (1827 genes). It is NOT a
# frequency threshold: it includes KRAS (52% of lines dependent) and some genes
# at 0%. Our frac_dependent >= 0.90 rule yields 974 genes - a strict SUBSET.
# We keep the stricter rule for the no-window verdict, because adopting the
# broader list would discard KRAS, and surface DepMap's call alongside it as
# `depmap_inferred_essential` for cross-reference.


def fuse_target_row(symbol, ot_score, effect_stats, lineage_row=None,
                    tractable_modalities=None, dep_threshold=-0.5):
    """Classify one target by how Open Targets evidence and DepMap dependency agree.

    effect_stats: dict from depmap_selectivity().
    lineage_row: optional dict with mean_effect / frac_dependent / p_bh for the
        disease-relevant lineage.
    Returns a dict with a `verdict` plus the numbers behind it.
    """
    mean_eff = effect_stats.get("mean_effect")
    frac = effect_stats.get("frac_dependent")
    frac = 0.0 if frac is None else frac
    cls = effect_stats.get("classification")
    lin_mean = (lineage_row or {}).get("mean_effect")
    lin_frac = (lineage_row or {}).get("frac_dependent")
    lin_p = (lineage_row or {}).get("p_bh")

    strong_evidence = ot_score is not None and ot_score >= 0.5
    has_evidence = ot_score is not None and ot_score >= EVIDENCE_FLOOR
    context = (lin_mean is not None and lin_p is not None
               and lin_mean <= dep_threshold and lin_p < 0.05
               and (mean_eff is None or lin_mean < mean_eff - 0.2))

    frac_pos = effect_stats.get("frac_positive")
    sd = effect_stats.get("sd_effect")

    suppressive = (mean_eff is not None and mean_eff > 0.1
                   and frac_pos is not None and sd is not None
                   and frac_pos >= SUPPRESSOR_MIN_FRAC_POS and sd >= SUPPRESSOR_MIN_SD)
    inert = (sd is not None and mean_eff is not None
             and sd < INERT_MAX_SD and abs(mean_eff) < INERT_MAX_ABS_MEAN
             and frac < 0.05)

    if sd is None or frac_pos is None or mean_eff is None:
        # The growth-suppressive and inert tests cannot run without dispersion.
        # This branch is FIRST deliberately: reporting any other verdict here
        # would assert a distinction that was never tested.
        verdict = "indeterminate"
    elif cls == "common essential":
        # Every cell needs it - killing tumour cells means killing normal tissue.
        verdict = "common-essential-no-window"
    elif suppressive:
        verdict = "growth-suppressive-mismatch"
    elif context:
        verdict = "context-restricted"
    elif strong_evidence and frac >= 0.10:
        verdict = "concordant-dependency"
    elif inert:
        verdict = "inert-in-panel"
    elif strong_evidence:
        verdict = "evidence-without-dependency"
    elif has_evidence and frac >= 0.10:
        # Evidence above EVIDENCE_FLOOR but below the strong bar, WITH a
        # dependency signal. Both halves are real, so neither
        # `dependency-without-evidence` nor `evidence-without-dependency`
        # describes the row: this branch must precede both.
        verdict = "dependency-with-thin-evidence"
    elif frac >= 0.10:
        verdict = "dependency-without-evidence"
    elif has_evidence:
        verdict = "evidence-without-dependency"
    else:
        verdict = "no-evidence-no-dependency"

    return {
        "gene": symbol,
        "depmap_inferred_essential": effect_stats.get("depmap_inferred_essential"),
        "ot_score": None if ot_score is None else round(float(ot_score), 3),
        "depmap_class": cls,
        "mean_effect": mean_eff,
        "frac_dependent": frac,
        "lineage_mean_effect": lin_mean,
        "lineage_frac_dependent": lin_frac,
        "lineage_p_bh": lin_p,
        "tractable_modalities": tractable_modalities or [],
        "verdict": verdict,
        "knockout_actionable": verdict in ("concordant-dependency",
                                          "context-restricted",
                                          "dependency-with-thin-evidence"),
    }


def fusion_notes():
    """Interpretation guide for the verdict vocabulary."""
    return {
        "concordant-dependency":
            "Human evidence and cell-line dependency agree. Strongest KO/degrader case.",
        "common-essential-no-window":
            "Nearly every line depends on it (frac_dependent >= 0.90). Real "
            "dependency, but no therapeutic window - normal tissue needs it too. "
            "Only viable with tumour-selective delivery or a genotype-restricted "
            "partner.",
        "growth-suppressive-mismatch":
            "High OT evidence but knockout FAVOURS growth across many lines "
            "(positive tail + high variance). Not a knockout target; consider "
            "synthetic-lethal partners instead.",
        "inert-in-panel":
            "Flat, low-variance effect near zero: the panel lacks the context that "
            "makes this gene matter (e.g. fusion-driven or ligand-driven oncogenes "
            "with no representative line). Absence of evidence, not evidence of "
            "absence - check lineage coverage before dismissing.",
        "context-restricted":
            "Dependency confined to a lineage or genotype. Target with a biomarker "
            "hypothesis, not broadly.",
        "evidence-without-dependency":
            "Genetic/clinical evidence without cell-autonomous dependency. May act "
            "non-cell-autonomously, be redundant, or need a model the panel lacks.",
        "dependency-with-thin-evidence":
            "Evidence above EVIDENCE_FLOOR but below the strong bar (ot_score "
            "in [0.10, 0.50)) together with a dependency signal "
            "(frac_dependent >= 0.10). Both halves are real, so the row is "
            "actionable, but on weaker evidence than `concordant-dependency` - "
            "read the ot_score alongside the verdict.",
        "dependency-without-evidence":
            "Cells need it but human evidence is thin - possible novel target or a "
            "general fitness gene with a poor therapeutic window.",
        "no-evidence-no-dependency":
            "Neither side supports this target: OT score below EVIDENCE_FLOOR (or "
            "absent) AND no cell-autonomous dependency. Nothing to pursue, and "
            "distinct from `evidence-without-dependency`, which requires real "
            "human evidence.",
        "indeterminate":
            "Required dispersion inputs (sd_effect / frac_positive / mean_effect) "
            "were missing, so the growth-suppressive and inert tests could not "
            "run. Pass the full depmap_selectivity() dict.",
    }
