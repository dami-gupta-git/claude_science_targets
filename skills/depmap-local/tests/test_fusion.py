"""depmap_fusion_models / depmap_fusion_contrast edge cases.

canonical_only crashed on a missing CanonicalFusionName (empty string reads
back as NaN, whose .str.split() also returns NaN -- a float, not a list --
so len(p) on the lambda's `p` raised TypeError rather than excluding the row).
"""
import csv
import os
import tempfile

import kernel


def _write_fusions(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ModelID", "IsDefaultEntryForModel", "CanonicalFusionName",
                    "Gene1", "Gene2"])
        for r in rows:
            w.writerow(r)


def test_canonical_only_skips_row_with_missing_fusion_name():
    d = tempfile.mkdtemp()
    _write_fusions(os.path.join(d, "OmicsFusionFiltered.csv"), [
        # CanonicalFusionName blank -> pandas reads NaN; must be excluded,
        # not crash.
        ["ACH-1", "Yes", "", "BCR (613)", "ABL1 (25)"],
        ["ACH-2", "Yes", "BCR--ABL1", "BCR (613)", "ABL1 (25)"],
    ])
    out = kernel.depmap_fusion_models("BCR", root=d, canonical_only=True)
    assert out == {"ACH-2"}


def test_canonical_only_keeps_named_partner_pairs():
    d = tempfile.mkdtemp()
    _write_fusions(os.path.join(d, "OmicsFusionFiltered.csv"), [
        ["ACH-1", "Yes", "BCR--ABL1", "BCR (613)", "ABL1 (25)"],
        ["ACH-2", "Yes", "BCR--BCR", "BCR (613)", "BCR (613)"],       # same-gene, dropped
        ["ACH-3", "Yes", "BCR--AC012345.1", "BCR (613)", "AC012345.1 (99)"],  # unnamed locus, dropped
    ])
    out = kernel.depmap_fusion_models("BCR", root=d, canonical_only=True)
    assert out == {"ACH-1"}


def test_fusion_contrast_untestable_when_too_few_fusion_positive_lines():
    d = tempfile.mkdtemp()
    with open(os.path.join(d, "CRISPRGeneEffect.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["", "GENE (1)"])
        for i in range(5):
            w.writerow(["ACH-%d" % i, -0.9])
    _write_fusions(os.path.join(d, "OmicsFusionFiltered.csv"), [
        ["ACH-0", "Yes", "GENE--X", "GENE (1)", "X (2)"],   # only 1 fusion-positive
        ["ACH-1", "Yes", "OTHER--Y", "OTHER (3)", "Y (4)"],
    ])
    out = kernel.depmap_fusion_contrast("GENE", root=d, min_n=3)
    assert out["status"] == "untestable"
    assert out["n_fusion_with_crispr"] == 1


def test_fusion_contrast_untestable_when_too_few_fusion_negative_lines():
    """Regression: an empty/undersized comparator must not fabricate a p-value.

    mannwhitneyu on a comparator this small returns NaN (SciPy). NaN > 0 is
    False in Python, so the old code's unconditional floor_pvalue(nan) floored
    it to sys.float_info.min -- a p-value that reads as maximally significant
    out of a comparison that was never actually computable. Almost every line
    is fusion-positive here, leaving too few fusion-negative comparator lines.
    """
    d = tempfile.mkdtemp()
    with open(os.path.join(d, "CRISPRGeneEffect.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["", "GENE (1)"])
        for i in range(5):
            w.writerow(["ACH-%d" % i, -0.9])
    rows = [["ACH-%d" % i, "Yes", "GENE--X", "GENE (1)", "X (2)"] for i in range(4)]
    rows.append(["ACH-4", "Yes", "OTHER--Y", "OTHER (3)", "Y (4)"])  # 1 fusion-negative
    _write_fusions(os.path.join(d, "OmicsFusionFiltered.csv"), rows)
    out = kernel.depmap_fusion_contrast("GENE", root=d, min_n=3)
    assert out["status"] == "untestable"
    assert out["n_fusion_with_crispr"] == 4
    assert out["n_fusion_negative"] == 1
    assert "p_fusion_more_dependent" not in out
