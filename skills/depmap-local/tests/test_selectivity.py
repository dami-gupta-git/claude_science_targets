"""depmap_selectivity's classification thresholds, using synthetic
CRISPRGeneEffect.csv fixtures with known frac_dependent values so the
0.90/0.02 boundaries can be pinned exactly (DEP_THRESHOLD = -0.5)."""
import csv

import kernel


def _write_effect_csv(path, gene, values):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["", "%s (1)" % gene])
        for i, v in enumerate(values):
            w.writerow(["ACH-%d" % i, v])


def test_classification_selective_midrange():
    values = ([-1.0] * 5) + ([0.0] * 5)  # 5/10 = 0.5 dependent
    stats = _classify(values)
    assert stats["classification"] == "selective"
    assert stats["frac_dependent"] == 0.5
    assert stats["n_dependent"] == 5
    assert stats["n_lines"] == 10
    assert stats["mean_effect"] == round(sum(values) / len(values), 3)
    assert stats["min_effect"] == -1.0


def test_classification_common_essential_at_boundary():
    # Exactly 0.90 -> common essential (boundary is inclusive, ">= 0.90").
    values = ([-1.0] * 9) + ([0.0] * 1)
    stats = _classify(values)
    assert stats["frac_dependent"] == 0.9
    assert stats["classification"] == "common essential"


def test_classification_selective_at_lower_boundary():
    # Exactly 0.02 -> selective (boundary is inclusive, ">= 0.02"), not
    # non-essential.
    values = ([-1.0] * 1) + ([0.0] * 49)  # 1/50 = 0.02
    stats = _classify(values)
    assert stats["frac_dependent"] == 0.02
    assert stats["classification"] == "selective"


def test_classification_non_essential_below_lower_boundary():
    values = ([-1.0] * 1) + ([0.0] * 199)  # 1/200 = 0.005 < 0.02
    stats = _classify(values)
    assert stats["frac_dependent"] == 0.005
    assert stats["classification"] == "non-essential"


def test_unknown_gene_raises_keyerror(tmp_path):
    _write_effect_csv(tmp_path / "CRISPRGeneEffect.csv", "KRAS", [-0.7, -0.2])
    import pytest
    with pytest.raises(KeyError):
        kernel.depmap_selectivity("NOTAGENE", root=str(tmp_path))


def _classify(values):
    import tempfile
    import os
    d = tempfile.mkdtemp()
    _write_effect_csv(os.path.join(d, "CRISPRGeneEffect.csv"), "GENE", values)
    return kernel.depmap_selectivity("GENE", root=d)
