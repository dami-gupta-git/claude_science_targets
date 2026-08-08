"""depmap_mutation_contrast's group-size guard.

Both arms of the mutant/WT comparison must be checked before a statistic is
computed. A one-sided guard let a too-small WT arm through: Cohen's d off a
single-element std() is NaN (Python truthiness doesn't catch NaN, so the
`if pooled else None` fallback does not fire), while a p-value is still
reported alongside it -- silently wrong rather than a clear "too few lines".
"""
import csv
import os
import tempfile

import kernel


def _write_effect_csv(path, gene, values):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["", "%s (1)" % gene])
        for i, v in enumerate(values):
            w.writerow(["ACH-%d" % i, v])


def _write_hotspot_csv(path, marker_gene, mutant_flags):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ModelID", "IsDefaultEntryForModel", "%s (1)" % marker_gene])
        for i, m in enumerate(mutant_flags):
            w.writerow(["ACH-%d" % i, "Yes", int(m)])


def _setup(effect_values, mutant_flags):
    d = tempfile.mkdtemp()
    _write_effect_csv(os.path.join(d, "CRISPRGeneEffect.csv"), "DEPGENE", effect_values)
    _write_hotspot_csv(os.path.join(d, "OmicsSomaticMutationsMatrixHotspot.csv"),
                       "MARKERGENE", mutant_flags)
    return d


def test_too_few_wt_lines_returns_note_not_nan():
    # 4 mutant, 1 WT: mutant arm clears the n>=3 floor, WT does not.
    root = _setup([-0.9, -0.8, -0.85, -0.1, -0.05], [1, 1, 1, 1, 0])
    out = kernel.depmap_mutation_contrast("DEPGENE", "MARKERGENE", root=root)
    assert out["note"] == "too few WT lines"
    assert out["n_mutant"] == 4
    assert out["n_wt"] == 1
    assert "cohens_d" not in out
    assert "p_mutant_more_dependent" not in out


def test_too_few_mutant_lines_returns_note():
    # Symmetric case, mutant arm short instead -- pins the existing branch.
    root = _setup([-0.9, -0.1, -0.05, -0.1, -0.05], [1, 0, 0, 0, 0])
    out = kernel.depmap_mutation_contrast("DEPGENE", "MARKERGENE", root=root)
    assert out["note"] == "too few mutant lines"
    assert out["n_mutant"] == 1


def test_both_arms_sufficient_returns_full_stats():
    root = _setup([-0.9, -0.8, -0.85, -0.1, -0.05, -0.15], [1, 1, 1, 0, 0, 0])
    out = kernel.depmap_mutation_contrast("DEPGENE", "MARKERGENE", root=root)
    assert out["n_mutant"] == 3
    assert out["n_wt"] == 3
    assert "note" not in out
    import math
    assert not math.isnan(out["cohens_d"])
    assert 0.0 <= out["p_mutant_more_dependent"] <= 1.0
