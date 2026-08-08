"""partial_corr, confounder_check, sensitivity_scan.

partial_corr is THE load-bearing control this skill exists to enforce
correctly (SKILL.md: "the step with the confounder"), so it's checked
against the classic closed-form first-order partial correlation formula

    r_xy.z = (r_xy - r_xz*r_yz) / sqrt((1 - r_xz^2)(1 - r_yz^2))

computed independently via scipy.stats.pearsonr on the raw pairwise
correlations -- not by re-deriving the residual-regression approach the
implementation itself uses, so the test isn't just checking the code
against itself.
"""
import numpy as np
import pandas as pd
from scipy import stats

from kernel import partial_corr, confounder_check, sensitivity_scan


def _closed_form_partial(x, y, z):
    r_xy, _ = stats.pearsonr(x, y)
    r_xz, _ = stats.pearsonr(x, z)
    r_yz, _ = stats.pearsonr(y, z)
    return (r_xy - r_xz * r_yz) / np.sqrt((1 - r_xz**2) * (1 - r_yz**2))


def test_partial_corr_matches_closed_form_formula():
    rng = np.random.default_rng(0)
    n = 60
    z = rng.normal(size=n)
    x = 0.9 * z + rng.normal(scale=0.5, size=n)
    y = 0.7 * z + 0.6 * rng.normal(size=n) + 0.3 * rng.normal(scale=1, size=n)

    r, p = partial_corr(x, y, z)
    expected = _closed_form_partial(x, y, z)
    assert abs(r - expected) < 1e-9


def test_partial_corr_removes_shared_confound():
    # x and y are BOTH pure linear functions of z plus independent noise, with
    # no direct relationship to each other beyond what z induces. Raw r(x,y)
    # should be large (driven by z); partial r after removing z should
    # collapse toward zero. (With n=200 even a small residual r can clear
    # p<0.05, so only the magnitude of r_partial is asserted, not p.)
    rng = np.random.default_rng(1)
    n = 200
    z = rng.normal(size=n)
    x = 2.0 * z + rng.normal(scale=0.3, size=n)
    y = 3.0 * z + rng.normal(scale=0.3, size=n)

    raw_r, _ = stats.pearsonr(x, y)
    partial_r, partial_p = partial_corr(x, y, z)

    assert raw_r > 0.9
    assert abs(partial_r) < 0.25


def test_partial_corr_preserves_confound_independent_signal():
    # y has a genuine z-independent relationship to x on top of the shared
    # confound; partialling out z should NOT destroy that signal.
    rng = np.random.default_rng(2)
    n = 200
    z = rng.normal(size=n)
    x_indep = rng.normal(size=n)
    x = z + x_indep
    y = z + 1.5 * x_indep + rng.normal(scale=0.2, size=n)

    partial_r, partial_p = partial_corr(x, y, z)
    assert partial_r > 0.7
    assert partial_p < 1e-6


def test_confounder_check_detects_panel_wide_confound():
    rng = np.random.default_rng(3)
    n = 200
    z = rng.normal(size=n)
    df = pd.DataFrame({"expr": z + rng.normal(scale=0.2, size=n)})
    drug_cols = []
    for i in range(5):
        name = "drug%d" % i
        drug_cols.append(name)
        df[name] = z + rng.normal(scale=0.2, size=n)

    out = confounder_check(df, "expr", drug_cols)
    assert out["r"] > 0.7
    assert out["p"] < 1e-6
    assert out["n"] == n


def test_confounder_check_drops_nan_rows_and_reports_n():
    # pandas .mean(axis=1) skips NaN per-row, so a row survives the panel-mean
    # step unless EVERY drug column is NaN; only expr-NaN (row idx 2) and
    # all-drugs-NaN (row idx 3) rows should actually drop.
    df = pd.DataFrame({"expr": [1.0, 2.0, None, 4.0, 5.0],
                       "d1": [1.0, 2.0, 3.0, None, 5.0],
                       "d2": [2.0, 3.0, 4.0, None, 6.0]})
    out = confounder_check(df, "expr", ["d1", "d2"])
    assert out["n"] == 3


def test_confounder_check_keeps_row_with_partial_nan_drug_coverage():
    # A row missing ONE drug column but not all of them must still contribute
    # to the panel mean (pandas skips the NaN element, not the whole row).
    df = pd.DataFrame({"expr": [1.0, 2.0, 3.0],
                       "d1": [1.0, None, 3.0],
                       "d2": [2.0, 3.0, 4.0]})
    out = confounder_check(df, "expr", ["d1", "d2"])
    assert out["n"] == 3


def test_sensitivity_scan_separates_true_hit_from_confound_only_drug():
    # Supply the TRUE confound directly (rather than relying on a noisy
    # leave-one-out proxy, which only ever approximates it) to isolate the
    # thing this function exists to get right: a confound-only drug's
    # r_partial should collapse, a drug with genuine independent signal
    # should keep most of its correlation.
    rng = np.random.default_rng(4)
    n = 60
    z = rng.normal(size=n)
    x_indep = rng.normal(size=n)
    df = pd.DataFrame({"expr": z + x_indep, "true_z": z})
    df["confound_only"] = 2.0 * z + rng.normal(scale=0.2, size=n)
    df["true_hit"] = z + 1.8 * x_indep + rng.normal(scale=0.2, size=n)

    res = sensitivity_scan(df, "expr", ["confound_only", "true_hit"],
                           panel_mean_col="true_z")

    hit = res[res.drug == "true_hit"].iloc[0]
    confound = res[res.drug == "confound_only"].iloc[0]
    assert hit.r_partial > 0.7
    assert abs(confound.r_partial) < 0.3
    assert "q_BH" in res.columns and "q_partial_BH" in res.columns


def test_sensitivity_scan_excludes_drugs_below_min_n():
    # leave_one_out=False so the min_n filter is tested in isolation, without
    # the leave-one-out confound (itself built from the other drug) also
    # losing rows to sparse_drug's own missingness.
    df = pd.DataFrame({
        "expr": [1.0, 2.0, 3.0, 4.0, 5.0],
        "sparse_drug": [1.0, None, None, None, None],
        "dense_drug": [5.0, 4.0, 3.0, 2.0, 1.0],
    })
    res = sensitivity_scan(df, "expr", ["sparse_drug", "dense_drug"], min_n=3,
                           leave_one_out=False)
    assert "sparse_drug" not in set(res.drug)
    assert "dense_drug" in set(res.drug)


def test_sensitivity_scan_leave_one_out_excludes_the_tested_drug_from_its_own_confound():
    # With leave_one_out=True and only 2 drugs, drug A's confound is drug B
    # ALONE (not the mean of A and B). If B is mostly missing, A's confound
    # is mostly missing too, and A drops out even though A itself is dense -
    # this is the intended leave-one-out behavior, not a min_n bug.
    df = pd.DataFrame({
        "expr": [1.0, 2.0, 3.0, 4.0, 5.0],
        "sparse_drug": [1.0, None, None, None, None],
        "dense_drug": [5.0, 4.0, 3.0, 2.0, 1.0],
    })
    res = sensitivity_scan(df, "expr", ["sparse_drug", "dense_drug"], min_n=3,
                           leave_one_out=True)
    assert len(res) == 0


def test_sensitivity_scan_supplied_panel_mean_disables_leave_one_out():
    # A caller-supplied panel_mean_col is a fixed external confound (its
    # constituent columns are unknown to sensitivity_scan), so leave_one_out
    # must be ignored rather than attempted.
    df = pd.DataFrame({
        "expr": [1.0, 2.0, 3.0, 4.0, 5.0],
        "sparse_drug": [1.0, None, None, None, None],
        "dense_drug": [5.0, 4.0, 3.0, 2.0, 1.0],
        "external_panel_mean": [3.0, 3.5, 3.0, 2.5, 2.0],
    })
    res = sensitivity_scan(df, "expr", ["sparse_drug", "dense_drug"], min_n=3,
                           panel_mean_col="external_panel_mean")
    assert "dense_drug" in set(res.drug)
    assert (res.confound == "whole-panel").all()


def test_sensitivity_scan_sorted_by_r_partial_ascending():
    rng = np.random.default_rng(5)
    n = 50
    df = pd.DataFrame({"expr": rng.normal(size=n)})
    for i in range(5):
        df["drug%d" % i] = rng.normal(size=n)
    res = sensitivity_scan(df, "expr", ["drug%d" % i for i in range(5)])
    assert list(res.r_partial) == sorted(res.r_partial)
