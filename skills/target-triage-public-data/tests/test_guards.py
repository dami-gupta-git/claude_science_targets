"""Guards and error paths.

These pin the behaviours a caller relies on when input is degenerate rather than
merely awkward: a partial correlation's p-value must account for the parameter
spent on the confound, and an input that makes the control meaningless must stop
the analysis instead of returning a number that looks like a controlled result.

The p-value case is checked against the t-distribution with n-3 degrees of
freedom computed here, rather than against scipy's pearsonr on the residuals,
because pearsonr's n-2 is precisely the value under test.
"""
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from kernel import (partial_corr, confounder_check, sensitivity_scan, load_prism,
                    load_prism_secondary, ot_essentiality_frame,
                    classify_gene_effect, DEPMAP_ESSENTIAL_THRESHOLD,
                    SCAN_COLUMNS)


def _write_csv(path, rows):
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        for r in rows:
            w.writerow(r)


# --- partial_corr p-value degrees of freedom ---------------------------------

@pytest.mark.parametrize("n", [20, 40, 200])
def test_partial_corr_p_uses_n_minus_3_degrees_of_freedom(n):
    rng = np.random.default_rng(7)
    z = rng.normal(size=n)
    x = z + rng.normal(size=n)
    y = 0.45 * x + z + rng.normal(size=n)

    r, p = partial_corr(x, y, z)
    t = r * np.sqrt((n - 3) / (1 - r ** 2))
    assert p == pytest.approx(2 * stats.t.sf(abs(t), n - 3), rel=1e-12)


def test_partial_corr_p_is_larger_than_the_naive_pearson_p_on_residuals():
    # The df error was anti-conservative: pearsonr on the residuals reports a
    # smaller p than the correct n-3 value. Asserting the direction means a
    # revert to `return stats.pearsonr(rx, ry)` fails here even if the
    # parametrised check above were somehow satisfied.
    n = 40
    rng = np.random.default_rng(7)
    z = rng.normal(size=n)
    x = z + rng.normal(size=n)
    y = 0.45 * x + z + rng.normal(size=n)

    r, p = partial_corr(x, y, z)
    rx = x - np.polyval(np.polyfit(z, x, 1), z)
    ry = y - np.polyval(np.polyfit(z, y, 1), z)
    naive_p = stats.pearsonr(rx, ry)[1]
    assert p > naive_p


def test_partial_corr_coefficient_unchanged_by_the_p_value_fix():
    # The fix must not move r: it is checked against the closed-form identity in
    # test_partial_corr.py, and a scan saved before the fix keeps its r_partial.
    rng = np.random.default_rng(0)
    n = 60
    z = rng.normal(size=n)
    x = 0.9 * z + rng.normal(scale=0.5, size=n)
    y = 0.7 * z + 0.6 * rng.normal(size=n)
    r, _ = partial_corr(x, y, z)
    rx = x - np.polyval(np.polyfit(z, x, 1), z)
    ry = y - np.polyval(np.polyfit(z, y, 1), z)
    assert r == pytest.approx(np.corrcoef(rx, ry)[0, 1], rel=1e-12)


# --- partial_corr degenerate input -------------------------------------------

def test_partial_corr_raises_on_constant_confound():
    # Previously returned the *uncontrolled* Pearson r (0.7746) as though z had
    # been removed, which is the one failure mode this skill must not have.
    with pytest.raises(ValueError, match="constant"):
        partial_corr([1., 2, 3, 4, 5], [2., 4, 5, 4, 5], [3., 3, 3, 3, 3])


def test_partial_corr_raises_on_non_finite_input():
    with pytest.raises(ValueError, match="finite"):
        partial_corr([1., 2, np.nan, 4], [1., 3, 2, 4], [1., 1, 2, 2])


def test_partial_corr_raises_below_four_observations():
    with pytest.raises(ValueError, match="n >= 4"):
        partial_corr([1., 2, 3], [1., 2, 4], [1., 3, 2])


def test_partial_corr_raises_on_mismatched_lengths():
    with pytest.raises(ValueError, match="same length"):
        partial_corr([1., 2, 3, 4], [1., 2, 3], [1., 2, 3, 4])


# --- sensitivity_scan guards --------------------------------------------------

def _scan_frame(n=30, seed=3, n_drugs=6):
    rng = np.random.default_rng(seed)
    df = pd.DataFrame({"expr": np.arange(1., n + 1)})
    for i in range(n_drugs):
        df["drug%d" % i] = rng.normal(size=n)
    return df


def test_sensitivity_scan_raises_on_reserved_column_name():
    # A drug column literally named _panel_mean used to be overwritten with the
    # panel mean and then correlated against expression as if it were a drug,
    # putting a quietly wrong row in the output table.
    df = _scan_frame(n_drugs=2)
    df["_panel_mean"] = np.arange(1., 31.)
    with pytest.raises(ValueError, match="reserved"):
        sensitivity_scan(df, "expr", ["drug0", "drug1", "_panel_mean"])


def test_sensitivity_scan_reserved_name_allowed_when_panel_mean_supplied():
    # With a caller-supplied confound the working columns are never written, so
    # the guard must not fire.
    df = _scan_frame(n_drugs=2)
    df["_panel_mean"] = np.arange(1., 31.)
    res = sensitivity_scan(df, "expr", ["drug0", "drug1"],
                           panel_mean_col="_panel_mean", min_n=5)
    assert len(res) == 2


def test_sensitivity_scan_raises_on_single_drug_panel():
    # With one drug column, "the whole-panel mean" IS that column, so the
    # computed confound would equal the drug itself. partial_corr would then
    # regress the drug on its own values and correlate expression against the
    # ~1e-16 floating-point residual noise of that perfect self-fit -- a
    # fabricated r_partial indistinguishable from a legitimate null result,
    # in the exact function this skill exists to keep from doing that.
    df = _scan_frame(n_drugs=1)
    with pytest.raises(ValueError, match=">=2 drug_cols"):
        sensitivity_scan(df, "expr", ["drug0"])


def test_sensitivity_scan_single_drug_panel_allowed_with_supplied_confound():
    # An externally supplied confound (e.g. from a wider panel) is not
    # self-referential, so a single-drug scan against it is legitimate.
    df = _scan_frame(n_drugs=1)
    df["_ext_confound"] = np.random.default_rng(1).normal(size=len(df))
    res = sensitivity_scan(df, "expr", ["drug0"],
                           panel_mean_col="_ext_confound", min_n=5)
    assert len(res) == 1


def test_sensitivity_scan_records_dropped_drugs():
    df = _scan_frame()
    df.loc[5:, "drug0"] = np.nan  # 5 surviving rows, below min_n
    res = sensitivity_scan(df, "expr", ["drug%d" % i for i in range(6)], min_n=15)

    assert "drug0" not in set(res.drug)
    dropped = res.attrs["dropped"]
    assert [d["drug"] for d in dropped] == ["drug0"]
    assert dropped[0]["n"] == 5
    assert "min_n" in dropped[0]["reason"]


def test_sensitivity_scan_empty_result_still_has_columns():
    # An all-dropped scan used to return a column-less frame, so res.r_partial
    # raised AttributeError instead of yielding an empty column. The BH columns
    # are asserted here too because they are added after the empty-result early
    # return, so they were the one part of the schema the guard did not cover.
    df = pd.DataFrame({"expr": [1., 2, 3], "drugA": [1., None, None],
                       "drugB": [2., None, None]})
    res = sensitivity_scan(df, "expr", ["drugA", "drugB"], min_n=15)
    assert len(res) == 0
    assert list(res.r_partial) == []
    assert list(res.q_BH) == []
    assert list(res.q_partial_BH) == []
    assert len(res.attrs["dropped"]) == 2


def test_sensitivity_scan_columns_match_between_empty_and_populated():
    # The two exit paths build the frame differently, so the schema a caller
    # filters on must be asserted identical rather than assumed.
    empty = sensitivity_scan(
        pd.DataFrame({"expr": [1., 2, 3], "drugA": [1., None, None],
                      "drugB": [2., None, None]}),
        "expr", ["drugA", "drugB"], min_n=15)
    populated = sensitivity_scan(_scan_frame(), "expr",
                                 ["drug%d" % i for i in range(6)])
    assert len(populated) > 0
    assert list(empty.columns) == list(populated.columns) == list(SCAN_COLUMNS)


def test_sensitivity_scan_raises_on_missing_drug_column():
    with pytest.raises(KeyError, match="absent"):
        sensitivity_scan(_scan_frame(), "expr", ["drug0", "nonexistent"])


def test_sensitivity_scan_raises_on_missing_supplied_panel_column():
    with pytest.raises(KeyError, match="panel_mean_col"):
        sensitivity_scan(_scan_frame(), "expr", ["drug0", "drug1"],
                         panel_mean_col="nope")


# --- confounder_check --------------------------------------------------------

def test_confounder_check_reports_thinnest_drug_coverage():
    # A panel mean resting on one measured drug is not a general-sensitivity
    # measure, and the returned r looks identical to a well-covered one, so the
    # coverage has to come back with it.
    n = 20
    df = pd.DataFrame({"expr": np.arange(1., n + 1)})
    for i in range(10):
        df["d%d" % i] = np.nan
    df["d0"] = np.arange(float(n), 0, -1)

    out = confounder_check(df, "expr", ["d%d" % i for i in range(10)])
    assert out["n"] == n
    assert out["n_drugs"] == 10
    assert out["n_drugs_min"] == 1
    assert out["r"] < -0.8


def test_confounder_check_raises_on_reserved_column():
    df = pd.DataFrame({"expr": [1., 2, 3], "d1": [1., 2, 3],
                       "_panel_mean": [1., 2, 3]})
    with pytest.raises(ValueError, match="reserved"):
        confounder_check(df, "expr", ["d1", "_panel_mean"])


def test_confounder_check_raises_when_too_few_complete_rows():
    df = pd.DataFrame({"expr": [1., None, None], "d1": [1., 2, 3]})
    with pytest.raises(ValueError, match="at least 3"):
        confounder_check(df, "expr", ["d1"])


# --- loaders -----------------------------------------------------------------

def test_load_prism_raises_on_duplicate_cell_line_ids(tmp_path):
    # A duplicated header id produced a duplicated index, so .loc returned a
    # Series where callers expect a scalar and merges multiplied rows.
    _write_csv(tmp_path / "c.csv", [["IDs", "Drug.Name"], ["BRD-A-1", "DRUGA"]])
    _write_csv(tmp_path / "m.csv", [["", "ACH-1", "ACH-1", "ACH-2"],
                                    ["BRD-A-1", "0.1", "0.5", "0.3"]])
    with pytest.raises(ValueError, match="repeats"):
        load_prism(str(tmp_path / "m.csv"), str(tmp_path / "c.csv"))


def test_load_prism_raises_on_truncated_row(tmp_path):
    _write_csv(tmp_path / "c.csv", [["IDs", "Drug.Name"], ["BRD-A-1", "DRUGA"]])
    _write_csv(tmp_path / "m.csv", [["", "ACH-1", "ACH-2", "ACH-3"],
                                    ["BRD-A-1", "0.1", "0.3"]])
    with pytest.raises(ValueError, match="truncated"):
        load_prism(str(tmp_path / "m.csv"), str(tmp_path / "c.csv"))


def test_load_prism_raises_when_no_brd_id_matches(tmp_path):
    _write_csv(tmp_path / "c.csv", [["IDs", "Drug.Name"], ["BRD-A-1", "DRUGA"]])
    _write_csv(tmp_path / "m.csv", [["", "ACH-1"], ["BRD-ZZZ", "0.1"]])
    with pytest.raises(ValueError, match="none of the"):
        load_prism(str(tmp_path / "m.csv"), str(tmp_path / "c.csv"))


def test_load_prism_raises_when_drug_names_select_nothing(tmp_path):
    _write_csv(tmp_path / "c.csv", [["IDs", "Drug.Name"], ["BRD-A-1", "DRUGA"]])
    _write_csv(tmp_path / "m.csv", [["", "ACH-1"], ["BRD-A-1", "0.1"]])
    with pytest.raises(ValueError, match="no compounds selected"):
        load_prism(str(tmp_path / "m.csv"), str(tmp_path / "c.csv"),
                   drug_names=["absent-drug"])


def test_load_prism_closes_the_compound_list_handle(tmp_path):
    import gc
    import io
    _write_csv(tmp_path / "c.csv", [["IDs", "Drug.Name"], ["BRD-A-1", "DRUGA"]])
    _write_csv(tmp_path / "m.csv", [["", "ACH-1"], ["BRD-A-1", "0.1"]])
    load_prism(str(tmp_path / "m.csv"), str(tmp_path / "c.csv"))
    gc.collect()
    open_paths = [o.name for o in gc.get_objects()
                  if isinstance(o, io.TextIOWrapper) and not o.closed
                  and isinstance(getattr(o, "name", None), str)]
    assert str(tmp_path / "c.csv") not in open_paths


def test_load_prism_secondary_loads_without_an_r2_column(tmp_path):
    # r2 is only needed for the min_r2 filter; requiring it unconditionally made
    # a file without that column fail on a call that never asked to filter.
    path = tmp_path / "s.csv"
    _write_csv(path, [["depmap_id", "name", "auc"], ["ACH-1", "drugx", "0.4"],
                      ["ACH-2", "drugx", "0.6"]])
    out = load_prism_secondary(str(path))
    assert out.loc["ACH-1", "drugx"] == 0.4


def test_load_prism_secondary_raises_when_drug_absent(tmp_path):
    path = tmp_path / "s.csv"
    _write_csv(path, [["depmap_id", "name", "auc", "r2"],
                      ["ACH-1", "trifluridine", "0.4", "0.9"]])
    with pytest.raises(ValueError, match="prism_secondary_has"):
        load_prism_secondary(str(path), drug_names=["decitabine"])


# --- Open Targets flattening -------------------------------------------------

def test_ot_essentiality_frame_tolerates_tissue_without_screens():
    df = ot_essentiality_frame({"approvedSymbol": "X",
                                "depMapEssentiality": [{"tissueName": "Lung"}]})
    assert len(df) == 0


def test_ot_essentiality_frame_empty_frame_keeps_columns():
    df = ot_essentiality_frame({"approvedSymbol": "X", "depMapEssentiality": None})
    assert list(df.columns) == ["gene", "tissue", "cell", "disease", "effect", "expr"]


def test_ot_essentiality_frame_rejects_non_dict():
    with pytest.raises(TypeError, match="ot_target"):
        ot_essentiality_frame(None)


# --- gene-effect classification ----------------------------------------------

def test_classify_gene_effect_boundaries_are_inclusive():
    # Exactly at a boundary must take the more essential label, and floating
    # point must not push the threshold value itself into the wrong bucket.
    assert classify_gene_effect(DEPMAP_ESSENTIAL_THRESHOLD) == "essential"
    assert classify_gene_effect(-1.0) == "essential"
    assert classify_gene_effect(-0.5) == "intermediate"
    assert classify_gene_effect(-0.499999) == "dispensable"


def test_classify_gene_effect_spans_the_documented_comparators():
    # Reference points from SKILL.md: RRM1/DUT around -2 to -3, KRAS -1.2 in CRC,
    # NUDT1 around 0.
    assert classify_gene_effect(-2.5) == "essential"
    assert classify_gene_effect(-1.2) == "essential"
    assert classify_gene_effect(-0.095) == "dispensable"
    assert classify_gene_effect(0.0) == "dispensable"


def test_classify_gene_effect_handles_missing_values():
    assert classify_gene_effect(None) == "unknown"
    assert classify_gene_effect(np.nan) == "unknown"
