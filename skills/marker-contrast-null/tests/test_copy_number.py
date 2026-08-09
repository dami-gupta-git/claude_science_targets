"""Tests for the copy-number deletion null.

All inputs are synthesised: a 60-sample panel where samples 0-19 carry a
deletion spanning FOCAL and four contiguous neighbours, plus one independent
marker deleted in a different, overlapping-but-distinct set of samples.
"""

import numpy as np
import pandas as pd
import pytest

from kernel import (
    CODELETION_MIN_OVERLAP,
    DELETION_CUT,
    codeletion_partners,
    deletion_marker_matrix,
    marker_null_scan,
    neighbourhood_check,
)

N = 60
DELETED = list(range(20))          # samples carrying the FOCAL-locus deletion
INDEPENDENT = list(range(40, 58))  # a separate locus, disjoint from DELETED
NEIGHBOURS = ["NB1", "NB2", "NB3", "NB4"]


@pytest.fixture
def cn():
    """Relative copy number: 1.0 diploid, 0.02 homozygously deleted."""
    idx = [f"S{i:02d}" for i in range(N)]
    frame = pd.DataFrame(1.0, index=idx, columns=["FOCAL"] + NEIGHBOURS + ["INDEP", "FLAT"])
    for col in ["FOCAL"] + NEIGHBOURS:
        frame.iloc[DELETED, frame.columns.get_loc(col)] = 0.02
    frame.iloc[INDEPENDENT, frame.columns.get_loc("INDEP")] = 0.02
    return frame


@pytest.fixture
def markers(cn):
    return deletion_marker_matrix(cn)


@pytest.fixture
def values(markers):
    """Dependency that really does track the FOCAL locus, plus noise."""
    rng = np.random.default_rng(0)
    base = pd.Series(rng.normal(-1.0, 0.25, N), index=markers.index)
    return base - 0.6 * markers["FOCAL"].astype(float)


class TestDeletionMarkerMatrix:
    def test_calls_deletions_below_the_cut(self, cn, markers):
        assert markers["FOCAL"].sum() == len(DELETED)
        assert markers["FLAT"].sum() == 0
        assert markers.dtypes.unique().tolist() == [np.dtype(bool)]

    def test_default_cut_is_the_documented_constant(self, cn):
        assert DELETION_CUT == 0.25
        loose = deletion_marker_matrix(cn.assign(FLAT=0.2))
        assert loose["FLAT"].sum() == len(cn)

    def test_log2_ratio_input_raises_by_calling_everything_deleted(self):
        # log2 ratios centre on 0, so every value sits below a 0.25 relative-CN
        # cut and the whole matrix would be called deleted.
        log2like = pd.DataFrame(0.0, index=[f"S{i}" for i in range(10)],
                                columns=["G1", "G2"])
        with pytest.raises(ValueError, match="(?i)relative copy number"):
            deletion_marker_matrix(log2like)

    def test_absolute_copy_number_input_raises_by_calling_nothing(self):
        # Absolute integer copy number never falls below a relative-scale cut.
        absolute = pd.DataFrame(2.0, index=[f"S{i}" for i in range(10)],
                                columns=["G1", "G2"])
        with pytest.raises(ValueError, match="(?i)relative copy number"):
            deletion_marker_matrix(absolute)

    def test_all_nan_frame_raises_its_own_error(self):
        with pytest.raises(ValueError, match="no numeric"):
            deletion_marker_matrix(pd.DataFrame({"G1": ["x", "y"]}))

    def test_non_numeric_columns_do_not_crash_the_comparison(self, cn):
        withtext = cn.assign(NOTE="hello")
        assert deletion_marker_matrix(withtext)["NOTE"].sum() == 0


class TestCodeletionPartners:
    def test_contiguous_neighbours_are_found(self, markers):
        found = set(codeletion_partners(markers, "FOCAL")["marker"])
        assert found == set(NEIGHBOURS)

    def test_independent_locus_is_excluded(self, markers):
        assert "INDEP" not in set(codeletion_partners(markers, "FOCAL")["marker"])

    def test_partial_overlap_below_threshold_is_excluded(self, cn):
        # HALF shares 10 of FOCAL's 20 samples and adds 10 of its own, so
        # containment is 10/20 = 0.5 in both directions — under the 0.7 default.
        cn = cn.copy()
        cn["HALF"] = 1.0
        cn.iloc[list(range(10, 30)), cn.columns.get_loc("HALF")] = 0.02
        markers = deletion_marker_matrix(cn)
        assert 0.5 < CODELETION_MIN_OVERLAP
        assert "HALF" not in set(codeletion_partners(markers, "FOCAL")["marker"])
        assert "HALF" in set(codeletion_partners(markers, "FOCAL", min_overlap=0.4)["marker"])

    def test_a_nested_narrow_deletion_is_still_a_partner(self, cn):
        # The case Jaccard gets wrong: NARROW is deleted in 6 of FOCAL's 20
        # samples and nowhere else. Containment = 6/6 = 1.0; Jaccard = 6/20 = 0.3.
        cn = cn.copy()
        cn["NARROW"] = 1.0
        cn.iloc[DELETED[:6], cn.columns.get_loc("NARROW")] = 0.02
        markers = deletion_marker_matrix(cn)
        partners = codeletion_partners(markers, "FOCAL").set_index("marker")
        assert "NARROW" in partners.index
        assert partners.loc["NARROW", "overlap"] == 1.0
        assert 6 / 20 < CODELETION_MIN_OVERLAP   # Jaccard would have excluded it

    def test_absent_marker_raises(self, markers):
        with pytest.raises(KeyError, match="NOPE"):
            codeletion_partners(markers, "NOPE")

    def test_marker_positive_in_no_sample_raises(self, markers):
        with pytest.raises(ValueError, match="positive in no sample"):
            codeletion_partners(markers, "FLAT")


class TestNeighbourhoodCheck:
    def test_codeleted_neighbours_are_not_counted_as_independent(self, cn):
        # Give the neighbours a slightly WIDER deletion than FOCAL so they score
        # marginally stronger and genuinely sort above it. Without this the
        # columns tie, FOCAL may rank first, and the assertion below is vacuous.
        cn = cn.copy()
        for col in NEIGHBOURS:
            cn.iloc[DELETED + [20, 21], cn.columns.get_loc(col)] = 0.02
        markers = deletion_marker_matrix(cn)
        rng = np.random.default_rng(0)
        vals = pd.Series(rng.normal(-1.0, 0.15, N), index=markers.index)
        vals -= 0.8 * markers["NB1"].astype(float)
        scan = marker_null_scan(vals, markers, min_n=5, max_n=40)
        out = neighbourhood_check(scan, markers, "FOCAL")
        assert out["n_above"] > 0, "test is vacuous unless something outranks FOCAL"
        assert set(out["codeleted_above"]) <= set(NEIGHBOURS)
        assert out["n_independent_above"] == 0
        assert out["n_above"] == out["n_codeleted_above"] + out["n_independent_above"]

    def test_an_independent_stronger_marker_is_counted(self, cn):
        # Make INDEP the true driver: it must show up as independent_above.
        markers = deletion_marker_matrix(cn)
        rng = np.random.default_rng(1)
        vals = pd.Series(rng.normal(-1.0, 0.1, N), index=markers.index)
        vals -= 0.3 * markers["FOCAL"].astype(float)
        vals -= 1.5 * markers["INDEP"].astype(float)
        scan = marker_null_scan(vals, markers, min_n=5, max_n=40)
        out = neighbourhood_check(scan, markers, "FOCAL")
        assert "INDEP" in out["independent_above"]
        assert out["n_independent_above"] >= 1

    def test_untested_marker_raises_through_rank_in_null(self, markers, values):
        scan = marker_null_scan(values, markers, min_n=5, max_n=40)
        with pytest.raises(KeyError):
            neighbourhood_check(scan, markers, "NOPE")
