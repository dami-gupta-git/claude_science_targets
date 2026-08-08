"""Tests for the marker-null control helpers.

Inputs are synthesised so the suite needs no DepMap release, no network and no
credentials. Two cases carry the skill's reason for existing: a planted marker
must rank first in its own scan (`test_planted_marker_ranks_first`), and a
stratum that is globally shifted must be caught by the global-shift control
rather than reported as a target-specific effect
(`test_global_shift_control_detects_a_globally_shifted_stratum`).
"""
import numpy as np
import pandas as pd
import pytest

from kernel import (MARKER_MAX_N, MARKER_MIN_N, bh_q, cohens_d,
                    gene_specificity_control, global_shift_control,
                    marker_null_scan, rank_in_null, sample_gene_means,
                    stratum_contrast)


def samples(n=200):
    return [f"S{i:03d}" for i in range(n)]


def marker_matrix(index, spec):
    """spec: {marker: number of positive samples} -- first k samples get a hit."""
    frame = pd.DataFrame(0, index=index, columns=list(spec))
    for marker, k in spec.items():
        frame.loc[index[:k], marker] = 1
    return frame


def test_cohens_d_sign_and_magnitude():
    # Both arms SD=1 (ddof=1) and means 3 apart -> d = -3.
    focal = pd.Series([-3.0, -4.0, -2.0])
    ref = pd.Series([0.0, -1.0, 1.0])
    assert cohens_d(focal, ref) == pytest.approx(-3.0)
    # order reverses the sign
    a = pd.Series([1.0, 2.0, 3.0, 4.0])
    b = pd.Series([2.0, 3.0, 4.0, 5.0])
    assert cohens_d(a, b) == pytest.approx(-cohens_d(b, a))


def test_cohens_d_raises_when_both_arms_are_constant():
    with pytest.raises(ValueError, match="pooled SD is zero"):
        cohens_d([1.0, 1.0, 1.0], [2.0, 2.0, 2.0])


def test_cohens_d_drops_nan_and_requires_two_per_arm():
    assert cohens_d([1.0, 2.0, np.nan, 3.0], [2.0, 3.0, 4.0]) == pytest.approx(
        cohens_d([1.0, 2.0, 3.0], [2.0, 3.0, 4.0]))
    with pytest.raises(ValueError, match="need >=2 per arm"):
        cohens_d([1.0], [1.0, 2.0, 3.0])


def test_bh_q_preserves_order_and_matches_hand_calculation():
    p = np.array([0.01, 0.5, 0.02, 0.9])
    q = bh_q(p)
    # sorted 0.01,0.02,0.5,0.9 -> raw 0.04,0.04,0.667,0.9, monotone from the top
    assert q[0] == pytest.approx(0.04)
    assert q[2] == pytest.approx(0.04)
    assert q[1] == pytest.approx(0.6667, abs=1e-4)
    assert q[3] == pytest.approx(0.9)
    assert bh_q([0.5]).tolist() == [0.5]
    assert (bh_q([0.9, 0.95]) <= 1.0).all()


def test_stratum_contrast_reports_both_arms():
    rng = np.random.default_rng(5)
    idx = samples(60)
    values = pd.Series(rng.normal(0, 0.2, 60), index=idx)
    values.iloc[:20] -= 1.0
    flags = pd.Series([True] * 20 + [False] * 40, index=idx)
    out = stratum_contrast(values, flags)
    assert (out["n_focal"], out["n_reference"]) == (20, 40)
    assert out["mean_focal"] == pytest.approx(-1.0, abs=0.15)
    assert out["d"] < 0 and out["p"] < 0.001


def test_stratum_contrast_excludes_samples_with_no_flag_entry():
    # A sample in `values` but absent from `flags`' index was never genotyped
    # for this marker -- it must be excluded, not silently counted as
    # confirmed wild-type. A prior version reindexed flags onto values and
    # filled missing entries with False, diluting the reference arm with
    # untested samples and reporting no signal that it had happened.
    idx = samples(20)
    # Only S000-S009 were genotyped for this marker (5 focal, 5 reference);
    # S010-S019 have no call at all (not: called and found wild-type), and
    # their values are an extreme outlier that must not leak into either
    # arm's mean.
    genotyped = [-2.5, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, 2.5]
    values = pd.Series(genotyped + [100.0] * 10, index=idx)
    flags = pd.Series([True] * 5 + [False] * 5, index=idx[:10])

    out = stratum_contrast(values, flags)
    assert out["n_focal"] == 5
    assert out["n_reference"] == 5
    assert out["n_excluded_no_flag"] == 10
    assert out["mean_focal"] == pytest.approx(-1.5)
    assert out["mean_reference"] == pytest.approx(1.5)

    # Fully-aligned callers (the normal case, e.g. marker_null_scan's
    # internal use) see no exclusions.
    aligned = stratum_contrast(values, pd.Series([True] * 10 + [False] * 10, index=idx))
    assert aligned["n_excluded_no_flag"] == 0


def test_stratum_contrast_enforces_min_arm_n():
    # A direct call (SKILL.md's documented workflow step 1) must not silently
    # return an underpowered d/p below MIN_ARM_N -- marker_null_scan already
    # pre-filters to this floor, but stratum_contrast itself did not enforce
    # it, so a caller bypassing the scan got no signal the result was too
    # small to interpret as a standardised mean difference.
    idx = samples(10)
    values = pd.Series(np.linspace(-2, 2, 10), index=idx)
    flags = pd.Series([True] * 2 + [False] * 8, index=idx)
    with pytest.raises(ValueError, match="MIN_ARM_N"):
        stratum_contrast(values, flags)


def test_stratum_contrast_alternative_flips_the_tail():
    idx = samples(60)
    values = pd.Series(np.linspace(-1, 1, 60), index=idx)
    flags = pd.Series([True] * 20 + [False] * 40, index=idx)
    low = stratum_contrast(values, flags, alternative="less")
    high = stratum_contrast(values, flags, alternative="greater")
    assert low["p"] < 0.05 < high["p"]
    assert low["d"] == pytest.approx(high["d"])


def test_planted_marker_ranks_first():
    # One marker's carriers are shifted; the rest are noise. The planted marker
    # must top the scan -- this is what the null is for.
    rng = np.random.default_rng(0)
    idx = samples(200)
    values = pd.Series(rng.normal(0, 1, 200), index=idx)
    values.iloc[:30] -= 2.0
    spec = {"PLANTED": 30}
    spec.update({f"NOISE{i}": 30 for i in range(1, 15)})
    matrix = marker_matrix(pd.Index(idx), spec)
    # give the noise markers different carriers so they are not aliases
    for i in range(1, 15):
        matrix[f"NOISE{i}"] = 0
        matrix.iloc[30 + i: 60 + i, matrix.columns.get_loc(f"NOISE{i}")] = 1

    scan = marker_null_scan(values, matrix)
    assert scan.iloc[0].marker == "PLANTED"
    hit = rank_in_null(scan, "PLANTED")
    assert hit["rank_by_d"] == 1
    assert hit["n_markers"] == len(scan)
    assert hit["fraction_more_extreme"] == 0.0
    assert hit["q"] < 0.05


def test_null_scan_finds_no_hit_when_nothing_is_planted():
    rng = np.random.default_rng(1)
    idx = samples(200)
    values = pd.Series(rng.normal(0, 1, 200), index=idx)
    matrix = pd.DataFrame(
        {f"M{i}": (pd.Index(idx).isin(idx[i * 3: i * 3 + 30])).astype(int)
         for i in range(20)}, index=idx)
    scan = marker_null_scan(values, matrix)
    assert len(scan) == 20
    assert (scan.q >= 0.05).all(), scan[scan.q < 0.05]


def test_null_scan_q_is_corrected_not_raw_p():
    # The scan's whole purpose is the multiplicity a single marker hides, so q
    # must differ from p. Built so at least one marker has raw p < 0.05 while
    # its BH q is above it -- returning p unchanged would pass the other tests.
    rng = np.random.default_rng(11)
    idx = samples(200)
    values = pd.Series(rng.normal(0, 1, 200), index=idx)
    matrix = pd.DataFrame(
        {f"M{i}": (pd.Index(idx).isin(idx[i * 4: i * 4 + 30])).astype(int)
         for i in range(30)}, index=idx)
    scan = marker_null_scan(values, matrix)
    assert (scan.q > scan.p).any(), "q never exceeds p - BH not applied"
    assert (scan.q >= scan.p - 1e-12).all(), "q below p is impossible under BH"
    assert (scan.p < 0.05).any() and (scan.q < 0.05).sum() < (scan.p < 0.05).sum()


def test_null_scan_respects_the_size_band():
    idx = samples(300)
    values = pd.Series(np.zeros(300), index=idx)
    matrix = marker_matrix(pd.Index(idx), {"TOO_RARE": MARKER_MIN_N - 1,
                                          "JUST_IN": MARKER_MIN_N,
                                          "TOO_COMMON": MARKER_MAX_N + 1,
                                          "AT_CEILING": MARKER_MAX_N})
    values.iloc[:] = np.linspace(-1, 1, 300)
    scan = marker_null_scan(values, matrix)
    assert set(scan.marker) == {"JUST_IN", "AT_CEILING"}


def test_null_scan_errors_when_nothing_is_eligible():
    idx = samples(50)
    values = pd.Series(np.zeros(50), index=idx)
    matrix = marker_matrix(pd.Index(idx), {"RARE": 2})
    with pytest.raises(ValueError, match="no marker had"):
        marker_null_scan(values, matrix)


def test_null_scan_errors_on_disjoint_sample_ids():
    values = pd.Series([0.0, 1.0, 2.0], index=["A", "B", "C"])
    matrix = pd.DataFrame({"M": [1, 0, 1]}, index=["X", "Y", "Z"])
    with pytest.raises(ValueError, match="share no sample ids"):
        marker_null_scan(values, matrix)


def test_rank_in_null_raises_for_an_untested_marker():
    idx = samples(120)
    values = pd.Series(np.linspace(-1, 1, 120), index=idx)
    matrix = marker_matrix(pd.Index(idx), {"TESTED": 30, "RARE": 3})
    scan = marker_null_scan(values, matrix)
    assert rank_in_null(scan, "TESTED")["rank_by_d"] == 1
    with pytest.raises(KeyError, match="outside the stratum-size band"):
        rank_in_null(scan, "RARE")


def test_global_shift_control_detects_a_globally_shifted_stratum():
    idx = samples(120)
    means = pd.Series(np.zeros(120), index=idx)
    means.iloc[:30] = -0.8                      # globally more sensitive
    flags = pd.Series(pd.Index(idx).isin(idx[:30]), index=idx)
    shifted = global_shift_control(means, flags)
    assert shifted["global_shift"] is True and shifted["p"] < 0.05

    flat = global_shift_control(pd.Series(np.linspace(-0.1, 0.1, 120), index=idx),
                               pd.Series(pd.Index(idx).isin(idx[45:75]), index=idx))
    assert flat["global_shift"] is False and flat["p"] > 0.05


def test_gene_specificity_control_ranks_unrelated_genes_alongside():
    rng = np.random.default_rng(3)
    idx = samples(120)
    flags = pd.Series(pd.Index(idx).isin(idx[:30]), index=idx)
    noise = lambda: rng.normal(0, 1, 120)
    frame = pd.DataFrame({"TARGET": noise(), "NEIGHBOUR": noise(),
                          "UNRELATED": noise()}, index=idx)
    frame.loc[idx[:30], "TARGET"] -= 1.0
    frame.loc[idx[:30], "NEIGHBOUR"] -= 2.0     # moves MORE than the target
    out = gene_specificity_control(frame, flags)
    assert list(out.gene[:2]) == ["NEIGHBOUR", "TARGET"]
    unrelated = out.set_index("gene").loc["UNRELATED"]
    assert abs(unrelated.d) < 0.5 and unrelated.p > 0.05


def test_sample_gene_means_is_seed_stable_and_checks_pool_size():
    idx = samples(40)
    pool = [f"G{i}" for i in range(50)]

    def read_matrix(kind, columns):
        assert kind == "effect"
        rng = np.random.default_rng(abs(hash(tuple(columns))) % 2**32)
        return pd.DataFrame(rng.normal(0, 1, (40, len(columns))),
                            index=idx, columns=columns), []

    first = sample_gene_means(read_matrix, pool, n_genes=10, seed=7)
    again = sample_gene_means(read_matrix, pool, n_genes=10, seed=7)
    pd.testing.assert_series_equal(first, again)
    assert len(first) == 40
    with pytest.raises(ValueError, match="need 999"):
        sample_gene_means(read_matrix, pool, n_genes=999)
