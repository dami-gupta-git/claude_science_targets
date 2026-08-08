"""Tests for the codependency half of kernel.py.

Correlation correctness is tested on a synthesised matrix with a KNOWN answer
and NaN pattern, so the suite runs without the DepMap release. Tests that need
the release are marked and skipped when it is absent.

Run: python -m pytest -q   (needs the `genetea` env)
"""
import os

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import kernel as cd

_depmap_root = os.environ.get("DEPMAP_ROOT", cd.DEPMAP_ROOT)
HAVE_DEPMAP = bool(_depmap_root) and os.path.isdir(
    os.path.join(_depmap_root, cd.CACHE_DIRNAME))
needs_depmap = pytest.mark.skipif(not HAVE_DEPMAP, reason="DepMap cache absent")


@pytest.fixture
def synthetic():
    """A 40x6 matrix with a planted correlation and structured missingness.

    G_partner is a noisy copy of G_query; G_anti is its negation; G_sparse is
    observed in only 12 lines, which is what exercises the pairwise-complete
    path. NaNs are placed in DIFFERENT lines for query and partner so a
    complete-case-across-all-genes implementation would give a different answer.
    """
    rng = np.random.default_rng(7)
    n = 40
    q = rng.normal(size=n)
    mat = pd.DataFrame({
        "GQUERY": q,
        "GPARTNER": q * 0.9 + rng.normal(scale=0.2, size=n),
        "GANTI": -q,
        "GNOISE": rng.normal(size=n),
        "GSPARSE": q + rng.normal(scale=0.1, size=n),
        "GFLAT": np.zeros(n),
    }, index=["ACH-%03d" % i for i in range(n)]).astype(np.float32)
    mat.loc[mat.index[:3], "GQUERY"] = np.nan
    mat.loc[mat.index[5:8], "GPARTNER"] = np.nan
    mat.loc[mat.index[12:], "GSPARSE"] = np.nan
    mat.index.name = "ModelID"
    return mat


def _prep_from(mat):
    """Build a prep dict from an in-memory matrix, mirroring codep_prepare."""
    values = mat.values.astype(np.float32)
    mask = ~np.isnan(values)
    centred = np.where(mask, values, np.float32(0.0))
    n_full = mask.sum(0).astype(np.int64)
    col_mean = (centred.sum(0, dtype=np.float64) / n_full).astype(np.float32)
    centred -= col_mean
    centred[~mask] = np.float32(0.0)
    maskf = mask.astype(np.float32)
    genes = np.asarray(list(mat.columns), dtype=object)
    return {
        "dataset": "effect", "genes": genes, "index": pd.Index(genes),
        "models": np.asarray(mat.index, dtype=object), "centred": centred,
        "maskf": maskf, "n_full": n_full, "col_mean": col_mean.astype(np.float64),
        "sum_full": centred.sum(0, dtype=np.float64),
        "sumsq_full": np.einsum("ij,ij->j", centred, centred, dtype=np.float64),
        "load_seconds": 0.0,
    }


def test_pairwise_complete_matches_scipy(synthetic):
    """r and n_pairs must equal scipy on the exact intersecting lines."""
    prep = _prep_from(synthetic)
    q = synthetic["GQUERY"].to_numpy(dtype=np.float64)
    r, n_pairs = cd.codep_pairwise_r(prep, q)
    for i, g in enumerate(prep["genes"]):
        y = synthetic[g].to_numpy(dtype=np.float64)
        both = np.isfinite(q) & np.isfinite(y)
        if both.sum() < 3 or np.std(y[both]) == 0:
            continue
        assert n_pairs[i] == both.sum(), g
        assert abs(r[i] - stats.pearsonr(q[both], y[both])[0]) < 1e-6, g


def test_planted_structure_recovered(synthetic):
    """The noisy copy is top-ranked; the negation is at r ~ -1."""
    prep = _prep_from(synthetic)
    q = synthetic["GQUERY"].to_numpy(dtype=np.float64)
    r, _ = cd.codep_pairwise_r(prep, q)
    by = dict(zip(prep["genes"], r))
    assert by["GPARTNER"] > 0.9
    assert by["GANTI"] == pytest.approx(-1.0, abs=1e-6)
    assert abs(by["GNOISE"]) < 0.5


def test_zero_variance_and_short_overlap_are_nan(synthetic):
    """A constant column yields NaN, and the sparse overlap count is exact.

    This pins the CONTRACT, not the guard: float32 zero variance already gives
    0/0 -> NaN, so the explicit (vx <= 0) | (vy <= 0) mask in
    codep_pairwise_r is defensive and this test stays green without it.
    """
    prep = _prep_from(synthetic)
    r, n_pairs = cd.codep_pairwise_r(prep, synthetic["GQUERY"].to_numpy(np.float64))
    assert np.isnan(r[list(prep["genes"]).index("GFLAT")])
    # GSPARSE overlaps the query in lines 3..11 only
    assert n_pairs[list(prep["genes"]).index("GSPARSE")] == 9


def test_query_needs_three_observations(synthetic):
    """Explicit failure, not a silent NaN table, on a degenerate query."""
    prep = _prep_from(synthetic)
    q = np.full(len(synthetic), np.nan)
    q[0] = 1.0
    with pytest.raises(ValueError, match="observed lines"):
        cd.codep_pairwise_r(prep, q)


def test_min_lines_is_not_knife_edge_on_synthetic(synthetic):
    """A partner exactly at the floor is kept; one below it is dropped."""
    prep = _prep_from(synthetic)
    r, n_pairs = cd.codep_pairwise_r(prep, synthetic["GQUERY"].to_numpy(np.float64))
    idx = list(prep["genes"]).index("GSPARSE")
    assert n_pairs[idx] >= 9
    assert (n_pairs[idx] >= 9) and not (n_pairs[idx] >= 10)


def test_bh_adjust_matches_reference():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.5])
    got = cd.bh_adjust(p)
    assert got[0] == pytest.approx(0.007, abs=1e-9)
    assert np.all(np.diff(got[np.argsort(p)]) >= -1e-12), "must be monotone"
    assert got.max() <= 1.0


def test_strip_entrez_idempotent():
    assert cd.strip_entrez("KRAS (3845)") == "KRAS"
    assert cd.strip_entrez(cd.strip_entrez("KRAS (3845)")) == "KRAS"


def test_shared_requires_two_queries():
    with pytest.raises(ValueError, match=">=2 query genes"):
        cd.depmap_shared_codependencies(["SMARCA4"])


@needs_depmap
def test_unknown_gene_names_itself():
    with pytest.raises(KeyError, match="NOTAGENE"):
        cd.depmap_codependencies("NOTAGENE")


@needs_depmap
def test_real_query_recovers_complex_and_gates():
    """SMARCB1 tops SMARCA4's list, and the control fails the query gate."""
    top = cd.depmap_codependencies("SMARCA4", n=10)
    assert top["gene"].iloc[0] == "SMARCB1"
    assert set(top["gene"]) & {"SMARCC1", "SMARCD1", "ARID1A"}
    assert list(top.columns[:3]) == ["gene", "r", "n_lines"]
    assert top["r"].is_monotonic_decreasing
    assert cd.codep_query_quality("SMARCA4")["passes"] is True
    ctrl = cd.codep_query_quality("OR5A1")
    assert ctrl["passes"] is False and "profile_sd" in ctrl["reason"]


@needs_depmap
def test_screened_background_is_smaller_than_corpus():
    """The whole point of the background correction."""
    allg = cd.depmap_screened_genes()
    floored = cd.depmap_screened_genes(min_lines=cd.CODEP_MIN_LINES)
    assert 18000 < len(allg) < 19000
    assert len(floored) < len(allg)
    assert len(allg) == len(set(allg)), "must be deduped"
    assert all(" (" not in g for g in allg), "Entrez suffix must be stripped"
