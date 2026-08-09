"""Coverage for the results/<topic>/<run>/ writer added to marker-contrast-null.

Uses a tmp_path root rather than the real results/ tree, and synthetic
stratum_contrast()/rank_in_null()-shaped dicts rather than a live DepMap panel.
"""
import os

import pandas as pd
import pytest

from kernel import (mcn_link_into, mcn_run_dir, mcn_run_readme, mcn_verdict,
                    mcn_write_run, results_root)


def test_results_root_raises_when_unconfigured(monkeypatch):
    monkeypatch.delenv("SCIENCE_RESULTS_ROOT", raising=False)
    with pytest.raises(FileNotFoundError, match="SCIENCE_RESULTS_ROOT"):
        results_root()


def test_run_dir_raises_and_creates_nothing_when_root_unconfigured(monkeypatch, tmp_path):
    monkeypatch.delenv("SCIENCE_RESULTS_ROOT", raising=False)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError):
        mcn_run_dir("USP1 in BRCA1-mutant lines")
    assert not (tmp_path / "results").exists()


def test_run_dir_honours_env_var(monkeypatch, tmp_path):
    monkeypatch.setenv("SCIENCE_RESULTS_ROOT", str(tmp_path))
    out_dir = mcn_run_dir("USP1 in BRCA1-mutant lines")
    assert out_dir == str(tmp_path / "marker_contrast_null" / "usp1_in_brca1_mutant_lines")


def contrast(n_focal=24, n_reference=1184, mean_focal=-0.424,
            mean_reference=-0.300, d=-0.478, p=0.007):
    return {"n_focal": n_focal, "n_reference": n_reference,
            "n_excluded_no_flag": 0, "mean_focal": mean_focal,
            "mean_reference": mean_reference, "d": d, "p": p}


def rank(marker="BRCA1", rank_by_d=67, n_markers=1719, percentile=3.9,
        d=-0.478, p=0.007, q=0.99, n_q_below=0):
    return {"marker": marker, "rank_by_d": rank_by_d, "n_markers": n_markers,
            "percentile": percentile, "d": d, "p": p, "q": q,
            "n_markers_q_below_05": n_q_below, "fraction_more_extreme": 0.04}


def test_run_dir_slugs_name_and_makes_scripts(tmp_path):
    out_dir = mcn_run_dir("USP1 in BRCA1-mutant lines", root=str(tmp_path))
    assert out_dir == str(tmp_path / "marker_contrast_null" / "usp1_in_brca1_mutant_lines")
    assert os.path.isdir(os.path.join(out_dir, "scripts"))


def test_run_dir_rejects_empty_name(tmp_path):
    with pytest.raises(ValueError):
        mcn_run_dir("   ", root=str(tmp_path))


def test_verdict_does_not_survive_when_q_above_floor():
    v = mcn_verdict(rank())
    assert v["survives"] is False
    assert v["reasons"]


def test_verdict_survives_when_q_below_floor_and_no_global_shift():
    v = mcn_verdict(rank(q=0.01, n_q_below=5),
                    global_shift={"global_shift": False, "p": 0.34})
    assert v["survives"] is True
    assert v["reasons"] == []


def test_verdict_fails_on_global_shift_even_with_low_q():
    v = mcn_verdict(rank(q=0.01), global_shift={"global_shift": True, "p": 0.001})
    assert v["survives"] is False


def test_write_run_does_not_survive_case(tmp_path):
    out_dir = mcn_run_dir("USP1 in BRCA1-mutant lines", root=str(tmp_path))
    scan = pd.DataFrame([
        {"marker": "BRCA1", "n_focal": 24, "n_reference": 1184,
         "mean_focal": -0.424, "mean_reference": -0.300, "d": -0.478,
         "p": 0.007, "q": 0.99},
        {"marker": "FANCD2", "n_focal": 12, "n_reference": 1196,
         "mean_focal": -0.6, "mean_reference": -0.3, "d": -0.643,
         "p": 0.001, "q": 0.4},
    ])
    written = mcn_write_run(
        out_dir, "USP1 in BRCA1-mutant lines", contrast(), rank(), scan=scan,
        global_shift={"global_shift": False, "p": 0.34},
        summary="BRCA1-damaging lines look more USP1-dependent, but the null "
                "scan and global-shift control both contradict it.",
        data_sources=["DepMap 24Q2"])
    assert os.path.isfile(written["scan"])
    assert os.path.isfile(written["readme"])
    text = open(written["readme"]).read()
    assert "does not survive" in text
    assert "## Limits" in text and "## Files" in text and "## Data sources" in text


def test_write_run_survives_case_has_no_reasons_block(tmp_path):
    out_dir = mcn_run_dir("survives case", root=str(tmp_path))
    written = mcn_write_run(
        out_dir, "survives case", contrast(), rank(q=0.01, n_q_below=5),
        global_shift={"global_shift": False, "p": 0.9},
        summary="Marker clears every control.")
    text = open(written["readme"]).read()
    assert "**BRCA1 survives" in text
    assert "does not survive because" not in text


def test_run_readme_rejects_bulleted_summary():
    with pytest.raises(ValueError):
        mcn_run_readme("name", contrast(), rank(), summary="- not prose")


def test_run_readme_requires_summary():
    with pytest.raises(ValueError):
        mcn_run_readme("name", contrast(), rank(), summary=None)


def test_link_into_creates_relative_symlink_to_canonical_dir(tmp_path):
    canonical = tmp_path / "results" / "marker_contrast_null" / "usp1_slug"
    canonical.mkdir(parents=True)
    link_path = tmp_path / "results" / "target_triage" / "usp1" / "marker_null"
    mcn_link_into(str(canonical), str(link_path))
    assert os.path.islink(str(link_path))
    assert os.path.realpath(str(link_path)) == os.path.realpath(str(canonical))
    assert not os.path.isabs(os.readlink(str(link_path)))


def test_link_into_is_a_noop_when_already_linked(tmp_path):
    canonical = tmp_path / "canonical"
    canonical.mkdir()
    link_path = tmp_path / "caller" / "marker_null"
    mcn_link_into(str(canonical), str(link_path))
    target_before = os.readlink(str(link_path))
    mcn_link_into(str(canonical), str(link_path))
    assert os.readlink(str(link_path)) == target_before


def test_link_into_replaces_a_stale_link(tmp_path):
    old_canonical = tmp_path / "old"
    old_canonical.mkdir()
    new_canonical = tmp_path / "new"
    new_canonical.mkdir()
    link_path = tmp_path / "caller" / "marker_null"
    mcn_link_into(str(old_canonical), str(link_path))
    mcn_link_into(str(new_canonical), str(link_path))
    assert os.path.realpath(str(link_path)) == os.path.realpath(str(new_canonical))


def test_link_into_refuses_to_overwrite_a_real_directory(tmp_path):
    canonical = tmp_path / "canonical"
    canonical.mkdir()
    link_path = tmp_path / "caller" / "marker_null"
    link_path.mkdir(parents=True)
    with pytest.raises(FileExistsError):
        mcn_link_into(str(canonical), str(link_path))
