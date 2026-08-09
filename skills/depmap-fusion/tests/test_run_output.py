"""Coverage for the results/<topic>/<run>/ writer added to depmap-fusion.

Uses a tmp_path root rather than the real results/ tree, and synthetic
fuse_target_row()-shaped rows rather than a live evidence/dependency join.
"""
import os

import pytest

from kernel import (fusion_link_into, fusion_run_dir, fusion_run_readme,
                    fusion_verdict_mix, fusion_write_run, results_root)


def test_results_root_raises_when_unconfigured(monkeypatch, tmp_path):
    monkeypatch.delenv("SCIENCE_RESULTS_ROOT", raising=False)
    with pytest.raises(FileNotFoundError, match="SCIENCE_RESULTS_ROOT"):
        results_root()


def test_run_dir_raises_and_creates_nothing_when_root_unconfigured(monkeypatch, tmp_path):
    monkeypatch.delenv("SCIENCE_RESULTS_ROOT", raising=False)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError):
        fusion_run_dir("TP53")
    # The bug this guards against: os.makedirs(exist_ok=True) running before
    # the root was validated, silently creating a "results/" directory
    # wherever the session's cwd happened to be.
    assert not (tmp_path / "results").exists()


def test_run_dir_honours_env_var(monkeypatch, tmp_path):
    monkeypatch.setenv("SCIENCE_RESULTS_ROOT", str(tmp_path))
    out_dir = fusion_run_dir("TP53")
    assert out_dir == str(tmp_path / "depmap_fusion" / "tp53")


def row(gene, verdict, knockout_actionable=True, ot_score=0.8,
        mean_effect=-0.9, frac_dependent=0.8):
    return {"gene": gene, "verdict": verdict,
            "knockout_actionable": knockout_actionable, "ot_score": ot_score,
            "depmap_class": "selective", "mean_effect": mean_effect,
            "frac_dependent": frac_dependent}


def test_run_dir_slugs_subject_and_makes_scripts(tmp_path):
    out_dir = fusion_run_dir("Lung Adenocarcinoma", root=str(tmp_path))
    assert out_dir == str(tmp_path / "depmap_fusion" / "lung_adenocarcinoma")
    assert os.path.isdir(os.path.join(out_dir, "scripts"))


def test_run_dir_rejects_empty_subject(tmp_path):
    with pytest.raises(ValueError):
        fusion_run_dir("   ", root=str(tmp_path))


def test_verdict_mix_counts_and_rejects_missing_verdict():
    rows = [row("EGFR", "concordant-dependency"), row("KRAS", "concordant-dependency"),
            row("TP53", "growth-suppressive-mismatch")]
    assert fusion_verdict_mix(rows) == {"concordant-dependency": 2,
                                        "growth-suppressive-mismatch": 1}
    with pytest.raises(ValueError):
        fusion_verdict_mix([{"gene": "X"}])


def test_write_run_dossier_writes_readme_and_no_table(tmp_path):
    out_dir = fusion_run_dir("TP53", root=str(tmp_path))
    written = fusion_write_run(
        out_dir, "TP53", [row("TP53", "growth-suppressive-mismatch",
                              knockout_actionable=False, mean_effect=0.42)],
        summary="TP53 is a growth-suppressive mismatch: knockout helps growth.",
        data_sources=["Open Targets Platform", "DepMap 24Q2"])
    assert os.path.isfile(written["readme"])
    text = open(written["readme"]).read()
    assert "growth-suppressive-mismatch" in text
    assert "## Limits" in text and "## Files" in text and "## Data sources" in text


def test_write_run_triage_writes_table_and_verdict_mix(tmp_path):
    out_dir = fusion_run_dir("lung adenocarcinoma", root=str(tmp_path))
    rows = [row("EGFR", "concordant-dependency"), row("KRAS", "concordant-dependency"),
            row("ALK", "inert-in-panel", knockout_actionable=False)]
    written = fusion_write_run(
        out_dir, "lung adenocarcinoma", rows,
        summary="Top OT targets for lung adenocarcinoma joined against 126 lines.",
        disease="MONDO_0005061", data_sources=["Open Targets Platform", "DepMap 24Q2"])
    assert os.path.isfile(written["table"])
    assert os.path.isfile(written["readme"])
    text = open(written["readme"]).read()
    assert "Joined n = 3" in text
    assert "concordant-dependency" in text and "inert-in-panel" in text


def test_write_run_rejects_row_without_verdict(tmp_path):
    out_dir = fusion_run_dir("BAD", root=str(tmp_path))
    with pytest.raises(ValueError):
        fusion_write_run(out_dir, "BAD", [{"gene": "BAD"}], summary="no verdict")


def test_run_readme_rejects_bulleted_summary():
    with pytest.raises(ValueError):
        fusion_run_readme("TP53", [row("TP53", "growth-suppressive-mismatch")],
                          summary="- not prose")


def test_run_readme_enforces_summary_word_cap():
    with pytest.raises(ValueError):
        fusion_run_readme("TP53", [row("TP53", "growth-suppressive-mismatch")],
                          summary=" ".join(["word"] * 200))


def test_link_into_creates_relative_symlink_to_canonical_dir(tmp_path):
    canonical = tmp_path / "results" / "depmap_fusion" / "egfr"
    canonical.mkdir(parents=True)
    link_path = tmp_path / "results" / "target_triage" / "egfr" / "fusion"
    fusion_link_into(str(canonical), str(link_path))
    assert os.path.islink(str(link_path))
    assert os.path.realpath(str(link_path)) == os.path.realpath(str(canonical))
    assert not os.path.isabs(os.readlink(str(link_path)))


def test_link_into_is_a_noop_when_already_linked(tmp_path):
    canonical = tmp_path / "canonical"
    canonical.mkdir()
    link_path = tmp_path / "caller" / "fusion"
    fusion_link_into(str(canonical), str(link_path))
    target_before = os.readlink(str(link_path))
    fusion_link_into(str(canonical), str(link_path))
    assert os.readlink(str(link_path)) == target_before


def test_link_into_replaces_a_stale_link(tmp_path):
    old_canonical = tmp_path / "old"
    old_canonical.mkdir()
    new_canonical = tmp_path / "new"
    new_canonical.mkdir()
    link_path = tmp_path / "caller" / "fusion"
    fusion_link_into(str(old_canonical), str(link_path))
    fusion_link_into(str(new_canonical), str(link_path))
    assert os.path.realpath(str(link_path)) == os.path.realpath(str(new_canonical))


def test_link_into_refuses_to_overwrite_a_real_directory(tmp_path):
    canonical = tmp_path / "canonical"
    canonical.mkdir()
    link_path = tmp_path / "caller" / "fusion"
    link_path.mkdir(parents=True)
    with pytest.raises(FileExistsError):
        fusion_link_into(str(canonical), str(link_path))
