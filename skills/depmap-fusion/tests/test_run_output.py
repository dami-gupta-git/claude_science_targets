"""Coverage for the results/<topic>/<run>/ writer added to depmap-fusion.

Uses a tmp_path root rather than the real results/ tree, and synthetic
fuse_target_row()-shaped rows rather than a live evidence/dependency join.
"""
import os

import pytest

from kernel import (fusion_run_dir, fusion_run_readme, fusion_verdict_mix,
                    fusion_write_run)


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
