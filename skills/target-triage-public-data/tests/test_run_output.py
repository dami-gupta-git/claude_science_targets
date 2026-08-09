"""Coverage for the results/<topic>/<run>/ writer added to target-triage-public-data.

Uses a tmp_path root and a synthetic step rather than a live DepMap/Open
Targets fetch.
"""
import os

import pytest

from kernel import results_root, triage_run_dir, triage_write_run

STEPS = [{"name": "Dependency", "finding": "Cells tolerate loss of the gene."}]


def test_results_root_raises_when_unconfigured(monkeypatch):
    monkeypatch.delenv("SCIENCE_RESULTS_ROOT", raising=False)
    with pytest.raises(FileNotFoundError, match="SCIENCE_RESULTS_ROOT"):
        results_root()


def test_run_dir_raises_and_creates_nothing_when_root_unconfigured(monkeypatch, tmp_path):
    monkeypatch.delenv("SCIENCE_RESULTS_ROOT", raising=False)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError):
        triage_run_dir("WRN")
    assert not (tmp_path / "results").exists()


def test_run_dir_honours_env_var(monkeypatch, tmp_path):
    monkeypatch.setenv("SCIENCE_RESULTS_ROOT", str(tmp_path))
    out_dir = triage_run_dir("WRN")
    assert out_dir == str(tmp_path / "target_triage" / "wrn")


def test_run_dir_slugs_gene_and_makes_scripts(tmp_path):
    out_dir = triage_run_dir("WRN", root=str(tmp_path))
    assert out_dir == str(tmp_path / "target_triage" / "wrn")
    assert os.path.isdir(os.path.join(out_dir, "scripts"))


def test_run_dir_rejects_empty_gene(tmp_path):
    with pytest.raises(ValueError):
        triage_run_dir("   ", root=str(tmp_path))


def test_write_run_lands_readme_in_its_own_run_dir(tmp_path):
    out_dir = triage_run_dir("WRN", root=str(tmp_path))
    written = triage_write_run(
        out_dir, gene="WRN",
        summary="Plain prose a non-specialist can follow.",
        steps=STEPS,
        data_sources=["DepMap 24Q2 CRISPRGeneEffect.csv"],
        limits=["Knockout is not pharmacological inhibition."])
    assert written["readme"] == os.path.join(out_dir, "README.md")
    assert os.path.isfile(written["readme"])
    text = open(written["readme"]).read()
    assert "# WRN target triage" in text
    assert "### Dependency" in text


def test_write_run_accepts_a_script_already_in_place(tmp_path):
    """Writing wiring straight to scripts/ is the normal path, not an error.

    shutil.copyfile raises SameFileError on a same-path copy, which would fail
    the run after the analysis had already completed.
    """
    out_dir = triage_run_dir("WRN", root=str(tmp_path))
    script = os.path.join(out_dir, "scripts", "build.py")
    with open(script, "w") as fh:
        fh.write("# wiring\n")
    written = triage_write_run(
        out_dir, gene="WRN", summary="Summary.", steps=STEPS,
        scripts=[script])
    assert written["scripts"] == [script]
    assert open(script).read() == "# wiring\n"


def test_write_run_copies_scripts(tmp_path):
    script = tmp_path / "build.py"
    script.write_text("# wiring\n")
    out_dir = triage_run_dir("WRN", root=str(tmp_path / "results"))
    written = triage_write_run(
        out_dir, gene="WRN", summary="Summary.", steps=STEPS,
        scripts=[str(script)])
    assert os.path.isfile(written["scripts"][0])
    assert written["scripts"][0].endswith("scripts/build.py")
