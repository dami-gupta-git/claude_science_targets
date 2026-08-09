"""Coverage for the results/<topic>/<run>/ writer added to get-protein-structure.

Uses synthetic get_structure()-shaped result dicts and plain placeholder files
rather than a live RCSB/AlphaFold fetch, so this needs no network access.
"""
import os

import pytest

from kernel import gps_run_dir, gps_run_readme, gps_write_run


def make_result(tmp_path, pocket=True, coverage_gaps=True):
    tmp_path.mkdir(parents=True, exist_ok=True)
    raw = tmp_path / "9IAY.cif"
    raw.write_text("data_fake\n")
    prep = tmp_path / "9IAY_prep.cif"
    prep.write_text("data_fake_prep\n")
    result = {
        "uniprot": {"accession": "P01116", "gene": "KRAS"},
        "candidates": [{"pdb_id": "9IAY", "method": "X-ray diffraction",
                        "resolution": 1.5, "r_free": 0.19, "state": "holo",
                        "coverage": 0.95, "n_mutations": 0, "score": 88.2}],
        "chosen": "9IAY", "raw_path": str(raw), "path": str(prep),
    }
    if pocket:
        result["pocket"] = [{"chain": "A", "seq_id": 12, "comp_id": "GTP",
                             "min_dist": 3.1}]
        result["pocket_ligand"] = "GTP"
        result["pocket_chain"] = "A"
        result["box"] = {"center": [1, 2, 3], "size": [20, 20, 20]}
    if coverage_gaps:
        result["coverage_gaps"] = {"gaps": [], "n_missing": 0}
    return result


def test_run_dir_slugs_target_and_makes_scripts(tmp_path):
    out_dir = gps_run_dir("KRAS", root=str(tmp_path))
    assert out_dir == str(tmp_path / "protein_structure" / "kras")
    assert os.path.isdir(os.path.join(out_dir, "scripts"))


def test_run_dir_rejects_empty_target(tmp_path):
    with pytest.raises(ValueError):
        gps_run_dir("   ", root=str(tmp_path))


def test_write_run_copies_structure_files_and_writes_tables(tmp_path):
    result = make_result(tmp_path / "src")
    out_dir = gps_run_dir("KRAS", root=str(tmp_path / "results"))
    written = gps_write_run(out_dir, "KRAS", result,
                            summary="KRAS 9IAY: 1.5 A holo structure bound to GTP.",
                            data_sources=["RCSB PDB"])
    assert os.path.isfile(written["path"])
    assert os.path.isfile(written["raw_path"])
    assert os.path.isfile(written["candidates"])
    assert os.path.isfile(written["pocket"])
    assert "coverage_gaps" not in written  # empty gaps list writes nothing
    assert os.path.isfile(written["readme"])
    text = open(written["readme"]).read()
    assert "9IAY" in text and "## Limits" in text


def test_write_run_rejects_error_result(tmp_path):
    out_dir = gps_run_dir("BAD", root=str(tmp_path))
    with pytest.raises(ValueError):
        gps_write_run(out_dir, "BAD", {"error": "no UniProt match"}, summary="x")


def test_write_run_requires_summary(tmp_path):
    result = make_result(tmp_path / "src")
    out_dir = gps_run_dir("KRAS", root=str(tmp_path / "results"))
    with pytest.raises(ValueError):
        gps_write_run(out_dir, "KRAS", result)


def test_run_readme_uses_structure_caveats_as_limits(tmp_path):
    result = make_result(tmp_path / "src")
    result["candidates"][0]["state"] = "apo"
    del result["pocket"]
    del result["pocket_ligand"]
    del result["pocket_chain"]
    del result["box"]
    text = gps_run_readme("KRAS", result, summary="Apo KRAS structure, no bound ligand.")
    assert "APO" in text or "apo" in text.lower()


def test_run_readme_rejects_bulleted_summary(tmp_path):
    result = make_result(tmp_path / "src")
    with pytest.raises(ValueError):
        gps_run_readme("KRAS", result, summary="- not prose")
