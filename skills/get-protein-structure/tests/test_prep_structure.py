"""Regression test: keep_ions must win over keep_ligands=False.

SKILL.md and prep_structure's own docstring promise structural ions are
kept "regardless" of keep_ligands. A prior version dropped every non-water
heteroatom unconditionally when keep_ligands=False, stripping ions too.
"""
import os
import kernel
from helpers import protein_residue, hetero_residue, write_structure


def _build(tmp_path):
    chains = {
        "A": [
            protein_residue("ALA", 1, (0, 0, 0)),
            protein_residue("GLY", 2, (1, 0, 0)),
            hetero_residue("ZN", 101, (5, 5, 5)),
            hetero_residue("SO4", 102, (6, 6, 6)),
        ],
    }
    return write_structure(chains, str(tmp_path / "in.cif"))


def test_keep_ions_survives_keep_ligands_false(tmp_path):
    in_path = _build(tmp_path)
    out_path = str(tmp_path / "out.cif")

    rep = kernel.prep_structure(in_path, out_path, keep_ligands=False, keep_ions=True)

    kept_ids = {h["comp_id"] for h in rep["hetero_kept"]}
    assert "ZN" in kept_ids
    assert "SO4" not in kept_ids


def test_keep_ions_false_and_keep_ligands_false_drops_ion_too(tmp_path):
    in_path = _build(tmp_path)
    out_path = str(tmp_path / "out2.cif")

    rep = kernel.prep_structure(in_path, out_path, keep_ligands=False, keep_ions=False)

    kept_ids = {h["comp_id"] for h in rep["hetero_kept"]}
    assert "ZN" not in kept_ids
    assert "SO4" not in kept_ids
