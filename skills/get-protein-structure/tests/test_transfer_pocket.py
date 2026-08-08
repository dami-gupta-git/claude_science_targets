"""Regression tests for transfer_pocket's chain handling.

Covers two bugs found in review:
  1. an invalid target_chain/donor_chain silently substitutes another chain
     instead of erroring, unlike residue_coverage()'s explicit-error
     discipline elsewhere in this module
  3. the injected ligand is always written to a hardcoded chain "L", which
     silently merges into any structure that already has a chain named "L"

No network access: fetch_pdb is monkeypatched to return a synthetic donor
structure built with the same helpers used elsewhere in this test suite.
"""
import math
import os

import gemmi
import kernel
from helpers import protein_residue, hetero_residue, write_structure


def _helix_pos(i, radius=2.3, rise=1.5, twist_deg=100.0):
    # Alpha-helix-like CA path (~3.8A CA-CA spacing). A straight line of CA
    # atoms is numerically degenerate for Kabsch superposition (SVD is
    # singular for collinear/coplanar point sets), so tests need real 3D
    # geometry, not just distinct coordinates.
    theta = math.radians(twist_deg) * i
    return (radius * math.cos(theta), radius * math.sin(theta), rise * i)


def _build_donor(tmp_path):
    # Chain A: 20-residue helix; LIG sits ~3A from residue 10's CA.
    chain_a = [protein_residue("ALA", i, _helix_pos(i)) for i in range(1, 21)]
    x, y, z = _helix_pos(10)
    chain_a.append(hetero_residue("LIG", 101, (x + 3, y, z)))
    return write_structure({"A": chain_a}, str(tmp_path / "donor.cif"), name="DONOR")


def _build_target(tmp_path, with_chain_l=False):
    chains = {"X": [protein_residue("ALA", i, _helix_pos(i)) for i in range(1, 21)]}
    if with_chain_l:
        # A pre-existing chain literally named "L", far from everything else,
        # so it can't be confused for the transferred ligand by position.
        chains["L"] = [protein_residue("GLY", i, (1000 + i, 0, 0)) for i in range(1, 4)]
    return write_structure(chains, str(tmp_path / "target.cif"), name="TARGET")


def _patch_fetch_pdb(monkeypatch, donor_path):
    monkeypatch.setattr(kernel, "fetch_pdb",
                         lambda pdb_id, out_dir="structures", fmt="cif": donor_path)


def test_invalid_target_chain_errors_instead_of_falling_back(monkeypatch, tmp_path):
    donor_path = _build_donor(tmp_path)
    target_path = _build_target(tmp_path)
    _patch_fetch_pdb(monkeypatch, donor_path)

    res = kernel.transfer_pocket(target_path, donor_pdb_id="TEST", ligand_comp_id="LIG",
                                 target_chain="NOSUCHCHAIN", out_dir=str(tmp_path / "out"),
                                 write_complex=False)

    assert "error" in res
    assert res.get("target_chain") != "X"  # must not have silently used chain X instead


def test_invalid_donor_chain_errors_instead_of_falling_back(monkeypatch, tmp_path):
    donor_path = _build_donor(tmp_path)
    target_path = _build_target(tmp_path)
    _patch_fetch_pdb(monkeypatch, donor_path)

    res = kernel.transfer_pocket(target_path, donor_pdb_id="TEST", ligand_comp_id="LIG",
                                 donor_chain="NOSUCHCHAIN", out_dir=str(tmp_path / "out"),
                                 write_complex=False)

    assert "error" in res
    assert res.get("donor_chain") != "A"  # must not have silently used chain A instead


def test_valid_transfer_succeeds_and_finds_pocket(monkeypatch, tmp_path):
    # Sanity check that the happy path still works once chain validation is added.
    donor_path = _build_donor(tmp_path)
    target_path = _build_target(tmp_path)
    _patch_fetch_pdb(monkeypatch, donor_path)

    res = kernel.transfer_pocket(target_path, donor_pdb_id="TEST", ligand_comp_id="LIG",
                                 target_chain="X", donor_chain="A",
                                 out_dir=str(tmp_path / "out"), write_complex=False)

    assert "error" not in res
    assert res["donor_chain"] == "A"
    assert res["target_chain"] == "X"
    assert res["rmsd"] < 0.5  # identical geometry -> near-perfect superposition
    assert any(p["seq_id"] == 10 for p in res["pocket"])


def test_written_complex_does_not_merge_ligand_into_existing_chain_L(monkeypatch, tmp_path):
    donor_path = _build_donor(tmp_path)
    target_path = _build_target(tmp_path, with_chain_l=True)
    _patch_fetch_pdb(monkeypatch, donor_path)

    res = kernel.transfer_pocket(target_path, donor_pdb_id="TEST", ligand_comp_id="LIG",
                                 target_chain="X", donor_chain="A",
                                 out_dir=str(tmp_path / "out"), write_complex=True)

    assert "error" not in res
    complex_path = res["complex_path"]
    st = gemmi.read_structure(complex_path)
    st.setup_entities()

    l_chain = next(c for c in st[0] if c.name == "L")
    l_residue_names = sorted(r.name for r in l_chain)
    # The pre-existing chain L's own residues (GLY x3) must survive untouched -
    # if the ligand got merged into this chain, LIG would show up here too.
    assert l_residue_names == ["GLY", "GLY", "GLY"]

    # The transferred ligand must still be present somewhere in the file.
    assert any(r.name == "LIG" for ch in st[0] for r in ch)
