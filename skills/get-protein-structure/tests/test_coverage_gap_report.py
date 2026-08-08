"""Regression test: authpos must be keyed by (chain, residue_num), not
residue_num alone, so a homodimer/heteromultimer's chains sharing residue
numbers don't overwrite each other's CA coordinates when no chain is given.
"""
import kernel
from helpers import protein_residue, hetero_residue, write_structure


def test_flank_distance_uses_correct_chain_not_last_written(monkeypatch, tmp_path):
    # Chain A: residues 1..20 near the origin. Chain B: same residue
    # numbers 1..20, but far away. The ligand sits close to chain A only.
    # Iteration order is A then B, so a residue-number-only dict would
    # have chain B silently overwrite chain A's coordinates for the same
    # numbers -- corrupting the flank-to-ligand distance.
    chain_a = [protein_residue("ALA", i, (i, 0, 0)) for i in range(1, 21)]
    chain_a.append(hetero_residue("LIG", 101, (9, 0, 5)))  # ~5A from residue 9
    chain_b = [protein_residue("ALA", i, (1000 + i, 0, 0)) for i in range(1, 21)]
    struct_path = write_structure({"A": chain_a, "B": chain_b}, str(tmp_path / "dimer.cif"))

    fake_cov = {
        "gaps": [{"kind": "internal", "unp_start": 110, "unp_end": 112,
                  "author_start": 10, "author_end": 12, "length": 3}],
    }
    monkeypatch.setattr(kernel, "residue_coverage",
                         lambda pdb_id, accession, chain=None: dict(fake_cov))

    out = kernel.coverage_gap_report("TEST", "ACC", chain=None,
                                      ligand_comp_id="LIG",
                                      struct_path=struct_path, plddt={})

    gap = out["gaps"][0]
    assert gap["near_site"] is True
    assert gap["flank_dist_to_ligand_A"] == 5.0
