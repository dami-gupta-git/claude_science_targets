"""Shared structure-building helpers for kernel.py tests. No network access."""
import gemmi


def make_atom(name, element, pos):
    a = gemmi.Atom()
    a.name = name
    a.element = gemmi.Element(element)
    a.pos = gemmi.Position(*pos)
    a.occ = 1.0
    a.b_iso = 20.0
    return a


def make_residue(name, seqnum, het_flag, atoms):
    r = gemmi.Residue()
    r.name = name
    r.seqid = gemmi.SeqId(seqnum, " ")
    r.het_flag = het_flag
    for a in atoms:
        r.add_atom(a)
    return r


def protein_residue(resname, seqnum, ca_pos):
    """A minimal amino-acid residue with just a CA, N and C atom."""
    x, y, z = ca_pos
    return make_residue(resname, seqnum, "A", [
        make_atom("CA", "C", (x, y, z)),
        make_atom("N", "N", (x - 1, y, z)),
        make_atom("C", "C", (x + 1, y, z)),
    ])


def hetero_residue(resname, seqnum, pos):
    return make_residue(resname, seqnum, "H", [make_atom(resname[:2].strip(), resname[:2].strip(), pos)])


def write_structure(chains, out_path, name="TEST"):
    """chains: dict of chain_id -> list of gemmi.Residue."""
    st = gemmi.Structure()
    st.name = name
    model = gemmi.Model("1")
    for chain_id, residues in chains.items():
        chain = gemmi.Chain(chain_id)
        for res in residues:
            chain.add_residue(res)
        model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()
    st.make_mmcif_document().write_file(out_path)
    return out_path
