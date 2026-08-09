"""Helpers for get-protein-structure: resolve identifier -> rank -> fetch -> prep."""
import csv
import json
import re
import shutil
import urllib.request
import urllib.parse
import gzip
import os

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
SIFTS_BEST = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/"
RCSB_GQL = "https://data.rcsb.org/graphql"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/"
RCSB_FILES = "https://files.rcsb.org/download/"

WATERS = ("HOH", "DOD", "WAT")
STRUCTURAL_IONS = ("ZN", "MG", "MN", "CA", "FE", "FE2", "CU", "CU1", "NI", "CO", "K", "NA")
ADDITIVES = (
    "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "FE", "FE2", "CU", "CU1", "NI", "CO",
    "CD", "HG", "AU", "AG", "PT", "BR", "IOD", "F", "CS", "RB", "SR", "BA", "LI", "AL", "TL",
    "GOL", "EDO", "PEG", "PGE", "PG4", "P6G", "1PE", "2PE", "MPD", "DMS", "DMF", "TRS", "MES",
    "EPE", "IMD", "ACT", "ACY", "ACE", "FMT", "OXL", "CIT", "FLC", "TAR", "MLA", "SCN", "AZI",
    "NO3", "CO3", "BCT", "NH4", "URE", "BME", "MRD", "BU3", "IPA", "SIN", "MOH", "EOH", "TFA",
    "PYR", "GLC", "NAG", "BMA", "MAN", "FUC", "GAL", "XYL", "SIA", "BEZ", "SGM", "LDA", "BOG",
    "C8E", "OCT", "HEZ", "D10", "D12", "UNX", "UNL", "UNK", "DTT", "TCE", "MYR", "PLM", "OLA",
    "STE", "EDT", "POP", "PPV", "VO4", "WO4", "MOO", "CAC", "BEF", "ALF", "SPD", "SPM", "PUT",
)
COFACTORS = (
    "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "ANP", "AGS", "GNP", "GSP", "GCP", "ACP", "APC",
    "NAD", "NAI", "NAP", "NDP", "FAD", "FMN", "SAM", "SAH", "COA", "ACO", "HEM", "HEC", "HEA",
    "SF4", "FES", "F3S", "CLA", "BCL", "PLP", "TPP", "B12", "BTN", "THF", "UDP", "UTP", "CTP",
    "TTP", "5GP", "IHP", "NGD", "GSH", "PPS",
)


def http_get(url, params=None, timeout=60, binary=False, retries=3):
    """GET with retry/backoff. EBI and RCSB both time out intermittently;
    a 404 is a real answer and is re-raised immediately."""
    import time
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "*/*"})
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                body = r.read()
            if url.endswith(".gz"):
                body = gzip.decompress(body)
            return body if binary else json.loads(body)
        except urllib.error.HTTPError as e:
            if e.code in (404, 422, 400):
                raise
            last = e
        except Exception as e:
            last = e
        if attempt < retries - 1:
            time.sleep(2 ** attempt)
    raise last


def classify_ligand(comp_id, mw=None):
    """Bucket a chem-comp id: water | additive | cofactor | drug_like.

    Heuristic — RCSB exposes no curated 'subject of investigation' flag.
    Used for RANKING (is this entry holo?), never for deciding what to keep
    during prep: structural ions land in 'additive' but must be retained.
    """
    cid = (comp_id or "").upper()
    if cid in WATERS:
        return "water"
    if cid in COFACTORS:
        return "cofactor"
    if cid in ADDITIVES:
        return "additive"
    if mw is not None and mw < 120:
        return "additive"
    return "drug_like"


def resolve_uniprot(query, organism_id=9606, reviewed=True):
    """Gene symbol / protein name / accession -> UniProt record(s).

    Returns list of dicts (accession, uniprot_id, name, gene, length).
    An input that already looks like an accession is looked up directly.
    """
    import re
    q = query.strip()
    # Real UniProt accession grammar (avoids misreading gene symbols like
    # TMEM238 or SLC25A6 as accessions, which makes the API 400).
    acc_re = r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}$"
    if re.match(acc_re, q.upper()):
        terms = ["accession:" + q.upper()]
    else:
        terms = ["(gene_exact:%s OR protein_name:%s)" % (q, q)]
        if organism_id:
            terms.append("organism_id:%d" % organism_id)
    if reviewed:
        terms.append("reviewed:true")
    d = http_get(UNIPROT_SEARCH, {
        "query": " AND ".join(terms), "format": "json", "size": 10,
        "fields": "accession,id,protein_name,length,gene_primary,organism_name"})
    out = []
    for r in d.get("results", []):
        desc = r.get("proteinDescription", {}).get("recommendedName", {})
        genes = r.get("genes") or [{}]
        out.append({
            "accession": r["primaryAccession"],
            "uniprot_id": r.get("uniProtkbId"),
            "name": desc.get("fullName", {}).get("value"),
            "gene": (genes[0].get("geneName") or {}).get("value"),
            "organism": (r.get("organism") or {}).get("scientificName"),
            "length": (r.get("sequence") or {}).get("length"),
        })
    return out


def sifts_chains(accession):
    """SIFTS best_structures: every PDB chain mapped to this accession.

    Each row: pdb_id, chain_id, experimental_method, resolution, coverage
    (fraction of the UniProt sequence observed), unp_start, unp_end.
    """
    try:
        d = http_get(SIFTS_BEST + accession.upper())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return []
        raise
    return d.get(accession.upper(), [])


def rcsb_entries(pdb_ids):
    """Batch entry metadata from the RCSB GraphQL API (ligands, R-free, mutations).

    One request for the whole list — do not loop per id.
    """
    q = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id
      struct{title}
      rcsb_entry_info{resolution_combined experimental_method
        deposited_polymer_entity_instance_count deposited_nonpolymer_entity_instance_count}
      refine{ls_R_factor_R_free}
      rcsb_accession_info{initial_release_date}
      nonpolymer_entities{
        rcsb_nonpolymer_entity_container_identifiers{auth_asym_ids}
        nonpolymer_comp{chem_comp{id name formula formula_weight}
          rcsb_chem_comp_descriptor{SMILES_stereo InChIKey}}}
      polymer_entities{
        rcsb_polymer_entity_container_identifiers{auth_asym_ids}
        entity_poly{rcsb_sample_sequence_length}
        rcsb_entity_source_organism{ncbi_scientific_name}
        rcsb_polymer_entity{pdbx_mutation}
        rcsb_polymer_entity_align{reference_database_accession}}}}"""
    out = {}
    ids = [p.upper() for p in pdb_ids]
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        body = json.dumps({"query": q, "variables": {"ids": chunk}}).encode()
        req = urllib.request.Request(RCSB_GQL, data=body, headers={
            "Content-Type": "application/json", "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=120) as r:
            res = json.loads(r.read())
        if "errors" in res and not res.get("data"):
            raise RuntimeError("RCSB GraphQL: " + res["errors"][0]["message"])
        for e in (res.get("data", {}).get("entries") or []):
            out[e["rcsb_id"]] = e
    return out


def summarize_entry(entry, accession=None):
    """Flatten one rcsb_entries record into a ranking-ready dict."""
    info = entry.get("rcsb_entry_info") or {}
    reso = (info.get("resolution_combined") or [None])[0]
    refine = entry.get("refine") or [{}]
    ligs = []
    for n in (entry.get("nonpolymer_entities") or []):
        comp = (n.get("nonpolymer_comp") or {})
        cc = comp.get("chem_comp") or {}
        desc = comp.get("rcsb_chem_comp_descriptor") or {}
        mw = cc.get("formula_weight")
        ligs.append({
            "comp_id": cc.get("id"), "name": cc.get("name"), "formula": cc.get("formula"),
            "mw": mw, "smiles": desc.get("SMILES_stereo"), "inchikey": desc.get("InChIKey"),
            "chains": (n.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get("auth_asym_ids"),
            "role": classify_ligand(cc.get("id"), mw),
        })
    muts, chains, organisms = [], [], []
    for p in (entry.get("polymer_entities") or []):
        aligns = p.get("rcsb_polymer_entity_align") or []
        accs = [a.get("reference_database_accession") for a in aligns]
        if accession and accession.upper() not in [str(a).upper() for a in accs]:
            continue
        m = (p.get("rcsb_polymer_entity") or {}).get("pdbx_mutation")
        if m:
            muts.append(m)
        chains += (p.get("rcsb_polymer_entity_container_identifiers") or {}).get("auth_asym_ids") or []
        organisms += [o.get("ncbi_scientific_name") for o in (p.get("rcsb_entity_source_organism") or [])]
    drug = [l for l in ligs if l["role"] == "drug_like"]
    return {
        "pdb_id": entry["rcsb_id"],
        "title": (entry.get("struct") or {}).get("title"),
        "method": info.get("experimental_method"),
        "resolution": reso,
        "r_free": refine[0].get("ls_R_factor_R_free") if refine else None,
        "released": (entry.get("rcsb_accession_info") or {}).get("initial_release_date", "")[:10],
        "n_protein_chains": info.get("deposited_polymer_entity_instance_count"),
        "target_chains": sorted(set(chains)),
        "organisms": sorted({o for o in organisms if o}),
        "mutations": "; ".join(muts) or None,
        "n_mutations": sum(len(m.split(",")) for m in muts),
        "state": "holo" if drug else "apo",
        "drug_like_ligands": [l["comp_id"] for l in drug],
        "cofactors": [l["comp_id"] for l in ligs if l["role"] == "cofactor"],
        "ions_additives": [l["comp_id"] for l in ligs if l["role"] == "additive"],
        "ligands": ligs,
    }


METHOD_ALIASES = (
    ("xray", ("x-ray", "x-ray diffraction", "xray", "xrd", "x ray")),
    ("em", ("em", "electron microscopy", "cryo-em", "cryoem",
            "electron crystallography", "single-particle em")),
    ("nmr", ("nmr", "solution nmr", "solid-state nmr")),
    ("neutron", ("neutron", "neutron diffraction")),
)


def method_family(value):
    """Normalise a method string to a family key: xray | em | nmr | neutron.

    Two APIs describe the same experiment differently — SIFTS returns
    'X-ray diffraction' / 'Electron Microscopy', RCSB returns 'X-ray' / 'EM' —
    so comparing raw strings across them silently fails to match.
    """
    v = (value or "").strip().lower()
    if not v:
        return None
    for key, forms in METHOD_ALIASES:
        if v in forms:
            return key
    for key, forms in METHOD_ALIASES:
        if any(f in v or v in f for f in forms):
            return key
    return v


def method_matches(value, wanted):
    """True when `value` and any entry of `wanted` denote the same method.

    Accepts either vocabulary in either position, so methods=["X-ray"] and
    methods=["X-ray diffraction"] both select X-ray entries.
    """
    fam = method_family(value)
    if fam is None:
        return False
    return any(method_family(w) == fam for w in (wanted or []))


def rank_structures(accession, prefer_ligand=True, prefer_wildtype=True,
                    min_coverage=0.0, methods=None, max_entries=60):
    """Rank every PDB entry for a UniProt accession. Returns list of dicts, best first.

    Score: resolution + UniProt coverage + method, with bonuses for a
    drug-like ligand (prefer_ligand) and penalty for engineered mutations
    (prefer_wildtype). `coverage` is the fraction of the FULL UniProt
    sequence observed — for multi-domain targets the useful entries often
    cover ~0.3, so read the coverage column rather than trusting rank alone.
    """
    rows = sifts_chains(accession)
    if not rows:
        return []
    best = {}
    for r in rows:
        pid = r["pdb_id"].upper()
        if r.get("coverage", 0) < min_coverage:
            continue
        if methods and not method_matches(r.get("experimental_method"), methods):
            continue
        prev = best.get(pid)
        if prev is None or (r.get("coverage") or 0) > (prev.get("coverage") or 0):
            best[pid] = r
    # Shortlist on a cheap resolution+coverage blend, NOT coverage alone:
    # a multi-domain target (EGFR: 531 chains at coverage ~0.3) would
    # otherwise fill the shortlist with full-length low-resolution EM and
    # starve the high-resolution single-domain entries that carry the pocket.
    def _pre(r):
        # Null resolution: EM gets the 3.5 A stand-in, everything else (NMR,
        # which has no resolution at all) must not be treated as mid-range or it
        # crowds high-resolution X-ray entries out of the shortlist.
        reso = r.get("resolution")
        if reso is None:
            reso = 3.5 if method_family(r.get("experimental_method")) == "em" else 99.0
        return -(30.0 * (r.get("coverage") or 0) + 40.0 * max(0.0, 1.0 - min(reso, 5.0) / 5.0))
    order = sorted(best.values(), key=lambda r: (_pre(r), r["pdb_id"]))
    picked = [r["pdb_id"].upper() for r in order[:max_entries]]
    meta = rcsb_entries(picked)
    out = []
    for pid in picked:
        e = meta.get(pid)
        if not e:
            continue
        s = summarize_entry(e, accession)
        sif = best[pid]
        s["coverage"] = round(sif.get("coverage") or 0, 3)
        s["unp_range"] = [sif.get("unp_start"), sif.get("unp_end")]
        s["sifts_chain"] = sif.get("chain_id")
        # Compare method FAMILIES, not raw strings: RCSB reports 'X-ray'/'EM'
        # while SIFTS reports 'X-ray diffraction'/'Electron Microscopy', and a
        # raw comparison silently stops matching if either changes.
        fam = method_family(s["method"])
        reso = s["resolution"] or (3.5 if fam == "em" else 99)
        score = 0.0
        score += max(0.0, 40.0 * (1.0 - min(reso, 5.0) / 5.0))
        score += 30.0 * s["coverage"]
        if prefer_ligand and s["state"] == "holo":
            score += 20.0
        if prefer_wildtype:
            score -= 3.0 * min(s["n_mutations"], 4)
        if fam == "xray":
            score += 5.0
        if s["r_free"] and s["r_free"] < 0.25:
            score += 3.0
        s["score"] = round(score, 1)
        out.append(s)
    return sorted(out, key=lambda s: -s["score"])


def fetch_pdb(pdb_id, out_dir="structures", fmt="cif"):
    """Download coordinates from RCSB. fmt: 'cif' (default) or 'pdb'."""
    os.makedirs(out_dir, exist_ok=True)
    pid = pdb_id.upper()
    try:
        body = http_get(RCSB_FILES + pid + "." + fmt + ".gz", binary=True)
    except urllib.error.HTTPError as e:
        # RCSB serves NO legacy .pdb file for entries too large for the format
        # (big assemblies, ribosomes). Fall back to mmCIF rather than dead-end.
        if e.code == 404 and fmt == "pdb":
            body = http_get(RCSB_FILES + pid + ".cif.gz", binary=True)
            fmt = "cif"
        else:
            raise
    path = os.path.join(out_dir, pid + "." + fmt)
    with open(path, "wb") as fh:
        fh.write(body)
    return path


def fetch_alphafold(accession, out_dir="structures", fmt="cif", with_pae=False):
    """Download the AlphaFold DB model for a UniProt accession.

    Returns dict with path, global pLDDT, per-confidence-bin residue
    fractions, version, and the model's UniProt range. has_model=False when
    AFDB has no prediction (not an error).
    """
    os.makedirs(out_dir, exist_ok=True)
    try:
        recs = http_get(AFDB_API + accession.upper())
    except urllib.error.HTTPError as e:
        if e.code in (404, 422):
            return {"has_model": False, "accession": accession}
        raise
    if not recs:
        return {"has_model": False, "accession": accession}
    m = recs[0]
    url = m["cifUrl"] if fmt == "cif" else m["pdbUrl"]
    path = os.path.join(out_dir, os.path.basename(url))
    with open(path, "wb") as fh:
        fh.write(http_get(url, binary=True))
    res = {
        "has_model": True, "accession": accession, "path": path,
        "entry_id": m.get("entryId"), "version": m.get("latestVersion"),
        "global_plddt": m.get("globalMetricValue"),
        "plddt_fractions": {
            "very_high_gt90": m.get("fractionPlddtVeryHigh"),
            "confident_70_90": m.get("fractionPlddtConfident"),
            "low_50_70": m.get("fractionPlddtLow"),
            "very_low_lt50": m.get("fractionPlddtVeryLow")},
        "uniprot_range": [m.get("uniprotStart"), m.get("uniprotEnd")],
        "n_models": len(recs), "gene": m.get("gene"),
        "organism": m.get("organismScientificName"),
    }
    if with_pae and m.get("paeDocUrl"):
        p = os.path.join(out_dir, os.path.basename(m["paeDocUrl"]))
        with open(p, "wb") as fh:
            fh.write(http_get(m["paeDocUrl"], binary=True))
        res["pae_path"] = p
    return res


def prep_structure(in_path, out_path, chains=None, keep_waters=False,
                   keep_ligands=True, drop_additives=True, keep_ions=True,
                   keep_hydrogens=False, model_index=0):
    """Clean a structure for docking / design. Returns a report dict.

    Drops waters, alt-locs, extra models and (by default) crystallisation
    additives; KEEPS structural ions and cofactors. Selects `chains` by
    author chain id. Output format follows out_path's extension.
    """
    import gemmi
    d = os.path.dirname(out_path)
    if d:
        os.makedirs(d, exist_ok=True)
    st = gemmi.read_structure(in_path)
    st.setup_entities()
    n_models = len(st)
    if model_index < 0 or model_index >= n_models:
        # Silently returning a different model would misreport which
        # conformer the output represents.
        raise IndexError("model_index %d out of range: %s has %d model(s) (0-%d)"
                         % (model_index, os.path.basename(in_path), n_models, n_models - 1))
    st.remove_alternative_conformations()
    if not keep_hydrogens:
        st.remove_hydrogens()
    while len(st) > model_index + 1:
        del st[len(st) - 1]
    while len(st) > 1:
        del st[0]
    ions = set(STRUCTURAL_IONS)
    model, kept, dropped = st[0], [], []
    for ch in model:
        doomed = []
        for i, res in enumerate(ch):
            info = gemmi.find_tabulated_residue(res.name)
            if res.name in WATERS or (info and info.is_water()):
                if not keep_waters:
                    doomed.append(i)
                    dropped.append({"chain": ch.name, "comp_id": res.name, "why": "water"})
                continue
            if info and (info.is_amino_acid() or info.is_nucleic_acid()):
                continue
            role = classify_ligand(res.name)
            # keep_ions wins over keep_ligands=False: structural ions (catalytic
            # Zn, Mg in a nucleotide site) are part of the protein's functional
            # state, and both SKILL.md and this docstring promise they survive.
            if keep_ions and res.name in ions:
                kept.append({"chain": ch.name, "comp_id": res.name,
                             "seq_id": res.seqid.num, "role": role})
            elif not keep_ligands:
                doomed.append(i)
                dropped.append({"chain": ch.name, "comp_id": res.name, "why": "keep_ligands=False"})
            elif drop_additives and role == "additive" and not (keep_ions and res.name in ions):
                doomed.append(i)
                dropped.append({"chain": ch.name, "comp_id": res.name, "why": "additive"})
            else:
                kept.append({"chain": ch.name, "comp_id": res.name,
                             "seq_id": res.seqid.num, "role": role})
        for i in reversed(doomed):
            del ch[i]
    if chains:
        want = set(chains)
        missing = want - {c.name for c in model}
        if missing:
            raise ValueError("chain(s) %s not present in %s (chains: %s)"
                             % (", ".join(sorted(missing)), os.path.basename(in_path),
                                ", ".join(c.name for c in model)))
        for i in reversed(range(len(model))):
            if model[i].name not in want:
                del model[i]
        # The hetero inventory was built while walking every chain, so entries
        # from chains just removed would describe residues absent from the
        # output file. Move them to `dropped` with the reason.
        kept, moved = [h for h in kept if h["chain"] in want], [h for h in kept if h["chain"] not in want]
        dropped.extend({"chain": h["chain"], "comp_id": h["comp_id"], "why": "chain not selected"}
                       for h in moved)
    st.remove_empty_chains()
    st.setup_entities()
    # Legacy PDB format has a 3-char resName field; the PDB has issued
    # 5-char CCD ids since 2023 (e.g. A1BEA), which get silently truncated
    # and break every downstream comp_id lookup. Coerce to mmCIF instead.
    long_ids = sorted({h["comp_id"] for h in kept if len(h["comp_id"]) > 3})
    coerced = None
    reasons = []
    if out_path.endswith(".pdb"):
        if long_ids:
            reasons.append("ligand id(s) %s exceed the 3-char PDB resName field"
                           % ", ".join(long_ids))
        # PDB format ceilings: 1 chain character and a 5-digit atom serial.
        # gemmi does NOT raise past these — it silently merges chains and
        # overflows serials, so guard rather than trust the writer.
        n_chains = len(st[0])
        n_atoms = sum(len(r) for c in st[0] for r in c)
        if n_chains > 62:
            reasons.append("%d chains exceed the 62 single-character chain ids"
                           % n_chains)
        if n_atoms > 99999:
            reasons.append("%d atoms exceed the 99,999 serial-number limit" % n_atoms)
        if any(len(c.name) > 1 for c in st[0]):
            reasons.append("multi-character chain id(s) %s do not fit the PDB "
                           "1-char chainID field"
                           % ", ".join(sorted({c.name for c in st[0] if len(c.name) > 1})))
    if reasons:
        coerced = out_path
        out_path = out_path[:-4] + ".cif"
    if out_path.endswith(".pdb"):
        with open(out_path, "w") as fh:
            fh.write(st.make_pdb_string())
    else:
        st.make_mmcif_document().write_file(out_path)
    n_res = sum(1 for c in st[0] for r in c
                if gemmi.find_tabulated_residue(r.name)
                and gemmi.find_tabulated_residue(r.name).is_amino_acid())
    rep = {"path": out_path, "chains": [c.name for c in st[0]],
           "n_residues": n_res,
           "n_atoms": sum(len(r) for c in st[0] for r in c),
           "hetero_kept": kept, "n_dropped": len(dropped),
           "dropped_kinds": sorted({d["why"] for d in dropped})}
    if coerced:
        rep["format_coerced"] = (
            "wrote mmCIF (%s) instead of %s — legacy PDB format would silently "
            "corrupt this structure: %s" % (out_path, coerced, "; ".join(reasons)))
        rep["format_coerced_reasons"] = reasons
    return rep


def ligand_copies(path_or_st, ligand_comp_id, model_index=0):
    """Every copy of a ligand in the structure: [{chain, seq_id, n_atoms}].

    A ligand present in several chains of a multimer yields several copies —
    pocket_residues/ligand_box operate on ONE copy, so check this first when
    the entry has more than one protein chain.
    """
    import gemmi
    st = gemmi.read_structure(path_or_st) if isinstance(path_or_st, str) else path_or_st
    st.setup_entities()
    return [{"chain": ch.name, "seq_id": r.seqid.num, "n_atoms": len(r)}
            for ch in st[model_index] for r in ch
            if r.name.upper() == ligand_comp_id.upper()]


def pocket_residues(path_or_st, ligand_comp_id, radius=5.0, model_index=0,
                    chain=None, all_copies=False):
    """Protein residues within `radius` A of a bound ligand, nearest-atom distance.

    Returns list of {chain, seq_id, comp_id, min_dist} sorted by sequence.
    Defaults to the FIRST copy of the ligand — in a multimer the same ligand
    sits in every chain, and pooling them would report one merged pseudo-site.
    Pass chain= to pick a copy or all_copies=True to pool deliberately.
    """
    import gemmi
    st = gemmi.read_structure(path_or_st) if isinstance(path_or_st, str) else path_or_st
    st.setup_entities()
    model = st[model_index]
    ns = gemmi.NeighborSearch(model, st.cell, radius).populate()
    targets = [r for ch in model for r in ch if r.name.upper() == ligand_comp_id.upper()
               and (chain is None or ch.name == chain)]
    if not targets:
        raise ValueError("ligand %s not present in %s%s" % (
            ligand_comp_id, getattr(st, "name", "structure"),
            "" if chain is None else " chain %s" % chain))
    if not all_copies:
        targets = targets[:1]
    found = {}
    for res in targets:
        for atom in res:
            for mark in ns.find_atoms(atom.pos, "\0", radius=radius):
                cra = mark.to_cra(model)
                info = gemmi.find_tabulated_residue(cra.residue.name)
                if not (info and info.is_amino_acid()):
                    continue
                d = cra.atom.pos.dist(atom.pos)
                key = (cra.chain.name, cra.residue.seqid.num, cra.residue.name)
                if key not in found or d < found[key]:
                    found[key] = d
    return [{"chain": c, "seq_id": n, "comp_id": nm, "min_dist": round(d, 2)}
            for (c, n, nm), d in sorted(found.items(), key=lambda kv: kv[0][1])]


def ligand_box(path_or_st, ligand_comp_id, padding=4.0, model_index=0, chain=None):
    """Docking box (center + size, Angstrom) around ONE copy of a bound ligand.

    Uses the first copy unless chain= is given. Pooling copies across a
    multimer would produce a box spanning both sites — never what docking
    wants, so it is not offered.
    """
    import gemmi
    st = gemmi.read_structure(path_or_st) if isinstance(path_or_st, str) else path_or_st
    st.setup_entities()
    hits = [r for ch in st[model_index] for r in ch
            if r.name.upper() == ligand_comp_id.upper() and (chain is None or ch.name == chain)]
    if not hits:
        raise ValueError("ligand %s not present%s" % (
            ligand_comp_id, "" if chain is None else " in chain %s" % chain))
    pos = [a.pos for a in hits[0]]
    xs = [p.x for p in pos]
    ys = [p.y for p in pos]
    zs = [p.z for p in pos]
    return {"center": [round(sum(v) / len(v), 2) for v in (xs, ys, zs)],
            "size": [round(max(v) - min(v) + 2 * padding, 2) for v in (xs, ys, zs)],
            "n_ligand_atoms": len(pos)}


PDBE_GRAPH_MAP = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/uniprot/"
PDBE_RES_LIST = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/"


def residue_coverage(pdb_id, accession, chain=None):
    """Per-residue observed/missing map for one PDB chain, in UniProt numbering.

    Answers the question resolution and coverage numbers hide: WHICH residues
    does this entry actually observe? Uses PDBe residue_listing
    (observed_ratio per residue) plus the SIFTS segment offsets.

    Returns dict with:
      chain, entity_id, unp_range, author_range, seq_offset
      n_observed / n_partial / n_missing / n_listed / n_unmapped_residues
        (n_unmapped_residues counts residues the chain contains but this
        accession does not claim — a fusion partner such as 2RH1's lysozyme.
        n_missing excludes expression_tag gaps, since those are construct
        additions rather than missing protein)
      residues: [{unp_num, author_num, comp_id, observed_ratio, status}]
      gaps: [{kind, unp_start, unp_end, length, author_start, author_end}]
        kind is 'internal' (a loop absent from an otherwise-present region —
        the dangerous one), 'terminal' (truncated N/C terminus), or
        'expression_tag' (construct addition numbered < 1 in UniProt space,
        e.g. a His-tag — not missing protein, when a SIFTS segment itself
        extends below residue 1).
    """
    pid = pdb_id.lower()
    acc = accession.upper()
    mp = http_get(PDBE_GRAPH_MAP + pid)
    segs = (((mp.get(pid) or {}).get("UniProt") or {}).get(acc) or {}).get("mappings") or []
    if not segs:
        return {"error": "no SIFTS mapping of %s onto %s" % (acc, pid.upper())}
    if chain:
        want = [s for s in segs if s.get("chain_id") == chain]
        if not want:
            # Never silently fall back to another chain's mapping: the offsets
            # would be wrong and the caller asked for a specific chain.
            return {"error": "chain %s of %s has no SIFTS mapping onto %s "
                             "(mapped chains: %s)"
                             % (chain, pid.upper(), acc,
                                ", ".join(sorted({str(s.get("chain_id")) for s in segs})))}
        segs = want
    # One accession can map to a chain in SEVERAL segments (fusion constructs,
    # e.g. 2RH1 = beta2-adrenergic receptor with a lysozyme insert: 1-230 and
    # 264-365). Each segment has its own author->UniProt offset, so keep them
    # all and resolve per residue rather than assuming segs[0] covers the chain.
    ch = segs[0].get("chain_id")
    ent = segs[0].get("entity_id")
    segs = [s for s in segs if s.get("chain_id") == ch and s.get("entity_id") == ent]

    usable = []
    for s in segs:
        # Either boundary can be null when that boundary residue is disordered
        # (2RH1 has null at both outer ends). Derive the offset from whichever
        # end is present; skip the segment only if both are missing.
        off = None
        if s.get("start", {}).get("author_residue_number") is not None:
            off = s["unp_start"] - s["start"]["author_residue_number"]
        elif s.get("end", {}).get("author_residue_number") is not None:
            off = s["unp_end"] - s["end"]["author_residue_number"]
        if off is None:
            continue
        usable.append({"unp_start": s["unp_start"], "unp_end": s["unp_end"], "offset": off})
    if not usable:
        return {"error": "no SIFTS segment of %s on %s chain %s has a usable "
                         "author-residue boundary" % (acc, pid.upper(), ch)}
    usable.sort(key=lambda s: s["unp_start"])
    offset = usable[0]["offset"]

    def unp_for(auth):
        """Map an author residue number through the segment that contains it."""
        for s in usable:
            if s["unp_start"] <= auth + s["offset"] <= s["unp_end"]:
                return auth + s["offset"], True
        return auth + offset, False

    rl = http_get(PDBE_RES_LIST + pid)
    residues, mols = [], (rl.get(pid) or {}).get("molecules") or []
    for mol in mols:
        if mol.get("entity_id") != ent:
            continue
        for c in mol.get("chains") or []:
            if c.get("chain_id") != ch and c.get("struct_asym_id") != ch:
                continue
            for r in c.get("residues") or []:
                ratio = r.get("observed_ratio")
                auth = r.get("author_residue_number")
                if auth is None:
                    continue
                status = ("missing" if not ratio else
                          "observed" if ratio >= 1.0 else "partial")
                unp, in_seg = unp_for(auth)
                residues.append({"unp_num": unp, "author_num": auth,
                                 "comp_id": r.get("residue_name"),
                                 "observed_ratio": ratio, "status": status,
                                 "in_mapped_segment": in_seg})
            break
        break
    if not residues:
        return {"error": "no residue listing for %s entity %s" % (pid.upper(), ent)}
    residues.sort(key=lambda r: (r["author_num"], r["unp_num"]))
    # Residues the chain contains but that this accession does not claim (a
    # fusion partner such as 2RH1's lysozyme insert) are not gaps in the target
    # protein — exclude them from gap analysis and from the counts.
    target = [r for r in residues if r["in_mapped_segment"]]
    n_foreign = len(residues) - len(target)
    target.sort(key=lambda r: r["unp_num"])
    obs = [r for r in target if r["status"] != "missing"]
    first_obs = obs[0]["unp_num"] if obs else None
    last_obs = obs[-1]["unp_num"] if obs else None

    # Walk the mapped UniProt range so a residue absent from the listing counts
    # as missing, and never bridge a gap across two separate mapped segments.
    gaps = []
    by_unp = {r["unp_num"]: r for r in target}
    for s in usable:
        run = []
        for u in range(s["unp_start"], s["unp_end"] + 1):
            r = by_unp.get(u)
            if r is None or r["status"] == "missing":
                run.append(r or {"unp_num": u, "author_num": u - s["offset"]})
                continue
            if run:
                gaps.append(run)
                run = []
        if run:
            gaps.append(run)
    out_gaps = []
    for g in gaps:
        s, e = g[0], g[-1]
        internal = (first_obs is not None and s["unp_num"] > first_obs
                    and e["unp_num"] < last_obs)
        # Residues numbered < 1 in UniProt space are construct additions
        # (His-tags, cloning scars), not missing protein.
        kind = ("expression_tag" if e["unp_num"] < 1 else
                "internal" if internal else "terminal")
        out_gaps.append({"kind": kind,
                         "unp_start": s["unp_num"], "unp_end": e["unp_num"],
                         "length": len(g), "author_start": s["author_num"],
                         "author_end": e["author_num"]})
    auths = [r["author_num"] for r in residues]
    return {"pdb_id": pid.upper(), "accession": acc, "chain": ch, "entity_id": ent,
            "seq_offset": offset,
            "segments": usable,
            "unp_range": [usable[0]["unp_start"], usable[-1]["unp_end"]],
            "author_range": [min(auths), max(auths)],
            "n_observed": sum(1 for r in target if r["status"] == "observed"),
            "n_partial": sum(1 for r in target if r["status"] == "partial"),
            # expression_tag gaps are construct additions, not missing protein
            # (see the classification above) -- excluded here for the same
            # reason, or a tagged construct would inflate this count with
            # residues that were never part of the target sequence.
            "n_missing": sum(g["length"] for g in out_gaps if g["kind"] != "expression_tag"),
            "n_listed": len(target),
            "n_unmapped_residues": n_foreign,
            "residues": target, "gaps": out_gaps}


def plddt_by_residue(accession_or_path, out_dir="structures"):
    """Per-residue pLDDT from an AlphaFold model (B-factor column).

    Accepts a UniProt accession (downloads the model) or a path to an
    AF cif/pdb. Returns {unp_num: plddt}.
    """
    import gemmi
    import re
    p = accession_or_path
    if not os.path.exists(p):
        af = fetch_alphafold(p, out_dir=out_dir)
        if not af.get("has_model"):
            return {}
        p = af["path"]
    st = gemmi.read_structure(p)
    st.setup_entities()
    out = {}
    for ch in st[0]:
        for res in ch:
            b = [a.b_iso for a in res]
            if b:
                out[res.seqid.num] = round(sum(b) / len(b), 2)
    return out


def coverage_gap_report(pdb_id, accession, chain=None, ligand_comp_id=None,
                        struct_path=None, plddt=None, near_site_angstrom=8.0):
    """Coverage gaps annotated with AlphaFold confidence and pocket proximity.

    For each gap: mean AF pLDDT over the missing stretch (can a predicted
    model credibly fill it, or is it disordered there too?) and — when a
    ligand and coordinate file are supplied — the distance from the gap's
    flanking residues to the ligand. A short internal gap whose flanks sit
    near the ligand is a loop that probably lines the pocket, which is the
    classic reason docking into an otherwise-good structure fails.

    `verdict` per gap: 'af_can_fill' (pLDDT >= 70), 'disordered_in_af'
    (< 70, likely genuinely flexible), or 'unknown' (no AF model).
    """
    cov = residue_coverage(pdb_id, accession, chain=chain)
    if "error" in cov:
        return cov
    if plddt is None:
        plddt = plddt_by_residue(accession)
    lig_pos = None
    if ligand_comp_id and struct_path and os.path.exists(struct_path):
        import gemmi
        st = gemmi.read_structure(struct_path)
        st.setup_entities()
        hits = [r for c in st[0] for r in c
                if r.name.upper() == ligand_comp_id.upper()
                and (chain is None or c.name == chain)]
        if hits:
            lig_pos = [a.pos for a in hits[0]]
            # Key by (chain, residue) — keying on residue number alone lets one
            # chain of a multimer overwrite another's coordinates, silently
            # corrupting every gap-to-ligand distance. When no chain was
            # requested, use the chain the chosen ligand copy actually sits in.
            lig_chain = chain
            if lig_chain is None:
                for c in st[0]:
                    if any(r.name.upper() == ligand_comp_id.upper() for r in c):
                        lig_chain = c.name
                        break
            authpos = {}
            for c in st[0]:
                for r in c:
                    ca = r.find_atom("CA", "*")
                    if ca:
                        authpos[(c.name, r.seqid.num)] = ca.pos
    for g in cov["gaps"]:
        vals = [plddt[u] for u in range(g["unp_start"], g["unp_end"] + 1) if u in plddt]
        g["af_mean_plddt"] = round(sum(vals) / len(vals), 1) if vals else None
        g["verdict"] = ("unknown" if g["af_mean_plddt"] is None else
                        "af_can_fill" if g["af_mean_plddt"] >= 70 else "disordered_in_af")
        g["near_site"] = None
        if lig_pos:
            ds = []
            for a in (g["author_start"] - 1, g["author_end"] + 1):
                pos = authpos.get((lig_chain, a))
                if pos is not None:
                    ds.append(min(pos.dist(lp) for lp in lig_pos))
            if ds:
                g["flank_dist_to_ligand_A"] = round(min(ds), 2)
                g["near_site"] = min(ds) <= near_site_angstrom
    internal = [g for g in cov["gaps"] if g["kind"] == "internal"]
    cov["n_internal_gaps"] = len(internal)
    cov["n_gaps_near_site"] = sum(1 for g in internal if g.get("near_site"))
    cov["fillable_by_af"] = [g for g in internal if g["verdict"] == "af_can_fill"]
    return cov


def copy_res_probe(st, comp_id, chain_name, model=0):
    """Atom positions of one ligand copy — used to find its surrounding shell."""
    for ch in st[model]:
        if chain_name and ch.name != chain_name:
            continue
        for r in ch:
            if r.name.upper() == comp_id.upper():
                return [a.pos for a in r]
    for ch in st[model]:
        for r in ch:
            if r.name.upper() == comp_id.upper():
                return [a.pos for a in r]
    return []


def aa_residues(chain):
    """Amino-acid residues of a chain, in order."""
    import gemmi
    out = []
    for r in chain:
        info = gemmi.find_tabulated_residue(r.name)
        if info and info.is_amino_acid():
            out.append(r)
    return out


def residue_correspondence(target_chain, donor_chain):
    """Map donor residue -> target residue by SEQUENCE alignment.

    Author numbering is a deposition choice, not a shared coordinate system: a
    donor crystallised as a construct fragment and a target taken from a
    predicted full-length model can describe the same residue with different
    numbers. Selecting a residue window by number therefore silently picks the
    wrong stretch whenever the two disagree (verified: shifting a WRN model by
    +1000 turns a 2.05 A transfer into 15.5 A with no error).

    Aligning the one-letter sequences removes that assumption. Returns
    (pairs, identity) where pairs is [(donor_residue, target_residue), ...] over
    aligned, non-gap positions and identity is the fraction of matched columns.
    """
    import gemmi
    t_res, d_res = aa_residues(target_chain), aa_residues(donor_chain)
    if not t_res or not d_res:
        return [], 0.0
    t_seq = gemmi.one_letter_code([r.name for r in t_res])
    d_seq = gemmi.one_letter_code([r.name for r in d_res])
    aln = gemmi.align_string_sequences(list(d_seq.upper()), list(t_seq.upper()), [])
    pairs, di, ti, matched = [], 0, 0, 0
    for op in aln.cigar_str().replace("M", "M ").replace("I", "I ").replace("D", "D ").split():
        n, kind = int(op[:-1] or 1), op[-1]
        if kind == "M":
            for _ in range(n):
                if di < len(d_res) and ti < len(t_res):
                    pairs.append((d_res[di], t_res[ti]))
                    if d_res[di].name == t_res[ti].name:
                        matched += 1
                di += 1
                ti += 1
        elif kind == "I":
            di += n
        else:
            ti += n
    return pairs, (matched / len(pairs) if pairs else 0.0)


def span_residues(chain, keep_residues):
    """(structure, polymer) for `chain` restricted to the given Residue objects.

    Selection is by object identity, not residue number, so it stays correct
    when donor and target number the same residue differently. Order follows
    keep_residues, which must be the aligned order — superposition pairs the two
    spans positionally. Returns the owning Structure too: the caller must hold it
    while using the polymer, since gemmi's ResidueSpan does not keep its parent
    alive.
    """
    import gemmi
    st = gemmi.Structure()
    st.add_model(gemmi.Model("1"))
    c = gemmi.Chain(chain.name)
    for r in keep_residues:
        c.add_residue(r)
    st[0].add_chain(c)
    st.setup_entities()
    return st, st[0][0].get_polymer()


def transfer_pocket(target_path, donor_pdb_id=None, accession=None, ligand_comp_id=None,
                    target_chain=None, donor_chain=None, out_dir="structures",
                    radius=5.0, padding=4.0, plddt=None, write_complex=True):
    """Give a ligand-free model a docking site by borrowing one from a holo structure.

    A predicted (AlphaFold) or apo structure has no bound ligand, so there is no
    ligand to build a box around. This superposes a LIGANDED donor entry onto the
    target and carries the ligand across, yielding the same outputs the
    experimental path gives — lining residues, centre, box — plus the two numbers
    that say whether to trust them: the superposition RMSD and the model's pLDDT
    *at the site* (global pLDDT can be high while the pocket itself is poor).

    Pass `donor_pdb_id` to choose the donor, or just `accession` to pick the
    best-scoring holo entry automatically. Returns a dict with `donor`,
    `rmsd`, `n_aligned`, `ligand`, `pocket`, `box`, `site_plddt` and
    `confidence`; `error` when no liganded donor exists.

    The transferred site is a HYPOTHESIS positioned by homology — good geometry
    for a pocket that exists in the donor, silent about a pocket that does not.
    Always report `rmsd` and `site_plddt` alongside any result derived from it.
    """
    import gemmi
    if donor_pdb_id is None:
        if not accession:
            return {"error": "supply donor_pdb_id or accession"}
        ranked = rank_structures(accession, prefer_ligand=True)
        holo = [c for c in ranked if c["drug_like_ligands"]]
        if not holo:
            return {"error": "no liganded (holo) entry exists for %s — nothing to "
                             "transfer; cavity detection or de-novo pocket "
                             "prediction would be needed" % accession}
        donor_pdb_id = holo[0]["pdb_id"]
        ligand_comp_id = ligand_comp_id or holo[0]["drug_like_ligands"][0]
    donor_raw = fetch_pdb(donor_pdb_id, out_dir=out_dir, fmt="cif")
    st_d = gemmi.read_structure(donor_raw)
    st_d.setup_entities()
    st_d.remove_alternative_conformations()
    st_d.remove_hydrogens()
    while len(st_d) > 1:
        del st_d[len(st_d) - 1]

    if ligand_comp_id is None:
        cands = {}
        for ch in st_d[0]:
            for r in ch:
                info = gemmi.find_tabulated_residue(r.name)
                if info and (info.is_amino_acid() or info.is_nucleic_acid() or info.is_water()):
                    continue
                if classify_ligand(r.name) == "drug_like":
                    cands.setdefault(r.name, len(r))
        if not cands:
            return {"error": "donor %s has no drug-like ligand" % donor_pdb_id}
        ligand_comp_id = max(cands, key=lambda k: cands[k])

    st_t = gemmi.read_structure(target_path)
    st_t.setup_entities()
    tchains = [c for c in st_t[0] if c.get_polymer()] or list(st_t[0])
    if target_chain is not None:
        want_t = [c for c in tchains if c.name == target_chain]
        if not want_t:
            # Never substitute another chain: the caller named one, and a
            # different chain would superpose and transfer into the wrong copy.
            return {"error": "target chain %s not found in %s (polymer chains: %s)"
                             % (target_chain, os.path.basename(target_path),
                                ", ".join(c.name for c in tchains) or "none")}
        tch = want_t[0]
    else:
        tch = tchains[0]
    # Donor chain must be one that actually contacts the ligand, else the
    # superposition aligns the wrong copy of a multimer.
    lig_res = [(ch.name, r) for ch in st_d[0] for r in ch
               if r.name.upper() == ligand_comp_id.upper()]
    if not lig_res:
        return {"error": "ligand %s not found in donor %s" % (ligand_comp_id, donor_pdb_id)}
    if donor_chain is None:
        ns = gemmi.NeighborSearch(st_d[0], st_d.cell, 6.0).populate()
        counts = {}
        for _, r in lig_res[:1]:
            for a in r:
                for m in ns.find_atoms(a.pos, "\0", radius=6.0):
                    cra = m.to_cra(st_d[0])
                    info = gemmi.find_tabulated_residue(cra.residue.name)
                    if info and info.is_amino_acid():
                        counts[cra.chain.name] = counts.get(cra.chain.name, 0) + 1
        donor_chain = max(counts, key=counts.get) if counts else lig_res[0][0]
    want_d = [c for c in st_d[0] if c.name == donor_chain]
    if not want_d:
        return {"error": "donor chain %s not found in %s (chains: %s)"
                         % (donor_chain, donor_pdb_id,
                            ", ".join(c.name for c in st_d[0]))}
    dch = want_d[0]

    pt, pd_ = tch.get_polymer(), dch.get_polymer()
    if not len(pt) or not len(pd_):
        return {"error": "no polymer to superpose (target chain %s / donor chain %s)"
                         % (tch.name, dch.name)}
    # Restrict the fit to the residues the DONOR observes. A donor that
    # crystallised one domain (WRN helicase, 516-946) against a full-length
    # model would otherwise be fitted globally: the spare domains drag the fit,
    # and a 1-2 A local match comes back as 16 A with a wrong-looking site.
    # The correspondence comes from SEQUENCE alignment, not residue numbers,
    # because author numbering is a deposition choice the two need not share.
    pairs, seq_identity = residue_correspondence(tch, dch)
    if not pairs:
        return {"error": "no sequence correspondence between target chain %s and "
                         "donor %s chain %s — is the donor the same protein?"
                         % (tch.name, donor_pdb_id, dch.name)}
    if seq_identity < 0.5:
        return {"error": "target chain %s and donor %s chain %s align at only %.0f%% "
                         "identity over %d residues — too divergent to transfer a site"
                         % (tch.name, donor_pdb_id, dch.name, 100 * seq_identity, len(pairs))}
    trimmed = False
    t_aa = aa_residues(tch)
    if len(pairs) < len(t_aa) - 5:
        keep_t_pairs, sub_t_poly = span_residues(tch, [t for _, t in pairs])
        keep_d_pairs, sub_d_poly = span_residues(dch, [d for d, _ in pairs])
        if len(sub_t_poly) >= 20 and len(sub_d_poly) >= 20:
            pt, pd_ = sub_t_poly, sub_d_poly
            trimmed = True
    sup = gemmi.calculate_superposition(pt, pd_, pt.check_polymer_type(),
                                        gemmi.SupSelect.CaP, trim_cycles=3)
    # A multi-domain protein can match locally and still refuse a global fit:
    # WRN's AlphaFold model has the helicase sub-domains hinged differently from
    # the crystal (78% of 20-residue windows agree under 2 A, yet no whole-domain
    # fit beats ~12 A). Docking only needs the pocket's own frame, so when the
    # global fit is poor, re-fit on the shell of residues that line the donor's
    # ligand — that is the geometry the transferred box depends on.
    site_local = False
    if sup.rmsd > 2.5:
        ns0 = gemmi.NeighborSearch(st_d[0], st_d.cell, 12.0).populate()
        shell = set()
        for a in copy_res_probe(st_d, ligand_comp_id, donor_chain):
            for m in ns0.find_atoms(a, "\0", radius=12.0):
                cra = m.to_cra(st_d[0])
                info = gemmi.find_tabulated_residue(cra.residue.name)
                if info and info.is_amino_acid() and cra.chain.name == donor_chain:
                    shell.add((cra.residue.seqid.num, cra.residue.name))
        # Carry the shell across via the sequence alignment, not residue numbers.
        shell_pairs = [(d, t) for d, t in pairs if (d.seqid.num, d.name) in shell]
        if len(shell_pairs) >= 20:
            keep_t, sub_t = span_residues(tch, [t for _, t in shell_pairs])
            keep_d, sub_d = span_residues(dch, [d for d, _ in shell_pairs])
            if len(sub_t) >= 20 and len(sub_d) >= 20:
                s2 = gemmi.calculate_superposition(sub_t, sub_d, sub_t.check_polymer_type(),
                                                   gemmi.SupSelect.CaP, trim_cycles=3)
                if s2.rmsd < sup.rmsd:
                    sup, site_local = s2, True
    # Move the chosen ligand copy into the target's frame.
    copy_res = [r for cn, r in lig_res if cn == donor_chain] or [lig_res[0][1]]
    lig_pos = [gemmi.Position(*sup.transform.apply(a.pos)) for a in copy_res[0]]

    found = {}
    for ch in st_t[0]:
        if target_chain and ch.name != target_chain:
            continue
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if not (info and info.is_amino_acid()):
                continue
            d = min((a.pos.dist(lp) for a in res for lp in lig_pos), default=None)
            if d is not None and d <= radius:
                found[(ch.name, res.seqid.num, res.name)] = round(d, 2)
    pocket = [{"chain": c, "seq_id": n, "comp_id": nm, "min_dist": d}
              for (c, n, nm), d in sorted(found.items(), key=lambda kv: kv[0][1])]

    xs = [p.x for p in lig_pos]
    ys = [p.y for p in lig_pos]
    zs = [p.z for p in lig_pos]
    box = {"center": [round(sum(v) / len(v), 2) for v in (xs, ys, zs)],
           "size": [round(max(v) - min(v) + 2 * padding, 2) for v in (xs, ys, zs)],
           "n_ligand_atoms": len(lig_pos)}

    res = {"target": target_path, "donor": donor_pdb_id, "donor_chain": dch.name,
           "target_chain": tch.name, "ligand": ligand_comp_id,
           "rmsd": round(sup.rmsd, 3), "n_aligned": sup.count,
           "superposed_on": ("site shell (12 A around the donor ligand) — the global "
                             "fit was worse, indicating different domain packing"
                             if site_local else
                             "%d residues aligned to the donor's observed span"
                             % len(pairs) if trimmed else "full chain"),
           "site_local_fit": site_local,
           "seq_identity": round(seq_identity, 3),
           "n_corresponding_residues": len(pairs),
           "pocket": pocket, "box": box}

    # pLDDT at the site — the number that matters for docking into a model.
    if plddt is None and accession:
        try:
            plddt = plddt_by_residue(accession, out_dir=out_dir)
        except Exception:
            plddt = None
    if plddt:
        vals = [plddt[p["seq_id"]] for p in pocket if p["seq_id"] in plddt]
        if vals:
            res["site_plddt"] = {"mean": round(sum(vals) / len(vals), 1),
                                 "min": round(min(vals), 1),
                                 "n_below_70": sum(1 for v in vals if v < 70)}
    # Confidence is calibrated on an 8-target benchmark (KRAS, EGFR, WRN, BRAF,
    # CDK2, HSP90AA1, PARP1, BTK): transfer onto each AlphaFold model from an
    # auto-selected donor, scored against that donor's crystallographic site.
    # Superposition RMSD is the signal that tracks recovery; site pLDDT is not.
    #   'high'   (6/8 cases): mean recovery 0.91, min 0.85
    #   'medium' (1/8 cases): mean recovery 0.56
    #   'low'    (1/8 cases): mean recovery 0.00
    # HSP90AA1 is the cautionary row: site pLDDT 89.8 with 33.8 A RMSD, 0.502
    # donor identity and 0.00 recovery, so a confident-looking model says
    # nothing about whether the donor's pocket transferred correctly. pLDDT
    # only downgrades — it answers "is the model reliable HERE", never "did
    # the transfer work".
    sp = res.get("site_plddt") or {}
    if seq_identity < 0.8:
        # A donor that aligns poorly is a different construct or a paralog:
        # HSP90AA1's auto-selected donor aligns at 0.50 identity and recovers
        # 0.00 of the site.
        res["confidence"] = "low"
    elif sup.rmsd > 3.0:
        res["confidence"] = "low"
    elif sup.rmsd <= 1.0:
        res["confidence"] = "low" if sp and sp["mean"] < 60 else "high"
    elif sup.rmsd <= 2.5:
        res["confidence"] = "low" if sp and sp["mean"] < 60 else "medium"
    else:
        res["confidence"] = "low"
    res["confidence_basis"] = ("superposition RMSD %.2f A; donor identity %.0f%%%s"
                              % (sup.rmsd, 100 * seq_identity,
                                 "" if not sp else
                                 "; site pLDDT mean %.1f" % sp["mean"]))

    if write_complex:
        st_out = gemmi.read_structure(target_path)
        st_out.setup_entities()
        # Pick a chain id not already present: a hardcoded "L" merges the ligand
        # into an existing chain L (antibody light chains) with no error.
        used = {c.name for c in st_out[0]}
        lig_chain_id = next(
            (n for n in ["LIG", "L", "X", "Z", "Y", "W"] + [chr(c) for c in range(65, 91)]
             if n not in used), None)
        if lig_chain_id is None:
            lig_chain_id = next("L%d" % i for i in range(1, 999) if "L%d" % i not in used)
        ch_out = gemmi.Chain(lig_chain_id)
        nr = gemmi.Residue()
        nr.name = ligand_comp_id
        nr.seqid = gemmi.SeqId(9001, " ")
        nr.het_flag = "H"
        for a, p in zip(copy_res[0], lig_pos):
            na = gemmi.Atom()
            na.name, na.element, na.pos = a.name, a.element, p
            na.occ, na.b_iso = a.occ, a.b_iso
            nr.add_atom(na)
        ch_out.add_residue(nr)
        st_out[0].add_chain(ch_out)
        st_out.setup_entities()
        base = os.path.splitext(os.path.basename(target_path))[0]
        out = os.path.join(out_dir, "%s_with_%s_from_%s.cif"
                           % (base, ligand_comp_id, donor_pdb_id.upper()))
        os.makedirs(out_dir, exist_ok=True)
        st_out.make_mmcif_document().write_file(out)
        res["complex_path"] = out
        res["complex_ligand_chain"] = lig_chain_id
    return res


def structure_caveats(result):
    """Read a get_structure() result and list what would bite downstream.

    Returns a list of human-readable warning strings — empty means nothing
    obvious. Always surface these to the user; a returned file is not the
    same as a file that is fit for the intended purpose.
    """
    w = []
    af = result.get("alphafold") or {}
    if af.get("has_model"):
        g = af.get("global_plddt")
        frac = af.get("plddt_fractions") or {}
        low = (frac.get("low_50_70") or 0) + (frac.get("very_low_lt50") or 0)
        w.append("PREDICTED model, not experimental — no ligand, no bound-state "
                 "conformation, side chains unreliable for docking.")
        if g is not None and g < 70:
            w.append("Global pLDDT %.1f is low; treat the fold as a hypothesis." % g)
        if low > 0.3:
            w.append("%.0f%% of residues are below pLDDT 70 (likely disordered "
                     "or flexible regions) — do not dock into them." % (100 * low))
    cand = (result.get("candidates") or [None])[0]
    if cand:
        if cand.get("state") == "apo":
            w.append("Chosen entry %s is APO (no drug-like ligand): the pocket may "
                     "be in a closed/unliganded conformation." % cand["pdb_id"])
        if cand.get("mutations"):
            w.append("Entry %s carries engineered mutations (%s) — confirm they do "
                     "not sit in your site of interest." % (cand["pdb_id"], cand["mutations"]))
        cov = cand.get("coverage")
        if cov is not None and cov < 0.8:
            r = cand.get("unp_range")
            w.append("Only %.0f%% of the UniProt sequence is observed (residues "
                     "%s-%s); the rest is absent from the file." % (100 * cov, r[0], r[1]))
        reso = cand.get("resolution")
        if reso and reso > 2.8:
            w.append("Resolution %.2f A is modest — side-chain placement is "
                     "uncertain for docking." % reso)
        if (cand.get("n_protein_chains") or 1) > 1:
            w.append("Entry has %d protein chains; a pocket may sit at an "
                     "interface or be duplicated across copies."
                     % cand["n_protein_chains"])
    tr = result.get("pocket_transfer") or {}
    if tr and not tr.get("error"):
        sp = tr.get("site_plddt") or {}
        msg = ("Binding site was NOT observed in this structure — it was transferred "
               "from %s (ligand %s) by superposition at %.2f A RMSD over %d residues. "
               "The site is a homology-positioned hypothesis: it locates a pocket that "
               "exists in the donor, and is silent about one that does not."
               % (tr["donor"], tr["ligand"], tr["rmsd"], tr["n_aligned"]))
        if sp:
            msg += (" Model confidence AT THE SITE: pLDDT mean %.1f, min %.1f"
                    % (sp["mean"], sp["min"]))
            if sp.get("n_below_70"):
                msg += " (%d site residues below 70)" % sp["n_below_70"]
            msg += "."
        w.append(msg)
        if tr.get("confidence") == "low":
            w.append("Transferred site is LOW confidence (RMSD %.2f A%s) — do not dock "
                     "into it without independent support."
                     % (tr["rmsd"], "" if not sp else ", site pLDDT %.1f" % sp["mean"]))
    elif tr.get("error"):
        w.append("No binding site available: %s" % tr["error"])
    cg = result.get("coverage_gaps") or {}
    for g in (cg.get("gaps") or []):
        if g["kind"] != "internal":
            continue
        where = "%d-%d" % (g["unp_start"], g["unp_end"])
        near = g.get("near_site")
        dist = g.get("flank_dist_to_ligand_A")
        msg = "Residues %s (%d aa) are UNRESOLVED inside the observed region" % (
            where, g["length"])
        if near:
            msg += (" and their flanks sit %.1f A from the ligand — this loop likely "
                    "lines the site; docking here may be unreliable" % dist)
        elif dist is not None:
            msg += " (flanks %.1f A from the ligand)" % dist
        if g.get("verdict") == "af_can_fill":
            msg += ". AlphaFold is confident there (pLDDT %.0f), so the loop can be modelled in." % g["af_mean_plddt"]
        elif g.get("verdict") == "disordered_in_af":
            msg += ". AlphaFold is also low-confidence there (pLDDT %.0f) — likely genuinely flexible." % g["af_mean_plddt"]
        w.append(msg + ("" if msg.endswith(".") else "."))
    if cg.get("n_missing") and not (cg.get("gaps") or []):
        w.append("%d residues are unresolved in the chosen chain." % cg["n_missing"])
    if len(result.get("ligand_copies") or []) > 1:
        w.append("Ligand %s occurs in %d copies — pocket/box use copy in chain %s only."
                 % (result.get("pocket_ligand"), len(result["ligand_copies"]),
                    result.get("pocket_chain")))
    if result.get("ambiguous_uniprot"):
        w.append("Query also matched %s — confirm the accession is the intended isoform/paralog."
                 % ", ".join(result["ambiguous_uniprot"]))
    if (result.get("prep") or {}).get("format_coerced"):
        w.append(result["prep"]["format_coerced"])
    return w


def get_pdb_entry(pdb_id, out_dir="structures", fmt="cif", prep=True, chains=None):
    """Fetch and prep ONE known PDB id, with its ligand inventory.

    Use when the user already named an entry; skips identifier resolution
    and ranking entirely.
    """
    pid = pdb_id.upper()
    meta = rcsb_entries([pid])
    if pid not in meta:
        return {"error": "PDB entry %s not found" % pid}
    s = summarize_entry(meta[pid])
    res = {"chosen": pid, "candidates": [s], "entry": s}
    raw = fetch_pdb(pid, out_dir=out_dir, fmt=fmt)
    res["raw_path"] = raw
    res["path"] = raw
    if prep:
        # Derive from the file actually fetched — fetch_pdb falls back to
        # mmCIF when RCSB has no legacy .pdb for an oversized entry.
        got = "pdb" if raw.endswith(".pdb") and fmt == "pdb" else "cif"
        out = os.path.join(out_dir, "%s_prep.%s" % (pid, got))
        res["prep"] = prep_structure(raw, out, chains=chains or (s["target_chains"] or None))
        res["path"] = res["prep"]["path"]
        if s["drug_like_ligands"]:
            lig = s["drug_like_ligands"][0]
            copies = ligand_copies(res["path"], lig)
            res["ligand_copies"] = copies
            ch = copies[0]["chain"] if copies else None
            res["pocket_ligand"], res["pocket_chain"] = lig, ch
            res["pocket"] = pocket_residues(res["path"], lig, chain=ch)
            res["box"] = ligand_box(res["path"], lig, chain=ch)
    return res


def get_structure(query, organism_id=9606, prefer_ligand=True, out_dir="structures",
                  fmt="cif", prep=True, top_n=10, min_coverage=0.0, fallback_alphafold=True,
                  coverage_gaps=True):
    """One-call entry: identifier -> ranked table -> best file on disk.

    Returns dict with `uniprot`, `candidates` (top_n ranked dicts),
    `chosen`, `path`, `prep` report and `alphafold`. Falls back to
    AlphaFold DB when no experimental structure passes the filters.
    """
    hits = resolve_uniprot(query, organism_id=organism_id)
    if not hits:
        return {"error": "no UniProt match for %r" % query, "uniprot": None}
    up = hits[0]
    acc = up["accession"]
    ranked = rank_structures(acc, prefer_ligand=prefer_ligand, min_coverage=min_coverage)
    res = {"uniprot": up, "n_pdb_entries": len(ranked),
           "candidates": ranked[:top_n], "ambiguous_uniprot": [h["accession"] for h in hits[1:4]]}
    if ranked:
        best = ranked[0]
        res["chosen"] = best["pdb_id"]
        raw = fetch_pdb(best["pdb_id"], out_dir=out_dir, fmt=fmt)
        res["raw_path"] = raw
        res["path"] = raw
        if prep:
            got = "pdb" if raw.endswith(".pdb") and fmt == "pdb" else "cif"
            out = os.path.join(out_dir, "%s_prep.%s" % (best["pdb_id"], got))
            res["prep"] = prep_structure(raw, out, chains=best["target_chains"] or None)
            out = res["prep"]["path"]          # may be coerced to .cif
            res["path"] = out
            if best["drug_like_ligands"]:
                lig = best["drug_like_ligands"][0]
                copies = ligand_copies(out, lig)
                res["ligand_copies"] = copies
                use_chain = copies[0]["chain"] if copies else None
                res["pocket"] = pocket_residues(out, lig, chain=use_chain)
                res["box"] = ligand_box(out, lig, chain=use_chain)
                res["pocket_ligand"] = lig
                res["pocket_chain"] = use_chain
        if coverage_gaps:
            try:
                res["coverage_gaps"] = coverage_gap_report(
                    best["pdb_id"], acc, chain=res.get("pocket_chain"),
                    ligand_comp_id=res.get("pocket_ligand"),
                    struct_path=res.get("path"))
            except Exception as exc:
                res["coverage_gaps"] = {"error": "%s: %s" % (type(exc).__name__, exc)}
    elif fallback_alphafold:
        res["alphafold"] = fetch_alphafold(acc, out_dir=out_dir, fmt=fmt)
        if res["alphafold"].get("has_model"):
            res["path"] = res["alphafold"]["path"]
            res["chosen"] = res["alphafold"]["entry_id"]
    return res


def get_structure_with_site(query, organism_id=9606, out_dir="structures",
                            donor_pdb_id=None, **kw):
    """get_structure(), and if the result has no bound ligand, borrow a site.

    The experimental path already returns a pocket and box when the chosen entry
    is holo. This wrapper covers the two cases where it does not: an AlphaFold
    fallback (no ligand by construction) and an apo experimental pick. In both it
    calls transfer_pocket() and attaches `pocket_transfer` to the result.

    Use when the downstream step needs a site regardless of what was available —
    e.g. handing a receptor plus box to a docking run.
    """
    res = get_structure(query, organism_id=organism_id, out_dir=out_dir, **kw)
    if res.get("error") or not res.get("path"):
        return res
    if res.get("pocket"):
        return res                     # experimental holo site already present
    acc = (res.get("uniprot") or {}).get("accession")
    if not acc:
        return res
    tr = transfer_pocket(res["path"], donor_pdb_id=donor_pdb_id, accession=acc,
                         out_dir=out_dir)
    res["pocket_transfer"] = tr
    if not tr.get("error"):
        res["pocket"] = tr["pocket"]
        res["box"] = tr["box"]
        res["pocket_ligand"] = tr["ligand"]
        res["pocket_source"] = "transferred from %s (RMSD %.2f A)" % (tr["donor"], tr["rmsd"])
    return res


# --------------------------------------------------------------------------
# Run output: results/<topic>/<run>/
# --------------------------------------------------------------------------
#
# get_structure()/get_pdb_entry()/fetch_pdb()/fetch_alphafold() all default
# out_dir="structures" for standalone use (e.g. handing a receptor straight
# to a docking tool without saving a run). get_structure_to_results() is the
# wrapper that routes a full run into results/<topic>/<run>/ instead, so a
# caller does not have to pass out_dir by hand and remember to write a README.
RESULTS_TOPIC = "protein_structure"
CANDIDATES_CSV = "ranked_structures.csv"
POCKET_CSV = "pocket_residues.csv"
COVERAGE_CSV = "coverage_gaps.csv"
RUN_README = "README.md"
SUMMARY_MAX_WORDS = 130

CANDIDATE_COLUMNS = ("pdb_id", "method", "resolution", "r_free", "state",
                     "coverage", "n_mutations", "score")


def gps_run_dir(target, root="results", topic=None, make=True):
    """Path for one structure run: <root>/<topic>/<slug>/, with scripts/ beside it.

    `target` is the query passed to get_structure() (gene symbol, accession,
    protein name); slugged to snake_case because result directories are named
    that way and a raw query may carry spaces or case.

    `topic` defaults to RESULTS_TOPIC through an explicit None check rather
    than in the signature: the kernel.py sidecar loader rejects a non-literal
    default, and a rejected file defines none of this module's helpers.
    """
    topic = RESULTS_TOPIC if topic is None else topic
    slug = re.sub(r"[^a-z0-9]+", "_", str(target).strip().lower()).strip("_")
    if not slug:
        raise ValueError("gps_run_dir got an empty target name; the run "
                         "directory is named after the query")
    out_dir = os.path.join(root, topic, slug)
    if make:
        os.makedirs(os.path.join(out_dir, "scripts"), exist_ok=True)
    return out_dir


def gps_write_table(path, rows, headers=None):
    """Write rows (list of dicts) to CSV. Returns the path, or None if empty."""
    rows = [dict(r) for r in (rows or [])]
    if not rows:
        return None
    headers = list(headers or [])
    extra = sorted({k for r in rows for k in r} - set(headers))
    cols = [c for c in headers if any(c in r for r in rows)] + extra
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: ("" if row.get(c) is None else row.get(c))
                             for c in cols})
    return path


def gps_check_words(text, cap, label):
    """Enforce a word cap on run-README prose, naming the count and the cap."""
    n = len(str(text or "").split())
    if n > cap:
        raise ValueError(
            f"{label} is {n} words; cap is {cap}. Move detail into the table "
            "rather than the prose.")
    return n


def gps_run_readme(name, result, summary, files=(), data_sources=(),
                   limits=(), title=None):
    """Render the run README for one get_structure() result. Returns markdown.

    `result` is whatever get_structure() / get_structure_with_site() /
    get_pdb_entry() returned. Sections follow `coding-standards`
    (Result, Files, Data sources, Limits). The Limits section is
    `structure_caveats(result)` verbatim, plus any run-specific ones passed
    in — those caveats are already computed by this skill and restating them
    by hand would risk drifting from what the function actually checks.
    """
    if result.get("error"):
        raise ValueError(f"gps_run_readme got a result with an error: "
                         f"{result['error']!r}; there is no structure to report")
    gps_check_words(summary, SUMMARY_MAX_WORDS, "summary")
    if str(summary or "").strip().startswith("-"):
        raise ValueError("summary must be plain prose, not bullets — it is "
                         "the paragraph a non-specialist reads first")

    up = result.get("uniprot") or {}
    default_title = f"{up.get('gene') or name} structure"
    parts = [f"# {title or default_title}", "", str(summary).strip(), "",
             "## Result", ""]

    cand = (result.get("candidates") or [None])[0]
    af = result.get("alphafold") or {}
    if cand:
        parts += [
            "| quantity | value |", "| --- | --- |",
            f"| chosen entry | {result.get('chosen')} |",
            f"| method | {cand.get('method')} |",
            f"| resolution | {cand.get('resolution')} |",
            f"| state | {cand.get('state')} |",
            f"| coverage | {cand.get('coverage')} |",
            f"| n_mutations | {cand.get('n_mutations')} |", ""]
    elif af.get("has_model"):
        parts += [
            "| quantity | value |", "| --- | --- |",
            f"| chosen entry | AlphaFold {result.get('chosen')} (predicted) |",
            f"| global pLDDT | {af.get('global_plddt')} |", ""]

    if result.get("pocket_transfer") and not result["pocket_transfer"].get("error"):
        tr = result["pocket_transfer"]
        parts += [f"Binding site transferred from {tr.get('donor')} "
                 f"(ligand {tr.get('ligand')}) at {tr.get('rmsd')} A RMSD over "
                 f"{tr.get('n_aligned')} residues.", ""]

    parts += ["## Files", ""]
    for entry in files:
        f_name, description = (entry["name"], entry["description"]) \
            if isinstance(entry, dict) else entry
        parts.append(f"- `{f_name}` — {description}")
    parts += ["", "## Data sources", ""] + [f"- {s}" for s in data_sources]
    parts += ["", "## Limits", ""]
    caveats = structure_caveats(result)
    if not caveats and not limits:
        parts.append("- None flagged by `structure_caveats()`.")
    else:
        parts += [f"- {s}" for s in caveats + list(limits)]
    return "\n".join(parts).rstrip() + "\n"


def gps_write_run(out_dir, name, result, summary=None, files=(),
                  data_sources=(), limits=(), scripts=()):
    """Write one structure run directory. Returns {name: path} for what was written.

    Copies the structure file(s) already on disk at `result["path"]` (and
    `result["raw_path"]` when different) into `out_dir`, writes the ranked
    candidates table, the pocket-residues table (when a pocket was found),
    the coverage-gap table (when computed), the run README, and copies
    `scripts` into `scripts/`.
    """
    if result.get("error"):
        raise ValueError(f"gps_write_run got a result with an error: "
                         f"{result['error']!r}; nothing to write")
    os.makedirs(os.path.join(out_dir, "scripts"), exist_ok=True)
    written = {}

    for key in ("path", "raw_path"):
        src = result.get(key)
        if not src or not os.path.isfile(src):
            continue
        dst = os.path.join(out_dir, os.path.basename(src))
        if os.path.abspath(src) != os.path.abspath(dst):
            shutil.copyfile(src, dst)
        written[key] = dst

    cand_path = gps_write_table(os.path.join(out_dir, CANDIDATES_CSV),
                                result.get("candidates"), CANDIDATE_COLUMNS)
    if cand_path:
        written["candidates"] = cand_path

    pocket_path = gps_write_table(os.path.join(out_dir, POCKET_CSV),
                                  result.get("pocket"),
                                  ["chain", "seq_id", "comp_id", "min_dist"])
    if pocket_path:
        written["pocket"] = pocket_path

    gaps = (result.get("coverage_gaps") or {}).get("gaps")
    gaps_path = gps_write_table(os.path.join(out_dir, COVERAGE_CSV), gaps)
    if gaps_path:
        written["coverage_gaps"] = gaps_path

    for src in scripts:
        dst = os.path.join(out_dir, "scripts", os.path.basename(src))
        shutil.copyfile(src, dst)
        written.setdefault("scripts", []).append(dst)

    listed = []
    if "path" in written:
        listed.append((os.path.basename(written["path"]),
                       "the prepared structure: cleaned, chain-selected, "
                       "ready for downstream use"))
    if "raw_path" in written and written["raw_path"] != written.get("path"):
        listed.append((os.path.basename(written["raw_path"]), "the unprocessed file as downloaded"))
    if "candidates" in written:
        listed.append((CANDIDATES_CSV, "every ranked candidate entry considered, "
                                       "not just the chosen one"))
    if "pocket" in written:
        listed.append((POCKET_CSV, "protein residues within the pocket radius "
                                   "of the bound or transferred ligand"))
    if "coverage_gaps" in written:
        listed.append((COVERAGE_CSV, "unresolved regions in the chosen chain, "
                                     "with AlphaFold fill confidence where "
                                     "relevant"))
    listed += list(files)

    if summary is None:
        raise ValueError("gps_write_run needs a summary — the plain-prose "
                         "paragraph a non-specialist reads first")
    readme = gps_run_readme(name, result, summary, files=listed,
                            data_sources=data_sources, limits=limits)
    readme_path = os.path.join(out_dir, RUN_README)
    with open(readme_path, "w", encoding="utf-8") as fh:
        fh.write(readme)
    written["readme"] = readme_path
    return written


def get_structure_to_results(query, name=None, root="results",
                             topic=None, summary=None, files=(),
                             data_sources=(), limits=(), scripts=(),
                             with_site=False, donor_pdb_id=None, **kwargs):
    """get_structure() (or get_structure_with_site()), saved as a full run.

    Routes the fetch into `<root>/<topic>/<slug>/` via `gps_run_dir` instead
    of the bare "structures" default, then writes the run (tables, README
    with `structure_caveats()` as its Limits) via `gps_write_run`. `name`
    defaults to `query`. Pass `with_site=True` to call
    `get_structure_with_site()` when the downstream step needs a pocket
    regardless of what the chosen entry has.

    `topic` defaults to RESULTS_TOPIC through an explicit None check rather
    than in the signature: the kernel.py sidecar loader rejects a non-literal
    default, and a rejected file defines none of this module's helpers.

    Returns the `get_structure()`-shaped result dict with a `run` key added:
    the `{name: path}` dict `gps_write_run` returned.
    """
    topic = RESULTS_TOPIC if topic is None else topic
    out_dir = gps_run_dir(name or query, root=root, topic=topic)
    if with_site:
        res = get_structure_with_site(query, out_dir=out_dir,
                                      donor_pdb_id=donor_pdb_id, **kwargs)
    else:
        res = get_structure(query, out_dir=out_dir, **kwargs)
    if res.get("error"):
        return res
    res["run"] = gps_write_run(out_dir, name or query, res, summary=summary,
                               files=files, data_sources=data_sources,
                               limits=limits, scripts=scripts)
    return res
