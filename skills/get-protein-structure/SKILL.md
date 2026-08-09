---
name: get-protein-structure
description: Get a protein 3D structure from the PDB or AlphaFold DB, starting from a gene symbol, UniProt accession, protein name, or PDB id. Ranks every experimental entry by resolution, UniProt coverage, apo/holo state and engineered mutations, then writes a cleaned chain-selected file ready for docking or design — with the bound-ligand inventory, binding-site residues, a docking box, and a residue-level coverage map flagging unresolved loops that line the pocket. Falls back to AlphaFold DB with pLDDT when no experimental structure exists, and can transfer a binding site onto a ligand-free or predicted model from a liganded relative with a calibrated confidence score. Use to fetch, download, find or prepare a structure for a target, to check whether a liganded (holo) structure of a pocket exists, to pick the best PDB entry, or to get a receptor plus box for DiffDock/Boltz/ProteinMPNN. For predicting a structure that does not exist, use alphafold2, boltz, chai1, or openfold3 instead.
---

# Get protein structure (PDB / AlphaFold)

Turn an identifier into a **structure file you can actually use**. The hard part
is never the download — it is choosing among hundreds of entries and knowing
what the chosen file is missing.

Helpers are preloaded in the Python kernel by `kernel.py` (no import needed).
Needs `gemmi`: `manage_packages(mode="install", environment=<env>, packages=["gemmi"])`.

## Decide the route first

| Situation | Call |
|---|---|
| Gene symbol / accession / protein name, want the best structure | `get_structure("KRAS")` |
| User named a specific entry | `get_pdb_entry("6OIM")` |
| Want the ranked table to choose from yourself | `rank_structures(acc)` |
| Only need the predicted model | `fetch_alphafold(acc)` |
| Structure does not exist anywhere | **wrong skill** — use `boltz` / `chai1` / `alphafold2` / `openfold3` |

## Workflow

```python
res = get_structure("KRAS")                 # human by default
res["uniprot"]        # {accession, gene, name, length, organism}
res["candidates"]     # ranked dicts, best first
res["chosen"]         # e.g. "9IAY"
res["path"]           # cleaned file ready for downstream use
res["pocket"]         # binding-site residues around the drug-like ligand
res["box"]            # {center, size} docking box

for c in structure_caveats(res):
    print(c)          # ALWAYS report these
```

Then hand `res["path"]` to `diffdock`, `boltz`, `proteinmpnn`, or a viewer
(save it as an artifact — `.pdb`/`.cif` render in an interactive 3D viewer).

To save the run instead of a one-off file (structure files, the ranked
candidates table, the pocket-residues table, and a README whose Limits
section is `structure_caveats()` verbatim), use `get_structure_to_results`,
not `get_structure` — it saves into `results/protein_structure/<target>/`
rather than the bare `structures/` default. Needs `$SCIENCE_RESULTS_ROOT` set
to this repo's `results/` directory, or an explicit `root=` — it raises
naming the variable rather than silently writing `results/` wherever the
session's cwd happens to be:

```python
res = get_structure_to_results("KRAS", summary="...",
                               data_sources=["RCSB PDB"])
res["run"]             # {name: path} for what get_structure_to_results wrote
```

Pass `with_site=True` to save a run built from `get_structure_with_site`
instead, when the downstream step needs a pocket regardless of what the
chosen entry has.

## Report the caveats — every time

`structure_caveats(res)` is not optional polish. A returned file is *not* the
same as a file fit for the purpose, and the failure modes are silent: docking
into an apo conformation, into a predicted disordered loop, or into a construct
whose engineered mutation sits in your pocket all produce confident nonsense.
Surface apo/holo state, coverage, mutations, resolution and pLDDT in your
answer to the user, not just the filename.

## How ranking works, and when to override it

Score blends resolution (40), UniProt coverage (30), a **+20 bonus for a
drug-like ligand**, a penalty per engineered mutation, and small bonuses for
X-ray and R-free < 0.25. Override via arguments:

- `prefer_ligand=False` — you want the apo/unliganded form.
- `min_coverage=0.9` — force near-full-length entries.
- `methods=["X-ray diffraction"]` — exclude EM/NMR. Either vocabulary is
  accepted (`"X-ray"` or `"X-ray diffraction"`, `"EM"` or
  `"Electron Microscopy"`): SIFTS and RCSB name the same experiment
  differently, so the filter compares method families rather than raw
  strings.
- `organism_id=10090` — mouse instead of human; `None` for any species.

**Coverage is a fraction of the FULL UniProt sequence, and for multi-domain
targets the useful entries score low on it.** EGFR has 531 chains at coverage
≈ 0.28 (the kinase domain, down to 1.07 Å) and only 12 at 1.00 (full-length EM
at ~3.1-3.6 Å). The default correctly returns the high-resolution kinase
domain; if the user needs the ectodomain or a full-length model, that is
`min_coverage=0.9` — a different question, not a better rank. Read the
`coverage` and `unp_range` columns before trusting the top row.

## Apo vs holo is a heuristic, so check it

The PDB exposes no curated "this ligand is the point of the structure" flag, so
`classify_ligand` buckets each chem-comp id as `water` / `additive` (ions,
glycerol, PEG, buffers) / `cofactor` (ATP, GDP, NAD, HEM...) / `drug_like`, and
`state` is holo iff something drug-like is present. Consequences:

- A **cofactor-only** entry (e.g. KRAS + GDP) reads as `apo` — correct for
  inhibitor docking, wrong if the nucleotide *is* your ligand.
- A novel ligand whose id happens to sit in the additive list would be missed.
  Inspect `candidates[i]["ligands"]` (each with `comp_id`, `name`, `mw`,
  `smiles`, `role`) when the call matters.
- Classification drives **ranking only**. Prep keeps structural ions and
  cofactors regardless.

## Prep: what it does and does not do

`prep_structure(in_path, out_path, chains=["A"])` removes waters, alt-locs,
extra NMR/EM models, hydrogens and crystallisation additives; **keeps**
structural ions (ZN, MG, MN, CA, FE, CU, NI, CO, K, NA), cofactors and
drug-like ligands. It does **not** protonate, add missing side chains, build
missing loops, or assign bond orders — use a dedicated preparation tool
(PDBFixer / Protoss / Schrödinger) if your downstream method needs those.
Missing residues stay missing; check `n_residues` against `uniprot.length`.

**Legacy PDB format corrupts silently — gemmi does not raise.** `prep_structure`
refuses to write `.pdb` when the structure would not survive it, writes mmCIF
instead, and explains why in `prep["format_coerced"]` /
`format_coerced_reasons`. Four triggers:

- **5-character CCD ids.** Since 2023 the PDB issues them (`A1BEA`, `A1I1R`);
  the 3-char resName field truncates `A1BEA` to `A1B`, breaking every
  downstream comp_id lookup.
- **> 62 chains** and **> 99,999 atoms** — the single-character chainID and
  5-digit serial ceilings. Verified: the 89-chain ribosome 4V6X written through
  `make_pdb_string()` collapses to **53 distinct chain characters** with serials
  overflowing to `A2YB1`, with no exception raised.
- **Multi-character chain ids** (`Aa`, `CB`), which do not fit the 1-char field.

**So always read the returned `prep["path"]` — never assume the extension you
asked for.** `fmt="pdb"` is honoured for ordinary structures and overridden only
when it would lose data. RCSB additionally serves no legacy `.pdb` at all for
oversized entries; `fetch_pdb` falls back to mmCIF rather than failing.

## Residue-level coverage: the gap that eats a docking run

Resolution and coverage are *summary* numbers and they hide the failure that
actually wastes a week: **a loop that forms part of the pocket is absent from
the crystal.** `get_structure` runs this automatically into
`res["coverage_gaps"]`; call it directly with
`coverage_gap_report(pdb_id, accession, ligand_comp_id=..., struct_path=...)`.

Each gap is classified `internal` (unresolved *inside* the observed region —
the dangerous kind), `terminal` (truncated N/C terminus, usually harmless), or
`expression_tag` (construct additions numbered < 1 in UniProt space — not
missing protein). Internal gaps are then annotated with two things that decide
what to do about them:

- **`flank_dist_to_ligand_A` / `near_site`** — how far the gap's flanking
  residues sit from the bound ligand (default threshold 8 Å).
- **`af_mean_plddt` / `verdict`** — mean AlphaFold pLDDT over the missing
  stretch. `af_can_fill` (≥ 70) means a predicted model can credibly supply the
  loop; `disordered_in_af` (< 70) means it is probably genuinely flexible and
  no model will rescue it.

Worked example — EGFR **8A2A**, a 1.43 Å entry that ranks 2nd on score, has an
8-residue gap at UniProt 868-875 whose flanks sit **7.5 Å from the inhibitor**,
with AF pLDDT 47 across it. Nothing in the resolution or coverage columns
reveals that; docking into it would be unreliable. `structure_caveats()` emits
this as a sentence automatically.

`residue_coverage()` alone returns the full per-residue table (`unp_num`,
`author_num`, `observed_ratio`, `status`) plus the `seq_offset` between author
and UniProt numbering. `plddt_by_residue(accession)` gives `{unp_num: pLDDT}`
from the AF model's B-factor column.

## No ligand? Borrow a site — with a confidence number

A predicted (AlphaFold) or apo structure has **no bound ligand, so there is no
box to compute**. The experimental path returns `pocket`/`box` only when the
chosen entry is holo. To get a docking site anyway:

```python
res = get_structure_with_site("KRAS")     # transfers a site if none was observed
res["pocket_transfer"]["confidence"]      # 'high' | 'medium' | 'low'
res["box"]                                # usable for docking
```

`transfer_pocket(target_path, accession=...)` superposes a **liganded donor**
entry onto the target and carries the ligand across, yielding lining residues, a
box, and a `_with_<lig>_from_<PDB>.cif` complex. Pass `donor_pdb_id=` to choose
the donor; otherwise the best-scoring holo entry is selected automatically.

**Read the confidence, not just the box.** Calibrated on 8 targets (KRAS, EGFR,
WRN, BRAF, CDK2, HSP90AA1, PARP1, BTK), scoring each transfer against the
donor's own crystallographic site:

| Reported confidence | n | mean recovery | min |
|---|---|---|---|
| `high` | 6 | 0.91 | 0.85 |
| `medium` | 1 | 0.56 | 0.56 |
| `low` | 1 | 0.00 | 0.00 |

Three inputs set the grade: superposition RMSD, donor sequence identity, and
site pLDDT. Correspondence between donor and target is established by sequence
alignment, so author numbering need not agree — `seq_identity` and
`n_corresponding_residues` report the alignment. Identity below 0.8 forces
`low`; below 0.5 the transfer is refused. HSP90AA1's auto-selected donor aligns
at 0.502 identity and recovers 0.00 of the site.

**Superposition RMSD is the signal; site pLDDT is not.** HSP90AA1 records site
pLDDT 89.8 with 33.8 Å RMSD and 0.00 recovery. Model confidence is independent
of whether the donor's pocket transferred correctly, so pLDDT can only
*downgrade* a grade, never raise it. It answers "is the model reliable here",
not "did the transfer work".

**Multi-domain targets get a site-local fit.** WRN's AlphaFold model has the
helicase sub-domains hinged differently from the crystal: 78% of 20-residue
windows agree under 2 Å, yet no whole-domain fit beats ~12 Å. When the global fit
exceeds 2.5 Å, the code re-fits on the 12 Å shell around the donor's ligand —
the only frame the box depends on. WRN goes 16.5 Å → 2.05 Å that way
(`site_local_fit: True`, and `superposed_on` says so).

**What this does not do.** It positions a pocket that exists *in the donor*; it
is silent about a pocket that does not. It requires a liganded relative — for
TMEM238 (no holo entry anywhere) it returns an explicit error rather than a
guess, and a genuinely novel pocket needs cavity detection or a de-novo
predictor. And it does not dock: hand `res["path"]` plus `res["box"]` to
`diffdock` or `boltz`.

## Multimers: one copy at a time

In a multimer the same ligand sits in every chain. `pocket_residues` and
`ligand_box` therefore use the **first copy** (or `chain=`) — pooling copies
would yield a box spanning both sites, which is never what docking wants. Call
`ligand_copies(path, comp_id)` to see them all.

## Function reference

- `resolve_uniprot(query, organism_id=9606, reviewed=True)` — identifier -> UniProt records. Gene symbols like `TMEM238` are not mistaken for accessions.
- `sifts_chains(accession)` — every PDB chain mapped to the accession, with per-chain resolution and coverage (SIFTS). `[]` when none.
- `rcsb_entries(pdb_ids)` — batch entry metadata in ONE GraphQL request. Never loop per id.
- `summarize_entry(entry, accession=None)` — flatten to a ranking-ready dict.
- `rank_structures(accession, ...)` — the ranked table.
- `fetch_pdb(pdb_id, fmt="cif"|"pdb")` / `fetch_alphafold(accession, with_pae=False)` — coordinates. AFDB returns `has_model=False` (not an error) when absent, plus global pLDDT and per-bin residue fractions.
- `prep_structure(...)` -> report incl. `path`, `n_residues`, `hetero_kept`, `format_coerced`.
- `pocket_residues(path, comp_id, radius=5.0, chain=None)` / `ligand_box(...)` / `ligand_copies(...)`.
- `residue_coverage(pdb_id, accession, chain=None)` — per-residue observed/missing table in UniProt numbering, with gaps classified internal / terminal / expression_tag.
- `coverage_gap_report(pdb_id, accession, ligand_comp_id=..., struct_path=...)` — gaps annotated with AF pLDDT and distance to the ligand.
- `plddt_by_residue(accession_or_path)` — `{unp_num: pLDDT}` from the AF B-factor column.
- `transfer_pocket(target_path, donor_pdb_id=None, accession=..., ...)` — borrow a binding site from a liganded donor; returns `rmsd`, `confidence`, `site_plddt`, `pocket`, `box`, `complex_path`.
- `get_structure_with_site(query, ...)` — `get_structure()` plus automatic transfer when the chosen structure has no ligand.
- `span_residues(chain, keep_residues)` / `copy_res_probe(st, comp_id, chain)` — superposition helpers.
- `get_structure(...)` / `get_pdb_entry(...)` — the two entry points.
- `structure_caveats(res)` — the warnings to report.

## Notes

- Prefer `fmt="cif"`: mmCIF has no 3-char resName limit and no 62-chain or
  99999-atom ceiling, so large assemblies survive intact.
- HTTP calls retry with backoff (EBI and RCSB time out intermittently); a 404
  is a real answer and is raised at once.
- Sequence numbering in the file is **author numbering**, which usually but not
  always matches UniProt — `unp_range` and SIFTS `start`/`end` give the mapping.
