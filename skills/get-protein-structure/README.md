# get-protein-structure

Retrieves and prepares protein 3D structures from the PDB and AlphaFold DB,
starting from a gene symbol, UniProt accession, protein name, or PDB id. A
well-studied target has hundreds of PDB entries that differ in resolution, in how
much of the sequence they observe, in whether a ligand is bound, and in what was
engineered into the construct. This skill ranks those entries, prepares the
selected one for downstream use, and reports the properties that determine whether
it is fit for purpose — bound ligands, binding-site residues, a docking box, and
which residues the entry does not observe.

## Functions

| Function | Description |
| --- | --- |
| `get_structure(query)` | identifier to prepared file: resolves the protein, ranks all mapped entries, downloads and cleans the best, and returns the pocket, docking box and coverage map. |
| `get_pdb_entry(pdb_id)` | same preparation for a specific entry, skipping resolution and ranking. |
| `rank_structures(accession)` | the ranked candidate table alone, scored on resolution, UniProt coverage, apo/holo state and engineered mutations. |
| `resolve_uniprot(query)` / `sifts_chains(accession)` / `rcsb_entries(pdb_ids)` | identifier resolution, per-chain resolution and coverage from SIFTS, and batched RCSB entry metadata. |
| `prep_structure(in_path, out_path, chains=...)` | removes waters, alt-locs, extra models and crystallisation additives; retains structural ions, cofactors and drug-like ligands. |
| `pocket_residues(path, comp_id)` / `ligand_box(...)` / `ligand_copies(...)` | binding-site residues, docking box, and the ligand copies present. |
| `residue_coverage(pdb_id, accession)` / `coverage_gap_report(...)` | per-residue observed/missing map in UniProt numbering, with gaps classified and annotated by distance to the ligand and AlphaFold confidence. |
| `fetch_alphafold(accession)` / `plddt_by_residue(...)` | predicted model with global and per-residue pLDDT. |
| `transfer_pocket(target_path, accession=...)` | assigns a binding site to a ligand-free or predicted model by superposing a liganded donor entry, with a calibrated confidence grade. |
| `get_structure_with_site(query)` | `get_structure()` plus automatic transfer when the selected structure has no bound ligand. |
| `structure_caveats(result)` | the properties of a result that affect downstream suitability, as reportable statements. |

## Run output

`get_structure_to_results(query)` is the entry point: it runs `get_structure()`
(or `get_structure_with_site()`) and writes a full run directory rather than the
bare `structures/` default the fetch functions use. The rest of the group are
the pieces it calls, exposed for a caller assembling a run by hand.

| Function | Description |
| --- | --- |
| `get_structure_to_results(query)` | one call from identifier to a written run: structure files, tables, README. |
| `results_root(root=None)` | resolves the repository `results/` directory from `$SCIENCE_RESULTS_ROOT` or an explicit `root=`; raises rather than defaulting to a cwd-relative path. |
| `gps_run_dir(target)` | the run directory for one target, `<root>/protein_structure/<slug>/`, with `scripts/` beside it. |
| `gps_write_run(out_dir, name, result)` | writes the run: structure files, ranked-candidate, pocket and coverage-gap tables, README, and any scripts. |
| `gps_run_readme(name, result, summary)` | renders the run README alone, with `structure_caveats()` supplying its Limits section. |
| `gps_write_table(path, rows)` / `gps_check_words(text, cap, label)` | CSV writer, and the word cap enforced on run-README prose. |

## Layout

`SKILL.md` holds the agent-facing instructions; `kernel.py` is a sidecar module
executed into the caller's Python kernel on skill load, so every top-level name
is exported (the loader rejects underscore-prefixed names, hence `http_get`).
Requires `gemmi`; all other dependencies are stdlib.

The tables above name the functions a caller invokes. `kernel.py` also defines
support functions those call — `http_get`, `fetch_pdb`, `summarize_entry`,
`method_family` / `method_matches`, `classify_ligand`, and the alignment
helpers `residue_correspondence`, `aa_residues`, `span_residues` and
`copy_res_probe`. They are exported because the sidecar loader exports every
top-level name, not because a caller needs them; their docstrings carry the
constraints that shaped them.

## Data sources

| Source | Used for |
|---|---|
| UniProt REST | identifier resolution |
| PDBe SIFTS `best_structures` | per-chain resolution and UniProt coverage |
| PDBe `residue_listing` + graph-api mappings | per-residue observed/missing, author-to-UniProt offset |
| RCSB GraphQL | batch entry metadata: ligands, R-free, mutations, chains |
| RCSB files | coordinates |
| AlphaFold DB | predicted models, global and per-residue pLDDT |

## Behaviour that affects callers

- Coverage is the fraction of the full UniProt sequence observed, so for a
  multi-domain target the best-resolved entry usually covers one domain while
  the full-length entries are low-resolution. Ranking blends resolution with
  coverage; `min_coverage` makes full-length an explicit request.
- `prep_structure` silently writes mmCIF when a requested `.pdb` cannot hold the
  structure — too many chains or atoms, multi-character chain ids, long CCD
  ligand ids. gemmi does not raise on any of them; it merges chains and
  overflows serials. The trigger is recorded in `format_coerced_reasons`.
- A ligand present in several chains of a multimer is several copies.
  `pocket_residues` and `ligand_box` operate on one copy; pooling them returns a
  box spanning both sites.
- Residue coordinates are keyed `(chain, residue)`, since residue numbers are
  unique only within a chain.

`SKILL.md` carries the per-parameter behaviour: prep options, coverage gaps,
multimer handling, and when to override the ranking.

## Pocket transfer confidence

`transfer_pocket` superposes a liganded donor onto a ligand-free target and
carries the ligand across. The reported grade was calibrated on eight targets
(KRAS, EGFR, WRN, BRAF, CDK2, HSP90AA1, PARP1, BTK), each transferred onto its
AlphaFold model from an auto-selected donor and scored against that donor's
crystallographic site. **That benchmark ships no script or result table with
this skill**, so the thresholds below cannot be re-derived as it stands;
rebuilding it is the prerequisite for changing any of them.

Donor-target correspondence is established by sequence alignment, not residue
numbers, since author numbering is a deposition choice the two need not share.
Superposition RMSD and donor sequence identity track recovery, while site pLDDT
does not — a model can be confident at the site and still transfer badly, which
is why pLDDT can only downgrade a grade and never raise one. Past a global RMSD
cut the fit is repeated on the shell around the donor ligand, which is the frame
the box depends on; `site_local_fit` and `superposed_on` record the result.

Transfer positions a pocket present in the donor and is uninformative about one
that is not; with no liganded relative it returns an error. A LIGSITE-style
geometric detector was evaluated and not adopted, its top-ranked cavity on apo
KRAS being a surface groove with no overlap with the switch-II site.

## Heuristics and defaults

- **Apo/holo classification.** RCSB exposes no curated "subject of
  investigation" annotation, so `classify_ligand` assigns each chem-comp id to
  `water` / `additive` / `cofactor` / `drug_like` from curated lists plus a
  molecular-weight floor. A cofactor-only entry such as KRAS with GDP is
  classified `apo` — correct for inhibitor docking, incorrect if the nucleotide
  is the ligand of interest. Classification affects ranking only.
- **Ranking weights.** Resolution 40, coverage 30, drug-like ligand +20, -3 per
  engineered mutation to a floor of four, small bonuses for X-ray and R-free
  below 0.25. An entry with no resolution scores at a 3.5 A stand-in when it is
  EM and is penalised otherwise, since NMR reports none at all. Override with
  `prefer_ligand`, `min_coverage`, `methods`.
- **Preparation is partial by design.** Waters, alt-locs, extra models and
  crystallisation additives are removed; protonation, side-chain and loop
  rebuilding, and bond-order assignment are not performed.

## Testing

Run `pytest tests/ -q` from the skill directory, using an interpreter that has
`gemmi` — the suite builds its structures with it, and `conftest.py` fails with
the environment to use when it is absent. Coverage is structure preparation,
residue mapping and pocket transfer against synthesised structures, so no PDB
or AlphaFold network access is needed. Ranking, format coercion and donor
selection depend on live API responses and are not covered; re-run a known
target end to end after changing them.

## Scope

Structure prediction (`boltz`, `chai1`, `alphafold2`, `openfold3`), docking
(`diffdock`, `boltz`), and protonation, loop rebuilding or bond-order assignment
(PDBFixer, Protoss) are out of scope and handled by those tools.
