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

- `get_structure(query)` — identifier to prepared file: resolves the protein,
  ranks all mapped entries, downloads and cleans the best, and returns the pocket,
  docking box and coverage map.
- `get_pdb_entry(pdb_id)` — same preparation for a specific entry, skipping
  resolution and ranking.
- `rank_structures(accession)` — the ranked candidate table alone, scored on
  resolution, UniProt coverage, apo/holo state and engineered mutations.
- `resolve_uniprot(query)` / `sifts_chains(accession)` / `rcsb_entries(pdb_ids)` —
  identifier resolution, per-chain resolution and coverage from SIFTS, and batched
  RCSB entry metadata.
- `prep_structure(in_path, out_path, chains=...)` — removes waters, alt-locs,
  extra models and crystallisation additives; retains structural ions, cofactors
  and drug-like ligands.
- `pocket_residues(path, comp_id)` / `ligand_box(...)` / `ligand_copies(...)` —
  binding-site residues, docking box, and the ligand copies present.
- `residue_coverage(pdb_id, accession)` / `coverage_gap_report(...)` — per-residue
  observed/missing map in UniProt numbering, with gaps classified and annotated by
  distance to the ligand and AlphaFold confidence.
- `fetch_alphafold(accession)` / `plddt_by_residue(...)` — predicted model with
  global and per-residue pLDDT.
- `transfer_pocket(target_path, accession=...)` — assigns a binding site to a
  ligand-free or predicted model by superposing a liganded donor entry, with a
  calibrated confidence grade.
- `get_structure_with_site(query)` — `get_structure()` plus automatic transfer when
  the selected structure has no bound ligand.
- `structure_caveats(result)` — the properties of a result that affect downstream
  suitability, as reportable statements.

### Saving a run

- `get_structure_to_results(query, name=None, root="results",
  topic="protein_structure", summary=None, with_site=False, ...)` —
  `get_structure()` (or `get_structure_with_site()` with `with_site=True`),
  routed into `results/protein_structure/<slug>/` instead of the bare
  `structures/` default `get_structure` otherwise writes to, then saved as a
  full run via `gps_write_run`. Returns the usual result dict with a `run`
  key added.
- `gps_run_dir(target, root="results", topic="protein_structure")` — the path
  for one run, `<root>/<topic>/<slug>/`, with `scripts/` created beside it.
- `gps_write_run(out_dir, name, result, summary=None, files=(),
  data_sources=(), limits=(), scripts=())` — copies the structure file(s)
  already on disk at `result["path"]`/`result["raw_path"]` into `out_dir`,
  writes `ranked_structures.csv`, `pocket_residues.csv` (when a pocket was
  found) and `coverage_gaps.csv` (when computed), the run README, and copies
  `scripts` into `scripts/`. Returns the paths written. Raises if `result`
  carries an `error`, before writing anything.
- `gps_run_readme(name, result, summary, files=(), data_sources=(),
  limits=(), title=None)` — renders the README text `gps_write_run` saves,
  per `coding-standards`' Result/Files/Data sources/Limits structure. The
  Limits section is `structure_caveats(result)` verbatim plus any
  run-specific ones passed in — reusing the function this skill already
  computes rather than restating its checks by hand.
- `gps_write_table(path, rows, headers=None)` — CSV writer used by
  `gps_write_run`, usable standalone.

Runs land in `results/protein_structure/<slug>/`, the topic already used by
the hand-written WRN and KRAS dossiers, so a new run is a sibling.

## Scope

Structure prediction (`boltz`, `chai1`, `alphafold2`, `openfold3`), docking
(`diffdock`, `boltz`), and protonation, loop rebuilding or bond-order assignment
(PDBFixer, Protoss) are out of scope and handled by those tools.

`SKILL.md` is the agent-facing guidance. This file documents the implementation
for a human reading the code.

## Layout

`SKILL.md` holds the agent-facing instructions; `kernel.py` is a sidecar module
executed into the caller's Python kernel on skill load, so every top-level name
is exported (the loader rejects underscore-prefixed names, hence `http_get`).
Requires `gemmi`; all other dependencies are stdlib.

## Data sources

| Source | Used for |
|---|---|
| UniProt REST | identifier resolution |
| PDBe SIFTS `best_structures` | per-chain resolution and UniProt coverage |
| PDBe `residue_listing` + graph-api mappings | per-residue observed/missing, author-to-UniProt offset |
| RCSB GraphQL | batch entry metadata: ligands, R-free, mutations, chains |
| RCSB files | coordinates |
| AlphaFold DB | predicted models, global and per-residue pLDDT |

RCSB metadata is fetched in one batched GraphQL request per 50 entries. Ranking
60 candidates costs two HTTP calls.

## Implementation notes

Each condition below yields a valid-looking file rather than an error, so the
code checks for it explicitly. Removing a check will not cause a test failure.

- **Ranking cannot be driven by coverage alone.** EGFR has 531 chains at
  coverage 0.28 (kinase domain, to 1.07 A) and 12 at 1.00 (full-length cryo-EM,
  3.1-3.6 A). The shortlist pre-filter blends resolution with coverage;
  `min_coverage` makes full-length an explicit request.
- **Legacy PDB format loses data without raising.** Writing the 89-chain
  ribosome 4V6X through `make_pdb_string()` yields 53 distinct chain characters
  and atom serials overflowing to `A2YB1`. `prep_structure` writes mmCIF instead
  above 62 chains, 99,999 atoms, multi-character chain ids, or a 5-character CCD
  ligand id (issued since 2023, e.g. `A1BEA`), recording the reason in
  `format_coerced_reasons`.
- **A ligand in a multimer is several ligands.** `A1BEA` occupies chains A and B
  of 9E3S; pooling the copies gives a box spanning both sites.
  `pocket_residues` and `ligand_box` operate on one copy.
- **An unresolved loop can line the pocket.** EGFR 8A2A is 1.43 A and ranks
  second, yet residues 868-875 are unresolved with their flanks 7.5 A from the
  inhibitor and AlphaFold pLDDT 47 across the gap.
- **One accession can map to a chain in several segments.** 2RH1 maps as 1-230
  and 264-365, each with its own offset, with `None` at boundaries whose residue
  is disordered. Every segment is retained, each residue resolves through the
  segment containing it, and residues the accession does not claim (168 lysozyme
  residues) are reported as `n_unmapped_residues`.
- **Residue numbers are unique only within a chain.** Coordinates are keyed
  `(chain, residue)`; with no chain requested, the chain of the selected ligand
  copy is used.
- **`keep_ions` outranks `keep_ligands=False`.** A catalytic metal is part of the
  functional state.
- **SIFTS and RCSB use different method vocabularies.** SIFTS returns
  `X-ray diffraction` / `Electron Microscopy`; RCSB returns `X-ray` / `EM`.
  `method_family` normalises both so the `methods=` filter and the scoring
  bonuses cannot drift apart. An entry with no resolution is scored at a 3.5 A
  stand-in when it is EM and penalised otherwise, since NMR has no resolution
  at all and would otherwise crowd out high-resolution X-ray entries.
- **`prep_structure` validates its selectors.** An out-of-range `model_index` or
  an absent `chains=` entry raises rather than substituting a different model or
  chain, and the `hetero_kept` inventory is reconciled after chain filtering so
  it never lists a ligand absent from the written file.
- **`fmt` applies to the prepared output, not only the download.** The extension
  follows the file `fetch_pdb` returned, since RCSB serves no legacy `.pdb` for
  oversized entries.

## Pocket transfer confidence

`transfer_pocket` superposes a liganded donor onto a ligand-free target and
carries the ligand across. Confidence is calibrated on 8 targets (KRAS, EGFR,
WRN, BRAF, CDK2, HSP90AA1, PARP1, BTK), each transferred onto its AlphaFold
model from an auto-selected donor and scored against that donor's
crystallographic site.

| Grade | n | mean recovery | min |
|---|---|---|---|
| `high` | 6 | 0.91 | 0.85 |
| `medium` | 1 | 0.56 | 0.56 |
| `low` | 1 | 0.00 | 0.00 |

Donor-target correspondence is established by sequence alignment, not residue
numbers, since author numbering is a deposition choice the two need not share.
Superposition RMSD and donor identity predict recovery; site pLDDT does not —
HSP90AA1 records site pLDDT 89.8 with 33.8 A RMSD, 0.502 donor identity and
0.00 recovery. Identity below 0.8 forces `low`; below 0.5 the transfer is
refused. pLDDT therefore only
downgrades a grade. Re-run the benchmark before changing these thresholds.

Above 2.5 A global RMSD the fit is repeated on the 12 A shell around the donor
ligand, which is the frame the box depends on: WRN's full-length model goes
16.5 A to 2.05 A, recovery 6/18 to 10/18, recorded in `site_local_fit` and
`superposed_on`.

Transfer positions a pocket present in the donor and is uninformative about one
that is not. With no liganded relative (TMEM238) it returns an error. A
LIGSITE-style geometric detector was evaluated and not adopted: on apo KRAS its
top-ranked cavity is a 101-residue surface groove with no overlap with the
switch-II site.

## Heuristics and defaults

- **Apo/holo classification.** RCSB exposes no curated "subject of
  investigation" annotation, so `classify_ligand` assigns each chem-comp id to
  `water` / `additive` / `cofactor` / `drug_like` from curated lists plus a
  molecular-weight floor. A cofactor-only entry such as KRAS with GDP is
  classified `apo` — correct for inhibitor docking, incorrect if the nucleotide
  is the ligand of interest. Classification affects ranking only.
- **Ranking weights.** Resolution 40, coverage 30, drug-like ligand +20, -3 per
  engineered mutation, small bonuses for X-ray and R-free below 0.25. Override
  with `prefer_ligand`, `min_coverage`, `methods`.
- **Preparation is partial by design.** Waters, alt-locs, extra models and
  crystallisation additives are removed; protonation, side-chain and loop
  rebuilding, and bond-order assignment are not performed.

## Worked examples

Reference outputs are in `../../protein_structure/` (KRAS, WRN).

| Target | Result | Exercises |
|---|---|---|
| KRAS | 9IAY, 0.95 A, switch-II inhibitor | Clean case, no caveats |
| EGFR | 8A27, 1.07 A kinase domain | Multi-domain ranking |
| EGFR 8A2A | 8-residue gap, 7.5 A from ligand | Near-site gap detection |
| WRN | 10AK, 1.37 A helicase core | Fragmented target, 33 inhibitors |
| WRN full-length | AF model + site from 10AK, 2.05 A | Site-local fit |
| TMEM238 | AlphaFold fallback, pLDDT 67 | No structure, no donor |
| 4V6X | 89 chains preserved as mmCIF | Format ceiling |
| 2RH1 | 2 segments, 168 unmapped residues | Fusion construct |
| 1CA2 | Zn retained with `keep_ligands=False` | Ion precedence |
| HSP90AA1 | donor identity 0.502, RMSD 33.8 A, recovery 0.00 | Confidence calibration |

## Testing

`tests/` (run with `pytest tests/ -q` from the skill directory) covers the
pure-logic paths with synthesised structures — chain-keyed ligand-distance
lookups, `keep_ions` precedence over `keep_ligands=False`, SIFTS multi-segment
offset resolution, and expression-tag exclusion from `n_missing` — so it needs
no PDB/AlphaFold network access. The worked-example table above exercises the
network-dependent behaviours (ranking, format coercion, pocket transfer) that
the suite can't synthesize; re-run those identifiers after changes to ranking,
preparation, or residue mapping. Re-run the 8-target benchmark after changing
`transfer_pocket` scoring and
confirm recovery still orders `high` > `medium` > `low`.
