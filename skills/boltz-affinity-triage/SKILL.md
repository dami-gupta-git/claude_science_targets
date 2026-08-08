---
name: boltz-affinity-triage
description: Rank candidate compounds against a protein target with Boltz-2 (or Chai-1) co-folding, calibrated against known actives and property-matched decoys pulled from ChEMBL, BindingDB and PubChem BioAssay. Use this whenever someone wants to score, rank, triage, prioritise or virtually screen small molecules against a target, asks for a binding affinity or Kd/Ki/IC50 estimate from a structure model, asks whether a compound "will bind" or "how tightly", or wants to sanity-check co-folding affinity output — even if they only name Boltz, Chai, DiffDock or "AI docking" and do not mention controls, calibration or knowledge bases. Boltz-2's affinity numbers are not calibrated potencies and its training set already contains most public bioactivity data, so this skill exists to make the ranking interpretable and to say when it is not.
---

# KB-grounded affinity triage

Boltz-2 will return an affinity number for any protein–ligand pair you hand it.
The number always looks plausible, which is the problem: on its own you cannot
tell a genuine 10 nM prediction from the model's output on a compound that
cannot bind at all. This skill wraps the co-folder in the two things that make
its output mean something — **internal controls** drawn from bioactivity
knowledge bases, and a **gate** that refuses to emit a ranking when the controls
show the run has no discriminative power on that target.

Use `boltz` (or `chai1`) for the mechanics of running the model. This skill is
the experimental design around it.

## Why controls, and why property-matched ones

Three facts about the affinity head drive the whole workflow:

1. **The output is not a potency.** `affinity_pred_value` is log10(IC50 in µM)
   on an assay-agnostic scale. Boltz-2 is trained on affinity records
   irrespective of activity type, so its own authors describe the prediction as
   a generic affinity estimate rather than a value corresponding to any single
   assay type — treat it as an ordering, never as a number to quote. Rank with
   `affinity_probability_binary`, which is what the authors intend for hit
   discovery.
2. **Public bioactivity data is mostly *in* the training set.** The affinity
   head was trained on continuous affinity data curated primarily from ChEMBL
   v34 and BindingDB, plus PubChem, as of June 2023. A control pulled from
   ChEMBL may be recalled rather
   than predicted, which inflates your apparent accuracy. This does not make
   controls useless — it means you must report them stratified by whether they
   plausibly predate the cutoff, and read a strong result on old compounds as a
   sanity check ("the model has not broken") rather than as evidence of
   prospective skill.
3. **Naive decoys make any scorer look good.** This is the trap that quietly
   invalidates most casual benchmarks. Potent actives are *bigger* than
   inactives: on EGFR, ChEMBL's sub-100 nM actives average 30.8 heavy atoms
   while a default weak-activity pull averages 16.0. A scorer using nothing but
   molecular size scores AUC **0.97** against that unmatched set. Property-match
   the decoys and the same size-only scorer drops to **0.49**, exactly where it
   belongs.

So: match your decoys on size and mass, and always report the size-only null
AUC beside the real one. `kba_null_auc` computes it; `kba_gate` enforces the
margin.

## Workflow

Helpers below ship in `kernel.py` and are already defined when this skill
loads. MCP calls go in the `repl` tool; analysis in `python`.

### 1. Resolve the target

Get a ChEMBL target id and a UniProt accession — you need both, because ChEMBL
keys on its own id and BindingDB keys on UniProt.

```python
# repl
t = host.mcp("chembl", "target_search", gene_symbol="EGFR",
             organism="Homo sapiens", target_type="SINGLE PROTEIN", limit=5)
# -> target_chembl_id CHEMBL203, components[0].accession P00533
```

Take the sequence for the YAML from UniProt (`mcp-genes-ontologies`) or from a
structure via `get-protein-structure`. One `protein` chain per YAML, and each
chain needs an MSA — pass `--use_msa_server` unless you already have an `.a3m`.

### 2. Pull actives from more than one knowledge base

Independent KBs matter here for a specific reason: agreement between ChEMBL and
BindingDB on a compound is evidence the measurement is real rather than a single
unreplicated assay, and BindingDB emphasises Ki/Kd from medicinal chemistry
papers where ChEMBL leans on curated IC50. Query both.

```python
# repl
act = host.mcp("chembl", "get_bioactivity", target_chembl_id="CHEMBL203",
               activity_type="IC50", min_pchembl=7, limit=200)
bdb = host.mcp("chemistry", "bindingdb_ligands_by_target",
               uniprot="P00533", affinity_cutoff_nm=100, max_rows=300)
```

Two retrieval hazards, both of which silently corrupt results:

- **BindingDB `affinity` is a string** that may carry `>` or `<`. A `>10000` is
  a genuine inactive; a `<1` is a censored ceiling whose true value is unknown.
  Parse with `kba_parse_nm`, which returns `(value, qualifier)` — never `float()`
  it directly. Boltz-2's own curation treats `>`-qualified values as decoys for
  hit discovery and as censored measurements for hit-to-lead, so preserving the
  qualifier keeps your controls consistent with how the model was trained.
- **BindingDB rows sort by (affinity_type, affinity)**, so a small `max_rows`
  can return one measurement type only. A 50-row EGFR pull came back 50/50
  EC50, zero Ki. Raise `max_rows` and check the type mix before trusting it.

For ChEMBL, prefer `assay_type == "B"` (direct binding) and drop rows carrying
a `data_validity_comment`. Use `kba_pactivity(row, source)` for a unified
pActivity — it falls back to `standard_value` when `pchembl_value` is absent,
which is not a nicety: 38 of 50 weak EGFR rows had no `pchembl_value`, so a
pipeline reading that field alone loses most of its negatives and then reports
a suspiciously clean AUC.

### 3. Build a property-matched decoy set

Pull the weak/inactive tail wide, then match. The pool must be large and
size-diverse or matching cannot work — a narrow pull has no big molecules to
pair with big actives, and the gap persists.

```python
# repl -- ask for a wide pool, not the default 20
weak = host.mcp("chembl", "get_bioactivity", target_chembl_id="CHEMBL203",
                activity_type="IC50", min_value=10000, unit="nM", limit=1000)
```

```python
# python
decoys = kba_match_decoys(actives, pool, per_active=1)
null = kba_null_auc(actives, decoys)   # must land near 0.5
```

Check `null` before spending GPU time. Near 0.5 means the controls are
size-balanced. Still high means your pool is too narrow — widen it, or loosen
`tol`, and if it will not balance, say so in the report rather than proceeding.

Presumed-inactive decoys from an unrelated-target pool are acceptable when a
measured weak set is unavailable, but they carry a real false-negative rate
(some will bind), which caps the AUC you can achieve. Prefer measured
inactives on the same target whenever they exist.

### 4. Co-fold controls and queries in one batch

Controls and queries must run under identical settings, or the comparison is
between configurations rather than compounds. Write one YAML per compound with
`kba_write_yaml`, which attaches the `properties: affinity:` block and skips
ligands at or over Boltz's 128-atom affinity cap instead of silently returning
structure-only output.

```python
# python
for cid, smi in compounds:
    kba_write_yaml(f"yaml/{cid}.yaml", target_seq, smi, ligand_id="L")
```

```bash
boltz predict yaml/ --use_msa_server --out_dir out/ \
      --diffusion_samples 5 --recycling_steps 3
```

Read results with `kba_read_result(pred_dir)`. It discovers score keys rather
than assuming them, because Boltz-2.1 replaced `affinity_pred_value` /
`affinity_probability_binary` with `binding_confidence` /
`optimization_score` — code that hardcodes the old names breaks silently on
upgrade. It normalises to `binder_score`, higher = more likely a binder.

`binding_confidence` is the Boltz-2.1 analogue of
`affinity_probability_binary` and ranks directly. `optimization_score` does
**not**: its sign convention is undocumented, and the published descriptions of
the analogous Boltz-2 field contradict each other — the model docs define
`affinity_pred_value` as log10(IC50 in µM), where lower binds tighter, while
third-party summaries describe the same field as a pIC50, where higher binds
tighter. Those are negatives of one another, so guessing would rank a screen
backwards while still producing a plausible-looking table. An output carrying
only `optimization_score` therefore comes back with `binder_score=None` and an
`unscored_reason`; the raw value is still reported. Once you have confirmed the
direction on a known potent/weak pair for your target, pass
`optimization_score_direction="higher_is_better"` (or `"lower_is_better"`) to
enable ranking on it.

Confidence is **per diffusion sample**, and this matters as soon as you pass
`--diffusion_samples`. Boltz writes one `confidence_<name>_model_K.json` per
sample and ranks them so `model_0` is its best; affinity is predicted once per
input. `kba_read_result` therefore reads confidence from a named sample
(`model_index=0` by default) rather than from whichever file the filesystem
happens to list first — reading an arbitrary sample can report an ipTM tens of
points away from the ranked pose, which inverts the check below. The returned
`per_sample` dict carries every sample's confidence if you want to see the
spread, and `model_index` records which one was used.

Also read the structural confidence: an affinity prediction sitting on a bad
pose is not trustworthy. Flag `iptm < 0.5` rows (use `ptm` for single-chain).
The published evaluations found affinity sometimes survives a wrong pose, so
treat low ipTM as a caveat to report rather than an automatic discard — but
never present it as a clean hit.

### 5. Gate, then rank

```python
# python
g = kba_gate(auc, n_pos, n_neg, null_auc=null)
```

If `g["pass"]` is False, the deliverable is the diagnostic — the AUC, the null,
the control counts, and what to try next. Do not also ship a ranked hit list;
a ranking from an uncalibrated run is noise wearing the costume of a result,
and it is worse than no answer because someone will order compounds from it.

When it passes, report for each query compound its `binder_score`, its
**percentile within the control distribution** (`kba_percentile`) — far more
interpretable than the raw number — its ipTM, and its nearest known active by
2D similarity, so the reader can see whether a hit is a genuine novel chemotype
or a near-analog of a training-set compound.

## Report format

Lead with whether the run is interpretable, because everything else is
conditional on it:

```
## Calibration
Target, controls (n actives / n decoys, KB sources), control AUC, size-only
null AUC, EF@10%, verdict from kba_gate. Leakage split: n likely_train /
n likely_novel, with the AUC on each if both are populated.

## Ranked queries        <-- only when the gate passes
compound | binder_score | control percentile | ipTM | nearest known active (Tanimoto) | leakage flag

## Caveats
Boltz-2 affinities are relative, not potencies. Name the compounds over the
128-atom cap, the low-ipTM rows, and any decoys that are presumed rather than
measured inactives.
```

## Things that will bite you

| Symptom | Cause / fix |
|---|---|
| Control AUC ~0.95 on the first try | Almost always unmatched decoys. Check `kba_null_auc`; if it is also high the separation is size, not binding. |
| Every active flagged `likely_train` | Normal for a well-studied target — a default ChEMBL pull on EGFR returned 50/50 papers from 1997–2003. Report it; do not claim prospective accuracy. |
| Most negatives vanish after parsing | Reading `pchembl_value` only. Use `kba_pactivity`. |
| No `affinity_*.json` in output | FASTA input, or the YAML lacks `properties:`. Affinity needs YAML plus a `ligand` binder chain. |
| `KeyError: 'iptm'` | Single-chain complex — read `ptm`. |
| BindingDB rows all one assay type | `max_rows` cap interacting with its sort order. Raise it. |
| No GPU available | Boltz-2 needs one. Check `list_compute`; if empty, tell the user and offer Colab or a smaller control set rather than starting a run that cannot finish. |

## Helpers in `kernel.py`

`kba_parse_nm` · `kba_pactivity` · `kba_heavy_atoms` · `kba_props` ·
`kba_match_decoys` · `kba_null_auc` · `kba_leakage_flag` · `kba_write_yaml` ·
`kba_read_result` · `kba_auc` · `kba_enrichment` · `kba_percentile` ·
`kba_gate`

RDKit gives exact atom counts and MW; without it `kba_heavy_atoms` falls back
to a regex estimate (mean error 0.36 atoms, max 2, on 100 ChEMBL structures)
and MW-based matching is skipped. Install rdkit when matching quality matters.
