---
name: target-triage-public-data
description: Assess whether a gene is a plausible drug target using public data — DepMap/Open Targets genetic dependency, tractability, expression-vs-drug-sensitivity in cell lines (GDSC, PRISM Repurposing), and TCGA survival/subtype correlates. Use when asked "is X a drug target for Y cancer", to sanity-check a target from the literature, or to test whether a gene's expression predicts response to a specific drug class. Covers the proliferation confounder that makes naive sensitivity correlations wrong.
---

# Public-data target triage

Answer "is this gene a drug target in this indication?" from public data, in an
order that spends the least effort before the cheapest disconfirming evidence.

**What this skill is not.** It does not decide whether a target is good — that is
judgment about mechanism, competition, and indication, and it belongs to you and
the user. This skill gets the evidence on the table correctly, and names the two
places where the obvious analysis gives the wrong answer.

## Run cheap disconfirming evidence first

Order matters. Step 1 is one API call and frequently settles the single-agent
question before you invest in anything else.

1. **Genetic dependency** — do cells in the indication need the gene?
2. **Tractability** — is it druggable, and does chemistry already exist?
3. **Expression vs drug sensitivity** — does expression predict response? (the
   step with the confounder)
4. **Clinical correlate** — does it stratify patients, and specifically
   treated patients?

Skip steps once an earlier one has answered the question being asked. A gene
with a flat dependency profile is not a single-agent target, and no amount of
prognostic correlation changes that.

## Step 1 — Genetic dependency

**If a local DepMap release is available, prefer the `depmap-local` skill** — it reads the
release directly, covers every screened line rather than a disease subset, and ships a
column-addressable cache:

```python
skill("depmap-local")
depmap_selectivity("DCTPP1")          # common-essential / selective / non-essential
depmap_lineage_enrichment("DCTPP1")   # which lineages, BH-adjusted
```

Without local files, Open Targets carries the same DepMap scores via GraphQL:

```python
t = ot_target("ENSG00000179958")             # Ensembl gene id
dep = ot_essentiality_frame(t)               # gene, tissue, cell, disease, effect, expr
crc = dep[dep.disease.str.contains("Colorectal|Colon|Rectal", case=False, na=False)]
print(crc.effect.median(), crc.effect.min(), len(crc))
```

The two agree — on one worked target the colorectal subsets matched at mean −0.095 and minimum
−0.433 from either source — but Open Targets returns only lines annotated to a disease, so the
local release is the one that supports a pan-cancer statement.

CRISPR gene effect (Chronos) is scaled so **0 = no effect and −1 = the median common-essential
gene**. A median near 0 means cells tolerate losing it.

**Always print comparators from the same cell lines** — an absolute number is uninterpretable
alone. Good nucleotide-metabolism reference points: RRM1 and DUT run around −2 to −3 (strongly
essential), KRAS around −1.2 in CRC (lineage dependency), NUDT1 around 0 (dispensable despite
being a drug target of interest). Fetch them the same way and report side by side.

**A stratified contrast needs a null before you believe it.** A genotype- or
subtype-restricted dependency — MSI vs MSS, BRCA1-mutant vs wild-type, deleted vs
intact — returns an uncorrected p-value for a hypothesis drawn from thousands of
candidate markers. Load `marker-contrast-null` and rank the chosen marker against
the same contrast over every other eligible marker. This is the step-1 analogue of
the proliferation confounder at step 3, and it has overturned a triage here: USP1's
BRCA1-mutant dependency (p = 0.007) ranks 67 of 1,719 markers with zero surviving
correction, which is where ~68 markers were always going to land.

Copy-number markers need the null built from deletions, not mutations: a gene lost
by deletion is absent from a damaging-mutation matrix entirely. Score the null over
every gene's deletion call, then read the markers above yours — for MTAP the eleven
stronger markers were all its own 9p21 neighbours, carried by the same deletion,
which is a different finding from eleven unrelated genes beating it.

**Check pan-cancer before concluding "not a target in indication X".** A gene that is flat
everywhere is a different finding from one that is flat in your indication but essential
elsewhere — and it extends the caution to every other indication's literature.

**Also test expression vs dependency.** If high expressors are no more dependent than low, the
prognostic literature's premise — that high-expressing tumours need the gene — is not supported.

## Step 2 — Tractability

`ot_target()` already returned it. Labels worth reporting: `Small Molecule
Binder`, `High-Quality Ligand`, `Structure with Ligand`, and the antibody/PROTAC
modality rows. This is usually the least surprising step; state it briefly and
move on.

## Step 3 — Expression vs drug sensitivity

### The confounder that will otherwise mislead you

Genes correlated with proliferation rate correlate with sensitivity to *every*
cytotoxic drug. A per-drug correlation therefore looks specific when it is not.

**Run `confounder_check` before believing any single-drug hit**, and report
`r_partial` (from `sensitivity_scan`) rather than `r`:

```python
print(confounder_check(bowel, "expr", drug_cols))   # r across the whole panel
res = sensitivity_scan(bowel, "expr", drug_cols)    # per-drug, raw + partial + BH
# The confound is a LEAVE-ONE-OUT panel mean: each drug is excluded from its
# own confound. Including it regresses part of a drug's signal out of itself
# and attenuates r_partial toward zero by roughly 1/len(drug_cols) — ~46% on
# a 3-drug panel, ~10% at 25. The `confound` column records which was used;
# pass leave_one_out=False to reproduce a whole-panel run.
```

A worked case — DCTPP1 in colorectal lines. **In GDSC** (40 COREAD lines, lnIC50):
raw 5-FU r = −0.27 (p = 0.10) and oxaliplatin r = −0.36 (p = 0.02), but panel-wide
r = −0.41 (p = 0.009), and after partialling 5-FU went to r = +0.04 (p = 0.80).
The apparent hit was entirely general sensitivity. **In PRISM 23Q2** (33 bowel lines,
single-dose LFC) the same drug gave raw r = −0.08 and partial r = +0.09 (p = 0.61).

Always state which screen a number came from. The two panels differ in cell lines,
readout (lnIC50 vs single-dose log-fold-change) and n, so their correlations are not
interchangeable even for the same gene and drug — quoting one as if it were the other
is an easy and invisible error.

Carried to completion, that target was tested across four screens — GDSC lnIC50, PRISM 23Q2 and
24Q2 single-dose, and the 20Q2 dose-response AUC panel — for 63 drug × screen combinations with
zero significant after correction. The dose-response screen is what retires the "single-dose
assay might miss a potency shift" objection, so run it whenever the compound is present.

Two more habits: correct across drugs tested (`sensitivity_scan` returns BH
q-values), and **read the sign of the whole scan**. If a resistance hypothesis
predicts positive correlations and every drug trends negative, that is
informative even when nothing is individually significant.

### Getting the data

**GDSC** (fitted dose-response, gives lnIC50 and AUC) — public, direct download:

```
https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.4/
  GDSC1_fitted_dose_response_24Jul22.csv     (~54 MB)
  GDSC2_fitted_dose_response_24Jul22.csv     (~40 MB)
  screened_compounds_rel_8.4.csv             (drug id -> name, target, pathway)
```

Use `curl`/`bash` rather than `urllib.request.urlretrieve` — the latter has been observed
stalling for 10+ minutes on this host where curl finishes in ~15 s. Filter by `TCGA_DESC` (e.g.
`COREAD`) and average replicate `DRUG_NAME` rows.

GDSC's oncology panel is broad but **lacks trifluridine, decitabine, azacitidine, floxuridine
and capecitabine** — check coverage before concluding a compound cannot be tested.

**PRISM Repurposing** covers those compounds. DepMap's portal sits behind a bot-verification
wall, so scripted download does not work — **ask the user to download it**; do not attempt to
circumvent the check. There are three distinct products and choosing the wrong one costs a
round trip:

| Release | Readout | Coverage | Use when |
|---|---|---|---|
| `Repurposing_Public_24Q2_Extended_Primary_*` | single-dose LFC (~2.5 uM) | ~6.5k compounds x 906 lines | Widest compound coverage; current primary |
| `Repurposing_Public_23Q2_Extended_Primary_*` | single-dose LFC | slightly fewer compounds | Superseded by 24Q2 |
| `prism-repurposing-20q2-secondary-screen-dose-response-curve-parameters.csv` | **fitted AUC / IC50 / EC50** | 1,489 compounds x ~500 lines | **Preferred when your compound is present** |

The **secondary dose-response file lives on the older "PRISM Repurposing 19Q4/20Q2" release
page, not on the 23Q2/24Q2 page** — those contain primary single-dose data only, whatever the
version number suggests. Only primary-screen hits were carried into the secondary screen, so
check membership before promising the analysis:

```python
sec = pd.read_csv(path, usecols=["depmap_id","name","auc","ic50","r2"], low_memory=False)
print(sorted(set(sec.name.str.lower()) & {"trifluridine","decitabine"}))
```

Column `name` is the compound (lowercase it), `depmap_id` is the ACH- id, `auc` is the fitted
area under the dose-response curve — **higher AUC = more resistant**, the opposite polarity to
the primary screen's log-fold-change. Average replicate rows per (line, compound) before
correlating. Prefer AUC over single-dose LFC: it is also less proliferation-confounded (on one
worked target, panel-wide r = −0.13 for AUC vs −0.41 for GDSC lnIC50).

For the single-dose matrices, `load_prism()` handles the BRD-id-to-name mapping:

```python
PR = load_prism(matrix_path, compound_list_path, drug_names=["TRIFLURIDINE", "DECITABINE"])
```

**Expression** — do NOT download DepMap's 500 MB expression matrix if cBioPortal will do. Its
CCLE study carries `DEPMAP_ID` alongside expression, so PRISM joins directly on ACH- ids:

```python
# repl tool (MCP/HTTP), then handoff via ./handoff/*.json
#   POST {CBIO_API}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL  body: ["DCTPP1"]
#   POST {CBIO_API}/molecular-profiles/ccle_broad_2025_rna_seq_mrna/molecular-data/fetch
#        body: {"entrezGeneIds":[...], "sampleListId":"ccle_broad_2025_all"}
#   GET  {CBIO_API}/studies/ccle_broad_2025/clinical-data?clinicalDataType=SAMPLE
#        -> keep DEPMAP_ID, ONCOTREE_LINEAGE, CANCER_TYPE_DETAILED, CELL_LINE_NAME
```

If the local release IS present, use `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`
instead — it keys on `ModelID` directly, avoiding name-matching entirely. The file is wide
(~19k gene columns as `SYMBOL (ENTREZ)`); read the one column you need with `csv.reader` rather
than loading the whole frame. Running the analysis through both sources is a cheap replication
that catches join errors.

Lineage lives in `ONCOTREE_LINEAGE` / `OncotreeLineage` (colorectal = `Bowel`). For GDSC, whose
ids are names not ACH-, join with `normalize_cell_name` on both sides and print how many lines
failed to match.

**Run every screen you have access to.** Agreement across independent panels with different
readouts is far more convincing than one null, and disagreement tells you the result is fragile.

**Power.** n is typically 30-45 lines per lineage, which powers only |r| >~ 0.45.
Say so when reporting a null.

## Step 4 — Clinical correlate

cBioPortal for expression and outcome; **GDC for treatment**, which cBioPortal
does not carry:

```
POST https://api.gdc.cancer.gov/cases
  filters: cases.project.project_id in ["TCGA-COAD","TCGA-READ"]
  expand:  diagnoses.treatments
  -> therapeutic_agents, treatment_type, treatment_or_therapy
```

The interesting split is **treated vs untreated**, not high vs low overall: a
survival difference confined to the treated arm is evidence the gene conditions
response. Fit the formal interaction term (`lifelines` CoxPHFitter with
`expression x treated`) rather than eyeballing two Kaplan-Meier curves — with
typical TCGA event counts (30-40 deaths in a treated arm) both arms are usually
underpowered, and the interaction p is the honest summary.

Report subtype checks (MSI/MSS, amplification status) with an **effect size**,
not just p. Cliff's delta is appropriate; a significant p on a 0.16 log2 median
difference is not a stratification hypothesis.

## Literature sweep (optional front end)

When the ask starts from "what does the literature say":

- **Expand synonyms** — genes carry legacy names that PubMed does not unify
  (DCTPP1 is also dCTPase, XTP3TPA, RS21-C6). Search the union; a synonym-only
  query typically surfaces the foundational enzymology the modern name misses.
- **Reconcile every batched fetch.** A search reporting 37 hits followed by a
  batched metadata call can silently yield 28 records, and the resulting counts
  look self-consistent. Always:

```python
missing = reconcile_fetch(search_pmids, [a["identifiers"]["pmid"] for a in arts])
```

- Tag records by evidence class (bioinformatic / functional / chemistry) to show
  what kind of evidence a field actually rests on — write the regexes per target
  rather than reusing another gene's.

## Reporting

Every run writes a `README.md` into its results directory. Build it with
`write_triage_readme()` rather than by hand, so runs stay comparable:

```python
write_triage_readme(
    f"{out_dir}/README.md", gene="WRN",
    summary="Plain prose a non-specialist can follow. No statistics, no bullets.",
    steps=[
        {"name": "Dependency", "finding": "...", "table": strata_df},
        {"name": "Tractability", "finding": "..."},
        {"name": "Sensitivity", "skipped": "No inhibitor in the screens on disk."},
        {"name": "Clinical", "finding": "..."},
    ],
    files=[("wrn_triage.png", "the figure")],
    data_sources=["DepMap 24Q2 CRISPRGeneEffect.csv"],
    limits=["Knockout is not pharmacological inhibition."])
```

The opening `summary` is the part a non-specialist reads: one plain paragraph,
capped at 130 words, no statistics — those go in the step findings, each capped
at 90. The caps are deliberate. A run README is an orientation document, not a
second copy of the analysis; the CSVs carry the detail. Raise the constants in
`kernel.py` if a whole class of runs needs more room, rather than splitting text
across fields to evade the cap.

A step that could not be run is reported with `skipped`, never omitted — "not
runnable on these data" and "not looked at" are different claims, and an absent
section reads as the second.

- Lead with the dependency number and its comparators.
- Report partial correlations; mention raw ones only to explain why they were
  discarded.
- State what the nulls do **not** cover. Baseline expression not predicting
  response is compatible with pharmacological inhibition still sensitising —
  no public dataset tests co-treatment, and mRNA is a poor proxy for enzyme
  activity.
- Separate the two theses a target can have: single-agent dependency, and
  combination/chemosensitiser. Public data speak to the first far better than
  the second.
