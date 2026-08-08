# target-triage-public-data

A Claude Science skill for assessing whether a gene is a plausible drug target
from public data alone: genetic dependency (DepMap via Open Targets),
tractability, expression-versus-drug-sensitivity in cell line panels (GDSC and
PRISM Repurposing), and TCGA survival and subtype correlates. It fixes an order
of operations that puts the cheapest disconfirming evidence first, and fixes the
controls and dataset choices at the three points where the obvious analysis
returns a wrong answer without raising an error. `SKILL.md` is the guidance
document an agent loads into context, and `kernel.py` is a sidecar of helper
functions that loads into the Python kernel alongside it.

## Functions

- `ot_target(ensembl_id)` — Open Targets tractability and DepMap essentiality
  for one Ensembl gene id, in a single GraphQL call.
- `ot_essentiality_frame(target_dict)` — flattens that response into a tidy
  DataFrame of gene, tissue, cell line, disease, gene effect and expression.
- `partial_corr(x, y, z)` — Pearson correlation of x and y after linearly
  removing z from both. The p-value uses n−3 degrees of freedom, since fitting z
  out of each variable consumes a parameter. Raises on a constant z, on
  non-finite input, and below n = 4.
- `confounder_check(df, expr_col, drug_cols)` — tests whether expression
  correlates with sensitivity across the entire drug panel. Returns `r`, `p`,
  `n`, `n_drugs` and `n_drugs_min`, the last being the thinnest drug coverage
  behind any surviving row.
- `sensitivity_scan(df, expr_col, drug_cols)` — per-drug raw and partial
  correlations against one expression column, Benjamini-Hochberg corrected.
  Drugs that could not be tested are listed in `.attrs["dropped"]` with their n
  and the reason, and summarised on stdout.
- `classify_gene_effect(effect)` — labels a Chronos gene effect `essential`,
  `intermediate`, `dispensable` or `unknown` against
  `DEPMAP_ESSENTIAL_THRESHOLD`.
- `load_prism(matrix_path, compound_list_path)` — PRISM primary matrix as a cell
  line by drug frame keyed by drug name rather than BRD id.
- `load_prism_secondary(path)` — PRISM 20Q2 dose-response AUC, indexed by
  ModelID.
- `prism_secondary_has(path, drug_names)` — which of the named compounds were
  carried into the dose-response screen.
- `normalize_cell_name(s)` — strips case and punctuation so GDSC and CCLE line
  names join.
- `reconcile_fetch(expected_ids, fetched_ids)` — asserts a batched fetch
  returned every id the search promised, and returns the missing ids.

## Thresholds

`DEPMAP_ESSENTIAL_THRESHOLD = -1.0` is not a fitted cutoff. Chronos gene effect
is scaled so that 0 is no effect and −1 is the median common-essential gene, so
the value is a property of the DepMap scaling and changes only if DepMap changes
that scaling; `classify_gene_effect` places `intermediate` at −0.5 by the same
convention. Neither is re-derivable from a screen, and neither substitutes for
comparators measured on the same cell lines — a gene effect of −0.6 means
something different in a panel where RRM1 reads −2.5 than in one where it reads
−0.9, which is why the guidance requires reference genes to be reported
alongside.

`sensitivity_scan` defaults to `min_n = 15` cell lines per drug. Drugs below it
are excluded from the returned rows and recorded in `.attrs["dropped"]`.

## Error paths

The helpers raise rather than returning a degraded number, because every one of
these inputs otherwise yields a value indistinguishable from a valid result:

- `partial_corr` raises on a constant confound. With a rank-deficient fit the
  residuals are the original variables, so the return value would be the
  uncontrolled Pearson correlation presented as a controlled one.
- `partial_corr` raises on non-finite input rather than returning `(nan, nan)`,
  and below n = 4, where no residual degrees of freedom remain.
- `sensitivity_scan` raises when built with a single drug column and no
  externally supplied `panel_mean_col`. The whole-panel mean of one column IS
  that column, so the confound would equal the drug being tested: `partial_corr`
  would then regress the drug on itself and correlate expression against the
  ~1e-16 floating-point residual of that perfect self-fit — a fabricated
  `r_partial` indistinguishable from a legitimate null result, rather than the
  "nothing to control for" error it actually is. Pass `panel_mean_col=` from a
  wider panel, or use `confounder_check` for a single-drug question.
- `sensitivity_scan` and `confounder_check` raise when a caller's frame already
  holds `_panel_mean` or `_panel_mean_loo`, which the panel-mean computation
  would otherwise overwrite, and on drug or expression columns absent from the
  frame.
- `load_prism` raises on duplicate cell-line ids in the matrix header, which
  produce a duplicated index that multiplies rows on any later merge, and on a
  data row whose length disagrees with the header, which indicates a truncated
  download.
- `ot_target` raises on a GraphQL `errors` payload and on a null target, both of
  which arrive with HTTP 200 for an unrecognised Ensembl id.

## Analysis constraints

- **The panel-mean confound is leave-one-out.** Each drug is excluded from the
  panel mean used as its own confound; including it regresses part of the
  drug's signal out of itself and attenuates `r_partial` toward zero by
  roughly `1/len(drug_cols)` (~46% on a 3-drug panel, ~10% at 25, ~3% at
  100). The bias is conservative and cannot manufacture a hit, but it can
  hide a borderline one. `confounder_check` is unaffected: there the panel
  mean is the quantity under test, not a confound.

**Proliferation confounding.** Genes that track proliferation rate correlate
with sensitivity to every cytotoxic drug, so per-drug correlations look specific
when they are not. Measured on DCTPP1 in colorectal lines, **GDSC** (40 COREAD
lines, lnIC50) gave a raw 5-FU correlation of r = −0.27, a panel-wide
r = −0.41 (p = 0.009), and a 5-FU association of r = +0.04 (p = 0.80) after partialling out each
line's mean response across all drugs. The independent **PRISM** panel (33 bowel
lines, single-dose log-fold-change) gave raw r = −0.08 and partial
r = +0.09 (p = 0.61) for the same drug — the same qualitative conclusion at
different numbers, because the two screens differ in cell lines, readout and n.
Both figures are from `dctpp1_drug_sensitivity_correlations.csv`;
`dctpp1_drug_sensitivity_all_screens.csv` reports the PRISM row separately per
release (23Q2 and 24Q2) and differs in the third decimal, while the GDSC row is
identical across both files.
Correlations from different panels are not interchangeable, so every reported
number should name the screen it came from. `confounder_check` and
`sensitivity_scan` make the partial correlation the default path.

**PRISM product selection.** DepMap ships three repurposing datasets whose
version numbers do not order them by content: `Repurposing_Public_23Q2`/`24Q2`
are single-dose (~2.5 uM) primary screens, while the 8-point dose-response data
with fitted AUC/IC50 sits on the older `19Q4`/`20Q2` release page. Only
primary-screen hits were carried into the dose-response screen, so membership
must be checked before promising the analysis (`prism_secondary_has`). AUC
polarity is inverted relative to the primary screen's log-fold-change: higher
AUC means more resistant.

**Silent losses in batched API fetches.** A literature search reporting 37 hits
followed by a batched metadata call can return 28 records with no error raised,
and the resulting counts remain self-consistent. `reconcile_fetch` asserts the
two agree.

## Data sources

- **Open Targets Platform** GraphQL API — tractability, DepMap CRISPR
  essentiality (Chronos).
- **GDSC** release 8.4, Wellcome Sanger Institute — fitted dose-response
  (lnIC50, AUC).
- **PRISM Repurposing Public** 23Q2 / 24Q2, Broad Institute — single-dose
  primary screen; covers nucleoside analogues GDSC lacks (trifluridine,
  decitabine, azacitidine, floxuridine, capecitabine).
- **PRISM Repurposing 20Q2 secondary screen**, Broad Institute — 8-point
  dose-response with fitted AUC/IC50; preferred readout when the compound is
  present. Screens described in Corsello et al. 2020,
  doi:10.1038/s43018-019-0018-6.
- **cBioPortal** — CCLE expression and cell-line annotation; TCGA PanCancer
  Atlas.
- **NCI GDC** API — TCGA treatment records, which cBioPortal does not carry.

DepMap's portal is behind bot verification, so PRISM files must be downloaded
manually rather than fetched from a script.

## Installation

Place the directory in your Claude Science skills folder:

```
~/.claude-science/orgs/<org-id>/skills/target-triage-public-data/
```

`kernel.py` requires `numpy`, `pandas`, `scipy` and `statsmodels`. `lifelines`
is needed for the survival step and is imported only where used.

## Tests

`python -m pytest` from the skill directory; requires pytest, numpy, pandas,
scipy and statsmodels. `pytest.ini` pins rootdir here so pytest does not walk up
to an unrelated config, and `conftest.py` fails with an actionable message when
scipy or statsmodels is missing. All inputs are synthesised, so no DepMap
release, GDSC download, network access or credentials are needed.

55 tests across four modules. `test_partial_corr.py` checks the partial
correlation coefficient against the closed-form first-order identity computed
from independent `pearsonr` calls, so the test cannot pass merely by agreeing
with the implementation's residual-regression route. `test_guards.py` checks the
p-value against the t-distribution at n−3 degrees of freedom, and asserts it
exceeds the value `pearsonr` would report on the residuals — the direction that
fails if the n−2 form is ever restored. Each guard was mutation-checked by
reverting it and confirming the intended test, and only that test, fails.

## Scope

The skill assembles and controls the evidence; it does not decide whether a
target is worth pursuing, which turns on mechanism, competition and indication
and stays with the user. Dependency analysis against a locally held DepMap
release is handled by the `depmap-local` skill, which the guidance defers to
instead of going over the network when such a release is present.
