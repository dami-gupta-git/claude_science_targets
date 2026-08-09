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

| Function | Description |
| --- | --- |
| `ot_target(ensembl_id)` | Open Targets tractability and DepMap essentiality for one Ensembl gene id, in a single GraphQL call. |
| `ot_essentiality_frame(target_dict)` | flattens that response into a tidy DataFrame of gene, tissue, cell line, disease, gene effect and expression. |
| `partial_corr(x, y, z)` | Pearson correlation of x and y after linearly removing z from both. The p-value uses n−3 degrees of freedom, since fitting z out of each variable consumes a parameter. Raises on a constant z, on non-finite input, and below n = 4. |
| `confounder_check(df, expr_col, drug_cols)` | tests whether expression correlates with sensitivity across the entire drug panel. Returns `r`, `p`, `n`, `n_drugs` and `n_drugs_min`, the last being the thinnest drug coverage behind any surviving row. |
| `sensitivity_scan(df, expr_col, drug_cols)` | per-drug raw and partial correlations against one expression column, Benjamini-Hochberg corrected. Columns are `SCAN_COLUMNS` whether or not any drug survived, so a scan in which every drug fell below `min_n` still exposes `r_partial` and `q_partial_BH` as empty columns. Drugs that could not be tested are listed in `.attrs["dropped"]` with their n and the reason, and summarised on stdout. |
| `classify_gene_effect(effect)` | labels a Chronos gene effect `essential`, `intermediate`, `dispensable` or `unknown` against `DEPMAP_ESSENTIAL_THRESHOLD`. |
| `load_prism(matrix_path, compound_list_path)` | PRISM primary matrix as a cell line by drug frame keyed by drug name rather than BRD id. |
| `load_prism_secondary(path)` | PRISM 20Q2 dose-response AUC, indexed by ModelID. |
| `prism_secondary_has(path, drug_names)` | which of the named compounds were carried into the dose-response screen. |
| `normalize_cell_name(s)` | strips case and punctuation so GDSC and CCLE line names join. |
| `reconcile_fetch(expected_ids, fetched_ids)` | asserts a batched fetch returned every id the search promised, and returns the missing ids. |
| `triage_readme(gene, summary, steps, files, data_sources, limits)` | renders a run README as markdown: an opening plain-prose paragraph, one section per triage step with an optional table, then Files, Data sources and Limits. |
| `write_triage_readme(path, ...)` | the same, written to `path`. |
| `markdown_table(rows, headers=None)` | a list of dicts or a DataFrame as a markdown table; p-value and q-value columns are formatted, missing values render as `—`. |
| `is_pvalue_key(key)` | whether a column name holds a p-value or an FDR-adjusted q-value, matching `p` or `q` as a whole underscore-separated component so that `q_partial_BH` is caught and `quantile` is not. |
| `format_p(p)` | renders a p-value for prose, in scientific notation below 0.001 and as a bound when it has underflowed to zero. |

## Where a run is written

A triage writes one directory, `results/target_triage/<gene>/`, holding every
table and figure the run produced under a `<gene>_` prefix, its wiring in
`scripts/`, and a `README.md`. `depmap-local` and `opentargets-evidence` have
no run-writer of their own, so their outputs are part of this run and are
written directly here.

`marker-contrast-null` and `depmap-fusion` each have their own canonical
output directory (`results/marker_contrast_null/<slug>/`,
`results/depmap_fusion/<slug>/`) that they always write to, standalone or not
— so a re-run of the same contrast or subject outside this triage lands on
the same files instead of a second, divergent copy. The triage links to those
directories with `mcn_link_into`/`fusion_link_into` rather than duplicating
their tables: `results/target_triage/chek2/marker_null` is a symlink, not a
copy.

| Function | Description |
| --- | --- |
| `results_root(root=None)` | resolves this repo's `results/` directory from `$SCIENCE_RESULTS_ROOT` or an explicit `root=`; raises naming the variable when neither is set, rather than silently creating a `results/` folder wherever the kernel session's cwd happens to be. A Claude Science kernel does not read the repo's `.env`, so the variable is set in the session. |
| `triage_run_dir(gene, root=None, topic=None)` | the path for one triage run, `<root>/<topic>/<slug>/`, with `scripts/` created beside it. Called before the analysis, so every step writes to a path that already exists. |
| `triage_write_run(out_dir, gene, summary, steps, files=(), data_sources=(), limits=(), title=None, scripts=())` | wraps `write_triage_readme` rather than replacing it: writes `README.md` into `out_dir` and copies `scripts` into `scripts/`. Returns the paths written. |

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

`SUMMARY_MAX_WORDS = 130` and `FINDING_MAX_WORDS = 90` cap the run README's
opening paragraph and each step finding, enforced by `check_words`. They are
editorial rather than derived, and exist so a run README stays an orientation
document rather than a second copy of the analysis. Raise them in `kernel.py`
rather than splitting text across fields to evade them.

## Error paths

The helpers raise rather than returning a degraded number, because every one of
these inputs otherwise yields a value indistinguishable from a valid result:

- `partial_corr` raises on a constant confound, where a rank-deficient fit
  leaves the residuals equal to the original variables and the return value
  would be an uncontrolled correlation presented as a controlled one. It also
  raises on non-finite input rather than returning `(nan, nan)`, and below
  n = 4.
- `sensitivity_scan` raises when built with a single drug column and no
  externally supplied `panel_mean_col`, since the panel mean of one column is
  that column: the drug would be regressed on itself and expression correlated
  against the floating-point residual, fabricating an `r_partial`
  indistinguishable from a legitimate null. Pass `panel_mean_col=` from a wider
  panel, or use `confounder_check` for a single-drug question.
- `sensitivity_scan` and `confounder_check` raise on absent drug or expression
  columns, and when the frame already holds the internal `_panel_mean` columns
  the computation would overwrite.
- `load_prism` raises on duplicate cell-line ids, which multiply rows on any
  later merge, and on a row whose length disagrees with the header, which
  indicates a truncated download.
- `ot_target` raises on a GraphQL `errors` payload and on a null target, both of
  which arrive with HTTP 200 for an unrecognised Ensembl id.

## Analysis constraints

- **The panel-mean confound is leave-one-out.** Each drug is excluded from the
  panel mean used as its own confound; including it regresses part of the drug's
  signal out of itself and attenuates `r_partial` toward zero, severely on a
  narrow panel and negligibly on a wide one. The bias is conservative and cannot
  manufacture a hit, but it can hide a borderline one. `confounder_check` is
  unaffected: there the panel mean is the quantity under test, not a confound.

**Proliferation confounding.** Genes that track proliferation rate correlate
with sensitivity to every cytotoxic drug, so per-drug correlations look specific
when they are not. On the worked DCTPP1 case a nominally significant raw
correlation with 5-FU disappears once each line's mean response across the whole
panel is partialled out, and an independent PRISM panel reproduces that
conclusion at different numbers, the two screens differing in cell lines,
readout and n. Correlations from different panels are not interchangeable, so
every reported number should name the screen it came from. `confounder_check`
and `sensitivity_scan` make the partial correlation the default path.

**PRISM product selection.** DepMap ships three repurposing datasets whose
version numbers do not order them by content: `Repurposing_Public_23Q2`/`24Q2`
are single-dose primary screens, while the 8-point dose-response data with
fitted AUC/IC50 sits on the older `19Q4`/`20Q2` release page. Only
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
scipy and statsmodels. All inputs are synthesised, so no DepMap release, GDSC
download, network access or credentials are needed. The suite covers the
run-README generator, the partial-correlation statistic, and the
degrees-of-freedom and column-consistency guards.

## Scope

The skill assembles and controls the evidence; it does not decide whether a
target is worth pursuing, which turns on mechanism, competition and indication
and stays with the user. Dependency analysis against a locally held DepMap
release is handled by the `depmap-local` skill, which the guidance defers to
instead of going over the network when such a release is present.
