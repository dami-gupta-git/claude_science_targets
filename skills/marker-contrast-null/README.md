# marker-contrast-null

Controls for a biomarker-restricted dependency claim. When one stratum is
compared against the rest — BRCA1-mutant versus wild-type cell lines, MSI versus
MSS, fusion-positive versus negative — the resulting p-value is uncorrected for
the fact that the stratum was chosen from thousands of candidates. This skill
computes the empirical null of that same contrast across every other eligible
marker, so the chosen marker can be reported as a rank rather than an isolated
p-value, and adds two confounder checks: whether the stratum is globally shifted
on unrelated genes, and whether unrelated genes shift as much as the target.

The helpers are screen-agnostic. They take a per-sample measurement vector and a
sample-by-marker matrix, so they apply to CRISPR gene effect, drug log-fold
change, or any other per-sample readout.

## Functions

| Function | Description |
| --- | --- |
| `stratum_contrast(values, flags, alternative="less")` | one contrast: n per arm, both means, Cohen's *d*, one-sided Mann-Whitney p, and `n_excluded_no_flag` — samples in `values` with no entry in `flags` (never assessed for the marker) are excluded from both arms rather than folded into the reference arm, and this field reports how many. Raises when either arm falls below `MIN_ARM_N`, since a contrast that small is not interpretable as a standardised mean difference. |
| `marker_null_scan(values, marker_matrix, min_n=None, max_n=None, alternative="less")` | that contrast for every eligible marker column, sorted by *d* with a BH q across markers tested. |
| `rank_in_null(scan, marker)` | a named marker's rank, percentile, *d*, p, q, the count of markers clearing q < 0.05, and the fraction more extreme. Raises `KeyError` when the marker was not tested. |
| `global_shift_control(sample_means, flags, alternative="less")` | is the focal stratum shifted across unrelated genes? Adds a boolean `global_shift`. |
| `gene_specificity_control(effect_frame, flags, alternative="less")` | the contrast applied to each column, for comparing the target against pathway neighbours and unrelated controls. |
| `sample_gene_means(read_matrix_fn, gene_names, kind="effect", n_genes=600, seed=1)` | per-sample mean over a random gene subset, feeding `global_shift_control`. |
| `deletion_marker_matrix(cn_frame, cut=None)` | boolean deletion calls from a relative copy-number matrix, for markers that a mutation matrix cannot carry. |
| `codeletion_partners(marker_matrix, marker, min_overlap=None)` | markers whose positive samples largely coincide with the named one, by asymmetric containment. |
| `neighbourhood_check(scan, marker_matrix, marker, min_overlap=None)` | splits the markers outranking yours into co-deleted partners and independent rivals. |
| `cohens_d(focal, reference)` and `bh_q(pvalues)` | the two statistics, exposed for reuse. |
| `mcn_verdict(rank, global_shift=None, q_floor=0.05)` | the rank and controls reduced to a single reportable verdict. |

**Run output**, writing to `results/marker_contrast_null/<slug>/`

| Function | Description |
| --- | --- |
| `mcn_run_dir(name, root=None, topic=None, make=True)` | resolves the canonical output directory. |
| `mcn_write_run(out_dir, name, contrast, rank, scan=None, global_shift=None, ...)` | writes the tables and run README together. |
| `mcn_run_readme(name, contrast, rank, global_shift=None, specificity=None, ...)` | the run README text alone. |
| `mcn_write_table(path, rows, headers=None)` | a single table. |
| `mcn_link_into(canonical_dir, link_path)` | links a run into a calling skill's output directory. |
| `results_root(root=None)` | resolves the results root. |

Module constants: `MARKER_MIN_N`, `MARKER_MAX_N`, `MIN_ARM_N`, `DELETION_CUT`,
`CODELETION_MIN_OVERLAP`, `SUMMARY_MAX_WORDS`. `mcn_check_words` enforces the
word cap on run prose and carries no leading underscore only because the skill
sidecar loader reserves that prefix; the writers above call it.

Runs always land in `results/marker_contrast_null/<slug>/`, via `mcn_run_dir`
and `mcn_write_run` — whether the check was asked standalone or triggered from
inside another skill's run (a triage, a fusion, a brief). There is exactly one
location for a given contrast, so a later re-run of the same contrast — standalone
or from a different caller — lands on the same directory and the same files,
rather than silently producing a second, possibly divergent, copy.

A caller that wants the marker-null evidence to read as part of its own run
directory links to it with `mcn_link_into(canonical_dir, link_path)`, e.g.
`mcn_link_into(out_dir, os.path.join(triage_out_dir, "marker_null"))`, giving
`target_triage/chek2/marker_null -> ../../marker_contrast_null/chek2_in_..._lines/`.
The caller's README references the link path in its own prose; the tables
themselves are never copied or written a second time.

## Usage

```python
scan = marker_null_scan(effect["USP1"], damaging_matrix)
rank_in_null(scan, "BRCA1")
# {'rank_by_d': ..., 'n_markers': ..., 'percentile': ..., 'd': ...,
#  'p': ..., 'q': ..., 'n_markers_q_below_05': ..., ...}
```

`marker_matrix` must be one row per sample. DepMap omics matrices key on
`SequencingID`, so filter `IsDefaultEntryForModel == "Yes"` and index by
`ModelID` before passing them in, or rows multiply per model.

## Thresholds

`MARKER_MIN_N = 10` and `MARKER_MAX_N = 200` bound eligible stratum size, with
`MIN_ARM_N = 5` per arm. The floor keeps a standardised mean difference
interpretable at small n; the ceiling excludes near-ubiquitous mutations whose
reference arm is the minority and which are not candidate biomarkers — TP53 is
damaged in a majority of DepMap lines. To re-derive, re-run with different
`min_n`/`max_n` and compare `rank_in_null(...)["percentile"]`: on the USP1 case a
narrow band over a sampled subset and the full band over every eligible marker
place BRCA1 in the same few percent.

Cohen's *d* raises on two constant arms instead of returning `±inf`, which would
otherwise sort a degenerate stratum to the top of a scan as the apparent
strongest effect.

`DELETION_CUT = 0.25` calls homozygous loss on DepMap relative copy number
(1.0 = diploid). Validate it per release against the deleted gene's own knockout
effect: if the gene is truly absent, knocking it out does nothing, which is what
MTAP shows in called-deleted lines against the rest of the panel.

`CODELETION_MIN_OVERLAP = 0.7` is the containment above which two deletion
markers are treated as one event. Calibrated on the MTAP/9p21 locus, where
genuine neighbours score well above the cut and nothing independent approaches
it; re-derive by listing containment for a known co-deleted locus and cutting
below its minimum. Containment rather than Jaccard is load-bearing — those same
neighbours score low on Jaccard, which would report every one of them as an
independent rival.

## What the null changes

The USP1/BRCA1 contrast is the case the skill was built for. BRCA1-damaging
lines are more USP1-dependent, at a nominally significant uncorrected p-value
that would ordinarily be reported as a biomarker hypothesis. Scanned against
every other eligible marker, that same contrast ranks in the top few percent of
the empirical null but carries a BH q near 1, and no marker in the scan clears
q < 0.05 — the effect is what marker selection produces from a panel this size,
not evidence for the stratum. The global-shift control is negative, so it is not
a general sensitivity difference in those lines either, and unrelated
DNA-damage-response genes give effects of comparable size in the same lines.

Reproduce it by running `marker_null_scan` on a USP1 effect vector against a
damaging-mutation matrix, then `rank_in_null(scan, "BRCA1")`; the scan table is
written to the run directory.

## Tests

`python -m pytest` from the skill directory; requires pytest, numpy, pandas and
scipy. Inputs are synthesised, so no DepMap release, network or credentials are
needed. The suite covers the marker-null ranking, the global-shift control and
the deletion null's containment logic.

## Scope

Controls only. Gene-effect statistics, lineage and mutation contrasts and the
DepMap file layout belong to `depmap-local`; human evidence and tractability to
`opentargets-evidence`; the verdict vocabulary to `depmap-fusion`; the
drug-sensitivity proliferation confounder to `target-triage-public-data`. No
expression, structure or clinical analysis.
