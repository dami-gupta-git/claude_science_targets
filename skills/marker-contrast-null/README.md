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

- `stratum_contrast(values, flags, alternative="less")` — one contrast: n per
  arm, both means, Cohen's *d*, one-sided Mann-Whitney p, and
  `n_excluded_no_flag` — samples in `values` with no entry in `flags` (never
  assessed for the marker) are excluded from both arms rather than folded into
  the reference arm, and this field reports how many. Raises when either arm
  falls below `MIN_ARM_N`, since a contrast that small is not interpretable as
  a standardised mean difference.
- `marker_null_scan(values, marker_matrix, min_n=None, max_n=None,
  alternative="less")` — that contrast for every eligible marker column, sorted
  by *d* with a BH q across markers tested.
- `rank_in_null(scan, marker)` — a named marker's rank, percentile, *d*, p, q,
  the count of markers clearing q < 0.05, and the fraction more extreme. Raises
  `KeyError` when the marker was not tested.
- `global_shift_control(sample_means, flags, alternative="less")` — is the focal
  stratum shifted across unrelated genes? Adds a boolean `global_shift`.
- `gene_specificity_control(effect_frame, flags, alternative="less")` — the
  contrast applied to each column, for comparing the target against pathway
  neighbours and unrelated controls.
- `sample_gene_means(read_matrix_fn, gene_names, kind="effect", n_genes=600,
  seed=1)` — per-sample mean over a random gene subset, feeding
  `global_shift_control`.
- `deletion_marker_matrix(cn_frame, cut=None)` — boolean deletion calls from a
  relative copy-number matrix, for markers that a mutation matrix cannot carry.
- `codeletion_partners(marker_matrix, marker, min_overlap=None)` — markers whose
  positive samples largely coincide with the named one, by asymmetric containment.
- `neighbourhood_check(scan, marker_matrix, marker, min_overlap=None)` — splits
  the markers outranking yours into co-deleted partners and independent rivals.
- `cohens_d(focal, reference)` and `bh_q(pvalues)` — the two statistics, exposed
  for reuse.

Module constants: `MARKER_MIN_N`, `MARKER_MAX_N`, `MIN_ARM_N`.

### Saving a run

- `mcn_run_dir(name, root="results", topic="marker_contrast_null")` — the path
  for one standalone run, `<root>/<topic>/<slug>/`, with `scripts/` created
  beside it. `name` identifies the contrast tested, not just the marker, since
  one marker can be scanned against several measurements.
- `mcn_write_run(out_dir, name, contrast, rank, scan=None, global_shift=None,
  specificity=None, neighbourhood=None, summary=None, files=(),
  data_sources=(), limits=(), scripts=())` — writes `null_scan.csv`,
  `gene_specificity.csv` (each only if the corresponding table is given), the
  run README, and copies `scripts` into `scripts/`. Returns the paths written.
- `mcn_run_readme(name, contrast, rank, global_shift=None, specificity=None,
  neighbourhood=None, summary=None, files=(), data_sources=(), limits=(),
  title=None)` — renders the README text `mcn_write_run` saves. Leads with
  `mcn_verdict()`'s plain survives / does-not-survive statement, since a ranked
  table with no stated conclusion leaves the reader to redo the judgement this
  skill exists to make. The Limits section always includes this skill's
  standing caveats (the null tests standing-out, not biology; a small arm only
  powers a large effect; mutation calls miss non-mutational loss; co-deleted
  markers are not independent tests) plus any run-specific ones passed in.
- `mcn_verdict(rank, global_shift=None, q_floor=0.05)` — does the contrast
  survive: BH q below `q_floor`, and no confirmed global shift on unrelated
  genes. Either failure alone is enough to explain the raw contrast without the
  target gene.
- `mcn_write_table(path, rows, headers=None)` — CSV writer for the null-scan or
  gene-specificity table (list of dicts or a DataFrame), used by
  `mcn_write_run` but usable standalone.

Runs land in `results/marker_contrast_null/<slug>/`. When this skill is invoked
from inside another skill's run instead, that run's own `scripts/` is the right
home for the wiring, per `coding-standards` — this writer is for the standalone
case.

## Usage

```python
scan = marker_null_scan(effect["USP1"], damaging_matrix)
rank_in_null(scan, "BRCA1")
# {'rank_by_d': 67, 'n_markers': 1719, 'percentile': 3.9, 'd': -0.4776,
#  'p': 0.0072, 'q': 0.9926, 'n_markers_q_below_05': 0, ...}
```

`marker_matrix` must be one row per sample. DepMap omics matrices key on
`SequencingID`, so filter `IsDefaultEntryForModel == "Yes"` and index by
`ModelID` before passing them in, or rows multiply per model.

## Thresholds

`MARKER_MIN_N = 10` and `MARKER_MAX_N = 200` bound eligible stratum size, with
`MIN_ARM_N = 5` per arm. The floor keeps a standardised mean difference
interpretable at small n; the ceiling excludes near-ubiquitous mutations whose
reference arm is the minority and which are not candidate biomarkers — TP53 is
damaged in 780 of 1,208 DepMap lines. To re-derive, re-run with different
`min_n`/`max_n` and compare `rank_in_null(...)["percentile"]`: on the USP1 case a
15–40 band over 400 sampled markers and the full 10–200 band over 1,719 markers
both place BRCA1 in the top 4%.

Cohen's *d* raises on two constant arms instead of returning `±inf`, which would
otherwise sort a degenerate stratum to the top of a scan as the apparent
strongest effect.

`DELETION_CUT = 0.25` calls homozygous loss on DepMap relative copy number
(1.0 = diploid). Validate it per release against the deleted gene's own knockout
effect: if the gene is truly absent, knocking it out does nothing — MTAP scores
+0.125 in called-deleted lines against −0.001 elsewhere.

`CODELETION_MIN_OVERLAP = 0.7` is the containment above which two deletion
markers are treated as one event. Calibrated on the MTAP/9p21 case, where genuine
neighbours score 0.82–1.00 and nothing independent comes close; re-derive by
listing containment for a known co-deleted locus and cutting below its minimum.
Containment rather than Jaccard is load-bearing — those same neighbours score
Jaccard 0.09–0.39, so a Jaccard rule reports every one of them as an independent
rival.

## Worked result

USP1 across 1,208 DepMap lines, from `usp1_marker_scan.csv` in the USP1 triage
outputs. BRCA1-damaging lines are more USP1-dependent (−0.424 vs −0.300,
*d* = −0.478, uncorrected *p* = 0.007, n = 24 vs 1,184), but BRCA1 ranks 67 of
1,719 markers with BH *q* = 0.99 and no marker clears q < 0.05. The
global-shift control is negative (−0.1424 vs −0.1416 over 600 random genes,
*p* = 0.34), so the signal is marker-selection noise rather than a general
sensitivity difference. The same contrast on FANCD2 gives *d* = −0.643 and on
PARP1 −0.514.

## Tests

`python -m pytest` from the skill directory; requires pytest, numpy, pandas and
scipy. Inputs are synthesised, so no DepMap release, network or credentials are
needed. Two cases carry the skill's purpose: a planted marker must rank first in
its own scan, and a globally shifted stratum must be caught by the global-shift
control. Each was mutation-checked — removing the BH correction, dropping the
size-band ceiling and inverting the *d* sign each fail a targeted test.

`test_copy_number.py` covers the deletion null on a synthesised panel: a nested
narrow deletion must still register as a co-deleted partner (the case Jaccard
gets wrong), a half-overlapping marker must not, and a genuinely independent
stronger marker must be counted as a rival. Reverting containment to Jaccard
fails two of these; the neighbour test was written against a fixture where
something actually outranks the focal marker, since a tie makes it vacuous.

## Scope

Controls only. Gene-effect statistics, lineage and mutation contrasts and the
DepMap file layout belong to `depmap-local`; human evidence and tractability to
`opentargets-evidence`; the verdict vocabulary to `depmap-fusion`; the
drug-sensitivity proliferation confounder to `target-triage-public-data`. No
expression, structure or clinical analysis.
