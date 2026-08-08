---
name: marker-contrast-null
description: Test whether a genotype/stratum-restricted dependency is real or a marker-selection artefact, by ranking the chosen marker against an empirical null of the same contrast over every other eligible marker, plus global-shift and gene-specificity controls. Use when a biomarker hypothesis is being checked against DepMap or any per-sample screen - "are BRCA1-mutant lines more dependent on X", "is this dependency MSI-restricted", "does mutation in Y sensitise to Z knockout" - or when a single mutant-vs-wild-type contrast has produced a promising uncorrected p-value that needs a multiplicity and confounder check before it is believed.
---

# Marker-contrast empirical null

A contrast of one pre-chosen stratum against the rest — BRCA1-mutant vs
wild-type, MSI vs MSS, fusion-positive vs negative — returns an uncorrected
p-value for a hypothesis drawn from a pool of thousands of candidate markers.
That p-value answers "would this stratum look this extreme by chance?" when the
question is "would *some* stratum look this extreme by chance?" The two differ by
orders of magnitude.

Three controls separate them. Run all three before reporting a
genotype-restricted dependency.

## Workflow

1. **The contrast itself** — `stratum_contrast(values, flags)`. `values` is one
   gene's measurement per sample (CRISPR effect, drug LFC); `flags` marks the
   focal stratum. Returns n per arm, both means, Cohen's *d* and a one-sided
   Mann-Whitney p.
2. **The null scan** — `marker_null_scan(values, marker_matrix)` runs that same
   contrast for every eligible column of a sample x marker matrix, and
   `rank_in_null(scan, "BRCA1")` reports where the chosen marker lands.
   **Interpret the rank, not the marker's own p.** A marker in the top 4% of
   1,700 candidates is where roughly 68 markers were always going to be.
3. **Global shift** — `sample_gene_means(read_matrix_fn, gene_names)` then
   `global_shift_control(means, flags)`. If the focal stratum is more sensitive
   across hundreds of unrelated genes, a single-gene difference carries no
   target-specific information. This is the genetic-screen analogue of the
   proliferation confounder in drug-sensitivity work.
4. **Gene specificity** — `gene_specificity_control(effect_frame, flags)` applies
   the contrast to the focal gene beside pathway neighbours and unrelated
   controls. When neighbours and non-neighbours shift together, the result is a
   property of those samples, not of the target.

`marker_matrix` must be collapsed to one row per sample first. DepMap's omics
matrices key on `SequencingID`, so filter `IsDefaultEntryForModel == "Yes"` and
index by `ModelID`, or rows multiply per model.

## Reading the output

The scan is sorted by *d*, most negative first, with a BH `q` across the markers
actually tested. Report, in this order: the marker's *d* and uncorrected p, its
rank out of n markers, its BH q, and how many markers cleared q < 0.05 at all.
A scan where **no** marker survives correction is the common outcome and is
itself the finding.

`rank_in_null` raises when the named marker was not tested — usually because it
fell outside the stratum-size band. That is deliberate: a rank silently computed
over a scan the marker is absent from would be meaningless.

## Calibration

`MARKER_MIN_N = 10` and `MARKER_MAX_N = 200` bound the stratum-size band, with
`MIN_ARM_N = 5` per arm. The floor keeps a standardised mean difference
interpretable at small n. The ceiling drops near-ubiquitous mutations whose
reference arm is the minority and which are not candidate biomarkers — TP53 is
damaged in 780 of 1,208 DepMap lines. Widening the band changes the size of the
null, not usually a marker's position in it; re-derive by re-running with
`min_n`/`max_n` and comparing `rank_in_null(...)["percentile"]`.

Cohen's *d* raises rather than returning `±inf` when both arms are constant, so a
degenerate stratum cannot sort to the top of a scan as the apparent strongest
effect.

## Worked result

USP1 in DepMap 1,208 lines. BRCA1-damaging lines are more USP1-dependent:
−0.424 vs −0.300, *d* = −0.478, uncorrected *p* = 0.007, n = 24 vs 1,184. Taken
alone that reads as confirmation of the HR-deficiency hypothesis behind USP1
inhibitors. All three controls contradict it:

- BRCA1 ranks **67 of 1,719** markers by effect size (3.9th percentile, i.e.
  top ~4%), BH *q* = 0.99, and **zero** markers clear q < 0.05.
- No global shift: BRCA1-mutant lines average −0.1424 over 600 random genes
  against −0.1416 for wild-type, *p* = 0.34. So the result is not a fitness
  artefact — it is marker-selection noise.
- The same contrast gives FANCD2 *d* = −0.643 and PARP1 *d* = −0.514, both
  stronger than USP1's.

## Limits

The null tests whether a marker stands out among markers; it does not test
whether the underlying biology is real. n = 24 in the worked case powers only
|*d*| ≳ 0.6, so "not distinguishable from selection noise" is a weaker claim than
"no effect" — say which one the data support. Damaging-mutation calls are not
functional pathway loss: promoter hypermethylation, reversion and structural
events are invisible to a mutation matrix, so a null result may reflect stratum
misassignment rather than absent dependency. Markers whose carriers largely
coincide are not independent tests, which makes BH conservative here.

## Scope

This skill supplies the controls only. Gene-effect statistics, lineage and
mutation contrasts and the DepMap file layout belong to `depmap-local`; human
evidence and tractability to `opentargets-evidence`; the verdict vocabulary to
`depmap-fusion`; and the drug-sensitivity confounder correction to
`target-triage-public-data`. It performs no expression, structure or clinical
analysis.
