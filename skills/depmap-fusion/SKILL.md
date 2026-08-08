---
name: depmap-fusion
description: Fuse Open Targets human evidence with local DepMap CRISPR dependency into a single target dossier or a ranked disease-level triage table, classifying each target as concordant-dependency, growth-suppressive-mismatch, context-restricted, inert-in-panel, or evidence-without-dependency. Use when asked whether a target is actually worth pursuing, to triage or prioritise targets for an indication, to sanity-check a target list against functional genomics, or when human genetic evidence and cell-line dependency need reconciling.
---

# Open Targets x DepMap fusion

Neither knowledge base answers the target question alone. Open Targets says
*this gene is linked to this disease in humans* but encodes **no direction** -
tumour suppressors rank at the top. DepMap says *cancer cells need this gene*
but has no human disease evidence. The interesting signal is where they
**disagree**.

Depends on `opentargets-evidence` (evidence) and `depmap-local` (dependency).
Load both; this skill only adds the join and the verdict.

## Workflow

1. **Evidence** - `open_targets_disease_targets(efo_id, size=N)` for a disease,
   or the dossier query for one target. Collect `ot_score` and tractability
   modalities. MCP calls run in the `repl` tool; write results to
   `handoff/*.json` for the analysis kernel.
2. **Dependency** - `depmap_selectivity(gene)` per target, plus the
   disease-relevant lineage slice from `Model.csv`.
3. **Fuse** - `fuse_target_row(symbol, ot_score, effect_stats, lineage_row,
   modalities)`. Feed it the full `depmap_selectivity()` dict: the verdict needs
   `sd_effect` and `frac_positive`, not just the mean.
4. Report the verdict mix, and always the joined n.

## Verdicts

| Verdict | Meaning |
|---|---|
| `concordant-dependency` | Human evidence and cell dependency agree. Strongest KO/degrader case. |
| `common-essential-no-window` | ~Every line depends on it. Real, but normal tissue needs it too - no therapeutic window without tumour-selective delivery. |
| `growth-suppressive-mismatch` | Knockout **helps** growth (positive tail, high variance). A suppressor - not a KO target. Pivot to synthetic lethality. |
| `context-restricted` | Dependency confined to a lineage/genotype. Needs a biomarker hypothesis. |
| `inert-in-panel` | Flat, low-variance, near zero. The panel lacks the relevant context. |
| `evidence-without-dependency` | Human evidence without cell-autonomous dependency - non-cell-autonomous, redundant, or wrong model. |
| `dependency-with-thin-evidence` | Evidence above `EVIDENCE_FLOOR` but below 0.5, plus a dependency signal. Both halves real; actionable on weaker evidence than concordant. |
| `dependency-without-evidence` | Cells need it, no human evidence at all (`ot_score` below `EVIDENCE_FLOOR` or absent). Novel target, or a fitness gene with a poor window. |
| `no-evidence-no-dependency` | Neither side supports it: OT score below `EVIDENCE_FLOOR` (0.10) or absent, AND no dependency. |
| `indeterminate` | `sd_effect` / `frac_positive` / `mean_effect` missing, so the suppressive and inert tests could not run. Pass the full `depmap_selectivity()` dict. |

## Calibration

Separating a real suppressor from an inert gene needs **dispersion**, not the
mean. Measured on this release:

| Gene | mean | sd | frac > +0.3 | truth |
|---|---|---|---|---|
| TP53 | +0.42 | 0.58 | 0.42 | suppressor |
| RB1 | +0.18 | 0.23 | 0.24 | suppressor |
| NRG1 | +0.19 | 0.08 | 0.08 | inert (ligand-driven) |
| ROS1 | +0.11 | 0.11 | 0.05 | inert (fusion-driven) |

Hence `SUPPRESSOR_MIN_FRAC_POS = 0.15` with `SUPPRESSOR_MIN_SD = 0.20`.
Mean-only logic misclassifies ROS1/NRG1 as suppressors - it was the first
version of this skill and it was wrong.

## Worked results

**Lung adenocarcinoma** (MONDO_0005061), top 30 OT targets x 126 lung lines:
concordant 11 (EGFR, KRAS, ERBB2, KEAP1) | inert-in-panel 11 (ALK, ROS1, RET,
EML4, NRG1 - fusion/ligand-driven, unrepresented) | growth-suppressive 3
(TP53, RB1) | evidence-without-dependency 4 | common-essential 1.

**Pancreatic** (MONDO_0009831), top 25 x 48 pancreas lines:
concordant 6 | common-essential-no-window 6 (POLD1, POLE, RRM1/2, PRIM2) |
inert 5 | evidence-without-dependency 4 (incl. SMAD4, a deleted suppressor) |
context-restricted 2 | growth-suppressive 2 (TP53, CDKN2A).
Actionable set: KRAS, BRCA1/2, PALB2, STK11, ARID1A, RBM10, TYMS - i.e. KRAS
plus the HR-deficiency genes behind PARP-inhibitor use in this indication.

**KRAS is `context-restricted` in pancreas but `concordant` in lung**, because
pancreatic lines are far more KRAS-dependent (-2.14 vs -0.73 global). The
verdict tracks the data, not a hardcoded per-gene rule.

Across both, **65% of the top OT-ranked targets are not knockout-actionable**
(36 of 55: lung 19/30, pancreatic 17/25). That is the finding, not a failure of
either source.

## Honest limits

- `inert-in-panel` is **absence of evidence**. Check lineage coverage before
  dismissing a target; a fusion-driven oncogene needs a fusion-positive line.
- DepMap is proliferation in 2D culture: immune, stromal, angiogenic and
  differentiation-dependent targets are invisible by construction.
- A `concordant` verdict is a prioritisation signal, not validation.
- Check `depmap_class` before trusting `knockout_actionable`: pan-essential
  genes are genuinely required but have no window. An earlier version of this
  skill listed POLD1/POLE/RRM1 as actionable - the guard exists because of it.
- **Each verdict asserts only what was supplied.** `evidence-without-dependency`
  and `dependency-with-thin-evidence` both require `ot_score >= EVIDENCE_FLOOR`
  (0.10); rows below the floor or `None` go to `no-evidence-no-dependency` or
  `dependency-without-evidence`. Evidence in [0.10, 0.50) with
  `frac_dependent >= 0.10` is its own tier rather than being folded into either
  neighbour, because both halves of the signal are real.
- Thresholds are conventions (`DEP_THRESHOLD = -0.5`), not decision boundaries;
  report the underlying numbers alongside any verdict.
