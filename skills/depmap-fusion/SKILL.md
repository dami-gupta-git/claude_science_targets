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
5. **Save the run** - `fusion_write_run(out_dir, subject, rows, summary, ...)`
   into `results/depmap_fusion/<subject>/`, not assembled by hand. A single row
   renders as a target dossier; several render as a triage table with the
   verdict mix. `out_dir` comes from `fusion_run_dir(subject)`, which needs
   `$SCIENCE_RESULTS_ROOT` set to this repo's `results/` directory (or an
   explicit `root=`) - it raises naming the variable rather than silently
   creating a `results/` folder wherever the session's cwd happens to be:

```python
rows = [fuse_target_row("EGFR", ot_score, effect_stats, lineage_row, modalities), ...]
out_dir = fusion_run_dir("lung adenocarcinoma")
fusion_write_run(out_dir, "lung adenocarcinoma", rows,
                 summary="...", disease="MONDO_0005061",
                 data_sources=["Open Targets Platform", "DepMap 24Q2"])
```

`fusion_run_dir`/`fusion_write_run` land in `results/depmap_fusion/<slug>/`
every time, whether the fusion was asked standalone or from inside a larger
run (a target-triage brief, say) — never at a path the caller picks. That way
a later standalone re-run of the same subject finds the existing directory
instead of producing a second, possibly divergent, copy. A caller that wants
the fusion evidence to read as part of its own run directory links to it with
`fusion_link_into(out_dir, os.path.join(caller_out_dir, "fusion"))` rather
than copying the tables in.

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

Hence `SUPPRESSOR_MIN_FRAC_POS = 0.15` with `SUPPRESSOR_MIN_SD = 0.20`. All
four genes have positive mean effects, so the mean alone does not separate
them: classifying on the mean labels ROS1 and NRG1 as suppressors when they are
untested in this panel.

## What a sweep returns

A top-N sweep of an OT disease ranking does not return a list of targets to
pursue. Tumour suppressors rank near the top and classify
`growth-suppressive-mismatch`; fusion- and ligand-driven oncogenes classify
`inert-in-panel` because the panel holds few or no lines carrying the driving
context; replication genes classify `common-essential-no-window`. The
actionable set is typically a minority of the list, and the discarded verdicts
carry the reason.

**The verdict tracks the lineage slice supplied, not a per-gene rule.** `KRAS`
has a global mean effect of -0.72 across 1208 lines, against -2.14 in the 48
pancreas lines and -0.88 in the 126 lung lines. Passed a pancreas
`lineage_row`, it clears the context test (lineage mean below the global mean
by more than 0.2, `p_bh < 0.05`) and returns `context-restricted`; passed lung,
it does not, and returns `concordant-dependency` on its evidence and
dependency alone.

## Honest limits

- `inert-in-panel` is **absence of evidence**. Check lineage coverage before
  dismissing a target; a fusion-driven oncogene needs a fusion-positive line.
- DepMap is proliferation in 2D culture: immune, stromal, angiogenic and
  differentiation-dependent targets are invisible by construction.
- A `concordant` verdict is a prioritisation signal, not validation.
- Check `depmap_class` before trusting `knockout_actionable`: pan-essential
  genes such as POLD1, POLE and RRM1 are genuinely required by cancer lines but
  have no therapeutic window, and are excluded from the actionable set for that
  reason.
- **Each verdict asserts only what was supplied.** `evidence-without-dependency`
  and `dependency-with-thin-evidence` both require `ot_score >= EVIDENCE_FLOOR`
  (0.10); rows below the floor or `None` go to `no-evidence-no-dependency` or
  `dependency-without-evidence`. Evidence in [0.10, 0.50) with
  `frac_dependent >= 0.10` is its own tier rather than being folded into either
  neighbour, because both halves of the signal are real.
- Thresholds are conventions (`DEP_THRESHOLD = -0.5`), not decision boundaries;
  report the underlying numbers alongside any verdict.
