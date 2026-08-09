# depmap-fusion

Answers whether a gene is worth pursuing as a drug target for a disease, by
combining two sources that are each misleading alone. Open Targets scores how
strongly a gene is linked to a disease in humans, but encodes no direction:
`TP53` scores 0.74 for lung adenocarcinoma while being a tumour suppressor
whose knockout accelerates growth. DepMap measures whether cell lines still
grow after CRISPR knockout of a gene, but carries no disease link, and the
genes every line needs are the ones normal tissue needs too. This skill crosses
the two and returns one verdict per target, plus a `knockout_actionable` flag.

Strong human evidence with a matching cell-line dependency gives
`concordant-dependency`. Strong evidence where knockout instead *helps* growth
gives `growth-suppressive-mismatch` — the case a ranking alone would hide. A
real dependency that every line shares gives `common-essential-no-window`. Ten
verdicts in total, of which three mean worth pursuing.

## Functions

**Classification**

| Function | Description |
| --- | --- |
| `fuse_target_row(symbol, ot_score, effect_stats, lineage_row=None, tractable_modalities=None, dep_threshold=-0.5)` | classifies one target and returns the verdict with the numbers behind it. `effect_stats` must be the **full** `depmap_selectivity()` dict: the verdict reads `sd_effect` and `frac_positive`, not just the mean. Pass `lineage_row` to enable `context-restricted`. |
| `fusion_notes()` | the interpretation guide for the verdict vocabulary, keyed by verdict. |

**Run output**, writing to `results/depmap_fusion/<slug>/`

| Function | Description |
| --- | --- |
| `fusion_run_dir(subject, root=None, topic=None, make=True)` | resolves the canonical output directory. Requires `$SCIENCE_RESULTS_ROOT` or an explicit `root=`, and raises naming the variable rather than creating a stray `results/` folder at the session cwd. |
| `fusion_write_run(out_dir, subject, rows, summary, disease=None, files=(), ...)` | writes the verdict table and run README together. One row renders as a target dossier, several as a triage table. |
| `fusion_verdict_mix(rows)` | verdict counts for a set of rows. |
| `fusion_write_table(path, rows, headers=None)` | the verdict table alone. |
| `fusion_markdown_table(rows, headers=None)` | the same rows as markdown. |
| `fusion_run_readme(subject, rows, summary, disease=None, files=(), ...)` | the run README text alone. |
| `fusion_link_into(canonical_dir, link_path)` | links a run into a calling skill's output directory, so a larger run can reference the fusion evidence without copying the tables. |
| `results_root(root=None)` | resolves the results root. |

Module constants: `VERDICTS`, `EVIDENCE_FLOOR`, `SUPPRESSOR_MIN_FRAC_POS`,
`SUPPRESSOR_MIN_SD`, `INERT_MAX_SD`, `INERT_MAX_ABS_MEAN`.

`fusion_check_words` carries no leading underscore only because the skill
sidecar loader reserves that prefix. It enforces the word cap on run-README
prose and is called by the writers above rather than directly.

## Companion skills

This skill adds only the join and the verdict. Evidence retrieval and
dependency statistics live elsewhere, and functions named `depmap_*` or
`open_targets_*` below belong to those skills, not this one.

| Skill | Answers | Source |
|---|---|---|
| `depmap-local` | Do cancer cells need this gene? | local DepMap files |
| `opentargets-evidence` | Is it linked to disease in humans, and druggable? | Open Targets GraphQL |
| `depmap-fusion` | Do those two agree? | both |

## Workflow

1. **Evidence** — `open_targets_disease_targets(efo_id, size=N)` for a disease,
   or the dossier query for one target (`opentargets-evidence`, in a `repl`
   cell). Collect `ot_score` and tractability modalities into `handoff/*.json`.
2. **Dependency** — `depmap_selectivity(gene)` per target, plus the
   disease-relevant lineage slice from `depmap_lineage_enrichment()`
   (`depmap-local`, in a `python` cell).
3. **Fuse** — `fuse_target_row(...)` per gene, then `fusion_write_run(...)`.
   Report the verdict mix and the joined n.

```python
st = depmap_selectivity("TP53")          # from depmap-local
fuse_target_row("TP53", ot_score=0.74, effect_stats=st)
# verdict='growth-suppressive-mismatch', knockout_actionable=False

st = depmap_selectivity("KRAS")
fuse_target_row("KRAS", ot_score=0.82, effect_stats=st)
# verdict='concordant-dependency', knockout_actionable=True
```

## The ten verdicts

| Verdict | Meaning | Actionable |
|---|---|---|
| `concordant-dependency` | Human evidence and cell dependency agree. Strongest KO/degrader case. | yes |
| `context-restricted` | Dependency confined to a lineage or genotype. Needs a biomarker. | yes |
| `dependency-with-thin-evidence` | `ot_score` in [0.10, 0.50) with `frac_dependent >= 0.10` — both halves real, evidence weaker than concordant. | yes |
| `common-essential-no-window` | Real dependency, but every cell has it — normal tissue included. | no |
| `growth-suppressive-mismatch` | Knockout **helps** growth (positive tail, high variance). A suppressor; pivot to synthetic lethality. | no |
| `inert-in-panel` | Flat near zero. The panel lacks the context that makes the gene matter. | no |
| `evidence-without-dependency` | Real human evidence (`ot_score >= EVIDENCE_FLOOR`) but no cell-autonomous dependency. | no |
| `dependency-without-evidence` | Cells need it, no human evidence at all. Possible novel target. | no |
| `no-evidence-no-dependency` | Neither side supports it. | no |
| `indeterminate` | Dispersion inputs missing, so the suppressive and inert tests could not run. | no |

Branch order is load-bearing: indeterminate, then common-essential, then
suppressive, then context, then concordant, then inert; first match wins.
Below the strong-evidence bar the thin-evidence tier is tested before either
single-sided verdict, so a row with both real evidence and a dependency signal
is never labelled as lacking one of them. A gene's verdict therefore depends on
which branch catches it first, and on the lineage slice supplied — `KRAS` is
`context-restricted` against pancreas lines and `concordant-dependency` against
lung.

## Reading the numbers

- Gene effect `0` means knockout changed nothing, `-1` that cells grow badly
  without the gene, `-2` or lower that they largely die. A **positive** effect
  means knockout *helps* growth; never read it as a weak dependency.
- `< -0.5` is the conventional dependent cut (`DEP_THRESHOLD`);
  `frac_dependent >= 0.90` marks a common essential.
- Separating a suppressor from an inert gene needs dispersion, not the mean.
  All four genes below have positive means, and the mean alone does not
  distinguish them:

| Gene | mean | sd | frac > +0.3 | reading |
|---|---|---|---|---|
| TP53 | +0.42 | 0.58 | 0.42 | suppressor |
| RB1 | +0.18 | 0.23 | 0.24 | suppressor |
| NRG1 | +0.19 | 0.08 | 0.08 | inert (ligand-driven) |
| ROS1 | +0.11 | 0.11 | 0.05 | inert (fusion-driven) |

  Hence `SUPPRESSOR_MIN_FRAC_POS = 0.15` and `SUPPRESSOR_MIN_SD = 0.20`. These
  are conventions, not decision boundaries — report the underlying numbers
  alongside any verdict.

## Limits

- `common-essential-no-window` uses this skill's `frac_dependent >= 0.90` rule
  (974 genes on the current release), not DepMap's shipped
  `CRISPRInferredCommonEssentials.csv` (1827 genes). The shipped list is
  broader and includes `KRAS`, which adopting it would discard from the
  actionable set; DepMap's call is surfaced as `depmap_inferred_essential` for
  cross-reference.
- `inert-in-panel` is absence of evidence, and is testable:
  `depmap_fusion_contrast()` (`depmap-local`) separates "no fusion-positive
  line exists" from "fusion-positive lines exist and still show no
  dependency". Check lineage coverage before dismissing a target.
- DepMap measures proliferation in 2D culture. Immune, stromal, angiogenic and
  differentiation-dependent targets are invisible by construction.
- Common essentials are excluded from `knockout_actionable`, so check
  `depmap_class` before trusting the flag.
- The thin-evidence tier does not appear in a top-25/30 sweep, where scores sit
  above the strong-evidence bar. It applies to deeper slices, where most
  targets fall in [0.10, 0.50).
- No verdict is validation. `concordant-dependency` is a prioritisation signal.

## Scope

This skill does not retrieve evidence or compute dependency statistics. Open
Targets queries belong to `opentargets-evidence`; gene-effect statistics,
lineage and mutation contrasts, fusion contrasts, PRISM compound summaries, the
file inventory and the parquet cache belong to `depmap-local`. It performs no
chemistry, structure or expression analysis.
