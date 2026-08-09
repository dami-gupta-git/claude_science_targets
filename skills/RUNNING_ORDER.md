# Running order for the skills

These are Claude Science skills, loaded with `skill("name")` inside a `repl`
session, not standalone scripts — "running" them means invoking them in one
conversation in the order below, not separate CLI calls. They do not all
belong to one pipeline: there is a target-judgment stack, a separate
molecule-triage pipeline that consumes the target stack's output, and a set of
skills that stand alone.

## Judging a drug target

1. **`target-triage-public-data`** — entry point. Internally sequences
   dependency → tractability → drug-sensitivity → survival, cheapest
   disconfirming evidence first. Its own step 1 already tells you to prefer
   `depmap-local` over its built-in fallback when a local DepMap release is
   available, and similarly points at Open Targets for tractability.
2. **`depmap-local`** and **`opentargets-evidence`** — invoked by step 1
   directly, or run separately first if you want full detail rather than
   `target-triage-public-data`'s reduced-coverage fallback. No ordering
   between the two; they read independent sources (local DepMap CRISPR data,
   Open Targets human evidence).
3. **`depmap-fusion`** — depends on both of the above by name in its own
   `SKILL.md` ("Depends on `opentargets-evidence` ... and `depmap-local`
   ... Load both"). Takes `ot_score` from Open Targets and the
   `depmap_selectivity()` dict from DepMap as direct arguments, so it cannot
   run before them. **`target-triage-public-data` does not call this skill —
   its `SKILL.md` never mentions `depmap-fusion`.** The triage skill reaches
   `depmap-local` and Open Targets directly and produces its own verdict from
   them; `depmap-fusion` is a separate, alternate way of combining the same
   two evidence sources into a fused verdict, not a step the triage skill
   invokes. Use one or the other, or both, but running
   `target-triage-public-data` does not run `depmap-fusion` for you.
4. **`marker-contrast-null`** — conditional, not part of the core stack. Run
   only when step 1 or step 3 surfaces a genotype/subgroup-restricted
   dependency claim ("only matters in mutation X") that needs checking against
   an empirical null before it's reported.
5. **`depmap-genetea`** — conditional, separate axis. Its own scope note says
   per-gene dependency is `depmap-local`'s and evidence joins are
   `depmap-fusion`'s; this skill answers a different question — what a gene's
   codependencies are, or what a hit list has in common — run after or
   independent of the target verdict, not as part of reaching it.
6. **`target-lit-brief`** — recent literature brief for one target. Useful
   alongside or after the quantitative stack: retrieves and screens recent
   PubMed papers, summarises them under a target-triage lens, and produces a
   themed synthesis. Independent of DepMap/Open Targets data; can be run at any
   point when a literature snapshot is needed.

## Molecules and structures

A second, separate pipeline, used once a target is chosen and the question
shifts to which molecules bind it:

1. **`get-protein-structure`** — turns a gene symbol/accession into a cleaned
   structure file plus a docking box.
2. **`boltz-affinity-triage`** — takes that structure and box, screens
   candidate compounds against it with calibrated controls (known actives and
   property-matched decoys) rather than raw co-folding scores.

## Standalone

Not part of either pipeline above:

- **`ai-origination-audit`** — audits what a company claims about AI's role in
  discovering a candidate against the public record. Unrelated question
  (evidence provenance, not target or molecule validation).
- **`code-review`, `doc-style`, `coding-standards`** — conventions, not
  analysis. Load as needed while writing code or docs, not as pipeline steps.
- **`morning`** — personal daily routine, runs only on explicit request.

## Summary

Not all skills run together in one order. Run the target-judgment stack
(`target-triage-public-data` → `depmap-local` / `opentargets-evidence` →
`depmap-fusion`, with `marker-contrast-null`, `depmap-genetea`, and
`target-lit-brief` as conditional or complementary add-ons) to decide whether
a gene is worth pursuing. Run `get-protein-structure` → `boltz-affinity-triage`
afterward, separately, to rank molecules against a chosen target. The remaining
skills are standalone: loaded individually for their specific purpose, not
sequenced with the others.
```