---
name: ai-origination-audit
description: Audit how much of a drug candidate's discovery was actually done by AI, and produce a tiered evidence map with provenance. Use when asked to verify, tier, or systematically map AI-designed / AI-discovered / AI-originated therapeutic candidates, to check company AI claims against clinical trial registries, or to build a candidate/trial/evidence table for AI-driven drug discovery programs.
---

# AI-origination audit

Separate what a company **claims** about AI's role in a drug program from what
the record **shows**. Output is three schema-locked tables plus a report, so
candidates audited months apart stay comparable and mergeable.

The core discipline: an AI-attribution claim is a claim about a *step*
(target identification, molecule generation, lead optimization, trial design),
made by a *source* of a given strength. Never record "AI-designed" as a
property of a candidate — record which step, from which source, at which grade.
Tier follows from the grid, and is defensible because of it.

`kernel.py` loads with this skill: `fetch_ctgov_study`, `search_ctgov`,
`extract_trial_row`, `scan_registry_for_ai`, `posted_results_summary`,
`blank_row`, `provenance_ref`, `audit_summary`, and the constants
`TIER_RUBRIC`, `EVIDENCE_GRADES`, `CLAIM_STEPS`, `AI_TERMS`,
`CANDIDATE_COLUMNS`, `TRIAL_COLUMNS`, `EVIDENCE_COLUMNS`.

## Tiers

    1  = both target and molecule primarily identified/designed by AI
    2  = AI primarily designed or optimized the molecule (target pre-validated)
    3  = significant but partial AI contribution
    ex = computational but not AI (FEP / structure-based only), or no AI role

Say "Tier N, qualified" and state the qualification whenever the tier rests on
different evidence at different steps. A candidate whose target was AI-ranked
but whose clinical structure emerged from human ADME optimization is Tier 1
under the rubric above and Tier 2 under a stricter one — write both readings
down rather than picking silently.

## Evidence grades

    A = peer-reviewed, methods-level (platform named, step described, params given)
    B = peer-reviewed assertion without methods detail (incl. conference abstracts)
    C = first-party company press release / website
    D = third-party media or unverified secondary source only

**Downgrade rule:** if no A or B evidence backs the AI-role claim, the tier is
press-release-only — say so explicitly and lower confidence.
`audit_summary()` returns `downgrade_recommended` for exactly this.

## Workflow

1. **Enumerate every trial, not just the ones you were given.** Run
   `search_ctgov(intervention=...)` for each synonym AND
   `search_ctgov(sponsor=...)`, then reconcile. Sponsor sweeps surface
   sibling molecules — exclude them by scope, but note the count so the
   reader knows what was screened.
2. **Pull full records and flatten.** `fetch_ctgov_study(nct)` then
   `extract_trial_row(study)`. Never work from search-result summaries: they
   omit endpoints, site counts, and the results-posting dates.
3. **Scan the registry for AI language.** `scan_registry_for_ai(study)`. An
   empty result across every trial is a *finding*, and a common one: the
   attribution exists only in literature and marketing, never in the
   regulatory record. Report it.
4. **Extract posted results where they exist.** `posted_results_summary(study)`
   returns primary and secondary outcome values, LS-mean analyses, participant
   flow, and adverse-event counts. Reconcile these against the publication and
   the press release digit by digit, and report the per-arm n actually analyzed
   — efficacy signals in early trials often rest on far fewer patients than
   were randomized.
5. **Cover the registries that are not ClinicalTrials.gov.** Programs routinely
   run phase 0/1 in ANZCTR (`ACTRN…`) or China CDE (`CTR…`), and Chinese
   phase 3 studies often carry a CDE ID paired with the NCT. Use `blank_row`
   for these and mark `verification_status="needs_verification"` unless you
   queried that registry directly. A count based on ClinicalTrials.gov alone
   undercounts.
6. **Sweep beyond PubMed.** Europe PMC and OpenAlex index conference abstracts
   (ATS, ERS, CHEST), book chapters, and DOI-registered commentary that PubMed
   misses — and those are frequently where the sharpest attribution wording and
   the discrepant numbers live. Reconstruct OpenAlex abstracts from
   `abstract_inverted_index` by sorting positions.
7. **Read the medicinal-chemistry paper specifically.** When a company
   publishes both a platform paper and a dedicated med-chem paper on the same
   molecule, compare them. A med-chem paper that describes only structure-based
   design while the platform paper claims generative AI is the most
   informative single discrepancy available, precisely because both are
   first-party and peer-reviewed.
8. **Build the three tables**, then `audit_summary(evidence_rows, trial_rows)`.
9. **Write the report** with methodology, tier rationale citing evidence IDs,
   discrepancy flags, a human-check queue, and limitations.
10. **Save the run** — `audit_write_run(out_dir, name, summary, candidate=...,
    evidence_rows=..., trial_rows=...)`, not assembled by hand, into
    `results/ai_origination_audit/<name>/`. `out_dir` comes from
    `audit_run_dir(name)`. Fails before writing anything if `candidate` has no
    `tier` set.

```python
out_dir = audit_run_dir("CandidateX")
audit_write_run(out_dir, "CandidateX", summary="...",
                candidate=candidate_row, evidence_rows=evidence_rows,
                trial_rows=trial_rows,
                data_sources=["ClinicalTrials.gov API v2", "PubMed", "company press releases"])
```

## Discrepancy classes to check every time

- **Announcement vs registry status** — "initiates Phase III" while the record
  reads NOT_YET_RECRUITING with a later estimated start.
- **Primary vs secondary endpoint framing** — the quoted efficacy number is
  often a secondary endpoint while the registry primary is a safety measure.
- **Publication before results posting** — note the gap in months.
- **Unstable comparators** — congress abstracts sometimes exclude outliers,
  shifting the placebo arm while the treatment value stays fixed. Always
  compare the *control* arm across sources, not just the treatment arm.
- **Missing results postings** on completed trials.
- **Stale records** — status RECRUITING with a primary completion date in
  the past.
- **Phase or indication inflation** relative to the registry.
- **Unregistered programs** promoted as clinical milestones.
- **Identifier drift** — one compound under several codes, including
  misspellings in published abstract titles. Search every variant; an audit
  keyed on one code undercounts.

## Provenance

Use `provenance_ref` for a uniform, auditable format:
`REG:CTGOV:NCT…@<date>`, `PMID:…|DOI:…`, `PREPRINT:<server>:<id>`,
`CO:<domain>:<slug>@<date>`, `MEDIA:<domain>@<date>`. Every claim in the
report resolves to one of these or to an evidence row that does. Mark anything
you did not verify against a machine-readable source as
`needs_verification` and put it in the human-check queue — do not let it
enter a headline count silently.

## Reporting

Lead with the tier and its qualification, then the trial table, then the
discrepancies. State what the evidence establishes and what it does not: an
unusually specific methods disclosure is evidence of disclosure quality, not
of counterfactual causation. No public source establishes what a conventional
program would have found for the same target — say so rather than implying
the AI claim has been independently validated.
