# ai-origination-audit

Audits how much of a drug candidate's discovery was performed by AI, separating
what a company claims from what the public record shows. An AI-attribution
claim is treated as a claim about a specific discovery step — target
identification, molecule generation, lead optimization, trial design — made by
a source of a given strength, rather than as a property of the candidate. The
output is three schema-locked tables (candidate, trial, evidence) plus a
report, so candidates audited months apart remain comparable and mergeable.
Tier follows from the claim-by-source grid, which is what makes it possible to
say which evidence a tier rests on.

`kernel.py` loads automatically with the skill and holds the registry access,
flattening and tallying helpers. `SKILL.md` holds the workflow, the tier
rubric, the evidence grades, and the list of discrepancy classes to check on
every audit.

## Functions

- `search_ctgov(intervention=None, sponsor=None, page_size=100)` — API v2 study
  records matching an intervention term, a sponsor, or both. Sponsor sweeps
  surface sibling molecules that an intervention-only search misses.
- `fetch_ctgov_study(nct_id)` — the full record for one NCT ID. Search results
  omit endpoints, site counts and results-posting dates, so trial rows are
  built from full records rather than summaries.
- `extract_trial_row(study, candidate_id="CAND-001", retrieval_date=None)` —
  flattens one API v2 record into `TRIAL_COLUMNS`.
- `scan_registry_for_ai(study, terms=None)` — the AI terms present in a study's
  title and description modules; defaults to `AI_TERMS`. An empty result across
  every trial of a program is itself a finding, and is reported as one.
- `posted_results_summary(study)` — primary and secondary outcome values,
  LS-mean analyses, participant flow and adverse-event counts; returns `{}`
  when no results are posted. Used to reconcile registry figures against the
  publication and the press release.
- `blank_row(columns, **values)` — a schema-complete row with empty strings for
  anything unset, for trials in registries other than ClinicalTrials.gov (CDE,
  ANZCTR, ChiCTR, jRCT, CTIS) or with no locatable registration.
- `provenance_ref(kind, ident, date=None, extra=None)` — one reference in the
  audit's fixed format, `kind` being one of `registry`, `pubmed`, `preprint`,
  `company`, `media`.
- `audit_summary(evidence_rows, trial_rows)` — tallies evidence grades, source
  mix and verification coverage, and returns `downgrade_recommended` when the
  AI-role evidence is press-release-only.

Module constants: `TIER_RUBRIC`, `EVIDENCE_GRADES`, `CLAIM_STEPS`, `AI_TERMS`,
`CANDIDATE_COLUMNS`, `TRIAL_COLUMNS`, `EVIDENCE_COLUMNS`, `CTGOV_API`.

### Saving a run

- `audit_run_dir(name, root="results", topic="ai_origination_audit")` — the
  path for one run, `<root>/<topic>/<slug>/`, with `scripts/` created beside
  it. `name` is the candidate or company audited.
- `audit_write_run(out_dir, name, summary, candidate=None, evidence_rows=(),
  trial_rows=(), files=(), data_sources=(), limits=(), scripts=())` — writes
  the three schema-locked tables (`candidate.csv`, `trials.csv`,
  `evidence.csv`), the run README, and copies `scripts` into `scripts/`.
  Returns the paths written. Fails before writing anything if `candidate` has
  no `tier` set.
- `audit_run_readme(name, summary, candidate=None, evidence_rows=(),
  trial_rows=(), files=(), data_sources=(), limits=(), title=None)` — renders
  the README text `audit_write_run` saves, following SKILL.md's own reporting
  order: tier and qualification first (with the downgrade note when
  `audit_summary()` recommends it), then the trial picture, then
  discrepancies and the human-check queue. The Limits section always
  includes this skill's standing caveats (methods disclosure is not
  counterfactual evidence, coverage stops at what was actually swept, an
  empty AI-language scan is not evidence of no AI role, a qualified tier is
  reported as such) plus any run-specific ones passed in.
- `audit_human_check_queue(trial_rows, evidence_rows)` — the rows either
  table marked `needs_verification`, as two lists.
- `audit_write_table(path, rows, headers=None)` — CSV writer for any of the
  three schema-locked tables, used by `audit_write_run` but usable standalone.

Runs land in `results/ai_origination_audit/<slug>/`.

## Grades and the downgrade rule

The four evidence grades run A (peer-reviewed with methods-level detail —
platform named, step described) through B (peer-reviewed assertion without
methods detail, including conference abstracts), C (first-party company press
release or website) to D (third-party media or unverified secondary source
only). These are editorial definitions set in `SKILL.md`, not values fitted to
data; the grade of a given source is re-derivable by reading that source
against those four definitions, and changing a definition means re-grading
every evidence row that used it.

The one automated threshold is the downgrade condition: `audit_summary()`
returns `downgrade_recommended` when no row in `evidence_rows` with
`supports_tier` set has grade A or B. Grade alone isn't enough — a grade-A
med-chem paper that contradicts the platform's AI claim (see workflow step 7)
must not count toward `peer`, so only rows that actually support the assigned
tier are tallied. A tier resting only on C/D evidence, or on A/B evidence that
argues against it, is reported as press-release-only with lowered confidence.

Tiers are 1 (target and molecule both primarily AI-derived), 2 (molecule
primarily AI-designed or optimized against a pre-validated target), 3 (partial
AI contribution) and `ex` (computational but not AI, or no AI role). Where a
candidate reads as a different tier under a stricter rubric, both readings are
written down and the tier is reported as qualified.

## Verification

Rows not checked against a machine-readable source are marked
`verification_status="needs_verification"` and enter the report's human-check
queue rather than a headline count. Counts drawn from ClinicalTrials.gov alone
undercount programs that run early phases in ANZCTR or China CDE, so a
registry-scope statement belongs with any trial count.

## Scope

Evidence provenance and tiering only. The skill does not assess whether a
candidate will succeed clinically, and it does not establish counterfactual
causation: an unusually specific methods disclosure is evidence of disclosure
quality, not evidence that a conventional program would have failed on the same
target. No public source supports that comparison, and audits say so. Target
biology and human genetic evidence belong to `opentargets-evidence`;
dependency and cell-line evidence to `depmap-local` and `depmap-fusion`;
structure and binding questions to `get-protein-structure` and
`boltz-affinity-triage`. No test suite ships with this skill; the ClinicalTrials.gov
helpers require network access to `clinicaltrials.gov`.
