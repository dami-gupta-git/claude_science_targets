"""Helpers for AI-origination evidence audits of clinical-stage drug candidates.

Loaded automatically when the ai-origination-audit skill is loaded.
See SKILL.md for the workflow these functions support.
"""
import csv
import json
import os
import re
import shutil
import urllib.parse
import urllib.request

CTGOV_API = "https://clinicaltrials.gov/api/v2/studies"

AI_TERMS = (
    "artificial intelligence", "machine learning", "deep learning", "generative",
    "neural network", "transformer", "in silico", "de novo design",
    "PandaOmics", "Chemistry42", "Pharma.AI", "inClinico", "AlphaFold",
)

TIER_RUBRIC = (
    "1  = both target and molecule primarily identified/designed by AI",
    "2  = AI primarily designed or optimized the molecule (target pre-validated)",
    "3  = significant but partial AI contribution",
    "ex = computational but not AI (e.g. FEP/structure-based only), or no AI role",
)

EVIDENCE_GRADES = (
    "A = peer-reviewed source with methods-level detail (platform named, step described)",
    "B = peer-reviewed assertion without methods detail (incl. conference abstract)",
    "C = first-party company press release / website claim",
    "D = third-party media or unverified secondary source only",
)

CANDIDATE_COLUMNS = (
    "candidate_id", "primary_name", "synonyms", "company", "company_type", "modality",
    "target", "target_novelty", "indication", "highest_phase_reached", "highest_phase_status",
    "first_in_human_year", "tier", "tier_rationale", "ai_role_target_id",
    "ai_role_molecule_design", "ai_platforms_named", "evidence_strength",
    "peer_reviewed_ai_claim", "registry_ids", "n_trials", "results_posted_any",
    "discrepancy_flags", "verification_status", "provenance_refs",
)

TRIAL_COLUMNS = (
    "candidate_id", "trial_id", "registry", "org_study_id", "title", "phase", "design",
    "sponsor", "collaborators", "status", "status_as_of", "start_date",
    "primary_completion_date", "enrollment", "enrollment_type", "countries", "n_sites",
    "population", "min_age", "primary_endpoint", "key_secondary_endpoints",
    "results_posted", "results_first_posted", "registry_mentions_ai",
    "discrepancy_flags", "verification_status", "provenance_ref",
)

EVIDENCE_COLUMNS = (
    "evidence_id", "candidate_id", "claim_step", "claim_text", "claim_type", "source_type",
    "source_title", "source_id", "source_date", "first_party", "peer_reviewed",
    "evidence_strength", "supports_tier", "notes",
)

CLAIM_STEPS = (
    "target_identification", "molecule_generation", "lead_optimization",
    "trial_design", "biomarker/patient_selection", "timeline_claim",
)


def provenance_ref(kind, ident, date=None, extra=None):
    """Format a provenance reference in the audit's fixed convention.

    kind: 'registry' | 'pubmed' | 'preprint' | 'company' | 'media'
    """
    if kind == "registry":
        reg = extra or "CTGOV"
        return "REG:%s:%s@%s" % (reg, ident, date or "unknown")
    if kind == "pubmed":
        return "PMID:%s|DOI:%s" % (ident, extra or "none")
    if kind == "preprint":
        return "PREPRINT:%s:%s" % (extra or "unknown", ident)
    if kind == "company":
        return "CO:%s:%s@%s" % (extra or "unknown", ident, date or "unknown")
    if kind == "media":
        return "MEDIA:%s@%s" % (ident, date or "unknown")
    raise ValueError("unknown provenance kind: %r" % (kind,))


def fetch_ctgov_study(nct_id):
    """Return the full ClinicalTrials.gov API v2 record for one NCT ID."""
    url = "%s/%s" % (CTGOV_API, nct_id)
    with urllib.request.urlopen(url, timeout=60) as fh:
        return json.load(fh)


def search_ctgov(intervention=None, sponsor=None, page_size=100):
    """Return every API v2 study record matching an intervention and/or sponsor term.

    Follows nextPageToken until the result set is exhausted, so callers never
    see a silently truncated sweep.
    """
    params = {"pageSize": str(page_size)}
    if intervention:
        params["query.intr"] = intervention
    if sponsor:
        params["query.spons"] = sponsor
    studies = []
    page_token = None
    while True:
        page_params = dict(params)
        if page_token:
            page_params["pageToken"] = page_token
        url = "%s?%s" % (CTGOV_API, urllib.parse.urlencode(page_params))
        with urllib.request.urlopen(url, timeout=60) as fh:
            payload = json.load(fh)
        studies.extend(payload.get("studies", []))
        page_token = payload.get("nextPageToken")
        if not page_token:
            break
    return studies


def extract_trial_row(study, candidate_id="CAND-001", retrieval_date=None):
    """Flatten one API v2 record into the audit's trial-table schema."""
    ps = study.get("protocolSection", {})
    idm = ps.get("identificationModule", {})
    sm = ps.get("statusModule", {})
    dm = ps.get("designModule", {})
    elig = ps.get("eligibilityModule", {})
    om = ps.get("outcomesModule", {})
    di = dm.get("designInfo", {})
    locs = ps.get("contactsLocationsModule", {}).get("locations", []) or []
    enroll = dm.get("enrollmentInfo") or {}
    nct = idm.get("nctId")
    countries = sorted({loc.get("country") for loc in locs if loc.get("country")})
    masking = (di.get("maskingInfo") or {}).get("masking")
    design_bits = [di.get("allocation"), di.get("interventionModel"),
                   di.get("primaryPurpose"), masking]
    row = {
        "candidate_id": candidate_id,
        "trial_id": nct,
        "registry": "ClinicalTrials.gov",
        "org_study_id": (idm.get("orgStudyIdInfo") or {}).get("id"),
        "title": idm.get("officialTitle"),
        "phase": "/".join(dm.get("phases", []) or []),
        "design": "; ".join([b for b in design_bits if b]),
        "sponsor": (ps.get("sponsorCollaboratorsModule", {})
                    .get("leadSponsor", {}).get("name")),
        "collaborators": "; ".join(
            c.get("name") for c in
            (ps.get("sponsorCollaboratorsModule", {}).get("collaborators") or [])
        ) or "none listed",
        "status": sm.get("overallStatus"),
        "status_as_of": (sm.get("lastUpdatePostDateStruct") or {}).get("date"),
        "start_date": (sm.get("startDateStruct") or {}).get("date"),
        "primary_completion_date": (sm.get("primaryCompletionDateStruct") or {}).get("date"),
        "enrollment": enroll.get("count"),
        "enrollment_type": enroll.get("type"),
        "countries": "; ".join(countries) or "not listed",
        "n_sites": len(locs),
        "population": ("healthy volunteers" if elig.get("healthyVolunteers")
                       else "patients"),
        "min_age": elig.get("minimumAge"),
        "primary_endpoint": " | ".join(
            o.get("measure", "") for o in (om.get("primaryOutcomes") or [])),
        "key_secondary_endpoints": " | ".join(
            o.get("measure", "") for o in (om.get("secondaryOutcomes") or [])[:4]),
        "results_posted": "yes" if study.get("hasResults") else "no",
        "results_first_posted": (sm.get("resultsFirstPostDateStruct") or {}).get("date") or "",
        "registry_mentions_ai": bool(scan_registry_for_ai(study)),
        "discrepancy_flags": "",
        "verification_status": "verified",
        "provenance_ref": provenance_ref("registry", nct, retrieval_date),
    }
    return row


def scan_registry_for_ai(study, terms=None):
    """Return AI terms present in a study's title/description modules.

    An empty result is a finding in its own right: the AI attribution exists
    only outside the regulatory record.
    """
    if terms is None:
        terms = AI_TERMS
    ps = study.get("protocolSection", {})
    blob = json.dumps(ps.get("descriptionModule", {}))
    blob += json.dumps(ps.get("identificationModule", {}))
    blob += json.dumps(ps.get("armsInterventionsModule", {}))
    hits = []
    for term in terms:
        if re.search(r"\b%s\b" % re.escape(term), blob, re.I):
            hits.append(term)
    if re.search(r"\bAI\b", blob):
        hits.append("AI")
    return sorted(set(hits))


def posted_results_summary(study):
    """Extract posted primary/secondary outcome values from an API v2 record.

    Returns {} when no results are posted. Use to reconcile registry numbers
    against publication and press figures.
    """
    rs = study.get("resultsSection")
    if not rs:
        return {}
    out = {"outcomes": [], "participant_flow": {}, "adverse_events": []}
    pf = rs.get("participantFlowModule", {})
    gmap = {g["id"]: g.get("title") for g in pf.get("groups", [])}
    for period in pf.get("periods", []):
        for milestone in period.get("milestones", []):
            out["participant_flow"][milestone.get("type")] = {
                gmap.get(a.get("groupId"), a.get("groupId")): a.get("numSubjects")
                for a in milestone.get("achievements", [])
            }
    for om in rs.get("outcomeMeasuresModule", {}).get("outcomeMeasures", []):
        ogmap = {g["id"]: g.get("title") for g in om.get("groups", [])}
        classes = []
        for cls in om.get("classes", []):
            for cat in cls.get("categories", []):
                classes.append({
                    "class": cls.get("title"),
                    "values": {ogmap.get(m.get("groupId"), m.get("groupId")): m.get("value")
                               for m in cat.get("measurements", [])},
                })
        analyses = [{
            "groups": [ogmap.get(g, g) for g in a.get("groupIds", [])],
            "param": a.get("paramType"),
            "value": a.get("paramValue"),
            "ci_low": a.get("ciLowerLimit"),
            "ci_high": a.get("ciUpperLimit"),
        } for a in om.get("analyses", []) or []]
        out["outcomes"].append({
            "type": om.get("type"),
            "title": om.get("title"),
            "unit": om.get("unitOfMeasure"),
            "param_type": om.get("paramType"),
            "classes": classes,
            "analyses": analyses,
        })
    for eg in rs.get("adverseEventsModule", {}).get("eventGroups", []) or []:
        out["adverse_events"].append({
            "group": eg.get("title"),
            "at_risk": eg.get("seriousNumAtRisk"),
            "serious": eg.get("seriousNumAffected"),
            "other": eg.get("otherNumAffected"),
            "deaths": eg.get("deathsNumAffected"),
        })
    return out


def blank_row(columns, **values):
    """Build a schema-complete row dict, empty strings for anything unset.

    Use for trials that live in non-ClinicalTrials.gov registries (CDE, ANZCTR,
    ChiCTR, jRCT, CTIS) or that have no locatable registration at all.
    """
    row = {c: "" for c in columns}
    unknown = [k for k in values if k not in row]
    if unknown:
        raise ValueError("not in schema: %r" % (sorted(unknown),))
    row.update(values)
    return row


def audit_summary(evidence_rows, trial_rows):
    """Tally evidence grades, source mix, and verification coverage.

    Flags the plan's downgrade condition: AI-role evidence that is
    press-release-only (no A or B grade behind it).
    """
    grades = {}
    sources = {}
    steps = {}
    for row in evidence_rows:
        grades[row.get("evidence_strength")] = grades.get(row.get("evidence_strength"), 0) + 1
        sources[row.get("source_type")] = sources.get(row.get("source_type"), 0) + 1
        steps[row.get("claim_step")] = steps.get(row.get("claim_step"), 0) + 1
    verif = {}
    posted = 0
    ai_in_registry = 0
    for row in trial_rows:
        verif[row.get("verification_status")] = verif.get(row.get("verification_status"), 0) + 1
        if row.get("results_posted") == "yes":
            posted += 1
        if row.get("registry_mentions_ai") is True:
            ai_in_registry += 1
    peer = sum(
        1 for row in evidence_rows
        if row.get("supports_tier") and row.get("evidence_strength") in ("A", "B")
    )
    return {
        "n_evidence": len(evidence_rows),
        "grades": grades,
        "source_types": sources,
        "claim_steps": steps,
        "n_trials": len(trial_rows),
        "verification": verif,
        "trials_with_posted_results": posted,
        "trials_whose_registry_mentions_ai": ai_in_registry,
        "press_release_only": peer == 0,
        "downgrade_recommended": peer == 0,
    }


# --------------------------------------------------------------------------
# Run output: results/<topic>/<run>/
# --------------------------------------------------------------------------
RESULTS_TOPIC = "ai_origination_audit"
CANDIDATE_CSV = "candidate.csv"
TRIALS_CSV = "trials.csv"
EVIDENCE_CSV = "evidence.csv"
RUN_README = "README.md"
SUMMARY_MAX_WORDS = 130

# Restated from this skill's SKILL.md so every run states the same standing
# caveats rather than a hand-picked subset that can drift from the skill as
# it is revised.
STANDING_LIMITS = (
    "An unusually specific methods disclosure is evidence of disclosure "
    "quality, not of counterfactual causation. No public source establishes "
    "what a conventional program would have found for the same target.",
    "Coverage stops at what was actually swept: sponsor and synonym sweeps "
    "of ClinicalTrials.gov, non-ClinicalTrials.gov registries queried "
    "directly, and the literature/press sources listed under Data sources. "
    "A registry not queried directly is not represented, even where a trial "
    "in it is listed via `blank_row`.",
    "An empty AI-language scan of the registry record is not evidence of no "
    "AI role — it means the attribution, if any, lives outside the "
    "regulatory record.",
    "Tier follows the evidence grid in SKILL.md; a candidate whose evidence "
    "differs by step is reported as 'Tier N, qualified', not silently "
    "rounded to one tier.",
)


def audit_run_dir(name, root="results", topic=None, make=True):
    """Path for one audit run: <root>/<topic>/<slug>/, with scripts/ beside it.

    `name` is the candidate or company audited; slugged to snake_case because
    result directories are named that way and a free-text name may carry
    spaces, hyphens or case.

    `topic` defaults to RESULTS_TOPIC through an explicit None check rather
    than in the signature: the kernel.py sidecar loader rejects a non-literal
    default, and a rejected file defines none of this module's helpers.
    """
    topic = RESULTS_TOPIC if topic is None else topic
    slug = re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")
    if not slug:
        raise ValueError("audit_run_dir got an empty name; the run directory "
                         "is named after the candidate or company audited")
    out_dir = os.path.join(root, topic, slug)
    if make:
        os.makedirs(os.path.join(out_dir, "scripts"), exist_ok=True)
    return out_dir


def audit_write_table(path, rows, headers=None):
    """Write rows (list of dicts) to CSV. Returns the path, or None if empty.

    Columns are `headers` first (the schema-locked column order for that
    table), then any further keys present in the rows, so a field a caller
    attaches beyond the fixed schema is written rather than dropped.
    """
    rows = [dict(r) for r in (rows or [])]
    if not rows:
        return None
    headers = list(headers or [])
    extra = sorted({k for r in rows for k in r} - set(headers))
    cols = [c for c in headers if any(c in r for r in rows)] + extra
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({c: ("" if row.get(c) is None else row.get(c))
                             for c in cols})
    return path


def audit_check_words(text, cap, label):
    """Enforce a word cap on run-README prose, naming the count and the cap."""
    n = len(str(text or "").split())
    if n > cap:
        raise ValueError(
            f"{label} is {n} words; cap is {cap}. Move detail into the table "
            "rather than the prose.")
    return n


def audit_human_check_queue(trial_rows):
    """Trial rows flagged `verification_status='needs_verification'`.

    The report's human-check queue: anything not verified against a
    machine-readable source, so it does not enter a headline count silently.
    Trial rows only — `EVIDENCE_COLUMNS` carries no `verification_status`
    field, since an evidence row's strength (its `evidence_strength` grade)
    is already how much it is trusted; there is no separate verified/unverified
    state to flag on it.
    """
    return [r for r in trial_rows if r.get("verification_status") == "needs_verification"]


def audit_run_readme(name, summary, candidate=None, evidence_rows=(),
                     trial_rows=(), files=(), data_sources=(), limits=(),
                     title=None):
    """Render the run README for one AI-origination audit. Returns markdown.

    `candidate` is a CANDIDATE_COLUMNS-shaped row carrying the tier and its
    rationale; `evidence_rows`/`trial_rows` are the other two schema-locked
    tables. Sections follow `coding-standards` (Result, Files, Data sources,
    Limits) and SKILL.md's own reporting order: tier and qualification first,
    then the trial picture, then discrepancies and the human-check queue.

    Raises if `candidate` has no `tier` — a report with no tier is not the
    audit's actual deliverable, whatever else it collected.
    """
    if not candidate or not candidate.get("tier"):
        raise ValueError("audit_run_readme needs a candidate row with a "
                         "'tier' set; a report with no tier assessed is not "
                         "the audit's result")
    audit_check_words(summary, SUMMARY_MAX_WORDS, "summary")
    if str(summary or "").strip().startswith("-"):
        raise ValueError("summary must be plain prose, not bullets — it is "
                         "the paragraph a non-specialist reads first")

    evidence_rows = list(evidence_rows)
    trial_rows = list(trial_rows)
    stats = audit_summary(evidence_rows, trial_rows)
    trial_needs = audit_human_check_queue(trial_rows)

    default_title = f"{candidate.get('primary_name') or name} — AI-origination audit"
    parts = [f"# {title or default_title}", "", str(summary).strip(), "",
             "## Result", ""]

    tier = candidate.get("tier")
    rationale = candidate.get("tier_rationale") or ""
    parts += [f"**Tier {tier}.**" + (f" {rationale}" if rationale else ""), ""]
    if stats["downgrade_recommended"]:
        parts += ["No A- or B-grade evidence backs the AI-role claim at any "
                 "step: this tier is press-release-only and confidence is "
                 "lowered accordingly.", ""]

    parts += [
        "| quantity | value |", "| --- | --- |",
        f"| trials covered | {stats['n_trials']} |",
        f"| trials with posted results | {stats['trials_with_posted_results']} |",
        f"| trials whose registry mentions AI language | "
        f"{stats['trials_whose_registry_mentions_ai']} |",
        f"| evidence rows | {stats['n_evidence']} |",
        f"| evidence grades | "
        f"{', '.join(f'{k}={v}' for k, v in sorted(stats['grades'].items()) if k is not None)} |",
        ""]

    discrepant = [r for r in trial_rows if r.get("discrepancy_flags")]
    if discrepant:
        parts += ["### Discrepancies", ""]
        for r in discrepant:
            parts.append(f"- `{r.get('trial_id')}`: {r['discrepancy_flags']}")
        parts += [""]

    if trial_needs or evidence_needs:
        parts += ["### Human-check queue", "",
                 f"{len(trial_needs)} trial row(s) and {len(evidence_needs)} "
                 "evidence row(s) are marked `needs_verification` and were "
                 "not queried against a machine-readable source; excluded "
                 "from headline counts above.", ""]

    parts += ["## Files", ""]
    for entry in files:
        f_name, description = (entry["name"], entry["description"]) \
            if isinstance(entry, dict) else entry
        parts.append(f"- `{f_name}` — {description}")
    parts += ["", "## Data sources", ""] + [f"- {s}" for s in data_sources]
    parts += ["", "## Limits", ""]
    parts += [f"- {s}" for s in list(STANDING_LIMITS) + list(limits)]
    return "\n".join(parts).rstrip() + "\n"


def audit_write_run(out_dir, name, summary, candidate=None, evidence_rows=(),
                    trial_rows=(), files=(), data_sources=(), limits=(),
                    scripts=()):
    """Write one audit run directory. Returns {name: path} for what was written.

    Writes the three schema-locked tables (`candidate.csv`, `trials.csv`,
    `evidence.csv`), the run README, and copies `scripts` into `scripts/`.
    `candidate` is validated (via `audit_run_readme`) before anything is
    written, so a run missing its tier fails before leaving a half-written
    directory.
    """
    candidate_row = [candidate] if candidate else []
    trial_rows = list(trial_rows)
    evidence_rows = list(evidence_rows)

    # Rendered up front, before any file touches disk, so a candidate missing
    # its tier raises here rather than after candidate.csv/trials.csv/
    # evidence.csv are already written -- a half-written directory is exactly
    # what this validation exists to prevent.
    readme = audit_run_readme(name, summary, candidate=candidate,
                              evidence_rows=evidence_rows,
                              trial_rows=trial_rows, files=(),
                              data_sources=data_sources, limits=limits)

    os.makedirs(os.path.join(out_dir, "scripts"), exist_ok=True)
    written = {}

    cand_path = audit_write_table(os.path.join(out_dir, CANDIDATE_CSV),
                                  candidate_row, CANDIDATE_COLUMNS)
    if cand_path:
        written["candidate"] = cand_path

    trials_path = audit_write_table(os.path.join(out_dir, TRIALS_CSV),
                                    trial_rows, TRIAL_COLUMNS)
    if trials_path:
        written["trials"] = trials_path

    evidence_path = audit_write_table(os.path.join(out_dir, EVIDENCE_CSV),
                                      evidence_rows, EVIDENCE_COLUMNS)
    if evidence_path:
        written["evidence"] = evidence_path

    for src in scripts:
        dst = os.path.join(out_dir, "scripts", os.path.basename(src))
        shutil.copyfile(src, dst)
        written.setdefault("scripts", []).append(dst)

    listed = []
    if "candidate" in written:
        listed.append((CANDIDATE_CSV, "one row: the candidate's tier, "
                                      "AI-role breakdown and evidence strength"))
    if "trials" in written:
        listed.append((TRIALS_CSV, "one row per trial: registry status, "
                                   "results-posting dates, discrepancy flags"))
    if "evidence" in written:
        listed.append((EVIDENCE_CSV, "one row per claim: which step, from "
                                     "which source, at which evidence grade"))
    listed += list(files)

    readme = audit_run_readme(name, summary, candidate=candidate,
                              evidence_rows=evidence_rows,
                              trial_rows=trial_rows, files=listed,
                              data_sources=data_sources, limits=limits)
    readme_path = os.path.join(out_dir, RUN_README)
    with open(readme_path, "w", encoding="utf-8") as fh:
        fh.write(readme)
    written["readme"] = readme_path
    return written
