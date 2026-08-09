"""Coverage for the results/<topic>/<run>/ writer added to ai-origination-audit.

Uses a tmp_path root and synthetic schema-shaped rows rather than a live
ClinicalTrials.gov fetch, so this needs no network access.
"""
import os

import pytest

import kernel as k


def candidate_row(tier="2", **overrides):
    row = {c: "" for c in k.CANDIDATE_COLUMNS}
    row.update({"candidate_id": "CAND-001", "primary_name": "CandidateX",
               "tier": tier,
               "tier_rationale": "Molecule was AI-optimized against a "
                                 "pre-validated target."})
    row.update(overrides)
    return row


def evidence_row(**overrides):
    row = {c: "" for c in k.EVIDENCE_COLUMNS}
    row.update({"evidence_id": "E1", "candidate_id": "CAND-001",
               "claim_step": "molecule_generation", "evidence_strength": "A",
               "supports_tier": True, "source_type": "peer_reviewed"})
    row.update(overrides)
    return row


def trial_row(**overrides):
    row = {c: "" for c in k.TRIAL_COLUMNS}
    row.update({"candidate_id": "CAND-001", "trial_id": "NCT00000001",
               "results_posted": "no", "registry_mentions_ai": False,
               "verification_status": "verified"})
    row.update(overrides)
    return row


def test_results_root_raises_when_unconfigured(monkeypatch):
    monkeypatch.delenv("SCIENCE_RESULTS_ROOT", raising=False)
    with pytest.raises(FileNotFoundError, match="SCIENCE_RESULTS_ROOT"):
        k.results_root()


def test_run_dir_raises_and_creates_nothing_when_root_unconfigured(monkeypatch, tmp_path):
    monkeypatch.delenv("SCIENCE_RESULTS_ROOT", raising=False)
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError):
        k.audit_run_dir("CandidateX")
    assert not (tmp_path / "results").exists()


def test_run_dir_honours_env_var(monkeypatch, tmp_path):
    monkeypatch.setenv("SCIENCE_RESULTS_ROOT", str(tmp_path))
    out_dir = k.audit_run_dir("CandidateX")
    assert out_dir == str(tmp_path / "ai_origination_audit" / "candidatex")


def test_write_run_needs_tier_before_writing_anything(tmp_path):
    out_dir = k.audit_run_dir("BadCase", root=str(tmp_path))
    with pytest.raises(ValueError):
        k.audit_write_run(out_dir, "BadCase", summary="x",
                          candidate=candidate_row(tier=""))
    written = []
    for root, _dirs, files in os.walk(out_dir):
        written.extend(files)
    assert not any(f.endswith(".csv") for f in written)


def test_write_run_happy_path(tmp_path):
    out_dir = k.audit_run_dir("CandidateX", root=str(tmp_path))
    written = k.audit_write_run(
        out_dir, "CandidateX",
        summary="CandidateX's molecule was AI-optimized; target ID was conventional.",
        candidate=candidate_row(), evidence_rows=[evidence_row()],
        trial_rows=[trial_row()], data_sources=["ClinicalTrials.gov API v2"])
    assert os.path.isfile(written["candidate"])
    assert os.path.isfile(written["trials"])
    assert os.path.isfile(written["evidence"])
    assert os.path.isfile(written["readme"])
    text = open(written["readme"]).read()
    assert "Tier 2" in text and "## Limits" in text


def test_human_check_queue_only_reads_trial_rows():
    needs = k.audit_human_check_queue(
        [trial_row(verification_status="needs_verification")])
    assert len(needs) == 1
