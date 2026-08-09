"""Cross-KB fusion: Open Targets human evidence x DepMap cell-line dependency."""
import csv
import os
import re
import shutil

VERDICTS = (
    "concordant-dependency",
    "common-essential-no-window",
    "growth-suppressive-mismatch",
    "context-restricted",
    "inert-in-panel",
    "evidence-without-dependency",
    "dependency-with-thin-evidence",
    "dependency-without-evidence",
    "no-evidence-no-dependency",
    "indeterminate",
)

# Minimum OT association score for a target to count as having ANY human
# evidence. Below this (or ot_score=None) the `evidence-without-dependency`
# label would assert evidence that was never supplied - such rows get
# `no-evidence-no-dependency` instead.
EVIDENCE_FLOOR = 0.10

# A real growth suppressor shows a POSITIVE tail across many lines; an inert
# gene sits flat near zero. Calibrated on TP53/RB1 (suppressors) vs
# ROS1/NRG1/ALK (inert in the DepMap panel) - see SKILL.md ## Calibration.
SUPPRESSOR_MIN_FRAC_POS = 0.15
SUPPRESSOR_MIN_SD = 0.20
INERT_MAX_SD = 0.20
INERT_MAX_ABS_MEAN = 0.25

# DepMap ships CRISPRInferredCommonEssentials.csv (1827 genes). It is NOT a
# frequency threshold: it includes KRAS (52% of lines dependent) and some genes
# at 0%. Our frac_dependent >= 0.90 rule yields 974 genes - a strict SUBSET.
# We keep the stricter rule for the no-window verdict, because adopting the
# broader list would discard KRAS, and surface DepMap's call alongside it as
# `depmap_inferred_essential` for cross-reference.


def fuse_target_row(symbol, ot_score, effect_stats, lineage_row=None,
                    tractable_modalities=None, dep_threshold=-0.5):
    """Classify one target by how Open Targets evidence and DepMap dependency agree.

    effect_stats: dict from depmap_selectivity().
    lineage_row: optional dict with mean_effect / frac_dependent / p_bh for the
        disease-relevant lineage.
    Returns a dict with a `verdict` plus the numbers behind it.
    """
    mean_eff = effect_stats.get("mean_effect")
    frac = effect_stats.get("frac_dependent")
    frac = 0.0 if frac is None else frac
    cls = effect_stats.get("classification")
    lin_mean = (lineage_row or {}).get("mean_effect")
    lin_frac = (lineage_row or {}).get("frac_dependent")
    lin_p = (lineage_row or {}).get("p_bh")

    strong_evidence = ot_score is not None and ot_score >= 0.5
    has_evidence = ot_score is not None and ot_score >= EVIDENCE_FLOOR
    context = (lin_mean is not None and lin_p is not None
               and lin_mean <= dep_threshold and lin_p < 0.05
               and (mean_eff is None or lin_mean < mean_eff - 0.2))

    frac_pos = effect_stats.get("frac_positive")
    sd = effect_stats.get("sd_effect")

    suppressive = (mean_eff is not None and mean_eff > 0.1
                   and frac_pos is not None and sd is not None
                   and frac_pos >= SUPPRESSOR_MIN_FRAC_POS and sd >= SUPPRESSOR_MIN_SD)
    inert = (sd is not None and mean_eff is not None
             and sd < INERT_MAX_SD and abs(mean_eff) < INERT_MAX_ABS_MEAN
             and frac < 0.05)

    if sd is None or frac_pos is None or mean_eff is None:
        # The growth-suppressive and inert tests cannot run without dispersion.
        # This branch is FIRST deliberately: reporting any other verdict here
        # would assert a distinction that was never tested.
        verdict = "indeterminate"
    elif cls == "common essential":
        # Every cell needs it - killing tumour cells means killing normal tissue.
        verdict = "common-essential-no-window"
    elif suppressive:
        verdict = "growth-suppressive-mismatch"
    elif context:
        verdict = "context-restricted"
    elif strong_evidence and frac >= 0.10:
        verdict = "concordant-dependency"
    elif inert:
        verdict = "inert-in-panel"
    elif strong_evidence:
        verdict = "evidence-without-dependency"
    elif has_evidence and frac >= 0.10:
        # Evidence above EVIDENCE_FLOOR but below the strong bar, WITH a
        # dependency signal. Both halves are real, so neither
        # `dependency-without-evidence` nor `evidence-without-dependency`
        # describes the row: this branch must precede both.
        verdict = "dependency-with-thin-evidence"
    elif frac >= 0.10:
        verdict = "dependency-without-evidence"
    elif has_evidence:
        verdict = "evidence-without-dependency"
    else:
        verdict = "no-evidence-no-dependency"

    return {
        "gene": symbol,
        "depmap_inferred_essential": effect_stats.get("depmap_inferred_essential"),
        "ot_score": None if ot_score is None else round(float(ot_score), 3),
        "depmap_class": cls,
        "mean_effect": mean_eff,
        "frac_dependent": frac,
        "lineage_mean_effect": lin_mean,
        "lineage_frac_dependent": lin_frac,
        "lineage_p_bh": lin_p,
        "tractable_modalities": tractable_modalities or [],
        "verdict": verdict,
        "knockout_actionable": verdict in ("concordant-dependency",
                                          "context-restricted",
                                          "dependency-with-thin-evidence"),
    }


def fusion_notes():
    """Interpretation guide for the verdict vocabulary."""
    return {
        "concordant-dependency":
            "Human evidence and cell-line dependency agree. Strongest KO/degrader case.",
        "common-essential-no-window":
            "Nearly every line depends on it (frac_dependent >= 0.90). Real "
            "dependency, but no therapeutic window - normal tissue needs it too. "
            "Only viable with tumour-selective delivery or a genotype-restricted "
            "partner.",
        "growth-suppressive-mismatch":
            "High OT evidence but knockout FAVOURS growth across many lines "
            "(positive tail + high variance). Not a knockout target; consider "
            "synthetic-lethal partners instead.",
        "inert-in-panel":
            "Flat, low-variance effect near zero: the panel lacks the context that "
            "makes this gene matter (e.g. fusion-driven or ligand-driven oncogenes "
            "with no representative line). Absence of evidence, not evidence of "
            "absence - check lineage coverage before dismissing.",
        "context-restricted":
            "Dependency confined to a lineage or genotype. Target with a biomarker "
            "hypothesis, not broadly.",
        "evidence-without-dependency":
            "Genetic/clinical evidence without cell-autonomous dependency. May act "
            "non-cell-autonomously, be redundant, or need a model the panel lacks.",
        "dependency-with-thin-evidence":
            "Evidence above EVIDENCE_FLOOR but below the strong bar (ot_score "
            "in [0.10, 0.50)) together with a dependency signal "
            "(frac_dependent >= 0.10). Both halves are real, so the row is "
            "actionable, but on weaker evidence than `concordant-dependency` - "
            "read the ot_score alongside the verdict.",
        "dependency-without-evidence":
            "Cells need it but human evidence is thin - possible novel target or a "
            "general fitness gene with a poor therapeutic window.",
        "no-evidence-no-dependency":
            "Neither side supports this target: OT score below EVIDENCE_FLOOR (or "
            "absent) AND no cell-autonomous dependency. Nothing to pursue, and "
            "distinct from `evidence-without-dependency`, which requires real "
            "human evidence.",
        "indeterminate":
            "Required dispersion inputs (sd_effect / frac_positive / mean_effect) "
            "were missing, so the growth-suppressive and inert tests could not "
            "run. Pass the full depmap_selectivity() dict.",
    }


# --------------------------------------------------------------------------
# Run output: results/<topic>/<run>/
# --------------------------------------------------------------------------
#
# Runs land under a topic named for this skill, so the directory a run appears
# in names the analysis that produced it — the convention `depmap_genetea`
# already follows. Every run of this skill is a sibling under it.
RESULTS_TOPIC = "depmap_fusion"
TABLE_CSV = "fusion_verdicts.csv"
RUN_README = "README.md"
SUMMARY_MAX_WORDS = 130  # the opening paragraph is an orientation, not a report

TABLE_COLUMNS = ("gene", "verdict", "knockout_actionable", "ot_score",
                 "depmap_class", "mean_effect", "frac_dependent",
                 "lineage_mean_effect", "lineage_frac_dependent",
                 "lineage_p_bh", "depmap_inferred_essential",
                 "tractable_modalities")

# Restated from SKILL.md's "Honest limits" so every run states the same
# standing caveats rather than a hand-picked subset that can drift from the
# skill as it is revised.
STANDING_LIMITS = (
    "`inert-in-panel` is absence of evidence, not evidence of absence: check "
    "lineage coverage before dismissing a target; a fusion-driven or "
    "ligand-driven oncogene needs a representative line.",
    "DepMap measures proliferation in 2D culture. Immune, stromal, "
    "angiogenic and differentiation-dependent targets are invisible by "
    "construction.",
    "No verdict is validation. `concordant-dependency` and "
    "`context-restricted` are prioritisation signals, not confirmation.",
    "Check `depmap_class` before trusting `knockout_actionable`: "
    "pan-essential genes are genuinely required but have no therapeutic "
    "window without tumour-selective delivery.",
    "Thresholds (`dep_threshold`, `EVIDENCE_FLOOR`, "
    "`SUPPRESSOR_MIN_FRAC_POS`, `SUPPRESSOR_MIN_SD`) are conventions, not "
    "decision boundaries; report the underlying numbers alongside any "
    "verdict.",
)


def results_root(root=None):
    """Resolve this repo's results/ directory. Requires $SCIENCE_RESULTS_ROOT or root=.

    No cwd-relative default (a bare "results" silently resolves against
    whatever directory the kernel session happens to be running in, which is
    not reliably this repo's checkout) and no isdir check (unlike
    depmap-local's depmap_root(), a results root is written to, not read
    from, so it may not exist yet on a first run — the topic/run makedirs
    call creates it).
    """
    if root is None:
        root = os.environ.get("SCIENCE_RESULTS_ROOT")
    if root is None:
        raise FileNotFoundError(
            "No results root configured. Set $SCIENCE_RESULTS_ROOT to this "
            "repo's results/ directory, or pass root= explicitly.")
    return root


def fusion_run_dir(subject, root=None, topic=None, make=True):
    """Path for one fusion run: <root>/<topic>/<slug>/, with scripts/ beside it.

    `subject` is a gene symbol for a single-target dossier or a disease or
    indication name for a disease-level triage table (e.g. "lung
    adenocarcinoma"); it is slugged to snake_case because result directories
    are named that way and a raw subject may carry spaces, hyphens or case.

    `root` resolves via `results_root()` — $SCIENCE_RESULTS_ROOT or an
    explicit `root=` — before anything is created, so a misconfigured or
    unset root raises here rather than silently creating a `results/`
    directory wherever the session's cwd happens to be.

    `topic` defaults to RESULTS_TOPIC through an explicit None check rather
    than in the signature: the kernel.py sidecar loader rejects a non-literal
    default, and a rejected file defines none of this module's helpers.
    """
    root = results_root(root)
    topic = RESULTS_TOPIC if topic is None else topic
    slug = re.sub(r"[^a-z0-9]+", "_", str(subject).strip().lower()).strip("_")
    if not slug:
        raise ValueError("fusion_run_dir got an empty subject name; the run "
                         "directory is named after it")
    out_dir = os.path.join(root, topic, slug)
    if make:
        os.makedirs(os.path.join(out_dir, "scripts"), exist_ok=True)
    return out_dir


def fusion_write_table(path, rows, headers=None):
    """Write verdict rows (list of dicts) to CSV. Returns the path, or None if empty.

    Columns are `headers` first, then any further keys present in the rows in
    sorted order, so a field a caller attaches beyond fuse_target_row()'s own
    is written rather than dropped. None becomes an empty cell; a list value
    (`tractable_modalities`) is joined with ';' so it survives a round trip
    through CSV.
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
            out = {}
            for c in cols:
                v = row.get(c)
                if isinstance(v, (list, tuple)):
                    v = ";".join(str(x) for x in v)
                out[c] = "" if v is None else v
            writer.writerow(out)
    return path


def fusion_verdict_mix(rows):
    """Count fuse_target_row() rows by verdict. Returns {verdict: count}.

    Every row must carry a 'verdict': a row without one is not a completed
    classification, and folding it into a count silently understates how many
    targets were actually joined. Raises, naming the gene, rather than
    skipping the row.
    """
    counts = {}
    for r in rows:
        v = r.get("verdict")
        if not v:
            raise ValueError(
                f"row for gene {r.get('gene')!r} has no 'verdict'; pass the dicts "
                "returned by fuse_target_row(), not a partial row")
        counts[v] = counts.get(v, 0) + 1
    return counts


def fusion_check_words(text, cap, label):
    """Enforce a word cap on run-README prose, naming the count and the cap."""
    n = len(str(text or "").split())
    if n > cap:
        raise ValueError(
            f"{label} is {n} words; cap is {cap}. Move detail into the table "
            "rather than the prose.")
    return n


def fusion_markdown_table(rows, headers=None):
    """Render fuse_target_row() rows (or any list of dicts) as a markdown table."""
    rows = [dict(r) for r in (rows or [])]
    if not rows:
        return ""
    headers = list(headers or [])
    cols = [c for c in headers if any(c in r for r in rows)]
    cols += sorted({k for r in rows for k in r} - set(cols))

    def fmt(v):
        if v is None:
            return ""
        if isinstance(v, (list, tuple)):
            return ";".join(str(x) for x in v)
        if isinstance(v, bool):
            return str(v)
        if isinstance(v, float):
            return "{:.3g}".format(v)
        return str(v)

    lines = ["| " + " | ".join(cols) + " |",
             "| " + " | ".join("---" for _ in cols) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(c)) for c in cols) + " |")
    return "\n".join(lines)


def fusion_run_readme(subject, rows, summary, disease=None, files=(),
                      data_sources=(), limits=(), title=None):
    """Render the run README for one fusion run. Returns markdown text.

    `rows` is a list of fuse_target_row() dicts: one row for a single-target
    dossier, many for a disease-level triage table. Every row must carry a
    'verdict' (fusion_verdict_mix raises otherwise, naming the gene), because a
    README that silently omits a target's verdict misrepresents how many
    targets the run actually covered.

    Sections follow `coding-standards` (Result, Files, Data sources, Limits).
    A single-row `rows` renders as a target dossier; more than one renders as
    a triage table with the verdict mix and the joined n. `summary` is plain
    prose capped at SUMMARY_MAX_WORDS. `disease` names the indication for a
    triage table and is not required for a dossier.
    """
    if not rows:
        raise ValueError("fusion_run_readme got no rows; a README with no "
                         "verdicts is not a result")
    fusion_verdict_mix(rows)  # validates every row before any prose is built
    fusion_check_words(summary, SUMMARY_MAX_WORDS, "summary")
    if str(summary or "").strip().startswith("-"):
        raise ValueError("summary must be plain prose, not bullets — it is the "
                         "paragraph a non-specialist reads first")

    is_dossier = len(rows) == 1
    default_title = (f"{subject} target dossier — Open Targets × DepMap"
                     if is_dossier else
                     f"{subject} target triage — Open Targets × DepMap")
    parts = [f"# {title or default_title}", "", str(summary).strip(), "",
             "## Result", ""]

    if is_dossier:
        row = rows[0]
        parts += [
            "| quantity | value |", "| --- | --- |",
            f"| gene | {row.get('gene')} |",
            f"| verdict | {row.get('verdict')} |",
            f"| knockout_actionable | {row.get('knockout_actionable')} |",
            f"| ot_score | {row.get('ot_score')} |",
            f"| depmap_class | {row.get('depmap_class')} |",
            f"| mean_effect | {row.get('mean_effect')} |",
            f"| frac_dependent | {row.get('frac_dependent')} |", ""]
        if row.get("lineage_mean_effect") is not None:
            parts += [
                f"| lineage_mean_effect | {row.get('lineage_mean_effect')} |",
                f"| lineage_frac_dependent | {row.get('lineage_frac_dependent')} |",
                f"| lineage_p_bh | {row.get('lineage_p_bh')} |", ""]
    else:
        mix = fusion_verdict_mix(rows)
        n = len(rows)
        n_actionable = sum(1 for r in rows if r.get("knockout_actionable"))
        parts += [f"Joined n = {n} targets"
                 + (f" for {disease}" if disease else "") + ".", ""]
        parts += ["| verdict | n | % |", "| --- | --- | --- |"]
        for v, c in sorted(mix.items(), key=lambda kv: -kv[1]):
            parts.append(f"| {v} | {c} | {100.0 * c / n:.1f}% |")
        parts += ["", f"{n_actionable} of {n} targets ({100.0 * n_actionable / n:.1f}%) "
                      "are `knockout_actionable`.", "",
                 "### Verdicts", "", fusion_markdown_table(rows, TABLE_COLUMNS), ""]

    parts += ["## Files", ""]
    for entry in files:
        name, description = (entry["name"], entry["description"]) \
            if isinstance(entry, dict) else entry
        parts.append(f"- `{name}` — {description}")
    parts += ["", "## Data sources", ""] + [f"- {s}" for s in data_sources]
    parts += ["", "## Limits", ""]
    parts += [f"- {s}" for s in list(STANDING_LIMITS) + list(limits)]
    return "\n".join(parts).rstrip() + "\n"


def fusion_write_run(out_dir, subject, rows, summary, disease=None, files=(),
                     data_sources=(), limits=(), scripts=()):
    """Write one fusion run directory. Returns {name: path} for what was written.

    Writes `fusion_verdicts.csv` (one row per target from fuse_target_row()),
    the run README, and copies `scripts` into `scripts/`. `rows` is validated
    up front through fusion_verdict_mix, so a run with a malformed row fails
    before anything is written rather than leaving a half-written directory.
    """
    rows = [dict(r) for r in (rows or [])]
    fusion_verdict_mix(rows)

    os.makedirs(os.path.join(out_dir, "scripts"), exist_ok=True)
    written = {}

    table_path = fusion_write_table(os.path.join(out_dir, TABLE_CSV), rows, TABLE_COLUMNS)
    if table_path:
        written["table"] = table_path

    for src in scripts:
        dst = os.path.join(out_dir, "scripts", os.path.basename(src))
        shutil.copyfile(src, dst)
        written.setdefault("scripts", []).append(dst)

    listed = []
    if "table" in written:
        listed.append((TABLE_CSV, "one row per target: verdict, knockout_actionable, "
                                  "ot_score, DepMap class and the underlying effect "
                                  "and lineage statistics"))
    listed += list(files)

    readme = fusion_run_readme(subject, rows, summary, disease=disease,
                               files=listed, data_sources=data_sources, limits=limits)
    readme_path = os.path.join(out_dir, RUN_README)
    with open(readme_path, "w", encoding="utf-8") as fh:
        fh.write(readme)
    written["readme"] = readme_path
    return written


def fusion_link_into(canonical_dir, link_path):
    """Symlink `link_path` to `canonical_dir` so a caller's run directory
    exposes this skill's output without a second copy of it.

    A fusion run triggered from inside another skill's run (a triage brief,
    say) still writes to its own canonical `results/depmap_fusion/<slug>/`
    via `fusion_run_dir`/`fusion_write_run` — never to a path the caller
    picks — so a later standalone re-run of the same subject lands on the
    same directory instead of silently diverging from it. This links that
    canonical directory into the caller's run so it still reads as local
    structure.

    The link is relative, so it survives the whole `results/` tree being
    moved or copied together. A no-op if `link_path` already points at
    `canonical_dir`; replaces a stale symlink pointing elsewhere. Raises if
    something that is not a symlink already exists at `link_path` —
    silently overwriting a real file or directory there is not this
    function's call to make.
    """
    canonical_dir = os.path.abspath(canonical_dir)
    link_path = os.path.abspath(link_path)
    if os.path.islink(link_path):
        if os.path.realpath(link_path) == os.path.realpath(canonical_dir):
            return link_path
        os.remove(link_path)
    elif os.path.exists(link_path):
        raise FileExistsError(
            f"{link_path} exists and is not a symlink; refusing to overwrite "
            "it with a link to the canonical fusion run")
    os.makedirs(os.path.dirname(link_path), exist_ok=True)
    rel_target = os.path.relpath(canonical_dir, os.path.dirname(link_path))
    os.symlink(rel_target, link_path, target_is_directory=True)
    return link_path
