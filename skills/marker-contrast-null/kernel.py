"""Empirical-null controls for a stratified dependency contrast.

A contrast of one pre-chosen stratum (BRCA1-mutant vs wild-type, MSI vs MSS)
returns an uncorrected p-value for a hypothesis that was selected from a large
candidate pool. These helpers supply the two things that p-value is missing: the
distribution of the SAME contrast over every other eligible stratum, and a check
that the stratum is not globally shifted on unrelated genes.
"""
import csv
import os
import re
import shutil

import numpy as np
import pandas as pd
from scipy import stats

# Stratum-size band for the null scan. The floor keeps a standardised mean
# difference interpretable; the ceiling drops near-ubiquitous strata whose
# "reference" arm is the minority and which are not candidate biomarkers.
MARKER_MIN_N = 10
MARKER_MAX_N = 200
MIN_ARM_N = 5


def cohens_d(focal, reference):
    """Pooled-SD standardised mean difference, focal minus reference.

    Sign follows the input order, so with CRISPR gene effect (negative = more
    dependent) a NEGATIVE d means the focal stratum is more dependent.
    """
    focal = np.asarray(focal, dtype=float)
    reference = np.asarray(reference, dtype=float)
    focal = focal[~np.isnan(focal)]
    reference = reference[~np.isnan(reference)]
    if len(focal) < 2 or len(reference) < 2:
        raise ValueError(f"need >=2 per arm, got {len(focal)} and {len(reference)}")
    pooled = np.sqrt(((len(focal) - 1) * focal.var(ddof=1)
                      + (len(reference) - 1) * reference.var(ddof=1))
                     / (len(focal) + len(reference) - 2))
    if not pooled > 0:
        # Both arms constant. Dividing would yield +/-inf and sort to the top of
        # a scan as if it were the strongest effect, so fail loudly instead.
        raise ValueError("pooled SD is zero (both arms constant); "
                         "Cohen's d is undefined")
    return float((focal.mean() - reference.mean()) / pooled)


def bh_q(pvalues):
    """Benjamini-Hochberg adjusted p-values, input order preserved."""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.minimum(ranked, 1.0)
    return out


def stratum_contrast(values, flags, alternative="less"):
    """One stratum's contrast: n, means, Cohen's d, one-sided Mann-Whitney p.

    values: per-sample measurement (e.g. one gene's effect column).
    flags: boolean/0-1 Series indexed by sample id; True = focal stratum. A
        sample in `values` that `flags` has no entry for was never assessed
        for this marker and is EXCLUDED, not folded into the reference arm --
        silently defaulting an unmeasured sample to False would count "not
        genotyped" as "confirmed wild-type" and dilute the reference mean
        with untested samples (the same distinction depmap-local's
        depmap_fusion_contrast draws between unprofiled and fusion-negative).
    alternative: 'less' tests focal < reference, which for gene effect means
        the focal stratum is MORE dependent. Use 'greater' for the reverse.
    """
    values = pd.Series(values).astype(float)
    flags = pd.Series(flags)
    shared = values.index.intersection(flags.index)
    n_excluded = int(len(values) - len(shared))
    values = values.loc[shared]
    flags = flags.loc[shared].astype(bool)
    focal = values[flags].dropna()
    reference = values[~flags].dropna()
    # cohens_d only requires >=2 per arm to avoid a divide-by-zero; MIN_ARM_N=5
    # is the size below which SKILL.md says a standardised mean difference
    # stops being interpretable. marker_null_scan enforces it before ever
    # calling stratum_contrast, but a direct call (the documented workflow's
    # own step 1) bypassed it entirely and would return an underpowered d/p
    # with no signal that the floor was violated.
    if len(focal) < MIN_ARM_N or len(reference) < MIN_ARM_N:
        raise ValueError(
            f"need >={MIN_ARM_N} per arm (MIN_ARM_N), got {len(focal)} focal "
            f"and {len(reference)} reference; a contrast this small is not "
            "interpretable as a standardised mean difference")
    return {"n_focal": int(len(focal)), "n_reference": int(len(reference)),
            "n_excluded_no_flag": n_excluded,
            "mean_focal": float(focal.mean()),
            "mean_reference": float(reference.mean()),
            "d": cohens_d(focal, reference),
            "p": float(stats.mannwhitneyu(focal, reference,
                                          alternative=alternative).pvalue)}


def marker_null_scan(values, marker_matrix, min_n=None, max_n=None,
                     alternative="less"):
    """Run `stratum_contrast` for every eligible column of `marker_matrix`.

    marker_matrix: samples x candidate markers, non-zero = marker present.
        Must already be collapsed to one row per sample.
    Returns a DataFrame sorted by d (most-negative first) with a BH q across the
    markers actually tested -- the multiplicity a single pre-chosen marker
    hides. Interpret a marker's RANK, not only its own p.
    """
    if min_n is None:
        min_n = MARKER_MIN_N
    if max_n is None:
        max_n = MARKER_MAX_N
    values = pd.Series(values).astype(float)
    shared = marker_matrix.index.intersection(values.index)
    if len(shared) == 0:
        raise ValueError("marker_matrix and values share no sample ids")
    matrix = marker_matrix.loc[shared]
    vals = values.loc[shared]
    present = (matrix != 0).sum()
    eligible = present[(present >= min_n) & (present <= max_n)].index
    rows = []
    for marker in eligible:
        flags = matrix[marker] != 0
        if int(flags.sum()) < MIN_ARM_N or int((~flags).sum()) < MIN_ARM_N:
            continue
        row = stratum_contrast(vals, flags, alternative=alternative)
        row["marker"] = marker
        rows.append(row)
    if not rows:
        raise ValueError(f"no marker had {min_n}-{max_n} positive samples "
                         f"(matrix has {matrix.shape[1]} columns, "
                         f"{len(shared)} shared samples)")
    scan = pd.DataFrame(rows)
    scan["q"] = bh_q(scan.p.values)
    cols = ["marker", "n_focal", "n_reference", "mean_focal", "mean_reference",
            "d", "p", "q"]
    return scan[cols].sort_values("d").reset_index(drop=True)


def rank_in_null(scan, marker):
    """Where a named marker sits in its own scan: rank, percentile, BH q.

    Raises when the marker was not tested -- a marker outside the size band is
    absent from the scan, and silently reporting a rank would be worse than an
    error.
    """
    if marker not in set(scan.marker):
        raise KeyError(f"{marker!r} not in scan; it likely fell outside the "
                       f"stratum-size band ({len(scan)} markers tested)")
    row = scan.reset_index(drop=True)
    idx = int(row.index[row.marker == marker][0])
    hit = row.loc[idx]
    return {"marker": marker, "rank_by_d": idx + 1, "n_markers": len(row),
            "percentile": round(100.0 * (idx + 1) / len(row), 2),
            "d": float(hit.d), "p": float(hit.p), "q": float(hit.q),
            "n_markers_q_below_05": int((row.q < 0.05).sum()),
            "fraction_more_extreme": float((row.d < hit.d).mean())}


DELETION_CUT = 0.25
# Containment above which two deletion markers are treated as one event rather
# than as rival hypotheses. Calibrated on the MTAP/9p21 case: genuine
# neighbours score 0.82-1.00 and the nearest independent marker scores well
# below, so 0.7 separates them with margin on both sides. Re-derive by listing
# containment for a known co-deleted locus and cutting below its minimum.
CODELETION_MIN_OVERLAP = 0.7


def deletion_marker_matrix(cn_frame, cut=None):
    """Boolean deletion calls from a relative copy-number matrix.

    A gene lost by deletion carries no damaging mutation, so it is absent from a
    mutation matrix entirely and `marker_null_scan` over mutations cannot rank
    it. Copy-number markers need their null built from deletions instead.

    `cut` defaults to DELETION_CUT (0.25), a homozygous-loss threshold on DepMap
    relative copy number. Validate it per release against the deleted gene's own
    knockout effect: if the gene is truly absent, knocking it out does nothing
    (MTAP scores +0.125 in called-deleted lines against -0.001 elsewhere).
    """
    if cut is None:
        cut = DELETION_CUT
    numeric = cn_frame.apply(pd.to_numeric, errors="coerce")
    called = numeric < cut
    measured = numeric.notna()
    n_measured = int(measured.to_numpy().sum())
    if n_measured == 0:
        raise ValueError("no numeric copy-number values in cn_frame")
    fraction = int((called & measured).to_numpy().sum()) / n_measured
    # Both degenerate outcomes indicate the wrong input scale rather than an
    # unusual panel: log2 ratios centre on 0 and call everything deleted, while
    # absolute integer copy number never falls below a relative-scale cut.
    if fraction == 0.0 or fraction > 0.5:
        raise ValueError(
            f"cut={cut} calls {fraction:.0%} of measured values deleted, which is "
            "degenerate; check that the frame holds RELATIVE copy number "
            "(1.0 = diploid) rather than log2 ratios or absolute copy number")
    return called


def codeletion_partners(marker_matrix, marker, min_overlap=None):
    """Markers whose positive samples largely coincide with `marker`'s.

    A deletion removes a contiguous stretch of chromosome, so genes in that
    stretch carry near-identical marker columns and score near-identical
    contrasts. Those neighbours are not independent tests and not competing
    hypotheses -- they are the same event seen through adjacent genes.

    Overlap is ASYMMETRIC CONTAINMENT -- the larger of the two conditional
    fractions -- not Jaccard and not correlation. Deletions at one locus vary in
    extent between samples, so a narrow deletion is typically NESTED inside a
    broader one rather than coextensive with it, and Jaccard punishes that size
    asymmetry until real neighbours drop out. On the worked MTAP case the eleven
    9p21 neighbours score containment 0.82-1.00 but Jaccard only 0.09-0.39, so a
    Jaccard rule at any usable threshold reports them as independent rivals --
    exactly the error this function exists to prevent. Correlation fails for the
    separate reason that the agreeing-negative majority dominates it.

    Returns a DataFrame sorted by descending overlap.
    """
    if min_overlap is None:
        min_overlap = CODELETION_MIN_OVERLAP
    if marker not in marker_matrix.columns:
        raise KeyError(f"{marker!r} is not a column of marker_matrix")
    focal = marker_matrix[marker] != 0
    n_focal = int(focal.sum())
    if n_focal == 0:
        raise ValueError(f"{marker!r} is positive in no sample; nothing to compare")
    others = marker_matrix.drop(columns=[marker]) != 0
    n_other = others.sum()
    intersection = others.loc[focal].sum()
    containment = pd.concat([intersection / n_other.replace(0, np.nan),
                             intersection / n_focal], axis=1).max(axis=1).dropna()
    partners = containment[containment >= min_overlap].sort_values(ascending=False)
    return pd.DataFrame({"marker": partners.index, "overlap": partners.values,
                         "n_shared": intersection.reindex(partners.index).values,
                         "n_marker": n_other.reindex(partners.index).values})


def neighbourhood_check(scan, marker_matrix, marker, min_overlap=None):
    """Are the markers outranking `marker` the same deletion event, or rivals?

    `rank_in_null` alone is misleading for a copy-number marker: eleven markers
    beating MTAP reads as a weak hypothesis, but all eleven were 9p21 genes
    carried by the MTAP deletion itself, so the rank reflects ONE locus rather
    than eleven competing ones. `independent_above` is the number that answers
    the multiplicity question; `rank_by_d` is not.

    Returns the counts plus the outranking markers split into co-deleted
    partners and independent ones.
    """
    hit = rank_in_null(scan, marker)
    above = scan.reset_index(drop=True).iloc[: hit["rank_by_d"] - 1]
    partners = set(codeletion_partners(marker_matrix, marker,
                                       min_overlap=min_overlap)["marker"])
    codeleted = [m for m in above.marker if m in partners]
    independent = [m for m in above.marker if m not in partners]
    return {"marker": marker, "rank_by_d": hit["rank_by_d"],
            "n_markers": hit["n_markers"], "n_above": len(above),
            "codeleted_above": codeleted, "independent_above": independent,
            "n_codeleted_above": len(codeleted),
            "n_independent_above": len(independent),
            "d": hit["d"], "q": hit["q"]}


def global_shift_control(sample_means, flags, alternative="less"):
    """Is the focal stratum shifted on unrelated genes too?

    sample_means: each sample's mean measurement over many random genes.
    A significant shift here means a single-gene difference carries no
    gene-specific information -- the analogue of the proliferation confounder
    in drug-sensitivity work.
    """
    out = stratum_contrast(sample_means, flags, alternative=alternative)
    out["global_shift"] = bool(out["p"] < 0.05)
    return out


def gene_specificity_control(effect_frame, flags, alternative="less"):
    """The same stratum contrast applied to each column of `effect_frame`.

    Pass the focal gene plus reference genes: pathway neighbours expected to
    move together, and genes with no mechanistic link. If unrelated genes shift
    as much, the contrast is a property of those samples, not of the target.
    """
    rows = []
    for gene in effect_frame.columns:
        row = stratum_contrast(effect_frame[gene], flags, alternative=alternative)
        row["gene"] = gene
        rows.append(row)
    cols = ["gene", "n_focal", "n_reference", "mean_focal", "mean_reference",
            "d", "p"]
    return pd.DataFrame(rows)[cols].sort_values("d").reset_index(drop=True)


def sample_gene_means(read_matrix_fn, gene_names, kind="effect", n_genes=600,
                      seed=1):
    """Mean measurement per sample over a random gene subset.

    Feeds `global_shift_control`. read_matrix_fn(kind, columns) must return
    (DataFrame, missing) as depmap-local's `depmap_read_matrix` does. `seed` is
    fixed so the control is reproducible; vary it to confirm stability.
    """
    import random
    pool = list(gene_names)
    if len(pool) < n_genes:
        raise ValueError(f"only {len(pool)} genes available, need {n_genes}")
    rng = random.Random(seed)
    block, _ = read_matrix_fn(kind, rng.sample(pool, n_genes))
    return block.mean(axis=1)


# --------------------------------------------------------------------------
# Run output: results/<topic>/<run>/
# --------------------------------------------------------------------------
#
# This skill is invoked both standalone (a fresh biomarker hypothesis to
# check) and from inside another skill's run (the USP1 triage run's controls
# stayed in that run's own results/target_triage/usp1/scripts/, per
# coding-standards). The writer below is for the standalone case, under its
# own topic, so a caller does not have to hand-build the path each time.
RESULTS_TOPIC = "marker_contrast_null"
SCAN_CSV = "null_scan.csv"
SPECIFICITY_CSV = "gene_specificity.csv"
RUN_README = "README.md"
SUMMARY_MAX_WORDS = 130

# Restated from this file's SKILL.md so every run states the same standing
# caveats rather than a hand-picked subset that can drift from the skill as
# it is revised.
STANDING_LIMITS = (
    "The null tests whether a marker stands out among markers; it does not "
    "test whether the underlying biology is real.",
    "A small focal arm powers only a large |d|, so 'not distinguishable from "
    "selection noise' is a weaker claim than 'no effect' — report which one "
    "the data support.",
    "Damaging-mutation calls are not functional pathway loss: promoter "
    "hypermethylation, reversion and structural events are invisible to a "
    "mutation matrix, so a null result may reflect stratum misassignment "
    "rather than absent dependency.",
    "Markers whose carriers largely coincide (e.g. co-deleted neighbours) are "
    "not independent tests, which makes the BH correction conservative here.",
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


def mcn_run_dir(name, root=None, topic=None, make=True):
    """Path for one contrast run: <root>/<topic>/<slug>/, with scripts/ beside it.

    `name` identifies the contrast tested (e.g. "USP1 in BRCA1-mutant lines"),
    not the marker alone, since the same marker can be scanned against several
    measurements. Slugged to snake_case because result directories are named
    that way and a free-text name may carry spaces or punctuation.

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
    slug = re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")
    if not slug:
        raise ValueError("mcn_run_dir got an empty name; the run directory is "
                         "named after the contrast tested")
    out_dir = os.path.join(root, topic, slug)
    if make:
        os.makedirs(os.path.join(out_dir, "scripts"), exist_ok=True)
    return out_dir


def mcn_write_table(path, rows, headers=None):
    """Write rows (list of dicts, or a DataFrame) to CSV. Returns the path, or
    None if empty.

    Columns are `headers` first, then any further keys present in the rows in
    sorted order, so a field beyond the usual scan columns is written rather
    than dropped.
    """
    if hasattr(rows, "to_dict"):
        rows = rows.to_dict("records")
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


def mcn_check_words(text, cap, label):
    """Enforce a word cap on run-README prose, naming the count and the cap."""
    n = len(str(text or "").split())
    if n > cap:
        raise ValueError(
            f"{label} is {n} words; cap is {cap}. Move detail into the table "
            "rather than the prose.")
    return n


def mcn_verdict(rank, global_shift=None, q_floor=0.05):
    """Does the contrast survive the null-scan and global-shift controls?

    A marker "survives" only when its BH q in the null scan clears `q_floor`
    AND (when a global-shift control was run) the stratum is not globally
    shifted on unrelated genes — either failure alone is enough to explain the
    raw contrast without invoking the target gene. Returns a dict with the
    verdict and the reason, so `mcn_run_readme` states a plain conclusion
    rather than leaving the reader to weigh four numbers themselves.
    """
    reasons = []
    if rank["q"] >= q_floor:
        reasons.append(
            f"{rank['marker']} ranks {rank['rank_by_d']} of {rank['n_markers']} "
            f"markers ({rank['percentile']:.1f}th percentile) with BH q = "
            f"{rank['q']:.3g}, at or above the {q_floor} floor — "
            f"{rank['n_markers_q_below_05']} marker(s) in the whole scan clear "
            "q < 0.05.")
    if global_shift is not None and global_shift.get("global_shift"):
        reasons.append(
            f"the focal stratum is globally shifted on unrelated genes "
            f"(p = {global_shift['p']:.3g}), so a single-gene difference "
            "carries no target-specific information.")
    survives = not reasons
    return {"survives": survives, "reasons": reasons}


def mcn_run_readme(name, contrast, rank, global_shift=None, specificity=None,
                   neighbourhood=None, summary=None, files=(),
                   data_sources=(), limits=(), title=None):
    """Render the run README for one contrast. Returns markdown text.

    `contrast` is a stratum_contrast() dict, `rank` a rank_in_null() dict —
    both required, since a rank with no contrast to explain (or vice versa) is
    half the finding. `global_shift`, `specificity` and `neighbourhood` are
    optional, matching which of the three controls the caller actually ran.

    Sections follow `coding-standards` (Result, Files, Data sources, Limits).
    The Result section leads with mcn_verdict()'s plain survives/does-not-
    survive statement, because that is this skill's actual deliverable — a
    ranked table without a stated conclusion leaves the reader to redo the
    judgement the skill exists to make.
    """
    if summary is None:
        raise ValueError("mcn_run_readme needs a summary — the plain-prose "
                         "paragraph a non-specialist reads first")
    mcn_check_words(summary, SUMMARY_MAX_WORDS, "summary")
    if str(summary).strip().startswith("-"):
        raise ValueError("summary must be plain prose, not bullets")

    verdict = mcn_verdict(rank, global_shift=global_shift)
    marker = rank["marker"]

    parts = [f"# {title or f'{name} — marker-contrast null'}", "",
             str(summary).strip(), "", "## Result", ""]

    headline = ("survives" if verdict["survives"] else "does not survive")
    parts += [f"**{marker} {headline} the empirical-null and confounder "
             "controls.**", ""]
    parts += [
        "| quantity | value |", "| --- | --- |",
        f"| contrast | n={contrast['n_focal']} focal vs "
        f"n={contrast['n_reference']} reference |",
        f"| mean (focal / reference) | {contrast['mean_focal']:.4g} / "
        f"{contrast['mean_reference']:.4g} |",
        f"| Cohen's d | {contrast['d']:.3g} |",
        f"| uncorrected p | {contrast['p']:.3g} |",
        f"| rank in null | {rank['rank_by_d']} of {rank['n_markers']} "
        f"({rank['percentile']:.1f}th percentile) |",
        f"| BH q | {rank['q']:.3g} |",
        f"| markers clearing q < 0.05 | {rank['n_markers_q_below_05']} |", ""]
    if global_shift is not None:
        parts += [f"| global shift on unrelated genes | "
                 f"{'yes' if global_shift.get('global_shift') else 'no'} "
                 f"(p = {global_shift['p']:.3g}) |", ""]
    if not verdict["survives"]:
        parts += ["The contrast does not survive because:"] + \
                 [f"- {r}" for r in verdict["reasons"]] + [""]
    if neighbourhood is not None:
        parts += [f"Of the {neighbourhood['n_above']} marker(s) outranking "
                 f"{marker}, {neighbourhood['n_codeleted_above']} are "
                 f"co-deleted neighbours (the same event) and "
                 f"{neighbourhood['n_independent_above']} are independent "
                 "rivals.", ""]
    if specificity is not None:
        rows = specificity.to_dict("records") if hasattr(specificity, "to_dict") \
            else list(specificity)
        if rows:
            parts += ["### Gene specificity", "",
                     "| gene | d | p |", "| --- | --- | --- |"]
            for r in rows:
                parts.append(f"| {r['gene']} | {r['d']:.3g} | {r['p']:.3g} |")
            parts += [""]

    parts += ["## Files", ""]
    for entry in files:
        f_name, description = (entry["name"], entry["description"]) \
            if isinstance(entry, dict) else entry
        parts.append(f"- `{f_name}` — {description}")
    parts += ["", "## Data sources", ""] + [f"- {s}" for s in data_sources]
    parts += ["", "## Limits", ""]
    parts += [f"- {s}" for s in list(STANDING_LIMITS) + list(limits)]
    return "\n".join(parts).rstrip() + "\n"


def mcn_write_run(out_dir, name, contrast, rank, scan=None, global_shift=None,
                  specificity=None, neighbourhood=None, summary=None,
                  files=(), data_sources=(), limits=(), scripts=()):
    """Write one contrast-check run directory. Returns {name: path} written.

    Writes `null_scan.csv` (the full marker_null_scan() table, when given),
    `gene_specificity.csv` (the gene_specificity_control() table, when given),
    the run README, and copies `scripts` into `scripts/`.
    """
    os.makedirs(os.path.join(out_dir, "scripts"), exist_ok=True)
    written = {}

    scan_cols = ["marker", "n_focal", "n_reference", "mean_focal",
                "mean_reference", "d", "p", "q"]
    scan_path = mcn_write_table(os.path.join(out_dir, SCAN_CSV), scan, scan_cols)
    if scan_path:
        written["scan"] = scan_path

    spec_cols = ["gene", "n_focal", "n_reference", "mean_focal",
                "mean_reference", "d", "p"]
    spec_path = mcn_write_table(os.path.join(out_dir, SPECIFICITY_CSV),
                                specificity, spec_cols)
    if spec_path:
        written["specificity"] = spec_path

    for src in scripts:
        dst = os.path.join(out_dir, "scripts", os.path.basename(src))
        shutil.copyfile(src, dst)
        written.setdefault("scripts", []).append(dst)

    listed = []
    if "scan" in written:
        listed.append((SCAN_CSV, "the full null scan: every eligible marker's "
                                 "n, means, Cohen's d, p and BH q"))
    if "specificity" in written:
        listed.append((SPECIFICITY_CSV, "gene_specificity_control() rows: the "
                                        "same contrast on pathway neighbours "
                                        "and unrelated controls"))
    listed += list(files)

    readme = mcn_run_readme(name, contrast, rank, global_shift=global_shift,
                            specificity=specificity, neighbourhood=neighbourhood,
                            summary=summary, files=listed,
                            data_sources=data_sources, limits=limits)
    readme_path = os.path.join(out_dir, RUN_README)
    with open(readme_path, "w", encoding="utf-8") as fh:
        fh.write(readme)
    written["readme"] = readme_path
    return written


def mcn_link_into(canonical_dir, link_path):
    """Symlink `link_path` to `canonical_dir` so a caller's run directory
    exposes this skill's output without a second copy of it.

    A marker-null check run from inside another skill's run (a triage, a
    fusion, a brief) still writes to its own canonical
    `results/marker_contrast_null/<slug>/` via `mcn_run_dir`/`mcn_write_run`
    — never to a path the caller picks — so a later standalone re-run of the
    same contrast lands on the same directory instead of silently diverging
    from it. This links that canonical directory into the caller's run so it
    still reads as local structure.

    The link is relative, so it survives the whole `results/` tree being
    moved or copied together. A no-op if `link_path` already points at
    `canonical_dir`; replaces a stale symlink pointing elsewhere. Raises if
    something that is not a symlink already exists at `link_path` — silently
    overwriting a real file or directory there is not this function's call
    to make.
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
            "it with a link to the canonical marker-null run")
    os.makedirs(os.path.dirname(link_path), exist_ok=True)
    rel_target = os.path.relpath(canonical_dir, os.path.dirname(link_path))
    os.symlink(rel_target, link_path, target_is_directory=True)
    return link_path
