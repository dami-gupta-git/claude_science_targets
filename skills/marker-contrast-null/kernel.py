"""Empirical-null controls for a stratified dependency contrast.

A contrast of one pre-chosen stratum (BRCA1-mutant vs wild-type, MSI vs MSS)
returns an uncorrected p-value for a hypothesis that was selected from a large
candidate pool. These helpers supply the two things that p-value is missing: the
distribution of the SAME contrast over every other eligible stratum, and a check
that the stratum is not globally shifted on unrelated genes.
"""
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
