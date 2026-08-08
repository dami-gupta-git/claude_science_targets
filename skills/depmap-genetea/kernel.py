"""Kernel helpers for the depmap-genetea skill.

Two halves that answer one question. depmap_codependencies() ranks every
screened gene by how closely its CRISPR gene-effect profile tracks a query
gene's, which proposes functional partners without naming them. The genetea_*
helpers then run overrepresentation analysis over gene DESCRIPTION text to name
what a gene set shares, and trace each term back to the sentence that produced
it. codep_enrich() runs the two in sequence with the background set correctly.

Requires the trained GeneTEA model (1.1 GB) and a built DepMap parquet cache;
genetea_load() and codep_prepare() each name the fix when their input is absent.

The primitives in the SHARED PRIMITIVES block are duplicated from depmap-local's
kernel.py so this module imports standalone. They are identical by intent: if
depmap-local's versions change, these should follow.
"""
import contextlib
import io
import os
import re
import time
import warnings

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- SHARED PRIMITIVES
# No repo-relative or machine-specific default, matching depmap-local's own
# depmap_root() (this is a separate copy, not a shared import): $DEPMAP_ROOT or
# an explicit root= is required.
DEPMAP_ROOT = None
CACHE_DIRNAME = "_cache"

MATRIX_SPECS = {
    "effect": {"stem": "crispr_gene_effect", "csv": "CRISPRGeneEffect.csv"},
    "dependency": {"stem": "crispr_gene_dependency", "csv": "CRISPRGeneDependency.csv"},
}


def depmap_root(root=None):
    """Resolve the DepMap data directory. Requires $DEPMAP_ROOT or root=."""
    if root is None:
        root = os.environ.get("DEPMAP_ROOT", DEPMAP_ROOT)
    if root is None:
        raise FileNotFoundError(
            "No DepMap root configured. Set $DEPMAP_ROOT to the directory "
            "holding CRISPRGeneEffect.csv, Model.csv, etc., or pass root=."
        )
    if not os.path.isdir(root):
        raise FileNotFoundError("DepMap root not found: %s" % root)
    return root


def strip_entrez(name):
    """'KRAS (3845)' -> 'KRAS'. Idempotent on already-clean names."""
    return re.sub(r" \(\d+\)$", "", str(name))


def bh_adjust(pvals):
    """Benjamini-Hochberg adjusted p-values, returned in the INPUT order."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(ranked, 0, 1)
    return out


# ------------------------------------------------------------------------ TUNABLES
CODEP_MIN_LINES = 200
"""Minimum pairwise-complete cell lines for a partner to be reported.

The per-gene valid-line count in 24Q2 is bimodal: 17,931 of 18,531 genes are
screened in >=100 lines (17,787 of them in >=1000) and the remaining 600 in
~70-80, with NO gene between 100 and 200. Any floor in [100, 200] therefore
selects the identical 17,931 genes, so this threshold is not knife-edge. Genes below it are dropped rather
than reported with a wide-CI r.
"""

CODEP_STRONG_R = 0.2
"""Convenience cutoff for the shared-partner view ONLY. NOT a significance test.

There is no value of r that separates real partners from background, and this
constant must not be read as one. Measured in calibrate_codependency.py over a
300-gene random null (codependency_calibration_null.csv): median top1 r for a
RANDOM screened gene is 0.270, and 88.7% of random genes reach top1 r >= 0.20.
Meanwhile true curated partners in the panel run down to r = 0.164 (SMC2/NCAPD3)
while non-partner rows in the same top-15 slices reach r = 0.550. The two
distributions overlap almost completely.

What does discriminate is the QUERY profile - see CODEP_MIN_PROFILE_SD and
codep_query_quality(). Rank and z within a query are interpretable; the absolute
r is not comparable across queries.
"""

CODEP_MIN_PROFILE_SD = 0.20
CODEP_MIN_EFFECT = -1.0
"""Query-side gate: a codependency list is only interpretable if the QUERY is
itself a dependency somewhere.

Derived in gate_codependency.py against a labelled panel (9 complex-forming
queries + the olfactory-receptor control OR5A1) and the 300-gene random null.
The gate `profile_sd >= 0.20 and min_effect <= -1.0` admits 9/9 real queries,
0/1 control, and 18.3% of random screened genes (codependency_gate_grid.csv).

Rationale: a gene that is never essential in any line has a profile made of
screen noise plus residual copy-number signal that Chronos does not fully
remove, and correlating that against 18,531 genes still returns large r. OR5A1
(profile_sd 0.128, min_effect -0.61) returned OR5AN1 at r = 0.345 - above
SMARCA4's best true partner - and 47% of its top 15 sit within 5 Mb of its own
locus on chr11, versus a null median of 6.7%. That is co-deleted neighbouring
sequence, not shared biology.

Re-derive both numbers by rerunning calibrate_codependency.py then
gate_codependency.py. Loosening min_effect to -1.5 drops real_pass to 0.889;
tightening profile_sd to 0.25 drops it to 0.778.
"""

CODEP_PROXIMITY_BP = 5_000_000
"""Window for flagging a partner as a possible copy-number artefact.

Panel median proximity fraction in the top 15 is 0.000 and the random-null
median 0.067, against 0.467 for the OR5A1 control - so a high fraction is a
usable warning sign rather than a filter.
"""


# ------------------------------------------------------------- MATRIX PREPARATION
PREP_CACHE = {}


def codep_matrix_path(dataset, root=None):
    """Parquet cache path for `dataset`; explicit error if the cache is unbuilt.

    Whole-matrix correlation needs every column, so the column-addressable CSV
    path depmap-local falls back to is not usable here - build the cache first
    with depmap-local's build_depmap_cache().
    """
    if dataset not in MATRIX_SPECS:
        raise ValueError("unknown dataset %r; have %s" % (dataset, sorted(MATRIX_SPECS)))
    root = depmap_root(root)
    path = os.path.join(root, CACHE_DIRNAME, MATRIX_SPECS[dataset]["stem"] + ".parquet")
    if not os.path.exists(path):
        raise FileNotFoundError(
            "codependency needs the parquet cache (all 18k columns at once); "
            "%s is missing. Run depmap-local's build_depmap_cache(root=%r)."
            % (path, root))
    return path


def codep_prepare(dataset="effect", root=None, rebuild=False):
    """Load and pre-centre the gene-effect matrix once, cached in the kernel.

    Returns a dict with the column-centred matrix (NaN filled with 0), the
    observation mask, and the full-matrix column sums that let a per-query
    pairwise-complete correlation be computed by subtracting only the rows the
    query is missing.

    Pre-centring each column over its own observed entries is what keeps the
    float32 matvecs accurate: it makes the sum-of-products terms small relative
    to the sum of squares, so the covariance numerator does not lose precision
    to cancellation. Verified against scipy.stats.pearsonr in the tests.

    Memory: three float32 (1208 x 18531) planes are NOT held - only the centred
    matrix and the mask (~90 MB each), plus float64 column sums (~0.15 MB).
    """
    key = (dataset, depmap_root(root))
    if not rebuild and key in PREP_CACHE:
        return PREP_CACHE[key]

    t0 = time.time()
    path = codep_matrix_path(dataset, root=root)
    frame = pd.read_parquet(path)
    genes = np.asarray([strip_entrez(c) for c in frame.columns], dtype=object)
    models = np.asarray(frame.index, dtype=object)
    values = frame.values
    del frame
    assert values.dtype == np.float32, "expected float32 cache, got %s" % values.dtype
    assert not pd.Index(genes).has_duplicates, "duplicate gene symbols in %s" % path

    mask = ~np.isnan(values)
    centred = np.where(mask, values, np.float32(0.0))
    del values
    n_full = mask.sum(0).astype(np.int64)
    assert n_full.min() > 0, "a gene column is entirely NaN"
    col_mean = (centred.sum(0, dtype=np.float64) / n_full).astype(np.float32)
    centred -= col_mean          # shifts the NaN cells too, so re-zero them
    centred[~mask] = np.float32(0.0)

    maskf = mask.astype(np.float32)
    del mask
    prep = {
        "dataset": dataset,
        "genes": genes,
        "index": pd.Index(genes),
        "models": models,
        "centred": centred,
        "maskf": maskf,
        "n_full": n_full,
        # Pre-centring per-gene mean. Required to recover ABSOLUTE Chronos
        # effects: sum_full below is computed after centring and is ~0, so it
        # cannot be used to restore the mean.
        "col_mean": col_mean.astype(np.float64),
        "sum_full": centred.sum(0, dtype=np.float64),
        "sumsq_full": np.einsum("ij,ij->j", centred, centred, dtype=np.float64),
        "load_seconds": round(time.time() - t0, 2),
    }
    PREP_CACHE[key] = prep
    return prep


# --------------------------------------------------------------------- CORRELATION
def codep_pairwise_r(prep, profile):
    """Exact pairwise-complete Pearson r of `profile` against every column.

    NaN handling: for each partner gene the correlation uses exactly the lines
    where BOTH the query and that partner were screened. Nothing is imputed and
    no line is dropped globally. This matters because DepMap missingness is
    structural, not random - 600 genes are screened in only ~70-80 lines because
    of library version, so a complete-case intersection across all 18,531 genes
    would collapse the usable line count, and mean-imputation would shrink r
    toward 0 for exactly those sparsely screened genes.

    Returns (r, n_pairs) as float64/int64 arrays over prep["genes"].
    """
    centred, maskf = prep["centred"], prep["maskf"]
    valid = np.isfinite(profile)
    n_valid = int(valid.sum())
    if n_valid < 3:
        raise ValueError("query profile has %d observed lines; need >=3" % n_valid)

    # Query vector centred over its own observed lines, zero elsewhere, so every
    # matvec below is implicitly restricted to those lines.
    q = np.zeros(centred.shape[0], dtype=np.float32)
    q[valid] = (profile[valid] - profile[valid].mean()).astype(np.float32)

    sxy = (q @ centred).astype(np.float64)
    sx = (q @ maskf).astype(np.float64)
    sxx = ((q * q) @ maskf).astype(np.float64)

    # Partner-side sums start from the full-matrix totals and subtract only the
    # rows the query is missing (typically 0-120 of 1208), which is O(|missing|*G).
    drop = np.flatnonzero(~valid)
    n = prep["n_full"].copy()
    sy = prep["sum_full"].copy()
    syy = prep["sumsq_full"].copy()
    if drop.size:
        sub_c = centred[drop]
        sub_m = maskf[drop]
        n -= sub_m.sum(0, dtype=np.int64)
        sy -= sub_c.sum(0, dtype=np.float64)
        syy -= np.einsum("ij,ij->j", sub_c, sub_c, dtype=np.float64)

    with np.errstate(invalid="ignore", divide="ignore"):
        safe_n = np.where(n > 0, n, 1)
        cov = sxy - sx * sy / safe_n
        vx = sxx - sx * sx / safe_n
        vy = syy - sy * sy / safe_n
        r = cov / np.sqrt(vx * vy)
    # The n < 3 term is load-bearing. The vx/vy terms are defensive only: a
    # zero-variance column already gives 0/0 -> NaN, and no synthesised input
    # reaches them (verified by mutation-testing the line - removing the
    # variance terms left the suite green). They are kept to make a
    # negative-variance cancellation an explicit NaN rather than a sqrt warning.
    r[(n < 3) | (vx <= 0) | (vy <= 0)] = np.nan
    return np.clip(r, -1.0, 1.0), n


def depmap_codependencies(gene, n=25, root=None, min_lines=None,
                          dataset="effect", include_self=False):
    """Genes whose CRISPR dependency profile correlates with `gene`'s.

    Ranks every screened gene by pairwise-complete Pearson r against `gene`'s
    gene-effect profile (Chronos: 0 = no effect, -1 = median common essential).
    One standardise-and-matmul pass over the whole matrix, not 18k pairwise loops.

    Args:
        gene: query symbol, Entrez suffix optional ('KRAS' or 'KRAS (3845)').
        n: rows to return; None or <=0 returns the full ranked table.
        min_lines: drop partners with fewer than this many pairwise-complete
            lines. See CODEP_MIN_LINES for why the default is not knife-edge.
        dataset: 'effect' (Chronos, the default and the right choice for
            correlation) or 'dependency' (probabilities, bounded 0-1 and so
            compressed at both ends).
        include_self: keep the query's own r=1.0 row.

    Returns:
        DataFrame ranked by r descending with columns
        [gene, r, n_lines, z, rank, p, fdr]. `z` is r standardised against the
        empirical distribution of r over all retained genes for THIS query,
        which is the comparable-across-queries statistic; `p`/`fdr` are the
        two-sided t-test on r and its BH adjustment over the retained genes.

    Raises:
        KeyError naming the query if it is not a screened gene.
    """
    if min_lines is None:
        min_lines = CODEP_MIN_LINES   # module default; loader forbids non-literal arg defaults

    prep = codep_prepare(dataset=dataset, root=root)
    query = strip_entrez(gene)
    if query not in prep["index"]:
        raise KeyError(
            "%r is not screened in the %s matrix (%d genes). Nearest symbols: %s"
            % (query, dataset, prep["genes"].size,
               sorted(g for g in prep["genes"] if g.startswith(query[:3]))[:8]))
    col = int(prep["index"].get_loc(query))

    profile = prep["centred"][:, col].astype(np.float64)
    profile[prep["maskf"][:, col] == 0] = np.nan
    r, n_pairs = codep_pairwise_r(prep, profile)

    keep = np.isfinite(r) & (n_pairs >= int(min_lines))
    if not include_self:
        keep[col] = False
    assert keep.sum() > 0, "no partner passed min_lines=%s for %s" % (min_lines, query)

    out = pd.DataFrame({
        "gene": prep["genes"][keep],
        "r": r[keep],
        "n_lines": n_pairs[keep].astype(int),
    })
    # z against this query's own r distribution: the background is the 18k
    # screened genes, so a z of 5 means 'far outside what any gene achieves',
    # independent of how correlated the query's profile is in general.
    out["z"] = (out["r"] - out["r"].mean()) / out["r"].std(ddof=1)
    out = out.sort_values("r", ascending=False, kind="mergesort").reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)

    from scipy import stats
    with np.errstate(invalid="ignore", divide="ignore"):
        rr = out["r"].to_numpy()
        dof = out["n_lines"].to_numpy() - 2
        tstat = rr * np.sqrt(dof / np.clip(1 - rr * rr, 1e-300, None))
    out["p"] = 2 * stats.t.sf(np.abs(tstat), dof)
    out["fdr"] = bh_adjust(out["p"].to_numpy())
    out.attrs["query"] = query
    out.attrs["dataset"] = dataset
    out.attrs["r_sd_background"] = float(rr.std(ddof=1))
    return out if not n or n <= 0 else out.head(int(n)).copy()


def depmap_shared_codependencies(genes, n=25, root=None,
                                 min_lines=None, dataset="effect",
                                 min_queries=2):
    """Partners shared by several query genes - the inverse of the single view.

    Cheap because each query is one matvec pass over an already-cached matrix:
    cost is len(genes) x the single-query time, with no extra memory.

    Args:
        genes: query symbols.
        n: rows to return, ranked by mean r across the queries.
        min_queries: a partner must clear CODEP_STRONG_R for at least this many
            queries to be reported. Set to len(genes) for a strict intersection.

    Returns:
        DataFrame indexed by partner gene with one r column per query plus
        [n_strong, mean_r], sorted by mean_r descending. Query genes themselves
        are excluded from the partner rows.
    """
    if min_lines is None:
        min_lines = CODEP_MIN_LINES   # module default; loader forbids non-literal arg defaults

    genes = [strip_entrez(g) for g in genes]
    if len(genes) < 2:
        raise ValueError("need >=2 query genes, got %r" % (genes,))
    cols = {}
    for g in genes:
        full = depmap_codependencies(g, n=None, root=root, min_lines=min_lines,
                                     dataset=dataset)
        cols[g] = full.set_index("gene")["r"]
    wide = pd.DataFrame(cols).drop(index=[g for g in genes], errors="ignore")
    wide["n_strong"] = (wide[genes] >= CODEP_STRONG_R).sum(axis=1)
    wide["mean_r"] = wide[genes].mean(axis=1)
    wide = wide[wide["n_strong"] >= int(min_queries)]
    wide = wide.sort_values("mean_r", ascending=False)
    return wide if not n or n <= 0 else wide.head(int(n)).copy()


# ------------------------------------------------------------------ QUALITY GATES
def codep_query_quality(gene, root=None, dataset="effect"):
    """Is a codependency list for `gene` worth interpreting at all?

    Checks the query's own Chronos profile, because r cannot do this job: 88.7%
    of random screened genes return a top partner at r >= 0.20 (see
    CODEP_STRONG_R). Call this before trusting any list, and before running term
    enrichment on one.

    Returns:
        dict with profile_sd, min_effect, n_dep_lines, n_lines, `passes`
        (bool), and `reason` naming which criterion failed.
    """
    prep = codep_prepare(dataset=dataset, root=root)
    query = strip_entrez(gene)
    if query not in prep["index"]:
        raise KeyError("%r is not screened in the %s matrix" % (query, dataset))
    col = int(prep["index"].get_loc(query))
    observed = prep["maskf"][:, col] > 0
    raw = prep["centred"][observed, col].astype(np.float64) + prep["col_mean"][col]

    out = {
        "gene": query,
        "profile_sd": float(raw.std(ddof=1)),
        "min_effect": float(raw.min()),
        "n_dep_lines": int((raw < -0.5).sum()),
        "n_lines": int(observed.sum()),
    }
    fails = []
    if out["profile_sd"] < CODEP_MIN_PROFILE_SD:
        fails.append("profile_sd %.3f < %.2f" % (out["profile_sd"], CODEP_MIN_PROFILE_SD))
    if out["min_effect"] > CODEP_MIN_EFFECT:
        fails.append("min_effect %.3f > %.2f (never a strong dependency)"
                     % (out["min_effect"], CODEP_MIN_EFFECT))
    out["passes"] = not fails
    out["reason"] = "ok" if not fails else "; ".join(fails)
    return out


def codep_gene_loci(root=None):
    """Gene -> (chrom, median guide start) from Chronos-used Brunello guides.

    Cached in the kernel. Used to flag partners that may reflect co-deleted
    neighbouring sequence rather than shared function.
    """
    key = ("loci", depmap_root(root))
    if key in PREP_CACHE:
        return PREP_CACHE[key]
    path = os.path.join(depmap_root(root), "BrunelloGuideMap.csv")
    guides = pd.read_csv(path, usecols=["chrom", "start", "Gene", "UsedByChronos"])
    guides = guides[guides["UsedByChronos"].astype(bool)].dropna(subset=["chrom"])
    guides["symbol"] = [strip_entrez(g) for g in guides["Gene"]]
    loci = guides.groupby("symbol").agg(
        chrom=("chrom", lambda s: s.mode().iloc[0]), start=("start", "median"))
    PREP_CACHE[key] = loci
    return loci


def codep_flag_proximity(table, query, root=None, window_bp=None):
    """Add a `near_query_locus` column flagging same-chromosome close partners.

    Chronos corrects copy number screen-wide but leaves residual correlation
    between neighbouring genes, so a partner within `window_bp` of the query may
    be co-deleted sequence. Flagged, not dropped - real tandem paralogue pairs
    are also neighbours.
    """
    if window_bp is None:
        window_bp = CODEP_PROXIMITY_BP   # module default; loader forbids non-literal arg defaults

    loci = codep_gene_loci(root=root)
    query = strip_entrez(query)
    out = table.copy()
    if query not in loci.index:
        out["near_query_locus"] = False
        return out
    qc, qs = loci.loc[query, "chrom"], loci.loc[query, "start"]

    def near(g):
        if g not in loci.index:
            return False
        return bool(loci.loc[g, "chrom"] == qc
                    and abs(loci.loc[g, "start"] - qs) <= window_bp)

    out["near_query_locus"] = [near(g) for g in out["gene"]]
    return out


# ------------------------------------------------------------------ SCREENED GENES
def depmap_screened_genes(root=None, dataset="effect", min_lines=None):
    """Symbols actually screened in this CRISPR release, deduped, Entrez stripped.

    This is the correct background for any term enrichment on a codependency
    list. GeneTEA's corpus carries ~35.8k genes but the release screens ~18.5k,
    so enriching against the whole corpus credits terms for being screenable -
    'protein coding', 'ribosome', 'kinase' - rather than for being shared by the
    hit list.

    Args:
        min_lines: if given, restrict to genes with at least this many observed
            lines, matching the filter depmap_codependencies applied. Pass
            CODEP_MIN_LINES to keep foreground and background consistent.

    Returns:
        Sorted list of symbols.
    """
    prep = codep_prepare(dataset=dataset, root=root)
    genes = prep["genes"]
    if min_lines is not None:
        genes = genes[prep["n_full"] >= int(min_lines)]
    return sorted(set(genes.tolist()))


# ===========================================================================
# GENETEA ENRICHMENT
#
# GeneTEA names what a gene set has in common using free-text gene descriptions
# (RefSeq, UniProt, Alliance, CIViC, Names/Aliases, Gene Location) rather than a
# curated ontology, so an enriched "term" is a word or phrase and every hit
# traces back to the sentence that produced it. These wrappers exist because the
# upstream API has three sharp edges a caller should not have to rediscover:
#
# 1. get_enriched_terms / get_enriched_continuous return None, not an empty
#    frame, when nothing survives the thresholds. That None is absorbed exactly
#    once, here, so downstream code can always call .empty / len().
# 2. Terms belonging to a synonym group carry a literal "~ " prefix in the model
#    vocabulary ("~ orotate", not "orotate") -- 17,393 of 24,831 terms, so this
#    is the common case. get_context and validate_terms require that prefix
#    verbatim and raise otherwise, so genetea_term_context resolves a bare term
#    to its synonym group for you.
# 3. The continuous path requires more than 10,000 values AFTER ID alignment to
#    the corpus, not as supplied. The upstream check raises a bare Exception;
#    genetea_correlated_terms pre-checks and raises a ValueError naming both
#    counts.
#
# Loading the trained model costs ~2.6 s and ~1.7 GB RSS, so genetea_load is
# memoised per resolved path: repeated calls return the same object.
# ===========================================================================

# No repo-relative or machine-specific default: the trained model is a 1.11 GB
# manual download (see SKILL.md), so its path is developer-specific. $GENETEA_MODEL
# or an explicit path= is required; genetea_load raises FileNotFoundError naming
# both when neither is supplied, rather than silently trying one developer's path.
DEFAULT_MODEL_PATH = None

# Aligned-length floor enforced by GeneTEA.get_correlated_terms, which raises
# when len(x_clean) < 10_000 -- so exactly 10,000 aligned values is ACCEPTED.
# Mirrored here so the wrapper can raise a ValueError naming the aligned count
# before the upstream bare Exception fires. The comparison must match upstream's
# exactly: `< MIN_CONTINUOUS_GENES`, not `<=`, or the wrapper rejects a vector
# the underlying function would have run.
MIN_CONTINUOUS_GENES = 10_000

# Default |Correlation| floor for continuous mode. GeneTEA's own default is 0.2,
# which returns zero terms on any DepMap codependency vector: term-level
# correlations are computed against a sparse TF-IDF column, so they are an order
# of magnitude smaller than gene-gene correlations. Measured over all 24,831
# terms against the genome-wide TP53 codependency vector (18,507 genes aligned),
# the largest attainable |r| was 0.1038 and the 99.9th percentile was 0.0635 --
# see continuous_term_correlations_TP53.csv. 0.03 admits 281 terms there, which
# after graph filtering yields 10 term groups. Re-derive by rerunning
# genetea_correlated_terms with corr_thresh=None and reading the |r| quantiles.
DEFAULT_CORR_THRESH = 0.03

# Column order returned by genetea_enrich. "Matching Genes in List" is a
# space-separated string of gene symbols, as GeneTEA emits it.
DISCRETE_COLUMNS = [
    "Term",
    "Synonyms",
    "p-val",
    "FDR",
    "Effect Size",
    "n Matching Genes in List",
    "n Matching Genes Overall",
    "Matching Genes in List",
    "Term Group",
]

# Column order returned by genetea_enrich_continuous. There is no
# "Matching Genes in List" here: a continuous query has no member list, and
# "Correlation" is signed while "Effect Size" is not.
CONTINUOUS_COLUMNS = [
    "Term",
    "Synonyms",
    "Correlation",
    "p-val",
    "FDR",
    "Effect Size",
    "Directed Effect Size",
    "n Matching Genes Overall",
    "Term Group",
]

MODEL_CACHE: dict[str, object] = {}


def genetea_load(path: str | None = None, force_reload: bool = False):
    """Load the trained GeneTEA model, memoised per resolved path.

    The pickle is 1.11 GB on disk and ~1.7 GB resident once loaded, so a second
    load in the same kernel would double peak memory on a 16 GiB machine for no
    benefit. Repeat calls therefore return the identical object; pass
    `force_reload=True` only when the file on disk has changed.

    Path resolution order: explicit `path`, then $GENETEA_MODEL, then
    DEFAULT_MODEL_PATH. A missing file raises FileNotFoundError naming the path
    that was tried rather than falling back to an untrained model.
    """
    resolved = path or os.environ.get("GENETEA_MODEL") or DEFAULT_MODEL_PATH
    if resolved is None:
        raise FileNotFoundError(
            "No GeneTEA model path configured. Set $GENETEA_MODEL to a trained "
            "GeneTEA.pkl, or pass path=."
        )
    resolved = os.path.abspath(os.path.expanduser(resolved))

    if not force_reload and resolved in MODEL_CACHE:
        return MODEL_CACHE[resolved]

    if not os.path.isfile(resolved):
        raise FileNotFoundError(
            f"GeneTEA model not found at {resolved!r}. Set $GENETEA_MODEL or pass "
            "path= to point at a trained GeneTEA.pkl (train with GeneTEA.train)."
        )

    with warnings.catch_warnings():
        # The pickle was written under this exact scikit-learn build; unpickling
        # emits version chatter that is not actionable for the caller.
        warnings.simplefilter("ignore")
        from GeneTEA import load as _load

        model = _load(resolved)

    MODEL_CACHE[resolved] = model
    return model


def genetea_empty_frame(columns: list[str]) -> pd.DataFrame:
    """Typed empty frame, so callers can index columns without a None check."""
    return pd.DataFrame({c: pd.Series(dtype="object") for c in columns})


def genetea_tidy(table, columns: list[str]) -> pd.DataFrame:
    """Absorb GeneTEA's None-means-nothing-enriched and fix the column order.

    The requested columns lead, in the documented order; GeneTEA's remaining
    annotation columns (Stopword, Total Info, IDF, ...) are kept after them
    rather than dropped, since discarding them would force a caller who wants
    one back to recompute the whole table. Only the LEADING columns are a stable
    contract. Columns GeneTEA did not produce are omitted rather than filled with
    NaN, so a schema change upstream surfaces as a missing column instead of a
    column of nulls.
    """
    if table is None or len(table) == 0:
        return genetea_empty_frame(columns)
    keep = [c for c in columns if c in table.columns]
    extra = [c for c in table.columns if c not in columns]
    return table[keep + extra].reset_index(drop=True)


def genetea_enrich(
    genes,
    n: int | None = 10,
    max_fdr: float = 0.05,
    background=None,
    sort_by: str = "Effect Size",
    group_subterms: bool = True,
    min_n_list: int = 2,
    n_range=None,
    exclude_terms=None,
    tea=None,
) -> pd.DataFrame:
    """Enriched description terms for a discrete gene list, as a tidy frame.

    Hypergeometric test per term against `background` (default: all 35,818
    corpus genes), Benjamini-Hochberg corrected. Returns an EMPTY DataFrame with
    DISCRETE_COLUMNS when nothing is enriched -- GeneTEA itself returns None
    there, and that is absorbed here so callers never branch on it.

    sort_by:
      "Effect Size"  (default) ranks by summed TF-IDF information the term
                     carries within the query, breaking ties on FDR. Favours
                     specific, description-dense terms.
      "Significance" ranks by FDR first, effect size second. Favours terms with
                     the strongest statistical support, which on large queries
                     skews toward broad terms matching many genes.

    group_subterms=True (default) does not merge rows: it caps the result at `n`
    TERM GROUPS instead of `n` terms, so a group like "pyrimidine|orotate++"
    contributes all of its member rows and the frame can exceed `n` rows. Rows
    always carry a "Term Group" label. Set False to cap at `n` rows.

    `n_range` bounds how many corpus genes a term may match, as [min, max];
    None uses GeneTEA's default of [0, model.max_n] (5373 for this model), which
    is what excludes near-universal words.
    """
    tea = tea if tea is not None else genetea_load()
    genes = list(genes)
    if len(genes) == 0:
        raise ValueError("genes is empty; genetea_enrich needs at least min_n_list genes")

    table = tea.get_enriched_terms(
        entities=genes,
        max_fdr=max_fdr,
        n=n,
        sort_by=sort_by,
        group_subterms=group_subterms,
        min_n_list=min_n_list,
        n_range=[None, "default"] if n_range is None else list(n_range),
        background=background,
        exclude_terms=[] if exclude_terms is None else list(exclude_terms),
        plot=False,
    )
    return genetea_tidy(table, DISCRETE_COLUMNS)


def genetea_permutation_corr_serial(tea, tfidfs, x_clean, terms, n_permutations, random_state):
    """In-process permutation null for term correlations.

    GeneTEA's own permutation path builds a multiprocessing.Pool and pickles the
    24 MB TF-IDF matrix into every task. Under the "spawn" start method -- the
    default on macOS, and unavoidable inside a kernel whose __main__ has no
    importable file -- that never completes: a 5-permutation call was killed at
    616 s wall having burned 1.0 s of CPU. The arithmetic itself is cheap
    (16 ms per permutation measured, so ~2.7 min for the full 10,000), so it is
    done serially here where it is deterministic and interruptible.
    """
    from sklearn.feature_selection import r_regression

    observed = pd.Series(r_regression(tfidfs, x_clean), index=terms)
    rng = np.random.default_rng(random_state)
    x_arr = np.asarray(x_clean, dtype=float)
    smaller = np.zeros(len(terms), dtype=np.int64)
    larger = np.zeros(len(terms), dtype=np.int64)
    obs = observed.to_numpy()

    for _ in range(n_permutations):
        null = r_regression(tfidfs, rng.permutation(x_arr))
        smaller += null <= obs
        larger += null >= obs

    pvals = np.clip(
        2 * np.minimum((smaller + 1) / (n_permutations + 1), (larger + 1) / (n_permutations + 1)),
        0,
        1,
    )
    return observed, pd.Series(pvals, index=terms)


def genetea_correlated_terms(
    x: pd.Series,
    method: str = "f-test",
    max_fdr: float | None = 0.05,
    corr_thresh: float | None | str = "default",
    handle_duplicates: str | bool = "mean",
    n_permutations: int = 10_000,
    n_jobs: int = 1,
    random_state: int = 27,
    tea=None,
) -> pd.DataFrame:
    """Per-term correlation against a whole-genome continuous vector.

    This is the unfiltered layer under `genetea_enrich_continuous`, exposed
    separately because `get_enriched_continuous` hardcodes the f-test and gives
    no way to request `method="permutation"`.

    The vector is aligned to corpus gene symbols first (`align_series` remaps
    aliases and drops genes absent from the corpus), and the >10,000 floor
    applies to the count AFTER alignment. A short vector therefore raises
    ValueError naming both counts, in place of the upstream bare Exception whose
    message does not say which number failed.

    `x` must be a float-dtype pandas Series indexed by gene symbol; GeneTEA
    asserts both. Duplicate symbols after alias remapping are averaged by
    default ("mean") because the downstream code asserts a unique index.

    method="f-test" is a closed-form F-test over the TF-IDF matrix and the
    default: it took 0.25 s end-to-end on an 18,507-gene vector. "permutation"
    shuffles the vector `n_permutations` times and runs SERIALLY here rather than
    through GeneTEA's Pool, which hangs under the spawn start method -- see
    genetea_permutation_corr_serial. Budget ~16 ms per permutation (~2.7 min at the
    10,000 default), and note that a permutation null cannot resolve p-values
    below 1/(n_permutations+1).

    `corr_thresh` defaults to DEFAULT_CORR_THRESH (0.03), not GeneTEA's 0.2,
    which returns nothing on a real DepMap vector. The literal "default"
    sentinel stands in for that constant because the skill sidecar loader
    rejects non-literal argument defaults; pass corr_thresh=None to disable
    correlation filtering entirely, which is a distinct request.
    """
    if corr_thresh == "default":
        corr_thresh = DEFAULT_CORR_THRESH   # module default; loader forbids non-literal arg defaults

    tea = tea if tea is not None else genetea_load()
    if not isinstance(x, pd.Series):
        raise TypeError(f"x must be a pandas Series indexed by gene symbol, got {type(x)!r}")
    if not pd.api.types.is_float_dtype(x):
        raise TypeError(f"x must be float dtype (GeneTEA asserts this), got dtype {x.dtype}")

    x_clean = tea.align_series(x, verbose=False, handle_duplicates=handle_duplicates)
    if len(x_clean) < MIN_CONTINUOUS_GENES:
        raise ValueError(
            f"GeneTEA continuous mode needs > {MIN_CONTINUOUS_GENES:,} values after "
            f"aligning gene IDs to the corpus; got {len(x_clean):,} aligned from "
            f"{len(x):,} supplied. This mode ranks the whole genome against the "
            "vector, so pass an unfiltered genome-wide vector (e.g. every gene's "
            "correlation to a query gene), not a shortlist -- use genetea_enrich "
            "for a discrete gene list instead."
        )

    if method == "f-test":
        corr = tea.get_correlated_terms(
            x_clean,
            verbose=False,
            handle_duplicates=handle_duplicates,
            method="f-test",
            max_fdr=max_fdr,
            corr_thresh=corr_thresh,
            random_state=random_state,
            plot=False,
        )
    elif method == "permutation":
        corr = genetea_correlated_terms_permutation(
            tea, x_clean, n_permutations, random_state, max_fdr, corr_thresh
        )
    else:
        raise ValueError(f"method must be 'f-test' or 'permutation', got {method!r}")

    if corr is None:
        return genetea_empty_frame(CONTINUOUS_COLUMNS + ["Stopword", "IDF", "Weighted Query Info"])
    return corr


def genetea_correlated_terms_permutation(
    tea, x_clean, n_permutations, random_state, max_fdr, corr_thresh
):
    """Assemble the same frame get_correlated_terms returns, permutation p-values.

    The annotation columns downstream of the test (Stopword, Synonyms, IDF,
    Weighted Query Info, Effect Size) are reproduced here because GeneTEA adds
    them inside the method whose test branch is being replaced; the column set
    must match so `corr_terms=` accepts either method's output.
    """
    import scipy.stats

    tfidfs, terms, entities = tea.get_subset_tfidfs(entities=x_clean.index, terms=None)
    x_clean = x_clean[entities]
    observed, pvals = genetea_permutation_corr_serial(
        tea, tfidfs, x_clean, terms, n_permutations, random_state
    )

    corr = pd.DataFrame({"p-val": pvals, "Correlation": observed})
    corr["FDR"] = scipy.stats.false_discovery_control(corr["p-val"], method="bh")
    corr.index.name = "Term"

    if max_fdr is not None:
        corr = corr.loc[corr["FDR"] < max_fdr]
    if corr_thresh is not None:
        corr = corr.loc[corr["Correlation"].abs() >= corr_thresh]
    if corr.shape[0] == 0:
        return None

    corr["Stopword"] = corr.index.isin(tea.stopwords)
    corr["Synonyms"] = tea.map_synonyms(corr.index)
    corr["Weighted Query Info"] = pd.DataFrame(
        tfidfs.T.multiply(x_clean).T.mean(axis=0), columns=terms
    ).iloc[0]
    corr["Total Info"] = tea.get_total_info(corr.index, entities=entities)
    corr["IDF"] = tea.idfs
    corr["n Matching Genes Overall"] = tea.n_matching_overall
    corr["Directed Effect Size"] = corr["Correlation"].abs() * corr["Weighted Query Info"]
    corr["Effect Size"] = corr["Directed Effect Size"].abs()
    return corr.reset_index()


def genetea_enrich_continuous(
    x: pd.Series | None = None,
    corr_terms: pd.DataFrame | None = None,
    n: int | None = 10,
    max_fdr: float = 0.05,
    corr_thresh: float | str = "default",
    sort_by: str = "Effect Size",
    group_subterms: bool = True,
    n_range=None,
    handle_duplicates: str | bool = "mean",
    method: str = "f-test",
    tea=None,
) -> pd.DataFrame:
    """Enriched terms for a whole-genome ranked vector, as a tidy frame.

    Use this when the query is a ranking rather than a set -- every gene's
    codependency correlation with a query gene, a differential-dependency
    t-statistic, a genome-wide effect column. Terms are scored by how strongly
    their TF-IDF weight tracks the vector, so both tails are informative:
    "Correlation" is signed, while "Effect Size" is its magnitude times the
    term's weighted query information.

    Returns an EMPTY DataFrame with CONTINUOUS_COLUMNS when nothing survives.
    Raises ValueError via genetea_correlated_terms when fewer than
    MIN_CONTINUOUS_GENES values align to the corpus.

    `sort_by` and `group_subterms` behave as in genetea_enrich. Pass a
    precomputed `corr_terms` (from genetea_correlated_terms) instead of `x` to
    re-threshold without recomputing the correlations, or to apply
    method="permutation" results. When doing so, pass the SAME `corr_thresh`
    used to build `corr_terms`: this stage re-applies the filter, so a looser
    value here than upstream silently returns an empty frame.
    """
    if corr_thresh == "default":
        corr_thresh = DEFAULT_CORR_THRESH   # module default; loader forbids non-literal arg defaults

    tea = tea if tea is not None else genetea_load()
    if (x is None) == (corr_terms is None):
        raise ValueError("provide exactly one of x= or corr_terms=")

    if corr_terms is None:
        corr_terms = genetea_correlated_terms(
            x,
            method=method,
            max_fdr=max_fdr,
            corr_thresh=corr_thresh,
            handle_duplicates=handle_duplicates,
            tea=tea,
        )
    if len(corr_terms) == 0:
        return genetea_empty_frame(CONTINUOUS_COLUMNS)

    table = tea.get_enriched_continuous(
        corr_terms=corr_terms,
        max_fdr=max_fdr,
        corr_thresh=corr_thresh,
        n=n,
        sort_by=sort_by,
        group_subterms=group_subterms,
        n_range=[None, "default"] if n_range is None else list(n_range),
        plot=False,
    )
    return genetea_tidy(table, CONTINUOUS_COLUMNS)


def codep_enrich(
    gene,
    n_partners: int = 50,
    n: int | None = 10,
    max_fdr: float = 0.05,
    root=None,
    dataset: str = "effect",
    min_lines=None,
    corpus_background: bool = False,
    check_query: bool = True,
    tea=None,
):
    """Codependency partners of `gene`, then the terms naming what they share.

    Runs depmap_codependencies() and passes the partner symbols to
    genetea_enrich() with the background set to the SCREENED genes intersected
    with the GeneTEA corpus. That background is the whole point of this wrapper:
    enriching a codependency list against all 35,818 corpus genes credits terms
    for being screenable rather than for being shared by the partners. Measured
    over six real queries (background_comparison_summary.csv), the whole-corpus
    background returns 427 terms against 276 for the screened background -- 185
    terms drop out, and the generic vocabulary that dominates the dropped set
    ("~ DNA", "~ transcription", "~ kinase", "~ cancer") is the artefact itself.
    On one of the six (SOX10) the top term changes, so this is not cosmetic.

    The correction costs real power and is not purely a cleanup: VPS4A loses
    "~ endosome" (corpus FDR 3.3e-3) and SOX10 loses "~ waardenburg", both
    relevant. Pass corpus_background=True to see the uncorrected run for
    context, but report the corrected one as primary.

    check_query=True refuses a query that fails codep_query_quality(), because a
    gene that is never a dependency in any line yields a partner list built from
    screen noise and residual copy number -- see CODEP_MIN_PROFILE_SD. Pass
    False to override deliberately.

    Returns:
        (partners, terms) -- the depmap_codependencies() frame (top n_partners,
        already sliced) and the genetea_enrich() frame. `terms` is empty when
        nothing is enriched; it is never None.
    """
    if min_lines is None:
        min_lines = CODEP_MIN_LINES   # module default; loader forbids non-literal arg defaults

    if check_query:
        q = codep_query_quality(gene, root=root, dataset=dataset)
        if not q["passes"]:
            raise ValueError(
                "%s fails the query-side gate (%s): profile_sd=%.3f, "
                "min_effect=%.3f, n_dep_lines=%d. A gene that is not a "
                "dependency anywhere gives a codependency list built from screen "
                "noise and residual copy number. Pass check_query=False to "
                "override." % (gene, q["reason"], q["profile_sd"],
                               q["min_effect"], q["n_dep_lines"]))

    partners = depmap_codependencies(
        gene, n=n_partners, root=root, dataset=dataset, min_lines=min_lines)

    tea = tea if tea is not None else genetea_load()
    if corpus_background:
        background = None
    else:
        background = codep_enrich_background(
            root=root, dataset=dataset, min_lines=min_lines, tea=tea)

    terms = genetea_enrich(
        partners["gene"].tolist(), n=n, max_fdr=max_fdr,
        background=background, tea=tea)
    return partners, terms


def codep_enrich_background(root=None, dataset: str = "effect",
                            min_lines=None, tea=None):
    """Screened genes intersected with the GeneTEA corpus, as a sorted list.

    The intersection is required, not defensive: 28 screened symbols are absent
    from the corpus (background_gene_counts.csv). Passing those to GeneTEA
    inflates the background denominator with genes that cannot contribute to any
    term. Of 18,531 screened symbols, 18,503 are in the 35,818-gene corpus.
    """
    if min_lines is None:
        min_lines = CODEP_MIN_LINES   # module default; loader forbids non-literal arg defaults

    tea = tea if tea is not None else genetea_load()
    screened = set(depmap_screened_genes(
        root=root, dataset=dataset, min_lines=min_lines))
    corpus = set(pd.Index(tea.entities).astype(str))
    shared = sorted(screened & corpus)
    if not shared:
        raise RuntimeError(
            "no overlap between %d screened genes and the %d-gene GeneTEA "
            "corpus -- check the model and release are both human"
            % (len(screened), len(corpus)))
    return shared


def genetea_resolve_term(term: str, tea=None) -> str:
    """Map a term as a human would type it onto the model's vocabulary entry.

    Terms that belong to a synonym group are stored with a literal ``"~ "``
    prefix naming the group -- ``"~ orotate"`` covers {orotate, orotic}. The
    prefix is not decoration: `validate_terms` and `get_context` match the
    vocabulary exactly and raise RuntimeError("no valid terms") on the bare
    form. This resolves, in order: exact vocabulary hit, the ``"~ "``-prefixed
    form, and membership in a synonym group.

    Raises KeyError listing the closest vocabulary entries when nothing matches,
    rather than silently substituting a near neighbour.
    """
    tea = tea if tea is not None else genetea_load()
    vocab = set(tea.terms)

    if term in vocab:
        return term
    prefixed = f"~ {term}"
    if prefixed in vocab:
        return prefixed
    if tea.all_syns is not None and term in tea.all_syns.index:
        group = tea.all_syns.loc[term]
        group = group.iloc[0] if isinstance(group, pd.Series) else group
        if group in vocab:
            return group

    nearest = list(tea.get_most_similar_terms(term, n=5).values)
    raise KeyError(
        f"{term!r} is not in the GeneTEA vocabulary (24,831 terms). Closest "
        f"entries: {nearest}. Synonym-group terms need the '~ ' prefix, which "
        "this function adds automatically when the group exists."
    )


def genetea_split_excerpt(raw: str, sources: set[str]) -> list[tuple[str, str]]:
    """Split GeneTEA's ``"Source: text | Source: text"`` excerpt string.

    Splitting on " | " alone is unsafe because description sentences can contain
    that sequence, so a fragment is only treated as a new source when its prefix
    before ": " is a known corpus source name; otherwise it is appended to the
    fragment before it.
    """
    out: list[list[str]] = []
    for frag in raw.split(" | "):
        head, sep, tail = frag.partition(": ")
        if sep and head in sources:
            out.append([head, tail])
        elif out:
            out[-1][1] += " | " + frag
        else:
            out.append(["", frag])
    return [(s, genetea_dedupe_sentences(t)) for s, t in out]


def genetea_dedupe_sentences(text: str) -> str:
    """Collapse repeated sentences in one source's excerpt, preserving order.

    A synonym-group term matches once per synonym, so get_context returns the
    same sentence once for each member that hit it and joins the copies with
    "...". Querying "~ orotate" (members orotate, orotidine) against UMPS
    returned its UniProt paragraph four times, at 1,357 characters for 337
    characters of distinct text. Deduplicating here rather than at the call site
    keeps the Excerpt column readable and quotable; the number of matches is
    already carried by the enriched table's "n Matching Genes" columns, so
    nothing is lost by dropping the repeats.
    """
    raw = [p.strip() for p in text.split("...")]
    # A sentence ending in "." immediately before the "..." separator loses
    # that period to the boundary: "X....Y" splits to ["X", ".Y"], so the
    # period lands as a stray leading "." on the FOLLOWING fragment rather
    # than staying on the sentence that owns it. Reattach it to the fragment
    # it belongs to before anything else -- otherwise every fragment but the
    # last silently loses its trailing period, not only duplicates.
    parts = []
    for p in raw:
        if p.startswith("."):
            p = p[1:].strip()
            if parts and not parts[-1].endswith("."):
                parts[-1] += "."
        parts.append(p)

    seen, kept = set(), []
    for p in parts:
        if not p:
            continue
        # Repeats differ only by internal whitespace, so exact string equality
        # can miss them; compare on a normalised key but keep the longest
        # surviving form of each distinct sentence.
        key = " ".join(p.split()).strip(".;, ").lower()
        if not key:
            continue
        if key in seen:
            for i, q in enumerate(kept):
                if " ".join(q.split()).strip(".;, ").lower() != key:
                    continue
                if len(p) > len(q):
                    kept[i] = p
            continue
        seen.add(key)
        kept.append(p)
    return "...".join(kept) if kept else text.strip()


def genetea_term_context(
    term: str,
    genes,
    tea=None,
    source: str | None = None,
    html: bool = False,
) -> pd.DataFrame:
    """Source sentences behind a term, one row per gene per source.

    This is the provenance layer: an enriched GeneTEA term is only interpretable
    if you can read the sentence that produced it, and which database said it.
    Returns columns Gene, Term, Source, Excerpt -- Source is one of Names/Aliases,
    Gene Location, RefSeq, UniProt, Alliance, CIViC.

    `term` is accepted with or without the ``"~ "`` synonym-group prefix and
    resolved via genetea_resolve_term, because the caller reading an enriched
    table should not have to know which terms are synonym groups.

    `source` restricts the output to one database and raises on an unrecognised
    name rather than returning an empty frame, since an empty frame is also the
    legitimate "no gene matched" answer and the two should be distinguishable.
    Note that Names/Aliases excerpts are gene-symbol synonym lists, not curated
    prose, so a term supported only by that source is weaker evidence than the
    same term from RefSeq or UniProt.

    An empty frame means no supplied gene has that term in its descriptions,
    which is a normal query outcome; the upstream call raises a bare Exception
    in that case and it is converted here. Genes absent from the corpus raise
    ValueError naming them rather than being dropped silently.
    """
    tea = tea if tea is not None else genetea_load()
    resolved = genetea_resolve_term(term, tea=tea)
    genes = [genes] if isinstance(genes, str) else list(genes)

    try:
        valid, invalid, _ = tea.validate_entities(genes, verbose=False)
    except RuntimeError as exc:
        # validate_entities raises RuntimeError("no valid entities") when EVERY
        # gene is unknown, and returns them in `invalid` when only some are; both
        # cases are the same caller error, so they get the same exception.
        raise ValueError(
            f"none of the {len(genes)} supplied gene(s) are in the GeneTEA corpus "
            f"and have no descriptions to quote: {sorted(genes)[:10]}"
        ) from exc
    if invalid:
        raise ValueError(
            f"{len(invalid)} gene(s) are not in the GeneTEA corpus and have no "
            f"descriptions to quote: {sorted(invalid)[:10]}"
        )

    sources = set(tea.corpus.sentences["Source"].unique())
    # get_context hardcodes verbose=True on its own validation call; keep the
    # returned frame as the only output.
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            excerpts = tea.get_context(resolved, valid, html=html)
    except Exception:
        return genetea_empty_frame(["Gene", "Term", "Source", "Excerpt"])

    if source is not None and source not in sources:
        raise ValueError(f"source must be one of {sorted(sources)}, got {source!r}")

    rows = [
        {"Gene": gene, "Term": resolved, "Source": src, "Excerpt": text}
        for gene, raw in excerpts.items()
        for src, text in genetea_split_excerpt(raw, sources)
        if source is None or src == source
    ]
    return pd.DataFrame(rows, columns=["Gene", "Term", "Source", "Excerpt"])
