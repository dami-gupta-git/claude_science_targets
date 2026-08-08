"""Kernel helpers for local DepMap release files."""
import os
import re

# No repo-relative or machine-specific default: the DepMap release is a manual
# download (depmap.org sits behind a bot-verification wall) that lands wherever
# a given machine put it. $DEPMAP_ROOT or an explicit root= is required; a
# missing value raises rather than silently trying one developer's path.
DEPMAP_ROOT = None
DEP_THRESHOLD = -0.5
CACHE_DIRNAME = "_cache"


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


def depmap_inventory(root=None):
    """List which known DepMap files are present, with sizes in MB."""
    root = depmap_root(root)
    out = {}
    for f in sorted(os.listdir(root)):
        p = os.path.join(root, f)
        if os.path.isfile(p):
            out[f] = round(os.path.getsize(p) / 1e6, 1)
    cache = os.path.join(root, CACHE_DIRNAME)
    out["_cache_built"] = os.path.isdir(cache) and bool(os.listdir(cache))
    return out


def depmap_models(root=None, columns=None):
    """Model.csv (cell-line metadata) indexed by ModelID."""
    import pandas as pd
    root = depmap_root(root)
    df = pd.read_csv(os.path.join(root, "Model.csv"), low_memory=False)
    df = df.set_index("ModelID")
    if columns is not None:
        df = df[columns]
    return df


# ---------------------------------------------------------------------------
# Wide-matrix reads: ONE contract, two interchangeable backends.
#
# Every DepMap wide matrix (genes as columns, models as rows) is read through
# depmap_read_matrix(). It resolves the parquet cache when present and the raw
# CSV otherwise, and both backends are held to the same contract:
#
#   returns (df, missing)
#     df      - index ModelID, columns = requested symbols that exist,
#               in REQUEST order, Entrez suffix stripped
#     missing - requested symbols absent from the matrix
#   raises KeyError        if NONE of the requested columns exist
#   raises FileNotFoundError if neither backend is available
#
# Keeping the two backends behind one dispatcher is deliberate: three separate
# dual-path implementations previously drifted apart, and the divergences
# (silent empty frame vs KeyError, leaked pyarrow errors) were real bugs.
# ---------------------------------------------------------------------------

MATRIX_SPECS = {
    "effect": {"stem": "crispr_gene_effect", "csv": "CRISPRGeneEffect.csv",
               "collapse_to_model": False},
    "dependency": {"stem": "crispr_gene_dependency", "csv": "CRISPRGeneDependency.csv",
                   "collapse_to_model": False},
    "hotspot": {"stem": "hotspot_mutations", "csv": "OmicsSomaticMutationsMatrixHotspot.csv",
                "collapse_to_model": True},
    "damaging": {"stem": None, "csv": "OmicsSomaticMutationsMatrixDamaging.csv",
                 "collapse_to_model": True},
}


def strip_entrez(name):
    """'KRAS (3845)' -> 'KRAS'. Idempotent on already-clean names."""
    return re.sub(r" \(\d+\)$", "", str(name))


def depmap_matrix_backend(kind, root=None):
    """Which backend serves `kind`: ('parquet', path) | ('csv', path)."""
    if kind not in MATRIX_SPECS:
        raise ValueError("unknown matrix %r; have %s" % (kind, sorted(MATRIX_SPECS)))
    spec = MATRIX_SPECS[kind]
    root = depmap_root(root)
    if spec["stem"]:
        pq_path = os.path.join(root, CACHE_DIRNAME, spec["stem"] + ".parquet")
        if os.path.exists(pq_path):
            return "parquet", pq_path
    csv_path = os.path.join(root, spec["csv"])
    if os.path.exists(csv_path):
        return "csv", csv_path
    raise FileNotFoundError(
        "no source for %r: neither %s cache nor %s in %s"
        % (kind, spec["stem"] or "(no cache)", spec["csv"], root))


def depmap_read_matrix(kind, columns, root=None):
    """Read named columns from a DepMap wide matrix. See contract above."""
    import pandas as pd
    if isinstance(columns, str):
        columns = [columns]
    columns = list(dict.fromkeys(columns))          # dedupe, keep order
    backend, path = depmap_matrix_backend(kind, root=root)
    collapse = MATRIX_SPECS[kind]["collapse_to_model"]

    if backend == "parquet":
        import pyarrow.parquet as pqmod
        have = {strip_entrez(c): c for c in pqmod.ParquetFile(path).schema.names}
        keep = [c for c in columns if c in have]
        missing = [c for c in columns if c not in have]
        if not keep:
            raise KeyError("none of %r found in %s (%s)" % (columns, kind, backend))
        df = pd.read_parquet(path, columns=[have[c] for c in keep])
        df.columns = keep
    else:
        import csv as csvmod
        with open(path) as fh:
            hdr = next(csvmod.reader(fh))
        pos, idcol = {}, 0
        for i, c in enumerate(hdr):
            clean = strip_entrez(c)
            if clean == "ModelID":
                idcol = i
            elif i or clean:
                pos.setdefault(clean, i)
        keep = [c for c in columns if c in pos]
        missing = [c for c in columns if c not in pos]
        if not keep:
            raise KeyError("none of %r found in %s (%s)" % (columns, kind, backend))
        usecols = [idcol] + [pos[c] for c in keep]
        if collapse and "IsDefaultEntryForModel" in pos:
            usecols.append(pos["IsDefaultEntryForModel"])
        df = pd.read_csv(path, usecols=sorted(set(usecols)))
        df.columns = [strip_entrez(c) for c in df.columns]
        if collapse and "IsDefaultEntryForModel" in df.columns:
            df = df[df.IsDefaultEntryForModel == "Yes"].drop(columns=["IsDefaultEntryForModel"])
        if "ModelID" in df.columns:
            df = df.set_index("ModelID")
        else:
            df = df.set_index(df.columns[0])
        df = df[keep]

    df.index.name = "ModelID"
    return df, missing


def depmap_read_cached(stem, genes, root=None):
    """Deprecated shim. Column read from the parquet cache; None when absent.

    Kept so existing callers keep working; new code should use
    depmap_read_matrix(kind, columns).
    """
    import pandas as pd
    kind = {v["stem"]: k for k, v in MATRIX_SPECS.items() if v["stem"]}.get(stem)
    if kind is None:
        return None
    root = depmap_root(root)
    if not os.path.exists(os.path.join(root, CACHE_DIRNAME, stem + ".parquet")):
        return None
    try:
        return depmap_read_matrix(kind, genes, root=root)
    except KeyError:
        return pd.DataFrame(), list(genes)


def depmap_gene_effect(genes, root=None, dataset="effect"):
    """Gene-effect (Chronos) or dependency-probability matrix for `genes`.

    dataset: 'effect' (Chronos score) or 'dependency' (probability 0-1) ONLY.
    'hotspot' and 'damaging' are mutation-count matrices, not gene-effect data,
    and are not selectable here -- read them via depmap_read_matrix(kind, ...)
    directly, as depmap_mutation_contrast does. Any other value raises rather
    than silently falling back to 'dependency', which a mistyped or
    out-of-scope dataset name would otherwise do unnoticed.
    Cell lines x genes. Missing symbols are listed in df.attrs['missing_genes'];
    a request where NO symbol exists raises KeyError on either backend.
    """
    if dataset not in ("effect", "dependency"):
        raise ValueError(
            "dataset must be 'effect' or 'dependency', got %r; 'hotspot' and "
            "'damaging' are mutation matrices, not gene-effect data -- read "
            "them via depmap_read_matrix(kind, columns) instead" % (dataset,))
    df, missing = depmap_read_matrix(dataset, genes, root=root)
    if missing:
        df.attrs["missing_genes"] = missing
    return df


def depmap_selectivity(gene, root=None):
    """Classify a gene's dependency profile: common-essential / selective / non-essential."""
    import numpy as np
    eff = depmap_gene_effect(gene, root=root)
    if gene not in eff.columns:
        raise KeyError("gene %s not in CRISPR matrix" % gene)
    v = eff[gene].dropna()
    frac = float((v < DEP_THRESHOLD).mean())
    if frac >= 0.90:
        cls = "common essential"
    elif frac >= 0.02:
        cls = "selective"
    else:
        cls = "non-essential"
    return {
        "gene": gene, "n_lines": int(v.size), "mean_effect": round(float(v.mean()), 3),
        "median_effect": round(float(v.median()), 3), "min_effect": round(float(v.min()), 3),
        "pct05_effect": round(float(np.percentile(v, 5)), 3),
        "frac_dependent": round(frac, 4), "n_dependent": int((v < DEP_THRESHOLD).sum()),
        "sd_effect": round(float(v.std()), 3),
        "frac_positive": round(float((v > 0.3).mean()), 4),
        "p95_effect": round(float(np.percentile(v, 95)), 3),
        "classification": cls,
    }


def depmap_lineage_enrichment(gene, root=None, min_lines=15):
    """Per-lineage mean gene effect, with a Mann-Whitney test vs all other lines."""
    import pandas as pd
    from scipy import stats
    effdf = depmap_gene_effect(gene, root=root)
    if gene not in effdf.columns:
        raise KeyError("gene %s not in CRISPR matrix" % gene)
    eff = effdf[gene]
    lin = depmap_models(root=root, columns=["OncotreeLineage"])["OncotreeLineage"]
    j = pd.concat([eff.rename("effect"), lin], axis=1, join="inner").dropna()
    rows = []
    for name, grp in j.groupby("OncotreeLineage"):
        if len(grp) < min_lines:
            continue
        rest = j[j.OncotreeLineage != name]["effect"]
        # A lineage covering every scored line leaves nothing to compare
        # against: mannwhitneyu on an empty sample silently returns NaN rather
        # than raising, so this must be excluded explicitly, not just skipped
        # by a size floor that doesn't apply to the comparator side.
        if len(rest) < 3:
            continue
        u = stats.mannwhitneyu(grp["effect"], rest, alternative="less")
        rows.append({
            "lineage": name, "n": len(grp),
            "mean_effect": round(float(grp["effect"].mean()), 3),
            "mean_effect_other": round(float(rest.mean()), 3),
            "frac_dependent": round(float((grp["effect"] < DEP_THRESHOLD).mean()), 3),
            "p_enriched": float(u.pvalue),
        })
    if not rows:
        # No lineage cleared both the min_lines floor and the empty-comparator
        # guard above. An empty DataFrame has no p_enriched column to sort by,
        # so sort_values would raise KeyError.
        return pd.DataFrame()
    out = pd.DataFrame(rows).sort_values("p_enriched").reset_index(drop=True)
    out["p_bh"] = bh_adjust(out["p_enriched"].to_numpy())
    return out


def bh_adjust(pvals):
    """Benjamini-Hochberg adjusted p-values. Input order is irrelevant -

    the function sorts internally and returns values aligned to the INPUT order.
    Verified against statsmodels multipletests(method="fdr_bh") to 1e-16.
    """
    import numpy as np
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    # enforce monotonicity from the largest p downward
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(ranked, 0, 1)
    return out


def depmap_mutation_contrast(dep_gene, marker_gene, root=None, kind="hotspot"):
    """Contrast dependency on `dep_gene` between marker_gene-mutant and WT lines."""
    import pandas as pd
    from scipy import stats
    try:
        mut_df, _ = depmap_read_matrix(kind, marker_gene, root=root)
    except KeyError:
        raise KeyError("marker gene %s not in the %s mutation matrix (only "
                       "recurrently mutated genes are present)" % (marker_gene, kind))
    mut = mut_df[marker_gene]
    eff = depmap_gene_effect(dep_gene, root=root)[dep_gene]
    j = pd.concat([eff.rename("effect"), mut.rename("mut")], axis=1, join="inner").dropna()
    a = j[j["mut"] > 0]["effect"]
    b = j[j["mut"] == 0]["effect"]
    if len(a) < 3 or len(b) < 3:
        short = "mutant" if len(a) < 3 else "WT"
        return {"dep_gene": dep_gene, "marker_gene": marker_gene,
                "n_mutant": int(len(a)), "n_wt": int(len(b)),
                "note": "too few %s lines" % short}
    u = stats.mannwhitneyu(a, b, alternative="less")
    pooled = ((a.std() ** 2 * (len(a) - 1) + b.std() ** 2 * (len(b) - 1)) / (len(a) + len(b) - 2)) ** 0.5
    return {
        "dep_gene": dep_gene, "marker_gene": marker_gene, "marker_kind": kind,
        "n_mutant": int(len(a)), "n_wt": int(len(b)),
        "mean_effect_mutant": round(float(a.mean()), 3),
        "mean_effect_wt": round(float(b.mean()), 3),
        "delta": round(float(a.mean() - b.mean()), 3),
        "cohens_d": round(float((a.mean() - b.mean()) / pooled), 3) if pooled else None,
        "p_mutant_more_dependent": floor_pvalue(u.pvalue),
    }


def floor_pvalue(p):
    """Clamp underflowed p-values to the float64 minimum so they stay reportable."""
    import sys
    p = float(p)
    return p if p > 0 else sys.float_info.min


def depmap_prism_releases(root=None):
    """PRISM releases present on disk, oldest first (e.g. ['23Q2', '24Q2'])."""
    root = depmap_root(root)
    rels = set()
    for f in os.listdir(root):
        m = re.match(r"Repurposing_Public_(\d+Q\d+)_Extended_Primary_Data_Matrix\.csv", f)
        if m:
            rels.add(m.group(1))
    return sorted(rels)


def depmap_prism_compounds(target, root=None, release=None, verbose=True):
    """PRISM Repurposing compounds annotated against `target`, with potency summary.

    release: '24Q2', '23Q2', ... Defaults to the newest present on disk.

    Join notes (verified on 24Q2): the compound list and the data matrix both key
    on the `IDs` column including its `BRD:` prefix - join directly, no string
    surgery. Two known one-to-many quirks are handled here:
      * a compound screened in more than one screen (e.g. AZ-628 in REP.1M and
        REP.300) has multiple list rows but ONE matrix row - the screens are
        collapsed into a `screens` field rather than emitting duplicate results;
      * a small number of matrix rows carry no compound-list annotation and are
        therefore unreachable by target (1 row in 24Q2).
    """
    import pandas as pd
    root = depmap_root(root)
    rels = depmap_prism_releases(root)
    if not rels:
        raise FileNotFoundError("no PRISM Repurposing files under %s" % root)
    if release is None:
        release = rels[-1]
    elif release not in rels:
        raise ValueError("release %s not on disk; have %s" % (release, rels))

    base = os.path.join(root, "Repurposing_Public_%s_" % release)
    cl = pd.read_csv(base + "Extended_Primary_Compound_List.csv")
    hit = cl[cl["repurposing_target"].fillna("").str.split(",").apply(
        lambda xs: target in [x.strip() for x in xs])]
    if not len(hit):
        return pd.DataFrame()

    # collapse multi-screen rows: one record per compound id
    agg = {"Drug.Name": "first", "MOA": "first", "dose": "first",
           "screen": lambda s: ",".join(sorted(set(s)))}
    hit = hit.groupby("IDs", as_index=False).agg(agg)

    mat = pd.read_csv(base + "Extended_Primary_Data_Matrix.csv", index_col=0)
    rows = []
    for _, r in hit.iterrows():
        if r["IDs"] not in mat.index:
            continue
        lfc = mat.loc[r["IDs"]].dropna()
        rows.append({
            "drug": r["Drug.Name"], "moa": r["MOA"], "dose_uM": r["dose"],
            "screens": r["screen"], "release": release,
            "n_lines": int(lfc.size),
            "median_lfc": round(float(lfc.median()), 3) if lfc.size else None,
            "min_lfc": round(float(lfc.min()), 3) if lfc.size else None,
            "frac_killed_lfc_lt_1": round(float((lfc < -1).mean()), 3) if lfc.size else None,
        })
    if not rows:
        # Every list-annotated compound id was absent from the data matrix
        # (the "small number of matrix rows carry no compound-list
        # annotation" quirk above, in reverse). An empty DataFrame has no
        # median_lfc column to sort by, so sort_values would raise KeyError.
        return pd.DataFrame()
    out = pd.DataFrame(rows).sort_values("median_lfc", na_position="last").reset_index(drop=True)
    if verbose and len(out):
        out.attrs["release"] = release
    return out


def depmap_common_essentials(root=None):
    """DepMap's own inferred common-essential gene set (symbols)."""
    import pandas as pd
    root = depmap_root(root)
    path = os.path.join(root, "CRISPRInferredCommonEssentials.csv")
    if not os.path.exists(path):
        raise FileNotFoundError("CRISPRInferredCommonEssentials.csv not in %s" % root)
    col = pd.read_csv(path).iloc[:, 0]
    return set(col.astype(str).str.replace(r" \(\d+\)$", "", regex=True))


def depmap_fusion_models(gene, root=None, canonical_only=False):
    """ModelIDs carrying a fusion involving `gene`.

    canonical_only keeps partners that are named genes, dropping readthrough and
    same-gene (`GENE--GENE`) calls plus unnamed/AC-prefixed loci.
    """
    import pandas as pd
    root = depmap_root(root)
    path = os.path.join(root, "OmicsFusionFiltered.csv")
    if not os.path.exists(path):
        raise FileNotFoundError("OmicsFusionFiltered.csv not in %s" % root)
    fu = pd.read_csv(path, low_memory=False,
                     usecols=["ModelID", "IsDefaultEntryForModel",
                              "CanonicalFusionName", "Gene1", "Gene2"])
    fu = fu[fu.IsDefaultEntryForModel == "Yes"]
    hit = fu[fu.Gene1.str.startswith(gene + " ", na=False)
             | fu.Gene2.str.startswith(gene + " ", na=False)]
    if canonical_only:
        parts = hit.CanonicalFusionName.str.split("--")
        keep = parts.apply(lambda p: isinstance(p, list) and len(p) == 2 and p[0] != p[1]
                           and not any(x.startswith(("AC0", "AC1", "RP11", "CTD-")) for x in p))
        hit = hit[keep]
    return set(hit.ModelID)


def depmap_fusion_contrast(gene, root=None, canonical_only=False, min_n=3):
    """Is `gene` a dependency specifically in its fusion-positive lines?

    Distinguishes 'no fusion-positive line exists here' (untestable) from
    'fusion-positive lines exist and are not dependent' (a real negative).
    """
    from scipy import stats
    pos = depmap_fusion_models(gene, root=root, canonical_only=canonical_only)
    profiled = depmap_fusion_profiled_models(root=root)
    eff = depmap_gene_effect(gene, root=root)
    if gene not in eff.columns:
        return {"gene": gene, "status": "gene absent from CRISPR matrix"}
    v = eff[gene].dropna()
    a = v[v.index.isin(pos)]
    # A line with NO fusion call is only a true negative if it was fusion-profiled
    # at all; unprofiled lines are UNKNOWN status and must not enter the contrast.
    b = v[(~v.index.isin(pos)) & (v.index.isin(profiled))]
    n_unprofiled = int(v[~v.index.isin(profiled)].size)
    out = {"gene": gene, "n_fusion_models": len(pos), "n_fusion_with_crispr": int(a.size),
           "n_fusion_negative": int(b.size), "n_unprofiled_excluded": n_unprofiled,
           "canonical_only": canonical_only}
    if a.size < min_n:
        out["status"] = "untestable"
        out["note"] = ("only %d fusion-positive line(s) have CRISPR data - absence of "
                       "evidence, not evidence of absence" % a.size)
        return out
    if b.size < min_n:
        # mannwhitneyu on a comparator this small returns NaN, and NaN > 0 is
        # False, so floor_pvalue(nan) would floor it to sys.float_info.min - a
        # fabricated near-zero p-value that reads as maximally significant.
        # Same failure mode depmap_lineage_enrichment already guards against;
        # this side of the contrast needs the identical guard.
        out["status"] = "untestable"
        out["note"] = ("only %d fusion-negative (profiled) line(s) available for "
                       "comparison - too few to test against" % b.size)
        return out
    p = floor_pvalue(stats.mannwhitneyu(a, b, alternative="less").pvalue)
    out.update({
        "mean_effect_fusion_pos": round(float(a.mean()), 3),
        "mean_effect_fusion_neg": round(float(b.mean()), 3),
        "frac_dependent_fusion_pos": round(float((a < DEP_THRESHOLD).mean()), 3),
        "p_fusion_more_dependent": p,
        "status": "dependent-in-fusion-positive" if (p < 0.05 and a.mean() < DEP_THRESHOLD)
                  else "not-dependent-despite-fusion",
    })
    return out


def depmap_msi_contrast(gene, root=None, msi_threshold=20.0):
    """Contrast dependency on `gene` between MSI-high and MSS lines.

    MSI status comes from OmicsGlobalSignatures.csv (`MSIScore`), not from a
    mutation matrix - `depmap_mutation_contrast()` cannot be used for it.
    Reproduces the WRN synthetic-lethality result (Chan 2019).
    """
    import pandas as pd
    from scipy import stats
    root = depmap_root(root)
    path = os.path.join(root, "OmicsGlobalSignatures.csv")
    if not os.path.exists(path):
        raise FileNotFoundError("OmicsGlobalSignatures.csv not in %s" % root)
    sig = pd.read_csv(path)
    if "IsDefaultEntryForModel" in sig.columns:
        sig = sig[sig.IsDefaultEntryForModel == "Yes"]
    msi = sig.set_index("ModelID")["MSIScore"]
    eff = depmap_gene_effect(gene, root=root)
    if gene not in eff.columns:
        raise KeyError("gene %s not in CRISPR matrix" % gene)
    j = pd.concat([eff[gene].rename("effect"), msi], axis=1, join="inner").dropna()
    hi = j[j.MSIScore >= msi_threshold]["effect"]
    lo = j[j.MSIScore < msi_threshold]["effect"]
    if min(hi.size, lo.size) < 3:
        return {"gene": gene, "status": "too few lines", "n_msi_high": int(hi.size)}
    u = stats.mannwhitneyu(hi, lo, alternative="less")
    return {
        "gene": gene, "msi_threshold": msi_threshold,
        "n_msi_high": int(hi.size), "n_mss": int(lo.size),
        "mean_effect_msi_high": round(float(hi.mean()), 3),
        "mean_effect_mss": round(float(lo.mean()), 3),
        "frac_dependent_msi_high": round(float((hi < DEP_THRESHOLD).mean()), 3),
        "frac_dependent_mss": round(float((lo < DEP_THRESHOLD).mean()), 3),
        "p_msi_more_dependent": floor_pvalue(u.pvalue),
    }


def build_depmap_cache(root=None, verbose=True):
    """Write column-addressable parquet copies of the wide matrices to _cache/.

    One-off, ~15 s. Cuts a single-gene read from ~1 s to ~0.2 s. Needs pyarrow.
    """
    import json
    import time
    import pandas as pd
    root = depmap_root(root)
    cache = os.path.join(root, CACHE_DIRNAME)
    os.makedirs(cache, exist_ok=True)
    built = {}
    for fn, stem in [("CRISPRGeneEffect.csv", "crispr_gene_effect"),
                     ("CRISPRGeneDependency.csv", "crispr_gene_dependency")]:
        src = os.path.join(root, fn)
        if not os.path.exists(src):
            continue
        t0 = time.time()
        df = pd.read_csv(src, index_col=0)
        df.index.name = "ModelID"
        df.columns = [re.sub(r" \(\d+\)$", "", c) for c in df.columns]
        df.astype("float32").to_parquet(os.path.join(cache, stem + ".parquet"),
                                        compression="zstd")
        built[stem] = {"shape": list(df.shape), "seconds": round(time.time() - t0, 1)}
        if verbose:
            print(stem, df.shape, "%.0fs" % (time.time() - t0), flush=True)
    hs_src = os.path.join(root, "OmicsSomaticMutationsMatrixHotspot.csv")
    if os.path.exists(hs_src):
        hs = pd.read_csv(hs_src, index_col=0)
        hs = hs[hs.IsDefaultEntryForModel == "Yes"].set_index("ModelID")
        hs = hs[[c for c in hs.columns if re.search(r" \(\d+\)$", c)]]
        hs.columns = [re.sub(r" \(\d+\)$", "", c) for c in hs.columns]
        hs.astype("float32").to_parquet(os.path.join(cache, "hotspot_mutations.parquet"),
                                        compression="zstd")
        built["hotspot_mutations"] = {"shape": list(hs.shape)}
    json.dump({"built": time.strftime("%Y-%m-%d %H:%M"), "format": "parquet",
               "datasets": built}, open(os.path.join(cache, "MANIFEST.json"), "w"), indent=1)
    return built


def depmap_fusion_profiled_models(root=None):
    """ModelIDs that underwent fusion calling at all (default entries).

    A model absent from this set has UNKNOWN fusion status, not negative status -
    excluding it from the negative group is what keeps `depmap_fusion_contrast()`
    from treating missing data as evidence of absence.
    """
    import pandas as pd
    root = depmap_root(root)
    path = os.path.join(root, "OmicsFusionFiltered.csv")
    if not os.path.exists(path):
        raise FileNotFoundError("OmicsFusionFiltered.csv not in %s" % root)
    fu = pd.read_csv(path, low_memory=False,
                     usecols=["ModelID", "IsDefaultEntryForModel"])
    return set(fu[fu.IsDefaultEntryForModel == "Yes"].ModelID)


def depmap_stratified_lineage_contrast(gene, stratum, root=None, min_per_group=3,
                                       threshold=None):
    """Test `gene` dependency across a binary stratum WITHIN every lineage.

    Returns EVERY lineage that meets `min_per_group` in both arms - including
    null results - with BH-adjusted p-values across the lineages tested. This
    exists because hand-picking the lineages that worked is a reporting error
    that is easy to make and hard to see: report the whole table.

    stratum: 'MSI' (MSIScore >= threshold, default 20) or the name of a column in
        OmicsGlobalSignatures.csv, or a boolean Series indexed by ModelID.
    """
    import pandas as pd
    from scipy import stats
    root = depmap_root(root)
    eff = depmap_gene_effect(gene, root=root)
    if gene not in eff.columns:
        raise KeyError("gene %s not in CRISPR matrix" % gene)
    lin = depmap_models(root=root, columns=["OncotreeLineage"])["OncotreeLineage"]

    if isinstance(stratum, pd.Series):
        flag = stratum.astype(bool).rename("stratum")
        stratum_name = stratum.name or "custom"
    else:
        sig_path = os.path.join(root, "OmicsGlobalSignatures.csv")
        if not os.path.exists(sig_path):
            raise FileNotFoundError("OmicsGlobalSignatures.csv not in %s" % root)
        sig = pd.read_csv(sig_path)
        if "IsDefaultEntryForModel" in sig.columns:
            sig = sig[sig.IsDefaultEntryForModel == "Yes"]
        sig = sig.set_index("ModelID")
        col = "MSIScore" if stratum.upper() == "MSI" else stratum
        if col not in sig.columns:
            raise KeyError("%s not in OmicsGlobalSignatures.csv" % col)
        thr = 20.0 if threshold is None else float(threshold)
        flag = (sig[col] >= thr).rename("stratum")
        stratum_name = "%s>=%g" % (col, thr)

    j = pd.concat([eff[gene].rename("effect"), lin, flag], axis=1,
                  join="inner").dropna(subset=["effect", "stratum"])
    rows = []
    for name, g in j.groupby("OncotreeLineage"):
        a = g[g["stratum"]]["effect"]
        b = g[~g["stratum"]]["effect"]
        if len(a) < min_per_group or len(b) < min_per_group:
            continue
        rows.append({
            "lineage": name, "n_positive": len(a), "n_negative": len(b),
            "mean_positive": round(float(a.mean()), 3),
            "mean_negative": round(float(b.mean()), 3),
            "delta": round(float(a.mean() - b.mean()), 3),
            "frac_dependent_positive": round(float((a < DEP_THRESHOLD).mean()), 3),
            "frac_dependent_negative": round(float((b < DEP_THRESHOLD).mean()), 3),
            "p_raw": floor_pvalue(stats.mannwhitneyu(a, b, alternative="less").pvalue),
        })
    out = pd.DataFrame(rows)
    if len(out):
        out = out.sort_values("p_raw").reset_index(drop=True)
        out["p_bh"] = bh_adjust(out["p_raw"].to_numpy())
        out["significant_bh_0.05"] = out["p_bh"] < 0.05
    out.attrs["stratum"] = stratum_name
    out.attrs["gene"] = gene
    out.attrs["n_lineages_tested"] = int(len(out))
    return out
