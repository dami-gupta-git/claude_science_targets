"""Calibrate the r threshold for depmap_codependencies.

The panel run showed that raw r is NOT comparable across queries: the negative
control OR5A1 returned a top partner at r=0.345, higher than SMARCA4's best true
complex partner (SMARCB1, r=0.303). This script establishes what does separate
signal from background:

  1. an empirical null for top1 r over random query genes, split by whether the
     query is ever a dependency in any line;
  2. the query-side signal metrics that predict an interpretable list;
  3. how much of a junk query's top list is explained by genomic proximity
     (shared local copy number), which Chronos does not fully remove.

Run: python calibrate_codependency.py   (writes codependency_calibration*.csv)
"""
import numpy as np
import pandas as pd

import kernel as cd
from validate_codependency import PANEL

N_NULL = 300
NULL_SEED = 0
PROXIMITY_BP = 5_000_000


def gene_loci(root=None):
    """Gene -> (chrom, median guide start) from the Brunello guide map.

    Chronos corrects copy number screen-wide but not perfectly, so neighbouring
    genes retain correlated residual gene effect. Coordinates let a proximity
    artefact be flagged rather than mistaken for a functional relationship.
    Uses only Chronos-used guides, so the coordinates match the scored matrix.
    """
    import os
    path = os.path.join(cd.depmap_root(root), "BrunelloGuideMap.csv")
    guides = pd.read_csv(path, usecols=["chrom", "start", "Gene", "UsedByChronos"])
    guides = guides[guides["UsedByChronos"].astype(bool)].dropna(subset=["chrom"])
    guides["symbol"] = [cd.strip_entrez(g) for g in guides["Gene"]]
    grouped = guides.groupby("symbol").agg(
        chrom=("chrom", lambda s: s.mode().iloc[0]), start=("start", "median"))
    return grouped


def build(root=None):
    prep = cd.codep_prepare(root=root)
    loci = gene_loci(root=root)
    effect = pd.read_parquet(
        cd.codep_matrix_path("effect", root=root))
    effect.columns = [cd.strip_entrez(c) for c in effect.columns]

    def signal(gene):
        v = effect[gene].to_numpy(dtype=np.float64)
        v = v[np.isfinite(v)]
        return {
            "profile_sd": float(v.std(ddof=1)),
            "min_effect": float(v.min()),
            "n_dep_lines": int((v < -0.5).sum()),
            "frac_dep_lines": float((v < -0.5).mean()),
            "n_lines": int(v.size),
        }

    def proximity_frac(query, top):
        """Fraction of top partners within PROXIMITY_BP of the query locus."""
        if query not in loci.index:
            return np.nan
        qc, qs = loci.loc[query, "chrom"], loci.loc[query, "start"]
        near = 0
        for g in top:
            if g not in loci.index:
                continue
            if loci.loc[g, "chrom"] == qc and abs(loci.loc[g, "start"] - qs) <= PROXIMITY_BP:
                near += 1
        return near / len(top)

    # ---- panel queries
    panel_rows = []
    for query, complex_name, expected in PANEL:
        full = cd.depmap_codependencies(query, n=None, root=root)
        top15 = full.head(15)
        row = {"query": query, "complex": complex_name, "kind": "panel"}
        row.update(signal(query))
        row.update({
            "top1_r": round(float(full["r"].iloc[0]), 4),
            "r_rank15": round(float(full["r"].iloc[14]), 4),
            "background_r_sd": round(float(full["r"].std(ddof=1)), 4),
            "top1_over_sd": round(float(full["r"].iloc[0] / full["r"].std(ddof=1)), 2),
            "n_expected_in_top15": int(top15["gene"].isin(expected).sum()),
            "prox_frac_top15": round(proximity_frac(query, top15["gene"].tolist()), 3),
        })
        panel_rows.append(row)

    # ---- empirical null over random queries
    rng = np.random.default_rng(NULL_SEED)
    eligible = np.asarray(cd.depmap_screened_genes(root=root,
                                                  min_lines=cd.CODEP_MIN_LINES))
    picks = rng.choice(eligible, N_NULL, replace=False)
    null_rows = []
    for query in picks:
        full = cd.depmap_codependencies(query, n=None, root=root)
        s = signal(query)
        null_rows.append({
            "query": query, "kind": "null",
            "top1_r": float(full["r"].iloc[0]),
            "r_rank15": float(full["r"].iloc[14]),
            "background_r_sd": float(full["r"].std(ddof=1)),
            "top1_over_sd": float(full["r"].iloc[0] / full["r"].std(ddof=1)),
            "prox_frac_top15": proximity_frac(query, full["gene"].head(15).tolist()),
            **s,
        })
    return pd.DataFrame(panel_rows), pd.DataFrame(null_rows)


if __name__ == "__main__":
    panel, null = build()
    panel.to_csv("codependency_calibration_panel.csv", index=False)
    null.to_csv("codependency_calibration_null.csv", index=False)
    pd.set_option("display.width", 220)

    print("PANEL")
    print(panel[["query", "profile_sd", "min_effect", "n_dep_lines", "top1_r",
                 "r_rank15", "background_r_sd", "top1_over_sd",
                 "n_expected_in_top15", "prox_frac_top15"]].to_string(index=False))

    strong = null[null["n_dep_lines"] >= 10]
    weak = null[null["n_dep_lines"] == 0]
    print("\nNULL n=%d  (>=10 dep lines: %d, zero dep lines: %d)"
          % (len(null), len(strong), len(weak)))
    for name, sub in [("all", null), ("has_dep>=10", strong), ("no_dep", weak)]:
        if not len(sub):
            continue
        q = np.percentile(sub["top1_r"], [50, 75, 90, 95, 99])
        print("  top1_r %-12s median %.3f p75 %.3f p90 %.3f p95 %.3f p99 %.3f max %.3f"
              % (name, *q, sub["top1_r"].max()))
    print("  frac of null queries with top1_r >= 0.20: %.3f"
          % (null["top1_r"] >= 0.20).mean())
    print("  frac of null queries with top1_r >= 0.30: %.3f"
          % (null["top1_r"] >= 0.30).mean())
    print("  frac of null queries with top1_r >= 0.40: %.3f"
          % (null["top1_r"] >= 0.40).mean())
    print("  proximity frac of top15, null median %.3f  panel median %.3f"
          % (null["prox_frac_top15"].median(), panel["prox_frac_top15"].median()))
    print("  null top1_over_sd: median %.2f p95 %.2f max %.2f"
          % (null["top1_over_sd"].median(),
             np.percentile(null["top1_over_sd"], 95), null["top1_over_sd"].max()))
