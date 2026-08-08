"""Validate depmap_codependencies against known complexes.

Runs a panel of query genes spanning complex types (chromatin remodeller,
cohesin, condensin, mitoribosome, paralog pair, lineage TF) plus a negative
control, and records where the last true partner sits relative to the background
r distribution. That gap is what calibrates CODEP_STRONG_R.

Run: python validate_codependency.py   (writes codependency_validation.csv)
"""
import time

import numpy as np
import pandas as pd

import kernel as cd

TOP_N = 15

# Curated expected partners. Each set is the textbook complex/pathway membership
# for the query, written down BEFORE looking at the output so recovery is a
# prediction and not a post-hoc relabelling of whatever came back.
PANEL = [
    ("SMARCA4", "SWI/SNF (BAF) remodeller", {
        "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1", "SMARCE1", "ARID1A",
        "ARID2", "DPF2", "PBRM1", "SMARCA2", "BRD9", "ACTL6A"}),
    ("ARID1B", "SWI/SNF (BAF) remodeller", {
        "SMARCB1", "SMARCC1", "SMARCC2", "SMARCD1", "SMARCE1", "ARID1A",
        "ARID2", "DPF2", "SMARCA4", "SMARCA2", "ACTL6A"}),
    ("RAD21", "cohesin", {
        "SMC1A", "SMC3", "STAG1", "STAG2", "NIPBL", "PDS5A", "PDS5B",
        "WAPL", "ESCO2", "SGO1", "MAU2"}),
    ("SMC2", "condensin", {
        "SMC4", "NCAPD2", "NCAPG", "NCAPH", "NCAPD3", "NCAPG2", "NCAPH2"}),
    ("MRPL11", "mitochondrial ribosome / translation", {
        "MRPL%d" % i for i in range(1, 60)} | {
        "MRPS%d" % i for i in range(1, 40)} | {
        "TUFM", "TSFM", "GFM1", "GFM2", "MRRF", "PTCD3", "MTIF2", "MTIF3",
        "ERAL1", "MTRF1L", "OXA1L", "MALSU1", "MTERF3", "MRPL57"}),
    ("BMS1", "90S pre-ribosome / rRNA processing", {
        "RCL1", "TSR1", "RRP9", "UTP4", "UTP6", "UTP14A", "UTP15", "UTP18",
        "UTP20", "WDR3", "WDR36", "WDR43", "WDR75", "PWP2", "HEATR1",
        "NOL10", "NOL11", "KRR1", "PNO1", "RRP7A", "RRP12", "NOB1", "LTV1",
        "TBL3", "IMP3", "IMP4", "NOP14", "NOP56", "NOP58", "FCF1", "ESF1"}),
    ("VPS4A", "ESCRT-III disassembly (paralog-buffered pair)", {
        "VPS4B", "VTA1", "CHMP2A", "CHMP4B", "CHMP6", "IST1", "VPS28",
        "VPS37A", "CHMP3", "SNF8", "VPS25", "VPS36"}),
    ("CTNNB1", "WNT / beta-catenin transcription", {
        "TCF7L2", "TCF7", "LEF1", "TCF7L1", "APC", "AXIN1", "CSNK1A1",
        "BCL9", "BCL9L", "PYGO2", "TLE3", "CREBBP", "EP300", "RNF43",
        "PORCN", "FZD5", "LGR5", "CTNNBIP1"}),
    ("SOX10", "melanocyte lineage transcription", {
        "MITF", "PAX3", "TFAP2A", "EDNRB", "DCT", "TYR", "MLANA", "PMEL",
        "IRF4", "ZEB2", "SOX9", "BRAF", "MAP2K1", "MAPK1"}),
    ("OR5A1", "olfactory receptor (negative control)", set()),
]


def _family_hits(top_genes, query):
    """Same-family partners by symbol prefix, as a looser recovery check.

    Complements the curated sets for families the sets cannot enumerate
    exhaustively (MRPL/MRPS, NCAP, CHMP), using the query's alphabetic stem.
    """
    stem = "".join(ch for ch in query if not ch.isdigit())[:4]
    return [g for g in top_genes if g.startswith(stem) and g != query]


def run_panel(root=None, top_n=TOP_N):
    """Score every panel query and return (per-partner table, per-query summary)."""
    rows, summary = [], []
    for query, complex_name, expected in PANEL:
        t0 = time.time()
        full = cd.depmap_codependencies(query, n=None, root=root)
        elapsed = time.time() - t0
        top = full.head(top_n)
        bg_sd = float(full["r"].std(ddof=1))
        bg_mean = float(full["r"].mean())

        hit_flags = [g in expected for g in top["gene"]]
        n_hit = int(sum(hit_flags))
        # r of the last expected partner inside the top slice - the boundary
        # between recovered biology and background for this query.
        hit_r = [float(r) for r, h in zip(top["r"], hit_flags) if h]
        last_hit_r = min(hit_r) if hit_r else np.nan
        first_miss_r = next(
            (float(r) for r, h in zip(top["r"], hit_flags) if not h), np.nan)

        for rank, (_, rec) in enumerate(top.iterrows(), start=1):
            rows.append({
                "query": query, "complex": complex_name, "rank": rank,
                "partner": rec["gene"], "r": round(float(rec["r"]), 4),
                "n_lines": int(rec["n_lines"]), "z": round(float(rec["z"]), 2),
                "fdr": float(rec["fdr"]),
                "expected_partner": bool(rec["gene"] in expected),
            })

        # Deepest rank at which any curated partner appears anywhere in the table.
        pos = full.reset_index(drop=True)
        found_all = pos[pos["gene"].isin(expected)]
        summary.append({
            "query": query, "complex": complex_name,
            "n_expected_curated": len(expected),
            "n_expected_in_top%d" % top_n: n_hit,
            "top1_partner": top["gene"].iloc[0],
            "top1_r": round(float(top["r"].iloc[0]), 4),
            "last_expected_r_in_top": (round(last_hit_r, 4)
                                       if np.isfinite(last_hit_r) else np.nan),
            "first_unexpected_r_in_top": (round(first_miss_r, 4)
                                          if np.isfinite(first_miss_r) else np.nan),
            "best_expected_rank": (int(found_all["rank"].min())
                                   if len(found_all) else np.nan),
            "n_expected_screened": int(len(found_all)),
            "background_r_mean": round(bg_mean, 4),
            "background_r_sd": round(bg_sd, 4),
            "top1_z": round(float(top["z"].iloc[0]), 2),
            "n_genes_ranked": int(len(full)),
            "family_prefix_hits_in_top": ";".join(_family_hits(top["gene"], query)),
            "seconds": round(elapsed, 3),
        })
    return pd.DataFrame(rows), pd.DataFrame(summary)


if __name__ == "__main__":
    cd.codep_prepare()
    partners, summary = run_panel()
    assert len(partners) == len(PANEL) * TOP_N, "panel table wrong shape"
    partners.to_csv("codependency_validation.csv", index=False)
    summary.to_csv("codependency_validation_summary.csv", index=False)
    pd.set_option("display.width", 200)
    print(summary.drop(columns=["complex", "family_prefix_hits_in_top"]).to_string(index=False))
    print()
    print(partners[partners["query"].isin(["OR5A1", "VPS4A"])].to_string(index=False))
