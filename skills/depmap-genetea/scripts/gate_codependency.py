"""Derive the query-side gate that raw r cannot provide.

calibrate_codependency.py measured the empirical null for top1 r over 300 random
screened genes and found NO usable threshold on r: median null top1_r = 0.270 and
88.7% of random query genes reach top1_r >= 0.20. The negative control OR5A1
returned a top partner at r=0.345, above SMARCA4's best true complex partner
(SMARCB1, 0.303). So the discriminator must be a property of the QUERY profile,
not of the correlation value.

This script tests candidate gates on the query's own Chronos profile against the
labelled panel (9 real complex queries, 1 junk control) and the 300-gene null,
and reports which separates them.

Run: python gate_codependency.py   (writes codependency_gate.csv)
"""
import numpy as np
import pandas as pd

CANDIDATES = ["profile_sd", "min_effect", "n_dep_lines"]


def load():
    panel = pd.read_csv("codependency_calibration_panel.csv")
    null = pd.read_csv("codependency_calibration_null.csv")
    panel["label"] = np.where(panel["query"] == "OR5A1", "control", "real")
    return panel, null


def separation(panel, null):
    """Per-metric separation between real panel queries, the control, and the null."""
    real = panel[panel["label"] == "real"]
    ctrl = panel[panel["label"] == "control"]
    rows = []
    for m in CANDIDATES:
        rows.append({
            "metric": m,
            "real_min": round(float(real[m].min()), 4),
            "real_median": round(float(real[m].median()), 4),
            "control_value": round(float(ctrl[m].iloc[0]), 4),
            "null_median": round(float(null[m].median()), 4),
            "null_p10": round(float(np.percentile(null[m], 10)), 4),
            "null_p90": round(float(np.percentile(null[m], 90)), 4),
        })
    return pd.DataFrame(rows)


def gate_performance(panel, null, sd_min, min_effect_max):
    """Fraction of each group passing a (profile_sd, min_effect) gate."""
    def frac(df):
        return float(((df["profile_sd"] >= sd_min)
                      & (df["min_effect"] <= min_effect_max)).mean())
    real = panel[panel["label"] == "real"]
    ctrl = panel[panel["label"] == "control"]
    return {
        "sd_min": sd_min, "min_effect_max": min_effect_max,
        "real_pass": round(frac(real), 3),
        "control_pass": round(frac(ctrl), 3),
        "null_pass": round(frac(null), 3),
    }


if __name__ == "__main__":
    panel, null = load()
    sep = separation(panel, null)
    pd.set_option("display.width", 220)
    print("METRIC SEPARATION")
    print(sep.to_string(index=False))

    grid = [gate_performance(panel, null, sd, me)
            for sd in (0.15, 0.20, 0.25)
            for me in (-0.75, -1.0, -1.5)]
    grid = pd.DataFrame(grid)
    print("\nGATE GRID (real_pass should be 1.0, control_pass 0.0)")
    print(grid.to_string(index=False))

    # Does the gate improve the r distribution it admits?
    passing = null[(null["profile_sd"] >= 0.20) & (null["min_effect"] <= -1.0)]
    failing = null.drop(passing.index)
    print("\nNULL top1_r by gate status")
    print("  passes gate n=%d  median top1_r %.3f" % (len(passing), passing["top1_r"].median()))
    print("  fails  gate n=%d  median top1_r %.3f" % (len(failing), failing["top1_r"].median()))
    print("  proximity frac top15: passes %.3f  fails %.3f"
          % (passing["prox_frac_top15"].median(), failing["prox_frac_top15"].median()))

    # Where does the LAST true partner sit, among gate-passing panel queries?
    partners = pd.read_csv("codependency_validation.csv")
    hits = partners[partners["expected_partner"]]
    per_query = hits.groupby("query")["r"].agg(["min", "max", "count"])
    print("\nCURATED PARTNER r RANGE within top15, per panel query")
    print(per_query.to_string())
    print("\nlowest true-partner r across all panel queries: %.4f" % hits["r"].min())
    print("highest r of a NON-expected top15 row: %.4f"
          % partners[~partners["expected_partner"] & (partners["query"] != "OR5A1")]["r"].max())

    sep.to_csv("codependency_gate.csv", index=False)
    grid.to_csv("codependency_gate_grid.csv", index=False)
