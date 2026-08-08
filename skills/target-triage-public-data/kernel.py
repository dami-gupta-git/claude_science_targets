"""Helpers for public-data target triage. See SKILL.md."""
import numpy as np
import pandas as pd

# Chronos gene effect is scaled so 0 = no effect and -1 = the median
# common-essential gene (SKILL.md, "Dependency"), which is what sets this value:
# it is a property of the DepMap scaling, not a cutoff fitted to any screen. Used
# by classify_gene_effect(); comparators on the same lines matter more than the
# threshold itself, since an absolute gene effect is uninterpretable alone.
DEPMAP_ESSENTIAL_THRESHOLD = -1.0

# Working columns sensitivity_scan writes into its own copy of the caller's
# frame. Named here so the collision guard and the confound label cannot drift
# apart from the strings actually used.
_PANEL_MEAN_COL = "_panel_mean"
_PANEL_MEAN_LOO_COL = "_panel_mean_loo"

# Declared so an all-dropped scan returns an empty frame with the real columns
# rather than a column-less one, on which res.r_partial raises AttributeError.
_SCAN_COLUMNS = ("drug", "n", "r", "p", "r_partial", "p_partial", "confound")
OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
GDSC_RELEASE_BASE = "https://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-8.4/"
CBIO_API = "https://www.cbioportal.org/api"


def ot_target(ensembl_id, url=None):
    """Open Targets tractability + DepMap essentiality for one Ensembl gene id.

    Returns the raw `target` dict: keys approvedSymbol, tractability,
    depMapEssentiality (tissueName -> screens with geneEffect, expression).
    """
    import json
    import urllib.request
    if url is None:
        url = OT_GRAPHQL
    q = ("query T($id:String!){target(ensemblId:$id){id approvedSymbol "
         "tractability{label modality value} "
         "depMapEssentiality{tissueName screens{cellLineName diseaseFromSource "
         "geneEffect expression}}}}")
    req = urllib.request.Request(
        url, data=json.dumps({"query": q, "variables": {"id": ensembl_id}}).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        payload = json.load(resp)
    # GraphQL reports query-level failures in `errors` with HTTP 200, and an
    # unrecognised Ensembl id comes back as a null target on an otherwise valid
    # response. Both must fail here, naming the id, rather than downstream as an
    # AttributeError inside ot_essentiality_frame.
    if payload.get("errors"):
        msgs = "; ".join(str(e.get("message", e)) for e in payload["errors"])
        raise RuntimeError(f"Open Targets GraphQL error for {ensembl_id}: {msgs}")
    target = (payload.get("data") or {}).get("target")
    if target is None:
        raise ValueError(
            f"Open Targets returned no target for {ensembl_id!r}; expected an Ensembl "
            "gene id such as 'ENSG00000162813'")
    return target


def ot_essentiality_frame(target_dict):
    """Flatten ot_target()'s depMapEssentiality into a tidy DataFrame.

    Columns: gene, tissue, cell, disease, effect, expr. Always those six, even
    when the target carries no essentiality data, so a caller can select columns
    on an empty frame.

    Tissue entries may omit `screens` just as the target may omit
    depMapEssentiality entirely; both are treated as empty rather than raising.
    """
    if not isinstance(target_dict, dict):
        raise TypeError(
            f"expected the target dict from ot_target(), got {type(target_dict).__name__}; "
            "ot_target() now raises on an unknown id rather than returning None")
    recs = []
    for tis in (target_dict.get("depMapEssentiality") or []):
        for s in (tis.get("screens") or []):
            recs.append(dict(gene=target_dict.get("approvedSymbol"),
                             tissue=tis.get("tissueName"),
                             cell=s.get("cellLineName"), disease=s.get("diseaseFromSource"),
                             effect=s.get("geneEffect"), expr=s.get("expression")))
    return pd.DataFrame(recs, columns=["gene", "tissue", "cell", "disease", "effect", "expr"])


def classify_gene_effect(effect):
    """Label a Chronos gene effect against the common-essential reference point.

    "essential" at or below DEPMAP_ESSENTIAL_THRESHOLD (-1.0, the median
    common-essential gene), "intermediate" below -0.5, "dispensable" above it.
    The boundaries come from the Chronos scaling described in SKILL.md rather
    than from a fit, so they are fixed rather than re-derivable from a screen.

    A label is not a substitute for comparators on the same cell lines: report
    the number and its reference genes alongside, since -0.6 in a panel where
    RRM1 reads -2.5 means something different than in a panel where it reads -0.9.
    Returns "unknown" for a missing value.
    """
    if effect is None:
        return "unknown"
    e = float(effect)
    if not np.isfinite(e):
        return "unknown"
    if e <= DEPMAP_ESSENTIAL_THRESHOLD:
        return "essential"
    if e <= -0.5:
        return "intermediate"
    return "dispensable"


def normalize_cell_name(s):
    """Strip punctuation/case so GDSC and CCLE line names join. 'NCI-H23' -> 'NCIH23'."""
    import re
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def partial_corr(x, y, z):
    """Pearson correlation of x,y after linearly removing z from both.

    THE load-bearing control for expression-vs-drug-sensitivity work: pass z =
    each cell line's mean response across the whole drug panel, so a gene that
    merely tracks proliferation rate does not read as a drug-specific effect.
    Returns (r, p).

    The p-value uses n-3 degrees of freedom, not the n-2 that scipy's pearsonr
    would apply to the residuals: fitting z out of each variable consumes a
    parameter. Borrowing pearsonr's p understates it in the anti-conservative
    direction (~0.81x at n=40), which is the wrong direction for a control whose
    purpose is to stop a confounded correlation being reported as real.

    Raises on a constant z rather than returning the uncontrolled correlation:
    polyfit is rank-deficient there, the residuals come back as the original
    variables, and the return value is then a plain Pearson r indistinguishable
    from a genuine partial one. A single-drug panel reaches this case.
    """
    from scipy import stats
    x, y, z = np.asarray(x, float), np.asarray(y, float), np.asarray(z, float)
    if not (x.shape == y.shape == z.shape):
        raise ValueError(f"x, y, z must be the same length; got {x.shape}, {y.shape}, {z.shape}")
    n = x.size
    if n < 4:
        raise ValueError(f"partial_corr needs n >= 4 to leave residual degrees of freedom; got {n}")
    n_bad = int((~np.isfinite(x) | ~np.isfinite(y) | ~np.isfinite(z)).sum())
    if n_bad:
        raise ValueError(
            f"partial_corr needs finite input; {n_bad} of {n} positions are NaN or inf. "
            "Drop them first and report the surviving n.")
    if np.ptp(z) == 0:
        raise ValueError(
            "confound z is constant, so there is nothing to partial out; the residuals "
            "would be the original variables and the result a plain Pearson r")
    rx = x - np.polyval(np.polyfit(z, x, 1), z)
    ry = y - np.polyval(np.polyfit(z, y, 1), z)
    r = float(np.corrcoef(rx, ry)[0, 1])
    if not np.isfinite(r):
        raise ValueError("residuals are constant after removing z; correlation is undefined")
    if abs(r) >= 1.0:
        return r, 0.0
    t = r * np.sqrt((n - 3) / (1 - r ** 2))
    return r, float(2 * stats.t.sf(abs(t), n - 3))


def sensitivity_scan(df, expr_col, drug_cols, panel_mean_col=None, min_n=15,
                     leave_one_out=True):
    """Correlate one expression column against many drug-response columns.

    df: one row per cell line, already subset to the lineage of interest.
    panel_mean_col: column holding each line's mean general-sensitivity level.
        If None it is computed from drug_cols. Always report r_partial, not r.
    leave_one_out: when the panel mean is computed here, exclude the drug being
        tested from its own confound. Including it regresses part of a drug's
        signal out of itself, attenuating r_partial toward zero by roughly
        1/len(drug_cols) — ~46% on a 3-drug panel, ~10% at 25, ~3% at 100.
        The bias is conservative (it cannot manufacture a hit) but it will hide
        a borderline one, so this defaults to True. Set False to reproduce a
        prior run, or supply panel_mean_col to use a fixed confound.
    Returns a DataFrame: drug, n, r, p, r_partial, p_partial, q_BH, q_partial_BH
    (BH-corrected across the drugs actually tested). Drugs that could not be
    tested are excluded from those rows but recorded in `.attrs["dropped"]` as
    (drug, n, reason) and summarised on stdout, so a scan that analysed three of
    twelve drugs says so rather than returning a short table silently.
    """
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    d = df.copy()
    drug_cols = list(drug_cols)
    supplied_panel = panel_mean_col is not None
    # Validate before any column is read or written, so a missing name fails with
    # the drug it names rather than with pandas' index message from the panel mean.
    missing_drugs = [c for c in drug_cols if c not in d.columns]
    if missing_drugs:
        raise KeyError(f"drug columns absent from df: {missing_drugs}")
    if expr_col not in d.columns:
        raise KeyError(f"expr_col {expr_col!r} is not a column of df")
    if not supplied_panel:
        if len(drug_cols) < 2:
            # With one drug column, "the whole-panel mean" IS that column: the
            # confound would equal the drug being tested, so partial_corr
            # would regress the drug on itself and correlate the expression
            # residual against ~1e-16 floating-point noise -- a fabricated
            # r_partial that looks like a legitimate null result instead of
            # the "nothing to control for" error it actually is.
            raise ValueError(
                "sensitivity_scan needs >=2 drug_cols to build a panel-mean "
                "confound from the data (a single drug's panel mean is itself); "
                "pass panel_mean_col= computed from a wider drug panel, or use "
                "confounder_check for a single-drug check.")
        for reserved in (_PANEL_MEAN_COL, _PANEL_MEAN_LOO_COL):
            if reserved in d.columns:
                raise ValueError(
                    f"{reserved!r} is reserved by sensitivity_scan for the panel-mean "
                    "confound and would overwrite the supplied column; rename it, or pass "
                    "panel_mean_col to use a fixed confound.")
        panel_mean_col = _PANEL_MEAN_COL
        d[panel_mean_col] = d[drug_cols].mean(axis=1)
    elif panel_mean_col not in d.columns:
        raise KeyError(
            f"panel_mean_col {panel_mean_col!r} is not a column of df; "
            f"available: {sorted(d.columns)[:20]}")
    out = []
    dropped = []
    for drug in drug_cols:
        # A leave-one-out confound needs its own column per drug; with a
        # user-supplied panel mean the constituent columns are unknown, so the
        # fixed column is used as given.
        if supplied_panel or not leave_one_out or len(drug_cols) < 2:
            conf_col = panel_mean_col
        else:
            conf_col = _PANEL_MEAN_LOO_COL
            d[conf_col] = d[[c for c in drug_cols if c != drug]].mean(axis=1)
        s = d[[drug, expr_col, conf_col]].dropna()
        if len(s) < min_n:
            dropped.append(dict(drug=drug, n=len(s), reason=f"n < min_n ({min_n})"))
            continue
        if s[conf_col].nunique() == 1:
            # A constant confound makes the partial correlation undefined; with a
            # leave-one-out confound on a small panel this is reachable whenever
            # the other drugs carry no variance on the surviving lines.
            dropped.append(dict(drug=drug, n=len(s), reason="confound constant on surviving rows"))
            continue
        r, p = stats.pearsonr(s[expr_col], s[drug])
        rp, pp = partial_corr(s[expr_col].values, s[drug].values, s[conf_col].values)
        out.append(dict(drug=drug, n=len(s), r=r, p=p, r_partial=rp, p_partial=pp,
                        confound="leave-one-out" if conf_col == _PANEL_MEAN_LOO_COL
                                 else "whole-panel"))
    if dropped:
        print(f"[sensitivity_scan] {len(dropped)} of {len(drug_cols)} drugs not tested: "
              + ", ".join(f"{r['drug']} (n={r['n']}, {r['reason']})" for r in dropped))
    res = pd.DataFrame(out, columns=_SCAN_COLUMNS)
    res.attrs["dropped"] = dropped
    if not len(res):
        return res
    res["q_BH"] = multipletests(res.p, method="fdr_bh")[1]
    res["q_partial_BH"] = multipletests(res.p_partial, method="fdr_bh")[1]
    res = res.sort_values("r_partial").reset_index(drop=True)
    # attrs propagation through sort/reset is version-dependent, so re-attach.
    res.attrs["dropped"] = dropped
    return res


def confounder_check(df, expr_col, drug_cols):
    """Is expression correlated with sensitivity to the ENTIRE panel?

    Run this BEFORE believing any single-drug hit. A strong result here means
    per-drug correlations are largely a proliferation/general-sensitivity
    artifact and only partial correlations are interpretable.
    Returns dict with r, p, n.

    Here the panel mean is the quantity under test rather than a confound, so
    the leave-one-out correction that sensitivity_scan applies does not arise.
    pandas skips NaN per row when averaging, so a line contributes as long as one
    drug column is measured; `n_drugs_min` reports the thinnest coverage behind
    any surviving row, since a panel mean resting on one of ten drugs is not the
    general-sensitivity measure the name implies.
    """
    from scipy import stats
    drug_cols = list(drug_cols)
    d = df.copy()
    missing = [c for c in drug_cols + [expr_col] if c not in d.columns]
    if missing:
        raise KeyError(f"columns absent from df: {missing}")
    if _PANEL_MEAN_COL in d.columns:
        raise ValueError(
            f"{_PANEL_MEAN_COL!r} is reserved by confounder_check for the panel mean "
            "and would overwrite the supplied column; rename it.")
    d[_PANEL_MEAN_COL] = d[drug_cols].mean(axis=1)
    d["_n_drugs"] = d[drug_cols].notna().sum(axis=1)
    s = d[[expr_col, _PANEL_MEAN_COL, "_n_drugs"]].dropna(subset=[expr_col, _PANEL_MEAN_COL])
    if len(s) < 3:
        raise ValueError(
            f"confounder_check needs at least 3 complete rows; got {len(s)} "
            f"from {len(d)} rows of {len(drug_cols)} drug columns")
    r, p = stats.pearsonr(s[expr_col], s[_PANEL_MEAN_COL])
    return {"r": float(r), "p": float(p), "n": int(len(s)),
            "n_drugs": len(drug_cols), "n_drugs_min": int(s["_n_drugs"].min())}


def load_prism(matrix_path, compound_list_path, drug_names=None):
    """Load PRISM Repurposing into a cell-line x drug DataFrame indexed by ACH- id.

    The matrix's row ids are BRD- broad ids and are meaningless without the
    compound list, which maps them to Drug.Name. Replicate rows for the same
    drug are averaged. drug_names: optional uppercase set to keep.

    The matrix header supplies the index, so duplicate cell-line ids there are
    rejected: a duplicated index turns `.loc[line, drug]` into a Series and
    multiplies rows in any downstream merge. Short data rows are rejected for the
    same reason a truncated file must not load quietly — DepMap's portal is behind
    bot verification, so these files are downloaded by hand and arriving
    incomplete is a realistic outcome.
    """
    import csv
    want = {}
    with open(compound_list_path) as fh:
        wanted_upper = {d.upper() for d in drug_names} if drug_names else None
        for r in csv.DictReader(fh):
            if wanted_upper is not None and r["Drug.Name"].upper() not in wanted_upper:
                continue
            for i in r["IDs"].split(","):
                want[i.strip()] = r["Drug.Name"]
    if not want:
        raise ValueError(
            f"no compounds selected from {compound_list_path}"
            + (f" for drug_names={sorted(drug_names)}" if drug_names else ""))
    with open(matrix_path) as f:
        rd = csv.reader(f)
        try:
            lines = next(rd)[1:]
        except StopIteration:
            raise ValueError(f"{matrix_path} is empty; expected a header of cell-line ids")
        dupes = sorted({c for c in lines if lines.count(c) > 1})
        if dupes:
            raise ValueError(
                f"{matrix_path} header repeats {len(dupes)} cell-line id(s): {dupes[:10]}. "
                "A duplicated index multiplies rows on merge; de-duplicate the matrix first.")
        cols = {}
        for lineno, row in enumerate(rd, start=2):
            if not row:
                continue
            if len(row) - 1 != len(lines):
                raise ValueError(
                    f"{matrix_path} line {lineno} (id {row[0]!r}) has {len(row) - 1} values "
                    f"but the header declares {len(lines)} cell lines; the file looks truncated")
            if row[0] in want:
                cols.setdefault(want[row[0]], []).append(
                    pd.to_numeric(pd.Series(row[1:], index=lines), errors="coerce"))
    if not cols:
        raise ValueError(
            f"none of the {len(want)} selected BRD ids appear in {matrix_path}; "
            "check that the compound list and matrix come from the same release")
    return pd.DataFrame({d: pd.concat(v, axis=1).mean(axis=1) for d, v in cols.items()})


def reconcile_fetch(expected_ids, fetched_ids, label="records"):
    """Assert a batched API fetch returned everything the search promised.

    Batched metadata fetches can silently drop ids. Call this after every
    search->fetch pair; returns the sorted missing ids so you can refetch.
    """
    missing = sorted(set(expected_ids) - set(fetched_ids))
    if missing:
        print(f"[reconcile] {label}: {len(missing)} of {len(set(expected_ids))} missing -> {missing}")
    return missing


def load_prism_secondary(path, drug_names=None, min_r2=None):
    """Load PRISM 20Q2 secondary dose-response into a cell-line x drug AUC frame.

    This is the 8-point dilution screen (fitted AUC/IC50), NOT the single-dose primary
    matrix - prefer it when your compound is present. Only primary-screen hits were
    carried forward, so check membership first. Note the polarity: higher AUC = MORE
    RESISTANT, opposite to the primary screen's log-fold-change.
    Returns a DataFrame indexed by ModelID (ACH- id), columns = lowercase drug names.

    `r2` is read only when min_r2 is given, so a file lacking that column still
    loads when no fit-quality filter was asked for.
    """
    cols = ["depmap_id", "name", "auc"] + (["r2"] if min_r2 is not None else [])
    try:
        df = pd.read_csv(path, usecols=cols, low_memory=False)
    except ValueError as e:
        raise ValueError(f"{path} is missing a required column ({cols}): {e}") from e
    df["name"] = df["name"].str.lower()
    if drug_names:
        df = df[df.name.isin({d.lower() for d in drug_names})]
        if df.empty:
            raise ValueError(
                f"none of drug_names={sorted(drug_names)} appear in {path}; "
                "check prism_secondary_has() first — only primary-screen hits were "
                "carried into the dose-response screen")
    if min_r2 is not None:
        df = df[df.r2 >= min_r2]
    out = df.groupby(["depmap_id", "name"])["auc"].mean().unstack()
    out.index.name = "ModelID"
    return out


def prism_secondary_has(path, drug_names):
    """Which of drug_names were carried into the PRISM secondary screen? Returns a dict."""
    import csv
    present = set()
    with open(path) as f:
        rd = csv.DictReader(f)
        for r in rd:
            n = (r.get("name") or "").lower()
            if n:
                present.add(n)
    return {d: (d.lower() in present) for d in drug_names}
