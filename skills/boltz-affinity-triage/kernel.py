
"""Helpers for KB-grounded Boltz-2 affinity triage. See SKILL.md."""
import re
import os
import csv
import json
import math
import shutil

AFFINITY_ATOM_CAP = 128
NM_TO_PACT = 9.0  # -log10(1e-9)
# Thresholds are compared with a small tolerance: a margin that exactly meets
# min_margin can land a few ulps low in binary floating point (0.60 - 0.50
# evaluates to 0.09999999999999998), which would reject a run that satisfies
# the documented requirement.
GATE_FLOAT_TOL = 1e-9


def kba_parse_nm(value):
    """BindingDB affinity is a STRING that may carry '>' / '<' qualifiers.

    Returns (nanomolar_float, qualifier) with qualifier in {'=', '>', '<'},
    or (None, None) if unparseable. Qualified values are censored
    measurements -- a '>10000' is a real inactive, a '<1' is a ceiling hit
    whose true value is unknown, so callers must not treat them as equalities.
    """
    if value is None:
        return None, None
    text = str(value).strip()
    if not text:
        return None, None
    qual = "="
    if text[0] in "><":
        qual = text[0]
        text = text[1:].strip()
    text = text.replace(",", "")
    try:
        return float(text), qual
    except ValueError:
        return None, None


def kba_pactivity(row, source="chembl"):
    """Unified pActivity (-log10 molar) from a ChEMBL or BindingDB row.

    ChEMBL's pchembl_value is absent whenever the measurement is qualified or
    non-nM, which in practice is most of the weak/inactive tail -- so falling
    back to standard_value is required or the negatives silently vanish.
    Returns (pact, qualifier, unit_source) or (None, None, None).
    """
    if source == "bindingdb":
        nm, qual = kba_parse_nm(row.get("affinity"))
        if nm is None or nm <= 0:
            return None, None, None
        return NM_TO_PACT - math.log10(nm), qual, "affinity"
    pch = row.get("pchembl_value")
    if pch not in (None, ""):
        try:
            return float(pch), row.get("standard_relation") or "=", "pchembl_value"
        except (TypeError, ValueError):
            pass
    val = row.get("standard_value")
    units = (row.get("standard_units") or "").strip()
    scale = {"pM": 1e-12, "nM": 1e-9, "uM": 1e-6, "M": 1.0, "mM": 1e-3}.get(units)
    if val in (None, "") or scale is None:
        return None, None, None
    try:
        molar = float(val) * scale
    except (TypeError, ValueError):
        return None, None, None
    if molar <= 0:
        return None, None, None
    return -math.log10(molar), row.get("standard_relation") or "=", "standard_value"


def kba_heavy_atoms(smiles):
    """Heavy-atom count. Uses RDKit when importable, else a regex estimate."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            return mol.GetNumHeavyAtoms()
    except Exception:
        pass
    stripped = re.sub(r"\[[^\]]*\]", "X", smiles or "")
    stripped = stripped.replace("Cl", "L").replace("Br", "R")
    return len(re.findall(r"[CNOSPFIBLRXcnospb]", stripped))


def kba_props(smiles):
    """Property vector for decoy matching: heavy_atoms, mw, clogp, rings.

    mw/clogp are None without RDKit; matching then falls back to heavy atoms
    alone, which still removes the dominant size confounder.
    """
    out = {"smiles": smiles, "heavy_atoms": kba_heavy_atoms(smiles),
           "mw": None, "clogp": None, "rings": None}
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            out["mw"] = Descriptors.MolWt(mol)
            out["clogp"] = Crippen.MolLogP(mol)
            out["rings"] = Chem.rdMolDescriptors.CalcNumRings(mol)
    except Exception:
        pass
    return out


def kba_match_decoys(actives, pool, per_active=1, tol=None):
    """Pick decoys whose bulk properties mirror the actives.

    Unmatched actives/decoys differ ~2x in heavy-atom count on real targets,
    so any size-correlated scorer separates them and the calibration passes
    for the wrong reason. Greedy nearest-neighbour on (heavy_atoms, mw),
    sampling without replacement.
    Both args: lists of dicts carrying 'smiles'. Returns chosen pool entries.
    """
    tol = {"heavy_atoms": 6, "mw": 75.0, **(tol or {})}
    act = [kba_props(a["smiles"]) for a in actives if a.get("smiles")]
    cand = [(i, kba_props(p["smiles"])) for i, p in enumerate(pool) if p.get("smiles")]
    used, chosen = set(), []
    for a in act:
        scored = []
        for i, c in cand:
            if i in used:
                continue
            dha = abs(c["heavy_atoms"] - a["heavy_atoms"])
            if dha > tol["heavy_atoms"]:
                continue
            dist = dha / max(tol["heavy_atoms"], 1)
            if a["mw"] is not None and c["mw"] is not None:
                dmw = abs(c["mw"] - a["mw"])
                if dmw > tol["mw"]:
                    continue
                dist += dmw / tol["mw"]
            scored.append((dist, i))
        scored.sort()
        for _, i in scored[:per_active]:
            used.add(i)
            chosen.append(pool[i])
    return chosen


def kba_leakage_flag(row, cutoff_year=2024):
    """Flag a ChEMBL row as likely-memorized by the affinity head.

    Boltz-2's affinity training drew on ChEMBL v34 / BindingDB / PubChem as of
    June 2023, so a measurement from an older publication is plausibly recalled
    rather than predicted. document_year is the cheapest available proxy.
    Returns 'likely_train', 'likely_novel', or 'unknown'.
    """
    year = row.get("document_year")
    try:
        year = int(year)
    except (TypeError, ValueError):
        return "unknown"
    return "likely_novel" if year >= cutoff_year else "likely_train"


def kba_write_yaml(path, target_sequence, smiles, ligand_id="L", protein_id="A",
                   msa=None, request_affinity=True):
    """Write one Boltz-2 input YAML. Returns the path, or None if over the cap.

    Affinity needs a ligand chain named in properties:, and Boltz refuses
    ligands at/above 128 atoms, so oversized queries are skipped rather than
    silently producing structure-only output.
    """
    n_heavy = kba_heavy_atoms(smiles)
    if request_affinity and n_heavy >= AFFINITY_ATOM_CAP:
        return None
    msa_line = "      msa: {}\n".format(msa) if msa else ""
    text = (
        "version: 1\nsequences:\n"
        "  - protein:\n      id: {pid}\n      sequence: {seq}\n{msa}"
        "  - ligand:\n      id: {lid}\n      smiles: '{smi}'\n"
    ).format(pid=protein_id, seq=target_sequence, msa=msa_line,
             lid=ligand_id, smi=smiles)
    if request_affinity:
        text += "properties:\n  - affinity:\n      binder: {}\n".format(ligand_id)
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def kba_read_result(pred_dir, model_index=0, optimization_score_direction=None):
    """Read Boltz affinity + confidence JSONs for ONE ranked diffusion sample.

    With --diffusion_samples N, Boltz writes confidence_<name>_model_K.json for
    each sample K, ranked so that model_0 is its best. Confidence must therefore
    be read from a named sample rather than from whichever file the filesystem
    lists first, or the reported ipTM belongs to an arbitrary sample and the
    ipTM < 0.5 pose check in SKILL.md step 4 can flag a good complex or clear a
    bad one. Affinity is predicted once per input, so affinity_*.json is not
    per-sample and is read whole.

    Returns 'binder_score' (higher = more likely a binder, normalised across
    Boltz versions), the selected sample's confidence, 'model_index' actually
    used, 'n_samples' found, and 'per_sample' mapping index -> confidence so a
    caller can inspect sampling spread.

    Boltz-2.1 replaced the affinity outputs with 'binding_confidence' and
    'optimization_score'. binding_confidence is the binder-vs-decoy analogue of
    affinity_probability_binary and is used for ranking directly.
    optimization_score is NOT used for ranking unless
    optimization_score_direction is given as 'higher_is_better' or
    'lower_is_better', because its sign convention is not documented and the
    published descriptions of the analogous Boltz-2 quantity disagree: the model
    docs define affinity_pred_value as log10(IC50 in uM), where lower binds
    tighter, while third-party summaries describe the same field as a pIC50,
    where higher binds tighter. Ranking a screen with the sign inverted would
    order every compound backwards while still producing a plausible-looking
    table, so an output carrying only optimization_score returns
    binder_score=None with 'unscored_reason' naming the problem. The raw value
    is always reported under 'optimization_score' so nothing is discarded.
    """
    out = {"pred_dir": pred_dir, "binder_score": None, "score_key": None,
           "affinity_pred_value": None, "optimization_score": None,
           "unscored_reason": None, "iptm": None, "ptm": None,
           "complex_plddt": None, "confidence_score": None,
           "model_index": None, "n_samples": 0, "per_sample": {}, "raw": {}}
    if not os.path.isdir(pred_dir):
        return out
    conf_re = re.compile(r"^confidence_.*_model_(\d+)\.json$")
    affinity = {}
    for name in sorted(os.listdir(pred_dir)):
        if not name.endswith(".json"):
            continue
        try:
            with open(os.path.join(pred_dir, name)) as fh:
                blob = json.load(fh)
        except (ValueError, OSError):
            continue
        if not isinstance(blob, dict):
            continue
        out["raw"][name] = blob
        hit = conf_re.match(name)
        if hit:
            out["per_sample"][int(hit.group(1))] = blob
        elif name.startswith("affinity"):
            affinity.update(blob)
    out["n_samples"] = len(out["per_sample"])
    chosen = {}
    if out["per_sample"]:
        pick = model_index if model_index in out["per_sample"] else min(out["per_sample"])
        out["model_index"] = pick
        chosen = out["per_sample"][pick]
    for key in ("ptm", "iptm", "complex_plddt", "confidence_score"):
        if key in chosen:
            out[key] = chosen[key]
    if "optimization_score" in affinity:
        out["optimization_score"] = affinity["optimization_score"]
    if "affinity_pred_value" in affinity:
        out["affinity_pred_value"] = affinity["affinity_pred_value"]
    if "affinity_probability_binary" in affinity:
        out["binder_score"] = affinity["affinity_probability_binary"]
        out["score_key"] = "affinity_probability_binary"
    elif "binding_confidence" in affinity:
        out["binder_score"] = affinity["binding_confidence"]
        out["score_key"] = "binding_confidence"
    elif "affinity_pred_value" in affinity:
        out["binder_score"] = -affinity["affinity_pred_value"]
        out["score_key"] = "neg_affinity_pred_value"
    elif "optimization_score" in affinity:
        if optimization_score_direction == "higher_is_better":
            out["binder_score"] = affinity["optimization_score"]
            out["score_key"] = "optimization_score"
        elif optimization_score_direction == "lower_is_better":
            out["binder_score"] = -affinity["optimization_score"]
            out["score_key"] = "neg_optimization_score"
        else:
            out["unscored_reason"] = (
                "only optimization_score is present and its sign convention is "
                "undocumented; pass optimization_score_direction="
                "'higher_is_better' or 'lower_is_better' after confirming it "
                "against a known potent/weak pair on this target")
    elif affinity:
        out["unscored_reason"] = (
            "affinity file present but carried no recognised score key: %s"
            % ", ".join(sorted(affinity)))
    return out


def kba_auc(pos, neg):
    """Mann-Whitney AUC with tie correction. Higher score = binder."""
    pos = [p for p in pos if p is not None]
    neg = [n for n in neg if n is not None]
    if not pos or not neg:
        return None
    merged = sorted([(v, 0) for v in neg] + [(v, 1) for v in pos])
    ranks, i = {}, 0
    while i < len(merged):
        j = i
        while j + 1 < len(merged) and merged[j + 1][0] == merged[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks.setdefault(k, avg)
        i = j + 1
    rsum = sum(ranks[k] for k, (_, lab) in enumerate(merged) if lab == 1)
    n1, n0 = len(pos), len(neg)
    return (rsum - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def kba_enrichment(scored, frac=0.1):
    """Enrichment factor at the top `frac` of a ranked list.

    scored: list of (score, is_active). Returns None when the slice is empty.
    """
    rows = [(s, a) for s, a in scored if s is not None]
    if not rows:
        return None
    rows.sort(key=lambda r: -r[0])
    n_top = max(1, int(round(len(rows) * frac)))
    base = sum(1 for _, a in rows if a) / float(len(rows))
    if base == 0:
        return None
    top = sum(1 for _, a in rows[:n_top] if a) / float(n_top)
    return top / base


def kba_percentile(value, reference):
    """Where a query score falls within the control distribution (0-100)."""
    ref = sorted(r for r in reference if r is not None)
    if value is None or not ref:
        return None
    below = sum(1 for r in ref if r < value)
    ties = sum(1 for r in ref if r == value)
    return 100.0 * (below + 0.5 * ties) / len(ref)


def kba_null_auc(actives, decoys):
    """AUC of a size-only scorer on the same controls -- the confounder check.

    If molecular size alone separates actives from decoys, then a scorer that
    merely tracks size scores well without knowing anything about binding.
    Report this next to the real AUC: a real AUC that does not clear the null
    by a margin is not evidence of binding discrimination.
    """
    pos = [kba_heavy_atoms(a.get("smiles")) for a in actives if a.get("smiles")]
    neg = [kba_heavy_atoms(d.get("smiles")) for d in decoys if d.get("smiles")]
    return kba_auc(pos, neg)


def kba_gate(auc, n_pos, n_neg, null_auc=None, min_auc=0.65,
             min_per_class=5, min_margin=0.1):
    """Decide whether this run's query rankings mean anything.

    Three ways a calibration fails, and all three must be reported rather
    than silently passed over:
      * too few controls -- AUC on 3 vs 3 compounds is noise;
      * AUC at chance -- no demonstrated discrimination on this target;
      * AUC not clearing the size-only null -- the separation is explained
        by a confounder, so it is not evidence about binding.
    Emitting a ranked hit list from a failed calibration would present noise
    as a result, so the honest output is the diagnostic, not the ranking.
    """
    reasons = []
    if n_pos < min_per_class or n_neg < min_per_class:
        reasons.append(
            "too few controls ({} active / {} decoy; need {} each)".format(
                n_pos, n_neg, min_per_class))
    if auc is None:
        reasons.append("control AUC not computable")
    else:
        if auc < min_auc - GATE_FLOAT_TOL:
            reasons.append("control AUC {:.2f} below {:.2f}".format(auc, min_auc))
        if null_auc is not None and (auc - null_auc) < min_margin - GATE_FLOAT_TOL:
            reasons.append(
                "control AUC {:.2f} does not clear size-only null {:.2f} "
                "by {:.2f} -- separation may be a size artifact".format(
                    auc, null_auc, min_margin))
    return {"pass": not reasons, "auc": auc, "null_auc": null_auc,
            "n_pos": n_pos, "n_neg": n_neg, "reasons": reasons,
            "verdict": "interpretable" if not reasons else "not interpretable"}


# --------------------------------------------------------------------------
# Run outputs: where a triage writes its files, and what it writes.
# --------------------------------------------------------------------------

RESULTS_TOPIC = "boltz_affinity_triage"
IPTM_FLOOR = 0.5          # below this the pose is reported as unreliable
SUMMARY_MAX_WORDS = 130   # the opening paragraph is an orientation, not a report

CALIBRATION_JSON = "calibration.json"
CONTROLS_CSV = "controls.csv"
QUERIES_CSV = "ranked_queries.csv"
RUN_README = "README.md"

CONTROL_COLUMNS = ("compound_id", "role", "smiles", "source", "pactivity",
                   "qualifier", "leakage_flag", "heavy_atoms", "binder_score",
                   "score_key", "iptm", "ptm", "model_index")
QUERY_COLUMNS = ("compound_id", "smiles", "binder_score", "score_key",
                 "control_percentile", "iptm", "ptm", "model_index",
                 "nearest_active", "nearest_tanimoto", "flags")


def kba_run_dir(target, root="results", topic=None, make=True):
    """Path for one triage run: <root>/<topic>/<target>/, with scripts/ beside it.

    Runs are siblings under one topic so that two targets stay comparable and
    the topic directory does not accumulate loose files. `target` is slugged to
    snake_case because result directories are named that way and a raw gene or
    ChEMBL label may carry spaces, hyphens or case.
    """
    topic = RESULTS_TOPIC if topic is None else topic
    slug = re.sub(r"[^a-z0-9]+", "_", str(target).strip().lower()).strip("_")
    if not slug:
        raise ValueError("kba_run_dir got an empty target name; the run "
                         "directory is named after the target")
    out_dir = os.path.join(root, topic, slug)
    if make:
        os.makedirs(os.path.join(out_dir, "scripts"), exist_ok=True)
    return out_dir


def kba_write_table(path, rows, headers=None):
    """Write rows (list of dicts) to CSV. Returns the path, or None if empty.

    Columns are `headers` first, then any further keys present in the rows in
    sorted order, so an extra field a caller attaches is written rather than
    dropped. None is written as an empty cell.
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
            writer.writerow({c: ("" if row.get(c) is None else row.get(c))
                             for c in cols})
    return path


def kba_tanimoto(smiles_a, smiles_b, radius=2, n_bits=2048):
    """Morgan-fingerprint Tanimoto between two SMILES, or None without RDKit.

    Returns None rather than a fallback similarity: a string-based substitute
    would order compounds by shared substring rather than by chemotype, and the
    report uses this number to say whether a hit is a novel scaffold or an
    analogue of a known active.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdFingerprintGenerator
        from rdkit import DataStructs
    except ImportError:
        return None
    mol_a, mol_b = Chem.MolFromSmiles(smiles_a or ""), Chem.MolFromSmiles(smiles_b or "")
    if mol_a is None or mol_b is None:
        return None
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    return DataStructs.TanimotoSimilarity(gen.GetFingerprint(mol_a),
                                          gen.GetFingerprint(mol_b))


def kba_nearest_active(smiles, actives):
    """Closest known active to a query compound: (compound_id, tanimoto).

    Returns (None, None) when RDKit is absent or nothing is comparable, so the
    report can say the column was not computed instead of implying the query is
    a novel chemotype.
    """
    best_id, best_sim = None, None
    for act in actives or []:
        sim = kba_tanimoto(smiles, act.get("smiles"))
        if sim is None:
            continue
        if best_sim is None or sim > best_sim:
            best_id, best_sim = act.get("compound_id") or act.get("id"), sim
    return best_id, best_sim


def kba_pose_flags(row, iptm_floor=None):
    """Caveat strings for one scored compound: low confidence, unscored, no pose.

    ipTM is read for a complex and ptm for a single chain, so whichever is
    present is compared against the floor. A low value is a caveat to report
    rather than grounds for discarding the row, because the published
    evaluations found affinity predictions that survive a poor pose.
    """
    iptm_floor = IPTM_FLOOR if iptm_floor is None else iptm_floor
    flags = []
    conf = row.get("iptm")
    conf = row.get("ptm") if conf is None else conf
    if conf is None:
        flags.append("no pose confidence")
    elif conf < iptm_floor:
        flags.append("low pose confidence ({:.2f} < {:.2f})".format(conf, iptm_floor))
    if row.get("binder_score") is None:
        reason = row.get("unscored_reason")
        flags.append("unscored" + (": " + reason if reason else ""))
    return flags


def kba_check_words(text, cap, label):
    n = len(str(text or "").split())
    if n > cap:
        raise ValueError("{} is {} words; cap is {}. Move detail into the "
                         "tables rather than the prose".format(label, n, cap))
    return n


def kba_markdown_table(rows, headers=None):
    """Render rows (list of dicts) as a markdown table."""
    rows = [dict(r) for r in (rows or [])]
    if not rows:
        return ""
    headers = list(headers or [])
    cols = [c for c in headers if any(c in r for r in rows)]
    cols += sorted({k for r in rows for k in r} - set(cols))
    fmt = lambda v: "" if v is None else ("{:.3g}".format(v) if isinstance(v, float) else str(v))
    lines = ["| " + " | ".join(cols) + " |",
             "| " + " | ".join("---" for _ in cols) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(fmt(row.get(c)) for c in cols) + " |")
    return "\n".join(lines)


def kba_run_readme(target, gate, summary, controls, queries=None, files=(),
                   data_sources=(), limits=(), enrichment=None,
                   leakage_split=None, skipped=(), title=None):
    """Render the run README. Returns markdown text.

    Sections follow `coding-standards` (Result, Files, Data sources, Limits)
    and the Result section leads with the calibration verdict, because every
    other statement in the document is conditional on it. When the gate fails
    the ranked table is replaced by the diagnostic and the reasons, matching
    what kba_write_run writes to disk.

    `gate` is a kba_gate() dict; `controls` and `queries` are the row lists
    written to the CSVs; `skipped` names compounds dropped at the 128-atom
    affinity cap. `summary` is plain prose capped at SUMMARY_MAX_WORDS.
    """
    if not isinstance(gate, dict) or "pass" not in gate:
        raise ValueError("kba_run_readme needs the dict returned by kba_gate; "
                         "got {!r}".format(gate))
    if not str(summary or "").strip():
        raise ValueError("kba_run_readme needs a summary — the plain-prose "
                         "paragraph a non-specialist reads first")
    kba_check_words(summary, SUMMARY_MAX_WORDS, "summary")
    if str(summary or "").strip().startswith("-"):
        raise ValueError("summary must be plain prose, not bullets")

    passed = bool(gate["pass"])
    auc, null_auc = gate.get("auc"), gate.get("null_auc")
    num = lambda v: "not computed" if v is None else "{:.2f}".format(v)

    parts = [f"# {title or f'{target} affinity triage'}", "",
             str(summary).strip(), "", "## Result", "", "### Calibration", ""]
    parts += [
        "| quantity | value |", "| --- | --- |",
        f"| target | {target} |",
        f"| controls | {gate.get('n_pos')} actives / {gate.get('n_neg')} decoys |",
        f"| control AUC | {num(auc)} |",
        f"| size-only null AUC | {num(null_auc)} |",
        f"| margin over null | {'not computed' if (auc is None or null_auc is None) else '%.2f' % (auc - null_auc)} |",
        f"| EF@10% | {'not computed' if enrichment is None else '%.2f' % enrichment} |",
        f"| verdict | {gate.get('verdict')} |", ""]
    if gate.get("reasons"):
        parts += ["The gate failed for the following reasons."] + \
                 [f"- {r}" for r in gate["reasons"]] + [""]
    if leakage_split:
        parts += ["Leakage split: " + ", ".join(
            f"{v} {k}" for k, v in sorted(leakage_split.items())) + ".", ""]

    parts += ["### Ranked queries", ""]
    if passed:
        if queries:
            parts += [kba_markdown_table(queries, QUERY_COLUMNS), ""]
        else:
            parts += ["No query compounds were submitted; this run is a "
                      "calibration only.", ""]
    else:
        parts += ["Not reported. The controls do not show discriminative power "
                  "on this target, so a ranking from this run would order "
                  "compounds by noise. The diagnostic above is the result.", ""]

    parts += ["## Files", ""]
    for entry in files:
        name, description = (entry["name"], entry["description"]) \
            if isinstance(entry, dict) else entry
        parts.append(f"- `{name}` — {description}")
    parts += ["", "## Data sources", ""] + [f"- {s}" for s in data_sources]
    parts += ["", "## Limits", ""]
    standing = ["Boltz-2 affinity output is an ordering, not a potency: "
                "`affinity_pred_value` is log10(IC50 in µM) on an assay-agnostic "
                "scale and is not quoted as a measured affinity."]
    low_conf = [r for r in (queries or []) + (controls or [])
                if any("low pose confidence" in f for f in kba_pose_flags(r))]
    if low_conf:
        standing.append("{} compound{} scored on a pose below ipTM {:.1f}: {}."
                        .format(len(low_conf), "" if len(low_conf) == 1 else "s",
                                IPTM_FLOOR,
                                ", ".join(str(r.get("compound_id")) for r in low_conf)))
    if queries and all(q.get("nearest_tanimoto") is None for q in queries):
        standing.append("Similarity to known actives was not computed; RDKit was "
                        "not available, so the nearest-active columns are empty "
                        "rather than indicating novel chemotypes.")
    if skipped:
        standing.append("Skipped at the {}-atom affinity cap: {}."
                        .format(AFFINITY_ATOM_CAP, ", ".join(str(s) for s in skipped)))
    parts += [f"- {s}" for s in standing + list(limits)]
    return "\n".join(parts).rstrip() + "\n"


def kba_write_run(out_dir, target, gate, controls, queries=None, summary=None,
                  files=(), data_sources=(), limits=(), enrichment=None,
                  leakage_split=None, skipped=(), scripts=()):
    """Write one run directory. Returns {name: path} for what was written.

    The ranked query table is written only when `gate["pass"]` is true. A
    ranking from a failed calibration is noise in the shape of a result, and a
    CSV on disk outlives the caveat in the chat that accompanied it, so the
    failing run's deliverable is the calibration and the controls alone.

    Query percentiles and pose flags are filled in here from the control
    distribution when the caller has not already set them, so that every run
    reports them the same way. `scripts` names files copied into scripts/.
    """
    if not isinstance(gate, dict) or "pass" not in gate:
        raise ValueError("kba_write_run needs the dict returned by kba_gate; "
                         "got {!r}".format(gate))
    controls = [dict(r) for r in (controls or [])]
    for row in controls:
        role = row.get("role")
        if role not in ("active", "decoy"):
            raise ValueError("control {!r} has role {!r}; every control row "
                             "must be 'active' or 'decoy'"
                             .format(row.get("compound_id"), role))
    queries = [dict(r) for r in (queries or [])]
    reference = [r.get("binder_score") for r in controls]
    actives = [r for r in controls if r["role"] == "active"]
    for row in queries:
        if row.get("control_percentile") is None:
            row["control_percentile"] = kba_percentile(row.get("binder_score"), reference)
        if row.get("nearest_active") is None and row.get("smiles"):
            near_id, near_sim = kba_nearest_active(row["smiles"], actives)
            row["nearest_active"], row["nearest_tanimoto"] = near_id, near_sim
        if not row.get("flags"):
            row["flags"] = "; ".join(kba_pose_flags(row))

    os.makedirs(os.path.join(out_dir, "scripts"), exist_ok=True)
    written = {}
    calibration = {"target": target, "gate": gate, "enrichment_at_10pct": enrichment,
                   "leakage_split": leakage_split or {},
                   "skipped_over_atom_cap": list(skipped),
                   "n_queries": len(queries),
                   "ranking_written": bool(gate["pass"]) and bool(queries)}
    cal_path = os.path.join(out_dir, CALIBRATION_JSON)
    with open(cal_path, "w", encoding="utf-8") as fh:
        json.dump(calibration, fh, indent=2, sort_keys=True)
    written["calibration"] = cal_path

    path = kba_write_table(os.path.join(out_dir, CONTROLS_CSV), controls, CONTROL_COLUMNS)
    if path:
        written["controls"] = path
    if gate["pass"]:
        path = kba_write_table(os.path.join(out_dir, QUERIES_CSV), queries, QUERY_COLUMNS)
        if path:
            written["queries"] = path

    for src in scripts:
        dst = os.path.join(out_dir, "scripts", os.path.basename(src))
        shutil.copyfile(src, dst)
        written.setdefault("scripts", []).append(dst)

    listed = [(CALIBRATION_JSON, "gate verdict, control counts, AUC, "
                                 "size-only null AUC and skipped compounds"),
              (CONTROLS_CSV, "one row per control compound with its knowledge-base "
                             "activity, leakage flag and co-folding scores")]
    if "queries" in written:
        listed.append((QUERIES_CSV, "query compounds ranked by binder_score, with "
                                    "control percentile, pose confidence and "
                                    "nearest known active"))
    listed += list(files)
    readme = kba_run_readme(target, gate, summary, controls,
                            queries=queries if "queries" in written else None,
                            files=listed, data_sources=data_sources, limits=limits,
                            enrichment=enrichment, leakage_split=leakage_split,
                            skipped=skipped)
    readme_path = os.path.join(out_dir, RUN_README)
    with open(readme_path, "w", encoding="utf-8") as fh:
        fh.write(readme)
    written["readme"] = readme_path
    return written
