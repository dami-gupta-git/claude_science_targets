"""Quantify the cost of using the wrong background for term enrichment.

GeneTEA's corpus carries ~35.8k genes; the CRISPR release screens ~18.5k. A
hypergeometric enrichment run with background=None tests the hit list against
the whole corpus, so a term is credited partly for being SCREENABLE - well
annotated, protein coding, in a druggable family - rather than for being shared
by the codependency list. The correct background is the screened set.

Run: python background_comparison.py   (writes background_comparison*.csv)
"""
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

import kernel as cd

QUERIES = ["SMARCA4", "RAD21", "MRPL11", "BMS1", "VPS4A", "SOX10"]
LIST_SIZE = 50
N_TERMS = 25


def load_tea(path=None):
    """Load the trained GeneTEA model via the skill's own resolver.

    Delegates to kernel.genetea_load so the path comes from $GENETEA_MODEL (or
    an explicit path=) and a missing model raises FileNotFoundError naming the
    variable, rather than this script carrying a path from one machine.
    """
    return cd.genetea_load(path)


def corpus_genes(tea):
    """Symbols in the GeneTEA corpus, as a set."""
    return set(pd.Index(tea.entities).astype(str))


def enrich(tea, genes, background, n=N_TERMS):
    """get_enriched_terms wrapper returning an empty frame instead of None.

    GeneTEA returns None when nothing clears max_fdr, which breaks a naive
    len() or merge; the comparison needs a frame either way.
    """
    out = tea.get_enriched_terms(entities=list(genes), background=background, n=n)
    if out is None:
        return pd.DataFrame(columns=["Term", "FDR", "Effect Size"])
    return out.reset_index() if out.index.name else out.copy()


def compare_one(tea, query, screened, root=None):
    """Enrichment on one codependency list, whole-corpus vs screened background."""
    table = cd.depmap_codependencies(query, n=LIST_SIZE, root=root)
    hits = table["gene"].tolist()

    wide = enrich(tea, hits, background=None)
    narrow = enrich(tea, hits, background=sorted(screened))

    tcol = "Term" if "Term" in wide.columns else wide.columns[0]
    w_terms = list(wide[tcol]) if len(wide) else []
    n_terms = list(narrow[tcol]) if len(narrow) else []

    dropped = [t for t in w_terms if t not in set(n_terms)]
    gained = [t for t in n_terms if t not in set(w_terms)]
    shared = [t for t in w_terms if t in set(n_terms)]

    # Rank movement among terms present in both.
    w_rank = {t: i + 1 for i, t in enumerate(w_terms)}
    n_rank = {t: i + 1 for i, t in enumerate(n_terms)}
    moves = [abs(w_rank[t] - n_rank[t]) for t in shared]

    fdr_rows = []
    if len(wide) and len(narrow) and "FDR" in wide.columns:
        wf = wide.set_index(tcol)["FDR"]
        nf = narrow.set_index(tcol)["FDR"]
        for t in shared:
            fdr_rows.append({
                "query": query, "term": t,
                "fdr_corpus_background": float(wf[t]),
                "fdr_screened_background": float(nf[t]),
                "log10_fdr_shift": float(np.log10(max(float(nf[t]), 1e-300))
                                         - np.log10(max(float(wf[t]), 1e-300))),
                "rank_corpus": w_rank[t], "rank_screened": n_rank[t],
            })

    summary = {
        "query": query,
        "list_size": len(hits),
        "n_terms_corpus_background": len(w_terms),
        "n_terms_screened_background": len(n_terms),
        "n_shared": len(shared),
        "n_dropped_when_corrected": len(dropped),
        "n_gained_when_corrected": len(gained),
        "n_rank_changed": int(sum(1 for m in moves if m > 0)),
        "max_rank_move": int(max(moves)) if moves else 0,
        "top1_corpus": w_terms[0] if w_terms else "",
        "top1_screened": n_terms[0] if n_terms else "",
        "top1_changed": bool(w_terms and n_terms and w_terms[0] != n_terms[0]),
        "dropped_terms": " | ".join(dropped[:6]),
        "gained_terms": " | ".join(gained[:6]),
    }
    return summary, pd.DataFrame(fdr_rows)


if __name__ == "__main__":
    cd.codep_prepare()
    screened = cd.depmap_screened_genes(min_lines=cd.CODEP_MIN_LINES)
    screened_all = cd.depmap_screened_genes()
    tea = load_tea()
    corpus = corpus_genes(tea)

    overlap = corpus & set(screened_all)
    counts = {
        "genetea_corpus": len(corpus),
        "screened_all": len(screened_all),
        "screened_min_lines_%d" % cd.CODEP_MIN_LINES: len(screened),
        "corpus_and_screened": len(overlap),
        "screened_not_in_corpus": len(set(screened_all) - corpus),
        "corpus_not_screened": len(corpus - set(screened_all)),
        "background_shrink_factor": round(len(corpus) / max(len(overlap), 1), 3),
    }
    print("GENE COUNTS")
    for k, v in counts.items():
        print("  %-32s %s" % (k, v))
    pd.DataFrame([counts]).to_csv("background_gene_counts.csv", index=False)

    # Enrichment background must be the intersection: genes outside the corpus
    # cannot contribute to any term, so leaving them in inflates the background.
    bg = sorted(overlap)

    summaries, fdrs = [], []
    for q in QUERIES:
        s, f = compare_one(tea, q, bg)
        summaries.append(s)
        fdrs.append(f)
    summary = pd.DataFrame(summaries)
    fdr = pd.concat(fdrs, ignore_index=True) if any(len(f) for f in fdrs) else pd.DataFrame()

    pd.set_option("display.width", 240)
    print("\nBACKGROUND COMPARISON (top %d terms per query)" % N_TERMS)
    print(summary[["query", "n_terms_corpus_background", "n_terms_screened_background",
                   "n_shared", "n_dropped_when_corrected", "n_gained_when_corrected",
                   "n_rank_changed", "max_rank_move", "top1_changed"]].to_string(index=False))
    if len(fdr):
        print("\nFDR SHIFT for terms present under both backgrounds")
        print("  n term-query pairs %d" % len(fdr))
        print("  median log10 FDR shift (screened - corpus): %+.3f"
              % fdr["log10_fdr_shift"].median())
        print("  pairs where screened background is LESS significant: %d / %d"
              % (int((fdr["log10_fdr_shift"] > 0).sum()), len(fdr)))
        print("  max log10 FDR worsening %+.3f   max improvement %+.3f"
              % (fdr["log10_fdr_shift"].max(), fdr["log10_fdr_shift"].min()))
        fdr.to_csv("background_comparison_fdr.csv", index=False)
    summary.to_csv("background_comparison_summary.csv", index=False)
    print("\nDROPPED / GAINED examples")
    for _, r in summary.iterrows():
        print("  %-8s dropped: %s" % (r["query"], r["dropped_terms"] or "(none)"))
        print("  %-8s gained : %s" % ("", r["gained_terms"] or "(none)"))
