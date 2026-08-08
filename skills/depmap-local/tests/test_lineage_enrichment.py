"""depmap_lineage_enrichment's per-lineage comparator group.

A lineage covering every scored line leaves nothing to compare against;
mannwhitneyu on an empty sample returns NaN silently (with only a scipy
warning) instead of raising, so the empty-comparator case must be excluded
explicitly -- the existing min_lines floor only guards the lineage's own
group size, not the "rest of the panel" side of the test.
"""
import csv
import os
import tempfile

import kernel


def _write(effect_values, lineages):
    d = tempfile.mkdtemp()
    with open(os.path.join(d, "CRISPRGeneEffect.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["", "GENE (1)"])
        for i, v in enumerate(effect_values):
            w.writerow(["ACH-%d" % i, v])
    with open(os.path.join(d, "Model.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ModelID", "OncotreeLineage"])
        for i, lin in enumerate(lineages):
            w.writerow(["ACH-%d" % i, lin])
    return d


def test_single_lineage_covering_every_line_is_excluded_not_nan():
    # All 20 lines are "Lung": the comparator ("rest") would be empty.
    d = _write([-0.5 - 0.01 * i for i in range(20)], ["Lung"] * 20)
    out = kernel.depmap_lineage_enrichment("GENE", root=d, min_lines=15)
    assert len(out) == 0


def test_two_lineages_both_meeting_floor_produces_real_pvalues():
    import math
    vals = [-0.9] * 15 + [-0.1] * 15
    lins = ["Lung"] * 15 + ["Breast"] * 15
    d = _write(vals, lins)
    out = kernel.depmap_lineage_enrichment("GENE", root=d, min_lines=15)
    assert set(out["lineage"]) == {"Lung", "Breast"}
    assert not out["p_enriched"].isna().any()
    assert not out["p_bh"].isna().any()
