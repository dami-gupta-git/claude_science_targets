"""depmap_prism_compounds's target match against the data matrix.

A compound can be annotated for `target` in the compound list but absent from
the data matrix (the converse of the documented "matrix row with no
compound-list annotation" quirk). When ALL matched compounds are of this kind,
`rows` is empty and `pd.DataFrame([]).sort_values("median_lfc")` raised
KeyError -- there is no such column on an empty frame.
"""
import csv
import os
import tempfile

import kernel


def _setup(list_rows, matrix_header, matrix_rows, release="24Q2"):
    d = tempfile.mkdtemp()
    with open(os.path.join(
            d, "Repurposing_Public_%s_Extended_Primary_Compound_List.csv" % release),
            "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["IDs", "Drug.Name", "MOA", "dose", "screen", "repurposing_target"])
        for r in list_rows:
            w.writerow(r)
    with open(os.path.join(
            d, "Repurposing_Public_%s_Extended_Primary_Data_Matrix.csv" % release),
            "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(matrix_header)
        for r in matrix_rows:
            w.writerow(r)
    return d


def test_no_annotated_compound_present_in_matrix_returns_empty_not_keyerror():
    d = _setup(
        list_rows=[["BRD:1", "DrugA", "inhibitor", "1", "REP.1M", "USP1"]],
        matrix_header=["", "ACH-1"],
        matrix_rows=[["BRD:99", "-1.5"]],  # BRD:1 never appears in the matrix
    )
    out = kernel.depmap_prism_compounds("USP1", root=d)
    assert len(out) == 0
    assert list(out.columns) == []  # untouched pd.DataFrame(), not a KeyError


def test_partial_matrix_match_returns_only_matched_rows():
    d = _setup(
        list_rows=[["BRD:1", "DrugA", "inhibitor", "1", "REP.1M", "USP1"],
                   ["BRD:2", "DrugB", "inhibitor", "1", "REP.1M", "USP1"]],
        matrix_header=["", "ACH-1", "ACH-2"],
        matrix_rows=[["BRD:1", "-1.5", "-0.2"]],  # BRD:2 absent
    )
    out = kernel.depmap_prism_compounds("USP1", root=d)
    assert len(out) == 1
    assert out.iloc[0]["drug"] == "DrugA"
