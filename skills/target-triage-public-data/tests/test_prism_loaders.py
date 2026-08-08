import csv

from kernel import load_prism, load_prism_secondary, prism_secondary_has


def _write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        for r in rows:
            w.writerow(r)


def test_load_prism_joins_brd_ids_and_averages_replicate_screens(tmp_path):
    # BRD-A-1 and BRD-A-2 are two screens of the SAME compound (AZ-628-like
    # one-to-many quirk documented in depmap-local's SKILL.md); they must be
    # averaged into one row per drug, not kept as duplicate/conflicting rows.
    compound_list = tmp_path / "compounds.csv"
    _write_csv(compound_list, [
        ["IDs", "Drug.Name"],
        ["BRD-A-1", "AZ-628"],
        ["BRD-A-2", "AZ-628"],
        ["BRD-B-1", "OTHER-DRUG"],
    ])
    matrix = tmp_path / "matrix.csv"
    _write_csv(matrix, [
        ["", "ACH-1", "ACH-2"],
        ["BRD-A-1", "0.1", "0.3"],
        ["BRD-A-2", "0.3", "0.5"],
        ["BRD-B-1", "1.0", "2.0"],
        ["BRD-C-UNKNOWN", "9.0", "9.0"],  # not in compound list -> must be dropped
    ])

    out = load_prism(str(matrix), str(compound_list))

    assert set(out.columns) == {"AZ-628", "OTHER-DRUG"}
    assert list(out.index) == ["ACH-1", "ACH-2"]
    assert out.loc["ACH-1", "AZ-628"] == 0.2   # mean(0.1, 0.3)
    assert out.loc["ACH-2", "AZ-628"] == 0.4   # mean(0.3, 0.5)
    assert out.loc["ACH-1", "OTHER-DRUG"] == 1.0


def test_load_prism_filters_to_requested_drug_names(tmp_path):
    compound_list = tmp_path / "compounds.csv"
    _write_csv(compound_list, [
        ["IDs", "Drug.Name"],
        ["BRD-A-1", "TRIFLURIDINE"],
        ["BRD-B-1", "DECITABINE"],
        ["BRD-C-1", "SOME-OTHER-DRUG"],
    ])
    matrix = tmp_path / "matrix.csv"
    _write_csv(matrix, [
        ["", "ACH-1"],
        ["BRD-A-1", "0.1"],
        ["BRD-B-1", "0.2"],
        ["BRD-C-1", "0.3"],
    ])

    out = load_prism(str(matrix), str(compound_list),
                     drug_names=["trifluridine", "decitabine"])

    assert set(out.columns) == {"TRIFLURIDINE", "DECITABINE"}


def test_load_prism_secondary_orientation_and_replicate_averaging(tmp_path):
    path = tmp_path / "secondary.csv"
    _write_csv(path, [
        ["depmap_id", "name", "auc", "r2", "ic50"],
        ["ACH-1", "trifluridine", "0.4", "0.9", "1.0"],
        ["ACH-1", "trifluridine", "0.6", "0.95", "1.1"],  # replicate -> averaged
        ["ACH-2", "trifluridine", "0.8", "0.99", "1.2"],
        ["ACH-1", "decitabine", "0.2", "0.8", "0.5"],
    ])

    out = load_prism_secondary(str(path))

    assert out.index.name == "ModelID"
    assert "trifluridine" in out.columns
    assert abs(out.loc["ACH-1", "trifluridine"] - 0.5) < 1e-9  # mean(0.4, 0.6)
    assert out.loc["ACH-2", "trifluridine"] == 0.8


def test_load_prism_secondary_min_r2_filter(tmp_path):
    path = tmp_path / "secondary.csv"
    _write_csv(path, [
        ["depmap_id", "name", "auc", "r2"],
        ["ACH-1", "drugx", "0.4", "0.5"],   # below threshold -> excluded
        ["ACH-2", "drugx", "0.6", "0.95"],
    ])

    out = load_prism_secondary(str(path), min_r2=0.9)

    assert "ACH-1" not in out.index
    assert out.loc["ACH-2", "drugx"] == 0.6


def test_load_prism_secondary_drug_names_filter_is_case_insensitive(tmp_path):
    path = tmp_path / "secondary.csv"
    _write_csv(path, [
        ["depmap_id", "name", "auc", "r2"],
        ["ACH-1", "Trifluridine", "0.4", "0.9"],
        ["ACH-1", "OtherDrug", "0.7", "0.9"],
    ])

    out = load_prism_secondary(str(path), drug_names=["TRIFLURIDINE"])

    assert list(out.columns) == ["trifluridine"]


def test_prism_secondary_has_reports_membership(tmp_path):
    path = tmp_path / "secondary.csv"
    _write_csv(path, [
        ["depmap_id", "name", "auc", "r2"],
        ["ACH-1", "trifluridine", "0.4", "0.9"],
    ])

    out = prism_secondary_has(str(path), ["trifluridine", "decitabine"])

    assert out == {"trifluridine": True, "decitabine": False}
