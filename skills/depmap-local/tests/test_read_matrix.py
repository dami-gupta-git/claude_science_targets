"""depmap_read_matrix's CSV backend: column-position independence and the
collapse_to_model (IsDefaultEntryForModel) filter.

Real DepMap files disagree on where ModelID sits: CRISPRGeneEffect.csv has an
unnamed first column acting as the id; OmicsSomaticMutationsMatrixHotspot.csv
carries an explicit "ModelID" column that is NOT first. The idcol-detection
logic must handle both without hardcoding a position.
"""
import csv
import os

import kernel


def _write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def test_effect_matrix_unnamed_first_column_is_the_id(tmp_path):
    # CRISPRGeneEffect.csv style: blank header for the id column, at position 0.
    _write_csv(tmp_path / "CRISPRGeneEffect.csv",
              ["", "KRAS (3845)", "TP53 (7157)"],
              [["ACH-1", "-0.7", "0.1"],
               ["ACH-2", "-0.2", "0.3"]])

    df, missing = kernel.depmap_read_matrix("effect", ["KRAS"], root=str(tmp_path))

    assert missing == []
    assert df.index.name == "ModelID"
    assert list(df.index) == ["ACH-1", "ACH-2"]
    assert df.loc["ACH-1", "KRAS"] == -0.7
    assert list(df.columns) == ["KRAS"]  # Entrez suffix stripped


def test_gene_effect_rejects_unsupported_dataset(tmp_path):
    # Regression: dataset values other than 'effect'/'dependency' used to fall
    # through to 'dependency' silently ("kind = 'effect' if dataset == 'effect'
    # else 'dependency'"), so depmap_gene_effect(genes, dataset='hotspot')
    # returned the dependency-probability matrix instead of hotspot mutation
    # counts, or raising. hotspot/damaging are mutation matrices, not
    # gene-effect data, and must be read via depmap_read_matrix directly.
    _write_csv(tmp_path / "CRISPRGeneEffect.csv",
              ["", "KRAS (3845)"], [["ACH-1", "-0.7"]])
    _write_csv(tmp_path / "CRISPRGeneDependency.csv",
              ["", "KRAS (3845)"], [["ACH-1", "0.9"]])
    import pytest
    with pytest.raises(ValueError, match="'hotspot' and 'damaging'"):
        kernel.depmap_gene_effect(["KRAS"], root=str(tmp_path), dataset="hotspot")
    with pytest.raises(ValueError, match="dataset must be"):
        kernel.depmap_gene_effect(["KRAS"], root=str(tmp_path), dataset="damaging")
    # The two genuinely supported values still work.
    eff, _ = kernel.depmap_read_matrix("effect", ["KRAS"], root=str(tmp_path))
    assert eff.loc["ACH-1", "KRAS"] == -0.7


def test_hotspot_matrix_modelid_not_in_first_column(tmp_path):
    # OmicsSomaticMutationsMatrixHotspot.csv style: ModelID is a named column
    # that is NOT at position 0, and rows must be collapsed to one per
    # ModelID via IsDefaultEntryForModel == "Yes".
    _write_csv(tmp_path / "OmicsSomaticMutationsMatrixHotspot.csv",
              ["ProfileID", "SomeOtherCol", "ModelID", "IsDefaultEntryForModel",
               "KRAS (3845)", "TP53 (7157)"],
              [["P1", "x", "ACH-1", "Yes", "1", "0"],
               ["P1-alt", "x", "ACH-1", "No", "0", "0"],   # non-default duplicate, must be dropped
               ["P2", "x", "ACH-2", "Yes", "0", "1"]])

    df, missing = kernel.depmap_read_matrix("hotspot", ["KRAS", "TP53"], root=str(tmp_path))

    assert missing == []
    assert df.index.name == "ModelID"
    assert sorted(df.index) == ["ACH-1", "ACH-2"]
    assert len(df) == 2  # the non-default duplicate row was dropped
    assert df.loc["ACH-1", "KRAS"] == 1
    assert df.loc["ACH-2", "TP53"] == 1


def test_partial_missing_columns_reported_not_raised(tmp_path):
    _write_csv(tmp_path / "CRISPRGeneEffect.csv",
              ["", "KRAS (3845)"],
              [["ACH-1", "-0.7"]])

    df, missing = kernel.depmap_read_matrix("effect", ["KRAS", "NOTAGENE"], root=str(tmp_path))

    assert missing == ["NOTAGENE"]
    assert list(df.columns) == ["KRAS"]


def test_all_columns_missing_raises_keyerror(tmp_path):
    _write_csv(tmp_path / "CRISPRGeneEffect.csv",
              ["", "KRAS (3845)"],
              [["ACH-1", "-0.7"]])

    import pytest
    with pytest.raises(KeyError):
        kernel.depmap_read_matrix("effect", ["NOTAGENE"], root=str(tmp_path))


def test_neither_backend_available_raises_filenotfound(tmp_path):
    import pytest
    with pytest.raises(FileNotFoundError):
        kernel.depmap_read_matrix("effect", ["KRAS"], root=str(tmp_path))


def test_requested_columns_deduped_keep_order(tmp_path):
    _write_csv(tmp_path / "CRISPRGeneEffect.csv",
              ["", "KRAS (3845)", "TP53 (7157)"],
              [["ACH-1", "-0.7", "0.1"]])

    df, missing = kernel.depmap_read_matrix("effect", ["TP53", "KRAS", "TP53"], root=str(tmp_path))

    assert list(df.columns) == ["TP53", "KRAS"]
