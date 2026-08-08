from kernel import normalize_cell_name, reconcile_fetch, ot_essentiality_frame


def test_normalize_cell_name_strips_punctuation_and_uppercases():
    assert normalize_cell_name("NCI-H23") == "NCIH23"
    assert normalize_cell_name("nci-h23") == "NCIH23"
    assert normalize_cell_name("HCT 116") == "HCT116"


def test_normalize_cell_name_idempotent():
    assert normalize_cell_name(normalize_cell_name("NCI-H23")) == "NCIH23"


def test_reconcile_fetch_reports_missing_ids_sorted():
    missing = reconcile_fetch(["3", "1", "2"], ["1", "3"])
    assert missing == ["2"]


def test_reconcile_fetch_empty_when_nothing_missing():
    missing = reconcile_fetch(["1", "2"], ["2", "1"])
    assert missing == []


def test_reconcile_fetch_ignores_extra_fetched_ids():
    # A fetch returning MORE than was searched for is not a reconciliation
    # failure -- only ids present in expected but absent from fetched count.
    missing = reconcile_fetch(["1"], ["1", "2", "3"])
    assert missing == []


def test_ot_essentiality_frame_flattens_multi_tissue_multi_screen():
    target_dict = {
        "approvedSymbol": "DCTPP1",
        "depMapEssentiality": [
            {"tissueName": "Large Intestine",
             "screens": [
                 {"cellLineName": "HCT116", "diseaseFromSource": "Colorectal Cancer",
                  "geneEffect": -0.1, "expression": 3.2},
                 {"cellLineName": "SW48", "diseaseFromSource": "Colorectal Cancer",
                  "geneEffect": -0.3, "expression": 2.1},
             ]},
            {"tissueName": "Lung",
             "screens": [
                 {"cellLineName": "A549", "diseaseFromSource": "Lung Cancer",
                  "geneEffect": 0.0, "expression": 4.0},
             ]},
        ],
    }
    df = ot_essentiality_frame(target_dict)

    assert list(df.columns) == ["gene", "tissue", "cell", "disease", "effect", "expr"]
    assert len(df) == 3
    assert set(df.gene) == {"DCTPP1"}
    hct = df[df.cell == "HCT116"].iloc[0]
    assert hct.tissue == "Large Intestine"
    assert hct.disease == "Colorectal Cancer"
    assert hct.effect == -0.1
    assert hct.expr == 3.2


def test_ot_essentiality_frame_empty_when_no_data():
    df = ot_essentiality_frame({"approvedSymbol": "X", "depMapEssentiality": None})
    assert len(df) == 0
