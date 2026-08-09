"""Unit tests for the boltz-affinity-triage kernel helpers.

Run from the skill directory:  pytest tests/ -q

The suite is written so that it does not need a GPU, a network connection, or a
real Boltz install: co-folding output is synthesised from the documented file
layout, and knowledge-base rows are literal dicts shaped like ChEMBL and
BindingDB responses. Tests that need RDKit skip cleanly without it, and the
regex fallback path is exercised either way.
"""
import csv
import importlib.util
import json
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
KERNEL = os.path.join(os.path.dirname(HERE), "kernel.py")


def _load():
    spec = importlib.util.spec_from_file_location("kba_kernel", KERNEL)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


k = _load()

try:
    import rdkit  # noqa: F401
    HAVE_RDKIT = True
except ImportError:
    HAVE_RDKIT = False

needs_rdkit = pytest.mark.skipif(not HAVE_RDKIT, reason="RDKit not installed")


# --------------------------------------------------------------------------
# fixtures: synthetic Boltz-2 prediction directories
# --------------------------------------------------------------------------

# model_0 is Boltz's best-ranked sample; the rest score materially worse, and
# straddle the ipTM 0.5 pose threshold so a mis-selected sample changes the
# verdict rather than just the digits.
RANKED = {
    0: {"iptm": 0.81, "ptm": 0.88, "complex_plddt": 0.86, "confidence_score": 0.84},
    1: {"iptm": 0.44, "ptm": 0.71, "complex_plddt": 0.70, "confidence_score": 0.62},
    2: {"iptm": 0.39, "ptm": 0.66, "complex_plddt": 0.65, "confidence_score": 0.55},
    3: {"iptm": 0.31, "ptm": 0.60, "complex_plddt": 0.61, "confidence_score": 0.48},
    4: {"iptm": 0.28, "ptm": 0.58, "complex_plddt": 0.59, "confidence_score": 0.44},
}


def write_pred(root, confidences, affinity=None, name="x", write_order=None,
               extra_files=True):
    """Materialise a predictions/ dir in the layout `boltz predict` produces."""
    pred = os.path.join(str(root), "boltz_results", "predictions", name)
    os.makedirs(pred, exist_ok=True)
    order = write_order if write_order is not None else sorted(confidences)
    for idx in order:
        with open(os.path.join(pred, "confidence_%s_model_%d.json" % (name, idx)), "w") as fh:
            json.dump(confidences[idx], fh)
        if extra_files:
            # pae/plddt siblings exist in real output and must not be mistaken
            # for confidence files
            with open(os.path.join(pred, "pae_%s_model_%d.json" % (name, idx)), "w") as fh:
                json.dump({"pae": [[0.1, 0.2]]}, fh)
    if affinity is not None:
        with open(os.path.join(pred, "affinity_%s.json" % name), "w") as fh:
            json.dump(affinity, fh)
    return pred


AFF_V2 = {"affinity_pred_value": -1.2, "affinity_probability_binary": 0.77}


# --------------------------------------------------------------------------
# kba_read_result -- ranked-sample selection
# --------------------------------------------------------------------------

class TestReadResultSampleSelection:
    """Confidence is per diffusion sample; affinity is per input.

    Regression: reading confidence from whichever file the filesystem listed
    first attributed an unranked sample's ipTM to the prediction. On the fixture
    above that reports 0.44 instead of 0.81 -- across the 0.5 threshold the
    workflow uses to flag unreliable poses.
    """

    def test_defaults_to_best_ranked_sample(self, tmp_path):
        pred = write_pred(tmp_path, RANKED, AFF_V2)
        got = k.kba_read_result(pred)
        assert got["model_index"] == 0
        assert got["iptm"] == RANKED[0]["iptm"]
        assert got["n_samples"] == 5

    @pytest.mark.parametrize("order", [
        [0, 1, 2, 3, 4],
        [4, 3, 2, 1, 0],
        [2, 0, 4, 1, 3],
    ])
    def test_result_is_independent_of_write_order(self, tmp_path, order):
        pred = write_pred(tmp_path, RANKED, AFF_V2, write_order=order)
        got = k.kba_read_result(pred)
        assert got["iptm"] == RANKED[0]["iptm"], (
            "confidence must come from model_0, not from directory order")

    def test_listdir_order_cannot_change_the_answer(self, tmp_path, monkeypatch):
        """Force a hostile traversal order to prove selection is by rank."""
        pred = write_pred(tmp_path, RANKED, AFF_V2)
        real = os.listdir
        monkeypatch.setattr(os, "listdir",
                            lambda p: list(reversed(real(p))))
        assert k.kba_read_result(pred)["iptm"] == RANKED[0]["iptm"]

    def test_explicit_sample_selection(self, tmp_path):
        pred = write_pred(tmp_path, RANKED, AFF_V2)
        got = k.kba_read_result(pred, model_index=3)
        assert got["model_index"] == 3
        assert got["iptm"] == RANKED[3]["iptm"]

    def test_missing_requested_sample_falls_back_to_lowest(self, tmp_path):
        pred = write_pred(tmp_path, RANKED, AFF_V2)
        assert k.kba_read_result(pred, model_index=99)["model_index"] == 0

    def test_non_contiguous_indices_pick_lowest_present(self, tmp_path):
        """A failed sample leaves gaps; selection must not assume index 0 exists."""
        pred = write_pred(tmp_path, {2: {"iptm": 0.50}, 7: {"iptm": 0.90}}, AFF_V2)
        got = k.kba_read_result(pred)
        assert (got["model_index"], got["iptm"]) == (2, 0.50)

    def test_per_sample_exposes_every_sample(self, tmp_path):
        pred = write_pred(tmp_path, RANKED, AFF_V2)
        per = k.kba_read_result(pred)["per_sample"]
        assert sorted(per) == [0, 1, 2, 3, 4]
        assert {i: per[i]["iptm"] for i in per} == {i: RANKED[i]["iptm"] for i in RANKED}

    def test_sibling_json_not_treated_as_confidence(self, tmp_path):
        """pae_*/plddt_* files share the directory and must not be counted."""
        pred = write_pred(tmp_path, RANKED, AFF_V2)
        assert k.kba_read_result(pred)["n_samples"] == 5


class TestReadResultScoreKeys:
    """Score-key names differ across Boltz versions and must be discovered."""

    def test_boltz_2_probability_preferred(self, tmp_path):
        pred = write_pred(tmp_path, {0: {"iptm": 0.7}}, AFF_V2)
        got = k.kba_read_result(pred)
        assert got["score_key"] == "affinity_probability_binary"
        assert got["binder_score"] == 0.77
        assert got["affinity_pred_value"] == -1.2

    def test_boltz_2_1_keys(self, tmp_path):
        pred = write_pred(tmp_path, {0: {"iptm": 0.7}},
                          {"binding_confidence": 0.66, "optimization_score": 2.1})
        got = k.kba_read_result(pred)
        assert got["score_key"] == "binding_confidence"
        assert got["binder_score"] == 0.66

    def test_log_ic50_only_is_sign_flipped(self, tmp_path):
        """affinity_pred_value is log10(IC50 uM): lower binds tighter, so the
        normalised binder_score must invert it to keep higher = better."""
        pred = write_pred(tmp_path, {0: {"iptm": 0.7}}, {"affinity_pred_value": -2.0})
        got = k.kba_read_result(pred)
        assert got["binder_score"] == 2.0
        assert got["score_key"] == "neg_affinity_pred_value"

    def test_binding_confidence_preferred_over_optimization_score(self, tmp_path):
        """binding_confidence is the binder-vs-decoy analogue and ranks directly."""
        pred = write_pred(tmp_path, {0: {"iptm": 0.7}},
                          {"binding_confidence": 0.66, "optimization_score": 2.1})
        got = k.kba_read_result(pred)
        assert got["score_key"] == "binding_confidence"
        assert got["optimization_score"] == 2.1

    def test_optimization_score_alone_is_not_ranked_by_default(self, tmp_path):
        """Its sign convention is undocumented, and the published descriptions of
        the analogous Boltz-2 field disagree (log10(IC50 uM) vs pIC50 -- negatives
        of each other). Guessing would order a screen backwards while still
        producing a plausible table, so the value is reported but not ranked."""
        pred = write_pred(tmp_path, {0: {"iptm": 0.7}}, {"optimization_score": 2.1})
        got = k.kba_read_result(pred)
        assert got["binder_score"] is None
        assert got["optimization_score"] == 2.1
        assert "sign convention" in got["unscored_reason"]

    @pytest.mark.parametrize("direction,expected,key", [
        ("higher_is_better", 2.1, "optimization_score"),
        ("lower_is_better", -2.1, "neg_optimization_score"),
    ])
    def test_optimization_score_ranked_when_direction_declared(self, tmp_path,
                                                              direction, expected, key):
        pred = write_pred(tmp_path, {0: {"iptm": 0.7}}, {"optimization_score": 2.1})
        got = k.kba_read_result(pred, optimization_score_direction=direction)
        assert got["binder_score"] == expected
        assert got["score_key"] == key
        assert got["unscored_reason"] is None

    def test_unrecognised_affinity_keys_are_reported(self, tmp_path):
        """A future rename must surface as a named reason, not a silent None."""
        pred = write_pred(tmp_path, {0: {"iptm": 0.7}}, {"some_future_key": 1.0})
        got = k.kba_read_result(pred)
        assert got["binder_score"] is None
        assert "some_future_key" in got["unscored_reason"]

    def test_no_affinity_file_has_no_unscored_reason(self, tmp_path):
        """A structure-only run is not a scoring failure."""
        pred = write_pred(tmp_path, {0: {"iptm": 0.7}}, affinity=None)
        got = k.kba_read_result(pred)
        assert got["binder_score"] is None and got["unscored_reason"] is None

    def test_tighter_binder_scores_higher(self, tmp_path):
        a = write_pred(tmp_path / "a", {0: {"iptm": 0.7}}, {"affinity_pred_value": -3.0})
        b = write_pred(tmp_path / "b", {0: {"iptm": 0.7}}, {"affinity_pred_value": 0.0})
        assert k.kba_read_result(a)["binder_score"] > k.kba_read_result(b)["binder_score"]


class TestReadResultDegradedInputs:

    def test_structure_only_run_has_no_binder_score(self, tmp_path):
        pred = write_pred(tmp_path, {0: {"iptm": 0.7}}, affinity=None)
        got = k.kba_read_result(pred)
        assert got["binder_score"] is None and got["score_key"] is None
        assert got["iptm"] == 0.7

    def test_single_chain_has_ptm_but_no_iptm(self, tmp_path):
        """ipTM is interface-only; a monomer run legitimately lacks it."""
        pred = write_pred(tmp_path, {0: {"ptm": 0.91, "complex_plddt": 0.88}}, AFF_V2)
        got = k.kba_read_result(pred)
        assert got["ptm"] == 0.91 and got["iptm"] is None

    def test_missing_directory_returns_empty_not_raises(self, tmp_path):
        got = k.kba_read_result(str(tmp_path / "absent"))
        assert got["binder_score"] is None and got["n_samples"] == 0

    def test_corrupt_json_is_skipped(self, tmp_path):
        pred = write_pred(tmp_path, {0: {"iptm": 0.7}}, AFF_V2)
        with open(os.path.join(pred, "confidence_x_model_1.json"), "w") as fh:
            fh.write("{not json")
        got = k.kba_read_result(pred)
        assert got["n_samples"] == 1 and got["iptm"] == 0.7


# --------------------------------------------------------------------------
# kba_parse_nm / kba_pactivity -- KB value parsing
# --------------------------------------------------------------------------

class TestParseNm:
    """BindingDB affinity arrives as a string that may carry a qualifier."""

    @pytest.mark.parametrize("raw,value,qual", [
        ("41.0", 41.0, "="),
        (">10000", 10000.0, ">"),
        ("<1", 1.0, "<"),
        ("  >  500 ", 500.0, ">"),
        ("1,200", 1200.0, "="),
    ])
    def test_parses_value_and_qualifier(self, raw, value, qual):
        assert k.kba_parse_nm(raw) == (value, qual)

    @pytest.mark.parametrize("raw", [None, "", "   ", "n/a", ">", "abc"])
    def test_unparseable_returns_none_pair(self, raw):
        assert k.kba_parse_nm(raw) == (None, None)

    def test_qualifier_is_preserved_for_censored_values(self):
        """A '>' is a real inactive and a '<' a ceiling hit; callers need the
        distinction, so it must not be collapsed into a bare float."""
        assert k.kba_parse_nm(">10000")[1] == ">"
        assert k.kba_parse_nm("<1")[1] == "<"


class TestPactivity:

    def test_prefers_pchembl_when_present(self):
        row = {"pchembl_value": 7.39, "standard_value": 41.0, "standard_units": "nM"}
        pact, qual, src = k.kba_pactivity(row)
        assert (round(pact, 2), src) == (7.39, "pchembl_value")

    def test_falls_back_to_standard_value(self):
        """ChEMBL omits pchembl_value on most weak/qualified rows; without this
        fallback the negative class is silently lost."""
        row = {"pchembl_value": None, "standard_value": 500000.0,
               "standard_units": "nM", "standard_relation": "="}
        pact, qual, src = k.kba_pactivity(row)
        assert src == "standard_value"
        assert round(pact, 3) == 3.301  # 500 uM

    @pytest.mark.parametrize("units,value,expected", [
        ("nM", 1.0, 9.0),
        ("uM", 1.0, 6.0),
        ("mM", 1.0, 3.0),
        ("pM", 1.0, 12.0),
        ("M", 1.0, 0.0),
    ])
    def test_unit_scaling(self, units, value, expected):
        row = {"standard_value": value, "standard_units": units}
        assert round(k.kba_pactivity(row)[0], 6) == expected

    def test_bindingdb_source(self):
        pact, qual, _ = k.kba_pactivity({"affinity": "1"}, source="bindingdb")
        assert (pact, qual) == (9.0, "=")

    def test_bindingdb_qualifier_surfaces(self):
        assert k.kba_pactivity({"affinity": "<1"}, source="bindingdb")[1] == "<"

    @pytest.mark.parametrize("row", [
        {"standard_value": 0, "standard_units": "nM"},
        {"standard_value": -5, "standard_units": "nM"},
        {"standard_value": 10, "standard_units": "ug.mL-1"},
        {"standard_value": None, "standard_units": "nM"},
        {},
    ])
    def test_unusable_rows_return_none(self, row):
        assert k.kba_pactivity(row)[0] is None

    def test_potent_scores_above_weak(self):
        potent = k.kba_pactivity({"standard_value": 41.0, "standard_units": "nM"})[0]
        weak = k.kba_pactivity({"standard_value": 500000.0, "standard_units": "nM"})[0]
        assert potent > weak


# --------------------------------------------------------------------------
# kba_leakage_flag
# --------------------------------------------------------------------------

class TestLeakageFlag:
    """The affinity head trained on public bioactivity data as of mid-2023, so
    publication year is the available proxy for whether a control is recalled."""

    @pytest.mark.parametrize("year,flag", [
        (1997, "likely_train"),
        (2003, "likely_train"),
        (2023, "likely_train"),
        (2024, "likely_novel"),
        (2026, "likely_novel"),
    ])
    def test_year_split(self, year, flag):
        assert k.kba_leakage_flag({"document_year": year}) == flag

    @pytest.mark.parametrize("row", [{}, {"document_year": None},
                                     {"document_year": "n/a"}])
    def test_unknown_year(self, row):
        assert k.kba_leakage_flag(row) == "unknown"

    def test_string_year_is_coerced(self):
        assert k.kba_leakage_flag({"document_year": "2025"}) == "likely_novel"

    def test_cutoff_is_configurable(self):
        assert k.kba_leakage_flag({"document_year": 2020}, cutoff_year=2019) == "likely_novel"


# --------------------------------------------------------------------------
# heavy atoms / properties
# --------------------------------------------------------------------------

ATOM_TRUTH = [
    ("CC(=O)OC1=CC=CC=C1C(=O)O", 13),   # aspirin
    ("c1ccccc1", 6),                     # benzene, aromatic lowercase
    ("CN1C=NC2=C1C(=O)N(C)C(=O)N2C", 14),  # caffeine
    ("CCO", 3),
    ("ClCCBr", 4),                       # two-letter elements
]


class TestHeavyAtoms:

    @needs_rdkit
    @pytest.mark.parametrize("smiles,n", ATOM_TRUTH)
    def test_exact_with_rdkit(self, smiles, n):
        assert k.kba_heavy_atoms(smiles) == n

    @pytest.mark.parametrize("smiles,n", ATOM_TRUTH)
    def test_fallback_within_tolerance(self, smiles, n, monkeypatch):
        """Without RDKit a regex estimate is used; it is only required to be
        close enough for property matching, not exact."""
        monkeypatch.setitem(sys.modules, "rdkit", None)
        assert abs(k.kba_heavy_atoms(smiles) - n) <= 2

    def test_aromatic_lowercase_counted(self, monkeypatch):
        """Regression: an earlier fallback stripped lowercase aromatic atoms and
        returned 2 for benzene."""
        monkeypatch.setitem(sys.modules, "rdkit", None)
        assert k.kba_heavy_atoms("c1ccccc1") == 6

    def test_bracketed_atoms_counted_once(self, monkeypatch):
        monkeypatch.setitem(sys.modules, "rdkit", None)
        assert k.kba_heavy_atoms("[Na+]") == 1

    def test_empty_smiles(self):
        assert k.kba_heavy_atoms("") == 0


class TestProps:

    def test_always_returns_heavy_atoms(self):
        p = k.kba_props("CCO")
        assert p["heavy_atoms"] == 3 and p["smiles"] == "CCO"

    @needs_rdkit
    def test_rdkit_fields_present(self):
        p = k.kba_props("CC(=O)OC1=CC=CC=C1C(=O)O")
        assert 179 < p["mw"] < 181 and p["rings"] == 1
        assert p["clogp"] is not None

    def test_invalid_smiles_does_not_raise(self):
        p = k.kba_props("not-a-molecule")
        assert p["mw"] is None


# --------------------------------------------------------------------------
# decoy matching and the size-only null
# --------------------------------------------------------------------------

BIG = ["CCCCCCCCCCCCCCCCCCCC", "CCCCCCCCCCCCCCCCCCCCC",
       "CCCCCCCCCCCCCCCCCCCCCC", "CCCCCCCCCCCCCCCCCCC"]
SMALL = ["CCO", "CCC", "CCN", "CCCC"]


class TestMatchDecoys:

    def test_prefers_size_matched_candidates(self):
        actives = [{"smiles": s} for s in BIG]
        pool = [{"smiles": s} for s in SMALL + BIG]
        chosen = k.kba_match_decoys(actives, pool, per_active=1)
        assert chosen, "expected at least one match"
        for c in chosen:
            assert k.kba_heavy_atoms(c["smiles"]) >= 15

    def test_sampling_is_without_replacement(self):
        actives = [{"smiles": s} for s in BIG]
        pool = [{"smiles": s} for s in BIG]
        chosen = k.kba_match_decoys(actives, pool, per_active=1)
        assert len(chosen) == len({id(c) for c in chosen})

    def test_no_candidate_within_tolerance_yields_nothing(self):
        actives = [{"smiles": s} for s in BIG]
        pool = [{"smiles": s} for s in SMALL]
        assert k.kba_match_decoys(actives, pool, per_active=1,
                                  tol={"heavy_atoms": 1, "mw": 5.0}) == []

    def test_empty_inputs(self):
        assert k.kba_match_decoys([], [{"smiles": "CCO"}]) == []
        assert k.kba_match_decoys([{"smiles": "CCO"}], []) == []

    def test_rows_without_smiles_are_ignored(self):
        chosen = k.kba_match_decoys([{"smiles": "CCO"}],
                                    [{"name": "no smiles"}, {"smiles": "CCC"}])
        assert len(chosen) == 1 and chosen[0]["smiles"] == "CCC"

    def test_partial_tol_override_does_not_raise(self):
        """SKILL.md tells the reader to 'loosen tol' when the null won't balance;
        a caller supplying only one key must not hit a KeyError on the other."""
        actives = [{"smiles": s} for s in BIG]
        pool = [{"smiles": s} for s in SMALL + BIG]
        chosen = k.kba_match_decoys(actives, pool, per_active=1,
                                    tol={"heavy_atoms": 10})
        assert chosen
        chosen2 = k.kba_match_decoys(actives, pool, per_active=1,
                                     tol={"mw": 50.0})
        assert chosen2


class TestNullAuc:
    """The size-only null is the confounder check: it reports what a scorer
    reading nothing but molecular size achieves on the same controls."""

    def test_size_separated_controls_give_high_null(self):
        actives = [{"smiles": s} for s in BIG]
        decoys = [{"smiles": s} for s in SMALL]
        assert k.kba_null_auc(actives, decoys) > 0.9

    def test_size_balanced_controls_give_chance_null(self):
        actives = [{"smiles": s} for s in BIG]
        decoys = [{"smiles": s} for s in BIG]
        assert k.kba_null_auc(actives, decoys) == pytest.approx(0.5, abs=0.05)

    def test_matching_reduces_the_null(self):
        """End-to-end property under test: matching against a size-diverse pool
        moves the size-only null toward chance."""
        actives = [{"smiles": s} for s in BIG]
        pool = [{"smiles": s} for s in SMALL + BIG]
        before = k.kba_null_auc(actives, [{"smiles": s} for s in SMALL])
        after = k.kba_null_auc(actives, k.kba_match_decoys(actives, pool, per_active=1))
        assert after < before


# --------------------------------------------------------------------------
# ranking statistics
# --------------------------------------------------------------------------

class TestAuc:

    def test_perfect_separation(self):
        assert k.kba_auc([3, 4, 5], [0, 1, 2]) == 1.0

    def test_inverted_separation(self):
        assert k.kba_auc([0, 1, 2], [3, 4, 5]) == 0.0

    def test_all_ties_is_chance(self):
        assert k.kba_auc([1, 1, 1], [1, 1, 1]) == 0.5

    def test_partial_ties_handled(self):
        assert 0.0 < k.kba_auc([2, 3], [1, 2]) < 1.0

    @pytest.mark.parametrize("pos,neg", [([], [1]), ([1], []), ([], []),
                                         ([None], [1]), ([1], [None])])
    def test_empty_class_returns_none(self, pos, neg):
        assert k.kba_auc(pos, neg) is None

    def test_none_values_are_dropped(self):
        assert k.kba_auc([3, None, 4], [0, 1]) == 1.0

    def test_higher_score_means_binder(self):
        """Directional contract the gate and ranking both depend on."""
        assert k.kba_auc([10, 9], [1, 2]) > 0.5


class TestEnrichment:

    def test_all_actives_on_top(self):
        scored = [(9, 1), (8, 1), (2, 0), (1, 0)]
        assert k.kba_enrichment(scored, frac=0.5) == 2.0

    def test_no_enrichment_when_interleaved(self):
        scored = [(9, 1), (8, 0), (7, 1), (6, 0)]
        assert k.kba_enrichment(scored, frac=0.5) == 1.0

    def test_actives_at_bottom_gives_zero(self):
        scored = [(9, 0), (8, 0), (2, 1), (1, 1)]
        assert k.kba_enrichment(scored, frac=0.5) == 0.0

    def test_no_actives_returns_none(self):
        assert k.kba_enrichment([(9, 0), (1, 0)]) is None

    def test_empty_returns_none(self):
        assert k.kba_enrichment([]) is None

    def test_top_slice_is_at_least_one_item(self):
        assert k.kba_enrichment([(9, 1), (1, 0)], frac=0.01) is not None


class TestPercentile:

    def test_midpoint(self):
        assert k.kba_percentile(5, [1, 2, 3, 4, 6, 7]) == pytest.approx(66.67, abs=0.01)

    def test_below_all(self):
        assert k.kba_percentile(0, [1, 2, 3]) == 0.0

    def test_above_all(self):
        assert k.kba_percentile(9, [1, 2, 3]) == 100.0

    def test_ties_count_half(self):
        assert k.kba_percentile(2, [1, 2, 3]) == pytest.approx(50.0)

    def test_none_inputs(self):
        assert k.kba_percentile(None, [1, 2]) is None
        assert k.kba_percentile(1, []) is None


# --------------------------------------------------------------------------
# the gate
# --------------------------------------------------------------------------

class TestGate:
    """A run whose known actives do not separate from matched decoys has no
    demonstrated discrimination, so the gate withholds the ranking."""

    def test_passes_on_strong_separation(self):
        g = k.kba_gate(0.85, 20, 20, null_auc=0.50)
        assert g["pass"] and g["verdict"] == "interpretable" and g["reasons"] == []

    def test_fails_at_chance(self):
        g = k.kba_gate(0.50, 20, 20, null_auc=0.50)
        assert not g["pass"]
        assert any("below" in r for r in g["reasons"])

    def test_fails_when_null_is_not_cleared(self):
        """A high AUC that merely matches the size-only null is a size artifact."""
        g = k.kba_gate(0.95, 20, 20, null_auc=0.94)
        assert not g["pass"]
        assert any("size" in r for r in g["reasons"])

    def test_high_auc_over_low_null_passes(self):
        assert k.kba_gate(0.95, 20, 20, null_auc=0.50)["pass"]

    def test_fails_on_too_few_controls(self):
        g = k.kba_gate(0.99, 2, 2, null_auc=0.5)
        assert not g["pass"]
        assert any("too few" in r for r in g["reasons"])

    def test_reports_every_failure_not_just_the_first(self):
        g = k.kba_gate(0.50, 1, 1, null_auc=0.49)
        assert len(g["reasons"]) >= 2

    def test_none_auc_fails(self):
        g = k.kba_gate(None, 20, 20)
        assert not g["pass"] and any("not computable" in r for r in g["reasons"])

    def test_null_is_optional(self):
        assert k.kba_gate(0.85, 20, 20)["pass"]

    def test_reported_values_are_echoed(self):
        g = k.kba_gate(0.85, 11, 12, null_auc=0.51)
        assert (g["auc"], g["null_auc"], g["n_pos"], g["n_neg"]) == (0.85, 0.51, 11, 12)

    def test_thresholds_are_configurable(self):
        assert k.kba_gate(0.60, 20, 20, null_auc=0.5, min_auc=0.55)["pass"]

    @pytest.mark.parametrize("auc,null", [
        (0.60, 0.50), (0.85, 0.75), (0.95, 0.85), (0.75, 0.65), (0.70, 0.60),
    ])
    def test_exactly_satisfied_margin_passes(self, auc, null):
        """Binary floating point makes several exact margins evaluate a few ulps
        low (0.60 - 0.50 is 0.09999999999999998), so a run that exactly meets
        the documented requirement must not be rejected."""
        g = k.kba_gate(auc, 20, 20, null_auc=null, min_auc=0.0)
        assert g["pass"], g["reasons"]

    def test_auc_exactly_at_min_passes(self):
        assert k.kba_gate(0.65, 20, 20, null_auc=0.30, min_auc=0.65)["pass"]

    @pytest.mark.parametrize("auc,null", [(0.84, 0.75), (0.59, 0.50)])
    def test_margin_just_below_still_fails(self, auc, null):
        """The tolerance must not swallow a genuine shortfall."""
        assert not k.kba_gate(auc, 20, 20, null_auc=null, min_auc=0.0)["pass"]


# --------------------------------------------------------------------------
# YAML writer
# --------------------------------------------------------------------------

class TestWriteYaml:

    def test_writes_affinity_block(self, tmp_path):
        p = k.kba_write_yaml(str(tmp_path / "a.yaml"), "MVTPEG", "CCO", ligand_id="L")
        text = open(p).read()
        assert "properties:" in text and "binder: L" in text

    def test_parses_as_expected_structure(self, tmp_path):
        yaml = pytest.importorskip("yaml")
        p = k.kba_write_yaml(str(tmp_path / "a.yaml"), "MVTPEG", "CCO")
        doc = yaml.safe_load(open(p))
        assert doc["version"] == 1
        assert doc["sequences"][0]["protein"]["sequence"] == "MVTPEG"
        assert doc["sequences"][1]["ligand"]["smiles"] == "CCO"
        assert doc["properties"][0]["affinity"]["binder"] == "L"

    def test_smiles_with_hash_survives_quoting(self, tmp_path):
        """SMILES contain characters YAML treats specially; the value is quoted
        so triple bonds and brackets round-trip."""
        yaml = pytest.importorskip("yaml")
        smi = "CC#CC1=CC=CC=C1[N+](=O)[O-]"
        p = k.kba_write_yaml(str(tmp_path / "b.yaml"), "MVT", smi)
        assert yaml.safe_load(open(p))["sequences"][1]["ligand"]["smiles"] == smi

    def test_structure_only_omits_properties(self, tmp_path):
        p = k.kba_write_yaml(str(tmp_path / "c.yaml"), "MVT", "CCO",
                             request_affinity=False)
        assert "properties:" not in open(p).read()

    def test_oversized_ligand_is_skipped(self, tmp_path):
        """Boltz refuses affinity ligands at or above 128 atoms; returning None
        keeps the caller from believing an affinity was requested."""
        big = "C" * (k.AFFINITY_ATOM_CAP + 10)
        target = tmp_path / "big.yaml"
        assert k.kba_write_yaml(str(target), "MVT", big) is None
        assert not target.exists()

    def test_oversized_ligand_allowed_without_affinity(self, tmp_path):
        big = "C" * (k.AFFINITY_ATOM_CAP + 10)
        assert k.kba_write_yaml(str(tmp_path / "ok.yaml"), "MVT", big,
                                request_affinity=False) is not None

    def test_msa_line_written_when_given(self, tmp_path):
        p = k.kba_write_yaml(str(tmp_path / "m.yaml"), "MVT", "CCO", msa="x.a3m")
        assert "msa: x.a3m" in open(p).read()

    def test_msa_omitted_by_default(self, tmp_path):
        p = k.kba_write_yaml(str(tmp_path / "n.yaml"), "MVT", "CCO")
        assert "msa:" not in open(p).read()

    def test_custom_chain_ids(self, tmp_path):
        p = k.kba_write_yaml(str(tmp_path / "i.yaml"), "MVT", "CCO",
                             ligand_id="Z", protein_id="B")
        text = open(p).read()
        assert "id: B" in text and "id: Z" in text and "binder: Z" in text

    def test_creates_missing_parent_directory(self, tmp_path):
        p = k.kba_write_yaml(str(tmp_path / "deep" / "nest" / "d.yaml"), "MVT", "CCO")
        assert os.path.exists(p)


# --------------------------------------------------------------------------
# integration: writer -> reader -> stats -> gate
# --------------------------------------------------------------------------

def test_round_trip_write_read(tmp_path):
    """A YAML written by the skill names the ligand chain the reader's affinity
    block refers to, and confidence comes back from the ranked sample."""
    k.kba_write_yaml(str(tmp_path / "q.yaml"), "MVTPEG", "CCO", ligand_id="L")
    pred = write_pred(tmp_path, RANKED, AFF_V2)
    got = k.kba_read_result(pred)
    assert got["binder_score"] == 0.77 and got["iptm"] == 0.81


def test_confounded_calibration_is_rejected(tmp_path):
    """A scorer reading only molecular size must not pass the gate once the
    size-only null is taken into account."""
    actives = [{"smiles": s} for s in BIG]
    decoys = [{"smiles": s} for s in SMALL]
    size_score = lambda r: k.kba_heavy_atoms(r["smiles"]) / 50.0
    auc = k.kba_auc([size_score(a) for a in actives], [size_score(d) for d in decoys])
    null = k.kba_null_auc(actives, decoys)
    assert not k.kba_gate(auc, len(actives), len(decoys), null_auc=null,
                          min_per_class=1)["pass"]


def test_genuine_separation_passes_the_gate():
    """A scorer tracking measured potency clears a chance-level null."""
    actives = [{"standard_value": v, "standard_units": "nM"} for v in (1, 5, 10, 20, 40)]
    decoys = [{"standard_value": v, "standard_units": "nM"}
              for v in (50000, 100000, 200000, 300000, 500000)]
    pos = [k.kba_pactivity(r)[0] for r in actives]
    neg = [k.kba_pactivity(r)[0] for r in decoys]
    g = k.kba_gate(k.kba_auc(pos, neg), len(pos), len(neg), null_auc=0.50)
    assert g["pass"] and g["auc"] == 1.0


# --------------------------------------------------------------------------
# Run outputs
# --------------------------------------------------------------------------

PASSING_GATE = {"pass": True, "auc": 0.82, "null_auc": 0.51, "n_pos": 8,
                "n_neg": 8, "reasons": [], "verdict": "interpretable"}
FAILING_GATE = {"pass": False, "auc": 0.97, "null_auc": 0.96, "n_pos": 8,
                "n_neg": 8, "verdict": "not interpretable",
                "reasons": ["control AUC 0.97 does not clear size-only null 0.96"]}


def _controls(n=4):
    rows = []
    for i in range(n):
        rows.append({"compound_id": f"A{i}", "role": "active", "smiles": "c1ccccc1O",
                     "binder_score": 0.7 + 0.01 * i, "iptm": 0.8, "source": "chembl"})
        rows.append({"compound_id": f"D{i}", "role": "decoy", "smiles": "c1ccccc1",
                     "binder_score": 0.2 + 0.01 * i, "iptm": 0.7, "source": "chembl"})
    return rows


def _queries():
    return [{"compound_id": "Q1", "smiles": "c1ccccc1N", "binder_score": 0.9, "iptm": 0.75},
            {"compound_id": "Q2", "smiles": "c1ccccc1C", "binder_score": 0.1, "iptm": 0.31}]


class TestRunDir:
    def test_slugs_target_to_snake_case(self, tmp_path):
        out = k.kba_run_dir("EGFR kinase-1", root=str(tmp_path))
        assert out.endswith(os.path.join("boltz_affinity_triage", "egfr_kinase_1"))

    def test_creates_scripts_subdirectory(self, tmp_path):
        out = k.kba_run_dir("EGFR", root=str(tmp_path))
        assert os.path.isdir(os.path.join(out, "scripts"))

    def test_make_false_does_not_create(self, tmp_path):
        out = k.kba_run_dir("EGFR", root=str(tmp_path), make=False)
        assert not os.path.exists(out)

    def test_empty_target_raises(self, tmp_path):
        with pytest.raises(ValueError):
            k.kba_run_dir("   ", root=str(tmp_path))

    def test_results_root_raises_when_unconfigured(self, monkeypatch):
        monkeypatch.delenv("SCIENCE_RESULTS_ROOT", raising=False)
        with pytest.raises(FileNotFoundError, match="SCIENCE_RESULTS_ROOT"):
            k.results_root()

    def test_raises_and_creates_nothing_when_root_unconfigured(self, monkeypatch, tmp_path):
        monkeypatch.delenv("SCIENCE_RESULTS_ROOT", raising=False)
        monkeypatch.chdir(tmp_path)
        with pytest.raises(FileNotFoundError):
            k.kba_run_dir("EGFR")
        assert not (tmp_path / "results").exists()

    def test_honours_env_var(self, monkeypatch, tmp_path):
        monkeypatch.setenv("SCIENCE_RESULTS_ROOT", str(tmp_path))
        out = k.kba_run_dir("EGFR")
        assert out == str(tmp_path / "boltz_affinity_triage" / "egfr")


class TestWriteTable:
    def test_empty_rows_write_nothing(self, tmp_path):
        assert k.kba_write_table(str(tmp_path / "x.csv"), []) is None
        assert not (tmp_path / "x.csv").exists()

    def test_none_becomes_empty_cell(self, tmp_path):
        path = k.kba_write_table(str(tmp_path / "x.csv"),
                                 [{"a": 1, "b": None}], ["a", "b"])
        assert open(path).read().splitlines()[1] == "1,"

    def test_extra_keys_are_kept_not_dropped(self, tmp_path):
        path = k.kba_write_table(str(tmp_path / "x.csv"),
                                 [{"a": 1, "zz": 2}], ["a"])
        header = open(path).read().splitlines()[0]
        assert header == "a,zz"

    def test_headers_absent_from_rows_are_omitted(self, tmp_path):
        path = k.kba_write_table(str(tmp_path / "x.csv"), [{"a": 1}], ["a", "b"])
        assert open(path).read().splitlines()[0] == "a"


class TestPoseFlags:
    def test_low_iptm_is_flagged(self):
        assert any("low pose confidence" in f
                   for f in k.kba_pose_flags({"iptm": 0.3, "binder_score": 0.5}))

    def test_iptm_at_floor_is_not_flagged(self):
        assert k.kba_pose_flags({"iptm": k.IPTM_FLOOR, "binder_score": 0.5}) == []

    def test_ptm_used_when_iptm_absent(self):
        assert any("low pose confidence" in f
                   for f in k.kba_pose_flags({"ptm": 0.2, "binder_score": 0.5}))

    def test_missing_confidence_is_its_own_flag(self):
        assert k.kba_pose_flags({"binder_score": 0.5}) == ["no pose confidence"]

    def test_unscored_reason_is_carried_through(self):
        flags = k.kba_pose_flags({"iptm": 0.9, "binder_score": None,
                                  "unscored_reason": "sign undocumented"})
        assert any("sign undocumented" in f for f in flags)


class TestWriteRunGating:
    def test_failing_gate_writes_no_ranked_table(self, tmp_path):
        written = k.kba_write_run(str(tmp_path), "EGFR", FAILING_GATE, _controls(),
                                  queries=_queries(), summary="Controls do not separate.")
        assert "queries" not in written
        assert not (tmp_path / k.QUERIES_CSV).exists()

    def test_failing_gate_still_writes_controls_and_calibration(self, tmp_path):
        k.kba_write_run(str(tmp_path), "EGFR", FAILING_GATE, _controls(),
                        queries=_queries(), summary="Controls do not separate.")
        assert (tmp_path / k.CONTROLS_CSV).exists()
        assert (tmp_path / k.CALIBRATION_JSON).exists()

    def test_failing_gate_readme_names_the_reasons(self, tmp_path):
        k.kba_write_run(str(tmp_path), "EGFR", FAILING_GATE, _controls(),
                        queries=_queries(), summary="Controls do not separate.")
        text = (tmp_path / k.RUN_README).read_text()
        assert "size-only null" in text and "Not reported" in text
        assert "Q1" not in text

    def test_passing_gate_writes_ranked_table(self, tmp_path):
        written = k.kba_write_run(str(tmp_path), "EGFR", PASSING_GATE, _controls(),
                                  queries=_queries(), summary="Controls separate.")
        assert written["queries"].endswith(k.QUERIES_CSV)
        assert "Q1" in (tmp_path / k.RUN_README).read_text()

    def test_calibration_json_records_whether_ranking_was_written(self, tmp_path):
        k.kba_write_run(str(tmp_path), "EGFR", FAILING_GATE, _controls(),
                        queries=_queries(), summary="No separation.")
        blob = json.loads((tmp_path / k.CALIBRATION_JSON).read_text())
        assert blob["ranking_written"] is False
        assert blob["n_queries"] == 2

    def test_gate_must_be_a_gate_dict(self, tmp_path):
        with pytest.raises(ValueError):
            k.kba_write_run(str(tmp_path), "EGFR", {"auc": 0.9}, _controls(),
                            summary="x")

    def test_bad_control_role_raises_naming_the_compound(self, tmp_path):
        rows = _controls()
        rows[0]["role"] = "positive"
        with pytest.raises(ValueError, match="A0"):
            k.kba_write_run(str(tmp_path), "EGFR", PASSING_GATE, rows, summary="x")


class TestWriteRunEnrichment:
    def test_percentile_filled_from_control_distribution(self, tmp_path):
        k.kba_write_run(str(tmp_path), "EGFR", PASSING_GATE, _controls(),
                        queries=_queries(), summary="Controls separate.")
        rows = list(csv.DictReader(open(tmp_path / k.QUERIES_CSV)))
        by_id = {r["compound_id"]: r for r in rows}
        assert float(by_id["Q1"]["control_percentile"]) == 100.0
        assert float(by_id["Q2"]["control_percentile"]) == 0.0

    def test_caller_supplied_percentile_is_not_overwritten(self, tmp_path):
        queries = _queries()
        queries[0]["control_percentile"] = 42.0
        k.kba_write_run(str(tmp_path), "EGFR", PASSING_GATE, _controls(),
                        queries=queries, summary="Controls separate.")
        rows = {r["compound_id"]: r for r in csv.DictReader(open(tmp_path / k.QUERIES_CSV))}
        assert float(rows["Q1"]["control_percentile"]) == 42.0

    def test_low_confidence_query_is_flagged_in_the_table(self, tmp_path):
        k.kba_write_run(str(tmp_path), "EGFR", PASSING_GATE, _controls(),
                        queries=_queries(), summary="Controls separate.")
        rows = {r["compound_id"]: r for r in csv.DictReader(open(tmp_path / k.QUERIES_CSV))}
        assert "low pose confidence" in rows["Q2"]["flags"]
        assert rows["Q1"]["flags"] == ""

    def test_skipped_compounds_reach_readme_and_calibration(self, tmp_path):
        k.kba_write_run(str(tmp_path), "EGFR", PASSING_GATE, _controls(),
                        queries=_queries(), summary="Controls separate.",
                        skipped=["BIG-1"])
        assert "BIG-1" in (tmp_path / k.RUN_README).read_text()
        blob = json.loads((tmp_path / k.CALIBRATION_JSON).read_text())
        assert blob["skipped_over_atom_cap"] == ["BIG-1"]

    def test_scripts_are_copied_into_the_run(self, tmp_path):
        src = tmp_path / "run_triage.py"
        src.write_text("# analysis wiring\n")
        out = tmp_path / "run"
        k.kba_write_run(str(out), "EGFR", PASSING_GATE, _controls(),
                        summary="Controls separate.", scripts=[str(src)])
        assert (out / "scripts" / "run_triage.py").exists()


class TestRunReadme:
    def test_summary_over_cap_raises(self):
        with pytest.raises(ValueError):
            k.kba_run_readme("EGFR", PASSING_GATE,
                             "word " * (k.SUMMARY_MAX_WORDS + 1), _controls())

    def test_summary_at_cap_is_accepted(self):
        text = k.kba_run_readme("EGFR", PASSING_GATE,
                                "word " * k.SUMMARY_MAX_WORDS, _controls())
        assert "## Result" in text

    def test_bulleted_summary_raises(self):
        with pytest.raises(ValueError):
            k.kba_run_readme("EGFR", PASSING_GATE, "- a bullet", _controls())

    def test_missing_summary_raises_rather_than_rendering_none(self):
        # kba_write_run defaults summary=None; omitting it must raise here
        # rather than the README literally saying "None" for its opening
        # paragraph.
        with pytest.raises(ValueError):
            k.kba_run_readme("EGFR", PASSING_GATE, None, _controls())

    def test_sections_are_present_and_ordered(self):
        text = k.kba_run_readme("EGFR", PASSING_GATE, "Prose.", _controls(),
                                queries=_queries(), data_sources=["ChEMBL v34"],
                                limits=["One target."])
        idx = [text.index(h) for h in ("## Result", "## Files", "## Data sources", "## Limits")]
        assert idx == sorted(idx)

    def test_potency_caveat_is_always_stated(self):
        text = k.kba_run_readme("EGFR", PASSING_GATE, "Prose.", _controls())
        assert "not a potency" in text or "not quoted as a measured affinity" in text

    def test_margin_reported_when_both_aucs_present(self):
        text = k.kba_run_readme("EGFR", PASSING_GATE, "Prose.", _controls())
        assert "0.31" in text

    def test_missing_null_auc_is_named_not_guessed(self):
        gate = dict(PASSING_GATE, null_auc=None)
        text = k.kba_run_readme("EGFR", gate, "Prose.", _controls())
        assert "not computed" in text

    def test_calibration_only_run_says_so(self):
        text = k.kba_run_readme("EGFR", PASSING_GATE, "Prose.", _controls(), queries=None)
        assert "calibration only" in text

    def test_single_low_confidence_row_is_not_pluralised(self):
        queries = [{"compound_id": "Q1", "binder_score": 0.5, "iptm": 0.2}]
        text = k.kba_run_readme("EGFR", PASSING_GATE, "Prose.", _controls(), queries=queries)
        assert "1 compound scored on a pose below" in text

    def test_uncomputed_similarity_is_stated_not_implied(self):
        queries = [{"compound_id": "Q1", "binder_score": 0.5, "iptm": 0.8,
                    "nearest_tanimoto": None}]
        text = k.kba_run_readme("EGFR", PASSING_GATE, "Prose.", _controls(), queries=queries)
        assert "not computed" in text and "RDKit" in text

    def test_similarity_note_absent_when_it_was_computed(self):
        queries = [{"compound_id": "Q1", "binder_score": 0.5, "iptm": 0.8,
                    "nearest_active": "A0", "nearest_tanimoto": 0.4}]
        text = k.kba_run_readme("EGFR", PASSING_GATE, "Prose.", _controls(), queries=queries)
        assert "RDKit" not in text


class TestSimilarity:
    def test_returns_none_without_rdkit_or_a_number_with_it(self):
        sim = k.kba_tanimoto("c1ccccc1", "c1ccccc1")
        assert sim is None if not HAVE_RDKIT else sim == pytest.approx(1.0)

    @needs_rdkit
    def test_nearest_active_picks_the_closest(self):
        actives = [{"compound_id": "A", "smiles": "c1ccccc1"},
                   {"compound_id": "B", "smiles": "CCCCCCCCCCN"}]
        best_id, sim = k.kba_nearest_active("c1ccccc1C", actives)
        assert best_id == "A" and sim > 0

    def test_unparseable_smiles_gives_no_similarity(self):
        assert k.kba_tanimoto("not-a-molecule", "c1ccccc1") is None

    def test_nearest_active_on_empty_set(self):
        assert k.kba_nearest_active("c1ccccc1", []) == (None, None)
