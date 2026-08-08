import sys

from kernel import strip_entrez, bh_adjust, floor_pvalue


def test_strip_entrez_removes_suffix():
    assert strip_entrez("KRAS (3845)") == "KRAS"


def test_strip_entrez_idempotent_on_clean_name():
    assert strip_entrez("KRAS") == "KRAS"


def test_strip_entrez_only_strips_trailing_paren_number():
    # A gene symbol that legitimately contains parentheses-like text elsewhere
    # must not be mangled; only a trailing " (<digits>)" is an Entrez suffix.
    assert strip_entrez("HLA-A (3105)") == "HLA-A"
    assert strip_entrez("C1orf100") == "C1orf100"


def test_bh_adjust_matches_hand_computed_example():
    # Classic textbook example: 5 p-values, alpha irrelevant here.
    # raw:            [0.01, 0.02, 0.03, 0.04, 0.50]
    # rank (1-indexed): 1     2     3     4     5
    # p*n/rank:        0.05  0.05  0.05  0.05  0.50
    # monotone from the top down: min-accumulate from the largest rank downward
    # -> [0.05, 0.05, 0.05, 0.05, 0.50]
    raw = [0.01, 0.02, 0.03, 0.04, 0.50]
    expected = [0.05, 0.05, 0.05, 0.05, 0.50]
    out = bh_adjust(raw)
    for got, exp in zip(out, expected):
        assert abs(got - exp) < 1e-12


def test_bh_adjust_returns_values_aligned_to_input_order_not_sorted_order():
    # Feed p-values in a shuffled, non-sorted order; the output must line up
    # index-for-index with the INPUT order, not come back pre-sorted.
    raw = [0.50, 0.01, 0.04, 0.02, 0.03]
    out = bh_adjust(raw)
    # index 1 (p=0.01, the smallest) must carry the smallest adjusted value.
    assert out[1] == min(out)
    # index 0 (p=0.50, the largest) must carry the largest adjusted value.
    assert out[0] == max(out)


def test_bh_adjust_never_exceeds_one():
    raw = [0.9, 0.95, 0.99, 0.999]
    out = bh_adjust(raw)
    assert all(0 <= v <= 1 for v in out)


def test_floor_pvalue_clamps_zero_to_smallest_representable():
    assert floor_pvalue(0.0) == sys.float_info.min


def test_floor_pvalue_passes_through_positive_values():
    assert floor_pvalue(0.037) == 0.037
