"""Regression tests for residue_coverage's SIFTS-segment handling.

Covers three bugs found in review:
  1. offset math crashing when a segment boundary's author_residue_number is None
  2. only segs[0] being used, dropping later UniProt mapping segments
  3. silently falling back to another chain's mapping when the requested
     chain has none
"""
import kernel


PID = "9xyz"
ACC = "P00000"


def _graph_map(mappings):
    return {PID: {"UniProt": {ACC: {"mappings": mappings}}}}


def _res_list(entity_id, chain_id, residues):
    return {
        PID: {
            "molecules": [{
                "entity_id": entity_id,
                "chains": [{
                    "chain_id": chain_id,
                    "struct_asym_id": chain_id,
                    "residues": residues,
                }],
            }]
        }
    }


def _residue(author_num, ratio=1.0, name="ALA"):
    return {"author_residue_number": author_num, "observed_ratio": ratio, "residue_name": name}


def _fake_http_get(mapping_resp, reslist_resp):
    def fn(url, params=None, timeout=60, binary=False, retries=3):
        if url == kernel.PDBE_GRAPH_MAP + PID:
            return mapping_resp
        if url == kernel.PDBE_RES_LIST + PID:
            return reslist_resp
        raise AssertionError("unexpected url: %s" % url)
    return fn


def test_null_start_boundary_falls_back_to_end_boundary(monkeypatch):
    # start.author_residue_number is None (the 2RH1 case: boundary residue
    # is disordered); end.author_residue_number is present, so the offset
    # must be derived from the end boundary instead of crashing.
    mappings = [{
        "chain_id": "A", "entity_id": "1",
        "unp_start": 1, "unp_end": 50,
        "start": {"author_residue_number": None},
        "end": {"author_residue_number": 50},
    }]
    residues = [_residue(i) for i in range(1, 51)]
    monkeypatch.setattr(kernel, "http_get",
                         _fake_http_get(_graph_map(mappings),
                                        _res_list("1", "A", residues)))

    cov = kernel.residue_coverage(PID, ACC, chain="A")

    assert "error" not in cov
    assert cov["seq_offset"] == 0
    assert cov["residues"][0]["unp_num"] == 1


def test_both_boundaries_null_segment_is_skipped_not_fabricated(monkeypatch):
    # Neither boundary is usable -> no offset can be derived. Per the
    # project rule against fabricating scientific values, this must
    # surface as an explicit error, never a guessed/defaulted offset.
    mappings = [{
        "chain_id": "A", "entity_id": "1",
        "unp_start": 1, "unp_end": 50,
        "start": {"author_residue_number": None},
        "end": {"author_residue_number": None},
    }]
    monkeypatch.setattr(kernel, "http_get",
                         _fake_http_get(_graph_map(mappings), {}))

    cov = kernel.residue_coverage(PID, ACC, chain="A")

    assert "error" in cov


def test_multiple_segments_each_resolve_their_own_residues(monkeypatch):
    # Fusion-construct case (e.g. 2RH1 = receptor + T4-lysozyme insert):
    # two UniProt segments map onto the same chain with different offsets.
    # A residue must be resolved through the segment it actually belongs
    # to, not through segs[0]'s offset applied globally.
    mappings = [
        {
            "chain_id": "A", "entity_id": "1",
            "unp_start": 1, "unp_end": 50,
            "start": {"author_residue_number": 1},
            "end": {"author_residue_number": 50},
        },
        {
            "chain_id": "A", "entity_id": "1",
            "unp_start": 260, "unp_end": 300,
            "start": {"author_residue_number": 60},
            "end": {"author_residue_number": 100},
        },
    ]
    residues = ([_residue(i) for i in range(1, 51)] +
                [_residue(i) for i in range(60, 101)])
    monkeypatch.setattr(kernel, "http_get",
                         _fake_http_get(_graph_map(mappings),
                                        _res_list("1", "A", residues)))

    cov = kernel.residue_coverage(PID, ACC, chain="A")

    assert "error" not in cov
    by_author = {r["author_num"]: r for r in cov["residues"]}
    # First segment: offset 0, unp == author.
    assert by_author[25]["unp_num"] == 25
    assert by_author[25]["in_mapped_segment"] is True
    # Second segment: offset 200, unp == author + 200. If segs[0] were
    # applied globally (the bug) this would incorrectly come out as 75.
    assert by_author[75]["unp_num"] == 275
    assert by_author[75]["in_mapped_segment"] is True


def test_expression_tag_gap_excluded_from_n_missing(monkeypatch):
    # Regression: n_missing summed every gap's length regardless of kind, so
    # an unresolved construct tag (unp_num < 1, classified expression_tag and
    # explicitly documented as "not missing protein") silently inflated the
    # count of actually-missing protein residues.
    mappings = [{
        "chain_id": "A", "entity_id": "1",
        "unp_start": -5, "unp_end": 20,
        "start": {"author_residue_number": 1},
        "end": {"author_residue_number": 26},
    }]
    # author 1..6 = unresolved tag (unp -5..0); author 7..26 = observed protein.
    residues = [_residue(a, ratio=(0.0 if a <= 6 else 1.0)) for a in range(1, 27)]
    monkeypatch.setattr(kernel, "http_get",
                         _fake_http_get(_graph_map(mappings),
                                        _res_list("1", "A", residues)))

    cov = kernel.residue_coverage(PID, ACC, chain="A")

    assert "error" not in cov
    assert [g["kind"] for g in cov["gaps"]] == ["expression_tag"]
    assert cov["n_missing"] == 0


def test_missing_chain_errors_instead_of_falling_back(monkeypatch):
    # Requesting a chain with no SIFTS mapping must error, never silently
    # substitute another chain's segments/offsets.
    mappings = [{
        "chain_id": "A", "entity_id": "1",
        "unp_start": 1, "unp_end": 50,
        "start": {"author_residue_number": 1},
        "end": {"author_residue_number": 50},
    }]
    monkeypatch.setattr(kernel, "http_get",
                         _fake_http_get(_graph_map(mappings), {}))

    cov = kernel.residue_coverage(PID, ACC, chain="Z")

    assert "error" in cov
    assert "Z" in cov["error"]
