"""
Tests for beckmann/dft/parse_cmo.py's runner-up-MO/antibond-lookup extraction
(all_weight_matches_for_target, coefficient_in_mo) and the byte-identical-regression
guarantee on max_weight_for_target().

Synthetic vir_mos fixtures are hand-built dicts matching parse_cmo_table()'s
shape ({"mo", "kind", "energy", "contribs": [(coeff, label), ...]}) rather than
parsed from a real .log, so these tests don't depend on any DFT output being
present on disk.
"""
from beckmann_nbo.parse_cmo import (
    all_weight_matches_for_target, coefficient_in_mo, max_weight_for_target,
)


def _mo(mo, energy, contribs):
    return {"mo": mo, "kind": "vir", "energy": energy, "contribs": contribs}


def test_all_weight_matches_sorted_descending_by_weight():
    """Two MOs both carry BD*(C1-N2) -- matches must come back winner-first."""
    vir_mos = [
        _mo(10, -0.0500, [(-0.660, "BD*(1) C1-N2"), (0.20, "BD*(1) C1-C3")]),
        _mo(11, -0.0412, [(0.30, "BD*(1) C1-N2")]),
    ]
    matches = all_weight_matches_for_target(vir_mos, 1, 2)
    assert [m["mo_index"] for m in matches] == [10, 11]
    assert matches[0]["weight"] > matches[1]["weight"]
    assert matches[0]["coefficient"] == -0.660
    assert matches[1]["coefficient"] == 0.30


def test_all_weight_matches_empty_when_antibond_absent():
    vir_mos = [_mo(10, -0.05, [(-0.660, "BD*(1) C1-N2")])]
    assert all_weight_matches_for_target(vir_mos, 5, 6) == []


def test_max_weight_for_target_is_byte_identical_to_first_match():
    """max_weight_for_target()'s return must stay exactly matches[0]'s fields
    (plus the delta_lumo/in_window derivation) after the refactor to build on
    all_weight_matches_for_target() -- other code unpacks this positionally."""
    vir_mos = [
        _mo(10, -0.0500, [(-0.660, "BD*(1) C1-N2")]),
        _mo(11, -0.0412, [(0.30, "BD*(1) C1-N2")]),
    ]
    matches = all_weight_matches_for_target(vir_mos, 1, 2)
    result = max_weight_for_target(vir_mos, 1, 2, lumo_e=-0.06)
    weight, mo_index, epsilon, coeff, delta_lumo, in_window = result
    assert weight == matches[0]["weight"]
    assert mo_index == matches[0]["mo_index"]
    assert epsilon == matches[0]["epsilon_i_star"]
    assert coeff == matches[0]["coefficient"]
    assert delta_lumo == epsilon - (-0.06)
    assert in_window is True


def test_max_weight_for_target_none_when_no_match():
    vir_mos = [_mo(10, -0.05, [(-0.660, "BD*(1) C1-N2")])]
    assert max_weight_for_target(vir_mos, 5, 6) == (None, None, None, None, None, None)


def test_tie_break_keeps_first_seen_mo_as_winner():
    """Equal weights: original running-max used strict '>', so the first-encountered
    MO wins ties -- the sort-based rewrite must preserve that via stability."""
    vir_mos = [
        _mo(10, -0.05, [(0.5, "BD*(1) C1-N2")]),
        _mo(11, -0.04, [(-0.5, "BD*(1) C1-N2")]),
    ]
    matches = all_weight_matches_for_target(vir_mos, 1, 2)
    assert matches[0]["mo_index"] == 10


def test_coefficient_in_mo_looks_up_specific_mo():
    vir_mos = [
        _mo(10, -0.05, [(-0.660, "BD*(1) C1-N2"), (0.20, "BD*(1) C1-C3")]),
        _mo(11, -0.04, [(0.30, "BD*(1) C1-N2")]),
    ]
    assert coefficient_in_mo(vir_mos, 1, 3, 10) == 0.20
    assert coefficient_in_mo(vir_mos, 1, 3, 11) is None  # MO11 doesn't carry the aryl antibond
    assert coefficient_in_mo(vir_mos, 1, 3, None) is None
