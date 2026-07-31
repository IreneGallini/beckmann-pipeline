"""
Unit tests for beckmann_core.wcnmax_rule against synthetic row data (no
filesystem dependency -- beckmann-core stays path-agnostic, see its own
module docstrings). Separately, during this session's migration, this same
logic was regression-checked against beckmann-pipeline's real committed
data/output/analysis/{cmo_channel_extraction,wcnmax_rule_results}.csv and
matched the trusted 'predicted' column on all 34 benchmark molecules
(0 mismatches) -- not part of this package's own test suite since it
depends on data this package doesn't ship, but recorded here as the
migration's validation evidence.
"""
from beckmann_core.wcnmax_rule import (
    find_wcnmax_extremum, find_wcnmax_minimum, predict_from_wcnmax, resolve_series,
)


def _row(mol, stage, r_no, weight, mo_index=10, eps=-0.05):
    return {
        "mol": mol, "stage": stage, "channel": "cn",
        "R_NO": r_no, "weight": weight, "MO_index": mo_index, "epsilon_i_star": eps,
    }


def test_resolve_series_orders_nbo_then_scan_numerically():
    by_stage = {
        "scan_2": {"tag": "scan_2"}, "nbo": {"tag": "nbo"},
        "scan_10": {"tag": "scan_10"}, "scan_1": {"tag": "scan_1"},
    }
    series = resolve_series(by_stage)
    assert [row["tag"] for row in series] == ["nbo", "scan_1", "scan_2", "scan_10"]


def test_find_wcnmax_minimum_detects_genuine_dip():
    mol = "test_mol"
    rows = [
        _row(mol, "nbo", 1.60, 0.40),
        _row(mol, "scan_1", 1.65, 0.30),
        _row(mol, "scan_2", 1.70, 0.10),  # interior minimum
        _row(mol, "scan_3", 1.75, 0.35),
        _row(mol, "scan_4", 1.80, 0.45),
    ]
    minimum = find_wcnmax_minimum(mol, rows)
    assert minimum is not None
    assert minimum["R_star"] == 1.70
    assert minimum["w_star"] == 0.10
    assert minimum["depth"] > 0
    assert predict_from_wcnmax(minimum) == "R"


def test_find_wcnmax_minimum_none_for_monotonic_series():
    mol = "test_mol"
    rows = [_row(mol, "nbo", 1.60, 0.10 + 0.05 * i) for i in range(5)]
    for i, row in enumerate(rows[1:], start=1):
        row["stage"] = f"scan_{i}"
        row["R_NO"] = 1.60 + 0.05 * i
    minimum = find_wcnmax_minimum(mol, rows)
    assert minimum is None
    assert predict_from_wcnmax(minimum) == "F"


def test_find_wcnmax_minimum_ignores_local_maximum():
    """A local MAXIMUM (depth < 0) must not be reported as a minimum."""
    mol = "test_mol"
    rows = [
        _row(mol, "nbo", 1.60, 0.10),
        _row(mol, "scan_1", 1.65, 0.40),  # interior maximum, not a minimum
        _row(mol, "scan_2", 1.70, 0.12),
    ]
    assert find_wcnmax_minimum(mol, rows) is None
    assert find_wcnmax_extremum(mol, rows) is not None  # extremum exists, just not a minimum


def test_find_wcnmax_minimum_requires_at_least_three_points():
    mol = "test_mol"
    rows = [_row(mol, "nbo", 1.60, 0.40), _row(mol, "scan_1", 1.70, 0.10)]
    assert find_wcnmax_minimum(mol, rows) is None


def test_find_wcnmax_minimum_picks_deepest_extremum_not_first():
    mol = "test_mol"
    rows = [
        _row(mol, "nbo", 1.60, 0.30),
        _row(mol, "scan_1", 1.65, 0.28),   # shallow wobble (depth ~0.02)
        _row(mol, "scan_2", 1.70, 0.31),
        _row(mol, "scan_3", 1.75, 0.05),   # deep genuine minimum (depth much larger)
        _row(mol, "scan_4", 1.80, 0.32),
    ]
    minimum = find_wcnmax_minimum(mol, rows)
    assert minimum is not None
    assert minimum["R_star"] == 1.75


def test_find_wcnmax_minimum_ignores_other_channels_and_other_molecules():
    rows = [
        _row("mol_A", "nbo", 1.60, 0.40),
        _row("mol_A", "scan_1", 1.65, 0.10),
        _row("mol_A", "scan_2", 1.70, 0.45),
        _row("mol_B", "nbo", 1.60, 0.10),
        _row("mol_B", "scan_1", 1.65, 0.05),
        _row("mol_B", "scan_2", 1.70, 0.02),
    ]
    other_channel = dict(rows[1]); other_channel["channel"] = "cc"; other_channel["weight"] = 0.99
    assert find_wcnmax_minimum("mol_A", rows + [other_channel]) is not None
    assert find_wcnmax_minimum("mol_B", rows) is None
