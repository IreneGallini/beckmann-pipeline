"""
Tests for scripts/analysis/classical_rule_benchmark.py output.

Validates: data/output/analysis/classical_rule_results.csv
This is not a pipeline step — it's a one-off analysis showing that the classical
anti-periplanar dihedral rule alone is insufficient to predict Beckmann outcomes,
which motivates the DFT/NBO approach.
"""
import csv
import pytest


REQUIRED_COLUMNS = {
    "mol_id", "beckmann_pred", "exp_outcome", "agreement",
    "dihedral_O_N_C_aryl", "dihedral_O_N_C_allyl",
    "Emin_E_eV", "Emin_Z_eV", "lowest_isomer",
}


@pytest.fixture(scope="module")
def classical_rule_csv(project_root):
    p = project_root / "data" / "output" / "analysis" / "classical_rule_results.csv"
    if not p.exists():
        pytest.skip(
            "classical_rule_results.csv not found — "
            "run scripts/analysis/classical_rule_benchmark.py first"
        )
    return p


def _rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def test_csv_has_required_columns(classical_rule_csv):
    with open(classical_rule_csv) as f:
        cols = set(csv.DictReader(f).fieldnames or [])
    missing = REQUIRED_COLUMNS - cols
    assert not missing, f"Missing columns: {missing}"


def test_csv_row_count(classical_rule_csv, benchmark_meta):
    rows = _rows(classical_rule_csv)
    assert len(rows) == len(benchmark_meta), (
        f"Expected {len(benchmark_meta)} rows, got {len(rows)}"
    )


def test_prediction_is_r_or_f(classical_rule_csv):
    for row in _rows(classical_rule_csv):
        assert row["beckmann_pred"] in ("R", "F", "inspect"), (
            f"{row['mol_id']}: beckmann_pred={row['beckmann_pred']!r}"
        )


def test_dihedral_range(classical_rule_csv):
    for row in _rows(classical_rule_csv):
        for col in ("dihedral_O_N_C_aryl", "dihedral_O_N_C_allyl"):
            if not row[col]:
                continue
            val = float(row[col])
            assert -180.0 <= val <= 180.0, (
                f"{row['mol_id']}: {col}={val} is outside [-180, 180]"
            )


def test_agreement_values(classical_rule_csv):
    for row in _rows(classical_rule_csv):
        assert row["agreement"] in ("yes", "no", "inspect", "unclear"), (
            f"{row['mol_id']}: agreement={row['agreement']!r}"
        )


def test_exp_outcome_matches_meta(classical_rule_csv, benchmark_meta):
    for row in _rows(classical_rule_csv):
        mol_id = row["mol_id"]
        if mol_id not in benchmark_meta or not row["exp_outcome"]:
            continue
        expected = benchmark_meta[mol_id]["exp_outcome"]
        assert row["exp_outcome"] == expected, (
            f"{mol_id}: CSV exp_outcome={row['exp_outcome']!r} != meta {expected!r}"
        )


def test_accuracy_is_not_perfect(classical_rule_csv):
    """The classical rule should not achieve 100% — that's the whole point."""
    rows = _rows(classical_rule_csv)
    agree = sum(1 for r in rows if r["agreement"] == "yes")
    assert agree < len(rows), (
        "Classical rule predicted all molecules correctly — "
        "check that experimental outcomes are loaded properly"
    )
