"""
Tests for script 04_extract_dihedrals_and_predict.py output.

Validates: data/output/week1_benchmark_results.csv
Checks completeness, value ranges, and consistency with benchmark_meta.json.
"""
import csv


REQUIRED_COLUMNS = {
    "mol_id", "beckmann_pred", "exp_outcome", "agreement",
    "dihedral_O_N_C_aryl", "dihedral_O_N_C_allyl",
    "Emin_E_eV", "Emin_Z_eV", "lowest_isomer",
}


def _rows(benchmark_csv_path):
    with open(benchmark_csv_path) as f:
        return list(csv.DictReader(f))


def test_benchmark_csv_exists(benchmark_csv_path):
    assert benchmark_csv_path.exists()


def test_benchmark_csv_has_required_columns(benchmark_csv_path):
    with open(benchmark_csv_path) as f:
        cols = set(csv.DictReader(f).fieldnames or [])
    missing = REQUIRED_COLUMNS - cols
    assert not missing, f"Missing columns in benchmark CSV: {missing}"


def test_benchmark_csv_row_count(benchmark_csv_path, benchmark_meta):
    rows = _rows(benchmark_csv_path)
    assert len(rows) == len(benchmark_meta), (
        f"Expected {len(benchmark_meta)} rows (one per molecule), got {len(rows)}"
    )


def test_prediction_is_r_or_f(benchmark_csv_path):
    for row in _rows(benchmark_csv_path):
        assert row["beckmann_pred"] in ("R", "F"), (
            f"{row['mol_id']}: beckmann_pred={row['beckmann_pred']!r}, expected 'R' or 'F'"
        )


def test_dihedral_range(benchmark_csv_path):
    for row in _rows(benchmark_csv_path):
        for col in ("dihedral_O_N_C_aryl", "dihedral_O_N_C_allyl"):
            val = float(row[col])
            assert -180.0 <= val <= 180.0, (
                f"{row['mol_id']}: {col}={val} is outside [-180, 180]"
            )


def test_agreement_is_yes_or_no(benchmark_csv_path):
    for row in _rows(benchmark_csv_path):
        assert row["agreement"] in ("yes", "no"), (
            f"{row['mol_id']}: agreement={row['agreement']!r}, expected 'yes' or 'no'"
        )


def test_lowest_isomer_is_e_or_z(benchmark_csv_path):
    for row in _rows(benchmark_csv_path):
        assert row["lowest_isomer"] in ("E", "Z"), (
            f"{row['mol_id']}: lowest_isomer={row['lowest_isomer']!r}"
        )


def test_exp_outcome_matches_meta(benchmark_csv_path, benchmark_meta):
    """exp_outcome in the CSV must agree with benchmark_meta.json (source of truth)."""
    for row in _rows(benchmark_csv_path):
        mol_id = row["mol_id"]
        if mol_id not in benchmark_meta:
            continue
        expected = benchmark_meta[mol_id]["exp_outcome"]
        assert row["exp_outcome"] == expected, (
            f"{mol_id}: CSV exp_outcome={row['exp_outcome']!r} "
            f"!= benchmark_meta {expected!r}"
        )


def test_all_meta_molecules_in_csv(benchmark_csv_path, benchmark_meta):
    csv_ids = {row["mol_id"] for row in _rows(benchmark_csv_path)}
    missing = set(benchmark_meta.keys()) - csv_ids
    assert not missing, f"Molecules in benchmark_meta missing from CSV: {sorted(missing)}"