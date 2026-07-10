"""
Tests for beckmann/ts_ml/neb.py output.

Pilot scope only: mol_002_E rearrangement TS (TS1_A1) via the AIMNet2-nse/
PySisyphus proxy route. Gated on the pilot output existing.
"""
import pytest
import numpy as np
import h5py

from beckmann.ts_ml.neb import verify_ts_ml


@pytest.fixture(scope="module")
def ts1_a1_hessian(project_root):
    p = (project_root / "data" / "output" / "ts_ml" / "mol_002_E_ts1_a1"
         / "ts_final_hessian.h5")
    if not p.exists():
        pytest.skip("mol_002_E_ts1_a1 AIMNet2-nse proxy run not found — "
                     "run beckmann.ts_ml.neb main() first")
    return p


def test_hessian_has_expected_fields(ts1_a1_hessian):
    with h5py.File(ts1_a1_hessian, "r") as f:
        for key in ("vibfreqs", "mw_cart_displs", "masses"):
            assert key in f


def test_verify_ts_ml_report_shape(ts1_a1_hessian):
    # mol_002_E atom map (0-based): ci=10, ni=11, oi=12, c_aryl=5, c_alkyl=9
    atoms_of_interest = {"ci": 10, "ni": 11, "oi": 12, "c_aryl": 5, "c_alkyl": 9}
    report = verify_ts_ml(ts1_a1_hessian, atoms_of_interest)
    assert report["n_imaginary_total"] >= 1
    assert report["dominant_mode_freq_cm1"] < 0
    assert len(report["top_displaced_atoms"]) <= 6
    assert set(report["atoms_of_interest_in_top6"]) == set(atoms_of_interest)


def test_dominant_mode_displaces_reaction_center(ts1_a1_hessian):
    """The dominant imaginary mode should implicate N, O, and the oxime carbon
    for ANY channel through this reaction center -- regardless of which specific
    migrating/breaking group is involved."""
    atoms_of_interest = {"ci": 10, "ni": 11, "oi": 12, "c_aryl": 5, "c_alkyl": 9}
    report = verify_ts_ml(ts1_a1_hessian, atoms_of_interest)
    top6 = report["atoms_of_interest_in_top6"]
    assert top6["ni"] and top6["oi"] and top6["ci"], (
        "dominant imaginary mode does not implicate the N-O reaction center at all"
    )
