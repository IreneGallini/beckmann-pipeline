"""
End-to-end pipeline tests on 3 known benchmark molecules (mol_002, mol_006,
mol_020 -- the ones Notes_pyscf_alt.md already has trusted PySCF-on-
Gaussian-geometry numbers for, in beckmann-pipeline's
data/output/analysis/wcnmax_rule_results_pyscf.csv).

Important: this pipeline's geometry comes from AIMNet2 alone (no Gaussian
step at all), while the existing PySCF benchmark ran on Gaussian's
DFT-converged geometry. Exact wCNmax value / R-F call match against that
reference is NOT the target -- see the beckmann-pyscf README's Validation
section. These tests hard-assert structural correctness (the pipeline runs
to completion and produces a well-shaped result) and print the reference
comparison for visual inspection (run with `-s`) rather than silently
passing or failing on an expected geometry-source divergence.
"""
import csv
from pathlib import Path

import pytest

from beckmann_pyscf import pipeline

# id -> (SMILES, exp_outcome) from beckmann-pipeline's data/input/benchmark.csv
BENCHMARK_CASES = {
    "mol_002": ("O=C1CCC2=C1C=CC(OC)=C2", "F"),
    "mol_006": ("O=C1CCC2=C1C=CC(C)=C2", "R"),
    "mol_020": ("O=C1CCCC2=C1C=CC(C)=C2", "R"),
}

# Sibling repo's reference CSV, if present on this machine -- not shipped in
# beckmann-pyscf itself (it's beckmann-pipeline's own benchmark artifact).
_REFERENCE_CSV = (
    Path(__file__).resolve().parent.parent.parent.parent
    / "beckmann-pipeline" / "data" / "output" / "analysis" / "wcnmax_rule_results_pyscf.csv"
)


def _load_reference() -> dict[str, dict]:
    if not _REFERENCE_CSV.exists():
        return {}
    with open(_REFERENCE_CSV) as f:
        return {row["mol"]: row for row in csv.DictReader(f)}


@pytest.mark.slow
@pytest.mark.parametrize("mol_id", sorted(BENCHMARK_CASES))
def test_predict_end_to_end(mol_id, tmp_path):
    smiles, exp_outcome = BENCHMARK_CASES[mol_id]
    result = pipeline.predict(smiles, workdir=tmp_path)

    # Structural correctness -- the actual thing this test hard-asserts.
    assert result["prediction"] in ("R", "F")
    assert result["sdf_block"]
    assert result["energy_ev"] is not None
    series = result["wcnmax_series"]
    assert len(series) == 7, f"expected 7 scan points (nbo + scan_1..6), got {len(series)}"
    r_values = [pt["R_NO"] for pt in series]
    assert r_values == sorted(r_values), "R(N-O) series must be strictly increasing"
    assert all(0.0 < pt["weight"] < 1.0 for pt in series), "wCNmax weight out of physical range"

    # Reference comparison -- reported, not asserted (see module docstring).
    reference = _load_reference().get(f"{mol_id}_E")
    print(f"\n--- {mol_id} ({smiles}) ---")
    print(f"experimental outcome:    {exp_outcome}")
    print(f"this pipeline (AIMNet2 geometry): {result['prediction']}")
    if result["minimum"]:
        m = result["minimum"]
        print(f"  interior minimum: R*={m['R_star']:.4f} A, w*={m['w_star']:.4f}, depth={m['depth']:.4f}")
    else:
        print("  no interior minimum found")
    if reference:
        print(f"beckmann-pipeline PySCF-on-Gaussian-geometry: {reference['PySCF']} "
              f"(R_star={reference['R_star'] or 'n/a'}, depth={reference['R_depth'] or 'n/a'})")
        print(f"beckmann-pipeline NBO7 (trusted):             {reference['NBO']}")
        if reference["PySCF"] != result["prediction"]:
            print("  NOTE: R/F call differs from the Gaussian-geometry PySCF reference "
                  "-- expected given the different geometry source, see README.")
    else:
        print("(no beckmann-pipeline reference CSV found on this machine -- skipping comparison)")
