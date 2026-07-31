"""
Benchmark/reference-case runners for the PySCF wCNmax validation work.
Originally part of beckmann_alt/pair_nbo.py itself; split out during the
monorepo restructuring since these are benchmark-harness convenience
wrappers (load a hardcoded case, call the engine), not part of the
validated wCNmax computation engine (build_local_iaos/compute_wcnmax/
run_from_case, all untouched, now living in
beckmann_pyscf.engine.pair_nbo). Keeping them there would have required the
product's own engine module to depend on this research package's
benchmark-log loaders -- backwards for a supposedly standalone product
package. research/ depending on the product's engine (the direction here)
is the correct one.
"""
from beckmann_pyscf.engine.pair_nbo import run_from_case

from pyscf_validation.geometry import (
    TEST_IDS, load_case, load_test_set_case, load_test_set_scan_series,
)


def run_case(name: str) -> dict:
    """One of the two hand-picked reference cases (mol_002, 5_s0_Me) -- see
    pyscf_validation.geometry.REFERENCE_CASES."""
    return run_from_case(load_case(name))


def run_test_set_case(mol_id: str) -> dict:
    """Any of the six main-pipeline test-set molecules (mol_002/006/014/020/021/029) --
    see pyscf_validation.geometry.load_test_set_case."""
    return run_from_case(load_test_set_case(mol_id))


def run_test_set_scan_series(mol_id: str, stages: list[str] | None = None) -> list[dict]:
    """wCNmax at every R(N-O) point of a test-set molecule's scan series
    (pyscf_validation.geometry.load_test_set_scan_series) -- one PySCF SCF +
    local per-atom-pair projection per point, not just the single
    equilibrium geometry run_test_set_case() computes. Returns one row per
    point shaped like a beckmann_nbo.parse_cmo 'cn'-channel extraction row
    (mol/stage/channel/R_NO/MO_index/weight/...) so
    beckmann_core.wcnmax_rule.find_wcnmax_minimum() can be called on the
    result directly, reusing the main pipeline's own interior-minimum
    criterion rather than reimplementing it here.

    stages, if given, restricts the run to just those stage labels (e.g.
    mol_034_E's STEP_SCAN_SOURCES-merged series has 12 points, ~2x every
    other molecule's 6 -- pass a 6-point subset to keep runtime comparable).
    """
    cases = load_test_set_scan_series(mol_id)
    if stages is not None:
        cases = [c for c in cases if c["stage"] in stages]
    rows = []
    for case in cases:
        result = run_from_case(case)
        cn = result["cn"]
        second = cn["second"]
        aryl_coeffs = result["aryl_coeffs"]
        rows.append({
            "mol": case["name"], "stage": case["stage"], "channel": "cn",
            "R_NO": case["r_no"], "MO_index": cn["mo_index"],
            "epsilon_i_star": cn["epsilon"], "coefficient": cn["coefficient"],
            "weight": cn["wmax"], "delta_lumo": None, "in_window": None,
            "MO_A": cn["mo_index"], "MO_B": second["mo_index"] if second else None,
            "eps_A": cn["epsilon"], "eps_B": second["epsilon"] if second else None,
            "CN_coeff_in_A": cn["coefficient"],
            "CN_coeff_in_B": second["coefficient"] if second else None,
            "arylCC_coeff_in_A": aryl_coeffs.get(cn["mo_index"]),
            "arylCC_coeff_in_B": aryl_coeffs.get(second["mo_index"]) if second else None,
        })
    return rows


def _print_wcnmax(result: dict) -> None:
    r, p = result["cn"], result["cn"]["pair"]
    print(f"\n=== {result['case']} ({result['basis_note']}) -- wCNmax, local per-atom-pair antibond ===")
    print(
        f"  wCNmax: wmax={r['wmax']:.4f}  MO={r['mo_index']}  "
        f"[local block: {p['n_local_iaos']} IAOs, {len(p['candidates_ao'])} antibond candidate(s), "
        f"winning_occ={r['antibond_occupation']:.4f}, bond_occ={p['bond_occupation']:.4f}]"
    )


if __name__ == "__main__":
    import sys
    name = sys.argv[1] if len(sys.argv) > 1 else "mol_002"
    result = run_test_set_case(name) if name in TEST_IDS else run_case(name)
    _print_wcnmax(result)
