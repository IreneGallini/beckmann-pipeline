"""
Validate beckmann.dft.diabatic_character against Tetiana's reference logs
(example_scans/5_s1_Me.log .. 5_s4_Me.log).

Gate 1 (reference case, hardcoded atom numbers -- the only path possible for
this molecule, which has no .gjf oxime label or best_per_substrate.sdf entry):
per-point argmax-w_CN / argmax-w_CC MO selection, checked against values
computed directly from extract_family_weights() (Gate-1-validated against
Detailed_Orbital_Character_Exchange_Handout.docx Section 6 -- see
beckmann/dft/diabatic_character.py's module docstring). Also demonstrates the
character-exchange pattern directly: the max-w_CC MO's own f_CN runs high
(N-side/mixed) at s1-s3 and drops to 0 (C-C-routed) at s4.

Gate 2 (atom-map generalization smoke test): confirms get_substituent_map()
resolves mol_002_E's real atom numbers correctly (cross-checked against the
known values already asserted in tests/test_descriptors.py) and that
track_diabatic_character() runs against its single-point _nbo.log without
error via the generalized (mol, mol_dir) path.

Gate 3 (real multi-point scan parsing): sanity-checks
extract_family_weights_series() against mol_002_E's actual 5-point
_scan.log (all points chained via --Link1-- into one log, unlike Gate 1's
reference case) -- confirms a sensible point count and that the first and
last scan points carry real, nonzero family-weight data rather than parsing
garbage. This is what track_diabatic_character_series() (see
beckmann/dft/diabatic_character.py) is built on for the benchmark-wide run.
"""
from pathlib import Path

from beckmann.config import DATA_OUTPUT
from beckmann.dft.descriptors import get_substituent_map
from beckmann.dft.diabatic_character import (
    extract_family_weights_series, track_diabatic_character,
)

EXAMPLE_DIR = Path(__file__).resolve().parents[2] / "example_scans"

# (point, mo_CN, E_CN, f_CN_CN, mo_CC, E_CC, f_CC_CC, f_CN_CC) -- computed directly
# from extract_family_weights()
GATE1_TARGETS = [
    ("s1", 33, 0.04779, 1.0, 35, 0.10259, 0.338325, 0.661675),
    ("s2", 32, -0.00951, 1.0, 35, 0.10039, 0.342862, 0.657138),
    ("s3", 32, -0.02490, 1.0, 35, 0.09966, 0.343994, 0.656006),
    ("s4", 32, -0.03881, 1.0, 38, 0.21060, 1.000000, 0.000000),
]

MOL_002E_DIR = DATA_OUTPUT / "dft_opt" / "mol_002_E"
MOL_002E_EXPECTED_SUBST = {"ci": 11, "ni": 12, "oi": 13, "c_aryl": 6, "c_alkyl": 10}


def run_gate1() -> bool:
    print("== Gate 1: reference-case argmax selection (hardcoded atom numbers) ==")
    logs = [EXAMPLE_DIR / f"5_{point}_Me.log" for point, *_ in GATE1_TARGETS]
    results = track_diabatic_character(logs)

    all_pass = True
    for (point, mo_cn, e_cn, f_cn_cn, mo_cc, e_cc, f_cc_cc, f_cn_cc), r in zip(GATE1_TARGETS, results):
        ok = (
            r["mo_CN"] == mo_cn and abs(r["E_CN"] - e_cn) < 1e-5 and abs(r["f_CN_CN"] - f_cn_cn) < 1e-6
            and r["mo_CC"] == mo_cc and abs(r["E_CC"] - e_cc) < 1e-5
            and abs(r["f_CC_CC"] - f_cc_cc) < 1e-6 and abs(r["f_CN_CC"] - f_cn_cc) < 1e-6
        )
        all_pass &= ok
        print(
            f"  {'PASS' if ok else 'FAIL'}: {point} -> "
            f"mo_CN=MO{r['mo_CN']} (E={r['E_CN']:.5f}) "
            f"mo_CC=MO{r['mo_CC']} (E={r['E_CC']:.5f}, f_CN={r['f_CN_CC']:.4f})"
        )
    print(f"Gate 1: {'PASS' if all_pass else 'FAIL'}\n")

    print("  Character-exchange pattern (max-w_CC MO's own f_CN, s1->s4):")
    print("   ", [round(r["f_CN_CC"], 4) for r in results])
    print("   N-side/mixed (s1-s3) -> C-C-routed (s4), as expected.\n")
    return all_pass


def run_gate2() -> bool:
    print("== Gate 2: atom-map generalization smoke test (mol_002_E) ==")
    subst = get_substituent_map("mol_002_E", MOL_002E_DIR)
    ok_subst = subst == MOL_002E_EXPECTED_SUBST
    print(f"  {'PASS' if ok_subst else 'FAIL'}: get_substituent_map -> {subst}")

    rows = track_diabatic_character(
        [MOL_002E_DIR / "mol_002_E_nbo.log"], mol="mol_002_E", mol_dir=MOL_002E_DIR,
    )
    row = rows[0]
    ok_row = row["mo_CN"] is not None and row["mo_CC"] is not None
    print(f"  {'PASS' if ok_row else 'FAIL'}: track_diabatic_character (generalized path) -> {row}")

    all_pass = ok_subst and ok_row
    print(f"Gate 2: {'PASS' if all_pass else 'FAIL'}\n")
    return all_pass


def run_gate3() -> bool:
    print("== Gate 3: real multi-point scan parsing (mol_002_E _scan.log) ==")
    subst = MOL_002E_EXPECTED_SUBST
    series = extract_family_weights_series(
        "mol_002_E", MOL_002E_DIR,
        c1_atom=subst["c_aryl"], cn_c_atom=subst["ci"], cn_n_atom=subst["ni"], ref_atom=subst["c_alkyl"],
    )
    print(f"  {len(series)} scan point(s) found (expect ~7: 1 nbo/R0 + 6 scan points)")

    if not series:
        print("Gate 3: FAIL (no points found)\n")
        return False

    ok = True
    for label, point in (("first", series[0]), ("last", series[-1])):
        rows = point["rows"]
        max_wcc = max((r["w_CC"] for r in rows), default=0.0)
        max_wcn = max((r["w_CN"] for r in rows), default=0.0)
        energies = [r["energy"] for r in rows]
        sane = (
            bool(rows) and (max_wcc > 0 or max_wcn > 0)
            and all(-1.0 < e < 1.0 for e in energies)
        )
        ok &= sane
        print(
            f"  {'PASS' if sane else 'FAIL'}: {label} point R(N-O)={point['r_no']:.4f} A -- "
            f"{len(rows)} virtual MO(s), max_w_CC={max_wcc:.4f}, max_w_CN={max_wcn:.4f}, "
            f"energy range=[{min(energies):.4f}, {max(energies):.4f}]"
        )
    print(f"Gate 3: {'PASS' if ok else 'FAIL'}\n")
    return ok


if __name__ == "__main__":
    run_gate1()
    run_gate2()
    run_gate3()
