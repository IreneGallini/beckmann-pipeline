"""
Validate beckmann.dft.branch_tracking against reference logs
(example_scans/5_s1_Me.log .. 5_s4_Me.log) and generate the branch diagrams.

Gate 1 (family-weight extraction, corrected/broadened per
Detailed_Orbital_Character_Exchange_Handout.docx) is checked against Section 6's
8 worked examples. Gate 2 (branch identity / candidate selection) reproduces
Section 5's table exactly at all 4 points now that track_branches()'s candidate
selection is restricted to positive-canonical-energy MOs -- see that function's
docstring for why (a negative-energy N-O 'activation coordinate' MO was
previously hijacking the ranking) and the caveat that this filter is an
empirical finding validated only against this one 4-point case, not yet
independently confirmed as Tetiana's actual rule. Both diagrams are still
generated so the two can be compared directly:
  - branch_tracking_reference_v2.png       : track_branches()'s ACTUAL output
  - branch_tracking_reference_v2_expected.png : Section 5's ground-truth table,
    plotted the same way, for direct visual comparison
"""
from pathlib import Path

from beckmann.config import DATA_OUTPUT
from beckmann.dft.branch_tracking import extract_family_weights, plot_branch_diagram, track_branches

EXAMPLE_DIR = Path(__file__).resolve().parents[2] / "example_scans"
PLOTS_DIR = DATA_OUTPUT / "analysis" / "plots"

# (log, mo, w_CC, w_CN, f_CC, f_CN) -- Detailed_Orbital_Character_Exchange_Handout.docx Section 6
GATE1_TARGETS = [
    ("5_s1_Me.log", 33, 0.000000, 0.779689, 0.000000, 1.000000),
    ("5_s1_Me.log", 35, 0.155236, 0.303601, 0.338325, 0.661675),
    ("5_s2_Me.log", 35, 0.158404, 0.303601, 0.342862, 0.657138),
    ("5_s2_Me.log", 38, 0.148225, 0.000000, 1.000000, 0.000000),
    ("5_s3_Me.log", 38, 0.146689, 0.000000, 1.000000, 0.000000),
    ("5_s3_Me.log", 35, 0.159201, 0.303601, 0.343994, 0.656006),
    ("5_s4_Me.log", 38, 0.220900, 0.000000, 1.000000, 0.000000),
    ("5_s4_Me.log", 35, 0.072900, 0.303601, 0.193625, 0.806375),
]

# Section 5's ground-truth table: (R_NO, MO1, E1, MO2, E2, MO_C7C8, E_C7C8)
GATE2_TARGETS = [
    (1.55, 33, 0.04779, 35, 0.10259, 39, 0.22551),
    (1.70, 35, 0.10039, 38, 0.21174, 39, 0.22475),
    (1.75, 38, 0.21117, 35, 0.09966, 39, 0.22440),
    (1.80, 38, 0.21060, 35, 0.09896, 39, 0.22402),
]


def run_gate1() -> bool:
    print("== Gate 1: family-weight extraction (corrected/broadened definitions) ==")
    all_pass = True
    for log_name, mo_index, w_cc, w_cn, f_cc, f_cn in GATE1_TARGETS:
        rows = extract_family_weights(EXAMPLE_DIR / log_name)
        row = next(r for r in rows if r["mo"] == mo_index)
        ok = (
            abs(row["w_CC"] - w_cc) < 1e-6 and abs(row["w_CN"] - w_cn) < 1e-6
            and abs(row["f_CC"] - f_cc) < 1e-6 and abs(row["f_CN"] - f_cn) < 1e-6
        )
        all_pass &= ok
        print(
            f"  {'PASS' if ok else 'FAIL'}: {log_name} MO{mo_index} -> "
            f"w_CC={row['w_CC']:.6f} w_CN={row['w_CN']:.6f} "
            f"f_CC={row['f_CC']:.6f} f_CN={row['f_CN']:.6f}"
        )
    print(f"Gate 1: {'PASS' if all_pass else 'FAIL'}\n")
    return all_pass


def run_gate2() -> bool:
    print("== Gate 2: branch tracking (current implementation, unmodified) ==")
    logs = [EXAMPLE_DIR / f"5_s{i}_Me.log" for i in range(1, 5)]
    results = track_branches(logs)

    all_pass = True
    for (r_no, mo1, e1, mo2, e2, mo_ref, e_ref), r in zip(GATE2_TARGETS, results):
        got = {r["mo_A"], r["mo_B"]}
        expected = {mo1, mo2}
        ok = got == expected
        all_pass &= ok
        print(
            f"  {'PASS' if ok else 'FAIL'}: R={r_no} -> got A=MO{r['mo_A']} ({r['E_A']:.5f} au) "
            f"B=MO{r['mo_B']} ({r['E_B']:.5f} au), expected MO{mo1}/MO{mo2}"
        )
    print(f"Gate 2: {'PASS' if all_pass else 'FAIL'}\n")
    return all_pass


def _ref_series() -> list[tuple[int, float]]:
    """Best single MO by w_ref (C7-C8 antibond weight) per point -- the separate
    third quantity, not part of the two-branch selection."""
    series = []
    for i in range(1, 5):
        rows = extract_family_weights(EXAMPLE_DIR / f"5_s{i}_Me.log")
        best = max(rows, key=lambda r: r["w_ref"])
        series.append((best["mo"], best["energy"]))
    return series


def plot_actual() -> None:
    r_values = [r_no for r_no, *_ in GATE2_TARGETS]
    logs = [EXAMPLE_DIR / f"5_s{i}_Me.log" for i in range(1, 5)]
    results = track_branches(logs)
    branch_a = [(r["mo_A"], r["E_A"]) for r in results]
    branch_b = [(r["mo_B"], r["E_B"]) for r in results]
    branch_ref = _ref_series()

    fig = plot_branch_diagram(
        r_values, branch_a, branch_b, branch_ref=branch_ref,
        title="Weighted orbital-character tracking (ACTUAL code output)",
        x_label="R(N-O) (Å)", y_label="Canonical MO eigenvalue (a.u.)",
        numeric_x=True,
    )
    out_path = PLOTS_DIR / "branch_tracking_reference_v2.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    print(f"-- wrote {out_path}")


def plot_expected() -> None:
    r_values = [r_no for r_no, *_ in GATE2_TARGETS]
    branch_a = [(mo1, e1) for _, mo1, e1, _, _, _, _ in GATE2_TARGETS]
    branch_b = [(mo2, e2) for _, _, _, mo2, e2, _, _ in GATE2_TARGETS]
    branch_ref = [(mo_ref, e_ref) for _, _, _, _, _, mo_ref, e_ref in GATE2_TARGETS]

    fig = plot_branch_diagram(
        r_values, branch_a, branch_b, branch_ref=branch_ref,
        title="Weighted orbital-character tracking (Section 5 ground truth)",
        x_label="R(N-O) (Å)", y_label="Canonical MO eigenvalue (a.u.)",
        numeric_x=True,
    )
    out_path = PLOTS_DIR / "branch_tracking_reference_v2_expected.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    print(f"-- wrote {out_path}")


if __name__ == "__main__":
    run_gate1()
    run_gate2()
    plot_actual()
    plot_expected()
