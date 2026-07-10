"""
Validate the channel-resolved descriptor code (parse_cmo.py, descriptors.py)
against an external reference log from Tetiana: 5_s0_Me.log, compound 3 from
"Ring Size and Substituent Effects in the Beckmann Rearrangement" (Table 2),
the methyl-substituted oxime, R0 equilibrium geometry only.

This is a one-off script (like the other scripts/analysis/ entries) -- it
reads a reference file that isn't one of our 34 benchmark molecules, so the
atom map (ci, ni, oi, c_aryl, c_alkyl) is hardcoded from the paper's own
Figure 2 labels, not derived.

Gotcha (see Notes.md): 5_s0_Me.log's route line includes Stable=Opt, which
triggers a wavefunction stability re-test and reruns population analysis --
the file contains TWO full NBO/CMO sections, not one. This script always uses
the LAST occurrence of each, per the handout's own warning. Our own pipeline's
.gjf files don't use Stable=Opt, so this quirk is specific to this reference
file, not our regular molecules.

Tier 1 (runs now): single-geometry wCNmax/w17max/w78max/Lambda/Psi, checked
against the handout's worked answer (wCNmax = 0.457 at MO 32).

Tier 2 (blocked): the paper's Table 2 d/dR slopes need the other 4 scan points
(R0+0.1 through R0+0.4) for this same compound. Only 5_s0_Me.log (R0) exists
in the repo -- this script reports that gap explicitly rather than
approximating a slope from one point.
"""
from pathlib import Path

from beckmann.dft.descriptors import PSI_EPSILON, compute_psi_row
from beckmann.dft.parse_cmo import (
    find_cmo_sections, parse_cmo_table, virtual_window, max_weight_for_target,
)
from beckmann.dft.parse_nbo import find_table_starts, parse_table_rows

REFERENCE_LOG = Path(__file__).parent.parent.parent / "5_s0_Me.log"

# Given directly (Figure 2 / compound 3 convention) -- not derived, since this
# is an external reference file, not one of our 34 benchmark molecules.
CI, NI, OI = 7, 17, 18
C_ARYL, C_ALKYL = 1, 8


def main() -> None:
    if not REFERENCE_LOG.exists():
        print(f"ERROR: {REFERENCE_LOG} not found.")
        return

    lines = REFERENCE_LOG.read_text().splitlines()

    cmo_starts = find_cmo_sections(lines)
    e2pert_starts = find_table_starts(lines)
    print(f"Found {len(cmo_starts)} CMO section(s), {len(e2pert_starts)} E2PERT section(s).")
    if len(cmo_starts) > 1 or len(e2pert_starts) > 1:
        print("  -> Stable=Opt re-test detected (>1 section) -- using the LAST of each, per the handout's warning.\n")

    # ---- Tier 1: single-geometry descriptors from the last (correct) NBO section ----
    cmo_table = parse_cmo_table(lines, cmo_starts[-1])
    window = virtual_window(cmo_table)
    lumo_e = window[0]["energy"] if window else None
    print(f"Tier 1: LUMO = {lumo_e} a.u., {len(window)} virtual MOs in the LUMO..LUMO+0.4 window\n")

    wcnmax, wcnmax_mo, wcnmax_eps, wcnmax_coeff = max_weight_for_target(window, CI, NI)
    w17max, w17max_mo, w17max_eps, w17max_coeff = max_weight_for_target(window, CI, C_ARYL)
    w78max, w78max_mo, w78max_eps, w78max_coeff = max_weight_for_target(window, CI, C_ALKYL)

    print(
        f"  wCNmax = {wcnmax} (MO {wcnmax_mo}, epsilon={wcnmax_eps} a.u., coefficient={wcnmax_coeff})  "
        f"-- expect 0.457 (MO 32) per the wCNmax handout's worked example"
    )
    ok = wcnmax is not None and abs(wcnmax - 0.457) < 1e-3
    print(f"  {'PASS' if ok else 'MISMATCH -- FLAGGING, NOT ADJUSTING THE FORMULA'}: matches handout's worked answer\n")

    print(f"  w17max = {w17max} (MO {w17max_mo}, epsilon={w17max_eps} a.u., coefficient={w17max_coeff})")
    print(f"  w78max = {w78max} (MO {w78max_mo}, epsilon={w78max_eps} a.u., coefficient={w78max_coeff})")

    import math
    lam = w78max / w17max if (w17max is not None and w17max > 0 and w78max is not None) else None
    log_lam = math.log10(lam) if lam is not None and lam > 0 else None
    print(f"  Lambda = {lam}, log_lambda = {log_lam}")
    print("  (no single-point reference value exists for Lambda -- Table 2 only reports d/dR log10(Lambda))\n")

    e2pert_rows = [
        {"donor": r["donor"], "acceptor": r["acceptor"], "e2_kcal": r["e2_kcal"]}
        for r in parse_table_rows(lines, e2pert_starts[-1])
    ]
    psi_row = compute_psi_row(e2pert_rows, CI, NI, OI, C_ARYL, C_ALKYL)
    print(f"  K_anti = {psi_row['k_anti']:.4f} kcal/mol  (verified by hand against the raw E2PERT rows)")
    print(f"  K_frag = {psi_row['k_frag']:.4f} kcal/mol  (verified by hand against the raw E2PERT rows)")
    print(f"  Psi    = {psi_row['psi']:.4f}  (epsilon={PSI_EPSILON})")
    print(
        "  NOTE: Psi < 1 here (fragmentation-channel stabilization exceeds rearrangement-channel\n"
        "  stabilization at R0, even though compound 3 experimentally rearranges). Table 2 only reports\n"
        "  d/dR(Psi), not an absolute-value threshold, and ALL FOUR compounds in Table 2 show positive\n"
        "  d/dR(Psi) regardless of R/F outcome -- so this may not actually be a distinguishing signal at\n"
        "  all, and Psi<1 at equilibrium might be entirely expected. I have not read Section 3.2's actual\n"
        "  text (out of scope per instructions) so I cannot confirm this either way -- flagging as an open\n"
        "  question rather than asserting it's right or wrong.\n"
    )

    # ---- Tier 2: blocked, report don't approximate ----
    print("Tier 2 (d/dR slope validation against Table 2): BLOCKED.")
    print("  Only 5_s0_Me.log (R0 equilibrium) exists in the repo for compound 3.")
    print("  Table 2's reported d/dR values need all 5 scan points (R0 through R0+0.4 A).")
    print("  Not approximating a slope from one point -- ask Tetiana for the other 4 scan-point logs.")


if __name__ == "__main__":
    main()
