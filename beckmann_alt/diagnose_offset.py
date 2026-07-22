"""
Diagnostic 1 of 2 (per the current investigation plan, see Notes_open_source_alt.md's
"wCNmax: PySCF vs NBO7" section) -- root-causing the fixed +5.2-6.2%/-1-MO-index offset
BEFORE touching pair_nbo.py's actual construction. One-off script, not part of the
regular wcnmax_scan_rule.py pipeline. Writes results to
data/output/analysis/wcnmax_offset_diagnostic.csv (new file -- does not touch
wcnmax_rule_results_opensource.csv or any file under data/output/dft_opt/).

Probes, in order:
  0. Sanity check: does PySCF's winning canonical MO for mol_002_E's 'nbo' point
     actually land on ITS OWN LUMO (mo_index == nocc), the same way NBO7's trusted
     MO48 lands on Gaussian's own LUMO (delta_lumo=0.0 in cmo_channel_extraction.csv,
     true for all 6 test molecules)? If so, the reported "MO index exactly 1 below"
     pattern is a 0-based (PySCF array index) vs 1-based (Gaussian/NBO7 orbital
     number) labeling mismatch, not a shifted virtual-orbital ordering.
  1. Same solvated Stage-1 geometry (the one everything else in this project uses),
     solvent removed entirely from the PySCF SCF (solvent_name=None) -- isolates the
     ddCOSMO-vs-(SMD or nothing) variable while holding geometry fixed.
  2. Same solvated geometry, density fitting turned off -- isolates the DF variable.
  3. PySCF gas-phase SCF on the REAL gas-phase-optimized Gaussian geometry (from
     data/output/dft_opt_gasphase_archive/mol_002_E/, a genuinely different, fully
     self-consistent gas-phase opt+NBO pair) -- lets the raw frontier virtual
     eigenvalues be compared directly against that archive's own gas-phase
     Alpha virt. eigenvalues (already confirmed to exist and match a real Gaussian
     gas-phase SP on that exact geometry).
"""
import csv

from beckmann.config import DATA_OUTPUT
from beckmann_alt.geometry import CHARGE, SPIN, final_geometry, load_case, pyscf_atom_spec
from beckmann_alt.pair_nbo import build_local_iaos, compute_wcnmax
from beckmann_alt.pyscf_livvo import build_mol, run_scf

OUT_CSV = DATA_OUTPUT / "analysis" / "wcnmax_offset_diagnostic.csv"
GAS_ARCHIVE_OPT_LOG = (
    DATA_OUTPUT / "dft_opt_gasphase_archive" / "mol_002_E" / "mol_002_E_opt.log"
)

# From data/output/dft_opt/mol_002_E/mol_002_E_nbo.log directly (grepped, not
# recomputed): 47 occupied orbitals, first ("Alpha virt.") eigenvalue -0.00029 a.u.
# = NBO7's own trusted MO48 (cmo_channel_extraction.csv: epsilon_i_star=-0.00029,
# delta_lumo=0.0) -- i.e. NBO7's winning MO IS its own LUMO, 1-based index 48.
TRUSTED_NBO7 = {"wcnmax": 0.4356, "mo_index_1based": 48, "epsilon": -0.00029, "nocc_1based_last_occ": 47}

# From the SAME log's gas-phase-archive twin, Alpha virt. eigenvalues (first 5):
# -0.14455 -0.13600 -0.08720 -0.07584 -0.03976 (39 occ orbitals, first virt at line 582)
GAS_GAUSSIAN_FRONTIER_VIRT = [-0.14455, -0.13600, -0.08720, -0.07584, -0.03976, -0.02810, -0.02613]


def run_one(label: str, case: dict, solvent_name, density_fit: bool) -> dict:
    mol = build_mol(case)
    mf = run_scf(mol, solvent_name=solvent_name, density_fit=density_fit)
    iaos_orth, s, atom_of_iao = build_local_iaos(mf)
    cn = compute_wcnmax(mf, s, iaos_orth, atom_of_iao, case["ci"], case["ni"])
    nocc = mol.nelectron // 2
    lumo_e = mf.mo_energy[nocc]
    mo_index = cn["mo_index"]
    return {
        "condition": label,
        "solvent": "ddCOSMO/water" if solvent_name else "gas (none)",
        "density_fit": density_fit,
        "wCNmax": round(cn["wmax"], 4),
        "MO_index_0based": mo_index,
        "MO_index_1based": mo_index + 1,
        "nocc": nocc,
        "epsilon": round(cn["epsilon"], 5),
        "lumo_epsilon": round(lumo_e, 5),
        "is_own_lumo": mo_index == nocc,
        "pct_diff_vs_NBO7": round(100 * (cn["wmax"] - TRUSTED_NBO7["wcnmax"]) / TRUSTED_NBO7["wcnmax"], 2),
        "frontier_virt_eigs_5": [round(e, 5) for e in mf.mo_energy[nocc:nocc + 5]],
    }


def main() -> None:
    solvated_case = load_case("mol_002")
    results = []

    print("[1/3] baseline: solvated Stage-1 geometry, ddCOSMO/water + density-fit (current default config)")
    results.append(run_one("baseline_ddcosmo_df", solvated_case, "water", True))

    print("[2/3] same geometry, solvent removed (gas phase), density-fit on")
    results.append(run_one("gas_same_geometry_df", solvated_case, None, True))

    print("[3/3] same geometry, ddCOSMO/water restored, density-fit OFF (slow)")
    results.append(run_one("ddcosmo_no_df", solvated_case, "water", False))

    print("[extra] PySCF gas-phase SCF on the REAL gas-phase-optimized Gaussian geometry")
    gas_atoms = final_geometry(GAS_ARCHIVE_OPT_LOG)
    gas_case = {
        **solvated_case,
        "name": "mol_002_gasphase_archive_geom",
        "atoms": gas_atoms,
        "atom_spec": pyscf_atom_spec(gas_atoms),
        "charge": CHARGE, "spin": SPIN,
    }
    gas_result = run_one("gas_archive_geometry_df", gas_case, None, True)
    gas_result["gaussian_gas_frontier_virt_5"] = GAS_GAUSSIAN_FRONTIER_VIRT[:5]
    results.append(gas_result)

    for r in results:
        print(r)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({k for r in results for k in r})
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    print(f"\nWrote {OUT_CSV}")


if __name__ == "__main__":
    main()
