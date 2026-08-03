"""
Build a PySCF-native R(N-O) scan series for one AIMNet2-optimized molecule --
the gap nothing in beckmann_alt/geometry.py covers, since its loaders all read
pre-existing Gaussian .gjf/.log files keyed to the fixed 34-molecule benchmark
set (see Notes_pyscf_alt.md). No Gaussian/NBO7/Citadel step anywhere here.

Two already-validated primitives (vendored unchanged) do the actual geometry
work -- this module is assembly, not new science:

  - beckmann.dft.geometry.displace_leaving_group() -- the same rigid
    pre-displacement the Gaussian rigid-scan architecture
    (beckmann.dft.inputs._scan_gjf_rigid) uses to build each scan point's
    starting geometry from a single fixed base structure. Handles negative
    delta (compression) the same as positive (stretching) -- pure signed
    vector math along the N->O unit vector.
  - beckmann.optimize.relax_geometry(..., restraints=[(ni, oi, target_R, k)])
    -- the harmonic-bond-restraint constrained-relaxation pattern already
    validated in beckmann.dft.ts_products for holding one bond near a target
    length while everything else relaxes. This is the AIMNet2/ASE analog of
    the Gaussian rigid-scan architecture's opt=(ModRedundant) frozen-bond
    step. RESTRAINT_K reuses ts_products.py's own value (15.0 eV/A^2,
    "comparable order-of-magnitude to a real bond force constant") rather
    than an invented number.

Scan window is a FIXED ABSOLUTE R(N-O) range (FIXED_R_MIN..FIXED_R_MAX),
the same for every molecule -- not anchored to each molecule's own AIMNet2-
optimized equilibrium R0. AIMNet2's unconstrained R0 is systematically
longer than the corresponding DFT-optimized R0 (confirmed directly: mol_006
AIMNet2 R0=1.664 A vs. Gaussian R0=1.511 A, a ~0.15 A gap), while NBO7's
trusted 34-molecule data shows the actual wCNmax-minimum crossing clustering
tightly at 1.594-1.703 A (mean 1.659, std 0.025) regardless of molecule. An
R0-anchored, outward-only scan starts already inside or past that window for
molecules like mol_006 and can never bracket the crossing as an interior
point. The fixed window here (1.50-1.80 A) pads ~0.09-0.10 A beyond the
observed NBO7 range on each side so the crossing lands as a genuine interior
scan point regardless of where any given molecule's own AIMNet2 R0 falls.

beckmann_alt.pair_nbo.run_from_case() itself is imported and called
unmodified -- nothing here changes its wCNmax computation logic. Its case
dict shape (atom_spec/charge/spin/ci/ni/c_aryl/name/basis_note) already
doesn't reference TEST_IDS/ALL_IDS anywhere in the function bodies, so no
interface change was needed to point it at a brand-new molecule instead of a
benchmark one.

Every run_from_case() call below goes through pyscf_subprocess.run_pyscf_isolated()
instead of calling it directly -- PyTorch (AIMNet2, used earlier in this same
loop) and PySCF's own BLAS/OpenMP runtime conflict when both are resident in
one process, confirmed to segfault PySCF's SCF reproducibly. See
pyscf_subprocess.py's module docstring for the full story.
"""
import time

from aimnet.calculators import AIMNet2Calculator
from ase import Atoms
from rdkit import Chem

from beckmann_core.geometry import displace_leaving_group, no_distance
from beckmann_core.optimize import relax_geometry
from beckmann_pyscf.pipeline.pyscf_subprocess import run_pyscf_isolated

FIXED_R_MIN = 1.50  # Angstrom -- see module docstring for why this is fixed/absolute, not R0-relative
FIXED_R_MAX = 1.80
FIXED_R_STEP = 0.05  # same point density as the original R0-relative scan (7 points total)
RESTRAINT_K = 15.0  # eV/A^2 -- beckmann.dft.ts_products.RESTRAINT_K's own value


def _mol_to_atom_tuples(mol: Chem.Mol) -> list[tuple]:
    """(sym, x, y, z) tuples, matching beckmann.dft.geometry's 1-based
    (atoms[i-1]) indexing convention when combined with RDKit's 0-based
    atom indices + 1."""
    conf = mol.GetConformer()
    return [
        (atom.GetSymbol(), *conf.GetAtomPosition(i))
        for i, atom in enumerate(mol.GetAtoms())
    ]


def _case_from_atoms(
    atoms: list[tuple], ci: int, ni: int, c_aryl: int, c_alkyl: int,
    charge: int, spin: int, name: str, stage: str, r_no: float,
) -> dict:
    """Same case-dict shape beckmann_alt.geometry's loaders build (atom_spec/
    charge/spin/ci/ni/c_aryl/c_alkyl/name/basis_note), sourced from a fresh
    AIMNet2 geometry instead of a parsed Gaussian log."""
    return {
        "name": name, "stage": stage, "r_no": r_no,
        "atom_spec": [[sym, (x, y, z)] for sym, x, y, z in atoms],
        "charge": charge, "spin": spin,
        "ci": ci, "ni": ni, "c_aryl": c_aryl, "c_alkyl": c_alkyl,
        "basis_note": "AIMNet2-optimized geometry, no Gaussian/DFT step -- see beckmann-pyscf README",
    }


def _row_from_case_result(case: dict, result: dict) -> dict:
    """Same row shape beckmann_alt.wcnmax_scan_rule's own PySCF adapter
    produces from run_test_set_scan_series() -- exactly what
    beckmann.dft.descriptors.find_wcnmax_minimum() expects
    (mol/stage/channel/R_NO/weight/MO_index/epsilon_i_star)."""
    cn = result["cn"]
    second = cn["second"]
    aryl_coeffs = result["aryl_coeffs"]
    return {
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
    }


def run_scan_series(
    mol: Chem.Mol, ci: int, ni: int, oi: int, c_aryl: int, c_alkyl: int, name: str,
    r_min: float = FIXED_R_MIN, r_max: float = FIXED_R_MAX, r_step: float = FIXED_R_STEP,
) -> list[dict]:
    """Fixed-absolute-window R(N-O) scan: one row per target in
    [r_min, r_max] stepped by r_step (7 points by default: r_min is the
    'nbo'-stage row, the rest are 'scan_1'..'scan_N'), each a fresh PySCF
    wCNmax single-point (beckmann_pyscf.engine.pair_nbo.run_from_case,
    unmodified, via pyscf_subprocess.run_pyscf_isolated) evaluated on its own
    AIMNet2-relaxed geometry, restrained toward that target R(N-O) the same
    way at every point -- see module docstring for why the window is fixed/
    absolute rather than anchored to this molecule's own AIMNet2 R0.
    ci/ni/oi/c_aryl/c_alkyl are 1-based (RDKit index + 1), typically from
    beckmann_core.classical.get_oxime_atoms(mol) + 1.

    Returns rows ready for beckmann.dft.descriptors.find_wcnmax_minimum().
    """
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    spin = 0  # closed-shell singlet -- matches beckmann_alt.geometry's SPIN (MULTIPLICITY=1) convention
    numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    base_atoms = _mol_to_atom_tuples(mol)
    r0 = no_distance(base_atoms, ni, oi)
    print(f"  [scan] AIMNet2 equilibrium R(N-O) = {r0:.4f} A "
          f"(fixed scan window {r_min:.2f}-{r_max:.2f} A, diagnostic only -- not the scan anchor)", flush=True)

    n_steps = round((r_max - r_min) / r_step)
    targets = [r_min + i * r_step for i in range(n_steps + 1)]

    rows = []
    base_calc = AIMNet2Calculator("aimnet2_2025")
    for pt, target_r in enumerate(targets, start=1):
        stage = "nbo" if pt == 1 else f"scan_{pt - 1}"
        delta = target_r - r0
        displaced = displace_leaving_group(base_atoms, ni, oi, delta)
        print(f"  [scan] point {pt}/{len(targets)}: target R(N-O) = {target_r:.4f} A, relaxing with AIMNet2...", flush=True)

        atoms_obj = Atoms(numbers=numbers, positions=[(x, y, z) for _, x, y, z in displaced])
        # HarmonicBondRestraint indexes directly into atoms.get_positions(),
        # a plain 0-based numpy array -- ni/oi are 1-based everywhere else in
        # this module (matching displace_leaving_group/no_distance), so
        # convert here only, same as beckmann.dft.ts_products does at its
        # own restraints-list call sites.
        t_relax = time.perf_counter()
        atoms_obj, _, converged = relax_geometry(
            atoms_obj, charge=charge, base_calc=base_calc,
            restraints=[(ni - 1, oi - 1, target_r, RESTRAINT_K)],
        )
        relax_elapsed = time.perf_counter() - t_relax
        relaxed_tuples = [
            (sym, *pos) for (sym, *_), pos in zip(displaced, atoms_obj.get_positions())
        ]
        actual_r = no_distance(relaxed_tuples, ni, oi)
        if not converged:
            print(f"    WARNING: relaxation did not converge within max_steps ({relax_elapsed:.1f}s)", flush=True)
        print(f"    relaxed: R(N-O) = {actual_r:.4f} A (target {target_r:.4f}, "
              f"delta {actual_r - target_r:+.4f}), {relax_elapsed:.1f}s", flush=True)

        case = _case_from_atoms(
            relaxed_tuples, ci, ni, c_aryl, c_alkyl, charge, spin, name, stage, actual_r,
        )
        t_pyscf = time.perf_counter()
        result = run_pyscf_isolated(case)
        pyscf_elapsed = time.perf_counter() - t_pyscf
        rows.append(_row_from_case_result(case, result))
        print(f"    wCNmax = {result['cn']['wmax']:.4f} at MO{result['cn']['mo_index']} "
              f"({pyscf_elapsed:.1f}s)", flush=True)

    return rows
