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
    starting geometry from a single fixed base structure.
  - beckmann.optimize.relax_geometry(..., restraints=[(ni, oi, target_R, k)])
    -- the harmonic-bond-restraint constrained-relaxation pattern already
    validated in beckmann.dft.ts_products for holding one bond near a target
    length while everything else relaxes. This is the AIMNet2/ASE analog of
    the Gaussian rigid-scan architecture's opt=(ModRedundant) frozen-bond
    step. RESTRAINT_K reuses ts_products.py's own value (15.0 eV/A^2,
    "comparable order-of-magnitude to a real bond force constant") rather
    than an invented number.

beckmann_alt.pair_nbo.run_from_case() itself is imported and called
unmodified -- nothing here changes its wCNmax computation logic. Its case
dict shape (atom_spec/charge/spin/ci/ni/c_aryl/name/basis_note) already
doesn't reference TEST_IDS/ALL_IDS anywhere in the function bodies, so no
interface change was needed to point it at a brand-new molecule instead of a
benchmark one.
"""
from aimnet.calculators import AIMNet2Calculator
from ase import Atoms
from rdkit import Chem

from beckmann_core.geometry import displace_leaving_group, no_distance
from beckmann_core.optimize import relax_geometry
from beckmann_pyscf.engine.pair_nbo import run_from_case

STEP = 0.05
N_POINTS = 6  # matches beckmann.dft.inputs._scan_gjf_rigid's current default
              # (R0+0.05 .. R0+0.30 A) for comparability with the trusted series
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
    step: float = STEP, n_points: int = N_POINTS,
) -> list[dict]:
    """R0 point ('nbo' stage) + n_points rigid-scan points ('scan_1'..
    'scan_N'), each a fresh PySCF wCNmax single-point
    (beckmann_pyscf.engine.pair_nbo.run_from_case, unmodified) evaluated on
    its own AIMNet2-relaxed geometry. ci/ni/oi/c_aryl/c_alkyl are 1-based
    (RDKit index + 1), typically from
    beckmann_core.classical.get_oxime_atoms(mol) + 1.

    Returns rows ready for beckmann.dft.descriptors.find_wcnmax_minimum().
    """
    charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    spin = 0  # closed-shell singlet -- matches beckmann_alt.geometry's SPIN (MULTIPLICITY=1) convention
    numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    base_atoms = _mol_to_atom_tuples(mol)
    r0 = no_distance(base_atoms, ni, oi)

    rows = []
    r0_case = _case_from_atoms(base_atoms, ci, ni, c_aryl, c_alkyl, charge, spin, name, "nbo", r0)
    rows.append(_row_from_case_result(r0_case, run_from_case(r0_case)))

    base_calc = AIMNet2Calculator("aimnet2_2025")
    for pt in range(1, n_points + 1):
        delta = pt * step
        displaced = displace_leaving_group(base_atoms, ni, oi, delta)
        target_r = r0 + delta

        atoms_obj = Atoms(numbers=numbers, positions=[(x, y, z) for _, x, y, z in displaced])
        # HarmonicBondRestraint indexes directly into atoms.get_positions(),
        # a plain 0-based numpy array -- ni/oi are 1-based everywhere else in
        # this module (matching displace_leaving_group/no_distance), so
        # convert here only, same as beckmann.dft.ts_products does at its
        # own restraints-list call sites.
        atoms_obj, _ = relax_geometry(
            atoms_obj, charge=charge, base_calc=base_calc,
            restraints=[(ni - 1, oi - 1, target_r, RESTRAINT_K)],
        )
        relaxed_tuples = [
            (sym, *pos) for (sym, *_), pos in zip(displaced, atoms_obj.get_positions())
        ]
        actual_r = no_distance(relaxed_tuples, ni, oi)

        case = _case_from_atoms(
            relaxed_tuples, ci, ni, c_aryl, c_alkyl, charge, spin, name, f"scan_{pt}", actual_r,
        )
        rows.append(_row_from_case_result(case, run_from_case(case)))

    return rows
