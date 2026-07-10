"""
Build product/intermediate geometries for the two competing TS channels
(rearrangement and stepwise fragmentation), for use as QST2/QST3 (Gaussian) and
NEB (AIMNet2/PySisyphus) endpoints.

Every stationary point in this mechanism -- reactant, rearrangement product,
fragmentation intermediate, fragmentation product -- has the SAME atom count and
formula as the starting protonated oxime (nothing dissociates to infinity in this
minimal model; see Notes.md). That means atom-index correspondence between
reactant and product, which QST2/QST3 and NEB both require, can be made exact by
construction: every builder here edits a COPY of the reactant's own RDKit Mol
(bonds/formal charges only, via RWMol -- same in-place-edit style as
beckmann.oximes.ketone_to_protonated_oximes) and displaces its existing 3D
coordinates, rather than building a fresh structure from SMILES/a reaction SMARTS
(which would not preserve atom order).

The RWMol bond edits below exist to document the mechanistic bond changes and to
produce a sensibly-bonded output Mol/SDF -- they are not required to pass strict
RDKit valence sanitization, since nothing downstream (AIMNet2, Gaussian) reads
RDKit bond orders. AIMNet2 sees only atomic numbers + total charge + coordinates;
Gaussian sees only atomic numbers + charge/multiplicity + coordinates.

Output: data/output/aimnet_optimized/{mol}_product_{rearr,frag_int,frag}.sdf
"""
import os
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

from pathlib import Path

import numpy as np
from rdkit import Chem
from ase import Atoms

from beckmann.config import DATA_OUTPUT, CHARGE
from beckmann.dft.descriptors import get_substituent_map, _load_mols
from beckmann.optimize import relax_geometry
from aimnet.calculators import AIMNet2Calculator

# A single soft nudge is not enough to cross the real barrier (~31 kcal/mol per
# the paper's TS1_A1) -- local optimization (LBFGS) just rolls back downhill into
# the reactant's own minimum (verified empirically: an early version of this
# module used a fractional nudge with no bias and relaxed to within 0.0001 eV of
# the unperturbed reactant energy). Instead: place the forming/breaking bond at a
# value that is topologically incompatible with the reactant (a real
# bonding/non-bonding distance), bias the relaxation toward holding it there with
# a smooth harmonic restraint (beckmann.optimize.HarmonicBondRestraint, via
# relax_geometry(..., restraints=...)) through a first relaxation (so the
# optimizer cannot slide back across the ridge), then relax again with no
# restraint for a short free "polish" from what is by then a genuinely different
# basin. (ASE's hard FixBondLength constraint was tried first and rejected: its
# iterative RATTLE-style solver failed to converge when two constraints share an
# atom, exactly the case here since both restraints involve N(ni).)
N_C_BOND_TARGET  = 1.55  # Angstrom, target N(ni)-C(migrating) distance (single C-N bond)
N_O_LEAVE_TARGET = 2.80  # Angstrom, target N(ni)-O(oi) distance (clearly non-bonding)
C_C_BREAK_TARGET = 2.60  # Angstrom, target C(ci)-C(alkyl) distance (clearly non-bonding)
RESTRAINT_K       = 15.0  # eV/Angstrom^2, comparable order-of-magnitude to a real bond force constant


def _mol_atoms(mol: Chem.Mol) -> Atoms:
    """RDKit Mol (with a conformer) -> ASE Atoms, same convention as optimize.py."""
    conf    = mol.GetConformer()
    numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    return Atoms(numbers=numbers, positions=conf.GetPositions())


def _hs_bonded_to(atoms: Atoms, center_idx: int, cutoff: float = 1.3) -> list[int]:
    """Indices of H atoms within `cutoff` Angstrom of atoms[center_idx].

    Geometric, not RDKit-bond-table-based, so it doesn't depend on the edited
    RWMol's bonds having survived sanitization.
    """
    pos = atoms.get_positions()
    center = pos[center_idx]
    return [
        i for i, num in enumerate(atoms.numbers)
        if num == 1 and np.linalg.norm(pos[i] - center) < cutoff
    ]


def _place_at_distance(pos: np.ndarray, anchor_idx: int, move_idx: int, target: float) -> None:
    """In place: move pos[move_idx] to exactly `target` Angstrom from pos[anchor_idx],
    along their current direction."""
    direction = pos[move_idx] - pos[anchor_idx]
    direction = direction / np.linalg.norm(direction)
    pos[move_idx] = pos[anchor_idx] + target * direction


def displace_for_rearrangement(atoms: Atoms, atom_map: dict) -> list[tuple[int, int, float, float]]:
    """Place N(ni)-C(c_aryl) at a bonding distance and N(ni)-O(oi) at a clearly
    non-bonding distance (dragging O's bonded H's along rigidly). Mutates `atoms`
    in place and returns the restraint list for the first-stage relax_geometry()
    call (see module docstring for why a restraint, not a hard constraint).
    """
    pos = atoms.get_positions().copy()
    ni, oi = atom_map["ni"] - 1, atom_map["oi"] - 1
    c_aryl = atom_map["c_aryl"] - 1

    _place_at_distance(pos, ni, c_aryl, N_C_BOND_TARGET)

    o_group = [oi] + _hs_bonded_to(atoms, oi)
    n_o_dir = pos[oi] - pos[ni]
    n_o_dir = n_o_dir / np.linalg.norm(n_o_dir)
    delta = (N_O_LEAVE_TARGET - np.linalg.norm(pos[oi] - pos[ni])) * n_o_dir
    for idx in o_group:
        pos[idx] = pos[idx] + delta

    atoms.set_positions(pos)
    return [(ni, c_aryl, N_C_BOND_TARGET, RESTRAINT_K), (ni, oi, N_O_LEAVE_TARGET, RESTRAINT_K)]


def displace_for_fragmentation(atoms: Atoms, atom_map: dict) -> list[tuple[int, int, float, float]]:
    """Place C(ci)-C(alkyl) and N(ni)-O(oi) at clearly non-bonding distances
    (ring-opening + leaving-group departure). Same mutate-in-place/return-restraints
    shape as displace_for_rearrangement, for the same reason.
    """
    pos = atoms.get_positions().copy()
    ci, ni, oi = atom_map["ci"] - 1, atom_map["ni"] - 1, atom_map["oi"] - 1
    c_alkyl = atom_map["c_alkyl"] - 1

    _place_at_distance(pos, ci, c_alkyl, C_C_BREAK_TARGET)

    o_group = [oi] + _hs_bonded_to(atoms, oi)
    n_o_dir = pos[oi] - pos[ni]
    n_o_dir = n_o_dir / np.linalg.norm(n_o_dir)
    delta = (N_O_LEAVE_TARGET - np.linalg.norm(pos[oi] - pos[ni])) * n_o_dir
    for idx in o_group:
        pos[idx] = pos[idx] + delta

    atoms.set_positions(pos)
    return [(ci, c_alkyl, C_C_BREAK_TARGET, RESTRAINT_K), (ni, oi, N_O_LEAVE_TARGET, RESTRAINT_K)]


def build_rearrangement_product(mol: Chem.Mol, atom_map: dict) -> Chem.Mol:
    """Nitrilium ion from migration of c_aryl: C(ci)-C(aryl) breaks, N(ni)-C(aryl)
    forms, C(ci)=N(ni) becomes a formal triple bond, the +1 charge moves from
    O(oi) to N(ni). O(oi)/its H's are left bond-free (untouched water fragment).
    Connectivity only -- does not touch the conformer/coordinates.
    """
    ci, ni, oi = atom_map["ci"] - 1, atom_map["ni"] - 1, atom_map["oi"] - 1
    c_aryl = atom_map["c_aryl"] - 1

    rw = Chem.RWMol(mol)
    rw.RemoveBond(ci, c_aryl)
    rw.AddBond(ni, c_aryl, Chem.BondType.SINGLE)
    rw.GetBondBetweenAtoms(ci, ni).SetBondType(Chem.BondType.TRIPLE)
    rw.RemoveBond(ni, oi)
    rw.GetAtomWithIdx(oi).SetFormalCharge(0)
    rw.GetAtomWithIdx(ni).SetFormalCharge(1)

    product = rw.GetMol()
    Chem.SanitizeMol(product, catchErrors=True)
    return product


def build_fragmentation_intermediate(mol: Chem.Mol, atom_map: dict) -> Chem.Mol:
    """Ring-opened nitrilium from cleavage of C(ci)-C(alkyl): that bond breaks,
    C(ci)=N(ni) becomes a formal triple bond, the +1 charge moves from O(oi) to
    N(ni). The alkyl fragment stays connected via whatever other ring bonds
    remain (ring-opening, not fragment separation). Connectivity only.
    """
    ci, ni, oi = atom_map["ci"] - 1, atom_map["ni"] - 1, atom_map["oi"] - 1
    c_alkyl = atom_map["c_alkyl"] - 1

    rw = Chem.RWMol(mol)
    rw.RemoveBond(ci, c_alkyl)
    rw.GetBondBetweenAtoms(ci, ni).SetBondType(Chem.BondType.TRIPLE)
    rw.RemoveBond(ni, oi)
    rw.GetAtomWithIdx(oi).SetFormalCharge(0)
    rw.GetAtomWithIdx(ni).SetFormalCharge(1)

    intermediate = rw.GetMol()
    Chem.SanitizeMol(intermediate, catchErrors=True)
    return intermediate


def build_fragmentation_product(intermediate: Chem.Mol, atom_map: dict) -> Chem.Mol:
    """P_B1: from the fragmentation intermediate, form a new N(ni)-C(alkyl) bond
    (TS2_B1's step), re-closing to a different connectivity. Connectivity only.
    """
    ni, c_alkyl = atom_map["ni"] - 1, atom_map["c_alkyl"] - 1

    rw = Chem.RWMol(intermediate)
    rw.AddBond(ni, c_alkyl, Chem.BondType.SINGLE)

    product = rw.GetMol()
    Chem.SanitizeMol(product, catchErrors=True)
    return product


def _write_sdf(mol_with_geometry: Chem.Mol, name: str, out_path: Path) -> Path:
    mol_with_geometry.SetProp("_Name", name)
    writer = Chem.SDWriter(str(out_path))
    writer.write(mol_with_geometry)
    writer.close()
    return out_path


def generate_rearrangement_product(mol: str, mol_dir: Path, output_dir: Path,
                                    base_calc: AIMNet2Calculator | None = None) -> Path:
    """Build, displace, and AIMNet2-relax the rearrangement product for `mol`
    (e.g. 'mol_002_E'). Writes {output_dir}/{mol}_product_rearr.sdf, same atom
    count/order as the reactant in best_per_substrate.sdf.
    """
    mols = _load_mols()
    if mol not in mols:
        raise ValueError(f"{mol}: not found in best_per_substrate.sdf")
    reactant = mols[mol]
    atom_map = get_substituent_map(mol, mol_dir)

    product_mol = build_rearrangement_product(reactant, atom_map)

    atoms = _mol_atoms(reactant)
    restraints = displace_for_rearrangement(atoms, atom_map)
    atoms, _ = relax_geometry(atoms, charge=CHARGE, base_calc=base_calc, restraints=restraints)
    atoms, energy_ev = relax_geometry(atoms, charge=CHARGE, base_calc=base_calc)  # free polish

    conf = product_mol.GetConformer()
    for i, pos in enumerate(atoms.get_positions()):
        conf.SetAtomPosition(i, pos.tolist())
    product_mol.SetProp("E_aimnet2_eV", f"{energy_ev:.6f}")

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / f"{mol}_product_rearr.sdf"
    return _write_sdf(product_mol, f"{mol}_product_rearr", out_path)


def main() -> None:
    """Pilot scope: mol_002_E rearrangement product only."""
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    output_dir  = DATA_OUTPUT / "aimnet_optimized"

    mol = "mol_002_E"
    mol_dir = dft_opt_dir / mol
    out_path = generate_rearrangement_product(mol, mol_dir, output_dir)
    print(f"-- {mol}: rearrangement product -> {out_path}")


if __name__ == "__main__":
    main()
