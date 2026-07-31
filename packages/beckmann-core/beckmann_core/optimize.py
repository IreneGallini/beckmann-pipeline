"""
Select lowest-energy conformer per molecule, run AIMNet2/ASE optimization.

select_and_optimize() writes two output files:
  best_aimnet_optimized.sdf  -- one optimized structure per isomer name (mol_XXX_E / mol_XXX_Z)
  best_per_substrate.sdf     -- one structure per substrate (lowest-energy isomer wins)
"""
import os
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

from pathlib import Path

import numpy as np
from rdkit import Chem

from aimnet.calculators import AIMNet2Calculator, AIMNet2ASE
from ase.calculators.calculator import Calculator as ASECalculator, all_changes
from ase.calculators.mixing import SumCalculator
from ase.optimize import LBFGS
from ase import Atoms


class HarmonicBondRestraint(ASECalculator):
    """Smooth harmonic bias E = 0.5*k*(r-r0)^2 on a set of atom-pair distances.

    Used to bias relax_geometry() toward a specific bond length (e.g. a forming
    or breaking bond) without a hard geometric constraint -- ASE's FixBondLength
    uses an iterative RATTLE-style solver that can fail to converge when two
    constraints share an atom (seen in practice restraining both a forming and a
    leaving bond at the same central atom). A smooth energy bias has no such
    convergence failure mode; it's just additional forces stacked via SumCalculator.
    """
    implemented_properties = ["energy", "forces"]

    def __init__(self, restraints: list[tuple[int, int, float, float]]):
        """restraints: list of (atom_i, atom_j, target_distance_angstrom, k_eV_per_A2)."""
        super().__init__()
        self.restraints = restraints

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        pos = atoms.get_positions()
        energy = 0.0
        forces = np.zeros_like(pos)
        for i, j, r0, k in self.restraints:
            vec = pos[j] - pos[i]
            r = np.linalg.norm(vec)
            dr = r - r0
            energy += 0.5 * k * dr ** 2
            direction = vec / r
            f = (-k * dr) * direction
            forces[j] += f
            forces[i] -= f
        self.results["energy"] = energy
        self.results["forces"] = forces


def relax_geometry(
    atoms: Atoms,
    charge: int,
    model: str = "aimnet2_2025",
    base_calc: AIMNet2Calculator | None = None,
    fmax: float = 0.05,
    restraints: list[tuple[int, int, float, float]] | None = None,
) -> tuple[Atoms, float]:
    """Relax an ASE Atoms object in place with AIMNet2/LBFGS.

    Returns (relaxed_atoms, energy_ev). Pass a pre-built base_calc to avoid reloading
    model weights across repeated calls. Pass `restraints` (atom_i, atom_j,
    target_distance, k) pairs to bias specific bond distances toward a target during
    this relaxation -- e.g. to hold a forming/breaking bond away from the reactant's
    own value so the optimizer can't just roll back downhill into it.
    energy_ev is the AIMNet2 energy alone (restraint bias excluded).
    """
    if base_calc is None:
        base_calc = AIMNet2Calculator(model)
    aimnet_calc = AIMNet2ASE(base_calc, charge=charge)
    if restraints:
        atoms.calc = SumCalculator([aimnet_calc, HarmonicBondRestraint(restraints)])
    else:
        atoms.calc = aimnet_calc
    opt = LBFGS(atoms, logfile=None)
    opt.run(fmax=fmax)
    energy_ev = aimnet_calc.results["energy"]
    return atoms, energy_ev


def select_and_optimize(
    conformers_sdf: Path,
    output_dir: Path,
) -> tuple[Path, Path]:
    """Pick best conformer per isomer, run AIMNet2/ASE optimization.

    Returns (best_aimnet_optimized.sdf, best_per_substrate.sdf).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading conformers from: {conformers_sdf}")
    suppl = Chem.SDMolSupplier(str(conformers_sdf), removeHs=False)
    molecules: dict[str, list[tuple[float, Chem.Mol]]] = {}

    for mol in suppl:
        if mol is None:
            continue
        name   = mol.GetProp("_Name") if mol.HasProp("_Name") else "unknown"
        energy = float(mol.GetProp("E_tot")) if mol.HasProp("E_tot") else float("inf")
        molecules.setdefault(name, []).append((energy, mol))

    print("Selecting lowest energy conformer per molecule:")
    best_conformers: dict[str, Chem.Mol] = {}
    for name, confs in molecules.items():
        confs.sort(key=lambda x: x[0])
        best_energy, best_mol = confs[0]
        best_conformers[name] = best_mol
        print(f"  {name}: selected E = {best_energy:.6f} Hartree "
              f"({best_energy * 627.509:.2f} kcal/mol) from {len(confs)} candidates")

    print(f"\nRunning AIMNet2 optimization on {len(best_conformers)} structure(s)...")
    base_calc = AIMNet2Calculator("aimnet2_2025")

    best_sdf_path = output_dir / "best_aimnet_optimized.sdf"
    writer = Chem.SDWriter(str(best_sdf_path))
    optimized: dict[str, tuple[float, Chem.Mol]] = {}

    for name, mol in best_conformers.items():
        print(f"  Optimizing {name}...")
        try:
            conf    = mol.GetConformer()
            coords  = conf.GetPositions()
            numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
            atoms   = Atoms(numbers=numbers, positions=coords)

            mol_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
            atoms, energy_ev = relax_geometry(atoms, charge=mol_charge, base_calc=base_calc)
            energy_kcal = energy_ev * 23.0605

            new_conf = mol.GetConformer()
            for i, pos in enumerate(atoms.get_positions()):
                new_conf.SetAtomPosition(i, pos.tolist())

            mol.SetProp("_Name",          name)
            mol.SetProp("E_aimnet2_eV",   f"{energy_ev:.6f}")
            mol.SetProp("E_aimnet2_kcal", f"{energy_kcal:.4f}")
            writer.write(mol)
            optimized[name] = (energy_ev, mol)
            print(f"    AIMNet2 E = {energy_ev:.6f} eV ({energy_kcal:.2f} kcal/mol)")

        except Exception as e:
            print(f"    WARNING: AIMNet2 optimization failed for {name}: {e}")

    writer.close()

    substrate_best: dict[str, tuple[float, Chem.Mol]] = {}
    for name, (energy_ev, mol) in optimized.items():
        base = name.rsplit("_", 1)[0]
        if base not in substrate_best or energy_ev < substrate_best[base][0]:
            substrate_best[base] = (energy_ev, mol)

    sub_sdf_path = output_dir / "best_per_substrate.sdf"
    writer2 = Chem.SDWriter(str(sub_sdf_path))
    for base in sorted(substrate_best):
        _, mol = substrate_best[base]
        writer2.write(mol)
    writer2.close()

    print(f"\nDone. Optimized structures saved to:")
    print(f"  {best_sdf_path}  ({len(optimized)} isomers)")
    print(f"  {sub_sdf_path}  ({len(substrate_best)} substrates, lowest-energy isomer only)")
    return best_sdf_path, sub_sdf_path
