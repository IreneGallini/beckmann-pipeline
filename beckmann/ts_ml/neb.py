"""
Fast AIMNet2-nse/PySisyphus proxy for TS location: climbing-image NEB between
reactant/product endpoints -> rsirfo TS refinement with a Hessian -> IRC
connectivity check. Same verification shape as beckmann.dft.ts.verify_ts (one
imaginary frequency after mass-weighting, displacement matches the expected
reaction coordinate) but on the much cheaper AIMNet2-nse surface, entirely
local -- no Citadel/Gaussian involvement.

Reuses geometries from beckmann.dft.ts_products (reactant from
best_per_substrate.sdf, product/intermediate SDFs built there) -- same
guaranteed atom-order correspondence NEB needs, same reason QST2 needs it.

Known PySisyphus/AIMNet2 limitation: for a multi-TS reaction (this project's
stepwise fragmentation channel), run a SEPARATE NEB for each consecutive
stationary-point pair, not one NEB across the whole path -- mirrors
beckmann.dft.ts's own TS1_B1/TS2_B1 split.

compile_model must stay False (the default) on any AIMNet2Calculator used here
-- the Hessian step (tsopt do_hess=True) is incompatible with compile_model=True
per the AIMNet2 docs.

Output: data/output/ts_ml/{mol}_{label}/  (pysisyphus working directory: cos
trajectory, ts_final_geometry.xyz, ts_final_hessian.h5, irc.trj, etc.)
"""
import os
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

from pathlib import Path

import h5py
import numpy as np
from rdkit import Chem

from aimnet.calculators.aimnet2pysis import AIMNet2Pysis
import pysisyphus.run as pysis_run

from beckmann.config import DATA_OUTPUT, CHARGE, MULTIPLICITY

MODEL = "aimnet2-nse"

pysis_run.CALC_DICT["aimnet"] = AIMNet2Pysis


def mol_to_xyz(mol: Chem.Mol, title: str) -> str:
    """Same simple XYZ format as beckmann.analysis.classical.mol_to_xyz -- not
    imported from there, to avoid a beckmann.analysis <-> beckmann.ts_ml
    dependency for a 5-line format."""
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), title]
    for atom in mol.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():<3}  {p.x:>12.6f}  {p.y:>12.6f}  {p.z:>12.6f}")
    return "\n".join(lines) + "\n"


def build_run_dict(reactant_xyz: Path, product_xyz: Path) -> dict:
    """NEB (climbing image) -> rsirfo TS refinement -> IRC, all on aimnet2-nse."""
    return {
        "geom": {
            "type": "cart",
            "fn": [str(reactant_xyz), str(product_xyz)],
        },
        # Without this, PySisyphus treats the 2 endpoints as a already-complete,
        # fully-fixed 2-image "path" with nothing to optimize (ZeroStepLength
        # crash, seen in practice) -- interpol generates the actual moving images
        # NEB needs. idpp (not plain linear Cartesian interpolation) avoids atoms
        # passing through each other across a bond-reorganization path.
        "interpol": {
            "type": "idpp",
            "between": 10,
        },
        "calc": {
            "type": "aimnet",
            "model": MODEL,
            "charge": CHARGE,
            "mult": MULTIPLICITY,
        },
        "cos": {
            "type": "neb",
            "climb": True,
        },
        "opt": {
            "type": "lbfgs",
            "align": True,
        },
        "tsopt": {
            "type": "rsirfo",
            "do_hess": True,
            "hessian_recalc": 5,  # recompute (not just Bofill-update) every 5 cycles --
                                  # AIMNet2's Hessian is cheap; first attempt without
                                  # this left 2 spurious small extra negative modes
                                  # (-97, -9 cm^-1) alongside the real one (-510 cm^-1).
        },
        "irc": {
            "type": "eulerpc",
        },
    }


def run_ts_search(mol: str, label: str, reactant: Chem.Mol, product: Chem.Mol) -> dict:
    """Run NEB -> TS refine -> IRC for one molecule/pathway. Writes XYZ endpoints
    and the PySisyphus working directory under data/output/ts_ml/{mol}_{label}/.
    Returns the raw run_result from pysisyphus.run.run_from_dict() for inspection
    -- this module does NOT auto-decide pass/fail, same precedent as
    beckmann.dft.ts.verify_ts.

    run_from_dict's `cwd` kwarg only affects logging setup -- internally,
    PySisyphus writes trajectories/hessians/logs relative to the *process*
    working directory regardless (confirmed empirically: an early version of
    this function passed cwd=out_dir and PySisyphus still dumped ~90 files into
    the repo root). We os.chdir() into out_dir ourselves and restore the
    original directory in a finally block, even on error.
    """
    out_dir = DATA_OUTPUT / "ts_ml" / f"{mol}_{label}"
    out_dir.mkdir(parents=True, exist_ok=True)

    reactant_xyz = out_dir / f"{mol}_reactant.xyz"
    product_xyz  = out_dir / f"{mol}_{label}_endpoint.xyz"
    reactant_xyz.write_text(mol_to_xyz(reactant, f"{mol} reactant"))
    product_xyz.write_text(mol_to_xyz(product, f"{mol} {label} endpoint"))

    run_dict = build_run_dict(reactant_xyz, product_xyz)

    original_cwd = Path.cwd()
    os.chdir(out_dir)
    try:
        return pysis_run.run_from_dict(run_dict, cwd=out_dir)
    finally:
        os.chdir(original_cwd)


def verify_ts_ml(hessian_h5_path: Path, atoms_of_interest: dict[str, int],
                  noise_thresh_cm1: float = 50.0) -> dict:
    """Human-review report for a completed AIMNet2-nse TS optimization: negative
    (imaginary) frequency count from the final Hessian, split into "significant"
    (|freq| >= noise_thresh_cm1) vs "possibly noise" (below it -- floppy/near-zero
    modes are common artifacts even in real DFT Hessians for flexible ring
    systems, per data/output/dft_opt/JOB_ISSUES.md's mol_020_E ring-pucker
    finding), and the dominant imaginary mode's per-atom displacement ranking.
    Does NOT auto-decide pass/fail, same precedent as beckmann.dft.ts.verify_ts.

    `atoms_of_interest` = {label: 0-based_atom_index} -- 0-based here, unlike
    verify_ts's 1-based Gaussian numbering, since this reads numpy arrays
    straight from PySisyphus's own atom ordering.
    """
    with h5py.File(hessian_h5_path, "r") as f:
        freqs = f["vibfreqs"][:]
        mw_displs = f["mw_cart_displs"][:]
        masses = f["masses"][:]

    negative_idx = np.where(freqs < 0)[0]
    significant = [(int(i), float(freqs[i])) for i in negative_idx if abs(freqs[i]) >= noise_thresh_cm1]
    possibly_noise = [(int(i), float(freqs[i])) for i in negative_idx if abs(freqs[i]) < noise_thresh_cm1]

    dominant_idx = int(negative_idx[np.argmin(freqs[negative_idx])]) if len(negative_idx) else None

    displacement_by_atom: dict[int, float] = {}
    if dominant_idx is not None:
        mode = mw_displs[:, dominant_idx].reshape(-1, 3)
        unweighted = mode / np.sqrt(masses)[:, None]  # mass-weighted -> plain Cartesian
        magnitudes = np.linalg.norm(unweighted, axis=1)
        displacement_by_atom = {i: float(m) for i, m in enumerate(magnitudes)}

    ranked = sorted(displacement_by_atom.items(), key=lambda t: -t[1])
    top_atoms = {idx for idx, _ in ranked[:6]}
    of_interest_ranks = {
        label: (idx, displacement_by_atom.get(idx))
        for label, idx in atoms_of_interest.items()
    }

    return {
        "n_imaginary_total": len(negative_idx),
        "n_imaginary_significant": len(significant),
        "significant_imaginary_freqs_cm1": significant,
        "possibly_noise_freqs_cm1": possibly_noise,
        "dominant_mode_freq_cm1": float(freqs[dominant_idx]) if dominant_idx is not None else None,
        "top_displaced_atoms": ranked[:6],
        "atoms_of_interest_displacement": of_interest_ranks,
        "atoms_of_interest_in_top6": {
            label: idx in top_atoms for label, (idx, _) in of_interest_ranks.items()
        },
    }


def print_verification_report_ml(mol: str, label: str, report: dict) -> None:
    print(f"-- {mol} {label} TS verification, AIMNet2-nse proxy (human review required) --")
    print(f"   Total imaginary (negative) frequencies: {report['n_imaginary_total']}")
    print(f"   Significant (|freq| >= 50 cm-1): {report['n_imaginary_significant']} "
          f"{report['significant_imaginary_freqs_cm1']}")
    if report["possibly_noise_freqs_cm1"]:
        print(f"   Possibly numerical noise (< 50 cm-1, common in floppy ring systems): "
              f"{report['possibly_noise_freqs_cm1']}")
    print(f"   Dominant imaginary mode: {report['dominant_mode_freq_cm1']} cm-1")
    print(f"   Top displaced atoms (0-based, |displacement|): {report['top_displaced_atoms']}")
    print(f"   Atoms of interest in top-6 displaced: {report['atoms_of_interest_in_top6']}")
    if report["n_imaginary_significant"] > 1:
        print("   -> WARNING: more than 1 significant imaginary frequency -- this is NOT a "
              "clean first-order TS on the AIMNet2-nse surface as converged. Do not treat "
              "as a verified TS without further investigation (mode-following, tighter "
              "optimization, or comparison against the DFT result).")


def main() -> None:
    """Pilot scope: mol_002_E rearrangement channel only."""
    from beckmann.dft.descriptors import _load_mols, get_substituent_map

    mol = "mol_002_E"
    label = "ts1_a1"
    reactant = _load_mols()[mol]
    product_sdf = DATA_OUTPUT / "aimnet_optimized" / f"{mol}_product_rearr.sdf"
    product = next(Chem.SDMolSupplier(str(product_sdf), removeHs=False))

    result = run_ts_search(mol, label, reactant, product)
    print(f"-- {mol} {label} (rearrangement, AIMNet2-nse proxy): {result}")

    atom_map = get_substituent_map(mol, DATA_OUTPUT / "dft_opt" / mol)  # 1-based
    atoms_of_interest_0based = {k: v - 1 for k, v in atom_map.items()}
    hessian_path = DATA_OUTPUT / "ts_ml" / f"{mol}_{label}" / "ts_final_hessian.h5"
    report = verify_ts_ml(hessian_path, atoms_of_interest_0based)
    print_verification_report_ml(mol, label, report)


if __name__ == "__main__":
    main()
