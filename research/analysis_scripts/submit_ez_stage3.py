"""
Step 3c helper: given a priority-molecule "other isomer" name whose Stage 1
(_opt.log) has already been downloaded to data/output/dft_opt_ez_other/,
compute R0 from the converged geometry, derive n_points to reach the
requested ~1.85 A ceiling at step=0.05 A, and generate the Stage 3
(_scan.gjf) rigid-scan input via beckmann_nbo.inputs.prepare_scan_rigid().

Usage: python research/analysis_scripts/submit_ez_stage3.py mol_003_Z
"""
import math
import sys
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from beckmann_core.geometry import no_distance
from beckmann_nbo.config import DATA_OUTPUT
from beckmann_nbo.geometry import parse_standard_orientations
from beckmann_nbo.inputs import prepare_scan_rigid
from beckmann_nbo.scan import oxime_atom_map_from_gjf

DFT_OPT_DIR = DATA_OUTPUT / "dft_opt_ez_other"
CEILING = 1.85
STEP = 0.05


def generate_stage3(name: str):
    mol_dir = DFT_OPT_DIR / name
    ci, ni, oi, label = oxime_atom_map_from_gjf(mol_dir / f"{name}_opt.gjf")
    lines = (mol_dir / f"{name}_opt.log").read_text().splitlines()
    base_atoms = parse_standard_orientations(lines)[-1][1]
    r0 = no_distance(base_atoms, ni, oi)
    n_points = max(1, math.floor((CEILING - r0) / STEP))
    out = prepare_scan_rigid(mol_dir, name, step=STEP, n_points=n_points)
    print(f"{name}: R0={r0:.4f} A, n_points={n_points} (reaches R0+{n_points*STEP:.2f}={r0+n_points*STEP:.4f} A) -> {out}")
    return out


if __name__ == "__main__":
    for name in sys.argv[1:]:
        generate_stage3(name)
