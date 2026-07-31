"""
Thin wrapper: Auto3D conformers -> AIMNet2/ASE optimization -> the single
winning (lowest-energy E/Z isomer) geometry. beckmann.optimize.
select_and_optimize() is vendored unchanged and already takes an explicit
SDF path + output dir with no benchmark-batch assumptions (only its own
main() hardcodes the benchmark's glob pattern -- not used here).
"""
from pathlib import Path

from rdkit import Chem

from beckmann_core.optimize import select_and_optimize


def run_optimize(conformers_sdf: Path, output_dir: Path) -> Chem.Mol:
    """Runs AIMNet2 optimization + E/Z substrate selection, returns the
    single winning RDKit Mol from best_per_substrate.sdf. Raises ValueError
    if optimization produced no usable structure."""
    _, best_per_substrate_sdf = select_and_optimize(conformers_sdf, output_dir)
    suppl = Chem.SDMolSupplier(str(best_per_substrate_sdf), removeHs=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        raise ValueError("AIMNet2 optimization produced no usable structure")
    return mols[0]
