"""
Run AIMNet2/ASE optimization on the full benchmark's Auto3D conformers.
Moved here (from beckmann/optimize.py's main()) since it hardcodes the
benchmark's own paths -- select_and_optimize() itself is beckmann_core's
reusable function, called unmodified.
"""
from beckmann_core.optimize import select_and_optimize
from beckmann_nbo.config import DATA_OUTPUT


def main() -> None:
    conformers_dir = DATA_OUTPUT / "conformers"
    output_dir     = DATA_OUTPUT / "aimnet_optimized"

    sdf_files = sorted(conformers_dir.glob("molecules_*/molecules_out.sdf"))
    if not sdf_files:
        raise FileNotFoundError(
            "No Auto3D output SDF found. Run 01_smiles_to_conformers.py first."
        )
    select_and_optimize(sdf_files[-1], output_dir)


if __name__ == '__main__':
    main()
