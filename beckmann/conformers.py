"""
Step 1: SMILES → 3D conformers via Auto3D + AIMNet2 pre-ranking.
"""
import os
import shutil
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"  # required on macOS: prevents libomp conflict

from pathlib import Path

from Auto3D.auto3D import options, main as _auto3d_main

from beckmann.config import DATA_INPUT, DATA_OUTPUT


def generate_conformers(
    smiles_file: Path,
    output_dir: Path,
    k: int = 5,
) -> Path:
    """Run Auto3D on smiles_file, return path to the *_out.sdf produced.

    Auto3D writes its output next to its input file, so we copy smiles_file
    into output_dir first and run from there.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    input_copy = output_dir / smiles_file.name
    shutil.copy(smiles_file, input_copy)

    args = options(
        path=str(input_copy),
        k=k,
        optimizing_engine="AIMNET",
        use_gpu=False,
        verbose=True,
    )
    _auto3d_main(args)

    stem = smiles_file.stem
    sdf_files = sorted(output_dir.glob(f"{stem}_*/{stem}_out.sdf"))
    if not sdf_files:
        raise FileNotFoundError(
            f"Auto3D produced no output SDF in {output_dir}. "
            "Check that Auto3D ran without errors."
        )
    return sdf_files[-1]


def main() -> None:
    input_src  = DATA_INPUT  / "molecules.smi"
    output_dir = DATA_OUTPUT / "conformers"

    sdf_path = generate_conformers(input_src, output_dir)
    print(f"\nDone. Output written to: {sdf_path}")


if __name__ == '__main__':
    main()
