"""
Step 1: SMILES → 3D conformers (k=5 per molecule, AIMNet2 for pre optimization)
Using Auto3D 2.3.1 - options/main API
"""

import os
import shutil
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")

from pathlib import Path
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"  # required on macOS: prevents libomp conflict

from Auto3D.auto3D import options, main

if __name__ == '__main__':
    PROJECT_ROOT = Path(__file__).parent.parent # Always finds the file relative to the script's location

    input_src  = PROJECT_ROOT / "data" / "input"  / "molecules.smi"
    output_dir = PROJECT_ROOT / "data" / "output" / "conformers"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Auto3D writes output next to wherever the input file lives
    # so we copy the input into the output folder first
    input_copy = output_dir / "molecules.smi"
    shutil.copy(input_src, input_copy)

    args = options(
        path=str(input_copy),
        k=5, #pick top 5 conformers
        optimizing_engine="AIMNET", # use AIMNet2 for ranking before optimization
        use_gpu=False,
        verbose=True,
    )

    out = main(args)
    print(f"\nDone. Output written to: {out}")