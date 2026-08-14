"""Stage 1+2 input generation: writes Gaussian {name}_opt.gjf (DFT geometry
optimization) and {name}_nbo.gjf (NBO7 single point) for this molecule.
Needs 01_optimize_aimnet2.py to have already run.

Edit MOL_NAME below to match what you used in 00/01, then:
    python 02_prepare_stage12_inputs.py
"""
from _common import sanitize_id, workdir_for

from beckmann_nbo.inputs import prepare_opt

MOL_NAME = "test1"


def main() -> None:
    mol_id = sanitize_id(MOL_NAME)
    workdir = workdir_for(mol_id)
    sub_sdf = workdir / "aimnet_optimized" / "best_per_substrate.sdf"
    dft_opt_dir = workdir / "dft_opt"

    prepare_opt(sub_sdf, dft_opt_dir, test_ids={mol_id})
    print(f"\nNext: edit MOL_NAME in 03_upload_submit_stage12.py to '{MOL_NAME}' and run it.")


if __name__ == "__main__":
    main()
