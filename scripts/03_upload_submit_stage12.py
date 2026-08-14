"""Upload Stage 1+2 input files to your cluster (over SSH, per .env) and
submit Stage 1 (geometry optimization). Needs 02_prepare_stage12_inputs.py
to have already run, and .env filled in (see /README.md).

Edit MOL_NAME below, then:
    python 03_upload_submit_stage12.py
Set DRY_RUN = True first to preview the commands without running them.
"""
from _common import sanitize_id, workdir_for

from beckmann_nbo.hpc import cmd_submit_opt, cmd_upload, load_config, require_config

MOL_NAME = "test1"
DRY_RUN = False


def main() -> None:
    mol_id = sanitize_id(MOL_NAME)
    dft_opt_dir = workdir_for(mol_id) / "dft_opt"

    config = load_config()
    require_config(config)

    cmd_upload(config, DRY_RUN, mol_id, dft_opt_dir)
    cmd_submit_opt(config, DRY_RUN, mol_id, dft_opt_dir)

    print(
        f"\nSubmitted. Poll with 04_check_status.py, then pull results with "
        f"05_download_results.py once Stage 1 shows Normal termination."
    )


if __name__ == "__main__":
    main()
