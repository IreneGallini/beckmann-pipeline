"""Upload the Stage 3 scan input and submit it. Needs
06_prepare_stage3_scan.py to have already run.

Edit MOL_NAME below, then:
    python 07_upload_submit_stage3.py
Then reuse 04_check_status.py / 05_download_results.py to poll/pull it.
"""
from _common import sanitize_id, workdir_for

from beckmann_nbo.hpc import cmd_submit_scan, cmd_upload, load_config, require_config

MOL_NAME = "test1"
DRY_RUN = False


def main() -> None:
    mol_id = sanitize_id(MOL_NAME)
    dft_opt_dir = workdir_for(mol_id) / "dft_opt"

    config = load_config()
    require_config(config)

    cmd_upload(config, DRY_RUN, mol_id, dft_opt_dir)
    cmd_submit_scan(config, DRY_RUN, mol_id, dft_opt_dir)

    print("\nSubmitted. Poll with 04_check_status.py.")


if __name__ == "__main__":
    main()
