"""Download whatever .log files exist on the cluster for this molecule.
Safe to re-run anytime reuse this same script after Stage 3 is submitted
(07_upload_submit_stage3.py) to pull that too.

Edit MOL_NAME below, then:
    python 05_download_results.py
"""
from _common import sanitize_id, workdir_for

from beckmann_nbo.hpc import cmd_download, load_config, require_config

MOL_NAME = "test1"
DRY_RUN = False


def main() -> None:
    mol_id = sanitize_id(MOL_NAME)
    dft_opt_dir = workdir_for(mol_id) / "dft_opt"

    config = load_config()
    require_config(config)

    cmd_download(config, DRY_RUN, mol_id, dft_opt_dir)

    print(
        "\nCheck each downloaded .log's last line for 'Normal termination of "
        "Gaussian 16' before moving on (run 04_check_status.py for an "
        "automated read of this)."
    )


if __name__ == "__main__":
    main()
