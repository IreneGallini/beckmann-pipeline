"""Optional: if 04_check_status.py flagged the Stage 3 scan as
OSCILLATING_DEGENERACY run CalcFC -> step=0.07 -> step=0.04 for this molecule. Safe
to re-run: each call just advances to the next untried rung, or does
nothing if already resolved or the ladder is exhausted. Upload+resubmit the
regenerated .gjf afterward

Edit MOL_NAME below, then:
    python 10_recover_if_crashed.py
"""
import sys

from _common import QUERY_PREFIX, sanitize_id, workdir_for

from beckmann_nbo.hpc import load_config, require_config
from beckmann_nbo.recovery import run_auto_recovery

MOL_NAME = "test1"
DRY_RUN = False


def main() -> None:
    mol_id = sanitize_id(MOL_NAME)
    dft_opt_dir = workdir_for(mol_id) / "dft_opt"

    matches = sorted(dft_opt_dir.glob(f"{QUERY_PREFIX}_{mol_id}_*"))
    if not matches:
        print(f"ERROR: no directory matching {QUERY_PREFIX}_{mol_id}_* under {dft_opt_dir}", file=sys.stderr)
        sys.exit(1)
    mol_name = matches[0].name

    config = load_config()
    require_config(config)

    run_auto_recovery([mol_name], dft_opt_dir, config, dry_run=DRY_RUN)
    print(
        "\nIf a new attempt was generated above, upload+resubmit it yourself "
        "(same cmd_upload/cmd_submit_scan calls as 07_upload_submit_stage3.py, "
        f"mol_id='{mol_id}'), then re-check with 04_check_status.py."
    )


if __name__ == "__main__":
    main()
