"""Check job status on your cluster for this molecule (per-stage: opt/nbo/
scan) and locally classify any downloaded logs (catches the ring-pucker
oscillation crash that hits roughly 1 in 5 substrates see
recovery.py/../README.md). Safe to re-run anytime; reuse this same script
after Stage 3 is submitted (07_upload_submit_stage3.py) to poll it too.

Edit MOL_NAME below, then:
    python 04_check_status.py
"""
from _common import sanitize_id, workdir_for

from beckmann_nbo.hpc import cmd_status, load_config, mol_dirs, require_config
from beckmann_nbo.log_diagnostics import FailureCategory, classify_scan

MOL_NAME = "test1"
DRY_RUN = False


def main() -> None:
    mol_id = sanitize_id(MOL_NAME)
    dft_opt_dir = workdir_for(mol_id) / "dft_opt"

    config = load_config()
    require_config(config)

    cmd_status(config, DRY_RUN, mol_id, dft_opt_dir)

    for d in mol_dirs(dft_opt_dir, mol_id):
        print(f"\n{d.name}")
        diagnoses = classify_scan(d, d.name)
        if not diagnoses:
            print("  no stage logs downloaded yet -- run 05_download_results.py")
            continue
        by_stage = {diag.stage: diag for diag in diagnoses}
        for stage in ("opt", "nbo", "scan"):
            if stage not in by_stage:
                continue
            diag = by_stage[stage]
            flag = "" if diag.category == FailureCategory.NORMAL else "  <-- needs attention"
            print(f"  Stage {stage:<5} {diag.category.value:<24}{flag}")


if __name__ == "__main__":
    main()
