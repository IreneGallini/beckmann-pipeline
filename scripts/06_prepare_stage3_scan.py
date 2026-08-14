"""Stage 3 input generation: the rigid N-O bond scan (rigid displacement ->
constrained re-optimization -> NBO7 single point, at each of 6 points by
default). Can only run after Stage 1's *converged* geometry is downloaded
(05_download_results.py) this script checks for "Normal termination" in
the Stage 1 log before generating anything.

Edit MOL_NAME below, then:
    python 06_prepare_stage3_scan.py
"""
import sys

from _common import QUERY_PREFIX, sanitize_id, workdir_for

from beckmann_nbo.inputs import prepare_scan_rigid

MOL_NAME = "test1"


def main() -> None:
    mol_id = sanitize_id(MOL_NAME)
    dft_opt_dir = workdir_for(mol_id) / "dft_opt"

    matches = sorted(dft_opt_dir.glob(f"{QUERY_PREFIX}_{mol_id}_*"))
    if not matches:
        print(f"ERROR: no directory matching {QUERY_PREFIX}_{mol_id}_* under {dft_opt_dir}", file=sys.stderr)
        sys.exit(1)
    mol_dir = matches[0]
    mol_name = mol_dir.name

    opt_log = mol_dir / f"{mol_name}_opt.log"
    if not opt_log.exists() or "Normal termination" not in opt_log.read_text()[-2000:]:
        print(
            f"ERROR: {opt_log} not found or Stage 1 hasn't reached Normal "
            f"termination yet. Check with 04_check_status.py.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Generating Stage 3 (N-O scan) for {mol_name}...")
    prepare_scan_rigid(mol_dir, mol_name)
    print(f"\nNext: edit MOL_NAME in 07_upload_submit_stage3.py to '{MOL_NAME}' and run it.")


if __name__ == "__main__":
    main()
