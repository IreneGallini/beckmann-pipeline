"""
Run one pass of the automated scan-failure recovery ladder (see
beckmann/dft/recovery.py): for each in-scope molecule, classify its
canonical Stage 3 scan and, if it's an oscillating-degeneracy crash,
automatically generate/upload/submit the next escalation rung
(CalcFC -> step=0.07 -> step=0.04).

Meant to be re-run periodically (by hand or via cron) -- a single
invocation does one "check everything, act on what needs it" round, not a
persistent background service, since Citadel has no job scheduler and this
is a research pipeline, not a production system.

Usage:
  python scripts/dft/auto_recover.py                # all molecules
  python scripts/dft/auto_recover.py --mol 020       # just mol_020_E
  python scripts/dft/auto_recover.py --dry-run       # preview only
"""
import argparse

from beckmann.config import DATA_OUTPUT
from beckmann.dft.hpc import load_config, mol_dirs, require_config
from beckmann.dft.recovery import run_auto_recovery


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mol", metavar="ID", help="Target a single molecule by ID (e.g. --mol 020).")
    parser.add_argument("--dry-run", action="store_true", help="Preview without uploading/submitting.")
    args = parser.parse_args()

    config = load_config()
    require_config(config)

    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    mols = [d.name for d in mol_dirs(dft_opt_dir, args.mol)]
    run_auto_recovery(mols, dft_opt_dir, config, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
