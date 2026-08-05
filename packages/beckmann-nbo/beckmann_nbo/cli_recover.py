"""`beckmann-nbo recover` -- exact wrapper around scripts/auto_recover.py's
body, exposed as a first-class subcommand instead of a script the user has
to go find. No new logic: same four calls auto_recover.py itself makes.
"""
from pathlib import Path

from beckmann_nbo.hpc import DEFAULT_LOCAL_DFT_DIR, load_config, mol_dirs, require_config
from beckmann_nbo.recovery import run_auto_recovery


def cmd_recover(args) -> None:
    config = load_config()
    require_config(config)

    dft_opt_dir = Path(args.dir) if args.dir else DEFAULT_LOCAL_DFT_DIR
    mols = [d.name for d in mol_dirs(dft_opt_dir, args.mol)]
    run_auto_recovery(mols, dft_opt_dir, config, dry_run=args.dry_run)
