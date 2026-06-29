#!/usr/bin/env python3
"""
HPC synchronisation tool for Gaussian 16 DFT jobs on Citadel.

Citadel (citadel.chem.cmu.edu) is a shared Ubuntu compute server — no SLURM.
Gaussian 16 is at /opt/g16/g16. Jobs run via nohup and survive disconnection.

Two job types are supported:
  Two-stage (dft_opt/):  submit-opt → submit-nbo
    Used for DFT geometry optimisation followed by NBO single-point.
    Input files: {name}_opt.gjf and {name}_nbo.gjf
  Single-point (dft_sp/): submit-sp
    Used for NBO single-point directly on AIMNet2 geometry.
    Input files: {name}.gjf

Typical two-stage workflow:
  python scripts/dft/hpc_sync.py --mol 002 upload
  python scripts/dft/hpc_sync.py --mol 002 submit-opt
  python scripts/dft/hpc_sync.py status
  python scripts/dft/hpc_sync.py --mol 002 submit-nbo   # after Stage 1 finishes
  python scripts/dft/hpc_sync.py --mol 002 download

For dft_sp/ (single-point jobs):
  python scripts/dft/hpc_sync.py --dir data/output/dft_sp upload
  python scripts/dft/hpc_sync.py --dir data/output/dft_sp submit-sp
  python scripts/dft/hpc_sync.py --dir data/output/dft_sp download

Flags (must come BEFORE the subcommand):
  --mol 002     target mol_002_E and mol_002_Z only
  --dir PATH    local directory to sync (default: data/output/dft_opt)
  --dry-run     print commands without executing

Configuration — copy .env.example to .env (gitignored) and fill in:
  HPC_HOST       igallini@citadel.chem.cmu.edu
  HPC_REMOTE_DIR ~/beckmann/dft_opt
  G16_PATH       /opt/g16/g16
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT          = Path(__file__).parent.parent.parent
DEFAULT_LOCAL_DFT_DIR = PROJECT_ROOT / "data" / "output" / "dft_opt"
ENV_FILE              = PROJECT_ROOT / ".env"


def load_config() -> dict:
    """Load HPC settings from .env, then override with os.environ."""
    config: dict[str, str] = {}
    if ENV_FILE.exists():
        for raw in ENV_FILE.read_text().splitlines():
            line = raw.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, value = line.partition("=")
            config[key.strip()] = value.strip()
    for key in ("HPC_HOST", "HPC_REMOTE_DIR", "G16_PATH"):
        if key in os.environ:
            config[key] = os.environ[key]
    return config


def require_config(config: dict) -> None:
    missing = [k for k in ("HPC_HOST", "HPC_REMOTE_DIR", "G16_PATH") if not config.get(k)]
    if missing:
        print(f"ERROR: missing config: {', '.join(missing)}", file=sys.stderr)
        print("  1. Copy .env.example to .env and fill in all values", file=sys.stderr)
        print("  2. HPC_HOST=username@hostname", file=sys.stderr)
        print("  3. HPC_REMOTE_DIR=~/beckmann/dft_opt", file=sys.stderr)
        print("  4. G16_PATH=/opt/g16/g16  (full path to g16 on the server)", file=sys.stderr)
        sys.exit(1)


def mol_dirs(local_dir: Path, mol: str | None) -> list[Path]:
    """Return local molecule directories to operate on."""
    if mol:
        dirs = sorted(local_dir.glob(f"mol_{mol.zfill(3)}_*/"))
        if not dirs:
            print(f"ERROR: no directories found for mol {mol} in {local_dir}", file=sys.stderr)
            sys.exit(1)
        return dirs
    return sorted(local_dir.glob("mol_*/"))


def run(cmd: list[str], dry_run: bool) -> None:
    """Print a command, then execute it unless dry_run is True."""
    print(f"$ {' '.join(str(c) for c in cmd)}")
    if dry_run:
        print("  [dry-run: not executed]")
        return
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"ERROR: command exited {result.returncode}", file=sys.stderr)
        sys.exit(result.returncode)


def cmd_upload(config: dict, dry_run: bool, mol: str | None, local_dir: Path) -> None:
    """Upload molecule directories to the cluster."""
    if not local_dir.exists():
        print(f"ERROR: {local_dir} does not exist.", file=sys.stderr)
        sys.exit(1)

    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    dirs       = mol_dirs(local_dir, mol)

    print(f"\n-- Upload: {', '.join(d.name for d in dirs)}  →  {host}:{remote_dir}/")
    run(["ssh", host, f"mkdir -p {remote_dir}"], dry_run)

    if mol:
        for d in dirs:
            run(["scp", "-r", str(d), f"{host}:{remote_dir}/"], dry_run)
    else:
        remote_parent = str(Path(remote_dir).parent)
        run(["scp", "-r", str(local_dir), f"{host}:{remote_parent}/"], dry_run)

    print("\nUpload complete.")


def cmd_submit_opt(config: dict, dry_run: bool, mol: str | None, local_dir: Path) -> None:
    """SSH into the server and run Stage 1 (opt) jobs via nohup g16."""
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    pattern    = f"mol_{mol.zfill(3)}_*" if mol else "*"
    g16          = config["G16_PATH"]
    gauss_exedir = str(Path(g16).parent)
    g16root      = str(Path(g16).parent.parent)
    # Export Gaussian env vars explicitly — non-interactive SSH shells do not
    # source ~/.bashrc so GAUSS_EXEDIR and g16root would otherwise be unset,
    # causing g16 to segfault on startup.
    submit_cmd = (
        f'export GAUSS_EXEDIR={gauss_exedir} && '
        f'export g16root={g16root} && '
        f'cd {remote_dir} && '
        f'for dir in {pattern}/; do '
        '  name="${dir%/}"; '
        f'  (cd "$dir" && nohup {g16} < "${{name}}_opt.gjf" > "${{name}}_opt.log" 2>&1 &); '
        'done'
    )
    print(f"\n-- Launching Stage 1 (opt) jobs on {host}:{remote_dir}")
    run(["ssh", host, submit_cmd], dry_run)
    print("\nStage 1 launched. Monitor with:\n  python scripts/dft/hpc_sync.py status")


def cmd_submit_nbo(config: dict, dry_run: bool, mol: str | None, local_dir: Path) -> None:
    """SSH into the server and run Stage 2 (NBO) jobs via nohup g16."""
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    print(
        "\nWARNING: NBO jobs read the .chk written by Stage 1.\n"
        "         Only proceed once ALL opt jobs have COMPLETED.\n"
        "         Check first: python scripts/dft/hpc_sync.py status"
    )
    pattern      = f"mol_{mol.zfill(3)}_*" if mol else "*"
    g16          = config["G16_PATH"]
    gauss_exedir = str(Path(g16).parent)
    g16root      = str(Path(g16).parent.parent)
    submit_cmd = (
        f'export GAUSS_EXEDIR={gauss_exedir} && '
        f'export g16root={g16root} && '
        f'cd {remote_dir} && '
        f'for dir in {pattern}/; do '
        '  name="${dir%/}"; '
        f'  (cd "$dir" && nohup {g16} < "${{name}}_nbo.gjf" > "${{name}}_nbo.log" 2>&1 &); '
        'done'
    )
    print(f"\n-- Launching Stage 2 (NBO) jobs on {host}:{remote_dir}")
    run(["ssh", host, submit_cmd], dry_run)
    print("\nStage 2 launched. Monitor with:\n  python scripts/dft/hpc_sync.py status")


def cmd_submit_sp(config: dict, dry_run: bool, mol: str | None, local_dir: Path) -> None:
    """SSH into the server and run single-point NBO jobs (dft_sp/ style)."""
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    pattern    = f"mol_{mol.zfill(3)}_*" if mol else "*"
    g16          = config["G16_PATH"]
    gauss_exedir = str(Path(g16).parent)
    g16root      = str(Path(g16).parent.parent)
    submit_cmd = (
        f'export GAUSS_EXEDIR={gauss_exedir} && '
        f'export g16root={g16root} && '
        f'cd {remote_dir} && '
        f'for dir in {pattern}/; do '
        '  name="${dir%/}"; '
        f'  (cd "$dir" && nohup {g16} < "${{name}}.gjf" > "${{name}}.log" 2>&1 &); '
        'done'
    )
    print(f"\n-- Launching single-point jobs on {host}:{remote_dir}")
    run(["ssh", host, submit_cmd], dry_run)
    print("\nJobs launched. Monitor with:\n  python scripts/dft/hpc_sync.py status")


def cmd_download(config: dict, dry_run: bool, mol: str | None, local_dir: Path) -> None:
    """Download *.log files from the cluster into the local directory."""
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    local_dir.mkdir(parents=True, exist_ok=True)

    pattern = f"mol_{mol.zfill(3)}_*" if mol else "*"
    print(f"\n-- Downloading *.log ({pattern}) from {host}:{remote_dir}/")
    run([
        "rsync", "-avz",
        f"--filter=+ {pattern}/",
        "--include=*/",
        "--include=*.log",
        "--exclude=*",
        f"{host}:{remote_dir}/",
        str(local_dir) + "/",
    ], dry_run)
    print(
        f"\nLogs written to {local_dir}/\n"
        "Note: *.log is gitignored — these files will not be committed."
    )


def cmd_status(config: dict, dry_run: bool, mol: str | None, local_dir: Path) -> None:
    """Show running g16 processes on the server."""
    host = config["HPC_HOST"]
    print(f"\n-- Running g16 processes on {host}")
    run(["ssh", host, "ps aux | grep '[g]16' | awk '{print $1, $2, $11, $12}'"], dry_run)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="HPC sync tool for Gaussian DFT jobs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples (two-stage opt→nbo, dft_opt/):\n"
            "  python scripts/dft/hpc_sync.py --mol 002 --dry-run upload\n"
            "  python scripts/dft/hpc_sync.py --mol 002 upload\n"
            "  python scripts/dft/hpc_sync.py --mol 002 submit-opt\n"
            "  python scripts/dft/hpc_sync.py status\n"
            "  python scripts/dft/hpc_sync.py --mol 002 submit-nbo\n"
            "  python scripts/dft/hpc_sync.py --mol 002 download\n"
            "\nexamples (single-point, dft_sp/):\n"
            "  python scripts/dft/hpc_sync.py --dir data/output/dft_sp upload\n"
            "  python scripts/dft/hpc_sync.py --dir data/output/dft_sp submit-sp\n"
            "  python scripts/dft/hpc_sync.py --dir data/output/dft_sp download\n"
            "\nnote: all flags (--mol, --dir, --dry-run) must come BEFORE the subcommand\n"
        ),
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print every command that would run, but do not execute it.",
    )
    parser.add_argument(
        "--mol", metavar="ID",
        help="Target a single molecule by ID (e.g. --mol 002). "
             "Operates on both E and Z isomers. Omit to process all molecules.",
    )
    parser.add_argument(
        "--dir", metavar="PATH", default=None,
        help="Local job directory to sync (default: data/output/dft_opt). "
             "Use data/output/dft_sp for single-point jobs.",
    )

    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True
    sub.add_parser("upload",      help="Upload molecule directories to cluster")
    sub.add_parser("submit-opt",  help="Submit Stage 1 geometry-opt jobs ({name}_opt.gjf)")
    sub.add_parser("submit-nbo",  help="Submit Stage 2 NBO jobs ({name}_nbo.gjf) — AFTER Stage 1 finishes")
    sub.add_parser("submit-sp",   help="Submit single-point NBO jobs ({name}.gjf, for dft_sp/)")
    sub.add_parser("download",    help="Download *.log files from cluster")
    sub.add_parser("status",      help="Show running g16 processes on server")

    args = parser.parse_args()

    local_dir = Path(args.dir) if args.dir else DEFAULT_LOCAL_DFT_DIR
    if not local_dir.is_absolute():
        local_dir = PROJECT_ROOT / local_dir

    config = load_config()
    require_config(config)

    dispatch = {
        "upload":     cmd_upload,
        "submit-opt": cmd_submit_opt,
        "submit-nbo": cmd_submit_nbo,
        "submit-sp":  cmd_submit_sp,
        "download":   cmd_download,
        "status":     cmd_status,
    }
    dispatch[args.command](config, args.dry_run, args.mol, local_dir)


if __name__ == "__main__":
    main()
