#!/usr/bin/env python3
"""
HPC synchronisation tool for Gaussian 16 DFT jobs.

Workflow (two-stage):
  1. python scripts/hpc_sync.py upload         # push input files to cluster
  2. python scripts/hpc_sync.py submit-opt     # Stage 1: geometry optimisation
     ... wait for Stage 1 jobs to finish ...
     python scripts/hpc_sync.py status         # check squeue --me
  3. python scripts/hpc_sync.py submit-nbo     # Stage 2: NBO single-point
     ... wait for Stage 2 jobs to finish ...
  4. python scripts/hpc_sync.py download       # pull *.log files back

Configuration:
  Copy .env.example to .env and fill in:
    HPC_HOST       e.g. igallini@submit.cluster.edu
                   or an alias from ~/.ssh/config
    HPC_REMOTE_DIR e.g. ~/beckmann/dft_opt_test

Credentials are never stored in committed code. All commands are printed
before execution. Use --dry-run to preview without executing.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT  = Path(__file__).parent.parent
LOCAL_DFT_DIR = PROJECT_ROOT / "data" / "output" / "dft_opt_test"
ENV_FILE      = PROJECT_ROOT / ".env"


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
    for key in ("HPC_HOST", "HPC_REMOTE_DIR"):
        if key in os.environ:
            config[key] = os.environ[key]
    return config


def require_config(config: dict) -> None:
    missing = [k for k in ("HPC_HOST", "HPC_REMOTE_DIR") if not config.get(k)]
    if missing:
        print(f"ERROR: missing config: {', '.join(missing)}", file=sys.stderr)
        print("  1. Copy .env.example to .env", file=sys.stderr)
        print("  2. Set HPC_HOST=username@cluster.hostname", file=sys.stderr)
        print("  3. Set HPC_REMOTE_DIR=~/beckmann/dft_opt_test", file=sys.stderr)
        print("  Or: export HPC_HOST=... before running this script.", file=sys.stderr)
        sys.exit(1)


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


def cmd_upload(config: dict, dry_run: bool) -> None:
    """Upload data/output/dft_opt_test/ to the cluster."""
    if not LOCAL_DFT_DIR.exists():
        print(
            f"ERROR: {LOCAL_DFT_DIR} does not exist.\n"
            "Run: python scripts/05_prepare_test_opt.py first.",
            file=sys.stderr,
        )
        sys.exit(1)

    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    # e.g. ~/beckmann/dft_opt_test → ~/beckmann
    remote_parent = str(Path(remote_dir).parent)

    print(f"\n-- Upload: {LOCAL_DFT_DIR}  →  {host}:{remote_dir}")
    run(["ssh", host, f"mkdir -p {remote_dir}"], dry_run)
    # scp -r copies the directory by name into remote_parent/,
    # creating remote_parent/dft_opt_test/ which matches HPC_REMOTE_DIR.
    run(["scp", "-r", str(LOCAL_DFT_DIR), f"{host}:{remote_parent}/"], dry_run)
    print("\nUpload complete.")


def cmd_submit_opt(config: dict, dry_run: bool) -> None:
    """SSH into the cluster and sbatch all Stage 1 (opt) jobs."""
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    # Parenthesized subshell: each molecule's cd is isolated, so a sbatch
    # failure does not corrupt the loop's working directory.
    submit_cmd = (
        f"cd {remote_dir} && "
        "for dir in */; do "
        '  name="${dir%/}"; '
        '  (cd "$dir" && sbatch "${name}_opt_submit.sh"); '
        "done"
    )
    print(f"\n-- Submitting Stage 1 (opt) jobs on {host}:{remote_dir}")
    run(["ssh", host, submit_cmd], dry_run)
    print("\nStage 1 submitted. Monitor with:\n  python scripts/hpc_sync.py status")


def cmd_submit_nbo(config: dict, dry_run: bool) -> None:
    """SSH into the cluster and sbatch all Stage 2 (NBO) jobs."""
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    print(
        "\nWARNING: NBO jobs read the .chk written by Stage 1.\n"
        "         Only proceed once ALL opt jobs have COMPLETED.\n"
        "         Check first: python scripts/hpc_sync.py status"
    )
    submit_cmd = (
        f"cd {remote_dir} && "
        "for dir in */; do "
        '  name="${dir%/}"; '
        '  (cd "$dir" && sbatch "${name}_nbo_submit.sh"); '
        "done"
    )
    print(f"\n-- Submitting Stage 2 (NBO) jobs on {host}:{remote_dir}")
    run(["ssh", host, submit_cmd], dry_run)
    print("\nStage 2 submitted. Monitor with:\n  python scripts/hpc_sync.py status")


def cmd_download(config: dict, dry_run: bool) -> None:
    """Download *.log files from the cluster into data/output/dft_opt_test/."""
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    LOCAL_DFT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"\n-- Downloading *.log from {host}:{remote_dir}/")
    # rsync filter order matters: allow dirs, allow *.log, exclude everything else.
    run([
        "rsync", "-avz",
        "--include=*/",
        "--include=*.log",
        "--exclude=*",
        f"{host}:{remote_dir}/",
        str(LOCAL_DFT_DIR) + "/",
    ], dry_run)
    print(
        f"\nLogs written to {LOCAL_DFT_DIR}/\n"
        "Note: *.log is gitignored — these files will not be committed."
    )


def cmd_status(config: dict, dry_run: bool) -> None:
    """Show the current SLURM queue for this user."""
    host = config["HPC_HOST"]
    print(f"\n-- SLURM queue on {host}")
    run(["ssh", host, "squeue --me"], dry_run)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="HPC sync tool for Gaussian DFT jobs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  python scripts/hpc_sync.py --dry-run upload\n"
            "  python scripts/hpc_sync.py upload\n"
            "  python scripts/hpc_sync.py submit-opt\n"
            "  python scripts/hpc_sync.py status\n"
            "  python scripts/hpc_sync.py submit-nbo\n"
            "  python scripts/hpc_sync.py download\n"
            "\nnote: --dry-run must come BEFORE the subcommand\n"
        ),
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print every command that would run, but do not execute it.",
    )

    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True
    sub.add_parser("upload",      help="Upload dft_opt_test/ to cluster")
    sub.add_parser("submit-opt",  help="Submit Stage 1 geometry-opt jobs")
    sub.add_parser("submit-nbo",  help="Submit Stage 2 NBO jobs (AFTER Stage 1 finishes)")
    sub.add_parser("download",    help="Download *.log files from cluster")
    sub.add_parser("status",      help="Show SLURM queue (squeue --me)")

    args = parser.parse_args()

    config = load_config()
    require_config(config)

    dispatch = {
        "upload":     cmd_upload,
        "submit-opt": cmd_submit_opt,
        "submit-nbo": cmd_submit_nbo,
        "download":   cmd_download,
        "status":     cmd_status,
    }
    dispatch[args.command](config, args.dry_run)


if __name__ == "__main__":
    main()
