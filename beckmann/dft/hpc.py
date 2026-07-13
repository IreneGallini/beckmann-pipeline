"""
HPC synchronisation tool for Gaussian 16 DFT jobs on Citadel.

Citadel (citadel.chem.cmu.edu) is a shared Ubuntu compute server — no SLURM.
Gaussian 16 is at /opt/g16/g16. Jobs run via nohup and survive disconnection.

Typical two-stage workflow:
  python scripts/dft/hpc_sync.py --mol 002 upload
  python scripts/dft/hpc_sync.py --mol 002 submit-opt
  python scripts/dft/hpc_sync.py status
  python scripts/dft/hpc_sync.py --mol 002 submit-scan  # after Stage 1 finishes
  python scripts/dft/hpc_sync.py --mol 002 submit-nbo   # optional equilibrium NBO
  python scripts/dft/hpc_sync.py --mol 002 download

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

from beckmann.config import ROOT

PROJECT_ROOT          = ROOT
DEFAULT_LOCAL_DFT_DIR = ROOT / "data" / "output" / "dft_opt"
ENV_FILE              = ROOT / ".env"


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
    for key in ("HPC_HOST", "HPC_REMOTE_DIR", "G16_PATH", "NBOEXE", "NBO_WRAPPER_DIR"):
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
    if not config.get("NBOEXE"):
        print("WARNING: NBOEXE not set — Gaussian will use bundled NBO 3.1 (no CMO support).", file=sys.stderr)
        print("         Set NBOEXE=/opt/nbo7/bin/nbo7.i8.exe in .env to use NBO7.", file=sys.stderr)
    if config.get("NBOEXE") and not config.get("NBO_WRAPPER_DIR"):
        print("WARNING: NBO_WRAPPER_DIR not set — pop=nbo7read needs the gaunbo7 script",
              file=sys.stderr)
        print("         on PATH with execute permission (it ships read-only under NBOEXE's",
              file=sys.stderr)
        print("         directory). Copy it somewhere you own and chmod +x, then set",
              file=sys.stderr)
        print("         NBO_WRAPPER_DIR=~/beckmann/nbo7_bin in .env.", file=sys.stderr)


def _gauss_exports(config: dict) -> str:
    """Build the export prefix for non-interactive SSH shells."""
    g16          = config["G16_PATH"]
    gauss_exedir = str(Path(g16).parent)
    g16root      = str(Path(g16).parent.parent)
    exports = f'export GAUSS_EXEDIR={gauss_exedir} && export g16root={g16root}'
    if config.get("NBOEXE"):
        nbo_bin = str(Path(config["NBOEXE"]).parent)
        exports += f' && export GAUSS_EXEDIR={nbo_bin}:{gauss_exedir}'
        exports += f' && export PATH={nbo_bin}:$PATH'
    if config.get("NBO_WRAPPER_DIR"):
        # gaunbo7/gaunbo6 under NBOEXE's directory are read-only (root-owned);
        # pop=nbo7read finds the executable copies here via a plain PATH search
        # (Gaussian's Link 612 external-program interface, not GAUSS_EXEDIR).
        exports += f' && export PATH={config["NBO_WRAPPER_DIR"]}:$PATH'
    return exports


def mol_dirs(local_dir: Path, mol: str | None) -> list[Path]:
    if mol:
        dirs = sorted(local_dir.glob(f"mol_{mol.zfill(3)}_*/"))
        if not dirs:
            print(f"ERROR: no directories found for mol {mol} in {local_dir}", file=sys.stderr)
            sys.exit(1)
        return dirs
    return sorted(local_dir.glob("mol_*/"))


def run(cmd: list[str], dry_run: bool) -> None:
    print(f"$ {' '.join(str(c) for c in cmd)}")
    if dry_run:
        print("  [dry-run: not executed]")
        return
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"ERROR: command exited {result.returncode}", file=sys.stderr)
        sys.exit(result.returncode)


def cmd_upload(config: dict, dry_run: bool, mol: str | None, local_dir: Path) -> None:
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
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    pattern    = f"mol_{mol.zfill(3)}_*" if mol else "*"
    g16        = config["G16_PATH"]
    submit_cmd = (
        f'{_gauss_exports(config)} && '
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
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    print(
        "\nWARNING: NBO jobs read the .chk written by Stage 1.\n"
        "         Only proceed once ALL opt jobs have COMPLETED.\n"
        "         Check first: python scripts/dft/hpc_sync.py status"
    )
    pattern = f"mol_{mol.zfill(3)}_*" if mol else "*"
    g16     = config["G16_PATH"]
    submit_cmd = (
        f'{_gauss_exports(config)} && '
        f'cd {remote_dir} && '
        f'for dir in {pattern}/; do '
        '  name="${dir%/}"; '
        f'  (cd "$dir" && nohup {g16} < "${{name}}_nbo.gjf" > "${{name}}_nbo.log" 2>&1 &); '
        'done'
    )
    print(f"\n-- Launching Stage 2 (NBO) jobs on {host}:{remote_dir}")
    run(["ssh", host, submit_cmd], dry_run)
    print("\nStage 2 launched. Monitor with:\n  python scripts/dft/hpc_sync.py status")


def cmd_submit_scan(config: dict, dry_run: bool, mol: str | None, local_dir: Path) -> None:
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    print(
        "\nWARNING: Scan jobs read _opt.chk from Stage 1.\n"
        "         Only proceed once ALL opt jobs show Normal termination.\n"
        "         Check first: python scripts/dft/hpc_sync.py status"
    )
    pattern = f"mol_{mol.zfill(3)}_*" if mol else "*"
    g16     = config["G16_PATH"]
    submit_cmd = (
        f'{_gauss_exports(config)} && '
        f'cd {remote_dir} && '
        f'for dir in {pattern}/; do '
        '  name="${dir%/}"; '
        f'  (cd "$dir" && nohup {g16} < "${{name}}_scan.gjf" > "${{name}}_scan.log" 2>&1 &); '
        'done'
    )
    print(f"\n-- Launching Stage 3 (scan) jobs on {host}:{remote_dir}")
    run(["ssh", host, submit_cmd], dry_run)
    print("\nStage 3 launched. Monitor with:\n  python scripts/dft/hpc_sync.py status")


def cmd_submit_scan_sp(config: dict, dry_run: bool, mol: str | None, local_dir: Path) -> None:
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    pattern    = f"mol_{mol.zfill(3)}_*" if mol else "*"
    g16        = config["G16_PATH"]
    submit_cmd = (
        f'{_gauss_exports(config)} && '
        f'cd {remote_dir} && '
        f'for dir in {pattern}/; do '
        '  name="${dir%/}"; '
        f'  for sp in 2 3 4 5; do '
        f'    gjf="${{name}}_sp${{sp}}.gjf"; '
        f'    [ -f "$dir/$gjf" ] && (cd "$dir" && nohup {g16} < "$gjf" > "${{gjf%.gjf}}.log" 2>&1 &); '
        '  done; '
        'done'
    )
    print(f"\n-- Launching intermediate scan SP jobs on {host}:{remote_dir}")
    run(["ssh", host, submit_cmd], dry_run)
    print("\nJobs launched. Monitor with:\n  python scripts/dft/hpc_sync.py status")


def cmd_submit_sp(config: dict, dry_run: bool, mol: str | None, local_dir: Path) -> None:
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    pattern    = f"mol_{mol.zfill(3)}_*" if mol else "*"
    g16        = config["G16_PATH"]
    submit_cmd = (
        f'{_gauss_exports(config)} && '
        f'cd {remote_dir} && '
        f'for dir in {pattern}/; do '
        '  name="${dir%/}"; '
        f'  (cd "$dir" && nohup {g16} < "${{name}}.gjf" > "${{name}}.log" 2>&1 &); '
        'done'
    )
    print(f"\n-- Launching single-point jobs on {host}:{remote_dir}")
    run(["ssh", host, submit_cmd], dry_run)
    print("\nJobs launched. Monitor with:\n  python scripts/dft/hpc_sync.py status")


def cmd_submit_ts(config: dict, dry_run: bool, mol: str | None, local_dir: Path) -> None:
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    print(
        "\nWARNING: TS (QST2/QST3) jobs carry an analytic CalcFC Hessian plus freq --\n"
        "         plausibly hours each, more expensive than the N-O scan jobs.\n"
        "         Submit ONE molecule at a time (both channels together is fine) --\n"
        "         do not batch across the test set. Confirm before every submission."
    )
    pattern = f"mol_{mol.zfill(3)}_*" if mol else "*"
    g16     = config["G16_PATH"]
    submit_cmd = (
        f'{_gauss_exports(config)} && '
        f'cd {remote_dir} && '
        f'for dir in {pattern}/; do '
        '  name="${dir%/}"; '
        '  (cd "$dir" && for gjf in "${name}"_ts*.gjf; do '
        '    case "$gjf" in *_irc.gjf) continue ;; esac; '
        '    [ -f "$gjf" ] || continue; '
        '    log="${gjf%.gjf}.log"; '
        f'    nohup {g16} < "$gjf" > "$log" 2>&1 & '
        '  done); '
        'done'
    )
    print(f"\n-- Launching TS (QST2/QST3) jobs on {host}:{remote_dir}")
    run(["ssh", host, submit_cmd], dry_run)
    print("\nTS jobs launched. Monitor with:\n  python scripts/dft/hpc_sync.py status")


def cmd_submit_irc(config: dict, dry_run: bool, mol: str | None, local_dir: Path) -> None:
    host       = config["HPC_HOST"]
    remote_dir = config["HPC_REMOTE_DIR"]
    print(
        "\nWARNING: IRC jobs read a TS job's .chk via %oldchk.\n"
        "         Only proceed once the corresponding TS job shows Normal termination\n"
        "         AND has passed manual verify_ts() review (exactly one imaginary\n"
        "         frequency, displacement vector matches the expected reaction\n"
        "         coordinate) -- do not IRC an unverified stationary point.\n"
        "         Check first: python scripts/dft/hpc_sync.py status"
    )
    pattern = f"mol_{mol.zfill(3)}_*" if mol else "*"
    g16     = config["G16_PATH"]
    submit_cmd = (
        f'{_gauss_exports(config)} && '
        f'cd {remote_dir} && '
        f'for dir in {pattern}/; do '
        '  name="${dir%/}"; '
        '  (cd "$dir" && for gjf in "${name}"_ts*_irc.gjf; do '
        '    [ -f "$gjf" ] || continue; '
        '    log="${gjf%.gjf}.log"; '
        f'    nohup {g16} < "$gjf" > "$log" 2>&1 & '
        '  done); '
        'done'
    )
    print(f"\n-- Launching IRC jobs on {host}:{remote_dir}")
    run(["ssh", host, submit_cmd], dry_run)
    print("\nIRC jobs launched. Monitor with:\n  python scripts/dft/hpc_sync.py status")


def cmd_download(config: dict, dry_run: bool, mol: str | None, local_dir: Path) -> None:
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
    parser.add_argument("--dry-run", action="store_true",
                        help="Print every command that would run, but do not execute it.")
    parser.add_argument("--mol", metavar="ID",
                        help="Target a single molecule by ID (e.g. --mol 002).")
    parser.add_argument("--dir", metavar="PATH", default=None,
                        help="Local job directory to sync (default: data/output/dft_opt).")

    sub = parser.add_subparsers(dest="command", metavar="COMMAND")
    sub.required = True
    sub.add_parser("upload",         help="Upload molecule directories to cluster")
    sub.add_parser("submit-opt",     help="Submit Stage 1 geometry-opt jobs ({name}_opt.gjf)")
    sub.add_parser("submit-nbo",     help="Submit Stage 2 NBO single-point jobs — AFTER Stage 1 finishes")
    sub.add_parser("submit-scan",    help="Submit Stage 3 N-O scan jobs — AFTER Stage 1 finishes")
    sub.add_parser("submit-scan-sp", help="Submit intermediate scan SP jobs — AFTER scan finishes")
    sub.add_parser("submit-sp",      help="Submit single-point NBO jobs (for dft_sp/)")
    sub.add_parser("submit-ts",      help="Submit QST2/QST3 TS jobs ({name}_ts*.gjf, excl. _irc)")
    sub.add_parser("submit-irc",     help="Submit IRC jobs — AFTER the matching TS job is verified")
    sub.add_parser("download",       help="Download *.log files from cluster")
    sub.add_parser("status",         help="Show running g16 processes on server")

    args = parser.parse_args()

    local_dir = Path(args.dir) if args.dir else DEFAULT_LOCAL_DFT_DIR
    if not local_dir.is_absolute():
        local_dir = PROJECT_ROOT / local_dir

    config = load_config()
    require_config(config)

    dispatch = {
        "upload":         cmd_upload,
        "submit-opt":     cmd_submit_opt,
        "submit-nbo":     cmd_submit_nbo,
        "submit-scan":    cmd_submit_scan,
        "submit-scan-sp": cmd_submit_scan_sp,
        "submit-sp":      cmd_submit_sp,
        "submit-ts":      cmd_submit_ts,
        "submit-irc":     cmd_submit_irc,
        "download":       cmd_download,
        "status":         cmd_status,
    }
    dispatch[args.command](config, args.dry_run, args.mol, local_dir)


if __name__ == "__main__":
    main()