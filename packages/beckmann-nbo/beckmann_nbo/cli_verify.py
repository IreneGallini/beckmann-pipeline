"""`beckmann-nbo verify` -- preflight checks against the configured cluster.

Runs three checks over SSH: connectivity, that G16_PATH is executable, and
that NBO_WRAPPER_DIR's gaunbo7/gaunbo6 are executable. The third check is
the one that catches the failure mode where a job reaches "Normal
termination" but silently used the bundled NBO 3.1 instead of NBO7 (no CMO
section in the log) -- otherwise only discoverable after the fact by
grepping a finished .log for "NBO 7.0".

Uses its own subprocess call (not beckmann_nbo.hpc.run) because hpc.run()
exits the process on a nonzero return code -- verify needs to run and
report all three checks, not stop at the first failure.
"""
import subprocess
import sys

from beckmann_nbo.hpc import load_config, require_config


def _ssh_check(host: str, remote_cmd: str, dry_run: bool) -> bool:
    cmd = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=10", host, remote_cmd]
    print(f"$ {' '.join(cmd)}")
    if dry_run:
        print("  [dry-run: not executed]")
        return True
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        if result.stderr.strip():
            print(f"  {result.stderr.strip()}")
        return False
    return True


def _fix_block(config: dict) -> str:
    nboexe_dir = "<NBOEXE's directory>"
    if config.get("NBOEXE"):
        nboexe_dir = config["NBOEXE"].rsplit("/", 1)[0]
    wrapper_dir = config.get("NBO_WRAPPER_DIR", "<NBO_WRAPPER_DIR>")
    host = config.get("HPC_HOST", "<HPC_HOST>")
    return (
        f'  ssh {host} "mkdir -p {wrapper_dir} && \\\n'
        f'    cp {nboexe_dir}/gaunbo7 {nboexe_dir}/gaunbo6 {wrapper_dir}/ && \\\n'
        f'    chmod +x {wrapper_dir}/gaunbo7 {wrapper_dir}/gaunbo6"'
    )


def cmd_verify(args) -> None:
    config = load_config()
    require_config(config)
    dry_run = args.dry_run
    host = config["HPC_HOST"]

    print(f"\n-- Check 1/3: SSH connectivity to {host}")
    ok_conn = _ssh_check(host, "echo ok", dry_run)
    print("  OK" if ok_conn else "  FAILED -- check HPC_HOST and SSH key auth (ssh-copy-id)")

    print(f"\n-- Check 2/3: G16_PATH executable ({config['G16_PATH']})")
    ok_g16 = _ssh_check(host, f"test -x {config['G16_PATH']}", dry_run) if ok_conn else False
    print("  OK" if ok_g16 else "  FAILED -- G16_PATH does not exist or is not executable on the cluster")

    print("\n-- Check 3/3: NBO7 wrapper (gaunbo7/gaunbo6) executable")
    wrapper_dir = config.get("NBO_WRAPPER_DIR")
    if not wrapper_dir:
        ok_nbo = False
        print("  FAILED -- NBO_WRAPPER_DIR not set in .env")
    elif ok_conn:
        ok_nbo = _ssh_check(
            host, f"test -x {wrapper_dir}/gaunbo7 && test -x {wrapper_dir}/gaunbo6", dry_run
        )
        print("  OK" if ok_nbo else "  FAILED -- gaunbo7/gaunbo6 missing or not executable")
    else:
        ok_nbo = False

    if not ok_nbo:
        print(
            "\nWithout an executable gaunbo7/gaunbo6 on PATH, Gaussian's pop=nbo7read\n"
            "silently falls back to the bundled NBO 3.1 (no CMO support) -- jobs will\n"
            "still reach 'Normal termination' but produce logs with no CMO/wCNmax data.\n"
            "Fix (copy the vendor binaries somewhere you own and chmod +x):\n" + _fix_block(config)
        )

    if not (ok_conn and ok_g16 and ok_nbo):
        sys.exit(1)
    print("\nAll checks passed.")
