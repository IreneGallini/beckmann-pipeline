"""
Automated escalation ladder for a crashed/oscillating Stage 3 scan: detect
via beckmann.dft.log_diagnostics, then automatically generate, upload, and
submit the next remediation rung on Citadel -- turning this week's manual
recovery (CalcFC test, then step07/step04 step-shift reruns, see
JOB_ISSUES.md) into a repeatable, unattended process.

Escalation order (per the user's explicit choice, even though this week's
data showed the step-shift rungs are individually more reliable than
CalcFC): CalcFC first, then step=0.07, then step=0.04. Three rungs, hard cap
-- if all three still show OSCILLATING_DEGENERACY, stop and flag for human
review rather than searching further parameter space unattended.

Submission here is FULLY AUTOMATIC (the user's explicit choice, despite
Citadel having no job scheduler) -- the one safeguard this module does
enforce is an in-flight check (_running_dirs()) so re-running this on a cron
schedule can't double-submit a rung that's already running from a previous
invocation.

This module intentionally does NOT touch a molecule that isn't currently
OSCILLATING_DEGENERACY (a SLOW_CONVERGENCE/NOISY_TRENDING/SEGFAULT/UNKNOWN
failure gets no automated suggestion at all) -- CalcFC/step-shift are only
validated against the ring-pucker oscillation signature; applying them to a
different, unvalidated failure mode would be exactly the "automating a guess
dressed up as a fix" risk flagged for the HPC product tier earlier.
"""
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from beckmann.dft.hpc import _gauss_exports, run
from beckmann.dft.inputs import prepare_scan_rigid, step_scan_dir
from beckmann.dft.log_diagnostics import FailureCategory, classify_log

RUNGS: list[tuple[str, dict]] = [
    ("calcfc", {"calcfc": True}),
    ("step07", {"step": 0.07, "calcfc": False}),
    ("step04", {"step": 0.04, "calcfc": False}),
]


@dataclass
class RemediationPlan:
    mol: str            # canonical name, e.g. 'mol_020_E'
    rung: str           # 'calcfc' | 'step07' | 'step04'
    name: str           # e.g. 'mol_020_E_calcfc'
    directory: Path
    kwargs: dict


def plan_next_attempt(mol: str, base_dir: Path, stepscan_root: Path | None = None) -> RemediationPlan | None:
    """Decide the next escalation rung for one molecule, or None if there's
    nothing automated to do (not oscillating, already resolved, or every
    rung already exhausted -- all three read as "leave it, needs a human").
    """
    stepscan_root = stepscan_root or step_scan_dir()
    diag = classify_log(base_dir / f"{mol}_scan.log", stage="scan")
    if diag.category != FailureCategory.OSCILLATING_DEGENERACY:
        return None

    for rung_name, kwargs in RUNGS:
        candidate_dir = stepscan_root / f"{mol}_{rung_name}"
        candidate_log = candidate_dir / f"{mol}_{rung_name}_scan.log"
        if candidate_log.exists():
            rung_diag = classify_log(candidate_log, stage="scan")
            if rung_diag.category == FailureCategory.NORMAL:
                return None  # already resolved by an earlier run of this ladder
            continue  # this rung was tried and is still failing -- escalate
        return RemediationPlan(
            mol=mol, rung=rung_name, name=f"{mol}_{rung_name}",
            directory=candidate_dir, kwargs=kwargs,
        )
    return None  # every rung exhausted, all still failing -- needs a human


def _running_dirs(config: dict) -> set[str]:
    """Directory basenames that currently have a live g16 process on
    Citadel, resolved via each PID's cwd -- 'ps aux' alone can't tell which
    molecule/rung a process belongs to, since Gaussian reads its input via
    shell redirection rather than a visible argv filename (confirmed
    manually earlier this project via readlink /proc/$pid/cwd)."""
    host = config["HPC_HOST"]
    g16  = config["G16_PATH"]
    pids_result = subprocess.run(["ssh", host, f"pgrep -f {g16}"], capture_output=True, text=True)
    pids = pids_result.stdout.split()
    if not pids:
        return set()
    cwd_cmd = "; ".join(f"readlink /proc/{pid}/cwd 2>/dev/null" for pid in pids)
    cwd_result = subprocess.run(["ssh", host, cwd_cmd], capture_output=True, text=True)
    return {Path(line).name for line in cwd_result.stdout.splitlines() if line.strip()}


def execute_plan(plan: RemediationPlan, base_dir: Path, config: dict, dry_run: bool = False) -> bool:
    """Generate the candidate .gjf (copying Stage 1 opt.gjf/opt.log from
    base_dir, same pattern used manually for this week's step07/step04
    reruns), then upload and submit it -- unless a job for this exact
    candidate is already running remotely. Returns True if a job was (or
    would be, under --dry-run) submitted, False if skipped as in-flight."""
    if plan.name in _running_dirs(config):
        print(f"-- {plan.name}: already running on Citadel, skipping (in-flight check)")
        return False

    host, remote_dir, g16 = config["HPC_HOST"], config["HPC_REMOTE_DIR"], config["G16_PATH"]
    print(f"\n-- Escalating {plan.mol}: rung '{plan.rung}' ({plan.name})")

    if dry_run:
        print(f"  [dry-run: would generate {plan.directory}/{plan.name}_scan.gjf "
              f"(kwargs={plan.kwargs}), then upload + submit]")
        return True

    plan.directory.mkdir(parents=True, exist_ok=True)
    shutil.copy(base_dir / f"{plan.mol}_opt.gjf", plan.directory / f"{plan.name}_opt.gjf")
    shutil.copy(base_dir / f"{plan.mol}_opt.log", plan.directory / f"{plan.name}_opt.log")
    prepare_scan_rigid(plan.directory, plan.name, **plan.kwargs)

    run(["ssh", host, f"mkdir -p {remote_dir}"], dry_run)
    run(["scp", "-r", str(plan.directory), f"{host}:{remote_dir}/"], dry_run)
    submit_cmd = (
        f'{_gauss_exports(config)} && cd {remote_dir} && '
        f'(cd "{plan.name}" && nohup {g16} < "{plan.name}_scan.gjf" > "{plan.name}_scan.log" 2>&1 &)'
    )
    run(["ssh", host, submit_cmd], dry_run)
    return True


def run_auto_recovery(mols: list[str], dft_opt_dir: Path, config: dict, dry_run: bool = False) -> None:
    """One pass over the given canonical mol names (e.g. ['mol_020_E']):
    for each, decide and (if anything to do) execute the next escalation
    rung. Meant to be re-run periodically (cron or by hand), not run as a
    persistent daemon -- see scripts/dft/auto_recover.py."""
    for mol in mols:
        base_dir = dft_opt_dir / mol
        if not (base_dir / f"{mol}_scan.log").exists():
            print(f"-- {mol}: no {mol}_scan.log downloaded yet, skipping")
            continue
        plan = plan_next_attempt(mol, base_dir)
        if plan is None:
            print(f"-- {mol}: nothing to do (not oscillating, already resolved, or ladder exhausted)")
            continue
        execute_plan(plan, base_dir, config, dry_run)
