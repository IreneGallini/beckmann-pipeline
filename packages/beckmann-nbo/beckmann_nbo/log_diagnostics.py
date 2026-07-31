"""
Turn the manual crash-diagnosis playbook in data/output/dft_opt/JOB_ISSUES.md
into code: classify a Gaussian .log's failure mode (if any) from its
'Maximum Force' trend, rather than a human eyeballing `grep`.

Built on top of log_terminated_normally() (parse_nbo.py) rather than
duplicating it -- that function already answers "did this converge," this
module only adds "if not, what kind of failure is it."

This is deliberately a DETECTION tool, not a remediation one -- see
beckmann/dft/recovery.py for what (if anything) automatically happens with a
classification. A classifier that force-fits every failure into a known
bucket is worse than useless (it hides genuinely novel failure modes), so
UNKNOWN is a legitimate, expected result, not a bug.
"""
import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

from beckmann_nbo.parse_nbo import log_terminated_normally

MAX_FORCE_RE   = re.compile(r"Maximum Force\s+(-?\d+\.\d+)")
STEP_NUMBER_RE = re.compile(r"Step number\s+(\d+)\s+out of a maximum of\s+(\d+)")


class FailureCategory(Enum):
    NORMAL                 = "normal"
    OSCILLATING_DEGENERACY = "oscillating_degeneracy"
    SLOW_CONVERGENCE       = "slow_convergence"
    NOISY_TRENDING         = "noisy_trending"
    SEGFAULT               = "segfault"
    ERROR_TERMINATION      = "error_termination"
    UNKNOWN                = "unknown"


@dataclass
class ScanPointDiagnosis:
    stage: str
    category: FailureCategory
    max_force_trace: list[float] = field(default_factory=list)
    step_count: int = 0
    step_cap: int | None = None


TAIL_WINDOW = 30  # steps -- the settled terminal behavior, not the initial transient


def _classify_trend(values: list[float]) -> FailureCategory:
    """Numeric version of JOB_ISSUES.md's three failure-trend shapes:
    - clean 2-state alternation, no trend  -> OSCILLATING_DEGENERACY
      (a genuine double-well degeneracy, commonly ring pucker -- more steps
      will not fix this)
    - clearly decaying                     -> SLOW_CONVERGENCE
      (monotonic but slow -- usually just needs more MaxCycles)
    - noisy but net-decreasing             -> NOISY_TRENDING
      (often self-resolves given enough steps)
    - anything else                        -> UNKNOWN
      (a judgment call even for a human -- don't force-fit it)

    Classified on the TAIL_WINDOW most recent steps, not the whole trace --
    a crashing point's optimization typically starts with a normal-looking
    descent before locking into a repeating 2-state bounce (see mol_003_E's
    step07 rerun: 0.026 -> 0.011 -> 0.002 ... before settling into
    0.025/0.016 alternation for 100+ steps), and the early transient would
    otherwise mask the terminal pattern that actually explains the crash.
    """
    tail = values[-TAIL_WINDOW:] if len(values) > TAIL_WINDOW else values
    if len(tail) < 4:
        return FailureCategory.UNKNOWN

    evens, odds = tail[0::2], tail[1::2]
    if len(evens) >= 2 and len(odds) >= 2:
        evens_spread = max(evens) - min(evens)
        odds_spread  = max(odds) - min(odds)
        band_gap     = abs(sum(evens) / len(evens) - sum(odds) / len(odds))
        # Each parity-band internally tight (values within it barely move)
        # but the two bands sit clearly apart -- the "same two values
        # repeating" signature from mol_020_E's original oscillation entry.
        if band_gap > 0.005 and band_gap > 3 * max(evens_spread, odds_spread, 1e-9):
            return FailureCategory.OSCILLATING_DEGENERACY

    n = max(1, len(tail) // 3)
    first_mean = sum(tail[:n]) / n
    last_mean  = sum(tail[-n:]) / n
    if last_mean < first_mean * 0.7:
        return FailureCategory.SLOW_CONVERGENCE
    if last_mean < first_mean:
        return FailureCategory.NOISY_TRENDING
    return FailureCategory.UNKNOWN


def classify_log(log_path: Path, stage: str | None = None) -> ScanPointDiagnosis:
    """Classify one .log file's outcome. stage defaults to the file's stem
    (e.g. 'mol_020_E_scan' -> use an explicit stage name via classify_scan()
    instead when you have one -- this default is only a fallback for
    standalone calls)."""
    stage = stage or log_path.stem

    if log_terminated_normally(log_path):
        return ScanPointDiagnosis(stage=stage, category=FailureCategory.NORMAL)

    text  = log_path.read_text()
    lines = text.splitlines()

    # A crash always sits after the last successfully-converged point (every
    # earlier point already printed its own "Stationary point found"); only
    # the trailing segment reflects the actual failing point's trajectory,
    # not a blend of every point's Maximum Force history in the file.
    last_good_idx = 0
    for i, line in enumerate(lines):
        if "Stationary point found" in line:
            last_good_idx = i
    tail = "\n".join(lines[last_good_idx:])

    max_force_trace = [float(m) for m in MAX_FORCE_RE.findall(tail)]
    step_matches = STEP_NUMBER_RE.findall(tail)
    step_count = len(step_matches)
    step_cap = int(step_matches[-1][1]) if step_matches else None

    if max_force_trace:
        category = _classify_trend(max_force_trace)
        if category == FailureCategory.UNKNOWN and ("Segmentation fault" in text or "Error termination" in text):
            category = FailureCategory.ERROR_TERMINATION
    elif "Segmentation fault" in text:
        category = FailureCategory.SEGFAULT
    elif "Error termination" in text:
        category = FailureCategory.ERROR_TERMINATION
    else:
        category = FailureCategory.UNKNOWN

    return ScanPointDiagnosis(
        stage=stage, category=category, max_force_trace=max_force_trace,
        step_count=step_count, step_cap=step_cap,
    )


def classify_scan(mol_dir: Path, mol: str, stages: list[str] = ("opt", "nbo", "scan")) -> list[ScanPointDiagnosis]:
    """Classify every present stage log for one molecule directory. Stages
    default to the three canonical Gaussian stages this project generates
    ({mol}_opt.log, {mol}_nbo.log, {mol}_scan.log) -- pass a different list
    for a step-scan side-experiment directory (e.g. just ["scan"], since
    those only ever hold {name}_scan.log)."""
    results = []
    for stage in stages:
        log_path = mol_dir / f"{mol}_{stage}.log"
        if not log_path.exists():
            continue
        results.append(classify_log(log_path, stage=stage))
    return results
