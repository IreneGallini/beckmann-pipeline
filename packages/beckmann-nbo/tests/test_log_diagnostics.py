"""
Tests for beckmann_nbo/log_diagnostics.py -- the classifier that turns
JOB_ISSUES.md's manual crash-diagnosis playbook into code. Uses synthetic
.log excerpts (via tmp_path) rather than real downloaded logs, so these run
in a fresh clone with no HPC output needed.
"""
from beckmann_nbo.log_diagnostics import FailureCategory, classify_log

NORMAL_LOG = """\
 SCF Done:  E(RwB97XD) =  -557.303090     A.U. after   12 cycles
 Optimization completed.
    -- Stationary point found.
 Elapsed time:       0 days  0 hours  6 minutes  3.4 seconds.
 Normal termination of Gaussian 16 at Fri Jul 17 17:07:58 2026.
"""

# Same signature documented for mol_020_E's original crash and this week's
# mol_003_E_step07 rerun: a clean 2-state Maximum Force alternation with no
# trend, then Error termination + Segmentation fault.
OSCILLATING_LOG = "\n".join(
    [" Maximum Force            " + ("0.025012     0.000450     NO " if i % 2 == 0
                                       else "0.015634     0.000450     NO ")
     for i in range(30)]
) + """
 Error termination via Lnk1e in /opt/g16/l9999.exe at Fri Jul 10 00:31:25 2026.
Segmentation fault (core dumped)
"""

SEGFAULT_NO_FORCE_LOG = """\
 Initial guess read from the checkpoint file:
Segmentation fault (core dumped)
"""

SLOW_CONVERGENCE_LOG = "\n".join(
    f" Maximum Force            {0.05 / (i + 1):.6f}     0.000450     NO "
    for i in range(15)
) + """
 Error termination via Lnk1e in /opt/g16/l9999.exe at Fri Jul 10 00:31:25 2026.
Segmentation fault (core dumped)
"""


def _write(tmp_path, name: str, content: str):
    path = tmp_path / name
    path.write_text(content)
    return path


def test_classify_log_normal_termination(tmp_path):
    log_path = _write(tmp_path, "mol_test_scan.log", NORMAL_LOG)
    diag = classify_log(log_path, stage="scan")
    assert diag.category == FailureCategory.NORMAL
    assert diag.max_force_trace == []


def test_classify_log_oscillating_degeneracy(tmp_path):
    log_path = _write(tmp_path, "mol_test_scan.log", OSCILLATING_LOG)
    diag = classify_log(log_path, stage="scan")
    assert diag.category == FailureCategory.OSCILLATING_DEGENERACY
    assert len(diag.max_force_trace) == 30


def test_classify_log_segfault_no_force_data(tmp_path):
    log_path = _write(tmp_path, "mol_test_scan.log", SEGFAULT_NO_FORCE_LOG)
    diag = classify_log(log_path, stage="scan")
    assert diag.category == FailureCategory.SEGFAULT
    assert diag.max_force_trace == []


def test_classify_log_slow_convergence(tmp_path):
    log_path = _write(tmp_path, "mol_test_scan.log", SLOW_CONVERGENCE_LOG)
    diag = classify_log(log_path, stage="scan")
    assert diag.category == FailureCategory.SLOW_CONVERGENCE


def test_classify_log_stage_label(tmp_path):
    """stage defaults to the file stem when not given explicitly."""
    log_path = _write(tmp_path, "mol_test_scan.log", NORMAL_LOG)
    diag = classify_log(log_path)
    assert diag.stage == "mol_test_scan"
