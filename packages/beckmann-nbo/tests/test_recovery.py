"""Tests for beckmann_nbo/recovery.py's describe_status() -- a read-only
report of a molecule's Stage 3 recovery state, shared internally by
plan_next_attempt() and the CLI's `status`/`report` commands. Synthetic
.log fixtures via tmp_path, matching test_log_diagnostics.py's pattern, so
these run in a fresh clone with no HPC output needed.
"""
from beckmann_nbo.recovery import RUNGS, describe_status

NORMAL_LOG = """\
 SCF Done:  E(RwB97XD) =  -557.303090     A.U. after   12 cycles
 Optimization completed.
    -- Stationary point found.
 Normal termination of Gaussian 16 at Fri Jul 17 17:07:58 2026.
"""

OSCILLATING_LOG = "\n".join(
    [" Maximum Force            " + ("0.025012     0.000450     NO " if i % 2 == 0
                                       else "0.015634     0.000450     NO ")
     for i in range(30)]
) + """
 Error termination via Lnk1e in /opt/g16/l9999.exe at Fri Jul 10 00:31:25 2026.
Segmentation fault (core dumped)
"""


def test_describe_status_clean_scan(tmp_path):
    base_dir = tmp_path / "mol_999_E"
    base_dir.mkdir()
    (base_dir / "mol_999_E_scan.log").write_text(NORMAL_LOG)

    status = describe_status("mol_999_E", base_dir, stepscan_root=tmp_path / "stepscan")
    assert status == {"needed_recovery": False, "rung": None, "resolved": True}


def test_describe_status_oscillating_no_rung_tried_yet(tmp_path):
    base_dir = tmp_path / "mol_999_E"
    base_dir.mkdir()
    (base_dir / "mol_999_E_scan.log").write_text(OSCILLATING_LOG)
    stepscan_root = tmp_path / "stepscan"

    status = describe_status("mol_999_E", base_dir, stepscan_root=stepscan_root)
    assert status == {"needed_recovery": True, "rung": None, "resolved": False}


def test_describe_status_resolved_at_step04(tmp_path):
    base_dir = tmp_path / "mol_999_E"
    base_dir.mkdir()
    (base_dir / "mol_999_E_scan.log").write_text(OSCILLATING_LOG)
    stepscan_root = tmp_path / "stepscan"

    # calcfc and step07 rungs were tried and are still failing...
    for rung in ("calcfc", "step07"):
        d = stepscan_root / f"mol_999_E_{rung}"
        d.mkdir(parents=True)
        (d / f"mol_999_E_{rung}_scan.log").write_text(OSCILLATING_LOG)
    # ...but step04 succeeded.
    d = stepscan_root / "mol_999_E_step04"
    d.mkdir(parents=True)
    (d / "mol_999_E_step04_scan.log").write_text(NORMAL_LOG)

    status = describe_status("mol_999_E", base_dir, stepscan_root=stepscan_root)
    assert status == {"needed_recovery": True, "rung": "step04", "resolved": True}


def test_describe_status_ladder_exhausted(tmp_path):
    base_dir = tmp_path / "mol_999_E"
    base_dir.mkdir()
    (base_dir / "mol_999_E_scan.log").write_text(OSCILLATING_LOG)
    stepscan_root = tmp_path / "stepscan"

    for rung_name, _ in RUNGS:
        d = stepscan_root / f"mol_999_E_{rung_name}"
        d.mkdir(parents=True)
        (d / f"mol_999_E_{rung_name}_scan.log").write_text(OSCILLATING_LOG)

    status = describe_status("mol_999_E", base_dir, stepscan_root=stepscan_root)
    assert status == {"needed_recovery": True, "rung": None, "resolved": False}
