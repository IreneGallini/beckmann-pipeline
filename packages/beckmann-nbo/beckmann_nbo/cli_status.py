"""`beckmann-nbo status` -- per-stage job status, recovery caveats, and (once
Stage 3 is clean) a live wCNmax R/F prediction, for one or all molecules
under a dft_opt/-shaped directory.
"""
from pathlib import Path

from beckmann_core.wcnmax_rule import find_wcnmax_minimum, predict_from_wcnmax
from beckmann_nbo.descriptors import get_substituent_map
from beckmann_nbo.hpc import (
    DEFAULT_LOCAL_DFT_DIR,
    cmd_download,
    cmd_status,
    load_config,
    mol_dirs,
    require_config,
)
from beckmann_nbo.inputs import STEP_SCAN_SOURCES
from beckmann_nbo.log_diagnostics import FailureCategory, classify_scan
from beckmann_nbo.parse_cmo import collect_molecule, collect_molecule_stepscan
from beckmann_nbo.recovery import describe_status


def _predict(mol: str, mol_dir: Path) -> None:
    try:
        c_map = get_substituent_map(mol, mol_dir)
    except ValueError as e:
        print(f"    prediction unavailable: {e}")
        return

    if mol in STEP_SCAN_SOURCES:
        _, channel_rows = collect_molecule_stepscan(mol, mol_dir, c_map["c_aryl"], c_map["c_alkyl"])
    else:
        _, channel_rows = collect_molecule(mol, mol_dir, c_map["c_aryl"], c_map["c_alkyl"])

    minimum = find_wcnmax_minimum(mol, channel_rows)
    cn_points = {
        r["stage"] for r in channel_rows
        if r["channel"] == "cn" and r["weight"] not in (None, "", "None")
    }
    if minimum is None and len(cn_points) < 3:
        print(f"    prediction unavailable: only {len(cn_points)} resolved scan point(s) (need >= 3)")
        return

    prediction = predict_from_wcnmax(minimum)
    if minimum is not None:
        print(
            f"    prediction: {prediction}  "
            f"(interior wCNmax minimum at R(N-O)={minimum['R_star']:.3f} A, "
            f"w={minimum['w_star']:.4f}, MO{minimum['MO_index']})"
        )
    else:
        print(f"    prediction: {prediction}  (no interior wCNmax minimum found)")


def cmd_status_command(args) -> None:
    config = load_config()
    require_config(config)
    local_dir = Path(args.dir) if args.dir else DEFAULT_LOCAL_DFT_DIR

    if not args.no_download:
        cmd_download(config, args.dry_run, args.mol, local_dir)

    cmd_status(config, args.dry_run, args.mol, local_dir)

    # Classification below reads only already-local logs (no network) --
    # still worth running under --dry-run, which only means "don't touch
    # the cluster," not "don't report local state."
    for d in mol_dirs(local_dir, args.mol):
        mol = d.name
        print(f"\n{mol}")
        diagnoses = classify_scan(d, mol)
        if not diagnoses:
            print("  no stage logs downloaded yet")
            continue

        by_stage = {diag.stage: diag for diag in diagnoses}
        for stage in ("opt", "nbo", "scan"):
            if stage not in by_stage:
                continue
            diag = by_stage[stage]
            flag = "" if diag.category == FailureCategory.NORMAL else "  <-- needs attention"
            print(f"  Stage {stage:<5} {diag.category.value:<24}{flag}")

        scan_diag = by_stage.get("scan")
        clean_scan = scan_diag is not None and scan_diag.category == FailureCategory.NORMAL

        if scan_diag is not None and scan_diag.category == FailureCategory.OSCILLATING_DEGENERACY:
            status = describe_status(mol, d)
            if status["resolved"]:
                print(
                    f"  Stage scan: recovered via '{status['rung']}' rerun -- "
                    f"NOT a clean single-run result, see JOB_ISSUES.md"
                )
                clean_scan = True
            elif status["needed_recovery"]:
                print("  Stage scan: still oscillating, recovery ladder in progress or exhausted")

        if clean_scan or mol in STEP_SCAN_SOURCES:
            _predict(mol, d)
