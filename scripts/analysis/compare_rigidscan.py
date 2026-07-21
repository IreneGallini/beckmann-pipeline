"""
Compare the rigid-scan architecture (data/output/dft_opt_rigidscan/) against
main's already-committed numbers (data/output/analysis/channel_descriptors.csv)
for the same 4 molecules, mol_002_E/006_E/020_E/021_E.

Standalone, throwaway comparison for the rigid-scan-architecture branch
experiment -- not part of the regular pipeline. Reuses parse_nbo.parse_log()/
parse_cmo.parse_log() (both file-path-agnostic) directly against the new
_scan.log files rather than duplicating their table-finding logic.

mol_006_E_rigidscan's point 4 (R0+0.4) crashed (see JOB_ISSUES.md) and was
excluded rather than fixed -- that molecule only has 4 of 5 points here.
"""
import csv
from pathlib import Path

from beckmann.config import DATA_OUTPUT
from beckmann.dft import parse_cmo, parse_nbo
from beckmann.dft.descriptors import compute_psi_row, get_substituent_map
from beckmann.dft.scan import oxime_atom_map_from_gjf

OLD_DIR = DATA_OUTPUT / "dft_opt"
NEW_DIR = DATA_OUTPUT / "dft_opt_rigidscan"
ANALYSIS_DIR = DATA_OUTPUT / "analysis"

MOLS = ["mol_002_E", "mol_006_E", "mol_020_E", "mol_021_E"]
DESCRIPTORS = ["psi", "log_lambda", "wcnmax", "w17max", "w78max"]


def _read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def group_by_r(rows: list[dict]) -> dict[float, list[dict]]:
    by_r: dict[float, list[dict]] = {}
    for row in rows:
        by_r.setdefault(row["r_no"], []).append(row)
    return by_r


def new_architecture_series(mol: str) -> list[dict]:
    """R0 (reused from the old, unchanged Stage 2 _nbo.log) + the rigid-scan
    points (new _scan.log), one row per R with psi/log_lambda/wcnmax/w17max/w78max."""
    old_mol_dir = OLD_DIR / mol
    new_name = f"{mol}_rigidscan"
    new_mol_dir = NEW_DIR / new_name

    ci, ni, oi, _ = oxime_atom_map_from_gjf(old_mol_dir / f"{mol}_opt.gjf")
    subst = get_substituent_map(mol, old_mol_dir)
    c_aryl, c_alkyl = subst["c_aryl"], subst["c_alkyl"]

    rows = []

    # R0 baseline -- Stage 2 methodology unchanged, reuse as-is.
    nbo_log = old_mol_dir / f"{mol}_nbo.log"
    e2pert_r0 = group_by_r(parse_nbo.parse_log(nbo_log, ni, oi))
    cmo_r0 = parse_cmo.parse_log(nbo_log, ci, ni, oi, c_aryl, c_alkyl)
    if cmo_r0:
        r0 = cmo_r0[0]["r_no"]
        psi_row = compute_psi_row(e2pert_r0.get(r0, []), ci, ni, oi, c_aryl, c_alkyl)
        rows.append({"r_no": r0, "stage": "nbo", **psi_row, **cmo_r0[0]})

    # Rigid-scan points -- new architecture.
    scan_log = new_mol_dir / f"{new_name}_scan.log"
    e2pert_by_r = group_by_r(parse_nbo.parse_log(scan_log, ni, oi))
    cmo_rows = parse_cmo.parse_log(scan_log, ci, ni, oi, c_aryl, c_alkyl)
    for point, cmo_row in enumerate(sorted(cmo_rows, key=lambda r: r["r_no"]), start=1):
        r_no = cmo_row["r_no"]
        psi_row = compute_psi_row(e2pert_by_r.get(r_no, []), ci, ni, oi, c_aryl, c_alkyl)
        rows.append({"r_no": r_no, "stage": f"rigid_pt{point}", **psi_row, **cmo_row})

    return rows


def main() -> None:
    old_channel_rows = _read_csv(ANALYSIS_DIR / "channel_descriptors.csv")

    for mol in MOLS:
        print(f"\n{'='*100}\n{mol}\n{'='*100}")
        old_rows = sorted(
            (r for r in old_channel_rows if r["mol"] == mol),
            key=lambda r: float(r["r_no"]),
        )
        new_rows = new_architecture_series(mol)

        print(f"{'R(N-O)':>8}  {'':>10} | {'psi':>10} {'log_lambda':>12} {'wcnmax':>10} {'w17max':>10} {'w78max':>10}")
        print("-" * 100)
        print("-- OLD (main, chained-walk / post-hoc extraction) --")
        for row in old_rows:
            print(
                f"{float(row['r_no']):>8.4f}  {row['stage']:>10} | "
                f"{float(row['psi']):>10.4f} {float(row['log_lambda']):>12.4f} "
                f"{float(row['wcnmax']):>10.4f} {float(row['w17max']):>10.4f} {float(row['w78max']):>10.4f}"
            )
        print("-- NEW (rigid-scan, independent per-point) --")
        for row in new_rows:
            print(
                f"{row['r_no']:>8.4f}  {row['stage']:>10} | "
                f"{row['psi']:>10.4f} {row['log_lambda']:>12.4f} "
                f"{row['wcnmax']:>10.4f} {row['w17max']:>10.4f} {row['w78max']:>10.4f}"
            )


if __name__ == "__main__":
    main()
