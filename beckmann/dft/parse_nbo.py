"""
Extract NBO second-order perturbation (E2PERT) donor->acceptor tables from
Gaussian .log files and combine them into one CSV across the DFT test set.

Each stage log (_nbo.log, _scan.log, _sp2.log, _sp3.log, _sp4.log) can contain
one or more "Second Order Perturbation Theory Analysis of Fock Matrix in NBO
Basis" tables. _scan.log has two (start and end of the relaxed scan); the rest
have exactly one. Each table is tagged with the N-O distance at that geometry,
computed from the nearest preceding "Standard orientation" block.

Output: data/output/analysis/nbo_e2pert.csv
"""
import csv
import re
from pathlib import Path

from beckmann.config import DATA_OUTPUT
from beckmann.dft.inputs import TEST_IDS
from beckmann.dft.scan import no_distance, oxime_indices_from_gjf, parse_standard_orientations

TABLE_HEADER = "Second Order Perturbation Theory Analysis"
STAGES       = ["nbo", "scan", "sp2", "sp3", "sp4"]

ROW_RE = re.compile(
    r"^\s*\d+\.\s+(?P<donor>.+?)\s*/\s*\d+\.\s+(?P<acceptor>.+?)\s+"
    r"(?P<e2>-?\d+\.\d+)\s+(?P<de>-?\d+\.\d+)\s+(?P<f>-?\d+\.\d+)\s*$"
)

FIELDS = ["mol", "stage", "r_no", "donor", "acceptor", "e2_kcal", "de_au", "f_au"]


def find_table_starts(lines: list[str]) -> list[int]:
    """Line indices of the '===...===' separator that opens each E2PERT table."""
    starts = []
    for i, line in enumerate(lines):
        if TABLE_HEADER in line:
            j = i
            while j < len(lines) and not re.match(r"^\s*=+\s*$", lines[j]):
                j += 1
            starts.append(j + 1)
    return starts


def parse_table_rows(lines: list[str], start: int) -> list[dict]:
    """Parse donor/acceptor rows starting just after a table's '===' line."""
    rows = []
    j = start
    while j < len(lines):
        line = lines[j]
        if line.strip() == "" or line.strip().lower().startswith("within unit"):
            j += 1
            continue
        match = ROW_RE.match(line)
        if not match:
            break
        rows.append({
            "donor":    " ".join(match.group("donor").split()),
            "acceptor": " ".join(match.group("acceptor").split()),
            "e2_kcal":  float(match.group("e2")),
            "de_au":    float(match.group("de")),
            "f_au":     float(match.group("f")),
        })
        j += 1
    return rows


def r_no_before(lines: list[str], idx: int, ni: int, oi: int) -> float | None:
    """N-O distance from the last 'Standard orientation' block before line idx."""
    so_blocks = [(i, atoms) for i, atoms in parse_standard_orientations(lines) if i < idx]
    if not so_blocks:
        return None
    _, atoms = max(so_blocks, key=lambda x: x[0])
    return no_distance(atoms, ni, oi)


def parse_log(log_path: Path, ni: int, oi: int) -> list[dict]:
    """Parse every E2PERT table in a .log file, each tagged with its R(N-O)."""
    lines  = log_path.read_text().splitlines()
    starts = find_table_starts(lines)
    rows   = []
    for start in starts:
        r_no = r_no_before(lines, start, ni, oi)
        for row in parse_table_rows(lines, start):
            row["r_no"] = round(r_no, 4) if r_no is not None else None
            rows.append(row)
    return rows


def collect_molecule(mol: str, mol_dir: Path) -> list[dict]:
    """Parse all available stage logs for one molecule, e.g. 'mol_002_E'."""
    ni, oi, _ = oxime_indices_from_gjf(mol_dir / f"{mol}_opt.gjf")

    all_rows = []
    for stage in STAGES:
        log_path = mol_dir / f"{mol}_{stage}.log"
        if not log_path.exists():
            continue
        table_rows_by_r = {}
        for row in parse_log(log_path, ni, oi):
            table_rows_by_r.setdefault(row["r_no"], []).append(row)

        # _scan.log has two tables (start/end of scan) — disambiguate by order.
        for point, (r_no, rows) in enumerate(sorted(
            table_rows_by_r.items(), key=lambda kv: (kv[0] is None, kv[0])
        ), start=1):
            stage_label = f"{stage}_{point}" if len(table_rows_by_r) > 1 else stage
            for row in rows:
                all_rows.append({"mol": mol, "stage": stage_label, "r_no": r_no, **row})

    return all_rows


def main() -> None:
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    out_path    = DATA_OUTPUT / "analysis" / "nbo_e2pert.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for mol_id in sorted(TEST_IDS):
        mol     = f"mol_{mol_id.zfill(3)}_E"
        mol_dir = dft_opt_dir / mol
        if not mol_dir.exists():
            print(f"-- {mol}: no directory, skipping")
            continue
        rows = collect_molecule(mol, mol_dir)
        print(f"-- {mol}: {len(rows)} E2PERT rows")
        all_rows.extend(rows)

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\n{len(all_rows)} total rows -> {out_path}")


if __name__ == "__main__":
    main()
