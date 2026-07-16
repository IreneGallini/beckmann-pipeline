"""
Extract Wiberg bond index (NAO basis) data from NBO7 Gaussian logs for the two
migrating C-C bonds -- central-C to aryl-C (rearrangement channel) and
central-C to alkyl-C (fragmentation channel) -- across the rigid N-O scan.

BNDIDX is already in the $NBO keylist of every _nbo.gjf/_scan.gjf this project
generates, and NBO7 prints a full "Wiberg bond index matrix in the NAO basis"
table (an NxN matrix over all atoms, split into 9-column blocks by Gaussian)
at every point -- this module is pure extraction, no new Gaussian jobs needed.

c_aryl/c_alkyl come from beckmann.dft.descriptors.get_substituent_map() (fresh
RDKit aromaticity check, not any pre-computed CSV), same as parse_cmo.py.

One row per (mol, point), point = 'nbo' (R0) or 'scan_N' -- same convention as
the 'stage' column in cmo_descriptors.csv/nbo_e2pert.csv, resolvable through
descriptors.resolve_series().

Output: data/output/analysis/bond_order_scan.csv
"""
import csv
import re
from pathlib import Path

from beckmann.config import DATA_OUTPUT
from beckmann.dft.descriptors import get_substituent_map
from beckmann.dft.inputs import TEST_IDS, resolve_mol_name
from beckmann.dft.parse_nbo import STAGES, log_terminated_normally, r_no_before
from beckmann.dft.scan import oxime_atom_map_from_gjf

WIBERG_HEADER = "wiberg bond index matrix"

ATOM_HEADER_RE = re.compile(r"^\s*Atom((?:\s+\d+)+)\s*$")
ROW_RE         = re.compile(r"^\s*(\d+)\.\s+\S+\s+(.+)$")

FIELDS = ["mol", "point", "R", "bond_order_aryl", "bond_order_alkyl"]


def find_wiberg_sections(lines: list[str]) -> list[int]:
    """Line indices right after each 'Wiberg bond index matrix' header."""
    return [i + 1 for i, line in enumerate(lines) if WIBERG_HEADER in line.lower()]


def parse_wiberg_table(lines: list[str], start: int) -> dict[tuple[int, int], float]:
    """Parse one full NxN Wiberg matrix starting right after its header line,
    spanning however many 9-column blocks Gaussian splits it into. Returns
    {(row_atom, col_atom): bond_order} for every entry seen (matrix is
    symmetric, so callers can look up either atom order)."""
    matrix: dict[tuple[int, int], float] = {}
    col_atoms: list[int] | None = None
    j = start
    while j < len(lines):
        line = lines[j]
        header = ATOM_HEADER_RE.match(line)
        if header:
            col_atoms = [int(n) for n in header.group(1).split()]
            j += 1
            continue
        if line.strip().startswith("-"):
            j += 1
            continue
        row = ROW_RE.match(line)
        if row and col_atoms is not None:
            row_atom = int(row.group(1))
            values = [float(v) for v in row.group(2).split()]
            for col_atom, val in zip(col_atoms, values):
                matrix[(row_atom, col_atom)] = val
            j += 1
            continue
        if line.strip() == "":
            # a blank line separates column-blocks within the same table --
            # only a real stop if the next non-blank line isn't another block.
            k = j + 1
            while k < len(lines) and lines[k].strip() == "":
                k += 1
            if k < len(lines) and ATOM_HEADER_RE.match(lines[k]):
                j = k
                continue
            break
        break
    return matrix


def parse_log(log_path: Path, ci: int, ni: int, oi: int, c_aryl: int, c_alkyl: int) -> list[dict]:
    """Extract bond_order_aryl/bond_order_alkyl for every Wiberg table in a
    .log file, tagged with R(N-O).

    When multiple tables share the same R (Stable=Opt's pre-optimization seed
    pass vs. the final post-optimization pass, see parse_cmo.py/Notes.md for
    the same dedup issue), only the LAST table at that R is kept.
    """
    lines  = log_path.read_text().splitlines()
    starts = find_wiberg_sections(lines)
    row_by_r: dict[float | None, dict] = {}
    for start in starts:
        matrix = parse_wiberg_table(lines, start)
        if not matrix:
            continue
        r_no = r_no_before(lines, start, ni, oi)
        r_no = round(r_no, 4) if r_no is not None else None
        row_by_r[r_no] = {
            "r_no": r_no,
            "bond_order_aryl": matrix.get((ci, c_aryl), matrix.get((c_aryl, ci))),
            "bond_order_alkyl": matrix.get((ci, c_alkyl), matrix.get((c_alkyl, ci))),
        }  # last table at this R wins
    return list(row_by_r.values())


def collect_molecule(mol: str, mol_dir: Path, c_aryl: int, c_alkyl: int) -> list[dict]:
    """Extract bond orders from all available stage logs of one molecule, e.g.
    'mol_002_E'. If any present stage log didn't reach Normal termination,
    the whole molecule is skipped -- same policy as parse_nbo.py/parse_cmo.py.
    """
    ci, ni, oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol}_opt.gjf")

    bad_logs = [
        p.name for stage in STAGES
        if (p := mol_dir / f"{mol}_{stage}.log").exists() and not log_terminated_normally(p)
    ]
    if bad_logs:
        print(f"   -- {mol}: {', '.join(bad_logs)} did not reach Normal termination "
              f"-- skipping whole molecule (see JOB_ISSUES.md)")
        return []

    rows = []
    for stage in STAGES:
        log_path = mol_dir / f"{mol}_{stage}.log"
        if not log_path.exists():
            continue
        stage_rows = parse_log(log_path, ci, ni, oi, c_aryl, c_alkyl)

        # _scan.log has one Wiberg table per rigid-scan point -- disambiguate by R(N-O) order.
        for point, row in enumerate(
            sorted(stage_rows, key=lambda r: (r["r_no"] is None, r["r_no"])), start=1
        ):
            point_label = f"{stage}_{point}" if len(stage_rows) > 1 else stage
            rows.append({"mol": mol, "point": point_label, **row})

    return rows


def main() -> None:
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    out_path    = DATA_OUTPUT / "analysis" / "bond_order_scan.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for mol_id in sorted(TEST_IDS):
        mol = resolve_mol_name(mol_id, dft_opt_dir)
        if mol is None:
            print(f"-- mol_{mol_id.zfill(3)}: no directory, skipping")
            continue
        mol_dir = dft_opt_dir / mol
        subst = get_substituent_map(mol, mol_dir)
        rows = collect_molecule(mol, mol_dir, subst["c_aryl"], subst["c_alkyl"])
        print(f"-- {mol} (aryl=C{subst['c_aryl']}, alkyl=C{subst['c_alkyl']}): {len(rows)} points")
        for row in sorted(rows, key=lambda r: r["point"]):
            print(
                f"     {row['point']:<8} R(N-O)={row['r_no']}  "
                f"aryl={row['bond_order_aryl']}  alkyl={row['bond_order_alkyl']}"
            )
        all_rows.extend(rows)

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for row in all_rows:
            writer.writerow({
                "mol": row["mol"], "point": row["point"], "R": row["r_no"],
                "bond_order_aryl": row["bond_order_aryl"],
                "bond_order_alkyl": row["bond_order_alkyl"],
            })
    print(f"\n{len(all_rows)} total rows -> {out_path}")


if __name__ == "__main__":
    main()
