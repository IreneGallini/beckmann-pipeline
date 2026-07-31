"""
Extract NBO second-order perturbation (E2PERT) donor->acceptor tables from
Gaussian .log files and combine them into one CSV across the DFT test set.

Each stage log (_nbo.log, _scan.log) can contain one or more "Second Order
Perturbation Theory Analysis of Fock Matrix in NBO Basis" tables. _nbo.log
has exactly one (the R0 baseline). _scan.log has one per rigid-scan point --
4 for the standard 0.1 A architecture, more for a finer scan (e.g. 8 for
mol_006_E's 0.05 A finescan, see Notes.md) -- since the rigid-scan
architecture (RIGID_SCAN_MIGRATION.md) natively runs NBO7 at every point,
unlike the old internal-walk architecture this superseded (which only ran it
at 2 of 5 points, needing 'sp2'/'sp3'/'sp4' extracted single-point
workarounds -- no longer used by any current test molecule). Each table is
tagged with the N-O distance at that geometry, computed from the nearest
preceding "Standard orientation" block.

Output: data/output/analysis/nbo_e2pert.csv
"""
import csv
import re
from pathlib import Path

from beckmann_core.geometry import no_distance
from beckmann_nbo.config import DATA_OUTPUT
from beckmann_nbo.geometry import parse_standard_orientations
from beckmann_nbo.inputs import (
    ALL_IDS, STEP_SCAN_SOURCES, build_stage_relabel_map, relabel_rows,
    resolve_mol_name, step_scan_dir,
)
from beckmann_nbo.scan import oxime_indices_from_gjf

TABLE_HEADER = "second order perturbation theory analysis"
STAGES       = ["nbo", "scan"]

NORMAL_TERMINATION = "Normal termination of Gaussian 16"


def log_terminated_normally(log_path: Path) -> bool:
    """True if the log's final non-blank line reports normal termination.

    A crashed/non-converged job (e.g. mol_020_E's Stage 3 scan see
    JOB_ISSUES.md) can still print complete-looking NBO/E2PERT tables upstream
    of the crash; that data is computed on a geometry that never converged and
    must not be trusted just because it parses cleanly.
    """
    for line in reversed(log_path.read_text().splitlines()):
        if line.strip():
            return NORMAL_TERMINATION in line
    return False

# NBO7 and NBO 3.1 print this table in different layouts: NBO 3.1 separates
# donor/acceptor with ' / '; NBO7 just column-aligns them with no delimiter
# (e.g. '14. LP ( 1) O  2            49. BD*( 1) C  1- H 14      3.47 ...').
# Locate the two '<int>.' index markers positionally instead of relying on a
# fixed separator -- works for both formats.
IDX_RE    = re.compile(r"(\d+)\.\s+")
FLOATS_RE = re.compile(r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")

FIELDS = ["mol", "stage", "r_no", "donor", "acceptor", "e2_kcal", "de_au", "f_au"]


def find_table_starts(lines: list[str]) -> list[int]:
    """Line indices of the '===...===' separator that opens each E2PERT table."""
    starts = []
    for i, line in enumerate(lines):
        if TABLE_HEADER in line.lower():
            j = i
            while j < len(lines) and not re.match(r"^\s*=+\s*$", lines[j]):
                j += 1
            starts.append(j + 1)
    return starts


def parse_e2pert_line(line: str) -> dict | None:
    """Parse one donor/acceptor row, tolerant of both NBO 3.1 and NBO7 layouts."""
    fm = FLOATS_RE.search(line)
    if not fm:
        return None
    idxs = list(IDX_RE.finditer(line[:fm.start()]))
    if len(idxs) != 2:
        return None
    donor    = line[idxs[0].end():idxs[1].start()].strip(" /")
    acceptor = line[idxs[1].end():fm.start()].strip(" /")
    return {
        "donor":    " ".join(donor.split()),
        "acceptor": " ".join(acceptor.split()),
        "e2_kcal":  float(fm.group(1)),
        "de_au":    float(fm.group(2)),
        "f_au":     float(fm.group(3)),
    }


def parse_table_rows(lines: list[str], start: int) -> list[dict]:
    """Parse donor/acceptor rows starting just after a table's '===' line."""
    rows = []
    j = start
    while j < len(lines):
        line = lines[j]
        if line.strip() == "" or line.strip().lower().startswith("within unit"):
            j += 1
            continue
        row = parse_e2pert_line(line)
        if row is None:
            break
        rows.append(row)
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
    """Parse every E2PERT table in a .log file, each tagged with its R(N-O).

    When multiple tables share the same R (e.g. Stable=Opt -- used in every
    rigid-scan NBO block -- prints a pre-optimization seed-geometry pass and
    a separate post-optimization pass, both at the same frozen scan-point R;
    see Notes.md), only the LAST table at that R is kept. The seed isn't a
    converged/trustworthy geometry -- keeping both would silently double-count
    that point's E2PERT contributions (inflating Psi's K_anti/K_frag sums).
    """
    lines  = log_path.read_text().splitlines()
    starts = find_table_starts(lines)

    start_by_r: dict[float | None, int] = {}
    for start in starts:
        r_no = r_no_before(lines, start, ni, oi)
        r_no = round(r_no, 4) if r_no is not None else None
        start_by_r[r_no] = start  # last table at this R wins

    rows = []
    for r_no, start in start_by_r.items():
        for row in parse_table_rows(lines, start):
            row["r_no"] = r_no
            rows.append(row)
    return rows


def collect_stage(mol: str, mol_dir: Path, ni: int, oi: int, stage: str) -> list[dict] | None:
    """Extract one stage's rows, requiring only that stage's own log (if
    present) to have converged -- unlike collect_molecule(), which requires
    every stage present in mol_dir to be clean. Returns [] if the log doesn't
    exist, None if it exists but didn't reach Normal termination (a hard
    failure the caller must decide how to handle)."""
    log_path = mol_dir / f"{mol}_{stage}.log"
    if not log_path.exists():
        return []
    if not log_terminated_normally(log_path):
        return None

    table_rows_by_r = {}
    for row in parse_log(log_path, ni, oi):
        table_rows_by_r.setdefault(row["r_no"], []).append(row)

    rows = []
    # _scan.log has one table per rigid-scan point disambiguate by order.
    for point, (r_no, trows) in enumerate(sorted(
        table_rows_by_r.items(), key=lambda kv: (kv[0] is None, kv[0])
    ), start=1):
        stage_label = f"{stage}_{point}" if len(table_rows_by_r) > 1 else stage
        for row in trows:
            rows.append({"mol": mol, "stage": stage_label, "r_no": r_no, **row})
    return rows


def collect_molecule(mol: str, mol_dir: Path) -> list[dict]:
    """Parse all available stage logs for one molecule, e.g. 'mol_002_E'.

    If any present stage log didn't reach Normal termination, the whole
    molecule is skipped rather than partially included -- a partial series
    (missing points from a crashed job) isn't comparable to another
    molecule's complete series (see JOB_ISSUES.md).
    """
    ni, oi, _ = oxime_indices_from_gjf(mol_dir / f"{mol}_opt.gjf")

    all_rows = []
    for stage in STAGES:
        rows = collect_stage(mol, mol_dir, ni, oi, stage)
        if rows is None:
            print(f"   -- {mol}: {mol}_{stage}.log did not reach Normal termination "
                  f"-- skipping whole molecule (see JOB_ISSUES.md)")
            return []
        all_rows.extend(rows)
    return all_rows


def collect_molecule_stepscan(mol: str, mol_dir: Path) -> list[dict]:
    """For a molecule whose canonical Stage 3 scan crashed but has one or
    more successful step-size reruns (STEP_SCAN_SOURCES, see inputs.py):
    'nbo' rows come from the canonical mol_dir (Stage 2 equilibrium NBO
    succeeded independently of the Stage 3 crash), 'scan' rows come from
    each listed dft_opt_stepscan/ source, merged into one R-ordered series
    under the canonical mol name -- see build_stage_relabel_map()."""
    ni, oi, _ = oxime_indices_from_gjf(mol_dir / f"{mol}_opt.gjf")

    nbo_rows = collect_stage(mol, mol_dir, ni, oi, "nbo")
    if nbo_rows is None:
        print(f"   -- {mol}: {mol}_nbo.log did not reach Normal termination -- skipping")
        return []

    all_scan_rows = []
    for source in STEP_SCAN_SOURCES[mol]:
        source_dir = step_scan_dir() / source
        s_ni, s_oi, _ = oxime_indices_from_gjf(source_dir / f"{source}_opt.gjf")
        rows = collect_stage(source, source_dir, s_ni, s_oi, "scan")
        if rows is None:
            print(f"   -- {mol}: {source}_scan.log did not reach Normal termination -- skipping this source")
            continue
        all_scan_rows.extend(rows)

    relabel = build_stage_relabel_map({r["r_no"] for r in all_scan_rows})
    scan_rows = relabel_rows(all_scan_rows, mol, relabel)
    return nbo_rows + scan_rows


def main() -> None:
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    out_path    = DATA_OUTPUT / "analysis" / "nbo_e2pert.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for mol_id in sorted(ALL_IDS):
        mol = resolve_mol_name(mol_id, dft_opt_dir)
        if mol is None:
            print(f"-- mol_{mol_id.zfill(3)}: no directory, skipping")
            continue
        mol_dir = dft_opt_dir / mol
        if mol in STEP_SCAN_SOURCES:
            rows = collect_molecule_stepscan(mol, mol_dir)
        else:
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
