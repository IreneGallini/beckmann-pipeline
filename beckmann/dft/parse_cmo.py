"""
Extract CMO (Canonical Molecular Orbital) data from NBO7 Gaussian logs and
compute the Lambda ("Frontier Dominance") and wCNmax ("CN-weighted Acceptor
Response") descriptors described in Notes.md.

Notes.md describes these verbally, not as a literal formula. The definitions
used here are the most direct reading of that text -- confirm against your
supervisor's source before treating them as final:

  Lambda     = squared leading-NBO-contribution coefficient (NBO's %-character
               "weight") of whichever virtual MO has the single largest such
               weight, among all virtual MOs from the LUMO up to LUMO+0.4 a.u.
               ("the maximum antibonding weights (w) found in the virtual
               manifold" -- not restricted to any one bond type).
  log_lambda = log10(Lambda).
  wCNmax     = the largest squared weight, in that same window, among
               contributions whose NBO label is specifically the developing
               C{ci}-N{ni} antibond (BD*, sigma or pi component) -- the bond
               common to both possible migrating groups (Notes.md's
               "nitrilium channel").

One row per (mol, stage, r_no) -- same grain as nbo_e2pert.csv.
Output: data/output/analysis/cmo_descriptors.csv
"""
import csv
import math
import re
from pathlib import Path

from beckmann.config import DATA_OUTPUT
from beckmann.dft.inputs import TEST_IDS
from beckmann.dft.parse_nbo import STAGES, r_no_before
from beckmann.dft.scan import oxime_atom_map_from_gjf

CMO_HEADER = "cmo: nbo analysis of canonical molecular orbitals"

MO_HEADER_RE = re.compile(r"MO\s+(\d+)\s+\((occ|vir)\):\s+orbital energy\s*=\s*(-?\d+\.\d+)\s*a\.u\.")
CONTRIB_RE   = re.compile(r"^\s*(-?\d+\.\d+)\*\[\s*\d+\]:\s+(.+?)\s*$")

FIELDS = [
    "mol", "stage", "r_no",
    "lambda", "log_lambda", "lambda_mo",
    "wcnmax", "wcnmax_mo",
    "n_virtual_mos_in_window",
]


def find_cmo_sections(lines: list[str]) -> list[int]:
    """Line indices right after each 'CMO: NBO Analysis...' header."""
    return [i + 1 for i, line in enumerate(lines) if CMO_HEADER in line.lower()]


def parse_cmo_table(lines: list[str], start: int) -> list[dict]:
    """Parse MO blocks (index, occ/vir, energy, leading NBO contributions)."""
    mos: list[dict] = []
    cur: dict | None = None
    j = start
    while j < len(lines):
        line = lines[j]
        header = MO_HEADER_RE.search(line)
        if header:
            cur = {
                "mo": int(header.group(1)),
                "kind": header.group(2),
                "energy": float(header.group(3)),
                "contribs": [],
            }
            mos.append(cur)
            j += 1
            continue
        contrib = CONTRIB_RE.match(line)
        if contrib and cur is not None:
            cur["contribs"].append((float(contrib.group(1)), contrib.group(2).strip()))
            j += 1
            continue
        if line.strip() == "":
            j += 1
            continue
        if mos:
            break
        j += 1
    return mos


def virtual_window(mo_table: list[dict], window_au: float = 0.4) -> list[dict]:
    """Virtual MOs from the LUMO up to LUMO + window_au."""
    vir = [m for m in mo_table if m["kind"] == "vir"]
    if not vir:
        return []
    lumo_e = vir[0]["energy"]
    return [m for m in vir if m["energy"] <= lumo_e + window_au]


def is_cn_antibond(label: str, ci: int, ni: int) -> bool:
    """True if label is an antibond (BD*) involving both atom C{ci} and N{ni}."""
    if "*" not in label:
        return False
    has_c = re.search(rf"C\s*{ci}(?!\d)", label) is not None
    has_n = re.search(rf"N\s*{ni}(?!\d)", label) is not None
    return has_c and has_n


def compute_descriptors(mo_table: list[dict], ci: int, ni: int) -> dict:
    window = virtual_window(mo_table)
    lambda_val = lambda_mo = wcn_val = wcn_mo = None

    for mo in window:
        if not mo["contribs"]:
            continue
        leading_coeff, _ = max(mo["contribs"], key=lambda t: abs(t[0]))
        w = leading_coeff ** 2
        if lambda_val is None or w > lambda_val:
            lambda_val, lambda_mo = w, mo["mo"]

        for coeff, label in mo["contribs"]:
            if is_cn_antibond(label, ci, ni):
                w_cn = coeff ** 2
                if wcn_val is None or w_cn > wcn_val:
                    wcn_val, wcn_mo = w_cn, mo["mo"]

    log_lambda = math.log10(lambda_val) if lambda_val else None
    return {
        "lambda": lambda_val, "log_lambda": log_lambda, "lambda_mo": lambda_mo,
        "wcnmax": wcn_val, "wcnmax_mo": wcn_mo,
        "n_virtual_mos_in_window": len(window),
    }


def parse_log(log_path: Path, ci: int, ni: int, oi: int) -> list[dict]:
    """Compute descriptors for every CMO table in a .log file, tagged with R(N-O)."""
    lines  = log_path.read_text().splitlines()
    starts = find_cmo_sections(lines)
    rows = []
    for start in starts:
        table = parse_cmo_table(lines, start)
        if not table:
            continue
        r_no = r_no_before(lines, start, ni, oi)
        row = compute_descriptors(table, ci, ni)
        row["r_no"] = round(r_no, 4) if r_no is not None else None
        rows.append(row)
    return rows


def collect_molecule(mol: str, mol_dir: Path) -> list[dict]:
    """Compute descriptors for all available stage logs of one molecule, e.g. 'mol_002_E'."""
    ci, ni, oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol}_opt.gjf")

    all_rows = []
    for stage in STAGES:
        log_path = mol_dir / f"{mol}_{stage}.log"
        if not log_path.exists():
            continue
        rows = parse_log(log_path, ci, ni, oi)

        # _scan.log has two CMO tables (start/end of scan) -- disambiguate by R(N-O) order.
        for point, row in enumerate(
            sorted(rows, key=lambda r: (r["r_no"] is None, r["r_no"])), start=1
        ):
            stage_label = f"{stage}_{point}" if len(rows) > 1 else stage
            all_rows.append({"mol": mol, "stage": stage_label, **row})

    return all_rows


def main() -> None:
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    out_path    = DATA_OUTPUT / "analysis" / "cmo_descriptors.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for mol_id in sorted(TEST_IDS):
        mol     = f"mol_{mol_id.zfill(3)}_E"
        mol_dir = dft_opt_dir / mol
        if not mol_dir.exists():
            print(f"-- {mol}: no directory, skipping")
            continue
        rows = collect_molecule(mol, mol_dir)
        print(f"-- {mol}: {len(rows)} stage points")
        for row in sorted(rows, key=lambda r: r["stage"]):
            print(
                f"     {row['stage']:<8} R(N-O)={row['r_no']}  "
                f"Lambda={row['lambda']}  logLambda={row['log_lambda']}  "
                f"wCNmax={row['wcnmax']}"
            )
        all_rows.extend(rows)

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\n{len(all_rows)} total rows -> {out_path}")


if __name__ == "__main__":
    main()
