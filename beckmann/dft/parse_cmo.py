"""
Extract CMO (Canonical Molecular Orbital) data from NBO7 Gaussian logs and
compute the channel-resolved wX^max descriptors and Lambda ("Frontier
Dominance") from "Ring Size and Substituent Effects in the Beckmann
Rearrangement" (Sections 2.2-2.4). See Notes.md for the full derivation.

  wX^max     = for a target antibond X (BD*, sigma or pi component -- whichever
               gives the larger squared coefficient wins), scan every virtual
               MO from the LUMO up to LUMO+0.4 a.u., take X's coefficient in
               that MO's CMO expansion (0 if X doesn't appear), square it, take
               the max across the window.
  w17max     = wX^max for X = BD*(C{ci}-C{c_aryl})  (rearrangement channel)
  w78max     = wX^max for X = BD*(C{ci}-C{c_alkyl}) (fragmentation channel)
  wcnmax     = wX^max for X = BD*(C{ci}-N{ni})      (nitrilium/routing channel)
  Lambda     = max(w78max) / max(max(w17max), 1e-3) -- fragmentation-channel
               dominance over rearrangement-channel dominance, NOT an
               unrestricted max over the whole window (that was the bug in an
               earlier version of this module -- see Notes.md).
  log_lambda = log10(Lambda).

c_aryl/c_alkyl come from beckmann.dft.descriptors.get_substituent_map() (fresh
RDKit aromaticity check, not any pre-computed CSV).

One row per (mol, stage, r_no) -- same grain as nbo_e2pert.csv.
Output: data/output/analysis/cmo_descriptors.csv
"""
import csv
import math
import re
from pathlib import Path

from beckmann.config import DATA_OUTPUT
from beckmann.dft.descriptors import get_substituent_map
from beckmann.dft.inputs import TEST_IDS
from beckmann.dft.parse_nbo import STAGES, r_no_before
from beckmann.dft.scan import oxime_atom_map_from_gjf

CMO_HEADER = "cmo: nbo analysis of canonical molecular orbitals"

MO_HEADER_RE = re.compile(r"MO\s+(\d+)\s+\((occ|vir)\):\s+orbital energy\s*=\s*(-?\d+\.\d+)\s*a\.u\.")
CONTRIB_RE   = re.compile(r"^\s*(-?\d+\.\d+)\*\[\s*\d+\]:\s+(.+?)\s*$")

LAMBDA_FLOOR = 1e-3

FIELDS = [
    "mol", "stage", "r_no",
    "lambda", "log_lambda",
    "w17max", "w17max_mo",
    "w78max", "w78max_mo",
    "wcnmax", "wcnmax_mo",
    "max_leading_weight", "max_leading_weight_mo",
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


def is_target_antibond(label: str, a: int, b: int) -> bool:
    """True if label is an antibond (BD*) involving atoms with numbers a and b."""
    if "*" not in label:
        return False
    has_a = re.search(rf"C\s*{a}(?!\d)", label) is not None
    has_b = re.search(rf"[CN]\s*{b}(?!\d)", label) is not None
    return has_a and has_b


def max_weight_for_target(window: list[dict], a: int, b: int) -> tuple[float | None, int | None]:
    """Max squared coefficient, across the window, of the BD* antibond between atoms a and b."""
    best_val = best_mo = None
    for mo in window:
        for coeff, label in mo["contribs"]:
            if is_target_antibond(label, a, b):
                w = coeff ** 2
                if best_val is None or w > best_val:
                    best_val, best_mo = w, mo["mo"]
    return best_val, best_mo


def compute_descriptors(mo_table: list[dict], ci: int, ni: int, c_aryl: int, c_alkyl: int) -> dict:
    window = virtual_window(mo_table)

    max_leading_val = max_leading_mo = None
    for mo in window:
        if not mo["contribs"]:
            continue
        leading_coeff, _ = max(mo["contribs"], key=lambda t: abs(t[0]))
        w = leading_coeff ** 2
        if max_leading_val is None or w > max_leading_val:
            max_leading_val, max_leading_mo = w, mo["mo"]

    w17max, w17max_mo = max_weight_for_target(window, ci, c_aryl)
    w78max, w78max_mo = max_weight_for_target(window, ci, c_alkyl)
    wcnmax, wcnmax_mo = max_weight_for_target(window, ci, ni)

    lambda_val  = (w78max or 0.0) / max(w17max or 0.0, LAMBDA_FLOOR)
    log_lambda  = math.log10(lambda_val) if lambda_val > 0 else None

    return {
        "lambda": lambda_val, "log_lambda": log_lambda,
        "w17max": w17max, "w17max_mo": w17max_mo,
        "w78max": w78max, "w78max_mo": w78max_mo,
        "wcnmax": wcnmax, "wcnmax_mo": wcnmax_mo,
        "max_leading_weight": max_leading_val, "max_leading_weight_mo": max_leading_mo,
        "n_virtual_mos_in_window": len(window),
    }


def parse_log(log_path: Path, ci: int, ni: int, oi: int, c_aryl: int, c_alkyl: int) -> list[dict]:
    """Compute descriptors for every CMO table in a .log file, tagged with R(N-O)."""
    lines  = log_path.read_text().splitlines()
    starts = find_cmo_sections(lines)
    rows = []
    for start in starts:
        table = parse_cmo_table(lines, start)
        if not table:
            continue
        r_no = r_no_before(lines, start, ni, oi)
        row = compute_descriptors(table, ci, ni, c_aryl, c_alkyl)
        row["r_no"] = round(r_no, 4) if r_no is not None else None
        rows.append(row)
    return rows


def collect_molecule(mol: str, mol_dir: Path, c_aryl: int, c_alkyl: int) -> list[dict]:
    """Compute descriptors for all available stage logs of one molecule, e.g. 'mol_002_E'."""
    ci, ni, oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol}_opt.gjf")

    all_rows = []
    for stage in STAGES:
        log_path = mol_dir / f"{mol}_{stage}.log"
        if not log_path.exists():
            continue
        rows = parse_log(log_path, ci, ni, oi, c_aryl, c_alkyl)

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
        subst = get_substituent_map(mol, mol_dir)
        rows = collect_molecule(mol, mol_dir, subst["c_aryl"], subst["c_alkyl"])
        print(f"-- {mol} (aryl=C{subst['c_aryl']}, alkyl=C{subst['c_alkyl']}): {len(rows)} stage points")
        for row in sorted(rows, key=lambda r: r["stage"]):
            print(
                f"     {row['stage']:<8} R(N-O)={row['r_no']}  "
                f"Lambda={row['lambda']:.4f}  logLambda={row['log_lambda']}  "
                f"w17max={row['w17max']}  w78max={row['w78max']}  wCNmax={row['wcnmax']}"
            )
        all_rows.extend(rows)

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\n{len(all_rows)} total rows -> {out_path}")


if __name__ == "__main__":
    main()
