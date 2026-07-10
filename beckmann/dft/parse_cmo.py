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
  Lambda     = w78max / w17max -- fragmentation-channel dominance over
               rearrangement-channel dominance, NOT an unrestricted max over
               the whole window (that was the bug in an earlier version of
               this module -- see Notes.md). Left undefined (None) when
               w17max wasn't found anywhere in the LUMO..LUMO+0.4 window --
               previously floored to w78max/1e-3 in that case, which silently
               produced huge Lambda values that were a division-floor
               artifact, not a real ratio (see JOB_ISSUES.md-adjacent
               discussion; caught because it inflated Lambda for two
               substrates where the rearrangement-channel antibond simply
               sits outside the window, not because the ratio was genuinely
               extreme).
  log_lambda = log10(Lambda), also None when Lambda is None.

c_aryl/c_alkyl come from beckmann.dft.descriptors.get_substituent_map() (fresh
RDKit aromaticity check, not any pre-computed CSV).

One row per (mol, stage, r_no) -- same grain as nbo_e2pert.csv.
Output: data/output/analysis/cmo_descriptors.csv (summary: w17max/w78max/wcnmax/Lambda)
        data/output/analysis/cmo_channel_extraction.csv (per-channel detail: which MO
        carried the max weight at each geometry, its orbital energy, and the raw signed
        coefficient before squaring -- the wX^max summary columns above are a max over
        this data, but which MO achieves it can shift identity between scan points
        (canonical MOs are energy-ordered and can change character along the scan), so
        that intermediate detail is kept here rather than discarded. This is what a
        later avoided-crossing check (small eigenvalue gap + character exchange between
        two nearby virtual MOs near the wCNmax extremum) would need -- not implemented
        yet, this module only preserves the data for it.
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

FIELDS = [
    "mol", "stage", "r_no",
    "lambda", "log_lambda",
    "w17max", "w17max_mo",
    "w78max", "w78max_mo",
    "wcnmax", "wcnmax_mo",
    "max_leading_weight", "max_leading_weight_mo",
    "n_virtual_mos_in_window",
]

# channel name -> which two atoms bound the target antibond BD*(C{ci}-{atom for b})
CHANNEL_TARGETS = {
    "cn": lambda ci, ni, c_aryl, c_alkyl: (ci, ni),
    "17": lambda ci, ni, c_aryl, c_alkyl: (ci, c_aryl),
    "78": lambda ci, ni, c_aryl, c_alkyl: (ci, c_alkyl),
}

EXTRACTION_FIELDS = [
    "mol", "stage", "channel", "R_NO", "MO_index", "epsilon_i_star", "coefficient", "weight",
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


def max_weight_for_target(window: list[dict], a: int, b: int) -> tuple[float | None, int | None, float | None, float | None]:
    """Max squared coefficient, across the window, of the BD* antibond between atoms a and b.

    Returns (weight, mo_index, epsilon_i_star, coefficient) -- the orbital energy and
    signed coefficient of the MO that achieves the max, not just the final weight.
    """
    best_val = best_mo = best_epsilon = best_coeff = None
    for mo in window:
        for coeff, label in mo["contribs"]:
            if is_target_antibond(label, a, b):
                w = coeff ** 2
                if best_val is None or w > best_val:
                    best_val, best_mo, best_epsilon, best_coeff = w, mo["mo"], mo["energy"], coeff
    return best_val, best_mo, best_epsilon, best_coeff


def compute_channel_weights(window: list[dict], ci: int, ni: int, c_aryl: int, c_alkyl: int) -> dict:
    """Raw (weight, mo_index, epsilon_i_star, coefficient) per channel: 'cn', '17', '78'."""
    return {
        name: max_weight_for_target(window, *target(ci, ni, c_aryl, c_alkyl))
        for name, target in CHANNEL_TARGETS.items()
    }


def compute_descriptors(mo_table: list[dict], ci: int, ni: int, c_aryl: int, c_alkyl: int) -> tuple[dict, dict]:
    """Returns (summary_dict, channels) -- channels is the raw per-channel detail
    (weight, mo_index, epsilon_i_star, coefficient) that the summary's wX^max/wX^max_mo
    columns are themselves derived from, kept around for the extraction table."""
    window = virtual_window(mo_table)
    channels = compute_channel_weights(window, ci, ni, c_aryl, c_alkyl)

    max_leading_val = max_leading_mo = None
    for mo in window:
        if not mo["contribs"]:
            continue
        leading_coeff, _ = max(mo["contribs"], key=lambda t: abs(t[0]))
        w = leading_coeff ** 2
        if max_leading_val is None or w > max_leading_val:
            max_leading_val, max_leading_mo = w, mo["mo"]

    wcnmax, wcnmax_mo, _, _ = channels["cn"]
    w17max, w17max_mo, _, _ = channels["17"]
    w78max, w78max_mo, _, _ = channels["78"]

    lambda_val  = w78max / w17max if (w17max is not None and w17max > 0 and w78max is not None) else None
    log_lambda  = math.log10(lambda_val) if lambda_val is not None and lambda_val > 0 else None

    summary = {
        "lambda": lambda_val, "log_lambda": log_lambda,
        "w17max": w17max, "w17max_mo": w17max_mo,
        "w78max": w78max, "w78max_mo": w78max_mo,
        "wcnmax": wcnmax, "wcnmax_mo": wcnmax_mo,
        "max_leading_weight": max_leading_val, "max_leading_weight_mo": max_leading_mo,
        "n_virtual_mos_in_window": len(window),
    }
    return summary, channels


def parse_log(log_path: Path, ci: int, ni: int, oi: int, c_aryl: int, c_alkyl: int) -> list[dict]:
    """Compute descriptors for every CMO table in a .log file, tagged with R(N-O).

    Each returned row carries the summary fields plus '_channels' (raw per-channel
    weight/mo/epsilon/coefficient detail, popped off by collect_molecule()).
    """
    lines  = log_path.read_text().splitlines()
    starts = find_cmo_sections(lines)
    rows = []
    for start in starts:
        table = parse_cmo_table(lines, start)
        if not table:
            continue
        r_no = r_no_before(lines, start, ni, oi)
        summary, channels = compute_descriptors(table, ci, ni, c_aryl, c_alkyl)
        summary["r_no"] = round(r_no, 4) if r_no is not None else None
        summary["_channels"] = channels
        rows.append(summary)
    return rows


def collect_molecule(mol: str, mol_dir: Path, c_aryl: int, c_alkyl: int) -> tuple[list[dict], list[dict]]:
    """Compute descriptors for all available stage logs of one molecule, e.g. 'mol_002_E'.

    Returns (summary_rows, channel_extraction_rows) -- the latter is the per-geometry,
    per-channel (cn/17/78) detail table: which MO carried the max weight, its orbital
    energy, and the signed coefficient before squaring.
    """
    ci, ni, oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol}_opt.gjf")

    summary_rows = []
    channel_rows = []
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
            channels = row.pop("_channels")
            summary_rows.append({"mol": mol, "stage": stage_label, **row})
            for channel_name, (weight, mo_index, epsilon, coeff) in channels.items():
                channel_rows.append({
                    "mol": mol, "stage": stage_label, "channel": channel_name,
                    "R_NO": row["r_no"], "MO_index": mo_index,
                    "epsilon_i_star": epsilon, "coefficient": coeff, "weight": weight,
                })

    return summary_rows, channel_rows


def main() -> None:
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    out_path        = DATA_OUTPUT / "analysis" / "cmo_descriptors.csv"
    extraction_path = DATA_OUTPUT / "analysis" / "cmo_channel_extraction.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    all_rows = []
    all_channel_rows = []
    for mol_id in sorted(TEST_IDS):
        mol     = f"mol_{mol_id.zfill(3)}_E"
        mol_dir = dft_opt_dir / mol
        if not mol_dir.exists():
            print(f"-- {mol}: no directory, skipping")
            continue
        subst = get_substituent_map(mol, mol_dir)
        rows, channel_rows = collect_molecule(mol, mol_dir, subst["c_aryl"], subst["c_alkyl"])
        print(f"-- {mol} (aryl=C{subst['c_aryl']}, alkyl=C{subst['c_alkyl']}): {len(rows)} stage points")
        for row in sorted(rows, key=lambda r: r["stage"]):
            lambda_str = f"{row['lambda']:.4f}" if row["lambda"] is not None else "None"
            print(
                f"     {row['stage']:<8} R(N-O)={row['r_no']}  "
                f"Lambda={lambda_str}  logLambda={row['log_lambda']}  "
                f"w17max={row['w17max']}  w78max={row['w78max']}  wCNmax={row['wcnmax']}"
            )
        all_rows.extend(rows)
        all_channel_rows.extend(channel_rows)

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\n{len(all_rows)} total rows -> {out_path}")

    with open(extraction_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=EXTRACTION_FIELDS)
        writer.writeheader()
        writer.writerows(all_channel_rows)
    print(f"{len(all_channel_rows)} total rows -> {extraction_path}")


if __name__ == "__main__":
    main()
