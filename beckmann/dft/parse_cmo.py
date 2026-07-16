"""
Extract CMO (Canonical Molecular Orbital) data from NBO7 Gaussian logs and
compute the channel-resolved wX^max descriptors and Lambda ("Frontier
Dominance") from "Ring Size and Substituent Effects in the Beckmann
Rearrangement" (Sections 2.2-2.4). See Notes.md for the full derivation.

  wX^max     = for a target antibond X (BD*, sigma or pi component -- whichever
               gives the larger squared coefficient wins), scan EVERY virtual MO
               (the full manifold, no energy cutoff), take X's coefficient in
               that MO's CMO expansion (0 if X doesn't appear), square it, take
               the max across all of them.
  w17max     = wX^max for X = BD*(C{ci}-C{c_aryl})  (rearrangement channel)
  w78max     = wX^max for X = BD*(C{ci}-C{c_alkyl}) (fragmentation channel)
  wcnmax     = wX^max for X = BD*(C{ci}-N{ni})      (nitrilium/routing channel)
  Lambda     = w78max / w17max -- fragmentation-channel dominance over
               rearrangement-channel dominance, NOT an unrestricted max over
               the whole window (that was the bug in an earlier version of
               this module -- see Notes.md). Left undefined (None) when
               w17max wasn't found anywhere in the full virtual manifold.

  Fixed bug (see Notes.md): this used to cap the search at LUMO+0.4 a.u. and
  silently return None if the target antibond's peak mixing fell outside that
  window. For mol_006_E/mol_021_E, BD*(C{ci}-C{c_aryl}) is a real NBO (confirmed
  in the NBOSUM table, occupancy ~0.032) that mixes strongly (25-38% coefficient,
  well above any print threshold) into virtual MOs sitting at LUMO+0.48 to
  LUMO+0.60 a.u. just past the old cutoff. That's not print-threshold noise,
  it's the antibond's character showing up in a higher-lying virtual MO than the
  paper's nominal window assumed, consistent with the project's own CN-handoff /
  avoided-crossing hypothesis (frontier character migrating between canonical
  MOs). The search now covers the whole virtual manifold unconditionally, and
  each wX^max carries a companion `*_delta_lumo` (energy above the LUMO, a.u.)
  and `*_in_window` (whether that fell inside the original 0.4 a.u. window) so
  it's visible whenever a descriptor is being driven by a MO outside the
  nominal frontier region.
  log_lambda = log10(Lambda), also None when Lambda is None.

c_aryl/c_alkyl come from beckmann.dft.descriptors.get_substituent_map() (fresh
RDKit aromaticity check, not any pre-computed CSV).

One row per (mol, stage, r_no) same grain as nbo_e2pert.csv.
Output: data/output/analysis/cmo_descriptors.csv (summary: w17max/w78max/wcnmax/Lambda)
        data/output/analysis/cmo_channel_extraction.csv (per-channel detail: which MO
        carried the max weight at each geometry, its orbital energy, and the raw signed
        coefficient before squaring the wX^max summary columns above are a max over
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
from beckmann.dft.inputs import TEST_IDS, resolve_mol_name
from beckmann.dft.parse_nbo import STAGES, log_terminated_normally, r_no_before
from beckmann.dft.scan import oxime_atom_map_from_gjf

CMO_HEADER = "cmo: nbo analysis of canonical molecular orbitals"

MO_HEADER_RE = re.compile(r"MO\s+(\d+)\s+\((occ|vir)\):\s+orbital energy\s*=\s*(-?\d+\.\d+)\s*a\.u\.")
CONTRIB_RE   = re.compile(r"^\s*(-?\d+\.\d+)\*\[\s*\d+\]:\s+(.+?)\s*$")

FIELDS = [
    "mol", "stage", "r_no",
    "lambda", "log_lambda",
    "w17max", "w17max_mo", "w17max_delta_lumo", "w17max_in_window",
    "w78max", "w78max_mo", "w78max_delta_lumo", "w78max_in_window",
    "wcnmax", "wcnmax_mo", "wcnmax_delta_lumo", "wcnmax_in_window",
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
    "delta_lumo", "in_window",
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


def max_weight_for_target(
    vir_mos: list[dict], a: int, b: int, lumo_e: float | None = None, window_au: float = 0.4,
) -> tuple[float | None, int | None, float | None, float | None, float | None, bool | None]:
    """Max squared coefficient, across the FULL virtual manifold (no energy cutoff), of the
    BD* antibond between atoms a and b.

    `vir_mos` should be every virtual MO in the log, not just a windowed subset -- capping the
    search at a fixed energy window was the bug: a target antibond's peak mixing can land above
    LUMO+0.4 a.u. for some substrates (see module docstring), and silently returning None there
    was indistinguishable from the antibond genuinely not existing.

    Returns (weight, mo_index, epsilon_i_star, coefficient, delta_lumo, in_window) -- the orbital
    energy and signed coefficient of the MO that achieves the max, plus how far above the LUMO
    it sits and whether that's inside the original 0.4 a.u. window (both None if lumo_e isn't
    given, or if the antibond was never found).
    """
    best_val = best_mo = best_epsilon = best_coeff = None
    for mo in vir_mos:
        for coeff, label in mo["contribs"]:
            if is_target_antibond(label, a, b):
                w = coeff ** 2
                if best_val is None or w > best_val:
                    best_val, best_mo, best_epsilon, best_coeff = w, mo["mo"], mo["energy"], coeff
    if best_val is None or lumo_e is None:
        return best_val, best_mo, best_epsilon, best_coeff, None, None
    delta_lumo = best_epsilon - lumo_e
    return best_val, best_mo, best_epsilon, best_coeff, delta_lumo, delta_lumo <= window_au + 1e-9


def compute_channel_weights(
    vir_mos: list[dict], ci: int, ni: int, c_aryl: int, c_alkyl: int, lumo_e: float | None = None,
) -> dict:
    """Raw (weight, mo_index, epsilon_i_star, coefficient, delta_lumo, in_window) per channel:
    'cn', '17', '78' -- searched across the full virtual manifold, see max_weight_for_target."""
    return {
        name: max_weight_for_target(vir_mos, *target(ci, ni, c_aryl, c_alkyl), lumo_e)
        for name, target in CHANNEL_TARGETS.items()
    }


def compute_descriptors(mo_table: list[dict], ci: int, ni: int, c_aryl: int, c_alkyl: int) -> tuple[dict, dict]:
    """Returns (summary_dict, channels) channels is the raw per-channel detail
    (weight, mo_index, epsilon_i_star, coefficient) that the summary's wX^max/wX^max_mo
    columns are themselves derived from, kept around for the extraction table."""
    window  = virtual_window(mo_table)
    vir_all = [m for m in mo_table if m["kind"] == "vir"]
    lumo_e  = vir_all[0]["energy"] if vir_all else None

    channels = compute_channel_weights(vir_all, ci, ni, c_aryl, c_alkyl, lumo_e)

    # max_leading_weight/n_virtual_mos_in_window are still scoped to the nominal
    # frontier window -- a different, generic "dominant character near the LUMO"
    # descriptor, not the per-channel target-antibond search that was buggy.
    max_leading_val = max_leading_mo = None
    for mo in window:
        if not mo["contribs"]:
            continue
        leading_coeff, _ = max(mo["contribs"], key=lambda t: abs(t[0]))
        w = leading_coeff ** 2
        if max_leading_val is None or w > max_leading_val:
            max_leading_val, max_leading_mo = w, mo["mo"]

    wcnmax, wcnmax_mo, _, _, wcnmax_delta, wcnmax_in_window = channels["cn"]
    w17max, w17max_mo, _, _, w17max_delta, w17max_in_window = channels["17"]
    w78max, w78max_mo, _, _, w78max_delta, w78max_in_window = channels["78"]

    lambda_val  = w78max / w17max if (w17max is not None and w17max > 0 and w78max is not None) else None
    log_lambda  = math.log10(lambda_val) if lambda_val is not None and lambda_val > 0 else None

    summary = {
        "lambda": lambda_val, "log_lambda": log_lambda,
        "w17max": w17max, "w17max_mo": w17max_mo,
        "w17max_delta_lumo": round(w17max_delta, 5) if w17max_delta is not None else None,
        "w17max_in_window": w17max_in_window,
        "w78max": w78max, "w78max_mo": w78max_mo,
        "w78max_delta_lumo": round(w78max_delta, 5) if w78max_delta is not None else None,
        "w78max_in_window": w78max_in_window,
        "wcnmax": wcnmax, "wcnmax_mo": wcnmax_mo,
        "wcnmax_delta_lumo": round(wcnmax_delta, 5) if wcnmax_delta is not None else None,
        "wcnmax_in_window": wcnmax_in_window,
        "max_leading_weight": max_leading_val, "max_leading_weight_mo": max_leading_mo,
        "n_virtual_mos_in_window": len(window),
    }
    return summary, channels


def parse_log(log_path: Path, ci: int, ni: int, oi: int, c_aryl: int, c_alkyl: int) -> list[dict]:
    """Compute descriptors for every CMO table in a .log file, tagged with R(N-O).

    Each returned row carries the summary fields plus '_channels' (raw per-channel
    weight/mo/epsilon/coefficient detail, popped off by collect_molecule()).

    When multiple tables share the same R (e.g. Stable=Opt's pre-optimization
    seed-geometry pass vs. the final post-optimization pass, both at the same
    frozen scan-point R -- see Notes.md), only the LAST table at that R is
    kept -- the seed isn't a converged/trustworthy geometry.
    """
    lines  = log_path.read_text().splitlines()
    starts = find_cmo_sections(lines)
    row_by_r: dict[float | None, dict] = {}
    for start in starts:
        table = parse_cmo_table(lines, start)
        if not table:
            continue
        r_no = r_no_before(lines, start, ni, oi)
        summary, channels = compute_descriptors(table, ci, ni, c_aryl, c_alkyl)
        summary["r_no"] = round(r_no, 4) if r_no is not None else None
        summary["_channels"] = channels
        row_by_r[summary["r_no"]] = summary  # last table at this R wins
    return list(row_by_r.values())


def collect_molecule(mol: str, mol_dir: Path, c_aryl: int, c_alkyl: int) -> tuple[list[dict], list[dict]]:
    """Compute descriptors for all available stage logs of one molecule, e.g. 'mol_002_E'.

    Returns (summary_rows, channel_extraction_rows) -- the latter is the per-geometry,
    per-channel (cn/17/78) detail table: which MO carried the max weight, its orbital
    energy, and the signed coefficient before squaring.

    If any present stage log didn't reach Normal termination, the whole molecule is
    skipped rather than partially included -- see parse_nbo.collect_molecule for why.
    """
    ci, ni, oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol}_opt.gjf")

    bad_logs = [
        p.name for stage in STAGES
        if (p := mol_dir / f"{mol}_{stage}.log").exists() and not log_terminated_normally(p)
    ]
    if bad_logs:
        print(f"   -- {mol}: {', '.join(bad_logs)} did not reach Normal termination "
              f"-- skipping whole molecule (see JOB_ISSUES.md)")
        return [], []

    summary_rows = []
    channel_rows = []
    for stage in STAGES:
        log_path = mol_dir / f"{mol}_{stage}.log"
        if not log_path.exists():
            continue
        rows = parse_log(log_path, ci, ni, oi, c_aryl, c_alkyl)

        # _scan.log has one CMO table per rigid-scan point -- disambiguate by R(N-O) order.
        for point, row in enumerate(
            sorted(rows, key=lambda r: (r["r_no"] is None, r["r_no"])), start=1
        ):
            stage_label = f"{stage}_{point}" if len(rows) > 1 else stage
            channels = row.pop("_channels")
            summary_rows.append({"mol": mol, "stage": stage_label, **row})
            for channel_name, (weight, mo_index, epsilon, coeff, delta_lumo, in_window) in channels.items():
                channel_rows.append({
                    "mol": mol, "stage": stage_label, "channel": channel_name,
                    "R_NO": row["r_no"], "MO_index": mo_index,
                    "epsilon_i_star": epsilon, "coefficient": coeff, "weight": weight,
                    "delta_lumo": round(delta_lumo, 5) if delta_lumo is not None else None,
                    "in_window": in_window,
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
        mol = resolve_mol_name(mol_id, dft_opt_dir)
        if mol is None:
            print(f"-- mol_{mol_id.zfill(3)}: no directory, skipping")
            continue
        mol_dir = dft_opt_dir / mol
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
