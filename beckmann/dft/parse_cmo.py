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
        that intermediate detail is kept here rather than discarded.

classify_crossing() (below) separately confirms whether a 'cn' channel MO handoff is a
real avoided crossing, by following the SPECIFIC pre/post-handoff MO pair across the
whole series (find_handoff_pair()/track_mo_pair(), against the raw log data via
collect_molecule_vir_mos()) rather than independently re-picking a 'runner-up' at each
geometry. See data/output/analysis/cn_crossing_report.csv /
data/output/analysis/cn_handoff_ledger/ (scripts/analysis/cn_crossing_report.py,
scripts/dft/inspect_cn_ledger.py) -- across all 34 benchmark molecules, the crossing
partner is consistently the N-O sigma*/sigma antibond (the breaking bond itself), never
the aryl-migrating C-C antibond, so classify_crossing() only checks the eigenvalue
(bracketed small gap) and CN-weight-conservation signatures, not any aryl-coefficient
swap.
"""
import csv
import json
import math
import re
from pathlib import Path

from beckmann.config import DATA_INPUT, DATA_OUTPUT
from beckmann.dft.descriptors import get_substituent_map, resolve_series
from beckmann.dft.inputs import (
    ALL_IDS, STEP_SCAN_SOURCES, build_stage_relabel_map, relabel_rows,
    resolve_mol_name, step_scan_dir,
)
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


def label_involves(label: str, atom_num: int) -> bool:
    """True if label mentions this atom number, for any element (any single/double
    letter symbol followed by the number, not followed by another digit)."""
    return re.search(rf"[A-Z][a-z]?\s*{atom_num}(?!\d)", label) is not None


def is_bond_between(label: str, a: int, b: int) -> bool:
    """Generalized is_target_antibond(): True if label mentions BOTH atom numbers a
    and b, any element -- unlike is_target_antibond(), doesn't assume atom a is
    carbon, so this also matches e.g. the N-O antibond (BD*(O{oi}-N{ni})), needed to
    identify the CN-channel handoff's crossing partner (see
    beckmann.dft.parse_cmo.build_cn_crossing_report(), scripts/analysis/
    cn_crossing_report.py)."""
    return label_involves(label, a) and label_involves(label, b)


def all_weight_matches_for_target(vir_mos: list[dict], a: int, b: int) -> list[dict]:
    """Every virtual MO containing the BD* antibond between atoms a and b, as
    {weight, mo_index, epsilon_i_star, coefficient} dicts, sorted descending by
    weight (squared coefficient) -- exposes the runner-up match(es), not just
    the single winner max_weight_for_target() reduces this to. A given MO
    contributes at most one entry (is_target_antibond matches whichever of
    that MO's own contribs is the target antibond). Empty list if the
    antibond never appears in vir_mos.
    """
    matches = [
        {
            "weight": coeff ** 2, "mo_index": mo["mo"],
            "epsilon_i_star": mo["energy"], "coefficient": coeff,
        }
        for mo in vir_mos
        for coeff, label in mo["contribs"]
        if is_target_antibond(label, a, b)
    ]
    matches.sort(key=lambda m: m["weight"], reverse=True)
    return matches


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

    Just matches[0] from all_weight_matches_for_target() -- kept as its own function with a
    byte-identical return shape/order since wX^max/Lambda/Psi and several scripts
    (validate_reference_descriptors.py, compare_wcnmax_window.py) unpack this positionally.
    """
    matches = all_weight_matches_for_target(vir_mos, a, b)
    if not matches:
        return None, None, None, None, None, None
    best = matches[0]
    best_val, best_mo = best["weight"], best["mo_index"]
    best_epsilon, best_coeff = best["epsilon_i_star"], best["coefficient"]
    if lumo_e is None:
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


def coefficient_in_mo(vir_mos: list[dict], a: int, b: int, mo_index: int | None) -> float | None:
    """Signed coefficient of the BD*(a,b) antibond in one SPECIFIC MO (looked up by
    index), or None if that MO doesn't carry it (or mo_index is None). Unlike
    max_weight_for_target()/all_weight_matches_for_target(), this doesn't search for
    the best match across the manifold -- it looks up one already-known MO, e.g. to
    ask "what's the aryl C-C antibond's own character in the MO that wins the CN
    channel" rather than wherever the aryl channel's own max happens to sit."""
    if mo_index is None:
        return None
    for mo in vir_mos:
        if mo["mo"] != mo_index:
            continue
        for coeff, label in mo["contribs"]:
            if is_target_antibond(label, a, b):
                return coeff
    return None


def mo_neighborhood(vir_mos: list[dict], target_mo: int, n: int = 3) -> list[dict]:
    """target_mo plus its n nearest-energy-ordered neighbor virtual MOs (vir_mos is
    already ascending-energy/ascending-index, straight from parse_cmo_table), each
    carrying its FULL leading-NBO-contribution list (mo['contribs']) -- not just one
    target antibond's coefficient. This is the human-inspectable ledger view: "what is
    this MO and its neighbors actually made of" (MO number + bond label), used to
    sanity-check a naive per-point-independent runner-up pick against what's really
    sitting next to the winner in energy. Returns [] if target_mo isn't in vir_mos.
    Each entry also carries is_target so the caller can highlight it."""
    idx = next((i for i, m in enumerate(vir_mos) if m["mo"] == target_mo), None)
    if idx is None:
        return []
    lo, hi = max(0, idx - n), min(len(vir_mos), idx + n + 1)
    return [
        {"mo": m["mo"], "energy": m["energy"], "contribs": m["contribs"], "is_target": m["mo"] == target_mo}
        for m in vir_mos[lo:hi]
    ]


def track_mo_pair(points: list[dict], mo_a: int, mo_b: int, ci: int, ni: int) -> list[dict]:
    """Follow two SPECIFIC, identity-fixed virtual MOs (typically a handoff's
    pre-winner and post-winner) across every geometry in `points`
    (collect_molecule_vir_mos()'s output) -- unlike naively re-picking, at each
    geometry independently, whichever OTHER MO has the 2nd-largest CN weight right
    there, which can latch onto a persistent, structurally unrelated third MO
    instead of the real handoff partner (see classify_crossing()'s docstring /
    Notes.md's mol_014_Z vs. mol_001_E cases, and
    data/output/analysis/cn_crossing_report.csv, which used this to confirm the
    real partner across all 34 benchmark molecules is consistently the N-O
    sigma*/sigma antibond, not the aryl-migrating C-C antibond -- so this no longer
    tracks aryl character at all).

    Returns one row per point: r_no, each MO's energy/CN-coefficient (None if that
    MO isn't even printed at that geometry -- below NBO7's 5% print threshold), and
    their energy gap (None if either MO is missing at that point). Reuses
    coefficient_in_mo() -- no new lookup logic."""
    rows = []
    for point in points:
        vir_mos = point["vir_mos"]
        energy = {m["mo"]: m["energy"] for m in vir_mos}
        e_a, e_b = energy.get(mo_a), energy.get(mo_b)
        rows.append({
            "r_no": point["r_no"],
            "mo_a": mo_a, "epsilon_a": e_a,
            "cn_coeff_a": coefficient_in_mo(vir_mos, ci, ni, mo_a) if e_a is not None else None,
            "mo_b": mo_b, "epsilon_b": e_b,
            "cn_coeff_b": coefficient_in_mo(vir_mos, ci, ni, mo_b) if e_b is not None else None,
            "gap": round(abs(e_a - e_b), 5) if e_a is not None and e_b is not None else None,
        })
    return rows


def parse_log_ledger(log_path: Path, ni: int, oi: int) -> dict:
    """Like parse_log() but keeps each geometry's raw virtual-MO table (mo/energy/
    contribs) instead of reducing it to compute_descriptors()'s per-channel summary
    -- for ledger/pair-tracking inspection (mo_neighborhood(), track_mo_pair()), not
    the main cmo_descriptors.csv/cmo_channel_extraction.csv pipeline. Same
    last-table-per-R dedup as parse_log() (Stable=Opt's seed-vs-final pass)."""
    lines = log_path.read_text().splitlines()
    vir_by_r: dict = {}
    for start in find_cmo_sections(lines):
        table = parse_cmo_table(lines, start)
        if not table:
            continue
        r_no = r_no_before(lines, start, ni, oi)
        r_key = round(r_no, 4) if r_no is not None else None
        vir_by_r[r_key] = [m for m in table if m["kind"] == "vir"]
    return vir_by_r


def collect_molecule_vir_mos(mol: str, mol_dir: Path) -> list[dict]:
    """Raw per-geometry virtual-MO data across all of a molecule's stage logs (nbo +
    scan, STEP_SCAN_SOURCES-merged the same way collect_molecule_stepscan() is for
    the main pipeline), sorted by R(N-O) ascending -- [{'r_no', 'vir_mos'}, ...].
    Bypasses compute_descriptors()'s reduction to summary/channel rows entirely; feeds
    mo_neighborhood()/track_mo_pair() for ledger inspection. Skips any log that
    doesn't exist or didn't reach Normal termination (silently, unlike the main
    pipeline's collect_molecule() -- this is a best-effort inspection tool, not a
    trusted-completeness gate)."""
    points: dict = {}
    if mol in STEP_SCAN_SOURCES:
        _, ni, oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol}_opt.gjf")
        nbo_log = mol_dir / f"{mol}_nbo.log"
        if nbo_log.exists() and log_terminated_normally(nbo_log):
            points.update(parse_log_ledger(nbo_log, ni, oi))
        for source in STEP_SCAN_SOURCES[mol]:
            source_dir = step_scan_dir() / source
            _, s_ni, s_oi, _ = oxime_atom_map_from_gjf(source_dir / f"{source}_opt.gjf")
            scan_log = source_dir / f"{source}_scan.log"
            if scan_log.exists() and log_terminated_normally(scan_log):
                points.update(parse_log_ledger(scan_log, s_ni, s_oi))
    else:
        _, ni, oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol}_opt.gjf")
        for stage in STAGES:
            log_path = mol_dir / f"{mol}_{stage}.log"
            if log_path.exists() and log_terminated_normally(log_path):
                points.update(parse_log_ledger(log_path, ni, oi))
    return [
        {"r_no": r, "vir_mos": vm}
        for r, vm in sorted(points.items(), key=lambda kv: (kv[0] is None, kv[0]))
    ]


def find_handoff_pair(mol: str, extraction_rows: list[dict]) -> dict:
    """The CN channel's winning-MO series (via resolve_series()) plus -- if the
    winning MO_index ever changes along the scan -- the specific pre-handoff and
    post-handoff MO indices and the R at which the switch happens. Shared by
    classify_crossing(), write_cn_ledger(), and build_cn_crossing_report(), which all
    need the exact same two MO indices to hand to track_mo_pair() -- the
    identity-tracked pair, as opposed to independently re-picking, at each geometry
    separately, whichever OTHER MO has the 2nd-largest CN weight right there (which
    can miss the real handoff partner -- see classify_crossing()'s docstring).

    Always returns {"pts", "handoff_idx", "pre_mo", "post_mo", "handoff_R"} -- "pts"
    is every point's {R, mo, weight}, sorted by R. "handoff_idx" (and the rest) are
    None when there are fewer than 2 usable points or the winning MO_index never
    changes along the series at all.
    """
    by_stage = {
        r["stage"]: r for r in extraction_rows
        if r["mol"] == mol and r["channel"] == "cn" and r["weight"] not in (None, "", "None")
    }
    rows = resolve_series(by_stage)
    pts = sorted(
        (
            {
                "R": float(r["R_NO"]),
                "mo": int(float(r["MO_index"])) if r.get("MO_index") not in (None, "", "None") else None,
                "weight": float(r["weight"]) if r.get("weight") not in (None, "", "None") else None,
            }
            for r in rows
        ),
        key=lambda p: p["R"],
    )
    handoff_idx = (
        next((i for i in range(1, len(pts)) if pts[i]["mo"] != pts[i - 1]["mo"]), None)
        if len(pts) >= 2 else None
    )
    if handoff_idx is None:
        return {"pts": pts, "handoff_idx": None, "pre_mo": None, "post_mo": None, "handoff_R": None}
    pre, post = pts[handoff_idx - 1], pts[handoff_idx]
    return {
        "pts": pts, "handoff_idx": handoff_idx,
        "pre_mo": pre["mo"], "post_mo": post["mo"], "handoff_R": post["R"],
    }


def classify_crossing(mol: str, extraction_rows: list[dict], mol_dir: Path, mo_gap_threshold: float = 0.03) -> dict:
    """Classify the CN channel's wCNmax MO handoff along one molecule's R(N-O)
    series as a CONFIRMED avoided crossing, an unconfirmed handoff, or no handoff at
    all -- a more specific signature than find_wcnmax_extremum()'s plain shape test.

    Unlike an earlier version of this function, which independently re-picked a
    'runner-up' MO at each geometry separately (whichever OTHER MO had the
    2nd-largest CN weight right there), this follows the SPECIFIC pre/post-handoff
    MO pair (find_handoff_pair()) across the WHOLE series via track_mo_pair(),
    against the raw log data (collect_molecule_vir_mos()) -- the earlier per-point
    approach could latch onto a persistent, structurally unrelated third MO instead
    of the real handoff partner (see Notes.md's mol_014_Z vs. mol_001_E write-up,
    data/output/analysis/cn_handoff_ledger/).

    "confirmed avoided crossing" requires BOTH:
      1. Eigenvalue signature -- the identity-tracked pair's gap(R) reaches a
         bracketed local minimum immediately at the handoff (not merely small
         somewhere unrelated in the series), that minimum is > 0 (a true
         intersection would mean 0, not avoided) and < mo_gap_threshold (default
         0.03 a.u., Notes.md's empirical calibration).
      2. Coefficient signature -- the CN-channel weight is roughly conserved across
         the pair at the bracket (their weights sum close to the single-orbital
         value away from the crossing -- character splits rather than vanishing).

    Does NOT check for any aryl-antibond coefficient swap: build_cn_crossing_report()
    established, across all 34 benchmark molecules, that the CN channel's crossing
    partner is consistently the N-O sigma*/sigma antibond (the breaking bond
    itself), never the aryl-migrating C-C antibond -- see
    data/output/analysis/cn_crossing_report.csv. That hypothesis (whether the
    migrating C-C antibond itself picks up CN character) would need to track the
    aryl channel's OWN winning MO (w17max_mo) for its own crossing, not check for
    aryl character inside the CN channel's pair -- not implemented here.

    Returns a dict: {"label": "confirmed avoided crossing" | "unconfirmed handoff" |
    "no handoff", "handoff_R", "gap_min", "R_gap_min", "reason"}.
    """
    no_result = {"handoff_R": None, "gap_min": None, "R_gap_min": None}
    handoff = find_handoff_pair(mol, extraction_rows)
    pts, handoff_idx = handoff["pts"], handoff["handoff_idx"]
    if handoff_idx is None:
        reason = "fewer than 2 scan points" if len(pts) < 2 else None
        return {"label": "no handoff", "reason": reason, **no_result}

    pre_mo, post_mo, handoff_R = handoff["pre_mo"], handoff["post_mo"], handoff["handoff_R"]
    pre_R = pts[handoff_idx - 1]["R"]

    ci, ni, oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol}_opt.gjf")
    points = collect_molecule_vir_mos(mol, mol_dir)
    track_rows = track_mo_pair(points, pre_mo, post_mo, ci, ni)

    gap_rows = [(r["r_no"], r["gap"]) for r in track_rows if r["gap"] is not None]
    if not gap_rows:
        return {
            "label": "unconfirmed handoff", "handoff_R": handoff_R, "gap_min": None, "R_gap_min": None,
            "reason": "neither handoff MO co-exists anywhere in the series",
        }
    R_gap_min, gap_min = min(gap_rows, key=lambda t: t[1])
    bracketed = R_gap_min in (pre_R, handoff_R)
    eigenvalue_ok = bracketed and 0 < gap_min < mo_gap_threshold

    if not eigenvalue_ok:
        reason = (
            f"gap minimum ({gap_min} a.u. at R={R_gap_min}) isn't at the handoff bracket (R={pre_R}/{handoff_R})"
            if not bracketed else
            f"gap minimum {gap_min} a.u. not below threshold {mo_gap_threshold}"
        )
        return {
            "label": "unconfirmed handoff", "handoff_R": handoff_R,
            "gap_min": gap_min, "R_gap_min": R_gap_min, "reason": reason,
        }

    bracket_row = next(r for r in track_rows if r["r_no"] == R_gap_min)
    cn_a, cn_b = bracket_row["cn_coeff_a"], bracket_row["cn_coeff_b"]
    bracket_weight = (cn_a or 0) ** 2 + (cn_b or 0) ** 2
    # single-orbital reference: the largest CN coefficient observed for EITHER MO
    # anywhere in the series -- what the pair's combined weight looks like away from
    # the crossing, when character sits fully on one orbital rather than splitting.
    all_coeffs = [
        abs(c) for r in track_rows for c in (r["cn_coeff_a"], r["cn_coeff_b"]) if c is not None
    ]
    other_weight = max(all_coeffs) ** 2 if all_coeffs else None
    conserved = (
        other_weight is not None and bracket_weight > 0
        and abs(bracket_weight - other_weight) <= 0.2 * other_weight
    )

    label = "confirmed avoided crossing" if conserved else "unconfirmed handoff"
    reason = None if conserved else "CN weight not roughly conserved across the pair at the bracket"
    return {"label": label, "handoff_R": handoff_R, "gap_min": gap_min, "R_gap_min": R_gap_min, "reason": reason}


def compute_descriptors(
    mo_table: list[dict], ci: int, ni: int, c_aryl: int, c_alkyl: int,
) -> tuple[dict, dict]:
    """Returns (summary_dict, channels). channels is the raw per-channel detail
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


def collect_stage(mol: str, mol_dir: Path, ci: int, ni: int, oi: int, c_aryl: int, c_alkyl: int,
                   stage: str) -> tuple[list[dict], list[dict]] | None:
    """Extract one stage's rows, requiring only that stage's own log (if
    present) to have converged -- see parse_nbo.collect_stage for the same
    pattern. Returns ([], []) if the log doesn't exist, None if it exists but
    didn't reach Normal termination."""
    log_path = mol_dir / f"{mol}_{stage}.log"
    if not log_path.exists():
        return [], []
    if not log_terminated_normally(log_path):
        return None

    rows = parse_log(log_path, ci, ni, oi, c_aryl, c_alkyl)
    summary_rows = []
    channel_rows = []
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


def collect_molecule(mol: str, mol_dir: Path, c_aryl: int, c_alkyl: int) -> tuple[list[dict], list[dict]]:
    """Compute descriptors for all available stage logs of one molecule, e.g. 'mol_002_E'.

    Returns (summary_rows, channel_extraction_rows) -- the latter is the per-geometry,
    per-channel (cn/17/78) detail table: which MO carried the max weight, its orbital
    energy, and the signed coefficient before squaring.

    If any present stage log didn't reach Normal termination, the whole molecule is
    skipped rather than partially included -- see parse_nbo.collect_molecule for why.
    """
    ci, ni, oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol}_opt.gjf")

    summary_rows = []
    channel_rows = []
    for stage in STAGES:
        result = collect_stage(mol, mol_dir, ci, ni, oi, c_aryl, c_alkyl, stage)
        if result is None:
            print(f"   -- {mol}: {mol}_{stage}.log did not reach Normal termination "
                  f"-- skipping whole molecule (see JOB_ISSUES.md)")
            return [], []
        s_rows, c_rows = result
        summary_rows.extend(s_rows)
        channel_rows.extend(c_rows)
    return summary_rows, channel_rows


def collect_molecule_stepscan(mol: str, mol_dir: Path, c_aryl: int, c_alkyl: int) -> tuple[list[dict], list[dict]]:
    """Merge one or more successful step-size reruns with the canonical
    equilibrium NBO -- see parse_nbo.collect_molecule_stepscan for the same
    pattern (this mirrors it for the CMO summary + channel-extraction rows,
    which must stay consistently renumbered with each other)."""
    ci, ni, oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol}_opt.gjf")

    nbo_result = collect_stage(mol, mol_dir, ci, ni, oi, c_aryl, c_alkyl, "nbo")
    if nbo_result is None:
        print(f"   -- {mol}: {mol}_nbo.log did not reach Normal termination -- skipping")
        return [], []
    nbo_summary, nbo_channel = nbo_result

    all_summary = []
    all_channel = []
    for source in STEP_SCAN_SOURCES[mol]:
        source_dir = step_scan_dir() / source
        s_ci, s_ni, s_oi, _ = oxime_atom_map_from_gjf(source_dir / f"{source}_opt.gjf")
        result = collect_stage(source, source_dir, s_ci, s_ni, s_oi, c_aryl, c_alkyl, "scan")
        if result is None:
            print(f"   -- {mol}: {source}_scan.log did not reach Normal termination -- skipping this source")
            continue
        s_rows, c_rows = result
        all_summary.extend(s_rows)
        all_channel.extend(c_rows)

    # Relabel using the SUMMARY rows' r_no as the authoritative point set --
    # channel rows must use the exact same mapping so a given point's summary
    # and channel-detail rows stay under the same renumbered stage label.
    relabel = build_stage_relabel_map({r["r_no"] for r in all_summary})
    summary_rows = relabel_rows(all_summary, mol, relabel)
    channel_rows  = relabel_rows(all_channel, mol, relabel, r_no_key="R_NO")
    return nbo_summary + summary_rows, nbo_channel + channel_rows


def _fmt(v, nd: int = 4) -> str:
    return f"{v:.{nd}f}" if isinstance(v, (int, float)) else "--"


def _fmt_contribs(contribs: list[tuple[float, str]], top_n: int = 4) -> str:
    ranked = sorted(contribs, key=lambda t: abs(t[0]), reverse=True)[:top_n]
    return "; ".join(f"{coeff:+.3f}*{label}" for coeff, label in ranked)


def write_cn_ledger(mol_id: str) -> Path | None:
    """Write a human-readable markdown ledger for one molecule's CN-channel handoff,
    straight from the raw NBO7 CMO output -- the "MO number and bond label" table
    requested to sanity-check a naive top-2-by-weight runner-up choice before
    trusting/extending classify_crossing(). Not part of the main
    cmo_descriptors.csv/cmo_channel_extraction.csv pipeline -- a one-off inspection
    tool, called by scripts/dft/inspect_cn_ledger.py.

    Two sections:
      1. The IDENTITY-TRACKED pre/post-handoff MO pair's full gap(R) trajectory
         (track_mo_pair()) -- unlike independently re-picking a 'runner-up' MO per
         point, a real avoided crossing shows a bracketed local minimum in `gap`
         near the handoff, not a monotonic drift or a persistently-large gap to an
         unrelated third MO (see Notes.md's mol_014_Z vs. mol_001_E write-up).
      2. Each of that pair's full leading-NBO-contribution neighborhood
         (mo_neighborhood()) at the two stages bracketing the handoff -- the raw
         MO-number + bond-label composition, not just the isolated CN/aryl
         coefficients.

    Returns the written path, or None if the molecule has no directory, no
    cmo_channel_extraction.csv data yet, or find_handoff_pair() finds no MO_index
    change at all (nothing to track)."""
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    mol = resolve_mol_name(mol_id, dft_opt_dir)
    if mol is None:
        print(f"-- mol_{mol_id.zfill(3)}: no directory, skipping")
        return None
    mol_dir = dft_opt_dir / mol

    extraction_path = DATA_OUTPUT / "analysis" / "cmo_channel_extraction.csv"
    with open(extraction_path) as f:
        extraction_rows = list(csv.DictReader(f))
    handoff = find_handoff_pair(mol, extraction_rows)
    if handoff["handoff_idx"] is None:
        print(f"-- {mol}: no CN-channel MO handoff found in cmo_channel_extraction.csv, skipping")
        return None
    pre_mo, post_mo, handoff_R = handoff["pre_mo"], handoff["post_mo"], handoff["handoff_R"]
    pts, idx = handoff["pts"], handoff["handoff_idx"]
    pre_R, post_R = pts[idx - 1]["R"], pts[idx]["R"]

    ci, ni, oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol}_opt.gjf")

    points = collect_molecule_vir_mos(mol, mol_dir)
    track_rows = track_mo_pair(points, pre_mo, post_mo, ci, ni)

    lines = [
        f"# CN-channel handoff ledger: {mol}",
        "",
        f"ci=C{ci}, ni=N{ni}. Handoff per cmo_channel_extraction.csv's "
        f"winning MO_index: MO{pre_mo} -> MO{post_mo} at R={handoff_R}.",
        "",
        "## 1. Identity-tracked pair trajectory",
        "",
        f"MO{pre_mo} and MO{post_mo} followed across every scan point (not "
        "re-picked per point) -- a real avoided crossing shows `gap` reaching a "
        "bracketed local minimum near the handoff row (marked below), not a "
        "monotonic drift or a persistently large value.",
        "",
        "| R(N-O) | MO A | eps_A | CN_A | MO B | eps_B | CN_B | gap |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in track_rows:
        marker = " **<- handoff**" if r["r_no"] == handoff_R else ""
        lines.append(
            f"| {r['r_no']}{marker} | {r['mo_a']} | {_fmt(r['epsilon_a'], 5)} | {_fmt(r['cn_coeff_a'], 3)} | "
            f"{r['mo_b']} | {_fmt(r['epsilon_b'], 5)} | {_fmt(r['cn_coeff_b'], 3)} | {_fmt(r['gap'], 5)} |"
        )

    lines += ["", "## 2. MO composition at the handoff bracket", ""]
    by_r = {p["r_no"]: p["vir_mos"] for p in points}
    for label, r_no, target in [
        (f"pre-handoff (R={pre_R}, winner=MO{pre_mo})", pre_R, pre_mo),
        (f"post-handoff (R={post_R}, winner=MO{post_mo})", post_R, post_mo),
    ]:
        lines.append(f"### {label}")
        lines.append("")
        lines.append("| MO | energy | winner? | leading contributions (bond label) |")
        lines.append("|---|---|---|---|")
        for m in mo_neighborhood(by_r.get(r_no, []), target, n=3):
            marker = "**yes**" if m["is_target"] else ""
            lines.append(f"| {m['mo']} | {_fmt(m['energy'], 5)} | {marker} | {_fmt_contribs(m['contribs'])} |")
        lines.append("")

    out_dir = DATA_OUTPUT / "analysis" / "cn_handoff_ledger"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{mol}.md"
    out_path.write_text("\n".join(lines) + "\n")
    print(f"-- {mol}: ledger -> {out_path}")
    return out_path


CROSSING_REPORT_FIELDS = [
    "mol", "exp", "pre_mo", "post_mo", "handoff_R",
    "gap_at_handoff", "gap_min", "R_at_gap_min", "bracketed_at_handoff",
    "eigenvalue_crossing_confirmed",
    "partner_label_pre", "partner_label_post", "partner_is_ON_sigma_star",
]


def build_cn_crossing_report(mol_ids: list[str] | None = None, mo_gap_threshold: float = 0.03) -> list[dict]:
    """For every benchmark molecule with a CN-channel MO handoff, document what the
    identity-tracked crossing partner actually IS -- straight from the raw NBO7 CMO
    output, using the SAME pre/post-handoff MO pair find_handoff_pair()/
    track_mo_pair() use (see write_cn_ledger()), not a naive per-point-repicked
    runner-up.

    Confirmed the finding from 6 hand-inspected priority molecules generalizes across
    the whole 34-molecule benchmark set: the CN channel's handoff partner is
    consistently the N-O sigma*/sigma antibond (BD*(O{oi}-N{ni}), the breaking bond
    itself), never the aryl-migrating C-C antibond -- see
    data/output/analysis/cn_crossing_report.csv (this function's own prior output,
    kept as the record) for the full 34/34 result. classify_crossing() (below)
    therefore doesn't check for any aryl-coefficient swap.

    One row per molecule with a handoff (skips molecules with no MO_index change,
    or missing cmo_channel_extraction.csv/log data). 'partner_label_pre'/'_post' is
    the dominant (largest |coefficient|) NBO contribution of whichever MO in the
    pair ISN'T the CN winner at that bracket stage -- e.g. 'BD*( 1) O 1- N 2*' --
    checked with is_bond_between() against (oi, ni) for 'partner_is_ON_sigma_star'.

    Output: data/output/analysis/cn_crossing_report.csv
    """
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    extraction_path = DATA_OUTPUT / "analysis" / "cmo_channel_extraction.csv"
    with open(extraction_path) as f:
        extraction_rows = list(csv.DictReader(f))
    outcomes = json.loads((DATA_INPUT / "benchmark_meta.json").read_text())

    report_rows = []
    for mol_id in sorted(mol_ids or ALL_IDS):
        mol = resolve_mol_name(mol_id, dft_opt_dir)
        if mol is None:
            print(f"-- mol_{mol_id.zfill(3)}: no directory, skipping")
            continue
        mol_dir = dft_opt_dir / mol

        handoff = find_handoff_pair(mol, extraction_rows)
        if handoff["handoff_idx"] is None:
            print(f"-- {mol}: no CN-channel MO handoff, skipping")
            continue
        pre_mo, post_mo, handoff_R = handoff["pre_mo"], handoff["post_mo"], handoff["handoff_R"]
        pts, idx = handoff["pts"], handoff["handoff_idx"]
        pre_R, post_R = pts[idx - 1]["R"], pts[idx]["R"]

        ci, ni, oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol}_opt.gjf")

        points = collect_molecule_vir_mos(mol, mol_dir)
        track_rows = track_mo_pair(points, pre_mo, post_mo, ci, ni)
        gap_rows = [(r["r_no"], r["gap"]) for r in track_rows if r["gap"] is not None]
        if not gap_rows:
            print(f"-- {mol}: neither handoff MO co-exists anywhere in the series, skipping")
            continue
        R_at_gap_min, gap_min = min(gap_rows, key=lambda t: t[1])
        gap_at_handoff = next((g for r, g in gap_rows if r == handoff_R), None)
        bracketed = R_at_gap_min in (pre_R, post_R)
        confirmed = bracketed and 0 < gap_min < mo_gap_threshold

        by_r = {p["r_no"]: p["vir_mos"] for p in points}

        def leading_label(r_no, mo_index):
            for m in by_r.get(r_no, []):
                if m["mo"] == mo_index and m["contribs"]:
                    _, label = max(m["contribs"], key=lambda t: abs(t[0]))
                    return label
            return None

        partner_label_pre = leading_label(pre_R, post_mo)   # not yet winner at pre_R
        partner_label_post = leading_label(post_R, pre_mo)  # no longer winner at post_R
        partner_is_on = all(
            lbl is not None and is_bond_between(lbl, oi, ni)
            for lbl in (partner_label_pre, partner_label_post)
        )

        mol_num = mol.split("_")[1]
        exp = outcomes.get(f"mol_{mol_num}", {}).get("exp_outcome", "")

        report_rows.append({
            "mol": mol, "exp": exp, "pre_mo": pre_mo, "post_mo": post_mo, "handoff_R": handoff_R,
            "gap_at_handoff": gap_at_handoff, "gap_min": gap_min, "R_at_gap_min": R_at_gap_min,
            "bracketed_at_handoff": bracketed, "eigenvalue_crossing_confirmed": confirmed,
            "partner_label_pre": partner_label_pre, "partner_label_post": partner_label_post,
            "partner_is_ON_sigma_star": partner_is_on,
        })
        print(
            f"-- {mol} ({exp}): MO{pre_mo}->MO{post_mo} @ R={handoff_R}, "
            f"gap_min={gap_min:.5f} @ R={R_at_gap_min}, confirmed={confirmed}, "
            f"N-O sigma* partner={partner_is_on}"
        )

    out_path = DATA_OUTPUT / "analysis" / "cn_crossing_report.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CROSSING_REPORT_FIELDS)
        writer.writeheader()
        writer.writerows(report_rows)
    print(f"\n{len(report_rows)} molecules -> {out_path}")
    return report_rows


def main() -> None:
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    out_path        = DATA_OUTPUT / "analysis" / "cmo_descriptors.csv"
    extraction_path = DATA_OUTPUT / "analysis" / "cmo_channel_extraction.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    all_rows = []
    all_channel_rows = []
    for mol_id in sorted(ALL_IDS):
        mol = resolve_mol_name(mol_id, dft_opt_dir)
        if mol is None:
            print(f"-- mol_{mol_id.zfill(3)}: no directory, skipping")
            continue
        mol_dir = dft_opt_dir / mol
        subst = get_substituent_map(mol, mol_dir)
        if mol in STEP_SCAN_SOURCES:
            rows, channel_rows = collect_molecule_stepscan(mol, mol_dir, subst["c_aryl"], subst["c_alkyl"])
        else:
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
