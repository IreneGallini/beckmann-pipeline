"""
Diabatic character-exchange tracking -- the production method for the
CN/C-C acceptor-character question, replacing an earlier continuity-tracked
"avoided crossing" approach that lived in beckmann.dft.branch_tracking
(deleted; see below for why). Per the PI's direct framing: acceptor
character starts N-side/C=N-dominant, becomes mixed, and ends C-C-routed
(migrating-bond-dominant) as R(N-O) increases. This is a diabatic,
per-point-independent quantity -- no branch continuity, no energy filtering.
If experimentally rearranging substrates consistently show this exchange
while fragmenting substrates do not, character exchange itself is the
predictive descriptor.

Why the old branch_tracking.py approach was retired: it tracked MO identity
across geometries via a signed-dot-product similarity metric, restricted to
candidates with positive canonical energy. The PI rejected that restriction
directly -- negative CMO eigenvalues are normal for low-lying virtuals in a
cationic system like this project's protonated oximes and don't imply
occupied/bonding character. Concrete evidence (from the reference case,
example_scans/5_s1_Me.log..5_s4_Me.log): the TRUE max-w_CN acceptor MO is
negative-energy at s2/s3/s4 (MO32 at -0.00951/-0.02490/-0.03881 a.u.) -- the
old filter would have discarded the correct acceptor state at 3 of 4 scan
points. That evidence, and the character-exchange pattern itself (max-w_CC
MO's own f_CN running 0.66 -> 0.66 -> 0.66 -> 0.00 across the same 4 points),
is plotted and tabulated in data/output/analysis/character_exchange_reference.csv
and data/output/analysis/plots/character_exchange_reference.png (produced by
scripts/analysis/plot_character_exchange.py) -- those artifacts are the
durable record now; branch_tracking.py's continuity-tracking code itself
added no further value once this module's simpler per-point approach was
validated against the same reference case, so it was deleted rather than
kept as unused legacy code.

extract_family_weights()/_is_family_member() (below) are the one part of the
old module that WAS still needed -- the per-virtual-MO family-weight
math (w_CC/w_CN/f_CC/f_CN), already Gate-1-validated against
Detailed_Orbital_Character_Exchange_Handout.docx Section 6 -- moved here
verbatim rather than re-derived.

Atom numbering: c1_atom/cn_c_atom/cn_n_atom/ref_atom default to the reference
molecule's own numbering (1/7/17/8). To run on a real benchmark substrate
instead, pass mol/mol_dir and the real per-molecule atom numbers are resolved
via beckmann_nbo.descriptors.get_substituent_map(), which returns
{ci, ni, oi, c_aryl, c_alkyl}. The mapping from that dict onto this module's
parameter names follows beckmann_nbo.parse_cmo's own channel definitions
(w17max = BD*(C{ci}-C{c_aryl}), w78max = BD*(C{ci}-C{c_alkyl})):

    c1_atom    (reference molecule's C1)  <-> c_aryl
    cn_c_atom  (reference molecule's C7)  <-> ci
    cn_n_atom  (reference molecule's N17) <-> ni
    ref_atom   (reference molecule's C8)  <-> c_alkyl

The reference case itself (example_scans/5_s*_Me.log) has no corresponding
.gjf oxime label or best_per_substrate.sdf entry, so it cannot be resolved
through get_substituent_map()/oxime_atom_map_from_gjf() -- mol/mol_dir must
stay unset and the hardcoded defaults used for that case.

Real benchmark molecules' scans live differently than the reference case:
one _scan.log per molecule with multiple CMO sections chained via --Link1--,
rather than one log per point. extract_family_weights_series()/
track_diabatic_character_series() (below) handle that by reusing
parse_cmo.collect_molecule_vir_mos() -- which already splits multi-section
logs, computes each section's actual R(N-O) (not the nominal scan step; see
JOB_ISSUES.md for a documented case where those diverged), dedups Stable=Opt
seed/final pairs, and redirects mol_003_E/mol_020_E/mol_034_E to their
working STEP_SCAN_SOURCES reruns -- rather than reimplementing any of that.
"""
import json
import re
from pathlib import Path

from beckmann_nbo.config import DATA_INPUT, DATA_OUTPUT
from beckmann_nbo.parse_cmo import (
    collect_molecule_vir_mos, find_cmo_sections, is_bond_between,
    parse_cmo_table, virtual_window,
)

HARTREE_TO_EV = 27.2114

_ATOM_MAP_KEYS = {
    "c1_atom": "c_aryl",
    "cn_c_atom": "ci",
    "cn_n_atom": "ni",
    "ref_atom": "c_alkyl",
}


def _is_family_member(label: str, element: str, atom_num: int) -> bool:
    """LP(n) {element}{atom_num} (not LP*), RY*(n) {element}{atom_num}, or BD*(n)
    involving {element}{atom_num} on either side -- the symmetric family-matching
    rule shared by both C-C-side (element='C', atom_num=c1_atom) and the
    broadened C=N-side (element='N', atom_num=cn_n_atom) per
    Detailed_Orbital_Character_Exchange_Handout.docx. Atom-mention regex is
    unanchored to position, so it matches both 'C1-X' and 'X-C1' orderings."""
    if not re.search(rf"{element}\s*{atom_num}(?!\d)", label):
        return False
    return bool(re.match(r"LP\s*\(", label) or re.match(r"RY\*\(", label) or re.match(r"BD\*\(", label))


def _family_weights_from_vir(
    vir: list[dict], c1_atom: int, cn_c_atom: int, cn_n_atom: int, ref_atom: int,
) -> list[dict]:
    """Per-MO family weights/fractions for one already virtual_window()-filtered
    MO list -- shared by extract_family_weights() (single-log/reference-case
    path) and extract_family_weights_series() (real-molecule multi-point path)
    so the coefficient math exists in exactly one place."""
    rows = []
    for mo in vir:
        w_cc = sum(c ** 2 for c, label in mo["contribs"] if _is_family_member(label, "C", c1_atom))
        w_cn = sum(c ** 2 for c, label in mo["contribs"] if _is_family_member(label, "N", cn_n_atom))
        w_ref = sum(
            c ** 2 for c, label in mo["contribs"]
            if re.match(r"BD\*\(", label) and is_bond_between(label, cn_c_atom, ref_atom)
        )
        w_target = w_cc + w_cn
        rows.append({
            "mo": mo["mo"], "energy": mo["energy"],
            "w_CC": w_cc, "w_CN": w_cn, "w_target": w_target,
            "f_CC": w_cc / w_target if w_target > 0 else None,
            "f_CN": w_cn / w_target if w_target > 0 else None,
            "w_ref": w_ref,
            "contribs": mo["contribs"],
        })
    return rows


def extract_family_weights(
    log_path: Path, c1_atom: int = 1, cn_c_atom: int = 7, cn_n_atom: int = 17,
    ref_atom: int = 8,
) -> list[dict]:
    """Per virtual MO in the 0.4 a.u.-above-LUMO window (virtual_window()), the
    C-C-side/C=N-side family weights and fractions, the separate w_ref
    (C{cn_c}-C{ref_atom} antibond) quantity, plus the MO's full raw contribs
    list. Uses the LAST CMO section in the log if there's more than one
    (Stable=Opt seed vs. final pass) -- correct for the reference case's
    one-log-per-scan-point layout (example_scans/5_s1_Me.log..5_s4_Me.log);
    for a real benchmark molecule's single multi-point _scan.log, use
    extract_family_weights_series() instead.

    Returns a list of {mo, energy, w_CC, w_CN, w_target, f_CC, f_CN, w_ref,
    contribs} dicts, one per virtual MO in the window, in the log's own
    (ascending energy) order. f_CC/f_CN are None when w_target == 0."""
    lines = Path(log_path).read_text().splitlines()
    starts = find_cmo_sections(lines)
    if not starts:
        return []
    table = parse_cmo_table(lines, starts[-1])
    vir = virtual_window(table)
    return _family_weights_from_vir(vir, c1_atom, cn_c_atom, cn_n_atom, ref_atom)


def extract_family_weights_series(
    mol: str, mol_dir: Path,
    c1_atom: int = 1, cn_c_atom: int = 7, cn_n_atom: int = 17, ref_atom: int = 8,
) -> list[dict]:
    """Per-scan-point family-weight tables for a real benchmark molecule's
    full scan (nbo/R0 stage + every scan_N stage), reusing
    parse_cmo.collect_molecule_vir_mos() for the actual multi-CMO-section
    splitting, Stable=Opt dedup, actual-R(N-O) computation, and
    STEP_SCAN_SOURCES redirection (mol_003_E/mol_020_E/mol_034_E) -- none of
    that is reimplemented here. Each point's raw vir_mos (all virtuals,
    unfiltered) is virtual_window()-filtered then run through the same
    per-MO weight math extract_family_weights() uses.

    Returns [{"r_no": float, "rows": list[dict]}, ...] sorted ascending by
    R(N-O) -- rows have the same shape extract_family_weights() returns."""
    points = collect_molecule_vir_mos(mol, mol_dir)
    return [
        {
            "r_no": point["r_no"],
            "rows": _family_weights_from_vir(
                virtual_window(point["vir_mos"]), c1_atom, cn_c_atom, cn_n_atom, ref_atom,
            ),
        }
        for point in points
        if point["r_no"] is not None
    ]


def _resolve_atom_numbers(
    mol: str | None, mol_dir: Path | None,
    c1_atom: int, cn_c_atom: int, cn_n_atom: int, ref_atom: int,
) -> tuple[int, int, int, int]:
    """Explicit params by default; if mol/mol_dir are both given, override with
    the real per-molecule atom numbers from get_substituent_map() instead --
    see module docstring for the c1/cn_c/cn_n/ref -> c_aryl/ci/ni/c_alkyl
    mapping. Imported lazily to avoid a hard dependency on best_per_substrate.sdf
    existing for callers that only ever use the explicit-atom-number path
    (e.g. the reference case, which has no sdf entry at all). Resolved once
    per molecule (not per scan point) -- get_substituent_map() depends only
    on mol/mol_dir, never on geometry, so there's nothing that could vary
    across a molecule's own scan points."""
    if mol is None and mol_dir is None:
        return c1_atom, cn_c_atom, cn_n_atom, ref_atom
    if mol is None or mol_dir is None:
        raise ValueError("mol and mol_dir must both be given, or both omitted")

    from beckmann_nbo.descriptors import get_substituent_map
    subst = get_substituent_map(mol, mol_dir)
    return subst["c_aryl"], subst["ci"], subst["ni"], subst["c_alkyl"]


def _diabatic_row(rows: list[dict]) -> dict | None:
    """Independently (no continuity, no energy filtering) pick the MO with
    max w_CN and the MO with max w_CC out of one scan point's family-weight
    rows. Shared by track_diabatic_character() (reference-case, per-log) and
    track_diabatic_character_series() (real molecule, per-point) so the
    selection rule exists in exactly one place. None if there are no rows
    (e.g. an empty virtual window at that point).

    Alongside each MO's f_CC/f_CN fraction, also carries through its raw
    (unnormalized) w_CC/w_CN weights (w_CC_CC/w_CN_CC for the max-w_CC MO,
    w_CC_CN/w_CN_CN for the max-w_CN MO) -- the fraction alone can't
    distinguish "strong absolute mixing" from "tiny w_target with an
    equally tiny w_CN dominating it," so both are exposed."""
    if not rows:
        return None
    mo_cn = max(rows, key=lambda m: m["w_CN"])
    mo_cc = max(rows, key=lambda m: m["w_CC"])
    return {
        "mo_CN": mo_cn["mo"], "E_CN": mo_cn["energy"], "f_CN_CN": mo_cn["f_CN"],
        "w_CC_CN": mo_cn["w_CC"], "w_CN_CN": mo_cn["w_CN"],
        "mo_CC": mo_cc["mo"], "E_CC": mo_cc["energy"],
        "f_CC_CC": mo_cc["f_CC"], "f_CN_CC": mo_cc["f_CN"],
        "w_CC_CC": mo_cc["w_CC"], "w_CN_CC": mo_cc["w_CN"],
    }


def track_diabatic_character(
    scan_logs: list[Path],
    mol: str | None = None,
    mol_dir: Path | None = None,
    c1_atom: int = 1, cn_c_atom: int = 7, cn_n_atom: int = 17, ref_atom: int = 8,
) -> list[dict]:
    """Per scan-point log, independently (no continuity across points, no
    energy filtering -- the full virtual_window() manifold as extract_family_
    weights() returns it):

      mo_CN = the MO with max w_CN (the N-side/C=N-family acceptor state)
      mo_CC = the MO with max w_CC (the C-C-side acceptor state)

    Returns one row per scan point: {mo_CN, E_CN, f_CN_CN, mo_CC, E_CC,
    f_CC_CC, f_CN_CC}. Energies are in Hartree (convert with this module's
    HARTREE_TO_EV for plotting). f_CN_CC -- the max-w_CC MO's OWN f_CN -- is
    exposed because it's the quantity the character-exchange framing tracks
    across scan points (N-side/mixed -> C-C-routed); this function does not
    itself classify the exchange, it just reports the numbers needed to do
    so. For a real benchmark molecule's single multi-point _scan.log, use
    track_diabatic_character_series() instead.
    """
    c1_atom, cn_c_atom, cn_n_atom, ref_atom = _resolve_atom_numbers(
        mol, mol_dir, c1_atom, cn_c_atom, cn_n_atom, ref_atom,
    )

    results = []
    for log_path in scan_logs:
        weights = extract_family_weights(log_path, c1_atom, cn_c_atom, cn_n_atom, ref_atom)
        row = _diabatic_row(weights)
        if row is not None:
            results.append(row)
    return results


def track_diabatic_character_series(
    mol: str, mol_dir: Path,
    c1_atom: int = 1, cn_c_atom: int = 7, cn_n_atom: int = 17, ref_atom: int = 8,
) -> list[dict]:
    """Whole-molecule diabatic character-exchange trace for a real benchmark
    substrate, using extract_family_weights_series() (handles the
    multi-CMO-section _scan.log and STEP_SCAN_SOURCES-merged molecules
    automatically) with atom numbers resolved ONCE via
    get_substituent_map(mol, mol_dir) -- see _resolve_atom_numbers().

    Returns one row per scan point, sorted ascending by R(N-O): {r_no, mo_CN,
    E_CN, f_CN_CN, mo_CC, E_CC, f_CC_CC, f_CN_CC} -- same per-point shape as
    track_diabatic_character(), plus r_no since these are real (non-uniform)
    geometries that need their own R(N-O) tag."""
    c1_atom, cn_c_atom, cn_n_atom, ref_atom = _resolve_atom_numbers(
        mol, mol_dir, c1_atom, cn_c_atom, cn_n_atom, ref_atom,
    )

    results = []
    for point in extract_family_weights_series(mol, mol_dir, c1_atom, cn_c_atom, cn_n_atom, ref_atom):
        row = _diabatic_row(point["rows"])
        if row is not None:
            row["r_no"] = point["r_no"]
            results.append(row)
    return results


# Side-experiment directories under data/output/dft_opt/ that are not part of
# the 34 canonical benchmark substrates (basis-set/step-size sensitivity
# tests, see CLAUDE.md) -- excluded from the benchmark-wide run below.
_NON_CANONICAL_DIRS = {
    "mol_002_E_rigidscan", "mol_006_E_631g", "mol_006_E_finescan",
    "mol_006_E_rigidscan", "mol_021_E_rigidscan",
}
_MOL_DIR_RE = re.compile(r"mol_\d{3}_[EZ]$")


def _canonical_mol_dirs(dft_opt_dir: Path) -> list[str]:
    """The 34 canonical benchmark molecule directory names under dft_opt/,
    e.g. 'mol_002_E' -- excludes the known side-experiment directories and
    anything not matching the mol_XXX_[EZ] naming pattern."""
    return sorted(
        d.name for d in dft_opt_dir.iterdir()
        if d.is_dir() and d.name not in _NON_CANONICAL_DIRS and _MOL_DIR_RE.fullmatch(d.name)
    )


def _exp_outcome(mol: str, meta: dict) -> str:
    """Substrate-level experimental R/F label for an isomer-specific molecule
    name (e.g. 'mol_002_E' -> meta['mol_002']['exp_outcome']) -- the same
    inline mol.split('_')[1] idiom every other script in this repo already
    uses for this join; benchmark_meta.json has no isomer suffix, only
    dft_opt/'s directory names do."""
    mol_id = mol.split("_")[1]
    return meta[f"mol_{mol_id}"]["exp_outcome"]


def build_benchmark_character_exchange(
    dft_opt_dir: Path | None = None,
) -> tuple[list[dict], list[dict], list[dict]]:
    """Run track_diabatic_character_series() over every canonical benchmark
    molecule, joining each against its substrate-level experimental R/F
    label. No classification threshold is computed here -- this only builds
    the raw per-point and per-molecule summary data the comparison plots
    need.

    Returns (detail_rows, summary_rows, failures):
      detail_rows  -- one row per (molecule, scan point), every metric
                       track_diabatic_character_series() computes at that
                       point: {mol, exp_outcome, point, r_no, delta_R,
                       mo_CC, E_CC, w_CC_CC, w_CN_CC, f_CC_CC, f_CN_CC,
                       mo_CN, E_CN, w_CC_CN, w_CN_CN, f_CN_CN}. r_no is the
                       point's actual R(N-O) in Angstroms; delta_R is
                       relative to that molecule's OWN first scan point
                       (typically its nbo/R0 stage, but see note below), not
                       an absolute R(N-O). w_CC_CC/w_CN_CC are the raw
                       (unnormalized) family weights of the max-w_CC MO --
                       f_CN_CC is just their ratio, and can't on its own
                       distinguish strong absolute mixing from a tiny
                       w_target dominated by an equally tiny w_CN.
      summary_rows -- one row per molecule: {mol, exp_outcome, f_CN_CC_start,
                       f_CN_CC_end, delta}, delta = f_CN_CC_start - f_CN_CC_end.
      failures     -- one row per molecule that didn't produce usable data:
                       {mol, category, detail}, category one of
                       'missing_scan_log' / 'atom_mapping' / 'no_usable_points'
                       / 'other'.

    Note on "first scan point": it is whichever point sorts first by actual
    R(N-O) in track_diabatic_character_series()'s output -- normally the nbo/
    R0 equilibrium stage, since scan points only ever stretch R(N-O) upward
    from there. If a molecule's nbo.log is missing/failed but its scan.log
    succeeded, the first scan point stands in as R0 instead -- delta_R stays
    internally consistent for that molecule (always relative to its own
    first available point), it just may not be the literal equilibrium
    geometry in that case.
    """
    dft_opt_dir = dft_opt_dir or (DATA_OUTPUT / "dft_opt")
    meta = json.loads((DATA_INPUT / "benchmark_meta.json").read_text())

    detail_rows: list[dict] = []
    summary_rows: list[dict] = []
    failures: list[dict] = []

    for mol in _canonical_mol_dirs(dft_opt_dir):
        mol_dir = dft_opt_dir / mol
        try:
            outcome = _exp_outcome(mol, meta)
        except KeyError as e:
            failures.append({"mol": mol, "category": "missing_scan_log", "detail": f"no benchmark_meta.json entry: {e}"})
            continue

        try:
            rows = track_diabatic_character_series(mol, mol_dir)
        except ValueError as e:
            failures.append({"mol": mol, "category": "atom_mapping", "detail": str(e)})
            continue
        except FileNotFoundError as e:
            failures.append({"mol": mol, "category": "missing_scan_log", "detail": str(e)})
            continue
        except Exception as e:
            failures.append({"mol": mol, "category": "other", "detail": f"{type(e).__name__}: {e}"})
            continue

        if not rows:
            failures.append({"mol": mol, "category": "no_usable_points", "detail": "track_diabatic_character_series() returned no rows"})
            continue

        r0 = rows[0]["r_no"]
        for point, row in enumerate(rows, start=1):
            detail_rows.append({
                "mol": mol, "exp_outcome": outcome, "point": point,
                "r_no": row["r_no"], "delta_R": row["r_no"] - r0,
                "mo_CC": row["mo_CC"], "E_CC": row["E_CC"],
                "w_CC_CC": row["w_CC_CC"], "w_CN_CC": row["w_CN_CC"],
                "f_CC_CC": row["f_CC_CC"], "f_CN_CC": row["f_CN_CC"],
                "mo_CN": row["mo_CN"], "E_CN": row["E_CN"],
                "w_CC_CN": row["w_CC_CN"], "w_CN_CN": row["w_CN_CN"],
                "f_CN_CN": row["f_CN_CN"],
            })
        summary_rows.append({
            "mol": mol, "exp_outcome": outcome,
            "f_CN_CC_start": rows[0]["f_CN_CC"], "f_CN_CC_end": rows[-1]["f_CN_CC"],
            "delta": rows[0]["f_CN_CC"] - rows[-1]["f_CN_CC"],
        })

    return detail_rows, summary_rows, failures
