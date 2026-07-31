"""
Cross-parser aggregation for the channel-resolved descriptors from "Ring Size
and Substituent Effects in the Beckmann Rearrangement" (Sections 2.2-2.5):
substituent role tagging (aryl vs alkyl), Psi, and least-squares d/dR slopes.
NBO-specific -- Psi/Lambda are E2PERT/CMO-derived quantities that only exist
on the Gaussian/NBO7 side. (The geometry-source-agnostic wCNmax-minimum rule
itself -- resolve_series()/find_wcnmax_minimum()/find_wcnmax_extremum() --
moved to beckmann_core.wcnmax_rule; import it from there, not from here.)

This sits above parse_nbo.py (E2PERT) and parse_cmo.py (CMO/Lambda/wX^max)
because Psi and the d/dR slopes need data from both, joined with a channel
role (aryl/alkyl) that neither of those parsers determines on its own.

Output: data/output/analysis/channel_descriptors.csv (per mol/stage/r_no)
        data/output/analysis/descriptor_slopes.csv   (one row per molecule)
"""
import csv
import re
from pathlib import Path

from rdkit import Chem

from beckmann_core.classical import get_oxime_atoms
from beckmann_core.wcnmax_rule import resolve_series
from beckmann_nbo.config import DATA_OUTPUT
from beckmann_nbo.inputs import ALL_IDS, resolve_mol_name
from beckmann_nbo.scan import oxime_atom_map_from_gjf

BEST_PER_SUBSTRATE_SDF = DATA_OUTPUT / "aimnet_optimized" / "best_per_substrate.sdf"

PSI_EPSILON = 1e-6  # unspecified in the paper/handouts only matters when K_frag ~ 0

_MOL_CACHE: dict[str, Chem.Mol] = {}


def _load_mols() -> dict[str, Chem.Mol]:
    """Lazily load and cache best_per_substrate.sdf, keyed by molecule name."""
    if not _MOL_CACHE:
        suppl = Chem.SDMolSupplier(str(BEST_PER_SUBSTRATE_SDF), removeHs=False)
        for mol in suppl:
            if mol is not None:
                _MOL_CACHE[mol.GetProp("_Name")] = mol
    return _MOL_CACHE


def get_substituent_map(mol: str, mol_dir: Path) -> dict:
    """Return {ci, ni, oi, c_aryl, c_alkyl} (1-based) for one substrate, e.g. 'mol_002_E'.

    Derived fresh via RDKit aromaticity (beckmann_core.classical.get_oxime_atoms)
    on a mol loaded live from best_per_substrate.sdf -- not read from any
    pre-computed CSV. Cross-validated against the independently-derived (ci, ni, oi)
    parsed from the molecule's own .gjf title line; a mismatch raises rather than
    silently trusting either source.
    """
    mols = _load_mols()
    if mol not in mols:
        raise ValueError(f"{mol}: not found in {BEST_PER_SUBSTRATE_SDF}")

    result = get_oxime_atoms(mols[mol])
    if result is None:
        raise ValueError(f"{mol}: oxime substructure or aryl/alkyl neighbors not found")
    cox, nox, oox, c_aryl, c_allyl = (idx + 1 for idx in result)

    gjf_ci, gjf_ni, gjf_oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol}_opt.gjf")
    if (cox, nox, oox) != (gjf_ci, gjf_ni, gjf_oi):
        raise ValueError(
            f"{mol}: atom map mismatch -- RDKit gives C{cox}=N{nox}-O{oox}, "
            f".gjf label gives C{gjf_ci}=N{gjf_ni}-O{gjf_oi}"
        )

    return {"ci": cox, "ni": nox, "oi": oox, "c_aryl": c_aryl, "c_alkyl": c_allyl}


# Same flexible atom-matching approach used in parse_nbo.py/parse_cmo.py (\s*
# tolerant of both spaced 'C 6' and compact 'C6' NBO label styles), reimplemented
# here rather than imported to avoid a circular import (parse_cmo.py already
# imports get_substituent_map from this module).
def _label_has_atom(label: str, num: int) -> bool:
    return re.search(rf"[A-Z]+\s*{num}(?!\d)", label) is not None


def compute_psi_row(e2pert_rows: list[dict], ci: int, ni: int, oi: int, c_aryl: int, c_alkyl: int) -> dict:
    """K_anti, K_frag, Psi for one (mol, stage) point's E2PERT rows."""
    k_anti = 0.0
    k_frag = 0.0
    for row in e2pert_rows:
        donor, acceptor, e2 = row["donor"], row["acceptor"], float(row["e2_kcal"])
        if (
            "*" in acceptor
            and _label_has_atom(donor, c_aryl) and _label_has_atom(donor, ci)
            and _label_has_atom(acceptor, ni) and _label_has_atom(acceptor, oi)
        ):
            k_anti += e2
        if (
            "*" in acceptor
            and _label_has_atom(acceptor, ci) and _label_has_atom(acceptor, c_alkyl)
        ):
            k_frag += e2
    psi = k_anti / (k_frag + PSI_EPSILON)
    return {"k_anti": k_anti, "k_frag": k_frag, "psi": psi}


def least_squares_slope(xs: list[float], ys: list[float]) -> float | None:
    """slope = sum((x-mean_x)*(y-mean_y)) / sum((x-mean_x)**2); None if underdetermined."""
    pts = [(x, y) for x, y in zip(xs, ys) if x is not None and y is not None]
    if len(pts) < 2:
        return None
    mean_x = sum(x for x, _ in pts) / len(pts)
    mean_y = sum(y for _, y in pts) / len(pts)
    denom = sum((x - mean_x) ** 2 for x, _ in pts)
    if denom == 0:
        return None
    return sum((x - mean_x) * (y - mean_y) for x, y in pts) / denom


def _read_csv(path: Path) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def build_channel_descriptors(mol: str, mol_dir: Path, e2pert_rows: list[dict], cmo_rows: list[dict]) -> list[dict]:
    """One row per (mol, stage) with k_anti/k_frag/psi joined to the CMO descriptors."""
    subst = get_substituent_map(mol, mol_dir)
    ci, ni, oi, c_aryl, c_alkyl = subst["ci"], subst["ni"], subst["oi"], subst["c_aryl"], subst["c_alkyl"]

    e2pert_by_stage: dict[str, list[dict]] = {}
    for row in e2pert_rows:
        if row["mol"] == mol:
            e2pert_by_stage.setdefault(row["stage"], []).append(row)

    cmo_by_stage = {row["stage"]: row for row in cmo_rows if row["mol"] == mol}

    out = []
    for stage, cmo_row in cmo_by_stage.items():
        psi_row = compute_psi_row(e2pert_by_stage.get(stage, []), ci, ni, oi, c_aryl, c_alkyl)
        out.append({
            "mol": mol, "stage": stage, "r_no": cmo_row["r_no"],
            "c_aryl": c_aryl, "c_alkyl": c_alkyl,
            **psi_row,
            "log_lambda": cmo_row["log_lambda"],
            "wcnmax": cmo_row["wcnmax"], "w17max": cmo_row["w17max"], "w78max": cmo_row["w78max"],
        })
    return out


def compute_slopes(mol: str, channel_rows: list[dict]) -> dict:
    """d/dR for psi, log_lambda, wcnmax, w17max, w78max over the 5-point series."""
    by_stage = {row["stage"]: row for row in channel_rows if row["mol"] == mol}
    series = resolve_series(by_stage)
    r_values = [float(row["r_no"]) for row in series]

    slopes = {"mol": mol, "n_points": len(series), "r_values": r_values}
    for descriptor in ["psi", "log_lambda", "wcnmax", "w17max", "w78max"]:
        y_values = [float(row[descriptor]) if row[descriptor] not in (None, "", "None") else None for row in series]
        slopes[f"d_{descriptor}_dR"] = least_squares_slope(r_values, y_values)
        slopes[f"{descriptor}_values"] = y_values
    return slopes


DESCRIPTORS = ["psi", "log_lambda", "wcnmax", "w17max", "w78max"]


def _float_or_none(v):
    if v in (None, "", "None"):
        return None
    return float(v)


def load_series(mol: str, channel_rows: list[dict]) -> tuple[list[float], dict[str, list[float | None]]]:
    """R(N-O) values and per-descriptor y-values for one molecule's series,
    R0 (if present) followed by scan_1..scan_N in order (see resolve_series()).
    Reused by beckmann_nbo.viz and its analysis scripts, not just one script."""
    by_stage = {row["stage"]: row for row in channel_rows if row["mol"] == mol}
    series = resolve_series(by_stage)
    r_values = [float(row["r_no"]) for row in series]
    y_by_descriptor = {d: [_float_or_none(row[d]) for row in series] for d in DESCRIPTORS}
    return r_values, y_by_descriptor


def main() -> None:
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    analysis_dir = DATA_OUTPUT / "analysis"
    e2pert_rows = _read_csv(analysis_dir / "nbo_e2pert.csv")
    cmo_rows    = _read_csv(analysis_dir / "cmo_descriptors.csv")

    all_channel_rows = []
    all_slopes = []
    for mol_id in sorted(ALL_IDS):
        mol = resolve_mol_name(mol_id, dft_opt_dir)
        if mol is None:
            print(f"-- mol_{mol_id.zfill(3)}: no directory, skipping")
            continue
        mol_dir = dft_opt_dir / mol
        channel_rows = build_channel_descriptors(mol, mol_dir, e2pert_rows, cmo_rows)
        if not channel_rows:
            # Upstream parse_nbo/parse_cmo excluded this molecule entirely (e.g. a
            # stage log that didn't reach Normal termination -- see JOB_ISSUES.md).
            print(f"-- {mol}: no rows in nbo_e2pert.csv/cmo_descriptors.csv, skipping")
            continue
        all_channel_rows.extend(channel_rows)
        slopes = compute_slopes(mol, channel_rows)
        all_slopes.append(slopes)
        print(f"-- {mol}: {len(channel_rows)} stage points, {slopes['n_points']} in the d/dR series")
        for descriptor in ["psi", "log_lambda", "wcnmax", "w17max", "w78max"]:
            print(f"     d({descriptor})/dR = {slopes[f'd_{descriptor}_dR']}")

    channel_path = analysis_dir / "channel_descriptors.csv"
    with open(channel_path, "w", newline="") as f:
        fields = ["mol", "stage", "r_no", "c_aryl", "c_alkyl", "k_anti", "k_frag", "psi",
                  "log_lambda", "wcnmax", "w17max", "w78max"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_channel_rows)
    print(f"\n{len(all_channel_rows)} rows -> {channel_path}")

    slopes_path = analysis_dir / "descriptor_slopes.csv"
    with open(slopes_path, "w", newline="") as f:
        fields = ["mol", "n_points", "r_values"] + [
            f"{prefix}{descriptor}{suffix}"
            for descriptor in ["psi", "log_lambda", "wcnmax", "w17max", "w78max"]
            for prefix, suffix in [("d_", "_dR"), ("", "_values")]
        ]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_slopes)
    print(f"{len(all_slopes)} rows -> {slopes_path}")


if __name__ == "__main__":
    main()
