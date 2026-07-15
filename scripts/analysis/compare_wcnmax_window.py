"""
Test whether restricting wCNmax's search to a LUMO-to-LUMO+0.4 a.u. window
(matching the paper's own method) reveals the interior minimum reported for
mol_006_E, which the current unrestricted (full virtual manifold) search does
not show.

Standalone re-analysis of already-downloaded log data -- no new Gaussian jobs,
no changes to the pipeline. Reuses parse_cmo.py's parse_cmo_table()/
virtual_window()/max_weight_for_target()/r_no_before() directly rather than
duplicating them.
"""
from pathlib import Path

from beckmann.config import DATA_OUTPUT
from beckmann.dft import parse_cmo
from beckmann.dft.descriptors import get_substituent_map
from beckmann.dft.parse_nbo import r_no_before
from beckmann.dft.scan import oxime_atom_map_from_gjf

OLD_DIR = DATA_OUTPUT / "dft_opt"
NEW_DIR = DATA_OUTPUT / "dft_opt_rigidscan"


def windowed_vs_unrestricted(log_path: Path, ci: int, ni: int) -> list[dict]:
    """For every CMO table in log_path, compute wCNmax both the current
    (unrestricted, full virtual manifold) way and restricted to LUMO..LUMO+0.4 a.u."""
    lines = log_path.read_text().splitlines()
    starts = parse_cmo.find_cmo_sections(lines)
    rows = []
    for start in starts:
        table = parse_cmo.parse_cmo_table(lines, start)
        if not table:
            continue
        vir_all = [m for m in table if m["kind"] == "vir"]
        if not vir_all:
            continue
        lumo_e = vir_all[0]["energy"]
        window = parse_cmo.virtual_window(table, window_au=0.4)

        unrestricted = parse_cmo.max_weight_for_target(vir_all, ci, ni, lumo_e)
        windowed = parse_cmo.max_weight_for_target(window, ci, ni, lumo_e)

        rows.append({
            "unrestricted_weight": unrestricted[0], "unrestricted_mo": unrestricted[1],
            "windowed_weight": windowed[0], "windowed_mo": windowed[1],
            "_start": start,
        })
    return rows, lines


def run(label: str, mol: str, log_paths: list[tuple[str, Path]]) -> None:
    """label is just for display; mol must be the real molecule dir name
    (e.g. 'mol_006_E') so its _opt.gjf can be found under OLD_DIR."""
    old_mol_dir = OLD_DIR / mol
    ci, ni, oi, _ = oxime_atom_map_from_gjf(old_mol_dir / f"{mol}_opt.gjf")

    print(f"\n{'='*90}\n{label}  (ci={ci}, ni={ni}, oi={oi})\n{'='*90}")
    print(f"{'stage':<12}{'R(N-O)':>10}  {'unrestricted':>14} {'(MO)':>6}   {'windowed':>10} {'(MO)':>6}   same-MO?")
    for label, log_path in log_paths:
        if not log_path.exists():
            print(f"{label:<12} -- missing: {log_path.name}")
            continue
        rows, lines = windowed_vs_unrestricted(log_path, ci, ni)
        for row in rows:
            r_no = r_no_before(lines, row["_start"], ni, oi)
            uw, um = row["unrestricted_weight"], row["unrestricted_mo"]
            ww, wm = row["windowed_weight"], row["windowed_mo"]
            same = "yes" if um == wm else "NO"
            uw_s = f"{uw:.4f}" if uw is not None else "None"
            ww_s = f"{ww:.4f}" if ww is not None else "None"
            r_s = f"{r_no:.4f}" if r_no is not None else "?"
            print(f"{label:<12}{r_s:>10}  {uw_s:>14} {str(um):>6}   {ww_s:>10} {str(wm):>6}   {same}")


def main() -> None:
    # Sanity check first: mol_002_E / mol_020_E, where the windowed and
    # unrestricted searches are already known to agree (in_window=True at
    # every point in the existing cmo_channel_extraction.csv) -- confirms this
    # re-implementation is correct before trusting what it says about mol_006.
    for mol in ("mol_002_E", "mol_020_E"):
        d = OLD_DIR / mol
        run(mol, mol, [
            ("nbo", d / f"{mol}_nbo.log"),
            ("sp2", d / f"{mol}_sp2.log"),
            ("sp3", d / f"{mol}_sp3.log"),
            ("sp4", d / f"{mol}_sp4.log"),
            ("scan", d / f"{mol}_scan.log"),
        ])

    # The actual question: mol_006_E, old architecture (full 4-stretched-point
    # series) and new rigid-scan architecture (3 of 4 points; pt4 crashed).
    d_old = OLD_DIR / "mol_006_E"
    run("mol_006_E [OLD architecture]", "mol_006_E", [
        ("nbo", d_old / "mol_006_E_nbo.log"),
        ("sp2", d_old / "mol_006_E_sp2.log"),
        ("sp3", d_old / "mol_006_E_sp3.log"),
        ("sp4", d_old / "mol_006_E_sp4.log"),
        ("scan", d_old / "mol_006_E_scan.log"),
    ])

    d_new = NEW_DIR / "mol_006_E_rigidscan"
    run("mol_006_E [NEW rigid-scan]", "mol_006_E", [
        ("nbo", d_old / "mol_006_E_nbo.log"),  # R0 baseline unchanged, reused
        ("rigidscan", d_new / "mol_006_E_rigidscan_scan.log"),
    ])


if __name__ == "__main__":
    main()
