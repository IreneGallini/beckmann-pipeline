"""
Step 1 of the E/Z re-run: for a given set of molecules, determine which
AIMNet2-optimized isomer (E or Z) is actually lower-energy, using energies
already computed by the benchmark batch pipeline -- no new conformer
generation or AIMNet2 optimization needed.

Flags any molecule where the isomer already carried through NBO/PySCF is
*not* the lower-energy one -- a separate finding worth surfacing on its own,
per Olexandr's ask.

ALREADY_RUN (mol_id -> isomer already run through NBO/PySCF) is derived
dynamically from wcnmax_rule_results_pyscf.csv's "mol" column -- the
canonical record of which isomer suffix was used for each of the full
34-molecule benchmark set -- rather than hand-maintained, so it covers any
subset of the benchmark without editing this file.
"""
import csv
import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated")
warnings.filterwarnings("ignore", category=DeprecationWarning)

from rdkit import Chem

from beckmann_nbo.config import DATA_OUTPUT

ENERGY_CSV = DATA_OUTPUT / "analysis" / "ez_energy_comparison.csv"
ENERGY_FIELDS = ["mol_id", "E_energy_eV", "Z_energy_eV", "diff_eV"]
EXISTING_RESULTS_CSV = DATA_OUTPUT / "analysis" / "wcnmax_rule_results_pyscf.csv"


def _already_run_map(csv_path=EXISTING_RESULTS_CSV) -> dict[str, str]:
    """{'mol_001': 'E', 'mol_014': 'Z', ...} for all 34 benchmark molecules."""
    out = {}
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            parent, ez = row["mol"].rsplit("_", 1)
            out[parent] = ez
    return out


ALREADY_RUN = _already_run_map()


def energy_table(sdf_path=DATA_OUTPUT / "aimnet_optimized" / "best_aimnet_optimized.sdf",
                  mol_ids=tuple(ALREADY_RUN)) -> list[dict]:
    """Returns one row per molecule:
    {mol_id, E_isomer_eV, Z_isomer_eV, lowest_energy_isomer, other_isomer,
     already_run_isomer, already_run_is_higher_energy}
    """
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    energies: dict[str, dict[str, float]] = {}
    for mol in suppl:
        if mol is None:
            continue
        name = mol.GetProp("_Name")
        parent, ez = name.rsplit("_", 1)
        if parent not in mol_ids or not mol.HasProp("E_aimnet2_eV"):
            continue
        energies.setdefault(parent, {})[ez] = float(mol.GetProp("E_aimnet2_eV"))

    rows = []
    for mol_id in mol_ids:
        e = energies.get(mol_id, {})
        if "E" not in e or "Z" not in e:
            rows.append({"mol_id": mol_id, "E_isomer_eV": e.get("E"), "Z_isomer_eV": e.get("Z"),
                         "lowest_energy_isomer": None, "other_isomer": None,
                         "already_run_isomer": ALREADY_RUN[mol_id],
                         "already_run_is_higher_energy": None})
            continue
        lowest = "E" if e["E"] < e["Z"] else "Z"
        other = "Z" if lowest == "E" else "E"
        already_run = ALREADY_RUN[mol_id]
        rows.append({
            "mol_id": mol_id,
            "E_isomer_eV": e["E"], "Z_isomer_eV": e["Z"],
            "lowest_energy_isomer": lowest, "other_isomer": other,
            "already_run_isomer": already_run,
            "already_run_is_higher_energy": already_run != lowest,
        })
    return rows


def write_energy_csv(rows: list[dict], out_path=ENERGY_CSV) -> None:
    """mol_id, E_energy_eV, Z_energy_eV, diff_eV (= E - Z) -- one line per
    molecule, the AIMNet2 energies auto3d's best-isomer selection compares."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=ENERGY_FIELDS)
        writer.writeheader()
        for r in rows:
            e, z = r["E_isomer_eV"], r["Z_isomer_eV"]
            writer.writerow({
                "mol_id": r["mol_id"],
                "E_energy_eV": e,
                "Z_energy_eV": z,
                "diff_eV": (e - z) if (e is not None and z is not None) else "",
            })
    print(f"Wrote {out_path}")


def main(mol_ids=tuple(ALREADY_RUN), out_path=ENERGY_CSV, write_csv=True) -> list[dict]:
    rows = energy_table(mol_ids=mol_ids)
    if write_csv:
        write_energy_csv(rows, out_path=out_path)
    print(f"{'mol_id':<10} {'E (eV)':>16} {'Z (eV)':>16}  {'lowest':>6}  {'already-run':>11}  flag")
    print("-" * 75)
    flagged = []
    for r in rows:
        if r["lowest_energy_isomer"] is None:
            print(f"{r['mol_id']:<10} missing E or Z energy -- skipping")
            continue
        flag = "<-- already-run isomer is HIGHER energy" if r["already_run_is_higher_energy"] else ""
        if flag:
            flagged.append(r["mol_id"])
        print(f"{r['mol_id']:<10} {r['E_isomer_eV']:>16.6f} {r['Z_isomer_eV']:>16.6f}  "
              f"{r['lowest_energy_isomer']:>6}  {r['already_run_isomer']:>11}  {flag}")
    print()
    if flagged:
        print(f"FINDING: already-run isomer is the higher-energy one for: {', '.join(flagged)}")
    else:
        print(f"FINDING: already-run isomer is the lower-energy one for all {len(rows)} molecules.")
    return rows


if __name__ == "__main__":
    main()
