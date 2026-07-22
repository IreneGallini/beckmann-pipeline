"""
Diagnostic 2 of 2 (Notes_open_source_alt.md, "wCNmax: PySCF vs NBO7"): does adding
deflation to the local per-atom-pair construction recover the dip depth at mol_006_E's
R=1.6608 A avoided-crossing point? One-off script, not part of the regular
wcnmax_scan_rule.py pipeline.

Geometry source: the REAL Citadel fine-resolution rerun,
data/output/dft_opt_finescan/mol_006_E_finescan/mol_006_E_finescan_scan.log -- NOT the
earlier interpolated Cartesian blends (removed from this branch; see
Notes_open_source_alt.md's "mol_006_E follow-up" section for why those were rejected).

Trusted NBO7 value at this point: 0.32604 (data/output/analysis/cmo_channel_extraction.csv,
mol_006_E / scan_3 / cn channel -- the same canonical R=1.6608 point, confirmed to match
the user-cited 0.3260 exactly).

Deflation pairs chosen: the two other bonds directly attached to the oxime carbon
(C-aryl, C-alkyl) plus the N-O bond directly attached to the oxime nitrogen -- the three
bonds most likely to leak density into the local C-N IAO block, approximating NBO's own
sequential whole-molecule deflation without building every bond in the molecule.

Writes data/output/analysis/wcnmax_deflation_diagnostic_mol006.csv (new file -- does not
touch wcnmax_rule_results_opensource.csv or anything under data/output/dft_opt/).
"""
import csv

from beckmann.config import DATA_OUTPUT
from beckmann.dft.descriptors import get_substituent_map
from beckmann.dft.scan import oxime_atom_map_from_gjf

import numpy as np

from beckmann_alt.geometry import CHARGE, SPIN, _stage_points_from_log, pyscf_atom_spec
from beckmann_alt.pair_nbo import (
    build_local_iaos, compute_wcnmax, compute_wcnmax_deflated, deflated_density_matrix, pair_density_matrix,
)
from beckmann_alt.pyscf_livvo import build_mol, run_scf

MOL = "mol_006_E"
MOL_DIR = DATA_OUTPUT / "dft_opt" / MOL
FINESCAN_LOG = DATA_OUTPUT / "dft_opt_finescan" / f"{MOL}_finescan" / f"{MOL}_finescan_scan.log"
TARGET_R = 1.6608
TRUSTED_NBO7_WEIGHT = 0.32604  # cmo_channel_extraction.csv, mol_006_E/scan_3/cn
OUT_CSV = DATA_OUTPUT / "analysis" / "wcnmax_deflation_diagnostic_mol006.csv"


def main() -> None:
    ci, ni, oi, _ = oxime_atom_map_from_gjf(MOL_DIR / f"{MOL}_opt.gjf")
    sub = get_substituent_map(MOL, MOL_DIR)
    c_aryl, c_alkyl = sub["c_aryl"], sub["c_alkyl"]
    print(f"{MOL} atom map: ci={ci} ni={ni} oi={oi} c_aryl={c_aryl} c_alkyl={c_alkyl}")

    by_r = _stage_points_from_log(FINESCAN_LOG, ni, oi)
    print(f"finescan R points available: {sorted(by_r)}")
    matches = [r for r in by_r if abs(r - TARGET_R) < 1e-3]
    if not matches:
        raise ValueError(f"no finescan point within 1e-3 of R={TARGET_R}; available: {sorted(by_r)}")
    r_actual = matches[0]
    atoms = by_r[r_actual]
    print(f"using finescan geometry at R={r_actual}")

    case = {
        "name": f"{MOL}_finescan_R{r_actual}",
        "atoms": atoms,
        "atom_spec": pyscf_atom_spec(atoms),
        "charge": CHARGE,
        "spin": SPIN,
        "ci": ci, "ni": ni,
    }

    mol = build_mol(case)
    mf = run_scf(mol)  # production defaults: ddCOSMO/water + density-fit
    iaos_orth, s, atom_of_iao = build_local_iaos(mf)

    before = compute_wcnmax(mf, s, iaos_orth, atom_of_iao, ci, ni)
    deflate_pairs = [(ci, c_aryl), (ci, c_alkyl), (ni, oi)]
    after = compute_wcnmax_deflated(mf, s, iaos_orth, atom_of_iao, ci, ni, deflate_pairs)

    dm_ao = pair_density_matrix(mf)
    dm_deflated = deflated_density_matrix(dm_ao, s, iaos_orth, atom_of_iao, deflate_pairs)
    print(f"\n[debug] full dm perturbation norm ||dm_deflated - dm_ao||_F = {np.linalg.norm(dm_deflated - dm_ao):.6f}")
    a0, b0 = ci - 1, ni - 1
    cn_idx = np.where((atom_of_iao == a0) | (atom_of_iao == b0))[0]
    cn_iaos = iaos_orth[:, cn_idx]
    block_before = cn_iaos.T @ s @ dm_ao @ s @ cn_iaos
    block_after = cn_iaos.T @ s @ dm_deflated @ s @ cn_iaos
    print(f"[debug] C-N local block perturbation norm = {np.linalg.norm(block_after - block_before):.8f}")
    print(f"[debug] before wmax (full precision) = {before['wmax']!r}, MO={before['mo_index']}")
    print(f"[debug] after  wmax (full precision) = {after['wmax']!r}, MO={after['mo_index']}")
    print(f"[debug] before candidate occupations = {before['pair']['candidate_occupations']}")
    print(f"[debug] after  candidate occupations = {after['pair']['candidate_occupations']}")

    rows = [
        {
            "condition": "before_deflation", "R_NO": r_actual,
            "wCNmax": round(before["wmax"], 4), "MO_index_0based": before["mo_index"],
            "pct_of_NBO7": round(100 * before["wmax"] / TRUSTED_NBO7_WEIGHT, 1),
        },
        {
            "condition": "after_deflation", "R_NO": r_actual,
            "wCNmax": round(after["wmax"], 4), "MO_index_0based": after["mo_index"],
            "pct_of_NBO7": round(100 * after["wmax"] / TRUSTED_NBO7_WEIGHT, 1),
        },
        {
            "condition": "NBO7_trusted", "R_NO": r_actual,
            "wCNmax": TRUSTED_NBO7_WEIGHT, "MO_index_0based": None, "pct_of_NBO7": 100.0,
        },
    ]
    for r in rows:
        print(r)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {OUT_CSV}")


if __name__ == "__main__":
    main()
