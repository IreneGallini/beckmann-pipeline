"""
One-off helper: compute beckmann_alt.pair_nbo.run_test_set_scan_series(mol_id)
and cache the result rows as JSON, so the expensive PySCF run (multiple
minutes per R(N-O) point) only has to happen once per molecule -- run in
parallel, one process per molecule, from the shell.

Usage: python -m beckmann_alt._compute_scan_cache <mol_id> <out_json_path>
"""
import json
import sys

from beckmann_alt.pair_nbo import run_test_set_scan_series


def main() -> None:
    mol_id, out_path = sys.argv[1], sys.argv[2]
    rows = run_test_set_scan_series(mol_id)
    with open(out_path, "w") as f:
        json.dump(rows, f, indent=2)
    print(f"{mol_id}: {len(rows)} points -> {out_path}")


if __name__ == "__main__":
    main()
