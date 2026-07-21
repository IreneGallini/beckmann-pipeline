"""
Compare beckmann_alt.pair_nbo's computed wCNmax against the trusted NBO7 wCNmax for
every main-pipeline test-set molecule (mol_002/006/014/020/021/029).

Trusted numbers are read directly from the already-generated
data/output/analysis/cmo_channel_extraction.csv (nbo/R0 stage, channel 'cn') --
never recomputed by a different method here. See Notes_open_source_alt.md,
"wCNmax across all 6 test-set molecules", for the write-up of the result this
reproduces.
"""
import csv

from beckmann.config import DATA_OUTPUT
from beckmann.dft.inputs import TEST_IDS

from beckmann_alt.pair_nbo import run_test_set_case


def trusted_wcnmax(mol_name: str) -> tuple[float, int]:
    """(wmax, MO_index) for mol_name's nbo/R0-stage cn channel from the trusted CSV."""
    path = DATA_OUTPUT / "analysis" / "cmo_channel_extraction.csv"
    for row in csv.DictReader(open(path)):
        if row["mol"] == mol_name and row["stage"] == "nbo" and row["channel"] == "cn":
            return float(row["weight"]), int(row["MO_index"])
    raise ValueError(f"{mol_name}: no nbo/cn row found in {path}")


def main() -> None:
    header = f"{'mol':<12} {'computed wCNmax':>16} {'MO':>4}   {'trusted wCNmax':>15} {'MO':>4}   {'% diff':>7} {'MO offset':>9}"
    print(header)
    print("-" * len(header))
    for mol_id in sorted(TEST_IDS):
        result = run_test_set_case(mol_id)
        mol_name = result["case"]
        computed_wmax, computed_mo = result["cn"]["wmax"], result["cn"]["mo_index"]
        trusted_wmax, trusted_mo = trusted_wcnmax(mol_name)

        pct_diff = 100 * (computed_wmax - trusted_wmax) / trusted_wmax
        mo_offset = computed_mo - trusted_mo

        print(
            f"{mol_name:<12} {computed_wmax:16.4f} {computed_mo:4d}   "
            f"{trusted_wmax:15.4f} {trusted_mo:4d}   {pct_diff:+6.1f}% {mo_offset:+9d}"
        )


if __name__ == "__main__":
    main()
