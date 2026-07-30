"""
Build the diabatic character-exchange dataset across the full benchmark set
(beckmann.dft.diabatic_character.build_benchmark_character_exchange()). Not
part of the regular prediction pipeline -- for supervisor discussion,
extending the reference-case result (scripts/analysis/plot_character_exchange.py)
to real substrates. No classification threshold is computed here; that's
deliberately deferred.

Produces:
  data/output/analysis/character_exchange_benchmark.csv         -- one row
      per (molecule, scan point)
  data/output/analysis/character_exchange_benchmark_summary.csv -- one row
      per molecule
"""
import csv

from beckmann.config import DATA_OUTPUT
from beckmann.dft.diabatic_character import build_benchmark_character_exchange

ANALYSIS_DIR = DATA_OUTPUT / "analysis"


def write_csv(path, rows: list[dict]) -> None:
    if not rows:
        print(f"-- skipped {path} (no rows)")
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"-- wrote {path} ({len(rows)} rows)")


def main() -> None:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
    detail_rows, summary_rows, failures = build_benchmark_character_exchange()

    write_csv(ANALYSIS_DIR / "character_exchange_benchmark.csv", detail_rows)
    write_csv(ANALYSIS_DIR / "character_exchange_benchmark_summary.csv", summary_rows)

    print(f"\n{len(summary_rows)} molecule(s) succeeded, {len(failures)} failed")
    if failures:
        by_category: dict[str, list[str]] = {}
        for f in failures:
            by_category.setdefault(f["category"], []).append(f["mol"])
        for category, mols in sorted(by_category.items()):
            print(f"  {category} ({len(mols)}): {', '.join(mols)}")
            for f in failures:
                if f["category"] == category:
                    print(f"    {f['mol']}: {f['detail']}")


if __name__ == "__main__":
    main()
