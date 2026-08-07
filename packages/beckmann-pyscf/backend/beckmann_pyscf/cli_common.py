"""Small output-writing helpers shared by cli_predict.py/cli_scan.py --
CLI-layer formatting only, no pipeline science lives here."""
import csv
from pathlib import Path

SERIES_FIELDS = ["stage", "R_NO", "weight", "MO_index"]


def write_series_csv(series: list[dict], path: Path) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SERIES_FIELDS)
        writer.writeheader()
        for row in series:
            writer.writerow({k: row.get(k) for k in SERIES_FIELDS})


def write_summary(path: Path, *, name: str, prediction: str, minimum: dict | None,
                   energy_ev: float | None = None) -> str:
    """Writes path and returns the same text (so callers can also print it)."""
    lines = [f"{name}: prediction = {prediction}"]
    if minimum is not None:
        lines.append(
            f"  interior wCNmax minimum: R* = {minimum['R_star']:.4f} A, "
            f"w* = {minimum['w_star']:.4f}, MO{minimum['MO_index']}"
        )
    else:
        lines.append("  no interior wCNmax minimum found")
    if energy_ev is not None:
        lines.append(f"  AIMNet2 energy: {energy_ev:.6f} eV")
    text = "\n".join(lines) + "\n"
    path.write_text(text)
    return text
