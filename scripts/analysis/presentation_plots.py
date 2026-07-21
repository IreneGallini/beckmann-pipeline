"""
Presentation-quality wCNmax figures for a small, curated subset of
substrates -- for a final deck, not data exploration. Reuses the same data
loaders as scripts/analysis/summarize_descriptors.py
(beckmann.dft.descriptors.load_series()/find_wcnmax_extremum()) and the same
reusable chart function (beckmann.dft.viz.plot_wcnmax_single()) so there is
one code path for "get a molecule's wCNmax curve," not a second copy tuned
for looks.

Kept separate from summarize_descriptors.py rather than a `--presentation`
flag on it: the input (a handful of hand-picked substrates vs. all 34) and
the concern (one-off deck polish vs. reproducible exploration) genuinely
differ.

Curated subset (chosen to tell the project's actual story, not arbitrarily):
  mol_006_E  -- the paper's central signature: a real interior wCNmax
                minimum (see Notes.md); also the project's own reference
                compound.
  mol_019_E  -- a decisive rearrangement case (100% product A).
  mol_002_E  -- a decisive fragmentation case (100% product B); the
                project's original/most-discussed reference substrate.
  mol_009_E  -- a borderline case (56/44 R/F) -- the kind of substrate the
                binary R/F label alone doesn't tell a clean story about.

Output: data/output/analysis/plots/presentation_{mol}.png
"""
import csv
import json

import matplotlib.pyplot as plt

from beckmann.config import DATA_INPUT, DATA_OUTPUT
from beckmann.dft.descriptors import find_wcnmax_extremum, load_series
from beckmann.dft.viz import plot_wcnmax_single

ANALYSIS_DIR = DATA_OUTPUT / "analysis"
PLOTS_DIR    = ANALYSIS_DIR / "plots"

CURATED_MOLS = ["mol_006_E", "mol_019_E", "mol_002_E", "mol_009_E"]


def _read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def main() -> None:
    channel_rows    = _read_csv(ANALYSIS_DIR / "channel_descriptors.csv")
    extraction_rows = _read_csv(ANALYSIS_DIR / "cmo_channel_extraction.csv")
    outcomes        = json.loads((DATA_INPUT / "benchmark_meta.json").read_text())

    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({
        "font.size": 13,
        "axes.titlesize": 14,
        "axes.labelsize": 13,
    })

    for mol in CURATED_MOLS:
        r_values, y_by_descriptor = load_series(mol, channel_rows)
        if not r_values:
            print(f"-- {mol}: no data in channel_descriptors.csv, skipping")
            continue
        extremum = find_wcnmax_extremum(mol, extraction_rows)
        mol_id = mol.split("_")[1]
        outcome = outcomes[f"mol_{mol_id}"]["exp_outcome"]

        fig = plot_wcnmax_single(mol, r_values, y_by_descriptor["wcnmax"], extremum, outcome)
        out_path = PLOTS_DIR / f"presentation_{mol}.png"
        fig.savefig(out_path, dpi=200)
        plt.close(fig)
        print(f"-- wrote {out_path}")


if __name__ == "__main__":
    main()
