"""
Gaussian/NBO7/Citadel-specific settings, plus this package's own data-path
convention. CHARGE/MULTIPLICITY moved to beckmann_core.constants (shared
with beckmann-pyscf) -- import them from there, not from here.

The official basis set for every DFT job in this project is
wB97XD/6-311+G(d,p) -- see the monorepo README for why.

ROOT resolution: unlike the old single-repo layout (where beckmann/config.py
sat directly under the repo root and `Path(__file__).parent.parent` just
worked), this package now lives under packages/beckmann-nbo/beckmann_nbo/,
two directories further from the shared data/ tree. Rather than hardcoding
a fixed number of `.parent` hops (brittle if this package is ever installed
or relocated independently), ROOT is resolved by searching upward from this
file for the nearest ancestor directory that actually contains a data/
directory -- overridable via the BECKMANN_DATA_ROOT environment variable for
any setup where that search wouldn't find the right place.
"""
import os
from pathlib import Path


def _find_data_root() -> Path:
    env_override = os.environ.get("BECKMANN_DATA_ROOT")
    if env_override:
        return Path(env_override).resolve()

    here = Path(__file__).resolve()
    for ancestor in [here.parent, *here.parents]:
        if (ancestor / "data").is_dir():
            return ancestor
    raise RuntimeError(
        "Could not locate a data/ directory above "
        f"{here} -- set BECKMANN_DATA_ROOT explicitly."
    )


ROOT = _find_data_root()

DATA_INPUT  = ROOT / "data" / "input"
DATA_OUTPUT = ROOT / "data" / "output"

FUNCTIONAL   = "wB97XD"
BASIS        = "6-311+G(d,p)"
NPROC        = 8
MEM_GB       = 16
NBO_KEYWORDS = "E2PERT BNDIDX NBOSUM CMO"
SOLVENT      = "scrf=(smd,solvent=water)"
