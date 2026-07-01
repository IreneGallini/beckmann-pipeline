from pathlib import Path

ROOT = Path(__file__).parent.parent

DATA_INPUT  = ROOT / "data" / "input"
DATA_OUTPUT = ROOT / "data" / "output"

# Gaussian / DFT settings — single source of truth for all prepare_*.py scripts
FUNCTIONAL   = "wB97XD"
BASIS        = "6-311+G(d,p)"
NPROC        = 8
MEM_GB       = 16
CHARGE       = 1        # protonated activated oxime (C=N-[OH2+])
MULTIPLICITY = 1
NBO_KEYWORDS_EQ   = "E2PERT BNDIDX NBOSUM"   # equilibrium NBO via NBO 3.1 bundled in g16
NBO_KEYWORDS_SCAN = "E2PERT BNDIDX NBOSUM"   # scan NBO via NBO 3.1
# CMO is NBO7-only; obtained post-hoc via gennbo on the .47 archive from _nbo.gjf
NBO_KEYWORDS_SP   = "E2PERT BNDIDX NBOSUM CMO"  # single-point (requires NBO7)
