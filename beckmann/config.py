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
# CMO (Lambda, wCNmax) requires NBO7 -- every stage below must use
# pop=nbo7read (not pop=nboread) in its route line, which routes through
# Gaussian's external-program interface (gaunbo7 -> g16nbo -> nbo7.i8.exe)
# instead of the bundled NBO 3.1, which doesn't support CMO at all.
NBO_KEYWORDS = "E2PERT BNDIDX NBOSUM CMO"
# Per supervisor: mimic experimental (aqueous) conditions in every calculation
# -- no gas-phase runs. SMD/water, added to every route line below.
SOLVENT = "scrf=(smd,solvent=water)"
