"""
Reference geometries and atom maps for the two validation cases, reusing the main
pipeline's own Gaussian-log geometry parser rather than reimplementing one.

mol_002_E: our own pipeline, wB97XD/6-311+G(d,p), SMD/water -- fully known level of
    theory. Geometry is the DFT-converged one from Stage 1 (mol_002_E_opt.log's final
    "Standard orientation"), the same geometry Stage 2/3's NBO7 analysis ran on.
5_s0_Me: Tetiana's external reference log (compound 3, "Ring Size and Substituent
    Effects in the Beckmann Rearrangement", Table 2). Its own route line
    (wb97xd/genecp scrf=(smd,solvent=water)) uses a custom hand-specified ECP basis
    applied to every C/N/O center (4 valence electrons per carbon, minimal
    split-valence primitives) -- NOT our all-electron 6-311+G(d,p). Reproducing that
    exact basis means transcribing every exponent/coefficient/ECP parameter from the
    log by hand; not attempted here. We run this case at OUR basis (6-311+G(d,p))
    instead, so any comparison against the paper's reported numbers is a
    different-basis check, not an apples-to-apples validation -- see
    Notes_open_source_alt.md.
"""
from pathlib import Path

from beckmann.dft.scan import ATOMIC_SYMBOLS, parse_standard_orientations

ROOT = Path(__file__).parent.parent

MOL_002_OPT_LOG = ROOT / "data" / "output" / "dft_opt" / "mol_002_E" / "mol_002_E_opt.log"
REFERENCE_LOG   = ROOT / "5_s0_Me.log"


def _parse_orientation_blocks(lines: list[str], header: str) -> list[tuple[int, list]]:
    """Same block-parsing loop as beckmann.dft.scan.parse_standard_orientations, just
    parameterized on the header string -- needed because 5_s0_Me.log's route line uses
    `nosymm`, which makes Gaussian print 'Input orientation:' instead of 'Standard
    orientation:'. Our own pipeline's .gjf files never use nosymm, so
    parse_standard_orientations (hardcoded to 'Standard orientation:') is correct and
    unmodified for every log the main pipeline produces -- this is a local fallback for
    this one external reference file's quirk, not a replacement for it."""
    blocks = []
    i = 0
    while i < len(lines):
        if header in lines[i]:
            j = i + 5
            atoms = []
            while j < len(lines) and "---" not in lines[j]:
                parts = lines[j].split()
                if len(parts) == 6:
                    atomic_num = int(parts[1])
                    x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                    sym = ATOMIC_SYMBOLS.get(atomic_num, f"X{atomic_num}")
                    atoms.append((sym, x, y, z))
                j += 1
            if atoms:
                blocks.append((i, atoms))
            i = j
        else:
            i += 1
    return blocks


def final_geometry(log_path: Path) -> list[tuple]:
    """Atoms (symbol, x, y, z) from the last orientation block in a Gaussian log.

    Tries 'Standard orientation:' first (via the main pipeline's own parser); falls
    back to 'Input orientation:' for logs run with nosymm (see _parse_orientation_blocks).
    """
    lines = log_path.read_text().splitlines()
    blocks = parse_standard_orientations(lines)
    if not blocks:
        blocks = _parse_orientation_blocks(lines, "Input orientation:")
    if not blocks:
        raise ValueError(f"{log_path}: no 'Standard orientation'/'Input orientation' block found")
    _, atoms = blocks[-1]
    return atoms


def pyscf_atom_spec(atoms: list[tuple]) -> list[list]:
    """Convert (symbol, x, y, z) tuples into PySCF's Mole.atom list format."""
    return [[sym, (x, y, z)] for sym, x, y, z in atoms]


# (charge, multiplicity) -- both cases are the protonated activated oxime, singlet.
CHARGE = 1
MULTIPLICITY = 1
SPIN = MULTIPLICITY - 1  # PySCF wants 2S, not 2S+1

REFERENCE_CASES = {
    "mol_002": {
        "log": MOL_002_OPT_LOG,
        "ci": 11, "ni": 12, "oi": 13, "c_aryl": 6, "c_alkyl": 10,
        "basis_note": "our own pipeline's 6-311+G(d,p), all-electron -- exact match to config.py",
    },
    "5_s0_Me": {
        "log": REFERENCE_LOG,
        # Hardcoded from the paper's own Figure 2 / compound-3 convention -- same
        # atom map validate_reference_descriptors.py uses for the NBO7 check.
        "ci": 7, "ni": 17, "oi": 18, "c_aryl": 1, "c_alkyl": 8,
        "basis_note": (
            "run at OUR 6-311+G(d,p), not the log's actual custom GenECP basis -- "
            "NOT an apples-to-apples reproduction of the paper's numbers, see module docstring"
        ),
    },
}


def load_case(name: str) -> dict:
    """Atoms, charge/spin, and atom map for one of the two reference cases."""
    case = REFERENCE_CASES[name]
    atoms = final_geometry(case["log"])
    return {
        "name": name,
        "atoms": atoms,
        "atom_spec": pyscf_atom_spec(atoms),
        "charge": CHARGE,
        "spin": SPIN,
        "ci": case["ci"], "ni": case["ni"], "oi": case["oi"],
        "c_aryl": case["c_aryl"], "c_alkyl": case["c_alkyl"],
        "basis_note": case["basis_note"],
    }
