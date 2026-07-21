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

from beckmann.dft.inputs import STEP_SCAN_SOURCES, TEST_IDS, resolve_mol_name, step_scan_dir
from beckmann.dft.parse_cmo import find_cmo_sections
from beckmann.dft.scan import (
    ATOMIC_SYMBOLS, no_distance, oxime_atom_map_from_gjf, parse_standard_orientations,
)
from beckmann.config import DATA_OUTPUT

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


def load_test_set_case(mol_id: str) -> dict:
    """Same shape as load_case(), but for any of the six main-pipeline test-set
    molecules (mol_002/006/014/020/021/029) instead of the two hand-picked reference
    cases above. Atom map (ci/ni) and geometry are both resolved fresh via the main
    pipeline's own utilities -- resolve_mol_name (beckmann.dft.inputs) and
    oxime_atom_map_from_gjf (beckmann.dft.scan) -- not hardcoded per molecule the way
    REFERENCE_CASES is, since that would mean re-transcribing six atom maps by hand
    with the usual risk of a transcription error.
    """
    if mol_id not in TEST_IDS:
        raise ValueError(f"{mol_id}: not in TEST_IDS ({sorted(TEST_IDS)})")
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    mol_name = resolve_mol_name(mol_id, dft_opt_dir)
    if mol_name is None:
        raise ValueError(f"mol_{mol_id}: no directory found under {dft_opt_dir}")

    mol_dir = dft_opt_dir / mol_name
    opt_gjf = mol_dir / f"{mol_name}_opt.gjf"
    opt_log = mol_dir / f"{mol_name}_opt.log"
    ci, ni, oi, _ = oxime_atom_map_from_gjf(opt_gjf)

    atoms = final_geometry(opt_log)
    return {
        "name": mol_name,
        "atoms": atoms,
        "atom_spec": pyscf_atom_spec(atoms),
        "charge": CHARGE,
        "spin": SPIN,
        "ci": ci, "ni": ni, "oi": oi,
        "basis_note": "our own pipeline's 6-311+G(d,p), all-electron -- exact match to config.py",
    }


def atoms_before(lines: list[str], idx: int) -> list[tuple]:
    """Atoms from the last 'Standard orientation' block before line idx --
    the same anchor beckmann.dft.parse_nbo.r_no_before uses to tag a CMO/
    E2PERT table with its R(N-O), but returning the full geometry instead of
    just the N-O distance."""
    so_blocks = [(i, atoms) for i, atoms in parse_standard_orientations(lines) if i < idx]
    if not so_blocks:
        return []
    _, atoms = max(so_blocks, key=lambda x: x[0])
    return atoms


def _stage_points_from_log(log_path: Path, ni: int, oi: int) -> dict[float, list]:
    """{r_no: atoms} for every CMO section in one stage log, keyed by the
    N-O distance of the geometry it was computed on. Last table at a given R
    wins (Stable=Opt prints a pre-optimization seed pass and a separate
    post-optimization pass at the same frozen scan-point R -- see
    beckmann.dft.parse_cmo.parse_log, which applies the identical rule to
    the NBO data itself)."""
    lines = log_path.read_text().splitlines()
    atoms_by_r: dict[float, list] = {}
    for start in find_cmo_sections(lines):
        atoms = atoms_before(lines, start)
        if not atoms:
            continue
        r_no = round(no_distance(atoms, ni, oi), 4)
        atoms_by_r[r_no] = atoms
    return atoms_by_r


def _mol_stage_points(mol: str, mol_dir: Path, ni: int, oi: int) -> dict[str, list]:
    """{stage_label: atoms} for the 'nbo' (R0) stage plus every 'scan_N'
    rigid-scan point of a normal (non-STEP_SCAN_SOURCES) molecule."""
    points: dict[str, list] = {}
    nbo_log = mol_dir / f"{mol}_nbo.log"
    if nbo_log.exists():
        by_r = _stage_points_from_log(nbo_log, ni, oi)
        if by_r:
            (atoms,) = by_r.values()  # exactly one CMO table in _nbo.log
            points["nbo"] = atoms

    scan_log = mol_dir / f"{mol}_scan.log"
    if scan_log.exists():
        by_r = _stage_points_from_log(scan_log, ni, oi)
        for point, r_no in enumerate(sorted(by_r), start=1):
            points[f"scan_{point}"] = by_r[r_no]
    return points


def _mol_stage_points_stepscan(mol: str, mol_dir: Path, ni: int, oi: int) -> dict[str, list]:
    """Same role as _mol_stage_points, but merges scan geometries from one or
    more dft_opt_stepscan/ reruns (STEP_SCAN_SOURCES) instead of the
    canonical, crashed _scan.log -- mirrors
    beckmann.dft.parse_cmo.collect_molecule_stepscan, applied to geometries
    instead of NBO data, so R(N-O) point identity matches the trusted series."""
    points: dict[str, list] = {}
    nbo_log = mol_dir / f"{mol}_nbo.log"
    (atoms,) = _stage_points_from_log(nbo_log, ni, oi).values()
    points["nbo"] = atoms

    all_by_r: dict[float, list] = {}
    for source in STEP_SCAN_SOURCES[mol]:
        source_dir = step_scan_dir() / source
        _, s_ni, s_oi, _ = oxime_atom_map_from_gjf(source_dir / f"{source}_opt.gjf")
        all_by_r.update(_stage_points_from_log(source_dir / f"{source}_scan.log", s_ni, s_oi))

    for point, r_no in enumerate(sorted(all_by_r), start=1):
        points[f"scan_{point}"] = all_by_r[r_no]
    return points


def load_test_set_scan_series(mol_id: str) -> list[dict]:
    """One 'case' dict (same shape as load_case()/load_test_set_case(), plus
    'stage'/'r_no') per available R(N-O) point -- 'nbo' followed by every
    'scan_N' -- for a main-pipeline test-set molecule. Lets the open-source
    method be run across a full scan instead of just the equilibrium
    geometry, using the same geometry-anchor convention and
    STEP_SCAN_SOURCES merging the main pipeline's own descriptor extraction
    (beckmann.dft.parse_cmo) uses, so point identity/R(N-O) values line up
    with the trusted NBO7 series.
    """
    if mol_id not in TEST_IDS:
        raise ValueError(f"{mol_id}: not in TEST_IDS ({sorted(TEST_IDS)})")
    dft_opt_dir = DATA_OUTPUT / "dft_opt"
    mol_name = resolve_mol_name(mol_id, dft_opt_dir)
    if mol_name is None:
        raise ValueError(f"mol_{mol_id}: no directory found under {dft_opt_dir}")

    mol_dir = dft_opt_dir / mol_name
    ci, ni, oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol_name}_opt.gjf")

    if mol_name in STEP_SCAN_SOURCES:
        points = _mol_stage_points_stepscan(mol_name, mol_dir, ni, oi)
    else:
        points = _mol_stage_points(mol_name, mol_dir, ni, oi)

    ordered_stages = ["nbo"] + sorted(
        (s for s in points if s.startswith("scan_")), key=lambda s: int(s.split("_")[1])
    )
    cases = []
    for stage in ordered_stages:
        if stage not in points:
            continue
        atoms = points[stage]
        cases.append({
            "name": mol_name,
            "stage": stage,
            "r_no": round(no_distance(atoms, ni, oi), 4),
            "atoms": atoms,
            "atom_spec": pyscf_atom_spec(atoms),
            "charge": CHARGE,
            "spin": SPIN,
            "ci": ci, "ni": ni, "oi": oi,
            "basis_note": "our own pipeline's 6-311+G(d,p), all-electron -- exact match to config.py",
        })
    return cases
