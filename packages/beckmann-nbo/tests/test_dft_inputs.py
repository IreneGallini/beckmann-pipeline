"""
Tests for scripts/dft/prepare_sp.py and scripts/dft/prepare_opt.py output.

Validates:
  data/output/dft_sp/   — single-point NBO inputs (all 34 molecules)
  data/output/dft_opt/  — two-stage opt+NBO inputs (test set: 002, 006, 020, 021)
"""
import re
import pytest
from rdkit import Chem

OXIME_PAT = Chem.MolFromSmarts("[C:1]=[N:2]-[O+:3]")
LABEL_RE = re.compile(r"\[oxime:\s*C(\d+)=N(\d+)-O(\d+)\]")
TEST_SET = {"002", "006", "020", "021"}


def _charge_mult(gjf_text: str) -> tuple[int | None, int | None]:
    """Parse (charge, multiplicity) from the charge/multiplicity line."""
    for line in gjf_text.splitlines():
        parts = line.split()
        if len(parts) == 2 and all(p.lstrip("-").isdigit() for p in parts):
            return int(parts[0]), int(parts[1])
    return None, None


# ── dft_sp/ fixtures and tests ─────────────────────────────────────────────────

@pytest.fixture(scope="module")
def dft_sp_dir(data_output):
    p = data_output / "dft_sp"
    if not p.exists() or not list(p.glob("**/*.gjf")):
        pytest.skip("dft_sp/ not found — run scripts/dft/prepare_sp.py first")
    return p


def test_sp_gjf_files_exist(dft_sp_dir):
    assert len(list(dft_sp_dir.glob("**/*.gjf"))) > 0


def test_sp_charge_plus1(dft_sp_dir):
    for gjf in dft_sp_dir.glob("**/*.gjf"):
        charge, mult = _charge_mult(gjf.read_text())
        assert charge == 1, f"{gjf.name}: charge={charge}, expected 1 (protonated oxime)"
        assert mult == 1,   f"{gjf.name}: multiplicity={mult}, expected 1"


def test_sp_route_is_single_point_with_nbo(dft_sp_dir):
    for gjf in dft_sp_dir.glob("**/*.gjf"):
        text = gjf.read_text()
        assert "sp" in text,           f"{gjf.name}: missing 'sp' in route"
        assert "pop=nbo7read" in text, f"{gjf.name}: missing 'pop=nbo7read' (needed for CMO)"


def test_sp_nbo_block_present(dft_sp_dir):
    for gjf in dft_sp_dir.glob("**/*.gjf"):
        text = gjf.read_text()
        assert "$NBO" in text and "$END" in text, f"{gjf.name}: missing $NBO...$END block"
        assert "E2PERT" in text, f"{gjf.name}: missing E2PERT keyword in $NBO block"


def test_sp_oxime_label_matches_sdf(dft_sp_dir, best_per_substrate_sdf_path):
    """[oxime: C=N-O] label in each .gjf must match the RDKit substructure match."""
    sdf_mols: dict[str, Chem.Mol] = {}
    for mol in Chem.SDMolSupplier(str(best_per_substrate_sdf_path), removeHs=False):
        if mol is not None:
            sdf_mols[mol.GetProp("_Name")] = mol

    for gjf in dft_sp_dir.glob("**/*.gjf"):
        name = gjf.stem
        if name not in sdf_mols:
            continue
        mol = sdf_mols[name]
        match = mol.GetSubstructMatch(OXIME_PAT)
        assert match, f"{name}: oxime pattern not found in SDF molecule"
        expected = tuple(idx + 1 for idx in match)

        m = LABEL_RE.search(gjf.read_text())
        assert m, f"{gjf.name}: no [oxime: C=N-O] label in title line"
        found = (int(m.group(1)), int(m.group(2)), int(m.group(3)))
        assert found == expected, (
            f"{name}: .gjf label C{found[0]}=N{found[1]}-O{found[2]} "
            f"!= expected C{expected[0]}=N{expected[1]}-O{expected[2]}"
        )


# ── dft_opt/ fixtures and tests ────────────────────────────────────────────────

@pytest.fixture(scope="module")
def dft_opt_dir(data_output):
    p = data_output / "dft_opt"
    if not p.exists() or not list(p.glob("**/*.gjf")):
        pytest.skip("dft_opt/ not found — run scripts/dft/prepare_opt.py first")
    return p


def test_opt_test_set_has_both_stages(dft_opt_dir):
    """Each test-set substrate must have exactly one mol directory with _opt.gjf and _nbo.gjf."""
    for mol_id in TEST_SET:
        dirs = list(dft_opt_dir.glob(f"mol_{mol_id}_*/"))
        assert len(dirs) == 1, (
            f"mol_{mol_id}: expected 1 directory (lowest-energy isomer only), found {len(dirs)}"
        )
        mol_dir = dirs[0]
        name = mol_dir.name
        assert (mol_dir / f"{name}_opt.gjf").exists(), f"Missing: {name}_opt.gjf"
        assert (mol_dir / f"{name}_nbo.gjf").exists(), f"Missing: {name}_nbo.gjf"


def test_opt_gjf_route_has_opt_not_nbo(dft_opt_dir):
    """Stage 1 files should request geometry optimisation but NOT NBO (that's Stage 2)."""
    for gjf in dft_opt_dir.glob("**/*_opt.gjf"):
        text = gjf.read_text().lower()
        assert "opt" in text, f"{gjf.name}: missing 'opt' keyword in route"
        assert "pop=nboread" not in text and "pop=nbo7read" not in text, (
            f"{gjf.name}: Stage 1 should not have pop=nbo(7)read — NBO belongs in Stage 2"
        )


def test_opt_gjf_charge(dft_opt_dir):
    for gjf in dft_opt_dir.glob("**/*_opt.gjf"):
        charge, mult = _charge_mult(gjf.read_text())
        assert charge == 1, f"{gjf.name}: charge={charge}, expected 1"
        assert mult == 1,   f"{gjf.name}: multiplicity={mult}, expected 1"


def test_nbo_gjf_reads_from_checkpoint(dft_opt_dir):
    """Stage 2 must specify geom=checkpoint and %oldchk to read the Stage 1 geometry."""
    for gjf in dft_opt_dir.glob("**/*_nbo.gjf"):
        text = gjf.read_text()
        assert "geom=checkpoint" in text, f"{gjf.name}: missing 'geom=checkpoint'"
        assert "%oldchk=" in text,        f"{gjf.name}: missing '%oldchk' line"


def test_nbo_gjf_nbo_block(dft_opt_dir):
    for gjf in dft_opt_dir.glob("**/*_nbo.gjf"):
        text = gjf.read_text()
        assert "$NBO" in text and "$END" in text, f"{gjf.name}: missing $NBO...$END block"
        assert "E2PERT" in text,      f"{gjf.name}: missing E2PERT keyword"
        assert "CMO" in text,         f"{gjf.name}: missing CMO keyword (needs NBO7)"
        assert "pop=nbo7read" in text, f"{gjf.name}: missing 'pop=nbo7read'"


def test_nbo_gjf_charge(dft_opt_dir):
    for gjf in dft_opt_dir.glob("**/*_nbo.gjf"):
        charge, mult = _charge_mult(gjf.read_text())
        assert charge == 1, f"{gjf.name}: charge={charge}, expected 1"
        assert mult == 1,   f"{gjf.name}: multiplicity={mult}, expected 1"


def test_scan_gjf_exists(dft_opt_dir):
    """Each test-set mol directory must have a _scan.gjf for the N-O bond scan."""
    for mol_id in TEST_SET:
        dirs = list(dft_opt_dir.glob(f"mol_{mol_id}_*/"))
        assert len(dirs) == 1, f"mol_{mol_id}: expected 1 directory, found {len(dirs)}"
        mol_dir = dirs[0]
        name = mol_dir.name
        assert (mol_dir / f"{name}_scan.gjf").exists(), f"Missing: {name}_scan.gjf"


def test_scan_gjf_has_ModRedundant_and_E2PERT(dft_opt_dir):
    """Scan file must request ModRedundant opt, have E2PERT+CMO in the NBO block, and a B scan line.

    pop=nbo7read (not pop=nboread) routes through Gaussian's external-program
    interface to NBO7, required for CMO-based descriptors (Lambda, wCNmax).
    """
    for gjf in dft_opt_dir.glob("**/*_scan.gjf"):
        text = gjf.read_text()
        assert "ModRedundant" in text,          f"{gjf.name}: missing 'ModRedundant' in route"
        assert "pop=nbo7read" in text.lower(),  f"{gjf.name}: missing 'pop=nbo7read'"
        assert "geom=checkpoint" in text, f"{gjf.name}: missing 'geom=checkpoint'"
        assert "%oldchk=" in text,        f"{gjf.name}: missing '%oldchk' line"
        assert "E2PERT" in text,          f"{gjf.name}: missing 'E2PERT' in NBO block"
        assert "CMO" in text,             f"{gjf.name}: missing 'CMO' in NBO block"
        lines = text.splitlines()
        scan_lines = [l for l in lines if l.startswith("B ") and " S " in l]
        assert scan_lines, f"{gjf.name}: missing 'B N O S 4 0.1' ModRedundant scan line"