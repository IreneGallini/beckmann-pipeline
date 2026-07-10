"""
Tests for beckmann/dft/ts.py and beckmann/dft/ts_products.py output.

Pilot scope only: mol_002_E rearrangement TS (TS1_A1). Fragmentation-channel
and other-molecule files don't exist yet, so tests here are gated on the
pilot files specifically rather than the full TEST_SET.
"""
import re
import pytest
from rdkit import Chem

LABEL_RE = re.compile(r"\[oxime:\s*C(\d+)=N(\d+)-O(\d+)\]")


def _charge_mult_blocks(gjf_text: str) -> list[tuple[int, int]]:
    """All (charge, multiplicity) pairs in a possibly-multi-structure .gjf."""
    pairs = []
    for line in gjf_text.splitlines():
        parts = line.split()
        if len(parts) == 2 and all(p.lstrip("-").isdigit() for p in parts):
            pairs.append((int(parts[0]), int(parts[1])))
    return pairs


@pytest.fixture(scope="module")
def ts1_a1_gjf(project_root):
    p = project_root / "data" / "output" / "dft_opt" / "mol_002_E" / "mol_002_E_ts1_a1.gjf"
    if not p.exists():
        pytest.skip("mol_002_E_ts1_a1.gjf not found — run beckmann.dft.ts main() first")
    return p


@pytest.fixture(scope="module")
def product_rearr_sdf(project_root):
    p = project_root / "data" / "output" / "aimnet_optimized" / "mol_002_E_product_rearr.sdf"
    if not p.exists():
        pytest.skip("mol_002_E_product_rearr.sdf not found — run beckmann.dft.ts_products main() first")
    return p


# ── ts_products.py: product geometry ────────────────────────────────────────

def test_product_atom_count_matches_reactant(product_rearr_sdf, project_root):
    reactant_sdf = project_root / "data" / "output" / "aimnet_optimized" / "best_per_substrate.sdf"
    reactant = next(m for m in Chem.SDMolSupplier(str(reactant_sdf), removeHs=False)
                     if m.GetProp("_Name") == "mol_002_E")
    product = next(Chem.SDMolSupplier(str(product_rearr_sdf), removeHs=False))
    assert product.GetNumAtoms() == reactant.GetNumAtoms()
    r_syms = [a.GetSymbol() for a in reactant.GetAtoms()]
    p_syms = [a.GetSymbol() for a in product.GetAtoms()]
    assert r_syms == p_syms, "atom order must match the reactant exactly for QST2/NEB"


def test_product_total_charge_is_plus1(product_rearr_sdf):
    product = next(Chem.SDMolSupplier(str(product_rearr_sdf), removeHs=False))
    total_charge = sum(a.GetFormalCharge() for a in product.GetAtoms())
    assert total_charge == 1


def test_product_has_finite_coordinates(product_rearr_sdf):
    product = next(Chem.SDMolSupplier(str(product_rearr_sdf), removeHs=False))
    conf = product.GetConformer()
    import math
    for i in range(product.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        assert all(math.isfinite(v) for v in (p.x, p.y, p.z))


# ── ts.py: QST2 input ───────────────────────────────────────────────────────

def test_ts_route_is_qst2_with_calcfc_and_freq(ts1_a1_gjf):
    text = ts1_a1_gjf.read_text()
    assert "opt=(qst2,calcfc)" in text
    assert " freq " in text or text.rstrip().endswith("freq")
    assert "scrf=(smd,solvent=water)" in text


def test_ts_has_two_charge_mult_blocks(ts1_a1_gjf):
    pairs = _charge_mult_blocks(ts1_a1_gjf.read_text())
    assert len(pairs) == 2, "QST2 needs exactly 2 molecule specifications"
    assert all(p == (1, 1) for p in pairs)


def test_ts_nbo_block_present(ts1_a1_gjf):
    text = ts1_a1_gjf.read_text()
    assert "$NBO" in text and "$END" in text
    assert "pop=nbo7read" in text


def test_ts_two_structures_have_matching_atom_order(ts1_a1_gjf):
    lines = ts1_a1_gjf.read_text().splitlines()
    coord_line_re = re.compile(r"^([A-Za-z]{1,2})\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+$")
    blocks: list[list[str]] = []
    current: list[str] = []
    for line in lines:
        m = coord_line_re.match(line.strip())
        if m:
            current.append(m.group(1))
        elif current:
            blocks.append(current)
            current = []
    if current:
        blocks.append(current)
    assert len(blocks) == 2, f"expected 2 coordinate blocks, found {len(blocks)}"
    assert blocks[0] == blocks[1], "element order must match between the two QST2 structures"
