"""
Tests for beckmann/dft/descriptors.py (aryl/alkyl role tagging, Psi, d/dR slopes).
"""
import pytest

from beckmann.dft.descriptors import get_substituent_map, least_squares_slope


@pytest.fixture(scope="module")
def dft_opt_dir(project_root):
    p = project_root / "data" / "output" / "dft_opt"
    if not p.exists():
        pytest.skip("dft_opt/ not found — run scripts/dft/prepare_opt.py first")
    return p


def test_substituent_map_mol_002(dft_opt_dir, best_per_substrate_sdf_path):
    """Known-good case: mol_002_E's aryl carbon is C6, alkyl carbon is C10."""
    subst = get_substituent_map("mol_002_E", dft_opt_dir / "mol_002_E")
    assert subst["ci"] == 11
    assert subst["ni"] == 12
    assert subst["oi"] == 13
    assert subst["c_aryl"] == 6
    assert subst["c_alkyl"] == 10


def test_substituent_map_mismatch_raises(dft_opt_dir, tmp_path):
    """A .gjf label that disagrees with RDKit's own atom map must raise, not be trusted silently."""
    bad_dir = tmp_path / "mol_002_E"
    bad_dir.mkdir()
    (bad_dir / "mol_002_E_opt.gjf").write_text("mol_002_E opt  [oxime: C1=N2-O3]\n")
    with pytest.raises(ValueError, match="atom map mismatch"):
        get_substituent_map("mol_002_E", bad_dir)


def test_least_squares_slope_matches_known_line():
    xs = [1.0, 2.0, 3.0, 4.0, 5.0]
    ys = [2.0, 4.0, 6.0, 8.0, 10.0]
    assert least_squares_slope(xs, ys) == pytest.approx(2.0)


def test_least_squares_slope_none_with_one_point():
    assert least_squares_slope([1.0], [2.0]) is None
