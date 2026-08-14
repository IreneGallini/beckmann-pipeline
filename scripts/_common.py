"""Shared helpers used by every numbered script in this directory. Not a
CLI, not an entry point on its own -- just the one bit of bookkeeping
(query-id sanitizing + default paths) that has to stay byte-identical
across all ten scripts, since beckmann_nbo.hpc's upload/submit/status/
download functions glob directories as "mol_{id.zfill(3)}_*" and would
silently find nothing if a script computed the id differently.
"""
from pathlib import Path

from rdkit import Chem

EXPORT_ROOT = Path(__file__).resolve().parent.parent
QUERY_PREFIX = "mol"  # required by beckmann_nbo.hpc.mol_dirs()'s "mol_*" glob -- not configurable


def sanitize_id(name: str) -> str:
    """'test1' -> 'qtest1' (mirrors beckmann_nbo.cli_predict._sanitize_id).
    The 'q' prefix keeps a query id from ever colliding with a 3-digit
    benchmark id (001-034, not present in this export anyway); underscores
    are collapsed because prepare_opt()'s test_ids filter does
    name.split('_')[1] on the 3-token 'mol_{id}_{E|Z}' convention."""
    return ("q" + name.replace("_", "-")).zfill(3)


def workdir_for(mol_id: str) -> Path:
    return EXPORT_ROOT / "data" / "output" / "query_predictions" / mol_id


def load_query_mol(mol_name: str, workdir: Path) -> Chem.Mol:
    """Load one molecule's optimized geometry from THIS query's own
    aimnet_optimized/best_per_substrate.sdf -- unlike
    beckmann_nbo.descriptors.get_substituent_map()/_load_mols(), which are
    hardcoded to the benchmark set's global data/output/aimnet_optimized/
    best_per_substrate.sdf and never see a query molecule's geometry at all
    (a known limitation of the beckmann-nbo CLI's `report` command for
    fresh molecules)."""
    sdf_path = workdir / "aimnet_optimized" / "best_per_substrate.sdf"
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    for mol in suppl:
        if mol is not None and mol.GetProp("_Name") == mol_name:
            return mol
    raise ValueError(f"{mol_name}: not found in {sdf_path}")


def local_substituent_map(mol_name: str, mol_dir: Path, workdir: Path) -> dict:
    """{ci, ni, oi, c_aryl, c_alkyl} (1-based) for one query molecule,
    cross-checked against its own {mol_name}_opt.gjf label -- the
    query-local equivalent of beckmann_nbo.descriptors.get_substituent_map()."""
    from beckmann_core.classical import get_oxime_atoms
    from beckmann_nbo.scan import oxime_atom_map_from_gjf

    rdkit_mol = load_query_mol(mol_name, workdir)
    result = get_oxime_atoms(rdkit_mol)
    if result is None:
        raise ValueError(f"{mol_name}: oxime substructure or aryl/alkyl neighbors not found")
    cox, nox, oox, c_aryl, c_allyl = (idx + 1 for idx in result)

    gjf_ci, gjf_ni, gjf_oi, _ = oxime_atom_map_from_gjf(mol_dir / f"{mol_name}_opt.gjf")
    if (cox, nox, oox) != (gjf_ci, gjf_ni, gjf_oi):
        raise ValueError(
            f"{mol_name}: atom map mismatch -- RDKit gives C{cox}=N{nox}-O{oox}, "
            f".gjf label gives C{gjf_ci}=N{gjf_ni}-O{gjf_oi}"
        )
    return {"ci": cox, "ni": nox, "oi": oox, "c_aryl": c_aryl, "c_alkyl": c_allyl}
