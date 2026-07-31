"""
Thin wrapper: ketone SMILES -> protonated activated oxime E/Z .smi file ->
Auto3D/AIMNet2-prefiltered conformers. Both underlying calls are vendored
unchanged from beckmann-pipeline:

  - beckmann.oximes.ketone_to_protonated_oximes()/enumerate_ez() -- the
    per-molecule primitives beckmann.oximes.benchmark_to_oximes() loops over
    for the 34-molecule CSV; not vendoring benchmark_to_oximes() itself,
    it's batch-CSV-only. This is the same call sequence
    beckmann.run_prediction() already uses for a single new SMILES.
  - beckmann.conformers.generate_conformers() -- already generic over any
    .smi file path, no benchmark assumptions; called unmodified.
"""
from pathlib import Path

from rdkit import Chem

from beckmann_core.conformers import generate_conformers
from beckmann_core.oximes import enumerate_ez, ketone_to_protonated_oximes


def smiles_to_oxime_smi(smiles: str, workdir: Path, mol_name: str = "query") -> Path:
    """Ketone SMILES -> a .smi file with one line per E/Z oxime isomer, named
    '{mol_name}_E'/'{mol_name}_Z' (matches beckmann.run_prediction()'s own
    naming for a single new molecule). Raises ValueError if RDKit can't
    parse the SMILES, or if no ketone/oxime-convertible group is found."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit cannot parse SMILES: {smiles!r}")

    oximes = ketone_to_protonated_oximes(mol)
    if not oximes:
        raise ValueError("No oxime products found -- is this a ketone SMILES?")

    smi_path = workdir / f"{mol_name}.smi"
    with open(smi_path, "w") as f:
        for ox in oximes:
            for iso, ez_label in enumerate_ez(ox):
                f.write(f"{Chem.MolToSmiles(iso)} {mol_name}_{ez_label}\n")
    return smi_path


def run_conformers(smi_path: Path, output_dir: Path) -> Path:
    """Call-through to the vendored Auto3D conformer generator, unchanged."""
    return generate_conformers(smi_path, output_dir)
