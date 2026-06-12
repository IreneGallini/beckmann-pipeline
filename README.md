# Beckmann reaction summer project

## Structure
```
beckmann-pipeline/
├── app.py           Flask server
├── templates/
│   └── index.html   UI
├── scripts/
│   └── 01_smiles_to_conformers.py
└── data/
    ├── inputs/
    └── outputs/
```
### /tests

Basic tests check that:
- imports don't generate an error
- script 1 and 2 generate output
- script 1 generates k conformers 
- script 2 generates lower energy structure + RMSD

### 01_smiles_to_conformers

Turn SMILES into top k conformers

### 02_select_and_optimize

Geometry optimization with AIMNet2 on lowest energy conformer

## Setting up python environment 

```bash
conda create -n beckmann python=3.10 
conda activate beckmann 
conda install -c conda-forge rdkit 
pip install torch --index-url https://download.pytorch.org/whl/cpu 
pip install Auto3D
conda env export > environment.yml
```
Anyone can run this command and have identical setup
```bash
conda env create -f environment.yml
```

Add:
```bash
pip install "aimnet[ase]"
pip install pytest
```

### Issues

UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
Should be replaced with importlib.metadata

**Warnings**
``` 
test_draft.py::test_imports
  /Users/irenegallini/miniconda3/envs/beckmann/lib/python3.11/site-packages/Auto3D/__init__.py:2: UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.
    from pkg_resources import get_distribution, DistributionNotFound

test_draft.py::test_imports
  /Users/irenegallini/miniconda3/envs/beckmann/lib/python3.11/site-packages/torch/jit/_script.py:1488: DeprecationWarning: `torch.jit.script` is deprecated. Please switch to `torch.compile` or `torch.export`.
    warnings.warn(
```