# Beckmann reaction summer project
### Structure
beckmann-pipeline/
├── app.py           Flask server
├── templates/
│   └── index.html   UI
├── scripts/
│   └── 01_smiles_to_conformers.py
└── data/
    ├── inputs/
    └── outputs/

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
pip install aimnet2calc
```

### Issues
UserWarning: pkg_resources is deprecated as an API. See https://setuptools.pypa.io/en/latest/pkg_resources.html. The pkg_resources package is slated for removal as early as 2025-11-30. Refrain from using this package or pin to Setuptools<81.

Should be replaced with importlib.metadata