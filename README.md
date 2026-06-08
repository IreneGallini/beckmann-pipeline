# Beckmann reaction summer project
### Structure

## Setting up python environment 
```bash
conda create -n beckmann python=3.10 conda activate beckmann conda install -c conda-forge rdkit pip install torch --index-url https://download.pytorch.org/whl/cpu 
pip install Auto3D
conda env export > environment.yml
```
Anyone in team can run this command and have identical setup
```bash
conda env create -f environment.yml
```