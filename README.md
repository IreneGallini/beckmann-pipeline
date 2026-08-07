# Beckmann Rearrangement Pipeline

Computational pipeline for predicting whether an oxime substrate undergoes Beckmann rearrangement or C–C fragmentation. Starts from a ketone SMILES, generates 3D conformers, optimises geometry with AIMNet2, and predicts the outcome via a wCNmax-minimum rule, computed either through Gaussian/NBO7 (validated against a 34-molecule benchmark) or through PySCF (open-source, no HPC required).

## Repository structure

This is a monorepo of three installable packages plus a `research/` directory for exploratory work:

```
beckmann-pipeline/
├── packages/
│   ├── beckmann-core/     shared library (oxime conversion, conformers, AIMNet2
│   │                      optimization, geometry primitives, the wCNmax-minimum
│   │                      R/F rule), method-agnostic, no filesystem conventions
│   ├── beckmann-nbo/      Gaussian/NBO7/Citadel product
│   └── beckmann-pyscf/    open-source, HPC-free product (Flask web app)
├── data/                  shared canonical benchmark data
├── research/              exploratory/investigation code (depends on the
│                          packages above, never the other way around)
└── environment.yml        single dev environment covering all three packages
```

See `CLAUDE.md` for the full architecture writeup, and `research/README.md` for what's in `research/` and how to run it.

## Environment setup

```bash
conda env create -f environment.yml
conda activate beckmann
```

Installs all three packages in editable mode. Python 3.11. `KMP_DUPLICATE_LIB_OK=TRUE` is set automatically by any module that imports Auto3D/AIMNet2 (required on macOS).

## Quick start

**Open-source, HPC-free prediction (no Gaussian/Citadel needed):**

```bash
cd packages/beckmann-pyscf/backend
python app.py   # http://localhost:5001
```

Or from Python directly:

```python
from beckmann_pyscf.pipeline import predict
result = predict("O=C1CCC2=CC=CC=C21")  # alpha-tetralone
# result["prediction"] is "R" or "F"
```

**Gaussian/NBO7 benchmark pipeline** (34-molecule set, requires Citadel access):

```bash
python research/benchmark_pipeline/00_benchmark_to_oximes.py
python research/benchmark_pipeline/01_smiles_to_conformers.py
python research/benchmark_pipeline/02_select_and_optimize.py
python packages/beckmann-nbo/scripts/prepare_opt.py
python packages/beckmann-nbo/scripts/hpc_sync.py --mol 002 upload
# ... see CLAUDE.md for the full Citadel workflow
```

## Tests

```bash
python -m pytest packages/beckmann-core/tests/ -m "not slow"
python -m pytest packages/beckmann-nbo/tests/
python -m pytest packages/beckmann-pyscf/tests/ packages/beckmann-pyscf/backend/tests/ -m "not slow" -c packages/beckmann-pyscf/backend/pytest.ini
python -m pytest research/tests/
```

## Baseline result

The classical anti-periplanar dihedral rule gives **20/34 correct predictions (59%)** on the benchmark set. The wCNmax-minimum rule gives **25/34 (74%)** via NBO7, which motivated both the DFT/NBO phase and, later, the PySCF work to make that rule available without a Gaussian/HPC dependency. See `CLAUDE.md` for the full writeup.
