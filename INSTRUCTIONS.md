# Handoff: how to use the pipeline

Written at the end of my internship as a practical "how do I actually run this"
guide, organized by the three ways to use the pipeline: the **web app**, the
**`beckmann-nbo` CLI**, and the **standalone scripts**. For architecture and
why things are structured this way, see `CLAUDE.md` (kept up to date, will
outlive this file).

## 0. One-time setup

```bash
conda env create -f environment.yml
conda activate beckmann
```

This installs all three packages (`beckmann-core`, `beckmann-nbo`,
`beckmann-pyscf`) in editable mode.

If you'll use the Gaussian/NBO7/Citadel path, copy `.env.example` to `.env`
at the repo root and fill in with SSH key.

## 1. Web app 

```bash
cd packages/beckmann-pyscf/backend
python app.py
# open http://localhost:5001
```

This is also callable directly from Python without the web UI, if you just
want a result in a script or notebook:

```python
from beckmann_pyscf.pipeline import predict
result = predict("O=C1CCC2=CC=CC=C21")  # example SMILES
result["prediction"]
```

## 2. `beckmann-nbo` CLI 

It needs Citadel (or any server with Gaussian16 + NBO7).

```bash
beckmann-nbo init                                   # writes .env SSH settings
beckmann-nbo verify                                 # SSH reachable? g16 executable? NBO7 wrapper set up?
beckmann-nbo predict --smiles "O=C1CCC2=C1C=CC=C2" --name test1   # SMILES -> conformers -> AIMNet2 opt -> submits Stage 1+2 to Citadel
# wait for Stage 1 to finish on the cluster 
beckmann-nbo predict --continue test1 --dir <workdir>/dft_opt     # generates + submits Stage 3 (the N-O bond scan)
beckmann-nbo status --mol test1 --dir <workdir>/dft_opt           # per-stage job status; once Stage 3 is done, a live R/F prediction
beckmann-nbo report --mol test1 --out <dir> --advanced            # wCNmax/bond-order plots + classical-vs-wCNmax comparison
```

Stage 1 = geometry optimization: Gaussian DFT optimization (wB97XD/6-311+G(d,p)) of the AIMNet2-optimized starting geometry, producing {name}_opt.gjf/.log. 

Stage 2 = NBO7 single-point: an NBO7 single-point calculation (pop=nbo7read) run on the Stage 1 converged geometry, producing {name}_nbo.gjf/.log. Result: E2PERT, BNDIDX, NBOSUM, and CMO data at the equilibrium geometry.

Stage 3 = rigid N–O bond scan: stretching the N–O bond away from the Stage 1 equilibrium (rigid displacement → constrained re-optimization → NBO7 single point, at each point), producing {name}_scan.gjf/.log. This generates the wCNmax-vs-R(N–O) series.

- `predict --csv path.csv` (with `id`/`SMILES` columns, same shape as
  `data/input/benchmark.csv`) submits Stage 1+2 for a whole batch of
  molecules at once.

## 3. Standalone scripts 

The CLI above is the recommended path for *new* molecules. The benchmark
set itself was processed with the underlying scripts directly.

**Stage 0: SMILES/AIMNet2 (no HPC):**

```bash
python research/benchmark_pipeline/00_benchmark_to_oximes.py   # benchmark.csv -> molecules.smi + benchmark_meta.json
python research/benchmark_pipeline/01_smiles_to_conformers.py  # 3D conformer generation (Auto3D)
python research/benchmark_pipeline/02_select_and_optimize.py   # AIMNet2 geometry optimization
```
Output lands in `data/output/conformers/` and
`data/output/aimnet_optimized/`.

**Stage 1+2: DFT/NBO7 input generation + Citadel submission:**

```bash
python packages/beckmann-nbo/scripts/prepare_opt.py    # two-stage opt+NBO gjf files for the test set
python packages/beckmann-nbo/scripts/hpc_sync.py --mol 002 upload
python packages/beckmann-nbo/scripts/hpc_sync.py --mol 002 submit-opt
python packages/beckmann-nbo/scripts/hpc_sync.py status
python packages/beckmann-nbo/scripts/hpc_sync.py --mol 002 download
```

**Stage 3: the N-O bond rigid scan** (needs Stage 1's converged geometry,
so it's generated as a separate step once that `.log` is downloaded):

```python
from beckmann_nbo.inputs import prepare_scan_rigid
from beckmann_nbo.config import DATA_OUTPUT
prepare_scan_rigid(DATA_OUTPUT / "dft_opt" / "mol_002_E", "mol_002_E")
```
then `hpc_sync.py --mol 002 upload` / `submit-scan` / `download` again.

**If a job crashes on Citadel** check
`data/output/dft_opt/JOB_ISSUES.md` first. About 1 in 5 benchmark
molecules hit a recurring non-convergence crash; there's an automated
recovery ladder for it:

```bash
python packages/beckmann-nbo/scripts/auto_recover.py            # all molecules
python packages/beckmann-nbo/scripts/auto_recover.py --mol 020  # one molecule
python packages/beckmann-nbo/scripts/auto_recover.py --dry-run  # preview only
```

**After `.log` files are downloaded, extract descriptors:**

```bash
python packages/beckmann-nbo/scripts/parse_nbo.py      # -> data/output/analysis/nbo_e2pert.csv
python packages/beckmann-nbo/scripts/parse_cmo.py       # -> data/output/analysis/cmo_descriptors.csv (Lambda, wCNmax)
python packages/beckmann-nbo/scripts/descriptors.py     # -> data/output/analysis/channel_descriptors.csv + descriptor_slopes.csv
python packages/beckmann-nbo/scripts/parse_wiberg.py    # -> data/output/analysis/bond_order_scan.csv
```

**Single molecule wCNmax plot**: wCNmax-vs-R(N-O) chart from already-parsed data (i.e. `parse_cmo.py`/
`descriptors.py` have already been run so `channel_descriptors.csv` and
`cmo_channel_extraction.csv` exist in `data/output/analysis/`).
Use `research/analysis_scripts/plot_single_wcnmax.py`:

```bash
PYTHONPATH=research python research/analysis_scripts/plot_single_wcnmax.py mol_002_E
# writes mol_002_E_wcnmax.png in the current directory; --out to choose a path
```

## Other documentation

- `packages/beckmann-nbo/README.md`, `packages/beckmann-pyscf/README.md`:
  per-package usage details, deployment notes.
- `research/README.md`: what's in `research/` and how to run it.
- `research/Notes.md`: debugging/investigation narrative (e.g. why some
  molecules needed a finer scan step size).
