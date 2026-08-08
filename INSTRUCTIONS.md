# Handoff: how to use the pipeline

four ways to use the pipeline:
- **`beckmann-pyscf` CLI**
- **`beckmann-nbo` CLI**
- **web app**
- **standalone scripts**. 

## 0. One-time setup

```bash
conda env create -f environment.yml
conda activate beckmann
```

This installs all three packages (`beckmann-core`, `beckmann-nbo`,
`beckmann-pyscf`) in editable mode.

If you'll use the Gaussian/NBO7/Citadel path, copy `.env.example` to `.env`
at the repo root and fill in with SSH key.

## 1. `beckmann-pyscf` CLI (AIMNet2 + PySCF, no HPC)

```bash
beckmann-pyscf predict --smiles "O=C1CCC2=C1C=CC=C2" --name test1 --plot
# [1/4]..[4/4] progress -> ./beckmann_pyscf_runs/test1/{optimized.sdf,wcnmax_series.csv,summary.txt,wcnmax_vs_rno.png}
```

That single command runs the whole pipeline: conformers → AIMNet2
optimization → PySCF wCNmax scan → R/F prediction, printing progress at
each stage and writing its results into `./beckmann_pyscf_runs/test1/` by
default (`--out` to choose a different directory).

**Batch mode, a whole CSV of molecules at once:**

```bash
beckmann-pyscf predict --csv molecules.csv --plot
```

`molecules.csv` needs `id`/`SMILES` columns, the same shape as
`data/input/benchmark.csv` 

**Stage-by-stage, to inspect each step's output on its own** (useful for
checking whether something looks wrong before trusting a full `predict` run):

```bash
beckmann-pyscf conformers --smiles "O=C1CCC2=C1C=CC=C2" --name test1
# -> prints the conformers SDF path

beckmann-pyscf optimize --conformers-sdf <path from above>
# -> prints AIMNet2 energy + the resolved oxime atom map (C/N/O/aryl/alkyl indices),
#    writes best.sdf. This is the checkpoint that catches a bad atom-map resolution
#    before it silently breaks the scan stage downstream

beckmann-pyscf scan --sdf <path from above> --plot
# -> runs the 7-point PySCF wCNmax scan, prints the R/F prediction,
#    writes wcnmax_series.csv/summary.txt/wcnmax_vs_rno.png
```

Each subcommand takes an explicit input path (copy whatever the previous
command printed) rather than any implicit state tracking, so any stage can
also be run standalone against an SDF you already have from somewhere
else: `optimize` doesn't require its input to have come from `conformers`,
and `scan` doesn't require its input to have come from `optimize`.

`scan` also takes `--ci`/`--ni`/`--oi`/`--c-aryl`/`--c-alkyl` (1-based atom
indices) to override auto-detection if `get_oxime_atoms()` picks the wrong
atom on an unusual substrate; `optimize`'s printed atom map gives you the
exact values to pass. And `--r-min`/`--r-max`/`--r-step` to adjust the scan
window/resolution.

**Callable directly from Python too**, if you want a result in a script or
notebook rather than the CLI:

```python
from beckmann_pyscf.pipeline import predict

if __name__ == "__main__":
    result = predict("O=C1CCC2=CC=CC=C21")  # example SMILES
    print(result["prediction"])
```

## 2. `beckmann-nbo` CLI

It needs Citadel (or any server with Gaussian16 + NBO7).

```bash
beckmann-nbo init                                   # writes .env SSH settings
beckmann-nbo verify                                 # SSH reachable? g16 executable? NBO7 wrapper set up?
beckmann-nbo predict --smiles "O=C1CCC2=C1C=CC=C2" --name test1   # SMILES -> conformers -> AIMNet2 opt -> submits Stage 1+2 to Citadel
```

**Read the directory path `predict` prints before doing anything else.**
Every command after this one needs a `--dir` pointing at the same job
folder, and the easiest way to get that right is to copy the path `predict`
itself prints rather than reconstruct it from memory.

### How the directories actually work

`--name test1` does **not** become the folder name verbatim. It's run
through an internal sanitizer (`_sanitize_id`) that prefixes a `q` and
strips underscores, so `--name test1` becomes the id `qtest1`. That id is
then used consistently everywhere: as the folder name, and as the value
you pass to `--mol`/`--continue` in every later command. You don't need to
compute this yourself: just re-use whatever `--name` you originally typed
(the CLI re-derives the same `qtest1` from it each time), or copy the
folder name straight out of the printed path.

Without `--workdir`, a fresh `predict --smiles` call creates:

```
data/output/query_predictions/qtest1/          <- workdir for this molecule
├── conformers/
├── aimnet_optimized/
└── dft_opt/                                    <- the "--dir" every later command needs
    └── qtest1_E/  (and/or qtest1_Z/)
        ├── qtest1_E_opt.gjf / .log             <- Stage 1
        ├── qtest1_E_nbo.gjf / .log              <- Stage 2
        └── qtest1_E_scan.gjf / .log             <- Stage 3 (after --continue)
```

After the first `predict` call finishes, it prints exactly this:

```
Submitted. Poll with:
  beckmann-nbo status --mol qtest1 --dir data/output/query_predictions/qtest1/dft_opt
Once Stage 1 shows Normal termination, continue to Stage 3 with:
  beckmann-nbo predict --continue qtest1 --dir data/output/query_predictions/qtest1/dft_opt
```

```bash
# wait for Stage 1 to finish on the cluster, then:
beckmann-nbo predict --continue qtest1 --dir data/output/query_predictions/qtest1/dft_opt   # generates + submits Stage 3 (the N-O bond scan)
beckmann-nbo status --mol qtest1 --dir data/output/query_predictions/qtest1/dft_opt          # per-stage job status; once Stage 3 is done, a live R/F prediction
beckmann-nbo report --mol qtest1 --dir data/output/query_predictions/qtest1/dft_opt --out /tmp/qtest1_report --advanced   # wCNmax/bond-order plots + classical-vs-wCNmax comparison
```

Note `report`'s `--out` where the generated PNG plots and `classical_vs_wcnmax.txt` get
written. It defaults to nothing (you must pass it).

Stage 1 = geometry optimization: Gaussian DFT optimization (wB97XD/6-311+G(d,p)) of the AIMNet2-optimized starting geometry, producing {name}_opt.gjf/.log.

Stage 2 = NBO7 single point: an NBO7 single point calculation (pop=nbo7read) run on the Stage 1 converged geometry, producing {name}_nbo.gjf/.log. Result: E2PERT, BNDIDX, NBOSUM, and CMO data at the equilibrium geometry.

Stage 3 = rigid N–O bond scan: stretching the N–O bond away from the Stage 1 equilibrium (rigid displacement → constrained re-optimization → NBO7 single point, at each point), producing {name}_scan.gjf/.log. This generates the wCNmax-vs-R(N–O) series.

- `predict --csv path.csv` (with `id`/`SMILES` columns, same shape as
  `data/input/benchmark.csv`) submits Stage 1+2 for a whole batch of
  molecules at once, each row gets its own `data/output/query_predictions/<sanitized-id>/`
  folder, following the same layout as above.

## 3. Web app (prototype, for external collaborators)

```bash
cd packages/beckmann-pyscf/backend
python app.py
# open http://localhost:5001
```

This is a hosted prototype of the same AIMNet2+PySCF pipeline as Section 1,
meant for external collaborators who don't have this repo or conda
environment set up locally, not the way to run this yourself. If you have
the repo cloned, use the `beckmann-pyscf` CLI (Section 1) instead: it's
faster to invoke, gives per-stage progress and stage-level inspection, and
doesn't require running a local Flask server.

## 4. Standalone scripts

These are the literal scripts used to build the 34-molecule benchmark set,
not general-purpose tools written to accept an arbitrary new molecule.
Before reusing any of them, check the two things below, both of which trip
people up if skipped.

### "Standalone" doesn't mean self-contained

Every script here still `import`s from `beckmann_core`/`beckmann_nbo`, so
downloading a single `.py` file on its own will not run: it needs the full
repo checked out and the `beckmann` conda env active with all three
packages installed (`conda env create -f environment.yml`, see Section 0).
Scripts under `research/` additionally need `PYTHONPATH=research` set for
their own local imports.

### Most of these scripts are hardwired to the benchmark set, not to "whatever `--dir` you point them at"

This is the part that matters most for testing a **new** molecule. Most
`packages/beckmann-nbo/scripts/*.py` files take **no arguments at all**:
they call a `main()` that hardcodes both the working directory
(`DATA_OUTPUT / "dft_opt"`, i.e. `data/output/dft_opt/`) and the set of
molecule IDs to loop over. A new molecule generated
outside that directory, under a different ID, is invisible to these
scripts: they'll just print `-- mol_XXX: no directory, skipping` for
every benchmark ID and never look at your molecule at all.

**If the molecule isn't one of the 34 benchmark substrates, use the
`beckmann-nbo` CLI (Section 2) instead of these scripts.** The CLI's
`predict`/`status`/`report` commands call the same underlying functions
these scripts call (`prepare_opt`/`prepare_scan_rigid`, `collect_molecule`,
`find_wcnmax_minimum`, etc.) but with an explicit `--mol`/`--dir` per
invocation rather than a hardcoded benchmark scope: that's the whole
reason the CLI exists as a separate layer on top of this package instead
of duplicating it.

The rest of this section documents what these scripts actually do, useful
for understanding/debugging the benchmark set itself or for adapting one
in place (e.g. temporarily editing `ALL_IDS`/`dft_opt_dir` in a copy of a
script) rather than for running unmodified against a new molecule.

**Stage 0: SMILES/AIMNet2 (no HPC):**

```bash
python research/benchmark_pipeline/00_benchmark_to_oximes.py   # benchmark.csv -> molecules.smi + benchmark_meta.json
python research/benchmark_pipeline/01_smiles_to_conformers.py  # 3D conformer generation (Auto3D)
python research/benchmark_pipeline/02_select_and_optimize.py   # AIMNet2 geometry optimization
```
Output lands in `data/output/conformers/` and
`data/output/aimnet_optimized/`.

Each of these three scripts' `main()` hardcodes
benchmark paths (`data/input/benchmark.csv`, `data/output/conformers/`,
etc.), but the actual work happens in a genuinely reusable, molecule-agnostic
`beckmann_core` function underneath (`generate_conformers()`,
`select_and_optimize()`) that takes explicit path arguments and knows
nothing about the benchmark set. For a single new molecule you don't need
these scripts at all: that's exactly what the web app and the
`beckmann-nbo predict`/beckmann-pyscf `predict()` calls in Sections 1-2 do
by calling those same `beckmann_core` functions directly. Reach for these
scripts only if you're regenerating conformers/geometries for the whole
benchmark set (or scripting something similar over a CSV of your own,
calling `generate_conformers()`/`select_and_optimize()` yourself with your
own paths).

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

