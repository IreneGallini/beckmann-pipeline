# beckmann-hpc-pipeline (export)

## Structure

```
exports/beckmann-hpc/
├── packages/
│   ├── beckmann-core/   oxime conversion, conformer generation, AIMNet2
│   │                    optimization, the wCNmax-minimum R/F rule
│   └── beckmann-nbo/    Gaussian/NBO7 DFT input generation, SSH job
│                        submission, NBO7 log parsing (library modules only)
├── scripts/              the numbered, run them one at a time pipeline
│                        (see "Running a molecule" below)
├── data/                empty on checkout; populated per-molecule as you run
│                        the scripts (gitignored see .gitignore)
├── environment.yml
└── .env.example         your own cluster's SSH/Gaussian/NBO7 paths
```

## One-time setup

```bash
cd exports/beckmann-hpc
conda env create -f environment.yml
conda activate beckmann

cp .env.example .env   # fill in your own cluster's SSH host, G16_PATH, etc.
```

`.env` needs, at minimum:
- `HPC_HOST` -- SSH destination (username@hostname), ideally an alias from `~/.ssh/config`
- `HPC_REMOTE_DIR` -- working directory on your cluster
- `G16_PATH` -- full path to the `g16` binary on your cluster
- `NBOEXE` / `NBO_WRAPPER_DIR` -- NBO7 interface binary + an *executable*
  wrapper directory for `gaunbo7`/`gaunbo6` (vendor installs are often
  root-owned and not executable -- see the "NBO7 setup" note below)

SSH key auth must already be set up (`ssh-copy-id`) so commands run without
password prompts.

## Running a molecule

Each script in `scripts/` has a `MOL_NAME`/`SMILES` edit it for your molecule, then run the script directly (`python scripts/00_smiles_to_conformers.py` from this
directory, or `cd scripts && python 00_smiles_to_conformers.py`). Every
downstream script re-derives the same `qtest1`-style job id from `MOL_NAME`
(see `scripts/_common.py`), so keep `MOL_NAME` consistent across all ten
scripts for one molecule.

Run in order:

| # | Script | What it does |
|---|---|---|
| 00 | `00_smiles_to_conformers.py` | SMILES -> protonated oxime isomers (E/Z) -> Auto3D conformers |
| 01 | `01_optimize_aimnet2.py` | lowest-energy conformer per isomer -> AIMNet2/ASE geometry optimization |
| 02 | `02_prepare_stage12_inputs.py` | writes Stage 1 (`_opt.gjf`) + Stage 2 (`_nbo.gjf`) Gaussian input |
| 03 | `03_upload_submit_stage12.py` | uploads to your cluster over SSH, submits Stage 1 |
| 04 | `04_check_status.py` | per-stage job status + local crash diagnosis; reusable, rerun anytime |
| 05 | `05_download_results.py` | pulls whatever `.log` files exist; reusable, rerun anytime |
| 06 | `06_prepare_stage3_scan.py` | (after Stage 1 converges) writes Stage 3 (`_scan.gjf`, the N-O rigid scan) |
| 07 | `07_upload_submit_stage3.py` | uploads + submits Stage 3 (then reuse 04/05 to poll/pull it) |
| 08 | `08_parse_descriptors.py` | (after Stage 3 downloads) parses `.log` files into E2PERT/CMO descriptor CSVs |
| 09 | `09_predict_rf.py` | the payoff: wCNmax R/F prediction + classical-baseline comparison + plot |
| 10 | `10_recover_if_crashed.py` | optional if 04 flags an oscillation crash |

This produces, per molecule:

```
data/output/query_predictions/qtest1/          <- workdir for this molecule
├── conformers/
├── aimnet_optimized/
├── analysis/                                    <- CSVs + wcnmax_vs_rno.png from 08/09
└── dft_opt/
    └── mol_qtest1_E/  (and/or mol_qtest1_Z/)
        ├── mol_qtest1_E_opt.gjf / .log          <- Stage 1: DFT geometry optimization
        ├── mol_qtest1_E_nbo.gjf / .log          <- Stage 2: NBO7 single point
        └── mol_qtest1_E_scan.gjf / .log         <- Stage 3: N-O bond rigid scan
```

## Stages, in brief

- **Stage 1** -- DFT geometry optimization (wB97XD/6-311+G(d,p)) of the
  AIMNet2-optimized starting geometry.
- **Stage 2** -- NBO7 single point (`pop=nbo7read`) on the Stage 1 converged
  geometry: E2PERT, BNDIDX, NBOSUM, CMO.
- **Stage 3** -- rigid N-O bond scan (rigid displacement -> constrained
  re-optimization -> NBO7 single point, at each point), generating the
  wCNmax-vs-R(N-O) series the R/F prediction is read off of.

## NBO7 setup, a real trap worth knowing up front

```bash
mkdir -p ~/beckmann/nbo7_bin
cp /path/to/nbo7/bin/gaunbo7 /path/to/nbo7/bin/gaunbo6 ~/beckmann/nbo7_bin/
chmod +x ~/beckmann/nbo7_bin/gaunbo7 ~/beckmann/nbo7_bin/gaunbo6
```

then set `NBO_WRAPPER_DIR` to that directory in `.env`. `03_upload_submit_stage12.py`/
`04_check_status.py` will fail with a clear message pointing here if `.env`
is incomplete. To confirm NBO7 (not 3.1) actually ran on a finished job,
check the `.log` for `NBO 7.0` and a `CMO: NBO Analysis of Canonical
Molecular Orbitals` section.
