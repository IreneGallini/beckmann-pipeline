# beckmann-nbo

Gaussian/NBO7 product for the Beckmann rearrangement pipeline: DFT input generation, HPC submission, NBO7 log parsing, and NBO-specific descriptors (Psi, Lambda). Depends on `beckmann-core` (pinned version) for the shared oxime/conformer/AIMNet2/geometry primitives and the wCNmax-minimum R/F rule. Not tied to Citadel specifically — any SSH-reachable server running Gaussian16 + NBO7 works.

```bash
pip install -e .   # pulls in beckmann-core automatically; installs the beckmann-nbo CLI
python -m pytest tests/ -v
```

## `beckmann-nbo` CLI

The recommended entry point for a new user — wraps this package's underlying modules so you can go from a SMILES string to a prediction without reading the full workflow below. Every subcommand is a thin wrapper: no pipeline logic lives in the CLI layer itself.

```bash
beckmann-nbo init                                    # write .env with your cluster's SSH/Gaussian/NBO7 settings
beckmann-nbo verify                                   # preflight: SSH reachable? G16_PATH executable? NBO7 wrapper set up?
beckmann-nbo predict --smiles "O=C1CCC2=C1C=CC=C2" --name test1   # SMILES -> conformers -> AIMNet2 opt -> submit Stage 1+2
beckmann-nbo predict --continue qtest1 --dir <workdir>/dft_opt    # once Stage 1 finishes, generate + submit Stage 3 (the N-O scan)
beckmann-nbo status --mol qtest1 --dir <workdir>/dft_opt          # per-stage status + live R/F prediction once Stage 3 is clean
beckmann-nbo recover --mol 020                        # run one pass of the automated oscillation-recovery ladder
beckmann-nbo report --mol 020 --out <dir> --advanced   # wCNmax/bond-order/E2PERT plots + classical-vs-wCNmax comparison
```

`--dry-run`, `--mol`, and `--dir` are global flags — pass them **before** the subcommand (e.g. `beckmann-nbo --dry-run --mol 002 status`), matching `hpc.py`'s own convention.

**`verify`'s NBO7 check is worth understanding**: Gaussian doesn't run NBO7 directly — it shells out to a helper script (`gaunbo7`/`gaunbo6`) which then calls the real NBO7 binary. On a shared server the vendor-installed copy is usually root-owned and not executable by your account. If that's not fixed, `pop=nbo7read` silently falls back to the old bundled NBO 3.1: the job still reaches "Normal termination," but the log has no CMO section, so wCNmax can't be computed — and nothing tells you this happened until you grep the log yourself. `verify`'s third check catches this *before* you submit anything, and on failure prints the exact fix (copy `gaunbo7`/`gaunbo6` into a directory you own, `chmod +x`, point `NBO_WRAPPER_DIR` at it) parameterized to your own `.env`, not hardcoded to Citadel.

**Multi-molecule prediction**: `predict --csv path.csv` accepts a file with `id`/`SMILES` columns (same shape as `data/input/benchmark.csv`) and submits Stage 1+2 for every row.

**Limitations**: `predict`/`status`/`report` for a fresh query molecule (not one of the 34 benchmark substrates) only take you through Stage 1+2 automatically and Stage 3 via `--continue`; `report`'s classical-baseline comparison currently only resolves molecules already present in `data/output/aimnet_optimized/best_per_substrate.sdf` (i.e. the benchmark set) — a query molecule's own `report` output will show `classical=unavailable` until that's extended.

For the full DFT/HPC workflow this CLI wraps (Citadel-specific setup notes, output file reference, NBO7 install details, the benchmark-batch scripts under `scripts/`), see the monorepo root `CLAUDE.md`.

Config (`.env`, cluster connection details) lives at the monorepo root — copy `.env.example` there and fill in values, or run `beckmann-nbo init`.
