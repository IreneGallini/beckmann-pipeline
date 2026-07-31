# research/

Exploratory and investigation code, kept separate from the two shipped
products (`packages/beckmann-nbo/`, `packages/beckmann-pyscf/`). Nothing
here is a dependency of either product; everything here depends on one or
both of them plus `beckmann-core`.

## Running scripts here

This directory is not an installable package. Put it on `PYTHONPATH` (or
`sys.path`) before running anything in it, so cross-directory imports like
`from pyscf_validation.geometry import ...` or `from viz import ...`
resolve:

```bash
PYTHONPATH=research python research/pyscf_validation/compare_wcnmax.py
```

## Contents

- `benchmark_pipeline/` -- the original scripts 00-02 that convert
  `data/input/benchmark.csv` into conformers and AIMNet2-optimized
  structures for all 34 molecules. Built on `beckmann_core`'s reusable
  functions; this is "how the benchmark set was actually processed," not
  part of either product's own API.
- `pyscf_validation/` -- the investigation comparing PySCF's wCNmax (the
  same engine `beckmann-pyscf` ships) against the trusted NBO7 numbers
  across the benchmark set. `geometry.py`'s loaders all read pre-existing
  Gaussian `.log` files, which is why this is validation tooling and not
  part of the product itself (see `beckmann_pyscf.engine.pair_nbo`, which
  this package imports from). `run_cases.py` holds the benchmark-harness
  runners (`run_case`/`run_test_set_case`/`run_test_set_scan_series`) that
  used to live inside `pair_nbo.py` itself before this restructuring.
- `diabatic_character/` -- the diabatic character-exchange descriptor
  work (unvalidated on PySCF; explicitly out of scope for either product).
- `analysis_scripts/` -- one-off comparison/plotting/validation scripts not
  covered above, plus `classical_benchmark.py`/`wcnmax_rule_benchmark.py`
  (the benchmark-CSV-driving harnesses around `beckmann_core.classical`'s
  and `beckmann_core.wcnmax_rule`'s own live product functions).
- `ts_ml/` -- transition-state location exploration (Gaussian QST2/QST3,
  AIMNet2-based product-geometry construction, and a NEB/PySisyphus ML
  proxy). Parked/exploratory -- no current callers in either product.
- `viz.py` -- generic wCNmax(R) plotting helpers used across the scripts
  above.
- `example_scans/` -- Tetiana's external reference logs (`5_s*_Me.log`)
  used throughout the PySCF validation work.
- `reference/` -- misc. reference artifacts (`oxime_001_scan.gjf`).
- `stray_logs/` -- accidental working-directory debris from local
  Gaussian/PySisyphus runs (`calculator.log`, `cos.log`, etc.) -- kept for
  now, safe to delete.
- `archived_ketone_pipeline/` -- the original (pre-benchmark) ChemDraw-based
  ketone-to-oxime pipeline, already marked archived before this
  restructuring.
- `data/` -- archived/side-experiment Gaussian output not part of the
  canonical 34-molecule benchmark set (basis-set sensitivity checks,
  alternate scan resolutions, etc.). The canonical data
  (`dft_opt/`, `dft_sp/`, `aimnet_optimized/`, `conformers/`, the
  non-archived parts of `analysis/`) stays in the shared top-level `data/`,
  since `beckmann-nbo`'s own scripts read from there too.
- `tests/` -- integration tests for `benchmark_pipeline/`'s generated
  artifacts and the `ts_ml/` exploration; not unit tests for
  `beckmann-core`'s own functions (see `packages/beckmann-core/tests/` for
  those).
