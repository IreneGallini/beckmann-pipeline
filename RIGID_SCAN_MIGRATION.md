# `rigid-scan-architecture` branch — merge notes

Written for a fresh session picking this up to merge into `main`. Self-contained
— don't assume access to the conversation that produced this branch.

## Why this branch exists

The supervisor reviewed `beckmann/dft/inputs.py`'s `_scan_gjf()` (Stage 3, the
N-O relaxed scan) and objected to its architecture: it's a single Gaussian job
using Gaussian's native multi-point scan (`opt=(ModRedundant) ... B {ni} {oi}
S 4 0.1`). Verified empirically: this genuinely relaxes the geometry at all 5
points (5x "Optimization completed" per log), but the external NBO7 population
analysis only actually fires at 2 of those 5 points (R0 and R0+0.4 — always
exactly 2 "NBO 7.0" banners, never 5, across every molecule checked). The
existing workaround (`beckmann/dft/scan.py::extract_scan_sp()`) reruns NBO on
the missing 3 points' already-converged geometries as bolted-on single-point
jobs scraped from the completed log. Numerically valid, but not what's
actually in the `.gjf` file — the supervisor's real objection.

She provided a reference file, `oxime_001_scan.gjf` (repo root, tracked on
this branch), showing her intended methodology: per target R, a **rigid
displacement** of just the leaving-group atoms from a fixed base geometry
(pure geometry, no SCF), then an **independent constrained optimization**
(bond frozen at the new length, everything else relaxes fresh from that rigid
guess — points do NOT chain from each other), then a same-checkpoint
`--Link1--` continuation into a full NBO7/CMO single point. This is `main`'s
architecture on this branch, replacing Stage 3 only.

## Code changes (all currently on the branch or staged, see "Commit status" below)

- **`beckmann/dft/geometry.py`** (new file) — `parse_standard_orientations()`,
  `no_distance()`, `find_leaving_group()`, `displace_leaving_group()`. Split
  out of `scan.py` into their own module so `inputs.py` can use
  `displace_leaving_group()` without a circular import (`scan.py` imports
  `TEST_IDS`/`resolve_mol_name` from `inputs.py`). `scan.py` re-imports these
  names from `geometry.py` for backward compatibility with existing callers
  (`parse_nbo.py`, `parse_cmo.py`, `descriptors.py` all still
  `from beckmann.dft.scan import ...` unchanged).
- **`beckmann/dft/inputs.py`**:
  - New `_scan_gjf_rigid(name, base_atoms, ni, oi, oxime_label, basis=BASIS,
    step=0.1, n_points=4)` — builds the Link1-chained architecture described
    above. `step`/`n_points` are parameterized (not hardcoded) so the same
    function covers both the standard 4-point/0.1 Å series and one-off finer
    scans (see mol_006 investigation below). NBO keywords are `CMO PRINT=2
    E2PERT=0.05 BNDIDX` — **deliberately omits** her literal `NBOMO=P120`
    print-window restriction, because `parse_cmo.py` was separately fixed
    earlier (see its own docstring) to search the *entire* virtual manifold
    with no MO-count/energy cutoff, after discovering real target antibonds
    mixing in above a narrower window for some substrates. A fixed narrow
    `NBOMO` range would reintroduce that bug at the Gaussian-printing level.
    Also adds `SCF=(Tight,XQC) NoSymm` to match her robustness settings.
  - `_opt_gjf()`, `_nbo_gjf()`, `_scan_gjf_rigid()` all gained an optional
    `basis: str = BASIS` parameter (defaults preserve existing behavior for
    every existing caller) — added to support the mol_006 basis-set
    investigation (see below) without touching the global `BASIS` constant.
  - **The old `_scan_gjf()` (native multi-point scan) is untouched** and
    still exists — nothing currently calls `_scan_gjf_rigid()` from
    `prepare_opt()`/`main_opt()`. This branch does NOT change default
    pipeline behavior; the new architecture was only used for the manual
    comparison run described below. **Deciding whether `prepare_opt()` should
    switch to the rigid architecture by default is an open decision, not yet
    made or implemented.**
- **`beckmann/dft/hpc.py`** (bug fix, unrelated to the architecture question
  but found and fixed while running these experiments): every `cmd_submit_*`
  function (opt/nbo/scan/scan-sp/sp/ts/irc) and `cmd_download` built its
  remote-side directory-matching from a bare glob
  (`mol_{id}_*` or `*`) run directly on the shared remote directory —
  **not scoped to the local `--dir` actually selected**. Since every local
  experiment directory (`dft_opt/`, `dft_opt_rigidscan/`, `dft_opt_631g/`,
  `dft_opt_finescan/`, ...) uploads into the *same* flat remote directory,
  this meant `--mol 006 submit-opt` from one experiment's `--dir` could (and
  twice did, during this session) also resubmit/overwrite a different
  experiment's already-completed job for the same numeric ID. Fixed via a new
  `_remote_dir_names(local_dir, mol)` helper that resolves exact directory
  names from what's actually present locally (reusing the existing
  `mol_dirs()`), and builds the remote loop (`for name in "dir1" "dir2"; do`)
  from that explicit list instead of a glob. Applies to all 7 submit
  commands. **This fix has no default-behavior risk** — it only narrows what
  gets matched, verified via `--dry-run` against both a `--mol`-scoped and a
  bare (no `--mol`) call before trusting it.

  `cmd_download` needed a **second, separate fix** — the first pass used
  `rsync --filter=+ {name}/` per selected directory, but paired it with a
  blanket `--include=*/` (needed so rsync would recurse into the selected
  directory's contents). That combination doesn't work: rsync include
  patterns without a leading `/` match at **any** depth, so `--include=*/`
  matched *every* top-level directory on the remote, not just the intended
  one — completely undermining the `--filter=+` scoping. This silently
  re-downloaded unrelated experiments' logs into the wrong local directory
  (confirmed happening: a `--dir data/output/dft_opt_rigidscan --mol 002
  download` pulled in `mol_006_E`, `mol_014_Z`, `mol_020_E`, `mol_021_E`,
  `mol_029_Z`, even `mol_006_E_631g` and `mol_006_E_finescan`, alongside the
  intended `mol_002_E_rigidscan`). Fixed by dropping the blanket
  `--include=*/` entirely and anchoring each directory's rules to the
  transfer root instead: `--include=/{name}/` + `--include=/{name}/*.log`
  per selected directory, then a single trailing `--exclude=*`. Verified via
  `--dry-run` that this matches only the intended directory/directories.
  **If you see stray non-experiment-matching directories appear in a
  `dft_opt_*` local folder after a download, this is the bug to check for
  regressions in** — confirm the current `cmd_download` still uses anchored
  `/name/` includes, not a bare `--include=*/`.

## Validation: does the new architecture change the science?

Ran the new architecture for the original 4 test molecules
(mol_002_E/006_E/020_E/021_E) in a parallel local+remote directory
(`data/output/dft_opt_rigidscan/`, remote dirs suffixed `_rigidscan` to avoid
colliding with the originals) via a throwaway comparison script,
`scripts/analysis/compare_rigidscan.py` (tracked, reuses
`parse_nbo.parse_log()`/`parse_cmo.parse_log()` directly rather than
duplicating their table-finding logic — no changes needed there, both already
handle "however many tables are in a log" generically).

**Result: wCNmax (the project's central descriptor) is essentially
unchanged** — differences ~0.000–0.004 at every point across all 4 molecules.
w78max similarly stable. w17max/Ψ show small-to-moderate systematic
differences (Ψ ~5–10% lower under the new architecture at normal points).

**The important finding**: mol_020's R0+0.4 point — which needed an ad-hoc
`sp5` crash-recovery restart earlier this session under the *old*
architecture (built by extracting/shifting a neighboring point's geometry) —
gave Ψ≈0 there, a sharp anomalous collapse. Under the new independent
rigid-scan architecture, the same point gives Ψ=1.35, smoothly continuing the
trend from the previous points. This is fairly direct evidence the old
ad-hoc restart geometry was subtly wrong, and the rigid-scan architecture's
independent, well-defined-from-R0 starting point avoided that failure mode.
This is the strongest concrete argument for adopting the new architecture,
not just a reproducibility/methods-section cleanup.

**mol_006_E_rigidscan's point 4 (R0+0.4) crashed** under the new architecture
too (same failure shape as mol_020 originally hit — ran out its step budget,
Force oscillating, then segfaulted). Excluded rather than fixed at the time
(user's call) — mol_006 only has 4 of 5 points in that comparison run.

## Unresolved: mol_006 doesn't show the wCNmax minimum the reference paper reports

Extensive investigation this session, all ruled out as the cause (see
`Notes.md` for the "open issue" section this connects to, and `JOB_ISSUES.md`
for the mol_020 crash precedent this pattern-matches against):

- **Not the scan architecture** — wCNmax trend is ~identical old vs. new for
  the points both have.
- **Not substrate misidentification** — mol_006_E rigorously confirmed
  (graph/ring-position analysis, not eyeballing) to be the same compound as
  the supervisor's reference log `5_s0_Me.log` (repo root, tracked): methyl
  group exactly 3 ring-bonds from the aryl-fusion carbon in both, vs. 2 and 1
  bonds for the other two methyl-indanone positional isomers in the benchmark
  (mol_009, mol_013).
- **Not the aryl/alkyl channel assignment** — verified correct via direct
  connectivity inspection (`get_substituent_map()` output cross-checked
  against the actual bonded atoms).
- **Not the LUMO-to-LUMO+0.4 a.u. energy window** — built
  `scripts/analysis/compare_wcnmax_window.py` (tracked) to test whether
  restricting the search to the paper's window (vs. this codebase's
  unrestricted full-virtual-manifold search) would reveal the minimum.
  Result: windowed and unrestricted searches give byte-identical results at
  every point for mol_006 (and, as a sanity check, for mol_002/mol_020 too) —
  the winning MO always already sits inside the window. Ruled out.
- **Not the basis set** — initially a strong lead (our R0 geometry gives
  wCNmax=0.4225 vs. 0.457 on the reference log's geometry; the reference log
  turned out to use a much smaller basis, ~120 basis functions vs. our
  336, consistent with `6-31G(d)`) — **but the supervisor confirmed directly
  that her `oxime_001_scan.gjf` file was only a basis/method sensitivity
  test, and results for this compound type don't depend on it.** Ruled out
  by her explicit statement, not further computation.

**Currently in progress at the time of writing**: a rerun of mol_006 with the
*same* `wB97XD/6-311+G(d,p)` basis but **8 scan points at 0.05 Å resolution**
instead of 4 points at 0.1 Å (same R0 to R0+0.4 Å range, double the density),
testing whether a narrow minimum is being stepped over by the coarser grid.
Directory: `data/output/dft_opt_finescan/mol_006_E_finescan/` (untracked,
uncommitted — see below). Built via the new `step`/`n_points` parameters on
`_scan_gjf_rigid()`. **Status as of this writing: 4 of 8 points complete,
job still running on Citadel** (`~/beckmann/dft_opt_test/mol_006_E_finescan/`)
— whoever picks this up should check `hpc_sync.py --dir
data/output/dft_opt_finescan status`, download once all 8 points show
"Optimization completed" and "NBO 7.0" banners, then parse and check the
wCNmax(R) trend for an interior minimum (reuse `compare_rigidscan.py`'s
`new_architecture_series()` pattern, or write a similar throwaway script
pointed at this directory).

## Commit status — action needed before merging

```
$ git log --oneline main..rigid-scan-architecture
b9ce92e test set ran again with rigid scan
faa5092 rigid replacement per R + optimization

$ git status --short
 M beckmann/dft/hpc.py       <- the remote-scoping bug fix, UNCOMMITTED
 M beckmann/dft/inputs.py    <- basis/step/n_points params, UNCOMMITTED
?? data/output/dft_opt_631g/     <- abandoned basis experiment, untracked
?? data/output/dft_opt_finescan/ <- in-progress fine-resolution rerun, untracked
```

**Before merging:**
1. Commit the `hpc.py`/`inputs.py` changes described above (two logically
   separate changes — the hpc.py bug fix could be its own commit, independent
   of the rigid-scan work, since it affects the whole pipeline not just this
   branch's experiment).
2. `data/output/dft_opt_631g/` — the abandoned basis-set experiment. Its
   Stage 1 finished (confirmed basis doesn't matter per the supervisor) but
   nothing downstream was run. Safe to delete; not needed for anything.
3. `data/output/dft_opt_finescan/` — decide once the mol_006 investigation
   concludes whether this is worth keeping/committing (the `.gjf` is small
   and useful to keep as a record; `.log` files are gitignored regardless).
4. **The core open decision this branch doesn't resolve**: should
   `prepare_opt()`/`main_opt()` switch to `_scan_gjf_rigid()` as the default
   for all future molecules, superseding `_scan_gjf()` and
   `extract_scan_sp()`/`submit-scan-sp`/`STAGES` in `parse_nbo.py`? The
   validation above supports it (wCNmax stable, fixed a real mol_020 error),
   but this wasn't executed — only the 4-molecule comparison run was done in
   a side directory. Making it the default requires: updating
   `prepare_opt()` to call `_scan_gjf_rigid()`, deciding whether to backfill
   mol_014_Z/mol_029_Z (untouched by this branch) and whether to regenerate
   mol_006/020's already-crashed-point-4 data under the new architecture, and
   updating `CLAUDE.md`'s pipeline documentation accordingly (currently still
   describes only the old `_scan_gjf()` architecture).

---

## MERGED (2026-07-16): rigid-scan architecture is now the default

All four open items above are resolved:

- `prepare_opt()` now writes only `_opt.gjf`/`_nbo.gjf` (Stages 1-2) upfront
  — it can't generate Stage 3 in the same call anymore, since
  `_scan_gjf_rigid()` needs Stage 1's *converged* geometry, which doesn't
  exist until Stage 1 has actually run and its log is downloaded. New
  function `prepare_scan_rigid(mol_dir, name, basis=BASIS, step=0.1,
  n_points=4)` generates Stage 3 as a separate step after that — call it
  once `{name}_opt.log` is on disk, then upload/`submit-scan`/download as
  before. `_scan_gjf()` (old architecture) stays in `inputs.py` for
  reference/rollback, just unused by default now.
- `STAGES` in `parse_nbo.py`/`parse_cmo.py` simplified to `["nbo", "scan"]`
  — the old `"sp2"/"sp3"/"sp4"/"sp5"` extracted-single-point workaround no
  longer applies to any current test molecule. `extract_scan_sp.py` and
  `hpc.py`'s `submit-scan-sp` are kept (harmless no-ops without matching
  `_sp*.gjf` files) but marked legacy in their docstrings/help text.
- `descriptors.py`'s `resolve_series()` no longer uses a fixed
  `SERIES_STAGES` list (which assumed exactly 4 stretched points via the old
  `sp2/sp3/sp4/scan_2` naming) — it now dynamically takes `"nbo"` plus every
  `"scan_N"` stage present, sorted numerically. This is what makes mol_006_E's
  9-point finescan series and the other 5 molecules' 5-point series both work
  through the same code path with no special-casing. `SERIES_FALLBACK` (the
  mol_020_E `sp5` patch) is gone — no longer needed now that a crashed point
  doesn't take the old internal-walk job's *other* points down with it (each
  rigid-scan point is independent by construction).
- mol_014_Z and mol_029_Z backfilled with the rigid-scan architecture
  (generated directly under their canonical names this time, not a
  side-experiment suffix — no promotion step needed). mol_002_E/mol_020_E/
  mol_021_E's already-completed `_rigidscan` data and mol_006_E's `_finescan`
  data were promoted into the canonical `data/output/dft_opt/mol_XXX_E/`
  locations (reused as-is, not re-run — their `.gjf`/`.log` internal
  `%chk`/title text still reads `..._rigidscan`/`..._finescan`, cosmetic
  only).

**A real bug was caught during this regeneration**, worth knowing about if
you're touching `parse_nbo.py`/`parse_cmo.py` again: every rigid-scan NBO
block uses `Stable=Opt` (a wavefunction stability re-check), which — exactly
as already documented for the supervisor's own reference log in `Notes.md`'s
Task-5 section — prints **two** full NBO/CMO tables at the same geometry (a
pre-optimization seed pass and the real post-optimization one), not one. The
old architecture never hit this because its two tables (R0 and R0+0.4) were
always at genuinely different R values, so the latent bug in both files'
point-disambiguation logic (which enumerated raw table order rather than
deduping by R first) never triggered. Under the rigid-scan architecture it
did, immediately and silently — first regeneration attempt produced 9
"points" for a 4-point molecule (mol_020_E/mol_021_E) with every R value
doubled, and doubled (summed) E2PERT rows feeding into Ψ. Fixed at the
source in both files' `parse_log()`: when multiple tables share the same R,
keep only the *last* one. Verified by re-running and confirming mol_020_E/
mol_021_E returned to exactly 5 points (not 9) and mol_006_E's real interior
minimum was still intact at exactly one point (not duplicated) afterward.
If you ever see a molecule's point count roughly double what you expect,
this is the first thing to check.

Old architecture's descriptor CSVs/plots and raw scan data (for the 4
molecules replaced) are archived at
`data/output/analysis/archive_pre_rigidscan_2026-07-15/` and
`data/output/dft_opt/_archive_pre_rigidscan/` respectively — not deleted.

**RESOLVED (2026-07-16)**: mol_014_Z and mol_029_Z's rigid-scan Citadel jobs
completed overnight (both 8/8 Normal termination, 0 errors) and were
downloaded, replacing their old-architecture `_sp2/3/4.gjf`/`.log` (archived
per above) with the new `_scan.log`. Full descriptor set (`parse_nbo.py` →
`parse_cmo.py` → `descriptors.py` → `summarize_descriptors.py` →
`parse_wiberg.py` → `plot_bond_orders.py`) regenerated for all 6 molecules.
All six now resolve through the same `resolve_series()` path: mol_002_E/
014_Z/020_E/021_E/029_Z each give a 5-point series (`nbo` + `scan_1..4`),
mol_006_E gives its 9-point finescan series (`nbo` + `scan_1..8`).

Notably, mol_014_Z now shows an interior wCNmax extremum (R=1.7034, depth
0.0586) at standard 0.1 Å resolution that wasn't visible in its stale
3-point data — despite being an F (fragmentation) outcome, unlike the
paper's pattern where the extremum was reported only for the rearranging
reference compound. Worth flagging at the supervisor meeting alongside the
0.05 Å resolution question, not something to over-interpret without her
input first.
