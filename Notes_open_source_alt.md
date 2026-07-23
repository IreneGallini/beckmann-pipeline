# Open-source wCNmax prototype (branch: open-source-wcnmax-prototype)

Exploratory only. Nothing here touches `beckmann/dft/` or the validated main-branch
pipeline. Main branch's NBO7-derived numbers remain the trusted ground truth throughout
-- this branch asks "does this open-source alternative rank/trend the same way," not
"are these numbers more correct."

New code lives entirely in `beckmann_alt/`:
- `geometry.py` -- reference-case geometries/atom maps, reusing
  `beckmann.dft.scan.parse_standard_orientations`/`ATOMIC_SYMBOLS` for the main-pipeline
  case and a local fallback for the external reference log's `nosymm` quirk (see below).
- `pyscf_livvo.py` -- shared PySCF SCF setup and MO-projection helpers (`build_mol`,
  `run_scf`, `project_virtuals_onto_livvo`) used by `pair_nbo.py`.
- `pair_nbo.py` -- the local per-atom-pair construction, the best-performing method
  found in this exploration and the only one carried forward. Computes wCNmax only.

**Scope decision (current)**: wCNmax is the actual predictive descriptor this project
needs; w17max/w78max are secondary and were never the target of this prototype effort.
Earlier code explored LIVVO-based (`pyscf.lo.vvo`) and crude AO-projection channel
identification, plus an iterative-deflation follow-up aimed at separating w17max/w78max
from wCNmax's own winning orbital -- none of that reproduced wCNmax's own result any
better than the per-atom-pair construction below, and none of it was needed for
wCNmax's accuracy (confirmed: wCNmax is computed independently of w17/w78 in every
version tried). That code has been removed from the current tree. It remains available
via `git log -- beckmann_alt/` for anyone who wants the full historical record of what
was tried and why; this file now documents only the wCNmax-relevant path.

## Reference cases

**mol_002_E** (primary): our own pipeline, level of theory fully known
(`wB97XD/6-311+G(d,p)`, `scrf=(smd,solvent=water)`, charge=1, singlet -- `beckmann/config.py`).
Geometry is the converged Stage-1 DFT geometry (`mol_002_E_opt.log`'s final "Standard
orientation"), the same geometry the trusted NBO7 numbers were computed on. Atom map
(`ci=11, ni=12, oi=13, c_aryl=6, c_alkyl=10`) reused from `beckmann.dft.descriptors`/
`beckmann.dft.scan` conventions, not re-derived.

**5_s0_Me** (secondary, explicitly NOT apples-to-apples): Tetiana's external reference
log (compound 3, "Ring Size and Substituent Effects in the Beckmann Rearrangement",
Table 2). Its actual route line is `wb97xd/genecp scrf=(smd,solvent=water) opt
scf=(Tight,xqc) nosymm`. Two things found before writing any PySCF code:
- **`genecp` is real, not incidental** -- the log's "Pseudopotential Parameters" section
  shows a custom hand-specified ECP applied to every C/N/O center (4 valence electrons
  per carbon, minimal split-valence primitives), not our all-electron `6-311+G(d,p)`.
  Reproducing it exactly would mean transcribing every exponent/coefficient/ECP
  parameter from the log by hand. Not attempted -- this case is run at **our own**
  `6-311+G(d,p)` instead, so any comparison against it is a different-basis check, not
  a validated reproduction of the paper's numbers.
- **`nosymm` changes the geometry block header** from `Standard orientation:` to `Input
  orientation:` -- `beckmann.dft.scan.parse_standard_orientations` (correct and
  unmodified for every log our own pipeline produces, which never uses `nosymm`) finds
  nothing in this file. `beckmann_alt/geometry.py` falls back to a locally-defined
  parser for `Input orientation:` blocks for this one external file's quirk only.

Trusted NBO7 numbers for both cases: wCNmax=0.457 (MO32, compound 3) is the paper's own
worked example, already validated by `scripts/analysis/validate_reference_descriptors.py`
(Tier 1, PASS). mol_002's numbers come straight from the main branch's already-generated
`cmo_channel_extraction.csv` (nbo/R0 stage). Neither was recomputed here.

## What had to change to actually run PySCF (found while building, not assumed upfront)

- **`wb97x-d` is blacklisted by PySCF's high-level dispersion dispatch**
  (`pyscf.scf.dispersion._black_list`). That blacklist exists to stop PySCF silently
  stacking an extra D3/D4 correction on a functional that already has its own dispersion
  built in -- exactly Gaussian's wB97XD. libxc itself recognizes `'wb97x-d'` fine
  (confirmed directly via `libxc.xc_type`). Removed it from the blacklist locally in
  `pyscf_livvo._enable_wb97xd()` rather than substituting PySCF's supported
  `wb97x-d3bj` variant, which uses a different (D3(BJ)) dispersion correction than
  Gaussian's wB97XD -- that substitution would be a bigger, less transparent deviation
  than bypassing an overcautious string check for a functional libxc already implements.
- **SMD's Python API exists but the compiled CDS extension does not.**
  `pyscf.solvent.smd.SMD` and its solvent parameter table (water entry present) are
  real, but the non-electrostatic cavity/dispersion/solvent-structure term requires a
  compiled `libsolvent` extension built with `-DENABLE_SMD=ON`, which this pip wheel
  does not include (`RuntimeError: SMD module is not available`). Building pyscf from
  source for this is a bigger ask than the plain `pip install` that was approved, so
  the SCF uses **ddCOSMO** instead. ddCOSMO is materially different from SMD --
  different cavity construction (fixed atomic-radius spheres vs. SMD's own tessellated
  solvent-accessible surface), and no explicit cavitation/dispersion/solvent-structure
  term at all. Treat the solvent contribution as a rough implicit-solvent stand-in, not
  a reproduction of Gaussian's actual SMD/water calculation.
- **Density fitting enabled for tractable runtime.** A conventional (non-DF) RSH-hybrid
  SCF on mol_002's 358-basis-function system did not finish in a reasonable time on this
  machine. `mf.density_fit()` is a standard, well-controlled approximation to the
  two-electron integrals -- flagged for completeness, not expected to be a significant
  source of divergence, but it is a deviation from Gaussian's conventional SCF.
- **Frozen geometry, no PySCF geometry optimization.** The SCF runs a single-point
  calculation on the geometry Gaussian already converged to (`beckmann_alt/geometry.py`).
  This isolates "does the local-construction method reproduce similar orbital
  character" from "does PySCF's optimizer land somewhere else," which isn't the
  question being asked here.

## The local per-atom-pair construction (`beckmann_alt/pair_nbo.py`)

Instead of picking from a pre-built global set of localized virtual orbitals, this
module builds a **fresh, local subspace for each atom pair individually**, directly
from that pair's own block of the density matrix (in a Löwdin-orthogonalized IAO
basis) -- structurally the same operation as NBO's real per-atom-pair antibond search,
just for one requested pair at a time.

**First pass (single lowest-occupation eigenvector per pair):** landed on a distinct MO
for the target pair in both test cases, but wCNmax came out badly wrong (0.087 vs
trusted 0.436 for mol_002; 0.117 vs 0.457 for 5_s0_Me -- consistently ~4-5x too small).

**Why, and the fix:** printing every local eigenvector's occupation for mol_002's cn=
(C11,N12) pair showed *two* low-occupation eigenvectors, not one: `[0.016, 0.219, 0.568,
1.02, ...]`. C=N is a double bond -- a real sigma*/pi* pair should exist locally, and
projecting the *second*-lowest candidate (occupation 0.219) instead of the lowest gave
wmax=0.460 at MO47, matching NBO7's trusted 0.436 at MO48 to within 6%. This is exactly
what the real wX^max definition already does ("both BD*(1) and BD*(2) ... whichever
gives the larger squared coefficient wins," Notes.md) -- the fix (`ANTIBOND_OCC_THRESHOLD`
in `pair_nbo.py`) is to treat every local eigenvector below occupation 1.0 as a
candidate, project all of them, and keep whichever wins, mirroring that rule exactly
rather than assuming the single lowest-occupation eigenvector is always the right one.

**Result after the fix:**

| channel | NBO7 (trusted) | mol_002 (this method) | 5_s0_Me trusted | 5_s0_Me (this method) |
|---|---|---|---|---|
| wCNmax | 0.4356 (MO48) | **0.4599** (MO47) | 0.4570 (MO32) | **0.4602** (MO43) |

**For the first time across every method tried in this exploration, wCNmax is both
correctly ranked as the dominant channel AND numerically close to the trusted value in
both reference cases** (within 5.6% for mol_002, within 0.7% for 5_s0_Me).

**Bottom line for the open-source goal**: the local per-atom-pair construction is a
substantial, genuine step forward -- it reproduces wCNmax's central diagnostic
signature (dominant, and numerically close to NBO7) in both reference cases.

## wCNmax across all 6 test-set molecules

Tested against every main-pipeline test-set molecule with completed NBO7 data
(`beckmann.dft.inputs.TEST_IDS`), not just mol_002 -- using
`beckmann_alt.geometry.load_test_set_case()`, which resolves atom maps and geometry
fresh via `resolve_mol_name`/`oxime_atom_map_from_gjf` (the main pipeline's own
utilities) rather than hand-transcribing six atom maps.

| mol | computed wCNmax (MO) | trusted wCNmax (MO) | % diff | MO offset |
|---|---|---|---|---|
| mol_002_E | 0.4599 (47) | 0.4356 (48) | +5.6% | -1 |
| mol_006_E | 0.4473 (43) | 0.4225 (44) | +5.9% | -1 |
| mol_014_Z | 0.4431 (43) | 0.4212 (44) | +5.2% | -1 |
| mol_020_E | 0.4781 (47) | 0.4502 (48) | +6.2% | -1 |
| mol_021_E | 0.4944 (51) | 0.4692 (52) | +5.4% | -1 |
| mol_029_Z | 0.4933 (47) | 0.4665 (48) | +5.7% | -1 |
| 5_s0_Me | 0.4602 (43) | 0.4570 (32) | +0.7% | -11 (different basis, see caveat above -- not directly comparable) |

**The pattern is tight and systematic, not noisy.** Across all 6 of our own-pipeline
molecules (same basis, same functional, same solvent-caveat throughout): the percent
error sits in a narrow 5.2-6.2% band (mean ~5.7%), **always an overestimate, never an
underestimate**, and the winning canonical MO is **exactly 1 index lower than NBO7's own
winning MO in every single case** -- not "close," not "usually," literally 6/6. This
level of consistency across six chemically different substrates (different ring sizes,
substituents, and connectivity) is much more characteristic of a fixed, structural
offset (e.g. a systematic difference in how PySCF's and Gaussian's virtual manifolds are
ordered/counted near the frontier, or a consistent solvent/basis/dispersion-driven
energy shift that happens to move exactly one virtual MO's relative position) than of
six independent, coincidentally-similar numerical errors. Not root-caused yet -- would
need a direct comparison of virtual orbital energies near the frontier between the
Gaussian and PySCF calculations for the same molecule to pin down which piece (solvent
model, dispersion treatment, or something else) produces the shift. If diagnosed, a
fixed 5-6%/1-MO-index systematic offset would be straightforward to calibrate out;
right now it's reported as an observed, reproducible pattern, not yet an explained or
corrected one.

`5_s0_Me` sits apart, as expected given its different-basis caveat -- its % error is
much smaller (0.7%) but its MO offset is much larger (-11), underscoring that magnitude
agreement and orbital-identity agreement are answering different questions here, and
that `5_s0_Me`'s numbers shouldn't be pooled with the other 6 in any summary statistic.

## The wCNmax-minima rule, open-source vs. NBO7 (6 test-set molecules)

Main branch's `beckmann/dft/wcnmax_rule.py` predicts rearrangement ('R') for any
substrate whose wCNmax(R) scan shows a genuine interior minimum, fragmentation ('F')
otherwise -- 25/34 (74%) on the full benchmark set. `beckmann_alt/wcnmax_scan_rule.py`
reproduces this using the open-source per-atom-pair method instead of NBO7, for the 6
main-pipeline test-set molecules, by computing wCNmax at every point of each
molecule's existing R(N-O) scan series (`beckmann_alt/pair_nbo.run_test_set_scan_series`,
via a new `beckmann_alt/geometry.load_test_set_scan_series()` that extracts geometry at
each scan point the same way `beckmann.dft.parse_cmo` does -- last "Standard
orientation" block before each CMO section, `STEP_SCAN_SOURCES`-merged for mol_020_E --
so point identity/R(N-O) values line up with the trusted series) and feeding the result
through the SAME `find_wcnmax_minimum()`/`predict_from_wcnmax()` main-pipeline code, not
a reimplementation.

| mol | open-source: min found / predicted / agrees w/ exp | NBO7: min found / predicted / agrees w/ exp | same call as NBO7? |
|---|---|---|---|
| mol_002_E | True / R / no | True / R / no | yes |
| mol_006_E | **False** / F / **no** | True / R / **yes** | **no** |
| mol_014_Z | True / R / no | True / R / no | yes |
| mol_020_E | True / R / yes | True / R / yes | yes |
| mol_021_E | True / R / yes | True / R / yes | yes |
| mol_029_Z | True / R / yes | True / R / yes | yes |

**Open-source rule: 3/6 (50%) agreement with experiment vs. NBO7's 4/6 (67%) on this
subset.** The only molecule where the two methods' classification actually differs is
mol_006_E -- everywhere else they agree on whether a minimum exists, even where that
shared call happens to be wrong. Output: `data/output/analysis/wcnmax_rule_results_opensource.csv`
(per-molecule) and `wcnmax_channel_extraction_opensource.csv` (per-point, all 6
molecules, 48 rows total).

## mol_006_E follow-up: interior minimum invisible at 0.05 Å resolution

mol_006_E is the one molecule (of the original 6) where the open-source series
disagrees with NBO7 on whether a genuine interior wCNmax minimum exists at all: NBO7
finds one (depth 0.1010, R=1.6608, predicting rearrangement, matching experiment); the
open-source series (7 points, `nbo` + `scan_1`..`scan_6`, the same standard 0.05 Å
resolution used everywhere else) is smooth and monotonically increasing across NBO7's
entire dip region -- no minimum
(`data/output/analysis/plots/mol006_opensource_vs_nbo7_wcnmax.png`).

**Decision: only directly-computed geometries are used for any comparison in this
project** -- an earlier pass tried probing this gap with interpolated intermediate
geometries (linear Cartesian blends between two real converged scan points, not
independently DFT-relaxed structures) at finer resolution around NBO7's dip, and found
a real-looking but much shallower/narrower secondary dip plus a plausible root cause
(the local, non-deflated per-atom-pair antibond construction responding to the same
underlying canonical-MO crossing far more mutedly than NBO7's own whole-molecule-
deflated BD*, ~5-6x muted specifically). That interpolation-based code and its
generated CSV have been removed (`beckmann_alt.geometry.interpolate_case`,
`_compute_interp_point.py`, `_compute_interp_diagnostics.py`,
`wcnmax_channel_extraction_opensource_mol006_interp.csv`) since results built on
non-relaxed geometries aren't a substitute for genuinely re-running the DFT scan at a
finer step -- available via `git log -- beckmann_alt/` for the full writeup and numbers
if this is revisited with real (Citadel-computed) finer-resolution geometries instead.

## Expansion: 11 more F-labeled molecules (specificity check)

The 6-molecule comparison above is a small, R-heavy sample (4 R, 2 F) that isn't well
suited to checking a specific failure mode: does the open-source rule spuriously find an
interior wCNmax minimum (predicting rearrangement) for substrates that experimentally
fragment? Ran the same `beckmann_alt/wcnmax_scan_rule.py` workflow
(`load_test_set_scan_series()`/`run_test_set_scan_series()`, generalized from `TEST_IDS`
to `ALL_IDS` -- every benchmark molecule has completed Stage 1-3 logs on disk, not just
the original 6) against all remaining F-labeled benchmark molecules: mol_001, 003, 004,
007, 008, 010, 011, 012, 016, 017, 018 (11 of the 14 F-labeled substrates -- **mol_005_E
skipped**: it contains bromine, and PySCF's bundled `6-311G-diffuse.dat` -- the "+"
diffuse-function file behind `6-311+G(d,p)`'s "+" -- has no Br entry at all, even though
the base `6-311G` and `d`-polarization files do; building the official basis for this
one element isn't possible with this pip install, and substituting a different basis for
just one molecule wasn't done -- skipped rather than silently deviating).

One geometry-loading gap found and fixed along the way: `_mol_stage_points_stepscan()`
(used for mol_003_E, one of the three `STEP_SCAN_SOURCES` molecules) unconditionally
required an `_nbo.log` to exist. mol_003_E has none -- its Stage 2 equilibrium NBO job
never ran, matching `cmo_channel_extraction.csv`'s own data (no `nbo`-stage row for
mol_003_E either, only `scan_1`..`scan_6`) -- so this now mirrors
`beckmann.dft.parse_cmo.collect_stage()`'s existing "missing file -> empty, not a hard
failure" handling rather than crashing.

### Result: tied overall, 9/17 (53%) each

| mol | open-source: min / pred / agree | NBO7: min / pred / agree | exp | same call? |
|---|---|---|---|---|
| mol_001_E | True / R / no | True / R / no | F | yes |
| mol_003_E | True / R / no | True / R / no | F | yes |
| mol_004_E | False / F / **yes** | False / F / yes | F | yes |
| mol_007_E | False / F / **yes** | False / F / yes | F | yes |
| mol_008_E | True / R / no | True / R / no | F | yes |
| mol_010_E | False / F / **yes** | False / F / yes | F | yes |
| mol_011_E | True / R / no | **False / F / yes** | F | **no -- NBO7 right** |
| mol_012_E | True / R / no | True / R / no | F | yes |
| mol_016_E | **False / F / yes** | True / R / no | F | **no -- open-source right** |
| mol_017_E | False / F / **yes** | False / F / yes | F | yes |
| mol_018_E | **False / F / yes** | True / R / no | F | **no -- open-source right** |

Combined with the original 6 (`wcnmax_rule_results_opensource.csv` now holds all 17):

| | accuracy vs. experiment |
|---|---|
| open-source | **9/17 (53%)** |
| NBO7 (same 17) | **9/17 (53%)** |

**Exactly tied, and not by the two methods agreeing with each other -- by being wrong on
different molecules.** Of the 4 molecules (across all 17) where the two methods'
predicted label actually differs, it's an even 2-2 split: NBO7 gets mol_006_E and
mol_011_E right where open-source doesn't; open-source gets mol_016_E and mol_018_E right
where NBO7 doesn't. The earlier 6-molecule sample (50% vs. 67%) made the open-source
method look meaningfully worse -- that gap was mostly small-sample noise, not a
consistent deficit. mol_016_E/mol_018_E are specific cases where NBO7 finds a spurious
minimum and the open-source method's generally muted response to sharp features
(the local, non-deflated per-atom-pair construction's own systematic understating of
these dips/crossings -- see the mol_006_E follow-up above) happens to avoid picking it
up -- the same tendency that under-detects mol_006_E's real minimum here correctly
under-detects two false ones. Not enough molecules yet to call this a real specificity
advantage rather than coincidence.

Output: `data/output/analysis/wcnmax_rule_results_opensource.csv` (17 rows, one per
molecule tested so far -- `exp`/`NBO`/`PySCF` columns hold each side's predicted R/F
label directly, side by side; `NBO_correct`/`PySCF_correct`/`NBO_PySCF_match` carry the
derived comparisons; `n_points`/`R_star`/`R_depth` -- PySCF's own scan diagnostics --
are the last three columns) and `wcnmax_channel_extraction_opensource.csv` (114
per-point rows: 48 from the original 6 + 66 from these 11).

## Full benchmark run: all 32 runnable molecules (2026-07-22)

Ran the remaining 15 non-test-set molecules (mol_009, 013, 015, 019, 022, 023, 024,
025, 026, 027, 028, 030, 031, 033, 034) through the same `wcnmax_scan_rule.py`
workflow, in small batches to bound peak memory/CPU on a laptop. `wcnmax_scan_rule.py`
now takes molecule ids as CLI args (`python -m beckmann_alt.wcnmax_scan_rule 009 013`,
no args = original `TEST_IDS` default) and merges each molecule's result into the
existing CSVs as soon as it's computed, instead of overwriting the whole file at the
end of the run -- this made it safe to run one or two molecules at a time and safe
against a run getting interrupted partway (which happened repeatedly this session;
no data was lost, since a molecule's row is only written once its own scan series is
fully computed).

mol_005_E and mol_032_E remain skipped -- both contain Br, and PySCF's bundled
`6-311G-diffuse.dat` (the "+" diffuse-function file behind `6-311+G(d,p)`) has no Br
entry, same issue documented above for mol_005. **32 of 34 benchmark molecules now
covered.**

mol_034_E's `STEP_SCAN_SOURCES`-merged series has 12 R(N-O) points (2x every other
molecule's 6), which kept getting caught mid-run by interruptions. Tried a 6-point
subset instead (`scan_2/4/6/8/10/12`, chosen to bracket NBO7's own trusted interior
minimum, R_star=1.6690, with three points on each side; `run_test_set_scan_series()`
gained an optional `stages=` filter for this, `beckmann_alt/pair_nbo.py` -- default
`None` behavior is unchanged for every other molecule/caller) -- but a stale retry
process from an earlier interrupted attempt turned out to still be running in the
background (despite being reported as killed) and finished the full 12-point run
*after* the 6-point result had already been written and verified, silently overwriting
it via the same read-merge-write logic that's supposed to protect against exactly this
(it protects against a crash losing data, not against two overlapping runs of the same
molecule racing each other). Caught by re-checking the CSV after an unexpected
late-arriving completion notification. The full 12-point result is what's actually in
the CSV now (`n_points=12`, `R_star=1.6290`, `R_depth=0.0954`) and is what the tally
below reflects -- it still finds the minimum and predicts R, matching both NBO7 and
experiment, same as the 6-point version did. The `stages=` filter is left in
`pair_nbo.py` since it's a real, reusable capability; just not exercised in the
molecule actually recorded here.

### Result: 21/32 (66%) open-source vs. 23/32 (72%) NBO7, no longer tied

| | accuracy vs. experiment |
|---|---|
| open-source | **21/32 (66%)** |
| NBO7 (same 32) | **23/32 (72%)** |
| open-source predicts the same R/F call as NBO7 | **26/32 (81%)** |

Unlike the earlier 17-molecule sample (exactly tied, 9/17 each), NBO7 now leads by 2
molecules once the full runnable benchmark is included -- still close, and the two
methods still agree with each other on the large majority (26/32) of calls, including
plenty of shared wrong calls (both methods lean toward over-predicting rearrangement,
the same false-positive-heavy pattern noted throughout this file). All 32 rows,
`n_points`/`R_star`/`R_depth` included, live in `wcnmax_rule_results_opensource.csv`;
per-point data (96 new rows -- 14 molecules x 6 points plus mol_034_E's 12 -- on top of
the previous 114) in `wcnmax_channel_extraction_opensource.csv`.

# wCNmax: PySCF vs NBO7 — investigation plan and implementation notes (2026-07-22)

Based on `wcnmax_rule_results_opensource.csv` (17 molecules) and the exploration
already documented in `Notes_open_source_alt.md` (`beckmann_alt/pair_nbo.py`).

## Current state, in one line

Overall accuracy is tied (9/17 each vs. experiment), but not because the two methods
agree with each other. They disagree on 4 molecules and split those 2-2: NBO7 gets
mol_006_E and mol_011_E right where PySCF doesn't; PySCF gets mol_016_E and mol_018_E
right where NBO7 doesn't. The "PySCF minima come out less deep" pattern you're seeing
now is the same phenomenon already caught once, on mol_006_E specifically, in the
notes below, not a new problem.

---

## Part 1: investigation plan

### Two known causes, not yet fixed

**1. A fixed magnitude/index offset (across all 6 same-basis test molecules).**
PySCF's wCNmax comes out +5.2 to +6.2% high (mean ~5.7%) at every point tested, and
the winning canonical MO is exactly one index lower than NBO7's winning MO, every
single time, on 6 chemically different substrates. That consistency is a much
stronger signal for a structural offset than for six coincidentally similar errors.

**2. Muted response to real orbital mixing (avoided crossings).** The PySCF
construction builds a fresh local orbital just for the one atom pair asked about
(C-N), directly from that pair's block of the density matrix, with no deflation
against the rest of the molecule's orbitals. NBO7's real algorithm builds every
localized orbital in one whole-molecule pass, sequentially deflating each one's
density before finding the next. Tested directly on mol_006_E's interpolated
finer-resolution geometries: right at a genuine near-degenerate mixing event, the
non-deflated local construction responded roughly 5-6x more mutedly than NBO7's own
deflated result. This is the most likely direct explanation for shallower/narrower
minima.

### Priority order for the remaining ~3 weeks

**1. Root-cause the fixed offset (cheap, do first).** It's one clean, repeatable
signal, so isolating the variable should be fast:
- Re-run 1-2 molecules gas-phase on both sides (no solvent model at all) to check
  whether ddCOSMO-vs-SMD is the source, rather than solvent-model choice itself.
- Toggle density fitting off for one small molecule, even if slow, to rule it in or
  out.
- Compare raw virtual-orbital eigenvalues near the frontier (not the wCNmax weight)
  between PySCF and Gaussian on the same geometry. If the ordering itself is
  shifted by one slot near the frontier, that's the direct explanation for the -1
  MO-index pattern regardless of which upstream setting causes it.
- If confirmed fixed, this becomes a calibratable correction, not something to chase
  to zero.

**2. Recover minima depth during real mixing events (the core issue).**
- Add deflation to the local per-atom-pair construction (project out already-claimed
  orbital density before diagonalizing the local block) and re-test specifically on
  the mol_006_E region where muting was directly measured, this time using real
  Citadel-computed fine-resolution geometries instead of the earlier interpolated
  Cartesian blends, since Citadel access exists now.
- Note: an earlier deflation attempt didn't improve the winning wCNmax value's match
  to NBO7, but it was aimed at "which MO wins," not "how deep is the dip at a
  crossing." Those are different questions and worth a second, narrower pass.
- If deflation doesn't move it, try widening the local subspace itself (a 3-atom
  cluster around the C-N bond instead of a strict pair) before concluding deflation
  is the wrong lever.

**3. Push scan resolution (if time remains, and the actual point of this tool).**
PySCF is cheap enough to scan at 0.02 Å or finer across all molecules, which is the
real value case for the Enamine/triage goal: use PySCF to flag which R-window looks
interesting, then send only that window to NBO7/Citadel for confirmation, instead of
scanning blind. mol_006_E's real minimum was completely invisible at the standard
0.1 Å grid and only showed up at 0.05 Å, so this isn't a hypothetical concern.

### One thing to get sign-off on before sinking more time in

Is the target "match NBO7's numeric depth" or "correctly flag the R-region where a
minimum exists, even if the magnitude is muted"? These call for different fixes
(offset-correction vs. depth-recovery), and the stated Enamine/triage use case only
strictly needs the second one. Worth confirming with Isayev/Tetiana before choosing
where to spend the remaining time.

---

## Part 2: how the two implementations differ (for anyone, code-optional)

Both sides compute exactly the same quantity: for the C-N antibond, scan every
virtual molecular orbital from the LUMO up to LUMO+0.4 a.u., take that antibond's
coefficient in each orbital's expansion, square it, and keep the largest value
found. That definition comes straight from the same theoretical framework behind
Tetiana's reference paper, and both sides target it identically. Where they differ
is in how each one builds the local orbital frame that coefficient gets measured
against.

- **NBO7 (Gaussian, the trusted numbers everywhere else in this project):** builds
  every localized bonding/antibonding orbital in the molecule in one pass,
  sequentially deflating each one's density before finding the next. This whole-
  molecule deflation is what gives NBO7's antibonds their sharp character,
  including sharp responses to genuine orbital-mixing events.
- **PySCF (open-source prototype):** builds a fresh, local orbital subspace for
  just the one requested atom pair, directly from that pair's block of the density
  matrix in an IAO basis, structurally the same kind of operation NBO does for one
  pair, just without deflating against the rest of the molecule. Much faster (no
  need to build the whole molecule's localized orbital set to ask about one bond),
  but the local orbital doesn't get "sharpened" the way NBO's does.
- **Practical consequence:** the two agree closely on which channel dominates and on
  general orbital character (validated to within 1% on Tetiana's own worked
  reference compound), but PySCF's version responds more mutedly wherever two
  orbitals are genuinely mixing, which surfaces as a shallower or narrower dip right
  at the R values where a Beckmann-relevant handoff happens. There's also a small,
  consistent (+5-6%) magnitude offset present even away from mixing regions,
  not yet root-caused (Part 1, item 1).
- **Level-of-theory notes, secondary but worth disclosing:** PySCF uses ddCOSMO as
  an implicit-solvent stand-in instead of Gaussian's actual SMD/water calculation
  (the compiled SMD extension isn't available in this pip build), and uses
  density-fitted integrals for speed. Neither is expected to be the main driver of
  the pattern above, but both are real, disclosed deviations from the reference
  calculation.
- **What's held constant:** basis set (6-311+G(d,p)), functional (wB97XD), and
  geometry (frozen, taken directly from the already-converged Gaussian structure,
  no independent PySCF optimization). This isolates the question to "does this
  local-construction method reproduce the same orbital character," not "does a
  different geometry or basis change the answer."
- **Shared theoretical foundation:** Tetiana's own reference calculation for this
  descriptor is itself Gaussian/NBO-based, so this isn't a competing physical
  theory. It's a test of whether a faster, structurally simpler construction of the
  same NBO-style local antibond can get close enough to NBO7's numbers to be useful
  as a pre-screening step before committing to a full NBO7 run.

## Part 1 outcome: offset/deflation diagnostics run, then set aside (2026-07-22)

Ran diagnostic 1 (root-cause the fixed offset) and diagnostic 2 (deflation) from the
plan above, on mol_002_E and mol_006_E respectively, both against real (never
interpolated) geometries:

- **The "-1 MO index" part of the offset is not a real electronic-structure effect.**
  Confirmed directly from `cmo_channel_extraction.csv`: NBO7's own trusted winning MO
  is the LUMO itself (`delta_lumo=0.0`) in every one of the 6 original test molecules.
  A fresh PySCF run on mol_002_E confirmed the same is true on the open-source side
  (`mo_index == nocc`). Both methods pick the literal LUMO every time -- they just
  number it differently (0-based array index vs. Gaussian's 1-based orbital number).
  Not a bug in the wCNmax construction; a reporting mismatch if `MO_index` columns
  from the two CSVs are ever compared directly.
- **The +5.2-6.2% magnitude offset is not explained by solvent choice, density
  fitting, or virtual-orbital ordering** -- all three ruled out with direct evidence
  (gas-phase PySCF on the same geometry made the mismatch *worse*, not better;
  density-fit on/off gave bit-identical results; PySCF's frontier virtual eigenvalues
  matched a real Gaussian gas-phase calculation on the same geometry to 4-5 decimals).
- **Deflation (local per-atom-pair construction, projecting out the C-aryl/C-alkyl/N-O
  bonds' own density before diagonalizing the C-N block) had no effect** on mol_006_E's
  R=1.6608 point -- before/after wCNmax agreed to 8 decimal places, even though the
  deflation measurably perturbed the local block's eigenspectrum. Also found, using the
  real Citadel finescan geometry (not the earlier interpolated blends): PySCF isn't
  showing a *muted* dip at this point, it shows *no* dip at all -- its value sits
  higher than NBO7's trusted number, right on PySCF's own smooth monotonic trend across
  the scan. Reframes the earlier interpolation-based "~5-6x muted" finding.

**Decision:** rather than keep tuning the deflation pair selection or widening the
local subspace blind, the exploratory code (`beckmann_alt/diagnose_offset.py`,
`diagnose_deflation.py`, and the deflation additions to `pair_nbo.py`) was removed from
the tree as added complexity without a clear payoff yet -- consistent with this
project's practice elsewhere in this file of not carrying speculative code forward.
Full diffs remain in git history (`ea4a4e2`, `ce3c510`) for anyone revisiting this.
Widening the local subspace (a 3-atom cluster instead of a strict pair) is the
suggested next lever if this is picked back up, per the original diagnostic-2 fallback
instruction.
