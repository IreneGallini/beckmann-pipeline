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
