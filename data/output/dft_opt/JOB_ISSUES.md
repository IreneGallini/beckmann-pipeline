# DFT job troubleshooting log

Add a new dated entry any time a Gaussian
job on Citadel crashes or fails to converge, don't just re-run and forget
why. 

---

## 2026-07-09/10 — mol_020_E `_scan.gjf` (Stage 3, N-O relaxed scan)

**Symptom:** job did not reach "Normal termination". Log ends with:
```
Error termination via Lnk1e in /opt/g16/l9999.exe at Fri Jul 10 00:31:25 2026.
Segmentation fault (core dumped)
```
Step count: 186 lines match `Step number` (scan point 5/5 alone used 162 of
its own step budget). Elapsed time: 6h37m.

**Root cause:** scan points 1-4 converged normally (few steps each, same as
mol_002/006/021). Point 5 (R = R0+0.4 Å, the most N-O-stretched geometry)
never converged — it oscillated between two nearly-degenerate ring-pucker
conformations of the saturated ring fused to the oxime carbon (the
benzosuberone-type CH2-CH2-CH2 bridge, atoms 8-10 in the `.gjf` numbering),
alternating every step:

| State | Energy (a.u.) | Max Force | C11-N12 bond |
|---|---|---|---|
| A | ≈ -557.30309 | ~0.031 | 1.2555 Å |
| B | ≈ -557.30361 | ~0.009 | 1.2820 Å |

Neither state ever approached the 0.00045 threshold; displacement was
essentially frozen at ~0.068 Å (target 0.0018 Å) for the last ~40+ steps
classic double-well oscillation, not slow monotonic convergence. The N-O
distance itself (the frozen scan coordinate) was correctly held fixed
throughout, so this is a ring-pucker degeneracy problem, unrelated to the
reaction coordinate.

After hitting the step cap, Gaussian proceeded anyway to the `pop=nbo7read`
single-point at that unconverged geometry — NBO7 actually completed
successfully (E2PERT, CMO, NBOSUM all present in the log) but the job's
exit/cleanup link (`l9999.exe`) then segfaulted, so the log ends in
`Error termination` rather than `Normal termination` even though usable-looking
NBO data is sitting in the file. **Do not trust that data** — it's computed on
an oscillating, non-converged geometry and was excluded from
`channel_descriptors.csv` / `cmo_descriptors.csv` / `nbo_e2pert.csv` /
`descriptor_slopes.csv` (see `scripts/dft/parse_nbo.py` /
`parse_cmo.py` — currently hardcoded to `TEST_IDS` in `beckmann/dft/inputs.py`
and don't filter by termination status, so this had to be stripped manually
post-hoc).

Also present throughout the whole scan (not just the failing point):
`Error on total polarization charges = 0.04625` — an SMD cavity-fitting
warning, repeated at every step. Didn't appear to be the direct cause (points
1-4 had the same warning and converged fine) but worth watching given it's
new since the solvent model was added (`SOLVENT = "scrf=(smd,solvent=water)"`
in `beckmann/config.py`, committed 2026-07-09).

**Fix plan (not yet executed):**
1. Pull the converged point-4 geometry from `mol_020_E_scan.log` (last
   "Standard orientation" block before point 5 starts) as the starting
   structure — don't restart from the crashed point-5 geometry.
2. Write a standalone single-point relaxed optimization for just R = R0+0.4 Å:
   freeze the N-O bond at the target distance with a `B ... F` ModRedundant
   line (not `S ... steps`, since we only need this one point, not a scan),
   `pop=nbo7read` as usual.
3. Add oscillation-damping options to the route line:
   `opt=(ModRedundant,MaxCycles=300,CalcFC,NoGDIIS,MaxStep=3)`
   - `CalcFC` computes an initial Hessian instead of a guessed/updated one —
     better curvature info in a flat region.
   - `NoGDIIS` disables Gaussian's default GDIIS step extrapolation, which is
     the most common cause of exactly this kind of two-state bounce in a flat
     PES.
   - `MaxStep=3` caps the step size (default is much larger) to stop the
     optimizer from jumping straight past the true minimum every time.
4. If it still oscillates: identify and freeze the specific ring dihedral
   that's flipping (candidates: dihedrals through C8-C9-C10-C11, the atoms
   with the largest coordinate deltas between states A/B) for this point only,
   converge, then release the constraint and re-verify with a short follow-up
   opt.
5. Re-run `scripts/dft/parse_nbo.py`, `parse_cmo.py`, `descriptors.py`,
   `scripts/analysis/summarize_descriptors.py` for all 4 test molecules once
   fixed (they currently only cover mol_002/006/021 — mol_020 rows were
   manually stripped from every output CSV on 2026-07-10).

**Resolution (2026-07-13):** Executed steps 1-2 with a lighter route line than
step 3 proposed — `opt=(ModRedundant,CalcFC,MaxCycles=300)` (no `NoGDIIS`, no
`MaxStep=3`) — as `mol_020_E_sp5.gjf`. Starting geometry: point 4's converged
Cartesian block from `mol_020_E_scan.log`, with the leaving group
(`O13`, `H26`, `H27` — the only atoms bonded to `O13` besides `N12`, per the
scan log's `$CHOOSE` block) rigidly translated +0.1 Å along the N12→O13 unit
vector to seed R(N-O) = 1.9029 Å, then `B 12 13 F` to hold that bond fixed
while everything else relaxes. Converged in 58 minutes (`Optimization
completed` / `Stationary point found`), clean NBO7 (98.01% Lewis structure,
no low-occupancy warnings) — `CalcFC` + a good initial guess was enough on its
own; didn't need step 4's dihedral-freeze fallback.

**Second issue found during recovery:** `mol_020_E_sp2/sp3/sp4.log` (dated
2026-07-07) turned out to be stale — their own title lines report R(N-O) =
1.7259/1.8259/1.9259 Å (a R0'=1.6259 Å baseline), not the 1.6029/1.7029/1.8029
Å expected from the current scan's R0 = 1.5029 Å. Stage 1 (`_opt.gjf`) must
have been re-optimized to a different converged geometry after these sp2-4
jobs were originally extracted, silently orphaning them — mixing them with
the current `nbo`/`scan`/`sp5` data would have spliced two different geometry
series into one non-monotonic "R(N-O) scan". Fixed by re-running
`scripts/dft/extract_scan_sp.py --mol 020` against the current
`mol_020_E_scan.log` (correctly reproduced R = 1.6029/1.7029/1.8029 Å) and
resubmitting just those three single-point jobs — sp5 was left alone since it
was already good. **Lesson: if a molecule's Stage 1 is ever re-run, its
sp2/sp3/sp4 (and any restart-point jobs like sp5) must be regenerated too —
nothing currently checks that they came from the same `_opt.chk` lineage.**

Pipeline code changes (`beckmann/dft/parse_nbo.py`, `parse_cmo.py`,
`descriptors.py`): added `sp5` as a recognized stage that supersedes
`_scan.log` entirely when present and clean (scan's point 1 duplicates `nbo`,
point 5 is what sp5 replaces) — `_scan.log` is no longer required to reach
Normal termination for a molecule that has a working `sp5.log`.
`SERIES_STAGES`/`SERIES_FALLBACK` in `descriptors.py` swap in `sp5` for
`scan_2` in the 5-point d/dR series when `scan_2` isn't available.

---

## General playbook: crashed or non-converging Gaussian jobs

**How to check status:**
- Last line of `.log` must be `Normal termination of Gaussian 16 at ...` —
  anything else (`Error termination`, no final line at all, a bare
  `Segmentation fault`) is a failure, even if earlier sections of the log
  look complete.
- `grep -c 'Step number' {name}_scan.log` vs the "out of a maximum of N"
  figure printed on each line tells you how close an optimization got to its
  step cap before stopping.
- `grep 'Maximum Force' {name}*.log` gives the convergence trend. Compare
  against a molecule that succeeded (similar size/system) — a healthy
  trajectory decays toward the YES/YES/YES/YES threshold within single-digit
  to ~30 steps per scan point; anything grinding past ~50 steps on one point
  is worth a closer look before assuming it'll "get there eventually."

**Distinguishing failure modes from the Max Force sequence:**
- **Monotonic but slow** (values steadily shrink, just not fast enough) →
  usually just needs more `MaxCycles`.
- **Noisy but trending down** (mol_006_E's scan point 4: spikes to 0.02-0.08
  but decaying) → often self-resolves given enough steps; not usually worth
  intervening unless it's burning through most of the step budget.
- **Clean 2-state alternation, same two Force/Energy values repeating with no
  trend** (mol_020_E's scan point 5) → a genuine double-well degeneracy
  (commonly ring pucker in fused non-aromatic rings, especially 7-membered
  ones — flatter pseudorotation landscape than 6-membered rings). More steps
  will not fix this. Needs `NoGDIIS`/`CalcFC`/reduced `MaxStep`, or freezing
  the offending internal coordinate.

**If the job crashes (not just fails to converge) after an apparently-successful
computation section:** check whether the crash is in the actual electronic
structure/optimization work or in a later link (e.g. `l9999.exe`, the exit/
archive step). A crash after "NBO analysis completed" or after an "SCF Done"
line means real data may be sitting in the log even though the job overall
reads as failed — but treat that data as suspect if it came from a geometry
that hadn't converged (check the preceding `Maximum Force` block), and exclude
it from downstream CSVs rather than silently including it.

**Diagnosing what's physically oscillating:** pull the "Standard orientation"
coordinate block from two steps on opposite sides of an oscillation cycle and
diff them per-atom (see the mol_020_E entry above for the method). The atoms
with the largest Cartesian displacement between the two states tell you which
internal coordinate to freeze or investigate — don't assume it's the
reaction coordinate itself without checking; here it was a ring bridge two
bonds away from N-O, not the scanned bond.

**Solvent-model-specific:** `Error on total polarization charges` is an SMD
cavity-fitting message, not necessarily fatal on its own (seen in converged
jobs too) — but new since `SOLVENT` was added to `beckmann/config.py`
(2026-07-09), so if a molecule that converged fine in gas phase starts
oscillating or crashing under `scrf=(smd,solvent=water)`, check whether the
implicit solvent cavity is doing something pathological at that geometry
before assuming it's an unrelated optimizer issue.

---

## 2026-07-16/18 — mol_020_E, mol_003_E, mol_016_E, mol_023_E, mol_030_E, mol_032_E — non-convergence crashes under the new 6-point/0.05 Å architecture — fix attempted, REVERTED pending PI input

**STATUS: OPEN.** A partial fix (CalcFC restart of only the crashed point,
spliced back into the other 5 unmodified points) was attempted for
mol_020_E and started for mol_003_E, then **reverted** — see "Why the fix
was reverted" below. All five molecules below are left in their crashed,
unresolved state (mol_023_E and mol_030_E crashed after the revert
decision, so neither was patched at all — straight to "leave it and
document it"). Do not re-attempt a per-point CalcFC patch without
re-reading this entry; the methodology question below needs an answer
first (from the PI), not another one-off fix.

Six crashes out of the 34 benchmark substrates processed so far (all 34
now have at least Stage 1 attempted; all 34 have reached the Stage 3 scan)
hit the same non-convergence signature — roughly 1 in 6. This is clearly a
recurring structural tendency (fused-
ring pucker oscillation, see the 2026-07-09/10 entry above for the original
mol_020_E case under the old architecture), not a one-off fluke — worth
raising as a general question,
not fixing molecule-by-molecule.

**mol_020_E** — points 1-2 converged normally (4/12 Normal termination).
Point 3's optimization (`R0+0.15 Å`) oscillated — `Maximum Force`
alternating between ~0.030 and ~0.012 for 30+ consecutive steps, never
trending down — hit its step budget (156 of 156), then errored via
`l9999.exe` and **segfaulted (core dumped)**, killing the whole job (points
4-6 never ran since it's one Gaussian process across all 12 Link1 blocks).
Crashed log (unmodified, as it happened):
`data/output/dft_opt/mol_020_E/mol_020_E_scan.log` (also duplicated at
`data/output/dft_opt/_archive_pre_6pt_scan/mol_020_E_crashed_attempt/mol_020_E_scan_crashed.log`).

**mol_003_E** — points 1-5 converged normally (10/12). Point 6's
optimization (`R0+0.30 Å`, the very last point) oscillated with the
identical signature (~0.0201/~0.0210 alternating for 20+ steps), then the
same `l9999.exe` error + segfault. Crashed log:
`data/output/dft_opt/_archive_pre_6pt_scan/mol_003_E_crashed_attempt/mol_003_E_scan_crashed.log`.

**mol_016_E** — this one crashed during **Stage 1** (the initial geometry
optimization, before the 6-point scan even starts), same oscillation
signature, caught proactively (still running, no convergence, after 4+
hours against the ~15-30 min every other Stage 1 job took) rather than
after an actual crash message. **The original crashed log for this one is
lost** — it was inspected live over SSH, then the job was killed and
resubmitted with the same output filename (`mol_016_E_opt.log`), which
overwrote it via shell redirect before a local copy was ever saved. Lesson
for next time: always `scp`/archive a crashed or suspect log to a local or
archived path *before* killing/resubmitting a job that writes to the same
filename — don't rely on being able to re-fetch it after the fact.

**mol_023_E** — crashed early: only point 1 converged
(2/12). Point 2's optimization (`R0+0.10 Å`) oscillated — `Maximum Force`
alternating between ~0.0091 and ~0.0348 for 20+ consecutive steps, never
trending down — same `l9999.exe` error + segfault. This one happened
*after* the revert decision below, so no fix was attempted at all — left
crashed and documented directly. Crashed log:
`data/output/dft_opt/_archive_pre_6pt_scan/mol_023_E_crashed_attempt/mol_023_E_scan_crashed.log`.

**mol_030_E** — same pattern again: only point 1 converged (2/12). Point
2's optimization (`R0+0.10 Å`, same point number as mol_023_E) oscillated
— `Maximum Force` alternating between ~0.0090 and ~0.0280 for 20+
consecutive steps — same `l9999.exe` error + segfault. Also happened after
the revert decision; left crashed, no fix attempted. Crashed log:
`data/output/dft_opt/_archive_pre_6pt_scan/mol_030_E_crashed_attempt/mol_030_E_scan_crashed.log`.

**mol_032_E** — same pattern, point 2 again: only point 1 converged
(2/12). `Maximum Force` alternating between ~0.0159 and ~0.0233 for 20+
consecutive steps — same `l9999.exe` error + segfault. Left crashed, no
fix attempted. Crashed log:
`data/output/dft_opt/_archive_pre_6pt_scan/mol_032_E_crashed_attempt/mol_032_E_scan_crashed.log`.
Notably its batch-6 sibling mol_033_E (also 24 atoms, both `[oxime:
C3=N2-O1]`) converged cleanly — this isn't simply "these two molecules
share a scaffold," so whatever's structurally triggering the oscillation
is more specific than atom count or oxime position alone.

**Why the fix was reverted:** the initial response (CalcFC + MaxCycles=300
on just the crashed point, chained through the other unmodified points,
spliced into one log) fixed mol_020_E's point 3 and was in progress for
mol_003_E's point 6 — but this leaves 5 of that molecule's 6 scan points
computed under one setting and 1 point under another, within the same
R(N-O) series that downstream descriptors (Ψ, Λ, wCNmax) treat as a single
continuous trend. That's a methodological inconsistency, not just a
technical fix — flagged directly by the user, who wants the PI's
opinion on the right general policy (e.g., "if any point in a molecule's
scan needs CalcFC, rerun all 6 points of that molecule with it" vs. some
other rule) before any more crashed points get patched one at a time.
Both in-flight fix jobs (mol_020_E's resume, mol_003_E's resume, mol_016_E's
CalcFC Stage 1 restart) were killed; mol_020_E's and mol_016_E's canonical
`.gjf` files were reverted to their original (non-CalcFC) route lines;
mol_020_E's canonical `.log` was restored to the original crashed content.
Partial/abandoned resume files kept in
`data/output/dft_opt/_archive_pre_6pt_scan/{mol}_crashed_attempt/` for
reference, not deleted.

**Open question for the PI:** what's the standard/expected way to
handle a non-converging point within a relaxed scan — rerun the whole
series with stronger optimizer settings uniformly, rerun just the failed
point with different settings (accepting the inconsistency), freeze the
specific oscillating internal coordinate, or something else? Send the PI the
five crashed logs above (mol_020_E, mol_003_E, mol_023_E, mol_030_E,
mol_032_E — mol_016_E's is unfortunately unrecoverable, see its entry
above).

**Takeaway for future batches:** mol_020_E-style fused-ring pucker
oscillation is a per-molecule structural tendency, not tied to any specific
scan point or architecture — expect it to keep recurring on other
substrates with similar non-aromatic fused rings (6 for 34 so far, roughly
1 in 6 — and not simply predictable from atom count or oxime label, see
the mol_032_E/mol_033_E comparison above). **Do
not apply the CalcFC-restart-one-point fix anymore** — that's the
methodology this whole entry exists to flag as inconsistent. Until the
the PI responds, the correct action on a new crash is: diagnose (confirm
it's the same oscillation signature, not something new), archive the
crashed log to `_archive_pre_6pt_scan/{mol}_crashed_attempt/`, add a
one-paragraph entry here, and leave the molecule crashed. Do not
patch, splice, or resubmit it.

**Update (2026-07-18/20) — uniform-CalcFC test results, still STATUS: OPEN.**
Despite the "do not patch" note above, a test was run anyway (commit
`b4f199f "test crashed mols with CalcFC"`) applying
`opt=(ModRedundant,CalcFC,MaxCycles=300)` uniformly across **all** scan
points of a molecule (not just the failing one, to avoid mixing optimizer
settings within a single R(N-O) series) for the five molecules above. This
was explicitly framed as a test pending the PI's answer, not a
policy change — logging results here for when they respond:

| Molecule | Original crash point | With uniform CalcFC |
|---|---|---|
| mol_023_E | pt2 (R0+0.10Å) | **Normal termination** |
| mol_030_E | pt2 (R0+0.10Å) | **Normal termination** |
| mol_032_E | pt2 (R0+0.10Å) | **Normal termination** |
| mol_020_E | pt3 (R0+0.15Å) | crashed pt6 (R0+0.30Å) — delayed, not fixed |
| mol_003_E | pt6 (R0+0.30Å) | crashed pt5 (R0+0.25Å) — no improvement (slightly worse) |

So the fix is substrate-dependent: it fully resolves the oscillation for 3
of 5, and for the other 2 it shifts which point fails without eliminating
the failure. This is useful data for the PI's methodology question
but is not evidence that uniform CalcFC is a general solution.

**mol_034_E — a 7th substrate with the same crash signature, found
2026-07-18/20, never CalcFC-tested.** Crashed at pt5 (R0+0.25Å),
`Maximum Force` alternating with no convergent trend, same
`l9999.exe`-then-segfault ending as every other case in this family.
Original (non-CalcFC) `.gjf`/route line, untouched. Total count is now
7 of 34 (~1 in 5), not 6 of 34.

**mol_016_E — Stage 1 resubmit succeeded.** The original Stage 1 crash
(entry above, log unfortunately unrecoverable) was resubmitted and reached
Normal termination on 2026-07-18 16:09. Stage 3 (`_scan.gjf`) has not been
generated yet — this one doesn't touch the methodology question since it's
the molecule's first-ever scan attempt at default settings, same as every
other substrate's initial try.

**PI response (2026-07-20) — step-size shift, not stronger optimizer
settings.** Rather than CalcFC/NoGDIIS/MaxStep, they suggested trying
different scan step increments for the molecules still failing at pt5 —
"another increments, say 0.7 or 0.4" (read as 0.07/0.04 Å, consistent with
our 0.05 Å scale) — the idea being to sample R(N-O) points that land off the
exact bond length where the ring-pucker double-well sits, rather than
fighting the optimizer at that geometry.

**Step-size test results (2026-07-20) — STATUS: RESOLVED for mol_020_E and
mol_034_E, PARTIALLY RESOLVED for mol_003_E.** Generated via
`prepare_scan_rigid(..., step=0.07/0.04, n_points=6)` in a new side-experiment
directory `data/output/dft_opt_stepscan/{mol}_step07/` and `..._step04/`
(same pattern as `dft_opt_631g/`, `dft_opt_finescan/` — copies of the
molecule's `_opt.gjf`/`_opt.log`, original settings otherwise, no CalcFC).
6 of 7 jobs succeeded:

| Molecule | 0.07 Å step | 0.04 Å step |
|---|---|---|
| mol_003_E | ❌ crashed — same oscillation at pt5 (R0+0.35Å), segfaulted after 3h12m stuck there | ✅ **Normal termination**, all 6 points |
| mol_020_E | ✅ **Normal termination**, all 6 points | ✅ **Normal termination**, all 6 points |
| mol_034_E | ✅ **Normal termination**, all 6 points | ✅ **Normal termination**, all 6 points |

So the ring-pucker degeneracy for mol_003_E sits close enough to the
default-step R(N-O) sampling that a 0.07Å offset still lands on it, while
0.04Å avoids it — a narrower miss than for mol_020_E/mol_034_E, where either
offset cleared the problem. All logs downloaded to
`data/output/dft_opt_stepscan/{mol}_step0{4,7}/`.

**Consistency question — RESOLVED (2026-07-21).** Per the user's follow-up
decision: mol_020_E and mol_034_E's two independently-successful scans
(step07 and step04) are **merged into one denser per-molecule series**, not
kept as two parallel rows — the scan's actual purpose is resolving whether
wCNmax is monotonic or has an interior minimum, and combining every
successfully-converged point is strictly better for that than picking one
series and discarding the other's points. Implemented via
`STEP_SCAN_SOURCES`/`build_stage_relabel_map()`/`relabel_rows()` in
`beckmann/dft/inputs.py`: `parse_nbo.py`/`parse_cmo.py`/`parse_wiberg.py`'s
`collect_molecule_stepscan()` pulls 'nbo' rows from the canonical
`dft_opt/{mol}/{mol}_nbo.log` (Stage 2 succeeded independently of the Stage 3
crash) and 'scan' rows from every listed `dft_opt_stepscan/` source, then
renumbers the combined scan points as one `scan_1..scan_N` sequence sorted by
actual R(N-O) — `resolve_series()`/`compute_slopes()` needed no changes,
since they already just sort by the `scan_N` stage suffix regardless of
count or spacing. Result: `mol_020_E` now has 13 points (1 nbo + 12 merged
scan), `mol_034_E` has 12 (no nbo log ever run for it — see below),
`mol_003_E` has 6 (step04 only, no merge needed).

The original crashed logs (`dft_opt/mol_003_E/`, `dft_opt/mol_020_E/`,
`dft_opt/mol_034_E/`, and the earlier uniform-CalcFC attempts) are left in
place, untouched — not archived, not deleted — since they're still "the
data," just not the series feeding the descriptor CSVs. Do not clean these up
without a deliberate decision.

**Unrelated data-completeness note found during this fix:** 28 of the 34
benchmark substrates (everything outside the original 6-molecule `TEST_IDS`
set) never had Stage 2 (`{mol}_nbo.log`, equilibrium NBO) run at all — only
Stage 1 (`opt`) and Stage 3 (`scan`). This is a pre-existing gap, not
introduced by anything above; `resolve_series()` already treats `nbo` as
optional, so these molecules' series simply start at `scan_1` with no R0
baseline point. Worth flagging for a future batch if the equilibrium point
specifically (not just the stretched series) is ever needed for these 28.

**Automated recovery now exists for this failure family** (see
`CLAUDE.md`'s "Automated failure detection + recovery" section) —
`beckmann/dft/log_diagnostics.py` classifies a log's failure mode
programmatically, and `beckmann/dft/recovery.py` /
`scripts/dft/auto_recover.py` automatically escalate an
`OSCILLATING_DEGENERACY` classification through CalcFC → step=0.07 →
step=0.04, fully automatically (including submission — an explicit product
decision, not an oversight). Verified working end-to-end against Citadel on
2026-07-21: `auto_recover.py --mol 020` and `--mol 003` both correctly
detected the still-oscillating canonical scans and launched fresh
`mol_020_E_calcfc`/`mol_003_E_calcfc` reruns; re-running immediately after
correctly skipped via the in-flight check instead of double-submitting.
Manual diagnosis via this file's playbook is still the right tool for any
failure category the classifier reports as `SLOW_CONVERGENCE`/`NOISY_TRENDING`/
`SEGFAULT`/`UNKNOWN` — those get no automated remediation attempt.
