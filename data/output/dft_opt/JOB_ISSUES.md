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
