# Agent Change Log

## 2026-05-28 (session 2)

### Scope
Verify that the `--fit-grid 40` geometry is genuinely regular at MPBES resolution,
diagnose why swirl alone does not trigger LO failure, and combine S-fold + swirl to
target angular violations.

### Finding 1: fit-grid 40 geometry confirmed genuinely regular
Running `poissonTHB_example.exe` on the `--fit-grid 40` geometry produced:

```
=== INITIAL MPBES DETERMINANT CHECK ===
[initial-mpbes] irregular points: 0 / 16000 (0%)
```

All patches show 0 initial irregular points. The geometry is **genuinely regular** at MPBES resolution. This closed the previous session's open question.

### Finding 2: S-fold alone → LO always succeeds

Over ~14 tested coarsening steps, `minusnumber(LO)=0` on every step. Root cause:

- S-fold violations appear at the **THB refinement boundary** (u ≈ 0.625 = 5/8 for level-3),
  clustered in patch 0, the same patch being coarsened.
- LO's Q_l + Q_u can reach these control points directly → LO converges in one iteration.
- **S-fold is a compression/accordion distortion** — LO's length penalty fixes it trivially.

### Finding 3: Swirl alone → LS never produces any violations

With `--swirl-alpha 3.0 --bell-sigma 0.25`:
- After LS on every coarsening step: `Total irregular points: 0`.
- Root cause: swirl creates a **smooth, low-frequency rotation**. The coarser THB basis can
  represent a smooth rotation without regularity loss. No violations → LO not even called.

### Proposed fix: S-fold first, then swirl on top
The S-fold creates high-frequency (cell-scale) structure that the coarser basis cannot
represent. The swirl then rotates that structure. Combined effect after coarsening:

- Violations appear (from S-fold's unrepresentable frequency — same mechanism as before).
- But the violations are now at the **patch center** (uv ≈ 0.46, 0.46) instead of the THB
  boundary — because the swirl shifted the violation locations.
- The violations have angular character (Jacobian with large off-diagonal component).
- LO (Q_l + Q_u) may not fix angular violations — this is the test in progress.

### Parameters used for combined S-fold + swirl geometry
```
fitting_mspline.exe mask_approximation_fine_L3.xml \
  --sfold-amplitude 5.0 --sfold-half-width 0.25 \
  --swirl-alpha 3.0 --bell-sigma 0.25 --focal-u 0.375 --focal-v 0.625 \
  --stress-strength 0.0 --stress-min-det 0.001 --fit-grid 40
```

Applied amplitudes:
- S-fold: p0=0.654, p1=0.771, p2=0.537, p3=0.078, p6=0.146, p7=0.068, p8=0.137, p9=0.088
- Swirl: p0=0.352 rad, p1=0.012, p2=1.400, p3=0.035, p4=1.547, p5=1.535, p6=0.334, p7=0.029, p8=0.023, p9=0.059 rad

Min oriented det after both distortions:
- Patch 0: **0.00384** (barely regular!), Patch 2: 0.072, Patch 3: 0.040, Patch 6: **0.029**

Initial MPBES check: **0 irregular points** ✓

Step 1 LS: **2 violations at uv=(0.462, 0.462)** and (0.487, 0.462) — patch center, not THB boundary.
LO outcome: **LO succeeded** (minusnumber(LO)=0 on every step).

### Finding 4: S-fold+swirl on ellipse_hole (4-patch ring domain)
The `ellipse_hole_approximation_fine_L3.xml` (4 patches in a ring around a hole) was tried next because more cross-patch coupling was hoped for (4 interfaces vs 1 for mask).

Distortion:
```
fitting_mspline.exe ellipse_hole_approximation_fine_L3.xml \
  --sfold-amplitude 5.0 --sfold-half-width 0.25 \
  --swirl-alpha 3.0 --bell-sigma 0.25 \
  --stress-strength 0.0 --stress-min-det 0.001 --fit-grid 40
```

Min oriented det: p0=0.033, p1=0.047, p2=0.322, p3=0.178. Initial MPBES: 0 violations.
Steps 1-3: 6, 24, 19 violations after LS (det magnitudes up to -154!) — all in patch 0 (mirrored).
LO: always minusnumber=0 in one iteration.

Total LO outcomes: **0 NLO triggers across all geometries tried.**

### Root Cause (Definitive)
For NLO to trigger, violations must appear in a DIFFERENT PATCH from the coarsened one.
LO can only adjust control points belonging to the coarsened MPBES basis. If violations are
in a remote patch, LO has no access to those control points → fails.

This only occurs when:
1. The initial geometry is highly stressed in a remote patch (min det << 0), AND
2. Coarsening in patch A creates a large enough MPBES rearrangement to push remote patch B negative.

With genuinely regular initial geometry (min det > 0 at all MPBES grid points), the rearrangement
is proportionally smaller — never large enough to push any remote patch below 0.

The -50 case worked because the initial geometry already had pre-existing violations in patch 7
(from the large distortion), which LO could not eliminate.

### `--input` flag added to poissonTHB_example.cpp
To support testing with different input files without recompiling, a `--input <path>` argument
was added to `main()` in `examples/poissonTHB_example.cpp` (near line 21494).

Usage: `poissonTHB_example.exe --input /path/to/file.xml`

### Files Changed
- `examples/poissonTHB_example.cpp` — added `--input` CLI flag (near line 21494)
- `filedata/generatedMPs/mask_approximation_fine_L3_NLO.xml` — overwritten with S-fold+swirl combined geometry.
- `filedata/generatedMPs/ellipse_hole_approximation_fine_L3_NLO.xml` — newly created distorted ellipse-hole.

### Finding 5: TV geometry (16 patches, undistorted) — LO never triggered

Run: `poissonTHB_example.exe --input ".../tv_approximation_fine_L3.xml"` (background task bens4pevv)

- Initial MPBES check: **0 / 25600 irregular points** (genuinely regular ✓)
- **23 coarsening steps completed; all `minusnumber: 0` after LS**
- LS alone was sufficient every step — LO not called, NLO not called
- `Success! iteration = 1, coarselevel = 2` on every step (LS converged in one pass)
- Confirms: undistorted regular geometry → proportionally small rearrangements → no violations → no LO/NLO

### Notes for Next Agent
- **Open question**: Is there ANY regular initial geometry + distortion that triggers NLO from scratch?
- Hypothesis: No, because regularity constraint limits distortion → limits rearrangement magnitude.
- Alternative: Accept the `-50` case as the paper's validation (shows NLO works even if input is irregular).
- Or: Explore geometries where LO is structurally infeasible (non-convex domains forcing impossible interface constraints).
- **Next geometry to try**: `joystick_approximation_fine_L3.xml` (30 patches) — highest cross-patch coupling.
- `epsilon_g=1e6` and `epsilon_f=0.1` remain HARDCODED — never add as CLI parameters.
- The `-50` geometry (valid NLO trigger but irregular initial) is archived in `Dokumentierung/2026-05-28_08-22-20/`.

---

## 2026-05-28

### Scope
Found first geometry that triggers NLO in `poissonTHB_example.cpp` — **the key paper validation milestone**.

### Background
The active open problem was: find a geometry where LO alone cannot restore regularity after a cell removal step, forcing NLO to be invoked and succeed.

### Root Cause of Earlier Failures
`fitting_mspline.cpp`'s bisection always found the maximum distortion keeping `minOrientedDet > tolerance`. With default `--stress-min-det 1e-3`, this produced barely-above-zero det values at the 20×20 check grid, so all MPBES irregular points after LS had det ≈ −0.001 to −0.5. The O(0.01 × Q_u) LO correction easily flipped those to positive.

### Solution
Run `fitting_mspline.exe` with:
```
--sfold-amplitude 5.0 --sfold-half-width 0.25 --stress-min-det -50.0 --fit-grid 20
```
The bisection now accepts geometries with `minOrientedDet > −50` (strongly negative), yielding patches with det ≈ −49 at many check points (18 + 3 + 10 + 3 + 2 + 3 + 3 + 2 points across 8 patches).

### Observed Result (first run)
After the **first coarsening step** (patch 0, level 2):
| Stage | minusnumber | Status |
|-------|-------------|--------|
| After LS | 119 | irregular |
| After LO | **11** | **LO FAILED** (minusnumber > 0) |
| NLO iter 1 | 3 | reducing |
| NLO iter 2 | **0** | **NLO converged ✓** |

This repeats on the **second coarsening step** as well (same LO failure, NLO reducing toward 0).

### Files Changed
- `filedata/generatedMPs/mask_approximation_fine_L3_NLO.xml` — regenerated with the extreme S-fold parameters above.

### How to Reproduce
```
fitting_mspline.exe filedata/generatedMPs/mask_approximation_fine_L3.xml \
  --sfold-amplitude 5.0 --sfold-half-width 0.25 --stress-min-det -50.0 --fit-grid 20 \
  --output filedata/generatedMPs/mask_approximation_fine_L3_NLO.xml
poissonTHB_example.exe   # uses mask_approximation_fine_L3_NLO.xml by default
```

### What to Watch in Log
- `[LO-step] ... minusnumber=11` — LO failed
- `LO candidate did not meet acceptance. Recomputing weights once for NLO fallback.`
- `[NLO-iter 1/10] ... minusnumber=3` (or 8) — NLO reducing
- `[NLO-iter 2/10] ... minusnumber=0, converged=true` — NLO succeeded
- `NLO worked!`

### Notes for Next Agent
- The 11 post-LO irregular points are always in **patch 7** (not the coarsened patch 0), suggesting a global MPBES coupling effect.
- The globalError after NLO ≈ 3.18 (up from 0.19 after LS), well within `epsilon_g=1e6`.
- If a cleaner geometry is needed (fewer initial irregular points in the loaded XML), reduce `--stress-min-det` (less negative) or reduce `--sfold-amplitude`.
- `epsilon_g=1e6` and `epsilon_f=0.1` remain HARDCODED — never add as CLI parameters.
- NLO triggered exclusively on patch 0 (first 6 steps); patches 1+ were clean LO throughout.

### Archived
`Dokumentierung/2026-05-28_08-22-20/` — executables, XMLs, and full README for reproducibility.

---

## 2026-05-27

### Scope
Gate INTERFACE_DIAGNOSTIC output behind a flag in `examples/poissonTHB_example.cpp`.

### Why
The `logSpecificInterface` calls emit ~60 lines per interface pair per coarsening step (5 pairs × many steps), inflating the console log by thousands of lines and making it hard to read progress output.

### Change
- Added `static bool g_enableInterfaceDiagnostics = false;` after the global file handles (near line 815).
- Wrapped all three call-sites of `logSpecificInterface` in `if (g_enableInterfaceDiagnostics) { ... }`:
  - Mesh-export function (formerly lines 4114–4118)
  - `main_postsolve` block (formerly lines 20455–20459)
  - LO/NLO postsolve block (formerly lines 20781–20785)

### How to re-enable
Set `g_enableInterfaceDiagnostics = true;` (one-liner) and rebuild — no other changes needed.

### Files Changed
- `examples/poissonTHB_example.cpp`

---

## 2026-05-07

### Scope
Stabilization and diagnostics updates in `examples/poissonTHB_example.cpp` for the THB unrefinement workflow.

### Why
Runs were repeatedly revisiting the same candidate states, producing loop-like behavior across outer attempt-selection cycles.

### Key Fixes
1. Added safe value-based candidate removal helper:
   - `removeCellIdsByValue(gsVector<int>&, const std::vector<int>&)`
2. Added per-attempt attempted-cell tracking:
   - `attemptedCellIds` stores the exact cell IDs attempted in current iteration.
3. Replaced index-based persistent-pool removal with ID-based removal:
   - Avoids removing the wrong element when `vectorS` and `nonCheckedCells` diverge.
4. Applied withdrawal on failure paths:
   - When LO and NLO both fail, attempted cell IDs are removed from `nonCheckedCells`.
5. Applied withdrawal on emergency no-spline path:
   - Prevents endless retries of non-constructible candidates.
6. Kept existing optimization policy from this session:
   - Single LO attempt, single NLO fallback attempt.
   - NLO early stop when residual increases for 2 consecutive iterations.
   - If both fail, restore `baseVectSol` and withdraw candidate.

### Files Changed
- `examples/poissonTHB_example.cpp`

### Important Locations
- Candidate selection loop starts near level/patch sweep around:
  - `examples/poissonTHB_example.cpp:18230`
- New helper:
  - `examples/poissonTHB_example.cpp:7468`
- Attempted-cell capture:
  - `examples/poissonTHB_example.cpp:18288`
  - `examples/poissonTHB_example.cpp:18400`
  - `examples/poissonTHB_example.cpp:18415`
- Persistent pool updates by value:
  - `examples/poissonTHB_example.cpp:20259`
  - `examples/poissonTHB_example.cpp:20633`
  - `examples/poissonTHB_example.cpp:20658`
  - `examples/poissonTHB_example.cpp:20697`

### Expected Runtime Effect
- `nonCheckedCells` should shrink consistently.
- Same single-cell or geo-cell-set proposals should not be repeatedly selected after they are attempted.
- Outer-loop repetition should reduce significantly.

### Verification Checklist
1. Run `adaptRefinementThb_example.exe`.
2. Confirm repeated patch/level attempt labels are reduced in main log.
3. Confirm `nonCheckedCells.size()` decreases over time.
4. Confirm LO/NLO still follow the configured one-pass policy.

### Notes For Next Agent
- If repetition persists, inspect proposal identity at selection time and log a normalized proposal signature (sorted cell IDs).
- If needed, add a `triedThisLevel` set keyed by signature to strictly prevent duplicates within each level.
