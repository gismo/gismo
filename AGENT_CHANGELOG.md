# Agent Change Log

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
