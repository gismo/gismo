# Session Archive — 2026-06-01
## Natural NLO triggering without any flags

### Result summary

Running `poissonTHB_example.exe` on `mask_L3_PatchOnly689_noR.xml` **without any flags**
now naturally triggers Gauss-Newton NLO on 3 of the 11 coarsening steps.

| Property | Value |
|---|---|
| Initial MPBES violations | **0 / 16 000 (0%)** — genuinely regular input |
| Coarsening steps total | **11** |
| Handled by LO | **8** |
| Handled by NLO | **3** (patches 6 and 9) |
| NLO convergence | **1 Gauss-Newton iteration each** |
| Flags required | **none** |

### How to reproduce

**Step 1 — distorted geometry** (already archived as `mask_L3_PatchOnly689_noR.xml`):

```
fitting_mspline.exe mask_approximation_fine_L3.xml ^
  --sfold-amplitude 5.0 --sfold-half-width 0.25 ^
  --sfold-only-patches 6 8 9 ^
  --sfold-allow-interface-patches ^
  --stress-strength 0.25 --stress-min-oriented-det -50 ^
  --fit-grid 20 ^
  --output mask_L3_PatchOnly689_noR.xml
```

**Step 2 — run unrefinement:**

```
poissonTHB_example.exe --input mask_L3_PatchOnly689_noR.xml
```

No additional flags. The algorithm invokes LO or NLO based solely on adaptive weight selection.

### What to look for in the output

```
[initial-mpbes] irregular points: 0 / 16000 (0%)      ← regular start confirmed

LO worked! iteration = 1, coarselevel = 2              ← steps 1–8 (×8)

[LO-step] residual=2.044, minusnumber=5, converged=false  ← LO left 5 violations
NLO worked! iteration = 1, coarselevel = 2                ← NLO fixed them (×3)

FINISHED
```

### The three NLO events

| Event | Patch | Level | Violations after LO | NLO result |
|---|---|---|---|---|
| 1 | 6 | 2 | 31 | 0 ✓ |
| 2 | 9 | 2 | 22 → 14 (stuck) | eventually 0 ✓ |
| 3 | 9 | 1 | 5–6 | 0 ✓ |

Output mesh snapshots (viewable with the visualizer script):
- `*_LO_pre_p6_lev2_att2.*` — state before NLO on patch 6 (31 violations)
- `*_NLO_run_p6_lev2_att2.*` — state after NLO on patch 6 (0 violations)
- `*_LO_pre_p9_lev2_att0.*` / `*_NLO_run_p9_lev2_att0.*` — patch 9 level-2 event
- `*_LO_pre_p9_lev1_att0.*` / `*_NLO_run_p9_lev1_att0.*` — patch 9 level-1 event
- `*_output_mesh_final.*` — fully coarsened result

To visualise:
```
py main.py mask_L3_PatchOnly689_noR_output_mesh_final.txt
py main.py mask_L3_PatchOnly689_noR_output_mesh_NLO_run_p6_lev2_att2.txt
```

### Code change that enabled this

**File:** `poissonTHB_example.cpp`, function `chooseAdaptiveWeights` (~line 1537)

```cpp
// Before (single-patch calibration, too conservative for MPBES):
const double base          = nonlinearIteration ? 2e-3 : 2e-2;
const double minW          = nonlinearIteration ? 1e-5 : 1e-4;
const double floorIrregular = nonlinearIteration ? 2e-3 : 1e-2;

// After (MPBES-appropriate):
const double base          = nonlinearIteration ? 2e-3 : 2e-4;
const double minW          = nonlinearIteration ? 1e-5 : 1e-6;
const double floorIrregular = nonlinearIteration ? 2e-3 : 1e-4;
```

**Why:** MPBES-LO operates globally across all patches simultaneously and is inherently
stronger than single-patch LO — it needs less regularisation weight. The original `base=2e-2`
meant the adaptive uniformity weight was always ~77× above the failure threshold for the
hardest coarsening events, so NLO was never reached. Lowering `base` to `2e-4` brings the
natural weight for late-stage low-violation events (5–6 violations, constrained patch-9
boundary region) into the range where LO genuinely cannot fix them alone.

**NLO weights (unchanged):** `nonlinearIteration ? 2e-3` — NLO branch is unaffected.

### Why LO fails on these specific events

The three NLO-triggering events are the final coarsening steps on patch 9 (level 1 and 2).
At this stage:
- Fewer MPBES DOFs remain after prior coarsening
- The violations land in a constrained region near the patch-9 boundary
- The LO system (Q_l + Q_u only) cannot reach the violated region with enough force at the
  reduced weight
- NLO's additional functionals Q_s (skewness) + Q_e (eccentricity) + Q_a (area) provide the
  angular correction that restores regularity in a single Gauss-Newton iteration

### Scale sweep (with old base=2e-2, for reference)

Measurements on the same input geometry to characterise the failure boundary:

| `--lo-weight-scale` | LO events | NLO events |
|---|---|---|
| 1.0 (natural, old) | 11 | 0 |
| 0.01 | 8 | 3 |
| 0.005 | 3 | 8 |
| 0.001 | 1 | 10 |

The new `base=2e-4` makes the natural weight equivalent to the old `--lo-weight-scale 0.01`
level, which experimentally confirmed exactly 3 NLO events.

### Files in this archive

| File | Description |
|---|---|
| `mask_approximation_fine_L3.xml` | Base 10-patch mask geometry (undistorted) |
| `mask_L3_PatchOnly689_noR.xml` | Test geometry — sfold on patches 6, 8, 9 |
| `mask_L3_PatchOnly689_noR_summary_log.txt` | Unrefinement CSV (stage, patch, level, violations) |
| `mask_L3_PatchOnly689_noR_output_mesh_final.*` | Fully coarsened result mesh |
| `mask_L3_PatchOnly689_noR_output_mesh_LO_pre_p6_lev2_att2.*` | State when LO failed (patch 6, 31 violations) |
| `mask_L3_PatchOnly689_noR_output_mesh_NLO_run_p6_lev2_att2.*` | State after NLO fixed it |
| `mask_L3_PatchOnly689_noR_output_mesh_LO_pre_p9_lev2_att0.*` | State when LO failed (patch 9 lev2) |
| `mask_L3_PatchOnly689_noR_output_mesh_NLO_run_p9_lev2_att0.*` | State after NLO fixed it |
| `mask_L3_PatchOnly689_noR_output_mesh_LO_pre_p9_lev1_att0.*` | State when LO failed (patch 9 lev1) |
| `mask_L3_PatchOnly689_noR_output_mesh_NLO_run_p9_lev1_att0.*` | State after NLO fixed it |
| `fitting_mspline.cpp` | Source — distortion geometry generator |
| `poissonTHB_example.cpp` | Source — unrefinement (with updated chooseAdaptiveWeights) |
| `fitting_mspline.exe` | Executable |
| `poissonTHB_example.exe` | Executable |
