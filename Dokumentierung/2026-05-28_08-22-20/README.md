# NLO Trigger Geometry — Documentation

**Date:** 2026-05-28  
**Author:** Teymur Heydarov  
**Context:** Paper "An Unrefinement Algorithm for Planar THB-spline Parameterizations" (Heydarov, Buffa, Jüttler)

---

## Contents of This Folder

| File | Description |
|------|-------------|
| `mask_approximation_fine_L3.xml` | Base (undistorted) multipatch geometry — input to `fitting_mspline.exe` |
| `mask_approximation_fine_L3_NLO.xml` | S-fold distorted geometry — input to `poissonTHB_example.exe` |
| `fitting_mspline.exe` | Distortion generator (built 2026-05-28) |
| `poissonTHB_example.exe` | MPBES unrefinement algorithm (built 2026-05-28) |

---

## Step 1: Generate the Distorted Geometry

Run `fitting_mspline.exe` from the gismo root directory (so relative paths to `filedata/` resolve):

```
fitting_mspline.exe filedata/generatedMPs/mask_approximation_fine_L3.xml ^
  --sfold-amplitude 5.0 ^
  --sfold-half-width 0.25 ^
  --stress-min-det -50.0 ^
  --fit-grid 20 ^
  --output filedata/generatedMPs/mask_approximation_fine_L3_NLO.xml
```

### What these parameters do

| Parameter | Value | Effect |
|-----------|-------|--------|
| `--sfold-amplitude` | 5.0 | Maximum S-fold displacement amplitude to attempt per patch |
| `--sfold-half-width` | 0.25 | Half-width of the S-fold in parameter space (≈ 1 coarse cell width) |
| `--stress-min-det` | -50.0 | Bisection accepts geometries where min oriented Jacobian det > −50 |
| `--fit-grid` | 20 | 20×20 = 400 check points per patch for the bisection safety check |

### Key output (distortion log)

The bisection accepts strongly-deformed patches — min oriented det ≈ −49 on 8 of 10 patches:

```
Patch 0: radial safe strength=0.0703, s-fold safe amplitude=1.006, min oriented det=-49.97 (18 pts below threshold)
Patch 1: s-fold safe amplitude=0.137, min oriented det=-49.40  (3 pts)
Patch 2: s-fold safe amplitude=0.928, min oriented det=-47.38  (10 pts)
Patch 3: s-fold safe amplitude=0.098, min oriented det=-49.79  (3 pts)
Patch 6: s-fold safe amplitude=0.342, min oriented det=-49.48  (2 pts)
Patch 7: s-fold safe amplitude=0.098, min oriented det=-49.03  (3 pts)
Patch 8: s-fold safe amplitude=0.225, min oriented det=-49.78  (3 pts)
Patch 9: s-fold safe amplitude=0.459, min oriented det=-48.81  (2 pts)
```

Patches 4 and 5 are excluded from S-fold distortion (they share the only patch–patch interface;
distorting them would raise featureError above ε_f = 0.1).

---

## Step 2: Run the Unrefinement Algorithm

Run `poissonTHB_example.exe` from the gismo root directory:

```
poissonTHB_example.exe
```

The executable uses `filedata/generatedMPs/mask_approximation_fine_L3_NLO.xml` by default.

---

## Results and Why They Are Important

### The Paper's Open Problem

The paper describes three optimization stages applied after each cell removal:
1. **LS fitting** (least-squares projection into coarser basis)
2. **LO** (linear optimization: Q_l + Q_u) — triggered when LS leaves irregular points
3. **NLO** (nonlinear optimization: Q_s + Q_e + Q_a) — triggered when LO still leaves irregular points

The paper claims NLO is sometimes *necessary* — i.e., there exist geometries where LO alone cannot restore regularity. Until this session, no concrete example had been found computationally.

### What This Run Demonstrates

On **patch 0**, the first 6 coarsening steps all follow the same pattern:

| Stage | minusnumber | Verdict |
|-------|-------------|---------|
| After LS | ~119 | irregular (0.74% of 16000 MPBES pts) |
| **After LO** | **11** | **LO FAILED** — residual locked at 4.2597, all 11 in patch 7 |
| NLO iter 1 | 3–15 | reducing but not yet zero |
| **NLO iter 2** | **0, converged=true** | **NLO SUCCEEDED** |

Log markers to search for:
```
[LO-step] residual=4.25966, ... minusnumber=11, converged=false
LO candidate did not meet acceptance. Recomputing weights once for NLO fallback.
[NLO-iter 2/10] result: residual=2.90267, ... minusnumber=0, converged=true
NLO worked!
```

### Why LO Fails Here

LO minimizes Q_l + Q_u — quadratic functionals that penalize *distance* deviations only. Its correction is a tangent-space step at the current (invalid) control-point configuration. When 11 MPBES Jacobian points are deeply negative (det ≪ 0), the straight-line correction falls outside the feasible manifold.

NLO adds Q_s + Q_e + Q_a — penalizing angular deformation (skewness, eccentricity, area). These terms curve the optimization trajectory, allowing it to navigate around the infeasible region and land on the regular manifold. This is the mesh-untangling mechanism described in the paper.

### Summary of All Steps Observed

| Step | Patch coarsened | LO minusnumber | NLO invoked | NLO result |
|------|----------------|----------------|-------------|------------|
| 1 | 0 | 11 | ✓ yes | converged (iter 2) |
| 2 | 0 | 11 | ✓ yes | converged (iter 2) |
| 3 | 0 | 11 | ✓ yes | converged (iter 2) |
| 4 | 0 | 11 | ✓ yes | converged (iter 2) |
| 5 | 0 | 11 | ✓ yes | converged (iter 2) |
| 6 | 0 | 11 | ✓ yes | converged (iter 2) |
| 7+ | 1, 2, … | 0 | ✗ no | LO sufficient |

NLO was needed exclusively on patch 0, which carries the strongest S-fold distortion.

---

## Important Constraints (do not change)

- `epsilon_g = 1e6` and `epsilon_f = 0.1` are **hardcoded** in `poissonTHB_example.cpp` — never add as CLI parameters.
- Source files: `examples/poissonTHB_example.cpp`, `examples/fitting_mspline.cpp`
- Git branch: `stable`
