# What to measure: body-fitted vs immersed FK on the torus (CMAME study)

**Claim to support:** the immersed discretization gives the same answer as
the body-fitted one, at comparable cost, without constructing a body-fitted
mesh. Each measurement below supports one part of that claim.

Codes:

- `examples/fischerkolmogrov_torus_example_bodyfitted.cpp`: exact single-patch
  NURBS torus. Quadratic solution space, **C1 everywhere** (since 2026-11-23):
  - single knots at 90/180/270°;
  - periodic B-splines across the seam at 0°. `--C0seam` restores the old
    glued C0 seam.
- `examples/fischerkolmogrov_torus_example_immersed.cpp`: exact signed-distance
  torus in [0,1]³, with Algoim quadrature on cut cells. p = 2, C1 everywhere.

Biomarker definitions follow Lorenzo et al. (2024), Eqs. 5–10; see
`fischerkolmogrov_torus_BIOMARKERS.md`. Status: ✅ implemented, ⚠️ partly or
runs needed, ❌ not implemented.

---

## 0. Output files (both codes, identical columns)

| file | content |
|---|---|
| `biomarkers.csv` | one row per time step, **t = 0 included**: `TimeStep,Time,NumDOFs,int_N,int_H,V_P,N_dash,N_T,A_p,N_max,N_min,V_exact,V_P_relerr,V_T,V_N,growth,dintN_dt`. Volumes in mm³, time in days, 12 digits. `V_T` = `int_H`; `growth` = ∫_Ω ρN(1−N) [mm³/day]; `dintN_dt` = (∫N^n − ∫N^(n−1))/Δt |
| `computational.csv` | one row per step: `TimeStep,Time,NewtonIts,LinIts,ResidualNorm,StepSolve_s,StepTotal_s,Wall_s`. `LinIts` = CG iterations summed over the Newton iterations of the step; `StepSolve_s` = Newton loop (assembly + linear solves); `StepTotal_s` includes biomarkers and output |
| `computational_summary.csv` | `quantity,value`: method, mesh, elements (interior/cut), **h_max/h_min/h_mean [mm]** (body-fitted also h_cross_max, h_tor_max), DoFs, active DoFs, quadrature points (interior/cut/total), quadrature cache [MB], quadrature build time, IC projection time, **IC_L2_relerr**, IC_Nmax/Nmin, V_P_relerr, dt, steps, Newton its (total, mean), **lin_solver, cg_tol, lin_its_total, lin_its_per_newton_mean**, setup/solve/loop/total times, and `cond_*` with `--cond` |

The log ends with a `Times:` line.

With `--plot`, the body-fitted code writes:
- `N.pvd`: N on a [0,1]³ grid of the geometry, one frame per step;
- `torus_geometry.pvd`: the exact surface;
- `mesh_surface.pvd` and `mesh_volume.vtp`: the computational mesh.

**Element size h:** for each element, its longest edge (chord between the mapped
corners), in mm.

| body-fitted | elements | h_max | h_min | h_mean | immersed | h |
|---|---|---|---|---|---|---|
| `-r 2 --refTor 1` | 4×4×32 | 4.64 | 3.14 | 3.97 | `-r 4` | 3.125 |
| `-r 3 --refTor 1` | 8×8×64 | 2.33 | 1.41 | 1.94 | `-r 5` | 1.56 |
| `-r 4 --refTor 1` | 16×16×128 | 1.16 | 0.68 | 0.95 | `-r 6` | 0.78 |

---

## 1. Linear solver and conditioning

### 1a. Solver option (both codes) ✅
`--solver cgdiag|cgic|direct`, **default `cgdiag`**. The helper is taken from
`adaptive_phase-field_fracture_example_DC_*.cpp`.

- `cgdiag`: Eigen CG with a diagonal (Jacobi) preconditioner.
- `cgic`: CG with incomplete Cholesky.
- `direct`: SimplicialLDLT (Pardiso if available). This was the old behaviour,
  so use `--solver direct` to go back to it.
- Options `--cgTol` (default 1e-8, relative to ‖rhs‖) and `--cgMaxIt`
  (default 20000). A CG solve that hits the iteration cap stops the run with
  an error.
- Body-fitted: the IC L² projection also uses the selected solver. Immersed:
  the IC projection uses gismo's `gsL2Projection` (see 5.4).

Check, 3 steps of Δt = 1 day (paper parameters, seed R = 2.5, σ = 1). The
biomarkers agree with `direct` to about 9 digits in every case.

| run | direct [s/step] | cgdiag [s/step] | CG its/step |
|---|---|---|---|
| body-fitted r3 (6.4 k DoFs) | 0.44 | 0.18 | 145 |
| body-fitted r4 (41 k DoFs) | 19.5 | **1.2** (16× faster) | 119 |
| immersed r4 (1.9 k active) | 0.54 | 0.16 | 225 |
| immersed r4, `cgic` | – | 0.32 | 185 |

The direct solver's cost grows superlinearly with the DoFs (3D fill-in). The
CG iteration count is roughly mesh-independent (see 1b), so CG is what makes
body-fitted r5 and immersed r6 affordable.

### 1b. Condition number (`--cond`) ✅
At the first Newton matrix, the code estimates κ = λ_max/λ_min:
- λ_max by power iteration;
- λ_min by inverse power iteration with an LDLT factorization.

It does this for K and for the Jacobi-scaled S = D^(−1/2) K D^(−1/2), which is
the matrix that CG with the diagonal preconditioner effectively sees. The
results go to the log and to `cond_*` in the summary. For the immersed code
they are also given for the **active** block (`cond_active_*`): the full matrix
contains the inactive background DoFs, which are only held by a 1e-10
regularization.

`--cond` costs one factorization plus up to 2 × 1000 iterations (in step 1).
Use it on a separate short run (`-N 1`), not in the production runs.

| run | κ(K) | κ(S), Jacobi-scaled |
|---|---|---|
| body-fitted r3 | 6.5e3 | 458 |
| body-fitted r4 | 4.9e3 | 416 |
| immersed r4, full (incl. inactive) | 9.8e12 | 1188 |
| immersed r4, **active DoFs** | **2.7e11** | **1178** |

- **Body-fitted:** κ is mesh-independent (r3 ≈ r4). The matrix is mass-dominated
  (D Δt/h² ≈ 0.001–0.006), and κ(S) ≈ 400–450 is the known conditioning of a
  3D quadratic B-spline mass matrix. Hence about 120–145 CG iterations.
- **Immersed:** κ(K) of the active block is 2.7e11. This is the small-cut-cell
  problem: basis functions with tiny support inside Ω. Diagonal scaling removes
  almost all of it, leaving κ(S) ≈ 1.2e3, only about 3× the body-fitted value.
  This is why `cgdiag` works for the immersed code too.
- **For the paper:** report κ(K) and κ(S) against h for both methods (run
  `--cond -N 1` at each level). Optional: κ against the cut position (F2).

---

## A. Geometry and initial data (t = 0): is the input the same?

### A1. Domain volume error ✅
e_V = |V_P − V_exact| / V_exact, with V_exact = 2π²Rr²L³ = **16 654.957 mm³**
(R = 15 mm, r = 7.5 mm). csv columns `V_exact` and `V_P_relerr`.

| method | mesh | `--quadDepth` | quad. points | e_V |
|---|---|---|---|---|
| immersed | `-r 2` | 2 | 9.5 k | 4.0e-5 |
| immersed | `-r 3` | 2 | 46 k | 4.3e-6 |
| immersed | `-r 4` | 0 | 38 k | 4.3e-6 |
| immersed | `-r 4` | 1, 2, 3 (identical) | 182 k | 3.1e-8 |
| immersed | `-r 5` | 0 | 200 k | 3.1e-8 |
| immersed | `-r 5` | 2 (default) | 850 k | 4.3e-10 |
| body-fitted | `-r 3 --refTor 1` | – | – | 4.0e-10 |

- The geometric error is far below the biomarker differences.
- `--quadDepth` ≥ 1 gives identical rules, because the `integralChange`
  indicator never subdivides further.
- `--quadDepth 0` halves the memory at r5 (1.3 GB instead of 2.4 GB), and the
  biomarkers change by about 1e-6.

### A2. Initial-condition projection error ✅
`IC_L2_relerr` = ‖N_h(0) − N₀‖_{L²(Ω)} / ‖N₀‖_{L²(Ω)}. It is unit-free and
uses the same definition in both codes.

Examples (R = 2.5, σ = 1):
- body-fitted r2: 0.43 (r2 cannot represent the seed; use it only as the
  coarsest convergence point);
- body-fitted r3 (R = 3, σ = 0.75): 0.11;
- immersed r4 (R = 3, σ = 0.75): 0.26.

### A3. Over- and undershoot of the initial condition ✅
`IC_Nmax` and `IC_Nmin`, against the exact peak (0.497 for a = 0.5, R = 2.5,
σ = 1). Rule of thumb: h ≤ σ.

| mesh | max N(0) | min N(0) |
|---|---|---|
| body-fitted `-r 2 --refTor 1` | 0.342 | −0.080 |
| body-fitted `-r 3 --refTor 1` | 0.566 | −0.022 |
| body-fitted `-r 4 --refTor 1` | 0.491 | −0.001 |
| immersed `-r 4` | 0.501 | −0.054 |
| immersed `-r 5` | 0.540 | −0.011 |

### A4. Initial tumour volume ✅ (t = 0 row of `biomarkers.csv`)
Ω_T(0) is a ball of radius R + σ·artanh(1 − 2N_th/a):

- R = 2.5, σ = 1: **V_T(0) = 104.7 mm³**.
- R = 3, σ = 0.75: **V_T(0) = 153.0 mm³**.

---

## B. The paper's biomarkers: are the clinical outputs the same?

### B1. Time histories ✅ (columns there, runs needed)
Plot V_T, V_N, N̄ (`N_dash`), N_T and A_p for both methods, plus a lower panel
with the relative difference to the reference solution.
- V_T and A_p are the sensitive ones.
- N̄ and N_T are averages and agree more closely.

### B2. Values at the final time ✅
These are the last row of `biomarkers.csv`. Post-processing: the relative
error at t = T, and max_t |q_h(t) − q_ref(t)| / |q_ref(t)|, per mesh level.

### B3. A_p normalization (open)
The code follows Eq. 10 as written (mm³/day). The paper reports A_p in 1/day,
which suggests it divides by V_T, so report `A_p / V_T` when comparing with the
paper's values. For the method comparison either form works.

---

## C. Local agreement: is the tumour in the same place with the same shape?

The paper's local metrics are the **DSC**, the **RMSE**, and the local PCC and
CCC between N̂ maps. We use the first two, in continuous (integral) form instead
of the paper's voxel sums.

### C1. Dice similarity coefficient ❌
DSC(a, b) = 2 |Ω_T^a ∩ Ω_T^b| / (|Ω_T^a| + |Ω_T^b|), with Ω_T = {N̂ ≥ 0.15}:

  DSC = 2 ∫_Ω H_a H_b dΩ / (∫_Ω H_a dΩ + ∫_Ω H_b dΩ),  H_x = H(N_x − 0.15).

DSC = 1 means identical tumour regions. For reference, the paper's model vs
MRI values are about 0.8. Use the same smoothed H as for V_T (k = 200), or the
sharp step.

### C2. RMSE and relative L² difference ❌
  RMSE(a, b) = sqrt( (1/V_P) ∫_Ω (N_a − N_b)² dΩ )   (N̂ is unit-free, so RMSE is too)
  e_L2 = ‖N_a − N_b‖_{L²(Ω)} / ‖N_b‖_{L²(Ω)}

The RMSE is the continuous version of the paper's voxel RMSE (paper: about
0.02, model vs MRI). Optional: the H¹ seminorm.

### How to compute C1 and C2 (plan, not yet implemented)
Everything is one quadrature loop over the **body-fitted** mesh (exact
geometry), accumulating ∫(N_a − N_b)², ∫N_b², ∫H_aH_b, ∫H_a, ∫H_b and V_P.

1. `--saveAt t1,t2,...` in both codes: write the solution at the chosen times
   with `gsWrite`.
   - Immersed: the background spline, a function of physical x (box units),
     because the background map is the identity.
   - Body-fitted: the solution spline in parameter space (u,v,w).
2. `--compareWith <file>` in the body-fitted code: at each quadrature point
   u, with x = G(u), evaluate
   - N_a = the own (body-fitted) solution at u;
   - N_b = a body-fitted file at u (same geometry map, so the same u), **or**
     an immersed file at x (`getCoeff(N_im, G)`).
3. Write DSC, RMSE and e_L2 per saved time to `comparison.csv`.

Uses: body-fitted vs body-fitted (convergence, reference = finest
body-fitted) and immersed vs body-fitted (method comparison). Effort: about
half a day.

### C3. Peak density and bounds ✅
`N_max`, `N_min` per step.

---

## D. Convergence: do both methods converge to the same solution?

### D1. Error against DoFs / h ⚠️ (runs needed)
- Body-fitted: `-r 2, 3, 4` with `--refTor 1`, and r5 as the reference (now
  affordable with `cgdiag`).
- Immersed: `-r 4, 5, 6` (r6 with `--quadDepth 0`).
- Plot, log-log against DoFs or `h_max_mm`, with both methods on one plot: the
  biomarker errors (B2), RMSE and e_L2 (C2), and 1 − DSC (C1).
- Expected for p = 2 and smooth N: L² ≈ h³.

### D2. Matched-resolution pairs
| body-fitted | h_mean [mm] | DoFs | immersed | h [mm] | active DoFs |
|---|---|---|---|---|---|
| `-r 3 --refTor 1` | 1.94 | 6 400 | `-r 5` | 1.56 | 8 848 |
| `-r 4 --refTor 1` | 0.95 | 41 k | `-r 6` | 0.78 | (pending) |

---

## E. Consistency checks: is the physics right?

### E1. Mass balance ✅ (columns `growth`, `dintN_dt`)
With zero flux, d/dt ∫N = ∫_Ω ρN(1−N). Compare `dintN_dt` with the trapezoidal
mean of `growth` over the step.

**Finding:** in both codes, step 1 has `dintN_dt` about 20% below `growth`;
from step 2 on they agree to within the time-discretization error. The cause
is the generalized-α start: dN/dt(0) is set to 0 instead of being computed from
the PDE. It is the same in both codes, so the comparison stays fair, but see
5.2.

### E2. Time-step independence ⚠️ (run included)
`bf_r3_dt5` (Δt = 5) against `bf_r3` (Δt = 1, the paper's value).

---

## F. Cost and robustness

### F1. Cost ✅
Everything is in `computational_summary.csv` and `computational.csv` (see 0).
Peak memory still comes from the scheduler (`sacct`) or `/usr/bin/time -l`.
Plot the error against DoFs and against wall-clock time.

Measured with `direct`, Δt = 1, 730 steps:
- body-fitted r2: Newton 43% of the loop;
- body-fitted r3: Newton 68%;
- the rest is biomarkers and output (see 5.3).

### F2. Small-cut robustness (optional) ❌ (no code needed)
Shift the torus by fractions of h (`--torusCx/Cy/Cz`) and show that the
biomarkers and κ(S) (`--cond -N 1`) vary much less than the discretization
error.

---

## Priorities

| priority | items |
|---|---|
| **Required** | A1, A2/A3, B1–B2, C1 (DSC), C2 (RMSE / L²), D1, F1 |
| **Strongly recommended** | E1, E2, 1b (κ against h, both methods) |
| **Optional** | F2, H¹ error, local PCC/CCC as in the paper |

---

## 5. What is still missing in the code

1. **DSC / RMSE / L² field comparison** (C1, C2): `--saveAt` and
   `--compareWith`, see the plan in C. This is the main gap.
2. **Consistent initial rate:** compute dN/dt(0) from M·dN₀ = −D K N₀ +
   ∫ρN₀(1−N₀)φ (one extra linear solve) instead of 0. This removes the step-1
   mass-balance error (E1). Same change in both codes.
3. **Body-fitted biomarker cost:** each step calls 6 separate `computeInt`, and
   V_P never changes. A single pass, like the immersed code does, saves most
   of the ~30% non-Newton time.
4. **The IC projection differs between the codes:** body-fitted projects onto
   the space over Ω, immersed uses `gsL2Projection` over the whole background
   box. Both are L²-accurate, but they are different operators near the
   boundary. For a strictly identical setup, project over Ω (the cut-cell
   quadrature) in the immersed code too.
5. **Output size:** `--plot` writes every step, 731 frames at Δt = 1. A
   `--plotEvery k` option would reduce that.
6. **Immersed r6:** memory and time not yet measured (first job on the server).

## 6. Run script

`run_fk_torus_cmame.sh` (SLURM): two lanes in parallel (body-fitted and
immersed), light → heavy, paper parameters (D = 1.26e-3, ρ = 1.86e-3, Δt = 1,
T = 730 days, R = 2.5, σ = 1, a = 0.5). All runs use the default `cgdiag`;
body-fitted r5 can now be added at the end of the BF lane.
