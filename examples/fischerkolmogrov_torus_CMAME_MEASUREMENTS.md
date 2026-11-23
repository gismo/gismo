# What to measure: body-fitted vs immersed FK on the torus (CMAME study)

**Claim to support:** the immersed discretization gives the same answer as
the body-fitted one, at comparable cost, without constructing a body-fitted
mesh. Each measurement below supports one part of that claim.

Codes:

- `examples/fischerkolmogrov_torus_example_bodyfitted.cpp`: exact single-patch
  NURBS torus. p = 2, C1 except at the glued seam (0°).
- `examples/fischerkolmogrov_torus_example_immersed.cpp`: exact signed-distance
  torus in [0,1]³, with Algoim quadrature on cut cells. p = 2, C1 everywhere.

Biomarker definitions follow Lorenzo et al. (2024), Eqs. 5–10; see
`fischerkolmogrov_torus_BIOMARKERS.md`. Status: ✅ in the csv/log, ⚠️ partly
there, ❌ not implemented.

---

## A. Geometry and initial data (t = 0): is the input the same?

### A1. Domain volume error ✅
e_V = |V_P − V_exact| / V_exact, with V_exact = 2π²Rr²L³.

With R = 15 mm, r = 7.5 mm (R = 0.30 and r = 0.15 box units, L = 50 mm),
V_exact = **16 654.957 mm³**.

- Body-fitted: the NURBS geometry is exact, so e_V is only the quadrature
  error of the rational Jacobian, about 1e-10.
- Immersed: e_V measures the geometric integration error of the cut-cell
  quadrature and should decrease with h.
- csv columns `V_exact` and `V_P_relerr` (12 digits). The immersed log also
  prints `Volume: relative error ...` right after building the quadrature.

Measured (seed R = 2.5 mm, σ = 1 mm):

| method | mesh | `--quadDepth` | quad. points | e_V |
|---|---|---|---|---|
| immersed | `-r 2` (h = 12.5 mm) | 2 | 9.5 k | 4.0e-5 |
| immersed | `-r 3` (h = 6.25 mm) | 2 | 46 k | 4.3e-6 |
| immersed | `-r 4` (h = 3.1 mm) | 0 | 38 k | 4.3e-6 |
| immersed | `-r 4` | 1, 2, 3 (identical) | 182 k | 3.1e-8 |
| immersed | `-r 5` (h = 1.6 mm) | 0 | 200 k | 3.1e-8 |
| immersed | `-r 5` | 2 (default) | 850 k | 4.3e-10 |
| body-fitted | `-r 3 --refTor 1` | – | – | 4.0e-10 |

Conclusions:

- The geometric error is many orders of magnitude below the differences in the
  biomarkers (about 1–5 % in V_T), so the geometry representation does not
  affect the comparison.
- `--quadDepth` ≥ 1 makes no difference, because the `integralChange`
  indicator never subdivides further. `--quadDepth 0` reduces memory by about
  2× (at r5: 1.3 GB instead of 2.4 GB) and changes the biomarkers by only
  ~1e-6 (relative).

### A2. Initial-condition projection error ⚠️
‖N_h(0) − N₀‖_{L²(Ω)} / ‖N₀‖_{L²(Ω)}, with N₀ the exact tanh seed.

- Both codes print `IC L2-projection error`. The body-fitted value is in
  physical units, the immersed one in box units, and neither is relative.
  ➜ Normalize both to a relative error in mm before comparing.

### A3. Over- and undershoot of the initial condition ✅ (log)
max N_h(0) and min N_h(0), against the exact peak (0.497 for a = 0.5,
R = 2.5, σ = 1). L² projection of a sharp tanh produces Gibbs oscillations, so
this shows whether the mesh resolves the seed. Rule of thumb: h ≤ σ.

Measured (R = 2.5 mm, σ = 1 mm):

| mesh | h [mm] | max N(0) | min N(0) |
|---|---|---|---|
| body-fitted `-r 3 --refTor 1` | 1.9 / 2.2 | 0.566 | −0.022 |
| body-fitted `-r 4 --refTor 1` | 0.9 / 1.1 | 0.491 | −0.001 |
| immersed `-r 4` | 3.1 | 0.501 | −0.054 |
| immersed `-r 5` | 1.6 | 0.540 | −0.011 |

### A4. Initial tumour volume and cell volume ✅
V_T(0) and V_N(0) against the exact seed values. For N₀ = a·½(1 − tanh((|x−x_c| − R)/σ)),
Ω_T(0) is a ball of radius R + σ·artanh(1 − 2N_th/a):

- R = 2.5 mm, σ = 1 mm, a = 0.5: radius 2.924 mm, so V_T(0) = **104.7 mm³**.
- R = 3 mm, σ = 0.75 mm, a = 0.5: radius 3.318 mm, so V_T(0) = **153.0 mm³**.

➜ To do: write a t = 0 row to the csv. Currently the first row is t = Δt.

---

## B. The paper's biomarkers: are the clinical outputs the same?

These are what a clinician would read off, so they are the headline results.

### B1. Time histories ⚠️
V_T(t), V_N(t), N̄(t), N_T(t), A_p(t) for both methods on one figure, plus a
lower panel with the relative difference to the reference solution.

- V_T and A_p are the sensitive ones: V_T depends on where the N̂ = 0.15
  contour lies, and A_p is the marker the paper uses for risk classification.
- N̄ and N_T are averages and will agree more closely.
- ➜ V_N is computed but **not written to the csv**: add a `V_N` column.

### B2. Values at the final time ✅
Relative error of each biomarker at t = T and the maximum over time,
max_t |q_h(t) − q_ref(t)| / |q_ref(t)|. Present as a table per mesh level.

### B3. A_p normalization (open)
The code follows Eq. 10 as written (mm³/day). The paper reports A_p in 1/day,
which suggests it divides by V_T. For a comparison between methods either form
works. Report A_p / V_T as well if the values are to be compared with the
paper.

---

## C. Local agreement: is the tumour in the same place with the same shape?

Agreeing volumes do not guarantee the same tumour shape.

### C1. Dice similarity coefficient of Ω_T ❌
DSC = 2|Ω_T^a ∩ Ω_T^b| / (|Ω_T^a| + |Ω_T^b|), between each solution and the
reference, at selected times (e.g. 0, T/2, T). The paper uses DSC as its local
metric, so it is a familiar measure for this application.

### C2. Field difference ❌
‖N_a − N_ref‖_{L²(Ω)} / ‖N_ref‖_{L²(Ω)}, and the H¹ seminorm, at selected
times.

- Route: the immersed solution is a plain spline in x (identity map). Write it
  with `gsWrite` at the chosen times.
- Load it in the body-fitted code and integrate over the body-fitted
  quadrature with `getCoeff(N_im, G)`, i.e. evaluated at physical points.
- The same loop gives DSC (C1) with H(N_a − 0.15)·H(N_b − 0.15).

### C3. Peak density and bounds ✅
N_max(t), N_min(t): peak density and violations of [0, 1] (undershoot).

---

## D. Convergence: do both methods converge to the same solution?

This is the central figure of the study. A single matched-mesh comparison is
not enough.

### D1. Error against DoFs ⚠️ (runs needed)
- Refine both methods:
  - body-fitted: `-r 2, 3, 4` with `--refTor 1`;
  - immersed: `-r 4, 5, 6` (r6 with `--quadDepth 0` for memory).
- Plot, log-log against DoFs (active DoFs for immersed) and with both methods
  on the same axes, the error of every B1 quantity and the C2 L² field error.
- Reference: the finest body-fitted run, or a Richardson extrapolation from
  the three finest levels.
- Expected for p = 2 and smooth N: L² rate ≈ h³ (DoFs^(−1)). A reduced rate in
  the immersed curves would point to cut-cell issues.
- What the figure should show: both curves tend to the same limit, with
  similar error per DoF.

### D2. Matched-resolution pairs
| body-fitted | h across / around [mm] | DoFs | immersed | h [mm] | active DoFs |
|---|---|---|---|---|---|
| `-r 3 --refTor 1` | 1.9 / 2.2 | 6 600 | `-r 5` | 1.6 | 8 848 |
| `-r 4 --refTor 1` | 0.9 / 1.1 | 42 120 | `-r 6` | 0.8 | (r6 test pending) |

---

## E. Consistency checks: is the physics right?

### E1. Mass balance ❌
With zero flux, d/dt ∫_Ω N̂ dΩ = ∫_Ω ρN̂(1−N̂) dΩ. Check, per step,

  ∫N̂(t_{n+1}) − ∫N̂(t_n)  ≈  Δt · ∫_Ω ρN̂(1−N̂) dΩ   (at the generalized-α
  intermediate state).

- In the immersed code this tests that the natural zero-flux condition holds
  on cut cells.
- ➜ Add a `growth = ∫_Ω ρN̂(1−N̂)` column. This is the integral over **Ω**,
  unlike A_p, which is over Ω_T.

### E2. Time-step independence ⚠️ (runs needed)
Run one mesh with Δt = 5 days and Δt = 1 day (the paper's value). Show the
time error is much smaller than the spatial error. This justifies Δt = 5 days
for the rest of the study.

---

## F. Cost and robustness: what does it cost?

### F1. Cost ⚠️
Per run, report:

- DoFs (active DoFs for immersed);
- number of quadrature points (interior and cut, for immersed);
- setup time: quadrature construction, IC projection, factorization;
- time per step and Newton iterations per step;
- peak memory.

Plot error against wall-clock time next to error against DoFs.

- Immersed: quadrature time and cache size are logged. Peak memory and total
  time currently come from `/usr/bin/time -l`.
- ➜ Add timing columns to the csv.

### F2. Small-cut robustness (optional, often asked by reviewers) ❌
- Shift the torus by fractions of h relative to the background grid
  (`--torusCx`, `--torusCy`, `--torusCz`) and show the biomarkers vary much
  less than the discretization error.
- Report the condition number of the system, or the iterative-solver iteration
  count, against the cut configuration.

---

## Priorities

| priority | items |
|---|---|
| **Required** | A1, A2/A3, B1 (V_T, N_T, A_p), C1 (DSC), C2 (L² difference), D1 (convergence against DoFs), F1 (cost) |
| **Strongly recommended** (cheap, heads off reviewer questions) | E1 (mass balance), E2 (Δt check), A4 (t = 0 row) |
| **Optional** | F2 (small-cut robustness), H¹ field error |

## Code changes needed

1. csv: add `V_N`, a t = 0 row, `growth` (E1), and timing columns
   (setup, time per step, Newton iterations).
2. Normalize the IC projection error identically in both codes (A2).
3. Field comparison: write the immersed solution at selected times; evaluate
   the L² difference and DSC in the body-fitted code (C1, C2).
