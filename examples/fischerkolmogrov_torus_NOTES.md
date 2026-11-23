# Fisher–Kolmogorov on a torus: body-fitted vs immersed

Work notes (2026-11-19). Goal: solve the same Fisher–Kolmogorov (FK) tumour-growth
problem on a genus-1 geometry (a solid torus) with a **body-fitted** single-patch
NURBS discretization and with the **immersed** (finite-cell, Algoim octree)
discretization. Then compare the two before going back to the prostate.

## 1. Files

| file | what it is |
|---|---|
| `examples/fischerkolmogrov_torus_example_bodyfitted.cpp` | FK on an exact NURBS solid torus (new) |
| `examples/fischerkolmogrov_torus_example_immersed.cpp` | FK on the same torus, immersed in `[0,1]^3` (new) |
| `examples/fischerkolmogrov_prostate_example_bodyfitted.cpp` | prostate body-fitted solver the torus one was copied from (unchanged) |
| `examples/fischerkolmogrov_prostate_example_immersed.cpp` | prostate immersed solver the immersed torus was copied from (unchanged) |

Both new examples are picked up automatically by `examples/CMakeLists.txt`, which globs all `.cpp` files there.
After adding a file, rerun `cmake .` in `build/`, then
`make fischerkolmogrov_torus_example_bodyfitted fischerkolmogrov_torus_example_immersed`.

## 2. Geometry

Default torus, the same in both codes, in box units: centre (0.5, 0.5, 0.5), axis parallel to z,
major radius R = 0.30, tube radius r = 0.15. It spans [0.05, 0.95]² × [0.35, 0.65].
Options `--torusR --torusr --torusCx --torusCy --torusCz` change it.
With `-L 50` (default) the torus has an outer diameter of 45 mm and a tube radius of 7.5 mm, roughly prostate-sized.
Exact volume: 2π²Rr² = 0.13324 box³ = **16 655 mm³**.

**Body-fitted** (`makeTorus()`): an exact single NURBS patch made by revolving
`gsNurbsCreator::NurbsDisk` (cross-section) along `gsNurbsCreator::NurbsCircle`
(toroidal direction):

    P_ijk = c + ((R + r a_ij) cx_k, (R + r a_ij) cy_k, r b_ij),   weight w_ij * v_k

- Parametric directions: u, v = disk cross-section; w = toroidal angle.
- **The seam.** The faces w = 0 and w = 1 coincide. They are glued as a patch interface:
  `mp.addInterface(0, front, 0, back)`, so N is continuous across the seam.
  This was verified: N at both faces is identical, and the DoF count drops by the glued layer.
- **The circle drawn in ParaView** is the seam face. ParaView draws all 6 faces of a structured block.
  It is a visual artifact only. To hide it, apply *Filters → Clean to Grid* on `N.pvd`.
- The disk patch has 4 corner points with a degenerate Jacobian. They are isolated lines along the
  torus, never quadrature points, and show up only as shading artifacts.

**Immersed**: Ω = {φ < 0} with the **exact signed distance** of the torus (class `gsTorusLevelSet`,
analytic value and gradient; no triangle mesh):

    phi(x) = sqrt( (sqrt((x-cx)^2+(y-cy)^2) - R)^2 + (z-cz)^2 ) - r

So both codes discretize *exactly the same domain*. The computed |Ω| is 16 655 mm³, which matches the exact value.

## 3. Parametric vs physical space (important)

- **Body-fitted:** the map is G : [0,1]³ (u,v,w) → torus.
  - Everything is pulled back: ∇ = J⁻ᵀ∇_u, dx = |det J| du.
  - Any function of x must be evaluated at G(u), via `A.getCoeff(f, G)` in gismo expressions.
- **Immersed:** the background map is the **identity** of [0,1]³.
  - Parameter space = physical space (box units), meas = 1, and basis derivatives are already physical.
  - The geometry enters only through φ, i.e. through which quadrature points exist.
- **Both:** 1 box unit = L mm.
  - All user inputs are physical (mm, days) and are converted once, in `fkToBox()`.
  - The csv reports volumes in mm³.
  - ParaView coordinates are box units in **both** codes, so the outputs overlay directly.

## 4. Units and parameters (both codes, identical flags and defaults)

| flag | meaning | unit | default |
|---|---|---|---|
| `-D` | diffusivity | mm²/day | 1e-3 |
| `-p` | proliferation rate ρ | 1/day | 5e-3 |
| `-t`, `-N` | time step, number of steps (T = t·N) | day | 5, 146 (T = 730 d = 2 y) |
| `-R`, `-s` | seed radius, tanh edge width | mm | 5, 1 |
| `-a` | seed peak density (≤ 1; < 1 so the core can densify) | – | 0.5 |
| `--Nth` | threshold that defines the tumour region Ω_T | – | 0.15 |
| `--xc --yc --zc` | seed centre | mm | on the core circle at 45°: (35.6, 35.6, 25) |
| `-L` | mm per box unit | mm | 50 |

Initial condition: N(x,0) = a · 0.5 · (1 − tanh((|x − x_c| − R)/σ)).
The seed must fit in the tube: |offset from the core circle| + R + 2σ ≤ r·L = 7.5 mm.
Both codes warn otherwise.

The xml file (`-f`) is read **only** for `rho_inf` and the Newton `tol`. Its model parameters
(id 3) are in box units, e.g. `radius = 0.15`, and are ignored on purpose.

### Diagnosis of the "tumour shrinks" problem (advisor's email)

Earlier runs used D in **box units**: `-D 5e-5` with 1 box unit = 61.9 mm (the prostate
normalization) means 0.19 mm²/day, which is about 190× the literature value of ~1e-3 mm²/day.
Together with T = 100 that gave a diffusion-dominated regime: the seed dilutes and both V_T and max N fall.
Both torus codes now print, before the time loop:

- **ρT:** growth factor e^(ρT). ρT < 1 means hardly any growth in the horizon.
- **D/R² vs ρ:** if D/R² > ρ, the tumour dilutes instead of growing (warning printed).
- **Front speed 2√(Dρ) and front width ≈ 4√(D/ρ):** the mesh h must resolve the front (warning printed).

With D = 1e-3, ρ = 5e-3, T = 730 d: ρT = 3.65, front ≈ 3.3 mm in 2 years, front width ≈ 1.8 mm.
The tumour expands and gets denser (max N 0.5 → 1.0, V_T ×3), which is the behaviour the advisor expects.
**To confirm with the advisor:** that "1e-3" means mm²/day and 1/day (Lorenzo et al. 2024).

## 5. Bugs found and fixed

1. **The initial condition was evaluated in parametric space** (body-fitted only).
   `gsL2Projection::project(dbasis, mp, f, coefs)` uses `A.getCoeff(f)` without G, so f is
   evaluated at (u,v,w), not at x. The seed came out deformed with the geometry: stretched about 1.9×
   around the torus and squeezed to 0.3× across the tube.
   - **Fix (torus body-fitted):** the example does its own L2 projection with
     `A.getCoeff(N0, G)` in the assembler's own space, which also respects the glued seam.
     Check: the initial V_T is 700 mm³ against the analytic 668 mm³; the old value was 165.
   - **Still present in:** `fischerkolmogrov_prostate_example_bodyfitted.cpp` and
     `fischerkolmogrov_prostate_example.cpp`. These must be fixed before any prostate
     body-fitted vs immersed comparison.
   - Not an issue in the immersed codes: their map is the identity.
   - gismo's `gsL2Projection` itself was left unchanged.
2. **Units:** see section 4.
3. **Update 2026-11-23: the solution space is now C1 everywhere** (quadratic, as in Lorenzo et al. 2024).
   - 90°/180°/270°: `knots(2).reduceMultiplicity(1)` before refinement.
   - Seam (0°): periodic B-splines in w. After refinement, the clamped w-knots are replaced by
     the unclamped uniform vector {−ph, …, 0, h, …, 1, …, 1+ph}, and the first p and last p
     functions are identified with `gsDofMapper::matchDof` and `w.setupMapper()`.
     gismo's own `setPeriodic()` is not used: on this clamped knot vector it loses 2 elements
     and V_P drops by 3%.
   - Check: a seed at the seam (0°) and at 90° (equivalent NURBS positions) give V_T, N_T, A_p
     equal to 12 digits. With the old C0 glued seam, V_T differs by 0.25%.
   - `--C0seam` restores the old glued interface for comparison.
   - `--plot` (or `--plotMesh` alone) also writes `torus_geometry.pvd` (the exact surface,
     300×300 points per side), `mesh_surface.pvd` (element edges on the surface) and
     `mesh_volume.vtp` (all element edges).

   Old text: **C0 planes of the circle.** `NurbsCircle` has double knots at 0.25, 0.5, 0.75, so the inherited
   solution basis was only C0 at 90°, 180° and 270°.
   - The solution basis now drops one copy (`knots(2).reduceMultiplicity(1)`), so N is C1
     everywhere except at the seam (0°), which stays C0.
   - Measured, this was *not* a visible problem: the projected seed matches the exact ball to within
     8% of its peak at 0°, 45° and 180°.
   - The "two halves" look came from the seam disk being drawn through a seed that sat on the seam.
     The default seed is now at 45°.

## 6. Discretization controls

**Body-fitted:**
- The degree is inherited from the geometry (`gsMultiBasis<> dbasis(mp, true)`): NurbsDisk and NurbsCircle
  are quadratic, so **p = 2**. There is no `-e` option yet. To add one, put
  `dbasis.setDegree(p)` right after the `dbasis` line, before the refinement loops.
- Continuity: C^(p−m) at a knot of multiplicity m, so C1 with single knots.

| flag | effect |
|---|---|
| `-r k` | uniform refinement in all directions: 8^k more elements |
| `--refTor t` | extra refinement in the toroidal direction only (elements there are ~2.4× longer) |
| `--plotPts n` | ParaView sample points per frame, spread evenly in parameter space (~n^(1/3) per direction); default 50 000 ≈ 1.6 MB/frame |

Element sizes with L = 50: across the tube h ≈ 15/2^k mm; around the outer rim h ≈ 141/(4·2^(k+t)) mm.

**Immersed:**

| flag | effect |
|---|---|
| `-r m` | background box refinement, h = L/2^m mm, isotropic |
| `-e` | degree elevation of the linear cube; default 1, so **p = 2**, same as body-fitted |
| `--quadDepth`, `--indicator`, `--indicatorTol` | cut-cell octree rule |
| `--plotRes` | N and φ sampled on a uniform grid (one `.vts` per frame); clip with φ < 0 in ParaView |
| `--plotPoints` | additionally dump N at the quadrature points |

Measured cost of the immersed code:

| `-r` | h (mm) | active DoFs | cut cells | quadrature points | cache |
|---|---|---|---|---|---|
| 4 | 3.1 | 1 912 | 664 | 182 k | 192 MB |
| 5 | 1.6 | 8 848 | 2 544 | 850 k | 901 MB |

**Caveat:** the cache of basis values at quadrature points grows about 4–5× per level, so `-r 6`
would need roughly 4+ GB. Nearly all of it comes from cut cells (about 300 points each with
`quadDepth 2`, `integralChange`). If memory becomes the limit, lower `--quadDepth` or `--indicatorTol`,
or stop caching the derivatives.

## 7. How to run

```bash
# body-fitted
./build/bin/fischerkolmogrov_torus_example_bodyfitted \
    -D 1e-3 -p 5e-3 -t 5 -N 146 -R 3 -s 0.75 -a 0.5 \
    -r 4 --refTor 1 --plotPts 200000 -L 50 \
    -o 00_output_FK_torus --plot

# immersed (same physics flags)
./build/bin/fischerkolmogrov_torus_example_immersed \
    -D 1e-3 -p 5e-3 -t 5 -N 146 -R 3 -s 0.75 -a 0.5 \
    -r 5 -L 50 \
    -o 00_output_FK_torus_immersed --plot
```

- Delete an output folder before rerunning (`rm -rf ...`). Otherwise stale frames from older runs stay
  next to the new ones. `N.pvd` only lists the current run.
- In zsh, a variable holding several flags must be expanded with `${=var}` to split into words.

## 8. How to compare the two (plan)

1. **Same physics, matched element sizes, both p = 2:**

   | body-fitted | h across / around (mm) | immersed | h (mm) |
   |---|---|---|---|
   | `-r 3 --refTor 1` | 1.9 / 2.2 | `-r 5` | 1.6 |
   | `-r 4 --refTor 1` | 0.9 / 1.1 | `-r 6` | 0.8 (memory!) |

   Compare at matched h, and above all check that each code **converges under refinement to the same curves**.
2. **Scalar comparison (csv, identical columns):**
   - V_P must equal 16 655 mm³ in both; it does.
   - Then V_T(t), int_N(t), N_dash(t), N_max(t), N_min(t). Plot both sets of curves on one figure.
   - `NumDOFs` is the active DoFs in the immersed code and all basis functions in the body-fitted code.
3. **Field comparison (not implemented yet):**
   - Easiest route: the immersed solution is a plain spline in x. Write it with `gsWrite` at chosen times.
   - Load it in the body-fitted code and integrate ‖N_bf − N_im‖² over the body-fitted quadrature
     with `getCoeff(N_im, G)`, i.e. evaluated at physical points.
   - That gives a relative L2 difference per time step.
4. **Visual:** both outputs are in box units, so open both `.pvd` files in one ParaView session.
   For the immersed one, use *Clip* by φ.

## 9. Open items / next steps

- [ ] Run the full T = 730 d comparison at matched meshes (section 8) and plot the csv curves.
- [ ] Implement the field L2-difference comparison (section 8, item 3).
- [ ] Port the IC fix (`getCoeff(N0, G)`) and physical units to `fischerkolmogrov_prostate_example_bodyfitted.cpp` / `fischerkolmogrov_prostate_example.cpp`.
- [ ] Port physical units (`-L`, `fkToBox`, regime diagnostics, `-a`) to `fischerkolmogrov_prostate_example_immersed.cpp`.
  There, L = 61.9 mm per box unit (fill 0.9) and the input is still in box units.
- [ ] Optional: an `-e` option in the body-fitted torus code; merged `.vtu` output without the seam disk; plot every k-th frame only.
- [ ] Confirm the parameter units (mm²/day, 1/day) with the advisor.
- [ ] The negative undershoot of N (N_min ≈ −0.03 to −0.08) comes from projecting a sharp tanh on a coarse mesh. It shrinks with refinement. Keep σ ≥ h.
