# Project Context for Claude Code

## Who I am

Teymur Heydarov — researcher at the Institute of Applied Geometry, Johannes Kepler University Linz, Austria. Co-author (with Annalisa Buffa, EPFL, and Bert Jüttler, JKU) of the paper "An Unrefinement Algorithm for Planar THB-spline Parameterizations."

Background: C++ (G+Smo / gismo library), isogeometric analysis (IGA), THB-splines, multipatch parameterizations, spline theory, optimization. This repo (`gismo`) is the library I work within.

---

## The Paper

**"An Unrefinement Algorithm for Planar THB-spline Parameterizations"**
Authors: Heydarov, Buffa, Jüttler. Supported by ERC CHANGE project (GA No. 694515).

**Goal:** Given a fine planar domain parameterization in THB-splines, reduce degrees of freedom (DoF) while keeping the parameterization regular (det J > 0) and geometrically close to the original.

**Three conditions:**
- C1: Regularity — no self-intersections, det J > 0 everywhere
- C2: Global approximation error < ε_global
- C3: Feature boundary error < ε_feature on subset F ⊂ ∂[0,1]²

**Three-stage Projection function:**
1. Local least-squares fitting (cheapest)
2. Linear optimization (LO) — quadratic functionals Q_l (length) + Q_u (uniformity)
3. Non-linear optimization (NLO) — Q_s (skewness) + Q_e (eccentricity) + Q_a (area) via Gauss-Newton

**Key result:** Indiana example: 1423 DoF → 112 DoF using NLO.

**Known weaknesses:** Greedy cell-by-cell (no global optimality), manual NLO weight selection, single-patch/planar only, no computational benchmarks, no simulation-driven validation.

---

## The Codebase

### `fitting_mspline.cpp`
**Role:** Test-case generator — creates distorted multipatch geometries as input for the main algorithm.

- Reads a G+Smo multipatch XML file
- Applies Gaussian bell-curve displacement to interior control points (`applyCompactCenterCompression`) — radial distortion from patch centroid
- Bisection to find max distortion strength keeping det J > 0 and C0 at interfaces
- Enforces C0 by copying interface control points patch1 → patch2
- Per-patch strengths are hardcoded as 1.1× values from a previous logfile
- Patch 7 skipped (no acceptable value found in previous run)
- Output: `*_NLO.xml` written to `filedata/generatedMPs/`

**Known issues:**
- `squareBoundaryIndices` assumes n×n control net, silently falls back for non-square patches
- Patch 7 gets zero distortion — inconsistent test geometry
- Distortion is always radial (centripetal) — does not create skewness or eccentricity

### `poissonTHB_example.cpp`
**Role:** Main algorithm — MPBES unrefinement (multipatch extension of the paper's algorithm). ~21,725 lines.

**Key components:**
- `MPBES` class (lines 143–801): Multi-patch THB-splines with C1/G1 continuity via twin identification (Kraft selection + truncation). Tracks `functionDescription[globalIdx][twinIdx] = {patch, level, tensorIdx}`
- `chooseAdaptiveWeights()` (lines 1482–1541): Automatic NLO weight selection. Normalizes residuals, boosts when Jacobian violations persist (`minusnumber > 0`). Addresses the paper's manual weight problem.
- `selectCellForCoarsening()` method `'g'` (line ~17430): Grenda's geometry-based cell selection. Calls `markCoarseningGeometryBased` with delta ∈ {0.2, 0.4, 0.6, 0.8, 1.0, 1.2}
- `analyzePartitionOfUnity()` (lines 874–917): Validates MPBES basis row sums ≈ 1
- `unrefinementAlgorithmHBJ()`: Main driver, called from `main()` at line 21531
- Default input: `mask_approximation_fine_L3_NLO.xml`

**Pipeline:** `fitting_mspline.cpp` generates distorted XML → `poissonTHB_example` runs MPBES unrefinement on it.

---

## Active Open Problem

**Finding a geometry that requires NLO (not just LO).**

Every test case so far is either trivially handled or fixable by LO alone. We need a geometry where:
1. It is regular at the start (det J > 0 everywhere)
2. It becomes irregular after one cell removal (det J ≤ 0 somewhere)
3. LO alone cannot restore regularity
4. NLO restores regularity

**Why LO fails vs NLO succeeds:**
- LO uses Q_l + Q_u — penalizes *distance* deviations only, blind to angles
- NLO adds Q_s + Q_e + Q_a — penalizes angular deformation (skewness, eccentricity)
- LO cannot fix isoparameter lines crossing at shallow angles or spiral/twist structure

**Proposed approaches (most practical first):**
1. **Swirl/tangential distortion** — add rotational displacement to `fitting_mspline.cpp` instead of radial (`--swirl-alpha` flag). Creates spiral isoparameter lines that LO can't untwist.
2. **Asymmetric focal bell** — center distortion near the cell Grenda selects first, not at patch centroid
3. **S-fold interior** — displacement ∝ sin(πu)·sin(πv)·sin(2π·phase), saddle with two regions pushed opposite directions
4. **High wave number on coarsening boundary** — wave number ≥ degree+1 in coarser space
5. **Anisotropic cells + shear** — long thin cells at 45° to parameter axes

**Validation criterion:** Check `minDet_afterLO < 0` but `minDet_afterNLO > 0` on at least one coarsening step (tracked via `minusnumber` in the summary log).

**Next step agreed on:** Add swirl distortion to `fitting_mspline.cpp` + use XML with level-3 refined cells adjacent to the outer boundary. XML input files are to be provided.
