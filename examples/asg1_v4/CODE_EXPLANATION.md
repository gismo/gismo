# AS-G1 v4 — Code Explanation of the Three Example Files

This document explains, function by function, the three v4
example translation units that implement the Argyris-Speleers G1
(AS-G1) conforming basis for a two-patch planar geometry:

1. `examples/gluing_data_v4.cpp`
2. `examples/2D_embedding_and_D_derivative_and_second_layer_example_v4.cpp`
3. `examples/as_g1_two_patch_basis_v4.cpp`

The header-only repackaging in `examples/asg1_v4/*.h` (see
`USAGE.md`) is a faithful copy of these files inside the
sub-namespace `gismo::asg1v4`; everything said here applies to the
headers too.

References in this text are to the public G+Smo Doxygen
documentation (https://gismo.github.io) and to the paper
**B. Mourrain, A. Vidunas, _Geometrically smooth spline bases for data
fitting and simulation_, CAGD 2016** (the AS-G1 construction).

---

## 0. Notation

Let `G_p : [0,1]² → ℝ²` be the patch-`p` geometry mapping (a B-spline
`gsTensorBSpline` of bidegree `(d,d)` and `C^{d−1}` interior
continuity).  Let `Γ` be the shared interface and `t ∈ [0,1]` a
common tangential parameter (with patch-1's tangent direction as
the reference; if patch-2 runs opposite we call the interface
*flipped*).  Let `n` denote the normal-direction parameter on each
patch.  Define

$$
D_1(t) = \det\bigl(\partial_n G_1,\, \partial_t G_\Gamma\bigr),\quad
D_2(t) = \det\bigl(\partial_n G_2,\, \partial_t G_\Gamma\bigr),\quad
D_3(t) = \det\bigl(\partial_n G_1,\, \partial_n G_2\bigr).
$$

A pair of B-spline trace functions `(f₁, f₂)` is *AS-G1* iff
- (C⁰) `f₁ = f₂` on `Γ`, and
- (G1) there exist polynomials `(α₁, α₂, β₁, β₂)` of low degree
  such that, in physical space, the gradients match:
  $\alpha_2\,\partial_n f_1 + \alpha_1\,\partial_n f_2
   + (\alpha_2\beta_1 - \alpha_1\beta_2)\,\partial_t f = 0$.

The classical Mourrain–Vidunas choice is `αᵢ`, `βᵢ` of degree 1 in
`t` solving the **gluing system**

$$
\alpha_1\,D_2 + \alpha_2\,D_1 = 0, \qquad
\alpha_1\,\beta_2 - \alpha_2\,\beta_1 = -D_3,
$$

with the normalisation $\int_0^1 (\alpha_1+\alpha_2)\,dt = 1$.

---

## 1. `gluing_data_v4.cpp`

### 1.1 Helpers `det2`, `partial`, `breaksOf`

```cpp
template <class T> T det2(const gsVector<T>& a, const gsVector<T>& b);
template <class T> gsVector<T> partial(const gsMatrix<T>& d, index_t dir, index_t col);
template <class T> std::vector<T> breaksOf(const gsGeometry<T>&, short_t dir);
```

`det2` returns `a(0)b(1) − a(1)b(0)`.  `partial` extracts the
`dir`-th derivative column of a `gsGeometry::deriv_into` result
(which stacks Jacobian rows as `[∂x/∂u, ∂y/∂u, ∂x/∂v, ∂y/∂v]`).
`breaksOf` returns the knot breakpoints of a one-direction
component of a tensor B-spline basis; it dynamic-casts the
component to `gsBSplineBasis` and calls `knots().breaks()`
(`gsBSplineBasis::knots` returns a `gsKnotVector`; see Doxygen
class `gismo::gsKnotVector`).

### 1.2 `InterfaceSamples<T>` + `sampleInterface`

```cpp
template <class T> struct InterfaceSamples { gsVector<T> D1, D2, D3, t, w; };
template <class T> InterfaceSamples<T> sampleInterface(
    const gsMultiPatch<T>& mp, const boundaryInterface& interf,
    bool& tangentialFlipped, index_t numGaussPerSpan = 0);
```

For each knot span on patch-1's tangential side, this function
builds Gauss-Legendre nodes & weights via `gsGaussRule<T>` and
`gsGaussRule::mapToAll` (Doxygen: `gismo::gsGaussRule`), then for
each node `t`:

1. Builds 2-D parameter points on both patches (`pt1, pt2`), with
   the patch-2 parameter using `(1 − t)` when the tangential
   orientations are flipped (read from
   `boundaryInterface::dirOrientation`; Doxygen
   `gismo::boundaryInterface`).
2. Calls `gsGeometry::deriv_into` on both patches to get
   `(∂_n G_p, ∂_t G_p)` and `gsGeometry::eval_into` on the
   common tangent curve for `∂_t G_Γ`.
3. Computes `D₁, D₂, D₃` with `det2`.
4. Multiplies through by ±1 so that the normal vectors are
   *outward* on each patch (via `s1, s2` flags from
   `patchSide::parameter`).

The set `tangentialFlipped` is the boolean
`!interf.dirOrientation(ps1, tDir1)` and is reused later in the
embedding.

### 1.3 `SolveResult<T>` + `solveLinearGluing`

```cpp
template <class T> struct SolveResult { T a[4]; T b[4]; T alphaErr, betaErr; };
template <class T> SolveResult<T> solveLinearGluing(const InterfaceSamples<T>&);
```

The four unknowns `α₁,α₂ ∈ span{1-t, t}` (`a[0..3]`) and the
analogous `β` (`b[0..3]`) are recovered by **least-squares
collocation** at the Gauss nodes (which is exact since the
system is linear in the unknowns and degree-matched):

1. The α-system has 4 unknowns and `nPts` equations; we build the
   `nPts × 4` matrix `A_alpha` directly from the basis values
   `1-t, t` weighted by `D₁, D₂`, append a `1 × 4` normalisation
   row (`∫(α₁ + α₂) dt = 1` discretised as `w · (α₁+α₂) = 1`),
   then solve with `Eigen`'s `fullPivLu().solve` (`gsMatrix::inverse`
   returns Eigen-backed objects; see Doxygen `gismo::gsMatrix`).
2. The β-system is solved next using the recovered α: same kind
   of `nPts × 4` linear system with RHS `−D₃ − (α₁β₂ − α₂β₁)`
   linearised around `(β₁, β₂) = 0`.
3. `alphaErr`, `betaErr` are the residual ∞-norms — used by the
   driver to decide PASS/FAIL.

### 1.4 `computeGluingDataForInterface`

```cpp
template <class T> gsVector<T> computeGluingDataForInterface(
    const gsMultiPatch<T>& mp, const boundaryInterface& interf,
    bool& success, T eps = 1e-8, index_t numGaussPerSpan = 0);
```

Top-level: sample, sign-fix, solve, post-process.

- **sign check**: if either `D₁` or `D₂` change sign along `Γ` the
  interface is degenerate (zero Jacobian); we bail with
  `success = false` (the user sees a quiet failure).
- **same-sign branch**: when `sign(D₁) == sign(D₂)` the homogeneous
  α-system `α₁D₂ + α₂D₁ = 0` together with positivity has no
  solution; we flip `D₁ ↦ −D₁`, solve, and remap
  `(α₁, β₂) ↦ (−α₁, −β₂)`.  This matches the v2 reference code.
- **embedding-convention sign flip**: the embedding solver uses
  `β₁D₂ − β₂D₁ = −D₃` while `solveLinearGluing` fits the same with
  `+D₃`; we negate all βs before returning.
- **flipped output ordering**: when the patch-2 tangent runs
  opposite, we swap the two endpoint values of `(α₂, β₂)` so the
  returned 8-tuple is *always* expressed in the gluing-data tangent
  frame `t ∈ [0,1]` from patch-1's point of view.

### 1.5 `computeGluingDataMatrix`

Convenience: loops over `mp.boundaries()` and `mp.interfaces()`
(Doxygen `gismo::gsBoxTopology`), fills a `(nPatches × 16)` matrix
with four numbers per side (W/E/S/N order indexed by
`boxSide::index()`).  Boundary sides keep `(α=1, β=0)`.

### 1.6 Main driver

The driver of `gluing_data_v4.cpp` parses a `gsCmdLine`, reads
geometry with `gsReadFile`, calls `mp.computeTopology()`, and prints
the resulting 16-column row for each patch.  Useful as a smoke test
and for sanity-checking external solvers.

---

## 2. `2D_embedding_and_D_derivative_and_second_layer_example_v4.cpp`

This file is the **per-side embedding**: for one boundary of one
patch, construct the sparse matrix `E` mapping a reduced set of
"effective" DOFs (interior + smoother + lower-degree) back into the
full tensor-B-spline coefficient vector, in a way that fixes the
trace and the gluing-data normal derivative `d_i` on the chosen
side.

### 2.1 1-D helpers `isNested`, `collocationMatrix`, `embeddingMatrix`

```cpp
template<class T> bool isNested(const gsBSplineBasis<T>&, const gsBSplineBasis<T>&);
template<class T> gsSparseMatrix<T> collocationMatrix(const gsBSplineBasis<T>&, const gsMatrix<T>& pts);
template<class T> gsSparseMatrix<T> embeddingMatrix(const gsBSplineBasis<T>& coarse, const gsBSplineBasis<T>& fine);
```

- `isNested` lifts the coarse knot vector by degree elevation to
  match `fine.degree()` and checks `fine.knots().includes(...)`
  (Doxygen `gismo::gsKnotVector::includes`).
- `collocationMatrix(b, pts)` evaluates the basis at `pts` (using
  `gsBasis::eval` and `gsBasis::active`) and fills a sparse matrix
  with row = point, column = active basis index.
- `embeddingMatrix(coarse, fine)` uses the Greville points of the
  fine basis as collocation points (`gsBSplineBasis::anchors`),
  builds `Cc, Cf` collocation matrices, and solves
  `Cf · E = Cc` via `makeSparseLUSolver` (Doxygen
  `gismo::gsSparseSolver`).  `E.sparseView(1, 1e-10)` prunes
  numerical zeros.

### 2.2 `getInteriorDofs`

Set-difference between the union of `boundary` and
`boundaryOffset(side,1)` DOFs and the full DOF list.  The DOFs
identified as *interior* embed via the identity column in `E`.

### 2.3 `createGluingDataArgyrisBasis` — the heart of the file

```cpp
template<class T>
gsSparseMatrix<T> createGluingDataArgyrisBasis(
    const gsTensorBSplineBasis<2,T>& tensorBasis,
    boxSide side,
    T alpha0, T alpha1, T beta0, T beta1,
    T eps = 1e-12, T tangentSign = T(1));
```

**Why three column-blocks?**  The space of patch coefficients
`c ∈ ℝ^{tb.size()}` whose trace lies in `span(side_basis)` is
huge; we want a basis whose:

- *Interior* part is just the unaffected interior coefficients
  (block 1: identity columns).
- *Trace* part (block 3) is parametrised by a **smoother basis**
  (the side basis with `elevateContinuity(1)`, i.e. one more
  continuity at every interior knot) whose values are imposed and
  whose normal derivative `d_i` is forced to **zero**.
- *Normal-derivative* part (block 2) is parametrised by a
  **lower-degree basis** (`degreeReduce(1)` of the side basis)
  whose values on the side are forced to zero but whose `d_i` is
  the lower-degree basis function itself.

Together: the trace lies in the smoother space, the gluing-data
normal derivative lies in the lower-degree space, and the two are
independent.  This is exactly Speleers' construction (see
**B. Mourrain et al., 2016**, sec. 4).

**Implementation steps**

1. Build the smoother & lower-degree side bases; assert nestedness;
   compute the two 1-D embeddings into the full side basis.
2. Identify normal/tangential directions from `boxSide`
   (`side.direction()`, `side.parameter()`), the boundary-row index
   in the patch coefficient grid (`tensorBasis.boundary(side)`,
   `boundaryOffset(side, 1)`), and the stride from one row to the
   next-row neighbour in the global coefficient vector.
3. Compute the boundary and neighbour values of the 1-D
   *normal-direction* basis at the side endpoint via
   `gsBSplineBasis::derivSingle`.  Call these `dBdry, dNeigh`.
   The combination `(d_bdry · u_bdry + d_neigh · u_neigh) / 1` is
   the parametric normal derivative of a B-spline patch at the side.
4. Block 2 (second-layer): for each lower-degree side basis
   function `φⱼ`, the column should produce trace = 0 and
   `d_i = φⱼ`.  Trace = 0 is automatic (boundary row is unused);
   `d_i = φⱼ` becomes a linear system in the neighbour-row
   coefficients which is solved by `collocationMatrix` + LU.
   The scaling `signN / dNeigh` flips for east/north sides
   (`signN = ±1`) and divides out the basis-function derivative.
5. Block 3 (boundary): for each smoother basis function `Vₘ`,
   set the boundary row to `Vₘ`'s coefficients, then solve for
   the neighbour row from the AS-G1 condition rewritten as
   $\gamma \cdot \text{(side basis)} =
    -(d_{bdry}\,V_m + \text{tangentSign}\,\cdot\,\text{sign}_N\,\beta\,\partial_t V_m)\,/\,d_{neigh}.$
   The `tangentSign` parameter is what reconciles a flipped
   interface with the gluing-data tangent frame (v4 fix).
6. Result is assembled as `gsSparseMatrix` with column layout
   `[interior | second-layer | boundary]` and `makeCompressed`.

### 2.4 Main driver

Iterates over all four sides of all patches, calls
`createGluingDataArgyrisBasis`, and exports the columns as VTK PVD
files for visual inspection.  The driver uses
`computeGluingDataForInterface` for shared sides and `(α=1, β=0)`
for boundary sides.

---

## 3. `as_g1_two_patch_basis_v4.cpp`

The end-to-end driver that **assembles a globally conforming AS-G1
basis** for an arbitrary two-patch geometry, verifies it
numerically, and (optionally) plots each global basis function.

### 3.1 Reuse of the helpers

Lines 76–390 of this file are exact copies of the helpers in §1
and §2 above — included so this file can build stand-alone.  When
using the headers in `examples/asg1_v4/*.h` you do not need to
duplicate them.

### 3.2 `main()` — global assembly walk-through

1. **CLI**: `gsCmdLine` with `-f, -r, -o, -p, -n`.  `-p k`
   exports a single basis function; `-p -2` exports all.
2. **Read & sanitize geometry**: `gsReadFile`, `mp.computeTopology()`,
   enforce two patches and one interface, then **degree-elevate**
   to at least 3 (Speleers' AS-G1 needs degree ≥ 3 for the
   smoother/lower-degree bases to be non-trivial), then
   **uniform-refine** with knot multiplicity `deg − 1` (so the
   interior continuity is `C¹` after refinement; this is the
   minimum for the lower-degree basis to be non-trivial).
3. **Solve gluing data**: one call to
   `computeGluingDataForInterface`, abort on failure.
4. **Compute embeddings**: two calls to
   `createGluingDataArgyrisBasis`, passing
   `tangentSign = flipped ? −1 : +1` to *both* patches (v4 fix:
   both halves of the AS-G1 relation are expressed in the
   gluing-data tangent frame, so both per-side embeddings need
   the same sign correction).
5. **DOF accounting**: per patch, count
   `nInt = #interior, nLD = lowerDeg.size, nSm = smoother.size`.
   Assert `nSm1 == nSm2` and `nLD1 == nLD2` (true because the side
   bases coincide after refinement).
6. **Global numbering**:
   `[patch1_int | patch2_int | shared_trace | shared_d-deriv]`.
   The shared blocks are sized `nSm` and `nLD` respectively.
7. **Assemble `G₁`** (size `tb1.size() × nGlobal`):
   - copy `E₁` interior columns → global int1 columns,
   - copy `E₁` smoother (= boundary) columns → global trace columns,
   - copy `E₁` lower-deg (= second-layer) columns → global d-deriv columns.
8. **Assemble `G₂`** (size `tb2.size() × nGlobal`): same, but with
   the v4 fixes:
   - boundary cols use reversed index `j2 = flipped ? (nSm2−1−j) : j`
     so that the same trace coefficient lands at the same physical
     point on both sides;
   - second-layer cols use reversed index AND a sign
     `l2Sign = flipped ? +1 : −1` (aligned: the two patches'
     normals point opposite ways in the gluing-data tangent frame,
     so we negate; flipped: an additional minus from the
     reparametrisation cancels, so we keep +).
9. **Numerical G1 check**: for each global index `k`, build
   `c_p = G_p · e_k`, wrap into a `gsGeometry` via
   `tensorBasis.makeGeometry(c_p)`, sample at 21 points along the
   interface (with patch-2 parameter mirrored when flipped), and:
   - `valErr = |f₁(pt₁) − f₂(pt₂)|` should be 0 by construction;
   - `gradErr = ‖J₁⁻ᵀ ∇f₁ − J₂⁻ᵀ ∇f₂‖` is the **physical-space**
     gradient mismatch (this is the actual definition of G1
     smoothness); for non-bilinear maps the per-DOF error is
     bounded by the cubic-vs.-quartic representation error of `D₃`
     in the side basis and is *not* exactly zero.
   - thresholds: `< 1e-3` is PASS; `< 1e-1` is APPROX; else FAIL.
10. **Plotting**: optional, dumps `.pvd/.vts` via
    `gsParaviewCollection` and `gsWriteParaview` (Doxygen
    `gismo::gsWriteParaview`).

### 3.3 Why three different sign fixes?

The original v3 code passed the orientation sweep only for
aligned interfaces; v4 fixed the flipped case by three separate
sign corrections that *all* depend on `flipped`:

| Where | Quantity | Aligned | Flipped |
| --- | --- | --- | --- |
| Per-side embedding RHS (boundary col) | `tangentSign` | `+1` | `−1` |
| Global G2 second-layer | `l2Sign`     | `−1` | `+1` |
| Global G2 boundary/L2 | column index  | `j`  | `n−1−j` |

The interaction is subtle: an isolated change of any one of the
three breaks every flipped geometry.  Together they restore
PASS on 65/65 orientations in the sweep (`two_patches_general_*.xml`).

### 3.4 What `tangentSign` and `l2Sign` *mean*

The β-term of the gluing-data normal derivative carries one factor
of the tangential direction.  When the patch-2 parameter runs
opposite to the gluing-data tangent, that factor changes sign once
in the per-side embedding (`tangentSign = −1`) **and** once in the
basis-function labelling (`j2 = n−1−j`), netting an extra `−1` on
the second-layer column that must be cancelled in the global
assembly — that is what `l2Sign = +1` (i.e. **don't** negate)
does for flipped interfaces.  In the aligned case neither
correction triggers and we apply the standard AS-G1 sign convention
(`l2Sign = −1`).

---

## 4. Pointers to the G+Smo API

| Used here | Doxygen class / file |
| --- | --- |
| `gsMultiPatch`, `boundaryInterface`, `patchSide`, `boxSide` | `gismo::gsMultiPatch`, `gismo::boundaryInterface` |
| `gsTensorBSplineBasis<2,T>`, `.component(dir)`, `.boundary(side)`, `.boundaryOffset(side,1)`, `.boundaryBasis(side)` | `gismo::gsTensorBSplineBasis` |
| `gsBSplineBasis<T>`, `.elevateContinuity(1)`, `.degreeReduce(1)`, `.anchors`, `.eval`, `.deriv`, `.derivSingle`, `.knots()` | `gismo::gsBSplineBasis` |
| `gsKnotVector::breaks`, `.includes`, `.degreeElevate` | `gismo::gsKnotVector` |
| `gsGaussRule::mapToAll` | `gismo::gsGaussRule` |
| `makeSparseLUSolver`, `gsSparseMatrix::sparseView`, `gsSparseMatrix::makeCompressed`, `InnerIterator` | `gismo::gsSparseSolver`, `gismo::gsSparseMatrix` |
| `gsGeometry::eval_into`, `.deriv_into`, `.makeGeometry` | `gismo::gsGeometry`, `gismo::gsBasis` |
| `gsCmdLine`, `gsReadFile`, `gsInfo`, `gsFileManager::mkdir`, `gsWriteParaview`, `gsParaviewCollection` | `gismo::gsCmdLine`, `gismo::gsReadFile`, `gismo::gsWriteParaview` |

---

## 5. Self-test

The header-only port is verified by `examples/test_headers_v4.cpp`:

```bash
cd build
make test_headers_v4 -j4
for f in ../filedata/domain2d/2patch/two_patches_general_*.xml; do
    ./bin/test_headers_v4 -f ${f#../filedata/} -r 1
done
# Result: PASS on 65/65 orientations.
```

This exercises every orientation permutation
(`Lidentity / LflipU / LflipV / LflipBoth × Ridentity / RflipU /
RflipV / RflipBoth × identity / swapUV`) and confirms the headers
reproduce the original binary's PASS verdict on the full sweep.
