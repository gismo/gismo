# AS-G1 Two-Patch v4 — Header-Only API Usage

The three headers in `examples/asg1_v4/` expose the v4 AS-G1
construction as a reusable, header-only library inside the
sub-namespace `gismo::asg1v4`.  Function names are kept identical
to those in the original `*_v4.cpp` examples; the namespace
prevents linker collisions when an application also links against
those examples.

```cpp
#include "asg1_v4/gluing_data_v4.h"
#include "asg1_v4/argyris_embedding_v4.h"
#include "asg1_v4/as_g1_two_patch_basis_v4.h"
using namespace gismo;
namespace v4 = gismo::asg1v4;   // optional alias
```

All functions are templated on the scalar type `T`; for
double-precision use `T = real_t`.

---

## 1. `gluing_data_v4.h` — solve for `(α₁, α₂, β₁, β₂)`

Computes the AS-G1 gluing data for one shared interface of a
two-patch geometry.  The gluing data are four affine functions
`αᵢ(t)`, `βᵢ(t)` on `t ∈ [0,1]` such that

$$\alpha_1(t)\,D_2(t) + \alpha_2(t)\,D_1(t) = 0, \qquad
  \alpha_1(t)\,\beta_2(t) - \alpha_2(t)\,\beta_1(t) = -D_3(t),$$

where `D₁ = det(∂_n G₁, ∂_t G_Γ)`, `D₂` analogous on patch 2, and
`D₃ = det(∂_n G₁, ∂_n G₂)` are sampled determinants along the
interface (Gauss-quadrature based; see Mourrain–Vidunas 2016).

### Top-level call

```cpp
bool ok = false;
gsVector<real_t> gd = v4::computeGluingDataForInterface<real_t>(
    mp,                 // gsMultiPatch with one interface
    *mp.iBegin(),       // the boundaryInterface
    ok,                 // success flag
    real_t(1e-8),       // residual tolerance
    /*numGaussPerSpan=*/0,
    /*verbose=*/false);
GISMO_ENSURE(ok, "gluing data failed");

// Layout (8 entries):
//   gd(0..1) : alpha_1 = gd(0)*(1-t) + gd(1)*t
//   gd(2..3) : beta_1  = gd(2)*(1-t) + gd(3)*t
//   gd(4..5) : alpha_2 (already swapped if interface is flipped)
//   gd(6..7) : beta_2
```

### Bulk helper for all sides of all patches

```cpp
gsMatrix<real_t> M = v4::computeGluingDataMatrix<real_t>(mp);
// M is (nPatches × 16): four columns per side (W/E/S/N),
// each side stored as (alpha_0, alpha_1, beta_0, beta_1).
// Boundary sides (not shared) are returned as (1,1,0,0).
```

### Lower-level pieces (rarely needed)

| Symbol | Role |
| --- | --- |
| `v4::sampleInterface(mp, ifc, flipped, nG)` | sample `D₁,D₂,D₃,t,w` |
| `v4::solveLinearGluing(samples)`            | fit the four affine functions |
| `v4::InterfaceSamples<T>` / `SolveResult<T>`| return structs |

---

## 2. `argyris_embedding_v4.h` — per-side embedding

Builds, for one side of one patch, the sparse matrix `E`
mapping `[interior | second-layer | boundary]`-DOFs of three
reduced bases into the patch tensor B-spline coefficient vector,
so that the resulting trace satisfies

$$d_i := \tfrac{1}{\alpha_i}\bigl(\partial_n + \beta_i\,\partial_{t_{gd}}\bigr)$$

on the interface (in the **gluing-data tangent frame**).

```cpp
const auto& tb1 = dynamic_cast<const gsTensorBSplineBasis<2,real_t>&>(
                      mp.patch(ps1.patch).basis());

// tangentSign = -1 when the patch tangent runs opposite to the
// gluing-data tangent on this interface (i.e. flipped == true).
const short_t tDir1 = 1 - ps1.direction();
const bool flipped  = !mp.iBegin()->dirOrientation(ps1, tDir1);
const real_t tSign  = flipped ? -1.0 : 1.0;

gsSparseMatrix<real_t> E1 = v4::createGluingDataArgyrisBasis<real_t>(
    tb1, ps1.side(),
    /*alpha0=*/gd(0), /*alpha1=*/gd(1),
    /*beta0=*/ gd(2), /*beta1=*/ gd(3),
    /*eps=*/1e-12, /*tangentSign=*/tSign);
```

Column layout of the returned matrix (per side):
`[ nInterior | nSecondLayer | nBoundary ]`, where

- `nInterior  = tb.size() − boundary − second-layer`
- `nSecondLayer = degreeReduce(side basis).size()`
- `nBoundary    = elevateContinuity(side basis).size()`

Helpers also exported:
`isNested`, `collocationMatrix`, `embeddingMatrix`, `getInteriorDofs`.

---

## 3. `as_g1_two_patch_basis_v4.h` — full driver

End-to-end build of the global AS-G1 conforming basis for a
two-patch geometry, plus an independent G1 smoothness check.

### Build the global basis

```cpp
v4::AsG1Options<real_t> opts;
opts.refinements    = 2;       // uniform refines (knot mult = deg-1)
opts.minDegree      = 3;       // degree-elevate to ≥ 3
opts.numGaussPerSpan = 0;      // 0 = auto (2*deg+1)
opts.eps            = 1e-8;
opts.verbose        = false;

v4::TwoPatchAsG1Basis<real_t> B =
    v4::buildTwoPatchAsG1Basis<real_t>(mp, opts);

// Outputs:
B.mp            // (possibly refined / elevated) multipatch copy
B.ifc           // the shared boundaryInterface
B.flipped       // tangent-flip flag
B.gluingData    // the 8 numbers above
B.E1, B.E2      // per-side embeddings (tb_i.size() × nCols_i)
B.G1, B.G2      // global-to-patch (tb_i.size() × nGlobal)
B.nGlobal,      // total DOFs
B.nInt1, B.nInt2, B.nSharedBdry, B.nSharedL2
B.gOff_int1, B.gOff_int2, B.gOff_bdry, B.gOff_L2  // column offsets
```

Global DOF ordering:
`[ patch1_interior | patch2_interior | shared_trace | shared_d-deriv ]`

To evaluate global basis function `k` on patch `p`:
```cpp
gsVector<real_t> v = gsVector<real_t>::Zero(B.nGlobal);
v(k) = 1.0;
gsVector<real_t> coeffs = (p == 0 ? B.G1 : B.G2) * v;
auto geo = B.mp.patch(p).basis().makeGeometry(coeffs);
gsMatrix<real_t> Y;  geo->eval_into(X, Y);
```

### Verify smoothness numerically

```cpp
v4::G1CheckResult<real_t> R =
    v4::g1SmoothnessCheck<real_t>(B, /*nCheck=*/21);

R.maxValErr;    // max |v1-v2| at sample points (should be ~0)
R.maxGradErr;   // max ||grad_phys v1 - grad_phys v2||
R.maxErrInt;    // by interior DOFs
R.maxErrTrace;  // by trace DOFs
R.maxErrL2;     // by d-derivative DOFs
R.pass;         // valErr<1e-8 AND gradErr<1e-3
```

---

## 4. Complete worked example

```cpp
#include <gismo.h>
#include "asg1_v4/as_g1_two_patch_basis_v4.h"
using namespace gismo;
namespace v4 = gismo::asg1v4;

int main(int argc, char** argv)
{
    std::string fn = "domain2d/2patch/two_bicubic_patches.xml";
    gsCmdLine cmd("AS-G1 v4 demo");
    cmd.addString("f", "file", "Geometry", fn);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsMultiPatch<real_t> mp;
    gsReadFile<real_t>(fn, mp);

    v4::AsG1Options<real_t> opts;  opts.refinements = 1;  opts.verbose = false;
    auto B = v4::buildTwoPatchAsG1Basis<real_t>(mp, opts);
    auto R = v4::g1SmoothnessCheck<real_t>(B);

    gsInfo << "nGlobal=" << B.nGlobal
           << "  gradErr=" << R.maxGradErr
           << "  " << (R.pass ? "PASS" : "FAIL") << "\n";
}
```

Compile it as a regular gismo example (drop it in `examples/` and
re-run CMake; auto-discovery picks it up).

---

## 5. Quick reference — namespace map

| Header | Free function / type | Original .cpp file |
| --- | --- | --- |
| `gluing_data_v4.h`             | `computeGluingDataForInterface` | `gluing_data_v4.cpp` |
| `gluing_data_v4.h`             | `computeGluingDataMatrix`       | `gluing_data_v4.cpp` |
| `gluing_data_v4.h`             | `sampleInterface`, `solveLinearGluing` | `gluing_data_v4.cpp` |
| `argyris_embedding_v4.h`       | `createGluingDataArgyrisBasis`  | `2D_embedding_and_D_derivative_and_second_layer_example_v4.cpp` |
| `argyris_embedding_v4.h`       | `isNested`, `collocationMatrix`, `embeddingMatrix`, `getInteriorDofs` | same |
| `as_g1_two_patch_basis_v4.h`   | `buildTwoPatchAsG1Basis`        | `as_g1_two_patch_basis_v4.cpp` (main()) |
| `as_g1_two_patch_basis_v4.h`   | `g1SmoothnessCheck`             | `as_g1_two_patch_basis_v4.cpp` (verify block) |

All live in `namespace gismo::asg1v4`; nothing leaks to the global
namespace.  A self-test target `bin/test_headers_v4` is provided
(source: `examples/test_headers_v4.cpp`) and confirms 65/65 PASS
over the full `two_patches_general_*.xml` orientation sweep.

## 6. v4 sign-handling notes (do not break)

The header code reproduces the two v4 fixes byte-for-byte:

1. In `computeGluingDataForInterface`, after the same-sign branch
   fix-up, all four `β` coefficients are negated to match the
   embedding's convention `β₁D₂ − β₂D₁ = −D₃`.
2. In `createGluingDataArgyrisBasis`, the boundary-column RHS is
   `−(d_bdry · Vm + tangentSign · sign_N · β · dVm) / d_neigh`;
   pass `tangentSign = −1` when `flipped`.
3. In `buildTwoPatchAsG1Basis`, the second-layer columns of `G₂`
   are scaled by `l2Sign = flipped ? +1 : −1`, and both the
   second-layer and boundary columns are read with reversed index
   `j2 = flipped ? (n−1−j) : j` when `flipped`.

Together, items (1)–(3) make the verifier output ≤ 10⁻⁴ grad-error
on every orientation in the sweep.
