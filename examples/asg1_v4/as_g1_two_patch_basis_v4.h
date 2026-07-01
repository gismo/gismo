/** @file as_g1_two_patch_basis_v4.h

    @brief Header-only AS-G1 conforming two-patch basis assembler (v4).

    Builds a global AS-G1 basis across the shared interface of a
    two-patch geometry by combining `gluing_data_v4.h` and
    `argyris_embedding_v4.h`.  Returns the global-to-patch matrices
    G1, G2 such that the global basis function `phi_k` has coefficients
        c1 = G1.col(k),   c2 = G2.col(k)
    on patch 1 and patch 2 respectively, and satisfies AS-G1 across
    the interface to machine precision (modulo cubic representability
    on curved interfaces).

    Public API (in namespace `gismo::asg1v4`):

      * struct `TwoPatchAsG1Basis<T>`           -- bundled output
      * `buildTwoPatchAsG1Basis(mp, opts)`      -- top-level driver
      * `g1SmoothnessCheck(basis, nCheck)`      -- diagnostic check

    Author(s): F. Hasanova, S. Takacs
*/

#pragma once

#include "gluing_data_v4.h"
#include "argyris_embedding_v4.h"

namespace gismo {
namespace asg1v4 {

/// Options bundle for `buildTwoPatchAsG1Basis`.
template <class T>
struct AsG1Options
{
    index_t refinements    = 0;        ///< uniform refinements before building
    short_t minDegree      = 3;        ///< elevate to at least this degree
    index_t numGaussPerSpan = 0;       ///< 0 = auto (2*deg+1)
    T       eps            = T(1e-8);  ///< gluing-data residual tolerance
    bool    verbose        = true;     ///< print diagnostics
};

/// Output bundle.  All matrices are global-to-patch:
///   G1: (tb1.size() x nGlobal),  G2: (tb2.size() x nGlobal).
template <class T>
struct TwoPatchAsG1Basis
{
    // Inputs (after refinement / elevation)
    gsMultiPatch<T>     mp;            ///< multipatch used (mutable copy)
    boundaryInterface   ifc;           ///< the single shared interface
    bool                flipped = false;

    // Per-patch trace embeddings (interior | second-layer | boundary)
    gsSparseMatrix<T>   E1, E2;

    // Global assembly
    gsSparseMatrix<T>   G1, G2;        ///< global-to-patch matrices
    index_t             nGlobal = 0;
    index_t             nInt1 = 0,  nInt2 = 0;
    index_t             nSharedBdry = 0, nSharedL2 = 0;

    // Convenience global-column offsets
    index_t             gOff_int1 = 0, gOff_int2 = 0;
    index_t             gOff_bdry = 0, gOff_L2   = 0;

    // Gluing data
    gsVector<T>         gluingData;    ///< 8 numbers: a1_0..b2_1
};

/// Build the AS-G1 conforming two-patch basis end-to-end.
template <class T>
TwoPatchAsG1Basis<T> buildTwoPatchAsG1Basis(
    const gsMultiPatch<T>& mpIn,
    const AsG1Options<T>& opts = AsG1Options<T>())
{
    TwoPatchAsG1Basis<T> out;
    out.mp = mpIn;
    gsMultiPatch<T>& mp = out.mp;
    mp.computeTopology();

    GISMO_ENSURE(mp.nPatches() == 2,
                 "buildTwoPatchAsG1Basis: expected exactly two patches.");
    GISMO_ENSURE(mp.nInterfaces() == 1,
                 "buildTwoPatchAsG1Basis: expected exactly one interface.");

    // Degree elevation
    const short_t inputDeg = mp.patch(0).basis().degree(0);
    if (inputDeg < opts.minDegree)
        mp.degreeElevate(opts.minDegree - inputDeg);

    // Refinement with knot multiplicity = deg-1
    if (opts.refinements > 0)
    {
        const short_t deg = mp.patch(0).basis().degree(0);
        const index_t mult = std::max<index_t>(deg - 1, 1);
        for (index_t i = 0; i < opts.refinements; ++i)
            mp.uniformRefine(1, mult);
    }

    out.ifc = *mp.iBegin();
    const patchSide ps1 = out.ifc.first(), ps2 = out.ifc.second();

    // Gluing data
    bool ok = false;
    out.gluingData = computeGluingDataForInterface<T>(
        mp, out.ifc, ok, opts.eps, opts.numGaussPerSpan, opts.verbose);
    GISMO_ENSURE(ok, "buildTwoPatchAsG1Basis: gluing data failed.");

    const T a1_0 = out.gluingData(0), a1_1 = out.gluingData(1);
    const T b1_0 = out.gluingData(2), b1_1 = out.gluingData(3);
    const T a2_0 = out.gluingData(4), a2_1 = out.gluingData(5);
    const T b2_0 = out.gluingData(6), b2_1 = out.gluingData(7);

    const short_t tdir1 = 1 - ps1.direction();
    out.flipped = !out.ifc.dirOrientation(ps1, tdir1);
    const T tSign = out.flipped ? T(-1) : T(1);

    const gsTensorBSplineBasis<2,T>& tb1 =
        dynamic_cast<const gsTensorBSplineBasis<2,T>&>(mp.patch(ps1.patch).basis());
    const gsTensorBSplineBasis<2,T>& tb2 =
        dynamic_cast<const gsTensorBSplineBasis<2,T>&>(mp.patch(ps2.patch).basis());

    out.E1 = createGluingDataArgyrisBasis<T>(
        tb1, ps1.side(), a1_0, a1_1, b1_0, b1_1, T(1e-12), tSign);
    out.E2 = createGluingDataArgyrisBasis<T>(
        tb2, ps2.side(), a2_0, a2_1, b2_0, b2_1, T(1e-12), tSign);

    // DOF counts
    gsBSplineBasis<T> sideBasis1 = *tb1.boundaryBasis(ps1.side());
    gsBSplineBasis<T> sideBasis2 = *tb2.boundaryBasis(ps2.side());
    gsBSplineBasis<T> smB1 = sideBasis1; smB1.elevateContinuity(1);
    gsBSplineBasis<T> ldB1 = sideBasis1; ldB1.degreeReduce(1);
    gsBSplineBasis<T> smB2 = sideBasis2; smB2.elevateContinuity(1);
    gsBSplineBasis<T> ldB2 = sideBasis2; ldB2.degreeReduce(1);

    gsMatrix<index_t> bdry1  = tb1.boundary(ps1.side());
    gsMatrix<index_t> neigh1 = tb1.boundaryOffset(ps1.side(), 1);
    std::vector<index_t> intDOFs1 = getInteriorDofs(tb1.size(), bdry1, neigh1);

    gsMatrix<index_t> bdry2  = tb2.boundary(ps2.side());
    gsMatrix<index_t> neigh2 = tb2.boundaryOffset(ps2.side(), 1);
    std::vector<index_t> intDOFs2 = getInteriorDofs(tb2.size(), bdry2, neigh2);

    out.nInt1 = static_cast<index_t>(intDOFs1.size());
    out.nInt2 = static_cast<index_t>(intDOFs2.size());
    const index_t nLD1 = ldB1.size(), nLD2 = ldB2.size();
    const index_t nSm1 = smB1.size(), nSm2 = smB2.size();

    GISMO_ENSURE(nSm1 == nSm2, "Smoother basis sizes differ across interface.");
    GISMO_ENSURE(nLD1 == nLD2, "Lower-degree basis sizes differ across interface.");
    out.nSharedBdry = nSm1;
    out.nSharedL2   = nLD1;

    out.gOff_int1 = 0;
    out.gOff_int2 = out.nInt1;
    out.gOff_bdry = out.nInt1 + out.nInt2;
    out.gOff_L2   = out.nInt1 + out.nInt2 + out.nSharedBdry;
    out.nGlobal   = out.gOff_L2 + out.nSharedL2;

    // Assemble G1
    out.G1 = gsSparseMatrix<T>(tb1.size(), out.nGlobal);
    for (index_t j = 0; j < out.nInt1; ++j)
        for (typename gsSparseMatrix<T>::InnerIterator it(out.E1, j); it; ++it)
            out.G1.insert(it.row(), out.gOff_int1 + j) = it.value();
    for (index_t j = 0; j < nLD1; ++j)
    {
        const index_t e1col = out.nInt1 + j;
        for (typename gsSparseMatrix<T>::InnerIterator it(out.E1, e1col); it; ++it)
            out.G1.insert(it.row(), out.gOff_L2 + j) = it.value();
    }
    for (index_t j = 0; j < nSm1; ++j)
    {
        const index_t e1col = out.nInt1 + nLD1 + j;
        for (typename gsSparseMatrix<T>::InnerIterator it(out.E1, e1col); it; ++it)
            out.G1.insert(it.row(), out.gOff_bdry + j) = it.value();
    }
    out.G1.makeCompressed();

    // Assemble G2 (with the two v4 conditional signs)
    out.G2 = gsSparseMatrix<T>(tb2.size(), out.nGlobal);
    for (index_t j = 0; j < out.nInt2; ++j)
        for (typename gsSparseMatrix<T>::InnerIterator it(out.E2, j); it; ++it)
            out.G2.insert(it.row(), out.gOff_int2 + j) = it.value();

    const T l2Sign = out.flipped ? T(1) : T(-1);
    for (index_t j = 0; j < nLD2; ++j)
    {
        const index_t j2 = out.flipped ? (nLD2 - 1 - j) : j;
        const index_t e2col = out.nInt2 + j2;
        for (typename gsSparseMatrix<T>::InnerIterator it(out.E2, e2col); it; ++it)
            out.G2.insert(it.row(), out.gOff_L2 + j) = l2Sign * it.value();
    }
    for (index_t j = 0; j < nSm2; ++j)
    {
        const index_t j2 = out.flipped ? (nSm2 - 1 - j) : j;
        const index_t e2col = out.nInt2 + nLD2 + j2;
        for (typename gsSparseMatrix<T>::InnerIterator it(out.E2, e2col); it; ++it)
            out.G2.insert(it.row(), out.gOff_bdry + j) = it.value();
    }
    out.G2.makeCompressed();

    return out;
}

/// Diagnostic result of `g1SmoothnessCheck`.
template <class T>
struct G1CheckResult
{
    T maxValErr    = T(0);   ///< max |v1-v2| at sample points
    T maxGradErr   = T(0);   ///< max ||grad v1 - grad v2|| in physical space
    T maxErrInt    = T(0);   ///< max grad err over interior-DOF basis funcs
    T maxErrTrace  = T(0);   ///< max grad err over trace-DOF basis funcs
    T maxErrL2     = T(0);   ///< max grad err over d-deriv-DOF basis funcs
    bool pass      = false;  ///< grad < 1e-3 and val < 1e-8
};

/// Verify AS-G1 of every global basis function at `nCheck` points along
/// the interface, in physical (gradient) space.
template <class T>
G1CheckResult<T> g1SmoothnessCheck(
    const TwoPatchAsG1Basis<T>& B,
    index_t nCheck = 21)
{
    G1CheckResult<T> R;
    const patchSide ps1 = B.ifc.first(), ps2 = B.ifc.second();
    const short_t normDir1 = ps1.direction(), tanDir1_ = 1 - normDir1;
    const short_t normDir2 = ps2.direction(), tanDir2_ = 1 - normDir2;
    const bool par1_ = ps1.parameter(), par2_ = ps2.parameter();

    const gsTensorBSplineBasis<2,T>& tb1 =
        dynamic_cast<const gsTensorBSplineBasis<2,T>&>(B.mp.patch(ps1.patch).basis());
    const gsTensorBSplineBasis<2,T>& tb2 =
        dynamic_cast<const gsTensorBSplineBasis<2,T>&>(B.mp.patch(ps2.patch).basis());

    gsMatrix<T> sup1 = B.mp.patch(ps1.patch).support();
    gsMatrix<T> sup2 = B.mp.patch(ps2.patch).support();
    const T ifcCoord1 = sup1(normDir1, par1_ ? 1 : 0);
    const T ifcCoord2 = sup2(normDir2, par2_ ? 1 : 0);
    const T t1a = sup1(tanDir1_, 0), t1b = sup1(tanDir1_, 1);
    const T t2a = sup2(tanDir2_, 0), t2b = sup2(tanDir2_, 1);

    for (index_t idx = 0; idx < B.nGlobal; ++idx)
    {
        gsVector<T> v = gsVector<T>::Zero(B.nGlobal);
        v(idx) = T(1);
        gsVector<T> c1 = B.G1 * v;
        gsVector<T> c2 = B.G2 * v;
        auto f1 = tb1.makeGeometry(c1);
        auto f2 = tb2.makeGeometry(c2);

        T thisMaxGrad = T(0);
        for (index_t i = 0; i < nCheck; ++i)
        {
            T s = T(i) / T(nCheck - 1);
            T t1 = t1a + s * (t1b - t1a);
            T s2v = B.flipped ? (T(1) - s) : s;
            T t2 = t2a + s2v * (t2b - t2a);

            gsMatrix<T> pt1(2, 1), pt2(2, 1);
            pt1(normDir1, 0) = ifcCoord1; pt1(tanDir1_, 0) = t1;
            pt2(normDir2, 0) = ifcCoord2; pt2(tanDir2_, 0) = t2;

            R.maxValErr = std::max(R.maxValErr,
                                   std::abs(f1->eval(pt1)(0,0) - f2->eval(pt2)(0,0)));

            gsMatrix<T> df1, df2, dG1m, dG2m;
            f1->deriv_into(pt1, df1);
            f2->deriv_into(pt2, df2);
            B.mp.patch(ps1.patch).deriv_into(pt1, dG1m);
            B.mp.patch(ps2.patch).deriv_into(pt2, dG2m);

            gsMatrix<T> J1(2,2), J2(2,2);
            J1(0,0) = dG1m(0,0); J1(0,1) = dG1m(1,0);
            J1(1,0) = dG1m(2,0); J1(1,1) = dG1m(3,0);
            J2(0,0) = dG2m(0,0); J2(0,1) = dG2m(1,0);
            J2(1,0) = dG2m(2,0); J2(1,1) = dG2m(3,0);

            gsVector<T> pg1(2), pg2(2);
            pg1(0) = df1(0,0); pg1(1) = df1(1,0);
            pg2(0) = df2(0,0); pg2(1) = df2(1,0);

            gsVector<T> phys1 = J1.inverse().transpose() * pg1;
            gsVector<T> phys2 = J2.inverse().transpose() * pg2;
            T ge = (phys1 - phys2).norm();
            R.maxGradErr = std::max(R.maxGradErr, ge);
            thisMaxGrad = std::max(thisMaxGrad, ge);
        }
        if (idx < B.gOff_bdry)      R.maxErrInt   = std::max(R.maxErrInt,   thisMaxGrad);
        else if (idx < B.gOff_L2)   R.maxErrTrace = std::max(R.maxErrTrace, thisMaxGrad);
        else                        R.maxErrL2    = std::max(R.maxErrL2,    thisMaxGrad);
    }
    R.pass = (R.maxValErr < T(1e-8)) && (R.maxGradErr < T(1e-3));
    return R;
}

} // namespace asg1v4
} // namespace gismo
