/** @file argyris_embedding_v4.h

    @brief Header-only per-side Argyris-type embedding with gluing-data
           directional derivative (v4).

    Wraps the `createGluingDataArgyrisBasis` family into the namespace
    `gismo::asg1v4` so that it may coexist with the example translation
    units (which define the same function names in the global namespace).

    Public API (in namespace `gismo::asg1v4`):

      * `isNested`, `collocationMatrix`, `embeddingMatrix`
                                                  -- 1D B-spline utilities
      * `getInteriorDofs(tbSize, bdry, second)`   -- DOF classifier
      * `createGluingDataArgyrisBasis(
            tensorBasis, side, a0, a1, b0, b1, eps, tangentSign)`
                                                  -- per-side embedding

    Column layout of `createGluingDataArgyrisBasis` result:
        [interior DOFs | second-layer cols | boundary cols]

    Pass `tangentSign = -1` when the patch's natural tangent runs
    opposite to the gluing-data tangent (`!ifc.dirOrientation(...)`);
    otherwise +1.  For boundary sides (alpha=1, beta=0), tangentSign
    is irrelevant.

    Author(s): F. Hasanova, S. Takacs
*/

#pragma once

#include <iostream>
#include <set>
#include <map>
#include <algorithm>
#include <numeric>
#include <gismo.h>

namespace gismo {
namespace asg1v4 {

template<typename T>
inline bool isNested(const gsBSplineBasis<T>& coarseBasis,
                     const gsBSplineBasis<T>& fineBasis)
{
    if (coarseBasis.dim() != 1 || fineBasis.dim() != 1) return false;
    const int coarseDegree = coarseBasis.degree();
    const int fineDegree   = fineBasis.degree();
    if (coarseDegree > fineDegree) return false;
    const gsKnotVector<T>& coarseKnots = coarseBasis.knots();
    const gsKnotVector<T>& fineKnots   = fineBasis.knots();
    gsKnotVector<T> elevated = coarseKnots;
    elevated.degreeElevate(fineDegree - coarseDegree);
    return fineKnots.includes(elevated);
}

template<typename T>
inline gsSparseMatrix<T> collocationMatrix(const gsBSplineBasis<T>& basis,
                                           const gsMatrix<T>& pts)
{
    gsMatrix<T>       vals = basis.eval(pts);
    gsMatrix<index_t> idx  = basis.active(pts);
    gsSparseEntries<T> entries;
    entries.reserve(vals.rows() * vals.cols());
    for (int i = 0; i < vals.cols(); ++i)
        for (int j = 0; j < vals.rows(); ++j)
            entries.add(i, idx(j, i), vals(j, i));
    gsSparseMatrix<T> C(pts.cols(), basis.size());
    C.setFrom(entries);
    return C;
}

template<typename T>
inline gsSparseMatrix<T> embeddingMatrix(const gsBSplineBasis<T>& coarse,
                                         const gsBSplineBasis<T>& fine)
{
    gsMatrix<T> greville = fine.anchors();
    gsSparseMatrix<T> Cc = collocationMatrix(coarse, greville);
    gsSparseMatrix<T> Cf = collocationMatrix(fine,   greville);
    gsMatrix<T> result;
    makeSparseLUSolver(Cf)->apply(Cc, result);
    return result.sparseView(1, 1e-10);
}

inline std::vector<index_t> getInteriorDofs(
    const index_t tbSize,
    const gsMatrix<index_t>& firstLayerDOFs,
    const gsMatrix<index_t>& secondLayerDOFs)
{
    std::vector<index_t> removedSet;
    removedSet.reserve(firstLayerDOFs.rows() + secondLayerDOFs.rows());
    std::merge(firstLayerDOFs.data(),
               firstLayerDOFs.data() + firstLayerDOFs.rows(),
               secondLayerDOFs.data(),
               secondLayerDOFs.data() + secondLayerDOFs.rows(),
               std::back_inserter(removedSet));
    std::sort(removedSet.begin(), removedSet.end());

    std::vector<index_t> all(tbSize);
    std::iota(all.begin(), all.end(), 0);

    std::vector<index_t> interior;
    interior.reserve(all.size() - removedSet.size());
    std::set_difference(all.begin(), all.end(),
                        removedSet.begin(), removedSet.end(),
                        std::back_inserter(interior));
    return interior;
}

/// Per-side Argyris-type embedding with gluing-data directional derivative.
///
///   d_i = (1/alpha_i) * (partial_n + beta_i * partial_t_gd)
///
/// Column layout of the result:  [interior | second-layer | boundary]
///
/// Pass `tangentSign = -1` when the patch tangent runs opposite to the
/// gluing-data tangent.
template<typename T>
gsSparseMatrix<T> createGluingDataArgyrisBasis(
    const gsTensorBSplineBasis<2,T>& tensorBasis,
    boxSide side,
    T alpha0, T alpha1,
    T beta0,  T beta1,
    T eps = T(1e-12),
    T tangentSign = T(1))
{
    gsBSplineBasis<T> sideBasis = *tensorBasis.boundaryBasis(side);

    GISMO_ENSURE(sideBasis.knots().minInteriorMultiplicity() > 1,
                 "Too small interior multiplicity.");

    gsBSplineBasis<T> sideSmootherBasis = sideBasis;
    sideSmootherBasis.elevateContinuity(1);
    GISMO_ASSERT(isNested(sideSmootherBasis, sideBasis),
                 "Computed bases not nested.");

    gsBSplineBasis<T> sideLowerDegreeBasis = sideBasis;
    sideLowerDegreeBasis.degreeReduce(1);
    GISMO_ASSERT(isNested(sideLowerDegreeBasis, sideBasis),
                 "Computed bases not nested.");

    gsSparseMatrix<T> embFirstLayer  = embeddingMatrix(sideSmootherBasis,     sideBasis);
    gsSparseMatrix<T> embSecondLayer = embeddingMatrix(sideLowerDegreeBasis,  sideBasis);

    const bool    isLow  = !side.parameter();
    const int     dir    = side.direction();
    const int     tanDir = 1 - dir;
    const index_t stride       = (dir == 0) ? 1 : tensorBasis.size(0);
    const index_t signedStride = isLow ? stride : -stride;

    const gsBSplineBasis<T> normalBasis = tensorBasis.component(dir);
    const index_t nNormal  = normalBasis.size();
    const index_t bdryIdx  = isLow ? 0 : nNormal - 1;
    const index_t neighIdx = isLow ? 1 : nNormal - 2;

    gsMatrix<T> bdryPt1D = normalBasis.support().col(isLow ? 0 : 1);
    const T dBdry  = normalBasis.derivSingle(bdryIdx,  bdryPt1D)(0, 0);
    const T dNeigh = normalBasis.derivSingle(neighIdx, bdryPt1D)(0, 0);
    GISMO_ENSURE(std::abs(dNeigh) > eps, "Neighbor derivative is zero.");

    const T signN = isLow ? T(-1) : T(1);

    gsMatrix<index_t> firstLayerDOFs  = tensorBasis.boundary(side);
    gsMatrix<index_t> secondLayerDOFs = tensorBasis.boundaryOffset(side, 1);
    std::vector<index_t> interiorDOFs =
        getInteriorDofs(tensorBasis.size(), firstLayerDOFs, secondLayerDOFs);

    const index_t nSide     = sideBasis.size();
    const index_t nSmoother = sideSmootherBasis.size();
    const index_t nLowerDeg = sideLowerDegreeBasis.size();
    const index_t nInterior = static_cast<index_t>(interiorDOFs.size());

    gsMatrix<T> greville = sideBasis.anchors();
    const index_t nPts = greville.cols();

    gsMatrix<T> supTan = tensorBasis.component(tanDir).support();
    const T tA = supTan(0, 0), tB = supTan(0, 1);

    gsMatrix<T> sideVals  = sideBasis.eval(greville);
    gsMatrix<T> sideDeriv = sideBasis.deriv(greville);
    gsMatrix<index_t> sideActives = sideBasis.active(greville);

    gsMatrix<T> Phi  = gsMatrix<T>::Zero(nSide, nPts);
    gsMatrix<T> dPhi = gsMatrix<T>::Zero(nSide, nPts);
    for (index_t pt = 0; pt < nPts; ++pt)
        for (index_t a = 0; a < sideActives.rows(); ++a)
        {
            const index_t k = sideActives(a, pt);
            Phi(k, pt)  = sideVals(a, pt);
            dPhi(k, pt) = sideDeriv(a, pt);
        }

    gsSparseMatrix<T> sideColloc = collocationMatrix(sideBasis, greville);

    gsVector<T> tNorm(nPts), betaVals(nPts), alphaVals(nPts);
    for (index_t pt = 0; pt < nPts; ++pt)
    {
        tNorm(pt) = (greville(0, pt) - tA) / (tB - tA);
        betaVals(pt)  = beta0  * (T(1) - tNorm(pt)) + beta1  * tNorm(pt);
        alphaVals(pt) = alpha0 * (T(1) - tNorm(pt)) + alpha1 * tNorm(pt);
    }

    const index_t numRows = tensorBasis.size();
    const index_t numCols = nInterior + nLowerDeg + nSmoother;
    gsSparseMatrix<T> result(numRows, numCols);

    // Block 1: identity for interior DOFs
    index_t col0 = 0;
    for (const index_t i : interiorDOFs) { result(i, col0) = T(1); ++col0; }

    // Block 2: second-layer columns (d_i derivative = phi_j)
    const index_t colOffsetL2 = nInterior;
    {
        gsMatrix<T> lowerVals(nLowerDeg, nPts);
        for (index_t j = 0; j < nLowerDeg; ++j)
        {
            gsMatrix<T> tmp = sideLowerDegreeBasis.evalSingle(j, greville);
            lowerVals.row(j) = tmp.row(0);
        }
        for (index_t j = 0; j < nLowerDeg; ++j)
        {
            gsMatrix<T> rhs(nPts, 1);
            for (index_t pt = 0; pt < nPts; ++pt)
                rhs(pt, 0) = alphaVals(pt) * lowerVals(j, pt);

            gsMatrix<T> coeffs;
            makeSparseLUSolver(sideColloc)->apply(rhs, coeffs);

            const T scale = signN / dNeigh;
            for (index_t k = 0; k < nSide; ++k)
                if (std::abs(coeffs(k, 0)) > eps)
                {
                    const index_t row2D = secondLayerDOFs(k, 0);
                    result(row2D, colOffsetL2 + j) = coeffs(k, 0) * scale;
                }
        }
    }

    // Block 3: boundary columns (value = smoother basis, d_i = 0)
    const index_t colOffsetBdry = colOffsetL2 + nLowerDeg;
    {
        gsMatrix<T> E_dense = embFirstLayer.toDense();
        for (index_t m = 0; m < nSmoother; ++m)
        {
            gsVector<T> Vm(nPts), dVm(nPts);
            Vm.setZero(); dVm.setZero();
            for (index_t k = 0; k < nSide; ++k)
            {
                const T ek = E_dense(k, m);
                if (std::abs(ek) < eps) continue;
                for (index_t pt = 0; pt < nPts; ++pt)
                {
                    Vm(pt)  += ek * Phi(k, pt);
                    dVm(pt) += ek * dPhi(k, pt);
                }
            }
            gsMatrix<T> rhs(nPts, 1);
            for (index_t pt = 0; pt < nPts; ++pt)
                rhs(pt, 0) = -(dBdry * Vm(pt)
                               + tangentSign * signN * betaVals(pt) * dVm(pt)) / dNeigh;

            gsMatrix<T> gamma;
            makeSparseLUSolver(sideColloc)->apply(rhs, gamma);

            const index_t col2D = colOffsetBdry + m;
            for (index_t k = 0; k < nSide; ++k)
            {
                const index_t bdryRow = firstLayerDOFs(k, 0);
                if (std::abs(E_dense(k, m)) > eps)
                    result(bdryRow, col2D) = E_dense(k, m);
                if (std::abs(gamma(k, 0)) > eps)
                {
                    const index_t neighRow = bdryRow + signedStride;
                    result(neighRow, col2D) = gamma(k, 0);
                }
            }
        }
    }

    result.makeCompressed();
    return result;
}

} // namespace asg1v4
} // namespace gismo
