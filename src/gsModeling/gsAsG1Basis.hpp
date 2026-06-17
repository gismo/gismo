/** @file gsAsG1Basis.hpp

    @brief Derive a C1 basis for an AS-G1 geometry.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

namespace gismo {


// ====================================================================
// Helper utilities (same as in the per-side example)
// ====================================================================

template<typename T>
bool isNested(const gsBSplineBasis<T>& coarseBasis, const gsBSplineBasis<T>& fineBasis)
{
    if (coarseBasis.dim() != 1 || fineBasis.dim() != 1)
        return false;
    const int coarseDegree = coarseBasis.degree();
    const int fineDegree   = fineBasis.degree();
    if (coarseDegree > fineDegree)
        return false;
    const gsKnotVector<T>& coarseKnots = coarseBasis.knots();
    const gsKnotVector<T>& fineKnots   = fineBasis.knots();
    gsKnotVector<T> elevated = coarseKnots;
    elevated.degreeElevate(fineDegree - coarseDegree);
    return fineKnots.includes(elevated);
}

template<typename T>
gsSparseMatrix<T> collocationMatrix(const gsBSplineBasis<T>& basis,
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
gsSparseMatrix<T> embeddingMatrix(const gsBSplineBasis<T>& coarse,
                                  const gsBSplineBasis<T>& fine)
{
    gsMatrix<T> greville = fine.anchors();
    gsSparseMatrix<T> Cc = collocationMatrix(coarse, greville);
    gsSparseMatrix<T> Cf = collocationMatrix(fine,   greville);
    gsMatrix<T> result;
    makeSparseLUSolver(Cf)->apply(Cc, result);
    return result.sparseView(1, 1e-10);
}

std::vector<index_t> getInteriorDofs(const index_t tbSize,
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

// ====================================================================
// Per-side Argyris embedding with gluing-data directional derivative
// ====================================================================

/// Creates the per-side embedding for ONE side of a tensor-product basis.
///
/// Column layout:  [interior DOFs | second-layer cols | boundary cols]
///
/// - Interior DOFs: identity
/// - Boundary cols: value spans smoother basis, d_i derivative = 0
/// - Second-layer cols: zero value, d_i derivative spans lower-deg basis
template<typename T>
gsSparseMatrix<T> createGluingDataArgyrisBasis(
    const gsTensorBSplineBasis<2,T>& tensorBasis,
    boxSide side,
    T alpha0, T alpha1,
    T beta0,  T beta1,
    T eps = 1e-12,
    T tangentSign = T(1))    // Multiplies the beta*dVm term in the
                              // trace gamma RHS.  Pass -1 when the
                              // patch's tangent runs opposite to the
                              // gluing-data tangent (i.e. `flipped`).
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

    gsSparseMatrix<T> embFirstLayer  = embeddingMatrix(sideSmootherBasis, sideBasis);
    gsSparseMatrix<T> embSecondLayer = embeddingMatrix(sideLowerDegreeBasis, sideBasis);

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
    GISMO_ENSURE(std::abs(dNeigh) > eps,
                 "Neighbor derivative is zero.");

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
        betaVals(pt)  = beta0  * (1 - tNorm(pt)) + beta1  * tNorm(pt);
        alphaVals(pt) = alpha0 * (1 - tNorm(pt)) + alpha1 * tNorm(pt);
    }

    const index_t numRows = tensorBasis.size();
    const index_t numCols = nInterior + nLowerDeg + nSmoother;
    gsSparseMatrix<T> result(numRows, numCols);

    // Block 1: interior DOFs → identity
    index_t col0 = 0;
    for (const index_t i : interiorDOFs)
    {
        result(i, col0) = T(1);
        ++col0;
    }

    // Block 2: second-layer columns  (d_i derivative = phi_j)
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

    // Block 3: boundary columns  (value = smoother basis, d_i derivative = 0)
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
                rhs(pt, 0) = -(dBdry * Vm(pt) + tangentSign * signN * betaVals(pt) * dVm(pt)) / dNeigh;

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




} // namespace gismo
