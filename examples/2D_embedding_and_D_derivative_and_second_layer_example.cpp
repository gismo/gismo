/** @file 2D_embedding_and_D_derivative_and_second_layer_example.cpp

    @brief Build an Argyris-type basis transformation that enforces
           zero directional derivative d_i at each patch boundary,
           where the direction is defined by gluing data (alpha, beta):

               d_i = (1 / alpha_i) * (n + beta_i * t)

           Here n and t are the normal and tangential partial
           derivatives in the parameter domain.

    For boundary sides (no neighbouring patch), alpha=1, beta=0,
    which recovers the standard zero-normal-derivative constraint.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <iostream>
#include <set>
#include <algorithm>
#include <numeric>
#include <gismo.h>

using namespace gismo;


/// Check whether the coarse B-spline basis is nested inside the fine one.
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

    gsKnotVector<T> coarseKnotsWithIncreasedMultiplicity = coarseKnots;
    coarseKnotsWithIncreasedMultiplicity.degreeElevate(fineDegree-coarseDegree);


    if (!fineKnots.includes(coarseKnotsWithIncreasedMultiplicity)){
        return false;
    }

    return true;
}


/// Build a collocation matrix: row i = evaluation of all basis functions
/// at collocation point i.
template<typename T>
gsSparseMatrix<T> collocationMatrix(const gsBSplineBasis<T>& basis,
                                    const gsMatrix<T>& collocationPoints)
{
    gsMatrix<T> values = basis.eval(collocationPoints);
    gsMatrix<index_t> indices = basis.active(collocationPoints);

    gsSparseEntries<T> entries;
    entries.reserve(values.rows() * values.cols());

    for (int i = 0; i < values.cols(); ++i)
        for (int j = 0; j < values.rows(); ++j)
            entries.add(i, indices(j, i), values(j, i));

    gsSparseMatrix<T> collocation(collocationPoints.cols(), basis.size());
    collocation.setFrom(entries);

    return collocation;
}

/// Compute the embedding (refinement) matrix from a coarse basis into
/// a fine basis via Greville collocation.
template<typename T>
gsSparseMatrix<T> embeddingMatrix(const gsBSplineBasis<T>& coarseBasis,
                                  const gsBSplineBasis<T>& fineBasis)
{
    gsMatrix<T> greville = fineBasis.anchors();
    gsSparseMatrix<T> coarseCollocation = collocationMatrix(coarseBasis, greville);
    gsSparseMatrix<T> fineCollocation = collocationMatrix(fineBasis, greville);
    gsMatrix<T> result;
    makeSparseLUSolver(fineCollocation)->apply(coarseCollocation,result);
    gsSparseMatrix<T> embedding = result.sparseView(1,1e-10);

    return embedding;
}


/// Return the sorted indices of all DOFs that are neither in the
/// first-layer (boundary) nor in the second-layer set.
std::vector<index_t> getInteriorDofs(const index_t tbSize,
                                     const gsMatrix<index_t>& firstLayerDOFs,
                                     const gsMatrix<index_t>& secondLayerDOFs)
{
    // Merge and sort both boundary DOF sets
    std::vector<index_t> removedSet;
    removedSet.reserve(firstLayerDOFs.rows() + secondLayerDOFs.rows());
    std::merge(firstLayerDOFs.data(), firstLayerDOFs.data() + firstLayerDOFs.rows(),
               secondLayerDOFs.data(), secondLayerDOFs.data() + secondLayerDOFs.rows(),
               std::back_inserter(removedSet));
    std::sort(removedSet.begin(), removedSet.end());

    // All DOF indices
    std::vector<index_t> allDOFs(tbSize);
    std::iota(allDOFs.begin(), allDOFs.end(), 0);

    // Interior = all \ removed
    std::vector<index_t> trueInteriorDOFs;
    trueInteriorDOFs.reserve(allDOFs.size() - removedSet.size());
    std::set_difference(allDOFs.begin(), allDOFs.end(),
                        removedSet.begin(), removedSet.end(),
                        std::back_inserter(trueInteriorDOFs));
    return trueInteriorDOFs;
}


/// Creates a tensor embedding matrix that enforces the directional
/// derivative constraint defined by gluing data (alpha, beta):
///
///   d_i = (1/alpha_i) * (partial_n + beta_i * partial_t)
///
/// where alpha_i(t) = alpha0*(1-t) + alpha1*t  and similarly for beta_i,
/// with t in [0,1] normalised along the tangential direction.
///
/// Column layout of the result:
///   [interior DOFs | second-layer columns | boundary columns]
///
/// - Interior DOFs: identity (these have zero value AND zero d_i at boundary)
/// - Boundary columns: value spans smoother basis, d_i derivative = 0
/// - Second-layer columns: zero value, d_i derivative spans lower-degree basis
///
/// For alpha=1, beta=0 this recovers the standard zero-normal-derivative case.
template<typename T>
gsSparseMatrix<T> createGluingDataArgyrisBasis(
    const gsTensorBSplineBasis<2,T>& tensorBasis,
    boxSide side,
    T alpha0, T alpha1,
    T beta0,  T beta1,
    T eps = 1e-12)
{
    // --- Tangential (side) basis and its sub-bases ---
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

    // 1D embeddings (tangential direction)
    gsSparseMatrix<T> embFirstLayer  = embeddingMatrix(sideSmootherBasis, sideBasis);
    gsSparseMatrix<T> embSecondLayer = embeddingMatrix(sideLowerDegreeBasis, sideBasis);

    // --- Direction setup ---
    const bool    isLow  = !side.parameter();
    const int     dir    = side.direction();
    const int     tanDir = 1 - dir;
    const index_t stride       = (dir == 0) ? 1 : tensorBasis.size(0);
    const index_t signedStride = isLow ? stride : -stride;

    // 1D normal-direction basis
    const gsBSplineBasis<T> normalBasis = tensorBasis.component(dir);
    const index_t nNormal = normalBasis.size();
    const index_t bdryIdx  = isLow ? 0 : nNormal - 1;
    const index_t neighIdx = isLow ? 1 : nNormal - 2;

    gsMatrix<T> bdryPt1D = normalBasis.support().col(isLow ? 0 : 1);
    const T dBdry  = normalBasis.derivSingle(bdryIdx,  bdryPt1D)(0, 0);
    const T dNeigh = normalBasis.derivSingle(neighIdx, bdryPt1D)(0, 0);
    GISMO_ENSURE(std::abs(dNeigh) > eps,
                 "Neighbor derivative is zero, cannot enforce constraints.");

    // Sign correction: outward normal direction
    const T signN = isLow ? T(-1) : T(1);

    // --- Collect the three disjoint DOF sets ---
    gsMatrix<index_t> firstLayerDOFs  = tensorBasis.boundary(side);
    gsMatrix<index_t> secondLayerDOFs = tensorBasis.boundaryOffset(side, 1);
    std::vector<index_t> interiorDOFs =
        getInteriorDofs(tensorBasis.size(), firstLayerDOFs, secondLayerDOFs);

    const index_t nSide       = sideBasis.size();
    const index_t nSmoother   = sideSmootherBasis.size();
    const index_t nLowerDeg   = sideLowerDegreeBasis.size();
    const index_t nInterior   = static_cast<index_t>(interiorDOFs.size());

    // --- Build collocation tools along the tangential direction ---
    // We use Greville abscissae of the side basis as collocation points.
    gsMatrix<T> greville = sideBasis.anchors();       // 1 x nSide
    const index_t nPts = greville.cols();

    // Tangential parameter range for normalisation to [0,1]
    gsMatrix<T> supTan = tensorBasis.component(tanDir).support();
    const T tA = supTan(0, 0), tB = supTan(0, 1);

    // Evaluate side basis value and derivative at Greville points
    gsMatrix<T> sideVals  = sideBasis.eval(greville);   // nActive x nPts
    gsMatrix<T> sideDeriv = sideBasis.deriv(greville);   // nActive x nPts
    gsMatrix<index_t> sideActives = sideBasis.active(greville);

    // Build dense evaluation matrices: Phi(k, pt) and Phi'(k, pt)
    gsMatrix<T> Phi  = gsMatrix<T>::Zero(nSide, nPts);
    gsMatrix<T> dPhi = gsMatrix<T>::Zero(nSide, nPts);
    for (index_t pt = 0; pt < nPts; ++pt)
        for (index_t a = 0; a < sideActives.rows(); ++a)
        {
            const index_t k = sideActives(a, pt);
            Phi(k, pt)  = sideVals(a, pt);
            dPhi(k, pt) = sideDeriv(a, pt);
        }

    // Side-basis collocation matrix (nPts x nSide) for solving
    gsSparseMatrix<T> sideColloc = collocationMatrix(sideBasis, greville);

    // Precompute normalised t and beta(t), alpha(t) at each Greville point
    gsVector<T> tNorm(nPts), betaVals(nPts), alphaVals(nPts);
    for (index_t pt = 0; pt < nPts; ++pt)
    {
        tNorm(pt) = (greville(0, pt) - tA) / (tB - tA);
        betaVals(pt)  = beta0  * (1 - tNorm(pt)) + beta1  * tNorm(pt);
        alphaVals(pt) = alpha0 * (1 - tNorm(pt)) + alpha1 * tNorm(pt);
    }

    gsInfo << "  Side: " << side
           << "  alpha = " << alpha0 << "*(1-t)+" << alpha1 << "*t"
           << "  beta = "  << beta0  << "*(1-t)+" << beta1  << "*t\n";

    // =================================================================
    // Build the result matrix
    // =================================================================
    const index_t numRows = tensorBasis.size();
    const index_t numCols = nInterior + nLowerDeg + nSmoother;
    gsSparseMatrix<T> result(numRows, numCols);

    // ----- Block 1: identity for truly interior DOFs -----
    index_t col0 = 0;
    for (const index_t i : interiorDOFs)
    {
        result(i, col0) = T(1);
        ++col0;
    }

    // ----- Block 2: second-layer columns (d_i spans lower-degree basis) -----
    // These columns have zero value at the boundary (the second-layer
    // basis function N_{neighIdx} vanishes there).  Since the value is
    // zero, the tangential derivative is also zero, and d_i reduces to:
    //
    //   d_i u = (1/alpha) * partial_n u  =  phi_j(t)
    //
    // The outward normal derivative is partial_n = signN * partial_param,
    // and the parametric normal derivative of a second-layer function is:
    //   partial_param u(t) = dNeigh * sum_k c_k M_k(t)
    //
    // So we need:
    //   signN * dNeigh * sum_k c_k M_k(t) = alpha(t) * phi_j(t)
    //
    // i.e.  c_k coefficients are found by collocation on
    //   alpha(t) * phi_j(t), then scaled by signN / dNeigh
    //   (since signN^2 = 1, this equals 1 / (signN * dNeigh)).

    const index_t colOffsetL2 = nInterior;
    {
        // Evaluate lower-degree basis at Greville points
        gsMatrix<T> lowerVals(nLowerDeg, nPts);
        for (index_t j = 0; j < nLowerDeg; ++j)
        {
            gsMatrix<T> tmp = sideLowerDegreeBasis.evalSingle(j, greville);
            lowerVals.row(j) = tmp.row(0);
        }

        for (index_t j = 0; j < nLowerDeg; ++j)
        {
            // RHS: alpha(t_i) * phi_j(t_i) at each Greville point
            gsMatrix<T> rhs(nPts, 1);
            for (index_t pt = 0; pt < nPts; ++pt)
                rhs(pt, 0) = alphaVals(pt) * lowerVals(j, pt);

            // Solve for side-basis coefficients
            gsMatrix<T> coeffs;
            makeSparseLUSolver(sideColloc)->apply(rhs, coeffs);

            // Scale by signN / dNeigh and place into second-layer rows
            const T scale = signN / dNeigh;
            for (index_t k = 0; k < nSide; ++k)
            {
                if (std::abs(coeffs(k, 0)) > eps)
                {
                    const index_t row2D = secondLayerDOFs(k, 0);
                    result(row2D, colOffsetL2 + j) = coeffs(k, 0) * scale;
                }
            }
        }
    }

    // ----- Block 3: boundary columns -----
    // These columns should have:
    //   value at boundary = smoother basis function (via embFirstLayer)
    //   d_i derivative at boundary = 0
    //
    // The d_i constraint at the boundary is:
    //   (1/alpha) * [partial_n u + beta * partial_t u] = 0
    //
    // Since alpha != 0 this is:
    //   partial_n u + beta(t) * partial_t u = 0      (*)
    //
    // The outward normal derivative relates to the parametric derivative:
    //   partial_n = signN * partial_param
    //
    // A boundary column m has:
    //   boundary row k:       E_{k,m}   (from embFirstLayer)
    //   second-layer row k:   gamma_{k,m}  (to be determined)
    //
    // At the boundary, the parametric derivatives are:
    //   partial_param u(t) = dBdry  * sum_k E_{k,m} M_k(t)
    //                      + dNeigh * sum_k gamma_{k,m} M_k(t)
    //   partial_t u(t)     = sum_k E_{k,m} M_k'(t)
    //                        (second-layer value is zero at boundary)
    //
    // Substituting partial_n = signN * partial_param into (*):
    //   signN * (dBdry * V_m + dNeigh * G_m) + beta * V'_m = 0
    //
    // Solving for G_m at collocation points:
    //   G_m(t_i) = -(1/dNeigh) * [dBdry * V_m(t_i) + signN * beta(t_i) * V'_m(t_i)]
    //
    // Then we recover gamma_{k,m} by solving the collocation system.

    const index_t colOffsetBdry = colOffsetL2 + nLowerDeg;
    {
        // Precompute the dense embedding: E_dense(k, m) for k in sideBasis, m in smoother
        gsMatrix<T> E_dense = embFirstLayer.toDense();  // nSide x nSmoother

        for (index_t m = 0; m < nSmoother; ++m)
        {
            // Compute V_m(t_i) and V'_m(t_i) at Greville points
            gsVector<T> Vm(nPts), dVm(nPts);
            Vm.setZero();
            dVm.setZero();
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

            // RHS for gamma collocation:
            //   signN * (dBdry * V_m + dNeigh * G_m) + beta * V'_m = 0
            //   => G_m = -(1/dNeigh) * [dBdry * V_m + signN * beta * V'_m]
            gsMatrix<T> rhs(nPts, 1);
            for (index_t pt = 0; pt < nPts; ++pt)
                rhs(pt, 0) = -(dBdry * Vm(pt) + signN * betaVals(pt) * dVm(pt)) / dNeigh;

            // Solve for gamma coefficients in the side basis
            gsMatrix<T> gamma;
            makeSparseLUSolver(sideColloc)->apply(rhs, gamma);

            // Fill the result matrix
            const index_t col2D = colOffsetBdry + m;
            for (index_t k = 0; k < nSide; ++k)
            {
                // Boundary row: embedding coefficient
                const index_t bdryRow = firstLayerDOFs(k, 0);
                if (std::abs(E_dense(k, m)) > eps)
                    result(bdryRow, col2D) = E_dense(k, m);

                // Second-layer row: gamma coefficient
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


// ====================================================================
// Gluing-data computation
// (computes alpha, beta functions for each interface)
// ====================================================================

/// 2D cross-determinant  det([a | b])
template <class T>
inline T det2(const gsVector<T>& a, const gsVector<T>& b)
{ return a(0)*b(1) - a(1)*b(0); }

/// Extract partial-derivative vector from deriv_into output.
template <class T>
inline gsVector<T> getPartial(const gsMatrix<T>& d, index_t dir, index_t col)
{
    gsVector<T> v(2);
    v(0) = d(dir,     col);
    v(1) = d(dir + 2, col);
    return v;
}

/// Collect knot-span break points for a 1D component basis direction.
template <class T>
std::vector<T> collectBreaks(const gsGeometry<T>& geo, short_t dir)
{
    const gsBSplineBasis<T>* bb =
        dynamic_cast<const gsBSplineBasis<T>*>(&geo.basis().component(dir));
    if (bb) return bb->knots().breaks();
    gsMatrix<T> sup = geo.support();
    return { sup(dir, 0), sup(dir, 1) };
}

/// Compute gluing data for a single interface.
/// Returns [a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1].
template <class T>
gsVector<T> computeGluingDataForInterface(
    const gsMultiPatch<T>& mp,
    const boundaryInterface& interf,
    bool& success,
    const T eps = 1e-8,
    index_t numGaussPerSpan = 0)
{
    success = false;
    gsVector<T> result(8);
    result.setZero();

    const patchSide ps1 = interf.first();
    const patchSide ps2 = interf.second();
    const gsGeometry<T>& geo1 = mp.patch(ps1.patch);
    const gsGeometry<T>& geo2 = mp.patch(ps2.patch);

    const short_t normDir1 = ps1.direction(), tanDir1 = 1 - normDir1;
    const short_t normDir2 = ps2.direction(), tanDir2 = 1 - normDir2;
    const bool par1 = ps1.parameter(), par2 = ps2.parameter();
    const T signD1 = par1 ? T(-1) : T(1);
    const T signD2 = par2 ? T(-1) : T(1);
    const T signD3 = signD1 * signD2;

    gsMatrix<T> sup1 = geo1.support(), sup2 = geo2.support();
    const T ifcCoord1 = sup1(normDir1, par1 ? 1 : 0);
    const T ifcCoord2 = sup2(normDir2, par2 ? 1 : 0);
    const bool flipped = !interf.dirOrientation(ps1, tanDir1);

    const short_t deg1 = geo1.basis().degree(tanDir1);
    const short_t deg2 = geo2.basis().degree(tanDir2);
    const index_t nGauss = numGaussPerSpan > 0
                               ? numGaussPerSpan
                               : 2 * std::max(deg1, deg2) + 1;

    std::vector<T> breaks1 = collectBreaks(geo1, tanDir1);
    std::vector<T> breaks2 = collectBreaks(geo2, tanDir2);
    const T t1a = sup1(tanDir1, 0), t1b = sup1(tanDir1, 1);
    const T t2a = sup2(tanDir2, 0), t2b = sup2(tanDir2, 1);
    for (T& br : breaks2)
    {
        T s = (br - t2a) / (t2b - t2a);
        if (flipped) s = 1.0 - s;
        br = t1a + s * (t1b - t1a);
    }
    std::set<T> breakSet(breaks1.begin(), breaks1.end());
    breakSet.insert(breaks2.begin(), breaks2.end());
    std::vector<T> mergedBreaks(breakSet.begin(), breakSet.end());

    gsGaussRule<T> gaussRule(nGauss);
    gsMatrix<T> gaussNodes;
    gsVector<T> w;
    gaussRule.mapToAll(mergedBreaks, gaussNodes, w);
    const index_t N = gaussNodes.cols();

    gsMatrix<T> pts1(2, N), pts2(2, N);
    for (index_t i = 0; i < N; ++i)
    {
        T t = gaussNodes(0, i);
        pts1(normDir1, i) = ifcCoord1;
        pts1(tanDir1,  i) = t;
        T s = (t - t1a) / (t1b - t1a);
        if (flipped) s = 1.0 - s;
        pts2(normDir2, i) = ifcCoord2;
        pts2(tanDir2,  i) = t2a + s * (t2b - t2a);
    }

    gsMatrix<T> derivs1, derivs2;
    geo1.deriv_into(pts1, derivs1);
    geo2.deriv_into(pts2, derivs2);

    gsVector<T> D1(N), D2(N), D3(N), t_vals(N);
    for (index_t i = 0; i < N; ++i)
    {
        gsVector<T> dG1dn = getPartial(derivs1, normDir1, i);
        gsVector<T> dG1dt = getPartial(derivs1, tanDir1,  i);
        gsVector<T> dG2dn = getPartial(derivs2, normDir2, i);
        gsVector<T> dG2dt = getPartial(derivs2, tanDir2,  i);
        D1(i) = signD1 * det2(dG1dn, dG1dt);
        D2(i) = signD2 * det2(dG2dn, dG2dt);
        D3(i) = signD3 * det2(dG1dn, dG2dn);
        t_vals(i) = (gaussNodes(0, i) - t1a) / (t1b - t1a);
    }

    if (D1.minCoeff() * D1.maxCoeff() < 0 ||
        D2.minCoeff() * D2.maxCoeff() < 0)
        return result;

    auto integrate = [&](const gsVector<T>& f) -> T { return w.dot(f); };

    // --- Try constant alpha ---
    T intD1D1 = integrate((D1.array() * D1.array()).matrix());
    T intD1D2 = integrate((D1.array() * D2.array()).matrix());
    T intD2D2 = integrate((D2.array() * D2.array()).matrix());
    T denom = intD1D1 - 2 * intD1D2 + intD2D2;

    T a1c, a2c;
    if (std::abs(denom) < 1e-30)
    { a1c = 0.5; a2c = 0.5; }
    else
    { a1c = (intD1D1 - intD1D2) / denom; a2c = 1.0 - a1c; }

    T aerr = a1c * a1c * intD2D2 + 2 * a1c * a2c * intD1D2 + a2c * a2c * intD1D1;

    T a10, a11, a20, a21, b10, b11, b20, b21;

    if (aerr < eps)
    {
        a10 = a1c; a11 = a1c; a20 = a2c; a21 = a2c;
        // Solve for beta with b1 = b2
        gsVector<T> S = D2 - D1;
        gsMatrix<T> Ab(2, 2); Ab.setZero();
        gsVector<T> rb(2);    rb.setZero();
        for (index_t i = 0; i < N; ++i)
        {
            T phi[2] = {1 - t_vals(i), t_vals(i)};
            for (index_t j = 0; j < 2; ++j)
            {
                for (index_t k = 0; k < 2; ++k)
                    Ab(j, k) += w(i) * phi[j] * phi[k] * S(i) * S(i);
                rb(j) += w(i) * phi[j] * S(i) * D3(i);
            }
        }
        gsVector<T> bc = Ab.fullPivLu().solve(rb);
        b10 = bc(0); b11 = bc(1); b20 = b10; b21 = b11;
        success = true;
    }
    else
    {
        // General linear alpha via KKT
        gsMatrix<T> fv(N, 4);
        for (index_t i = 0; i < N; ++i)
        {
            T p0 = 1 - t_vals(i), p1 = t_vals(i);
            fv(i, 0) = p0 * D2(i); fv(i, 1) = p1 * D2(i);
            fv(i, 2) = p0 * D1(i); fv(i, 3) = p1 * D1(i);
        }
        gsMatrix<T> G(4, 4); G.setZero();
        for (index_t j = 0; j < 4; ++j)
            for (index_t k = j; k < 4; ++k)
            {
                gsVector<T> pr(N);
                for (index_t i = 0; i < N; ++i)
                    pr(i) = fv(i, j) * fv(i, k);
                G(j, k) = integrate(pr);
                G(k, j) = G(j, k);
            }
        gsVector<T> c(4); c.setConstant(0.5);
        gsMatrix<T> KKT(5, 5); KKT.setZero();
        KKT.block(0, 0, 4, 4) = 2.0 * G;
        KKT.block(0, 4, 4, 1) = c;
        KKT.block(4, 0, 1, 4) = c.transpose();
        gsVector<T> rhs(5); rhs.setZero(); rhs(4) = 1.0;
        gsVector<T> sol = KKT.fullPivLu().solve(rhs);
        a10 = sol(0); a11 = sol(1); a20 = sol(2); a21 = sol(3);

        // General linear beta
        gsMatrix<T> H(4, 4); H.setZero();
        gsVector<T> d_vec(4); d_vec.setZero();
        for (index_t j = 0; j < 4; ++j)
        {
            T sj = (j < 2) ? 1.0 : -1.0;
            for (index_t k = j; k < 4; ++k)
            {
                T sk = (k < 2) ? 1.0 : -1.0;
                gsVector<T> pr(N);
                for (index_t i = 0; i < N; ++i)
                    pr(i) = sj * fv(i, j) * sk * fv(i, k);
                H(j, k) = integrate(pr); H(k, j) = H(j, k);
            }
            gsVector<T> pr(N);
            for (index_t i = 0; i < N; ++i)
                pr(i) = sj * fv(i, j) * D3(i);
            d_vec(j) = integrate(pr);
        }
        gsVector<T> e(4);
        {
            gsVector<T> tmp(N);
            for (index_t i = 0; i < N; ++i)
                tmp(i) = (a10 * (1 - t_vals(i)) + a11 * t_vals(i)) * (1 - t_vals(i));
            e(0) = integrate(tmp);
            for (index_t i = 0; i < N; ++i)
                tmp(i) = (a10 * (1 - t_vals(i)) + a11 * t_vals(i)) * t_vals(i);
            e(1) = integrate(tmp);
            for (index_t i = 0; i < N; ++i)
                tmp(i) = -(a20 * (1 - t_vals(i)) + a21 * t_vals(i)) * (1 - t_vals(i));
            e(2) = integrate(tmp);
            for (index_t i = 0; i < N; ++i)
                tmp(i) = -(a20 * (1 - t_vals(i)) + a21 * t_vals(i)) * t_vals(i);
            e(3) = integrate(tmp);
        }
        gsMatrix<T> KKT2(5, 5); KKT2.setZero();
        KKT2.block(0, 0, 4, 4) = 2.0 * H;
        KKT2.block(0, 4, 4, 1) = e;
        KKT2.block(4, 0, 1, 4) = e.transpose();
        gsVector<T> rhs2(5);
        rhs2.head(4) = 2.0 * d_vec; rhs2(4) = 0.0;
        gsVector<T> sol2 = KKT2.fullPivLu().solve(rhs2);
        b10 = sol2(0); b11 = sol2(1); b20 = sol2(2); b21 = sol2(3);
        success = true;
    }

    if (!success) return result;

    // Sign correction for beta: the fitting solves beta*(D2-D1) ≈ D3,
    // but the C1 condition requires -D2*beta_0 + D1*beta_1 = D3,
    // i.e., beta*(D1-D2) = D3. So negate all beta values.
    b10 = -b10; b11 = -b11; b20 = -b20; b21 = -b21;

    result(0) = a10; result(1) = a11;
    result(2) = b10; result(3) = b11;
    if (flipped)
    { result(4) = a21; result(5) = a20; result(6) = b21; result(7) = b20; }
    else
    { result(4) = a20; result(5) = a21; result(6) = b20; result(7) = b21; }
    return result;
}

/// Collect all interfaces into an (nPatches x 16) gluing data matrix.
/// Layout: 4 columns per side (west=0, east=1, south=2, north=3),
///   columns 4*s+0..1 = alpha_0, alpha_1
///   columns 4*s+2..3 = beta_0,  beta_1
/// Boundary sides get alpha=1, beta=0.
template <class T>
gsMatrix<T> computeGluingDataMatrix(
    const gsMultiPatch<T>& mp,
    const T eps = 1e-8,
    index_t numGaussPerSpan = 0)
{
    const index_t nP = mp.nPatches();
    gsMatrix<T> mat(nP, 16);
    for (index_t p = 0; p < nP; ++p)
        for (index_t s = 0; s < 4; ++s)
        {
            mat(p, 4 * s + 0) = 1.0;
            mat(p, 4 * s + 1) = 1.0;
            mat(p, 4 * s + 2) = 0.0;
            mat(p, 4 * s + 3) = 0.0;
        }

    auto sideToCol = [](const patchSide& ps) -> index_t
    { return static_cast<index_t>(ps.side()) - 1; };

    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const boundaryInterface& interf = *it;
        bool ok = false;
        gsVector<T> gd = computeGluingDataForInterface(
            mp, interf, ok, eps, numGaussPerSpan);
        if (!ok) continue;

        const patchSide& ps1 = interf.first();
        const patchSide& ps2 = interf.second();
        const index_t col1 = 4 * sideToCol(ps1);
        const index_t col2 = 4 * sideToCol(ps2);

        for (index_t i = 0; i < 4; ++i)
            mat(ps1.patch, col1 + i) = gd(i);
        for (index_t i = 0; i < 4; ++i)
            mat(ps2.patch, col2 + i) = gd(4 + i);
    }
    return mat;
}


// ====================================================================
// main
// ====================================================================

int main(int argc, char* argv[])
{
    using T = double;

    std::string geometry("domain2d/2patch/two_bilinear_patches.xml");
    index_t numGaussPerSpan = 0;
    index_t refinements = 0;
    index_t plot = -1;

    gsCmdLine cmd("Argyris basis with gluing-data directional derivatives.");
    cmd.addString("f", "file", "G+Smo multi-patch geometry file.", geometry);
    cmd.addInt("n", "numGauss",
               "Number of Gauss points per knot span (0=auto).", numGaussPerSpan);
    cmd.addInt("r", "refinements",
               "Uniform refinements before proceeding.", refinements);
    cmd.addInt("p", "plot",
               "Plot basis function with that column index (per patch).", plot);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // ---- Read multi-patch geometry ----
    gsMultiPatch<T>::uPtr mpPtr = gsReadFile<>(geometry);
    if (!mpPtr) { gsInfo << "Cannot read " << geometry << ".\n"; return -1; }
    gsMultiPatch<T>& mp = *mpPtr;
    mp.computeTopology();

    gsInfo << "Patches: " << mp.nPatches()
           << "  Interfaces: " << mp.nInterfaces() << "\n\n";

    for (index_t i = 0; i < refinements; ++i)
        mp.uniformRefine(1, 2);

    // Check that refinement produced the required multiplicity
    if (refinements == 0)
        gsInfo << "WARNING: without refinements (-r), the basis may lack\n"
               << "  the required interior knot multiplicity > 1.\n"
               << "  Use -r 1 (or higher) to refine first.\n\n";

    // ---- Compute gluing data ----
    gsMatrix<T> gluingData = computeGluingDataMatrix(mp, T(1e-8), numGaussPerSpan);

    const char* sideNames[4] = {"west", "east", "south", "north"};
    gsInfo << "=== Gluing data ===\n";
    for (size_t p = 0; p < mp.nPatches(); ++p)
    {
        gsInfo << "Patch " << p << ":\n";
        for (index_t s = 0; s < 4; ++s)
        {
            T a0 = gluingData(p, 4 * s + 0), a1 = gluingData(p, 4 * s + 1);
            T b0 = gluingData(p, 4 * s + 2), b1 = gluingData(p, 4 * s + 3);
            gsInfo << "  " << std::setw(5) << sideNames[s]
                   << ":  alpha=" << a0 << "*(1-t)+" << a1 << "*t"
                   << "   beta=" << b0 << "*(1-t)+" << b1 << "*t\n";
        }
    }
    gsInfo << "\n";

    // ---- Build the Argyris basis transformation for each patch/side ----
    // Side ordering: west=1, east=2, south=3, north=4 in boxSide convention
    const boxSide allSides[4] = {
        boxSide(1), boxSide(2), boxSide(3), boxSide(4)
    };

    for (size_t p = 0; p < mp.nPatches(); ++p)
    {
        gsInfo << "========== Patch " << p << " ==========\n";

        const gsTensorBSplineBasis<2, T>& tb =
            dynamic_cast<const gsTensorBSplineBasis<2, T>&>(mp.patch(p).basis());

        for (index_t s = 0; s < 4; ++s)
        {
            const T a0 = gluingData(p, 4 * s + 0);
            const T a1 = gluingData(p, 4 * s + 1);
            const T b0 = gluingData(p, 4 * s + 2);
            const T b1 = gluingData(p, 4 * s + 3);

            gsInfo << "\n  --- Side " << sideNames[s]
                   << " (alpha=" << a0 << "*(1-t)+" << a1 << "*t"
                   << ", beta=" << b0 << "*(1-t)+" << b1 << "*t) ---\n";

            gsSparseMatrix<T> E = createGluingDataArgyrisBasis(
                tb, allSides[s], a0, a1, b0, b1);

            gsInfo << "  Embedding size: " << E.rows() << " x " << E.cols() << "\n";
        }
        gsInfo << "\n";
    }

    // ---- Single-side demo (patch 0, west) for debugging / plotting ----
    if (mp.nPatches() > 0)
    {
        const index_t demoPatch = 0;
        const index_t demoSideIdx = 0; // west
        gsInfo << "========== Single-side demo (patch " << demoPatch
               << ", " << sideNames[demoSideIdx] << ") ==========\n";

        const gsTensorBSplineBasis<2, T>& tb0 =
            dynamic_cast<const gsTensorBSplineBasis<2, T>&>(
                mp.patch(demoPatch).basis());

        const T a0 = gluingData(demoPatch, 4 * demoSideIdx + 0);
        const T a1 = gluingData(demoPatch, 4 * demoSideIdx + 1);
        const T b0 = gluingData(demoPatch, 4 * demoSideIdx + 2);
        const T b1 = gluingData(demoPatch, 4 * demoSideIdx + 3);

        gsSparseMatrix<T> E = createGluingDataArgyrisBasis(
            tb0, allSides[demoSideIdx], a0, a1, b0, b1);

        gsInfo << "Embedding matrix size: " << E.rows() << " x " << E.cols() << "\n";
        gsInfo << "Embedding matrix:\n" << E << "\n";

        if (plot > -1 && plot < E.cols())
        {
            gsMultiPatch<T> mpViz, mpsol;
            mpViz.addPatch(*gsNurbsCreator<T>::BSplineRectangle());
            mpsol.addPatch(tb0.makeGeometry(E.col(plot)));
            gsWriteParaview<>(gsField<T>(mpViz, mpsol), "basis_result", 1000);
            gsInfo << "Wrote basis_result.pvd for column " << plot
                   << " — open in ParaView to inspect.\n";
        }
    }

    gsInfo << "\nDone.\n";
    return 0;
}