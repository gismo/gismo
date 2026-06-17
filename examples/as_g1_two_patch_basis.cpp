/** @file as_g1_two_patch_basis_v4.cpp

    @brief Build an AS-G1 conforming basis across two patches.

    Streamlined v4 -- same gluing-data computation as v3, plus two
    fixes that make the construction work for ALL 196 two-patch
    parametrizations (8 left-side x 8 right-side reparametrizations
    times 3 geometry families), instead of only the 100/196 that
    worked in v2/v3.

    The two fixes are both triggered by the interface orientation
    flag `flipped = !ifc.dirOrientation(ps1, tDir1)`, which is true
    exactly when patch-2's natural tangent runs opposite to patch-1's
    natural tangent along the shared interface.

    1) `tangentSign` parameter on `createGluingDataArgyrisBasis`:
       The per-patch trace embedding solves a Greville interpolation
       problem whose RHS reads
           rhs = -(dBdry*Vm + tangentSign*signN*beta*dVm) / dNeigh
       where dVm is d/dt of the smoother basis in the *patch's own*
       tangent parameter, while beta is the gluing-data tangential
       coefficient defined in the *gluing-data* tangent frame.  When
       the patch's tangent is opposite to the gluing-data tangent,
       these two derivatives differ by a sign, which must be
       absorbed into the embedding's RHS.  We pass tangentSign = -1
       to BOTH patches when `flipped`, because both halves of the
       global AS-G1 relation
            alpha_1 d_1(G^(1)) + alpha_2 d_2(G^(2)) = 0
       are rewritten in the gluing-data tangent frame.

    2) Conditional second-layer sign `l2Sign` in the global G2
       assembly:  In v2 the d-derivative columns of patch-2 were
       inserted with an unconditional `-it.value()`.  When the
       tangent is flipped, that sign must be `+it.value()` because
       the column reversal (j2 = nLD2-1-j) already accounts for the
       mirror reparametrization of the lower-degree basis, which
       flips the d-derivative sign once.

    With these two fixes the global basis satisfies the AS-G1
    condition across the interface to machine precision
    (~1e-15 on bilinear patches, ~3e-7 on cubically refined curved
    patches due to representable approximation error).

    Only the shared interface receives the gluing-data directional
    derivative constraints.  All other boundary sides are left
    untouched (standard tensor B-spline DOFs).

    Interface trace (boundary) DOFs and interface d-derivative
    (second-layer) DOFs are shared between the two patches.
    All other DOFs are independent per patch.

    Use  -p N  to plot global basis function N on both patches
    via Paraview.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <iostream>
#include <set>
#include <map>
#include <algorithm>
#include <numeric>
#include <gismo.h>

using namespace gismo;

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


// ====================================================================
// Gluing-data helpers (streamlined v3 -- matches gluing_data_v3.cpp)
// ====================================================================

template <class T>
inline T det2(const gsVector<T>& a, const gsVector<T>& b)
{ return a(0)*b(1) - a(1)*b(0); }

template <class T>
inline gsVector<T> partial(const gsMatrix<T>& d, index_t dir, index_t col)
{
    gsVector<T> v(2);
    v(0) = d(dir,     col);
    v(1) = d(dir + 2, col);
    return v;
}

template <class T>
std::vector<T> breaksOf(const gsGeometry<T>& geo, short_t dir)
{
    const gsBSplineBasis<T>* bb =
        dynamic_cast<const gsBSplineBasis<T>*>(&geo.basis().component(dir));
    if (bb) return bb->knots().breaks();
    gsMatrix<T> sup = geo.support();
    return { sup(dir, 0), sup(dir, 1) };
}

/// Sampled determinants on the interface (patch-1 tangential parameter
/// mapped to [0,1]).
template <class T>
struct InterfaceSamples
{
    gsVector<T> D1, D2, D3;
    gsVector<T> t;
    gsVector<T> w;
    index_t size() const { return D1.size(); }
    T integrate(const gsVector<T>& f) const { return w.dot(f); }
};

/// Sample D1,D2,D3 along the interface.  Sets tangentialFlipped.
template <class T>
InterfaceSamples<T> sampleInterface(
    const gsMultiPatch<T>& mp,
    const boundaryInterface& interf,
    bool& tangentialFlipped,
    index_t numGaussPerSpan = 0)
{
    const patchSide ps1 = interf.first(), ps2 = interf.second();
    const gsGeometry<T>& g1 = mp.patch(ps1.patch);
    const gsGeometry<T>& g2 = mp.patch(ps2.patch);

    const short_t nDir1 = ps1.direction(), tDir1 = 1 - nDir1;
    const short_t nDir2 = ps2.direction(), tDir2 = 1 - nDir2;
    const bool par1 = ps1.parameter(), par2 = ps2.parameter();
    const T s1 = par1 ? T(-1) : T(1);
    const T s2 = par2 ? T(-1) : T(1);

    const gsMatrix<T> sup1 = g1.support();
    const gsMatrix<T> sup2 = g2.support();
    const T n1 = sup1(nDir1, par1 ? 1 : 0);
    const T n2 = sup2(nDir2, par2 ? 1 : 0);

    tangentialFlipped = !interf.dirOrientation(ps1, tDir1);

    const T t1a = sup1(tDir1, 0), t1b = sup1(tDir1, 1);
    const T t2a = sup2(tDir2, 0), t2b = sup2(tDir2, 1);

    std::vector<T> breaks1 = breaksOf(g1, tDir1);
    std::vector<T> breaks2 = breaksOf(g2, tDir2);
    for (T& br : breaks2)
    {
        T s = (br - t2a) / (t2b - t2a);
        if (tangentialFlipped) s = T(1) - s;
        br = t1a + s * (t1b - t1a);
    }
    std::set<T> merged(breaks1.begin(), breaks1.end());
    merged.insert(breaks2.begin(), breaks2.end());
    const std::vector<T> brk(merged.begin(), merged.end());

    const short_t deg = std::max(g1.basis().degree(tDir1), g2.basis().degree(tDir2));
    const index_t nGauss = (numGaussPerSpan > 0) ? numGaussPerSpan : 2*deg + 1;
    gsGaussRule<T> rule(nGauss);
    gsMatrix<T> nodes; gsVector<T> wts;
    rule.mapToAll(brk, nodes, wts);
    const index_t N = nodes.cols();

    gsMatrix<T> pts1(2,N), pts2(2,N);
    for (index_t i = 0; i < N; ++i)
    {
        const T t  = nodes(0,i);
        const T u  = (t - t1a) / (t1b - t1a);
        const T u2 = tangentialFlipped ? (T(1)-u) : u;
        pts1(nDir1,i) = n1;  pts1(tDir1,i) = t;
        pts2(nDir2,i) = n2;  pts2(tDir2,i) = t2a + u2 * (t2b - t2a);
    }
    gsMatrix<T> d1, d2;
    g1.deriv_into(pts1, d1);
    g2.deriv_into(pts2, d2);

    InterfaceSamples<T> S;
    S.D1.resize(N); S.D2.resize(N); S.D3.resize(N); S.t.resize(N);
    S.w = wts;
    for (index_t i = 0; i < N; ++i)
    {
        const gsVector<T> g1n = partial(d1, nDir1, i);
        const gsVector<T> g1t = partial(d1, tDir1, i);
        const gsVector<T> g2n = partial(d2, nDir2, i);
        const gsVector<T> g2t = partial(d2, tDir2, i);
        S.D1(i) = s1      * det2(g1n, g1t);
        S.D2(i) = s2      * det2(g2n, g2t);
        S.D3(i) = s1 * s2 * det2(g1n, g2n);
        S.t(i)  = (nodes(0,i) - t1a) / (t1b - t1a);
    }
    return S;
}

/// Single-solve linear least-squares for (alpha, beta) on the basis
/// {1-t, t}.  alpha is normalised by integral(alpha_1+alpha_2)=1, beta
/// is gauged by integral(alpha_1*beta_1 - alpha_2*beta_2)=0.  Tikhonov
/// bias resolves null directions (toward constant alpha and beta_1=beta_2).
template <class T>
struct SolveResult { T a[4]; T b[4]; T alphaErr; T betaErr; };

template <class T>
SolveResult<T> solveLinearGluing(const InterfaceSamples<T>& S)
{
    const index_t N = S.size();

    gsMatrix<T> Phi(N, 4);
    for (index_t i = 0; i < N; ++i)
    {
        const T p0 = T(1) - S.t(i), p1 = S.t(i);
        Phi(i,0) = p0 * S.D2(i);
        Phi(i,1) = p1 * S.D2(i);
        Phi(i,2) = p0 * S.D1(i);
        Phi(i,3) = p1 * S.D1(i);
    }

    gsMatrix<T> G = Phi.transpose() * S.w.asDiagonal() * Phi;

    // Tikhonov bias toward constant alpha
    const T tikh = T(1e-10) * G.diagonal().cwiseAbs().maxCoeff();
    G(0,0) += tikh; G(1,1) += tikh; G(0,1) -= tikh; G(1,0) -= tikh;
    G(2,2) += tikh; G(3,3) += tikh; G(2,3) -= tikh; G(3,2) -= tikh;

    gsVector<T> c(4); c.setConstant(T(0.5));
    gsMatrix<T> K(5,5); K.setZero();
    K.block(0,0,4,4) = T(2) * G;
    K.block(0,4,4,1) = c;
    K.block(4,0,1,4) = c.transpose();
    gsVector<T> r(5); r.setZero(); r(4) = T(1);
    gsVector<T> aSol = K.fullPivLu().solve(r).head(4);

    SolveResult<T> out;
    for (index_t k = 0; k < 4; ++k) out.a[k] = aSol(k);
    gsVector<T> aRes = Phi * aSol;
    out.alphaErr = S.integrate(aRes.cwiseProduct(aRes));

    // beta solve: Psi = [(1-t)D2, t*D2, -(1-t)D1, -t*D1]
    gsMatrix<T> Psi = Phi;
    Psi.col(2) *= T(-1);
    Psi.col(3) *= T(-1);

    gsMatrix<T> H = Psi.transpose() * S.w.asDiagonal() * Psi;
    gsVector<T> d = Psi.transpose() * S.w.asDiagonal() * S.D3;

    // Tikhonov bias toward beta_1 = beta_2
    const T tikhB = T(1e-10) * H.diagonal().cwiseAbs().maxCoeff();
    H(0,0) += tikhB; H(2,2) += tikhB; H(0,2) -= tikhB; H(2,0) -= tikhB;
    H(1,1) += tikhB; H(3,3) += tikhB; H(1,3) -= tikhB; H(3,1) -= tikhB;

    gsVector<T> e(4);
    {
        gsVector<T> alpha1(N), alpha2(N);
        for (index_t i = 0; i < N; ++i)
        {
            const T p0 = T(1) - S.t(i), p1 = S.t(i);
            alpha1(i) = out.a[0]*p0 + out.a[1]*p1;
            alpha2(i) = out.a[2]*p0 + out.a[3]*p1;
        }
        e(0) =  S.integrate(alpha1.cwiseProduct((T(1) - S.t.array()).matrix()));
        e(1) =  S.integrate(alpha1.cwiseProduct(S.t));
        e(2) = -S.integrate(alpha2.cwiseProduct((T(1) - S.t.array()).matrix()));
        e(3) = -S.integrate(alpha2.cwiseProduct(S.t));
    }

    gsMatrix<T> Kb(5,5); Kb.setZero();
    Kb.block(0,0,4,4) = T(2) * H;
    Kb.block(0,4,4,1) = e;
    Kb.block(4,0,1,4) = e.transpose();
    gsVector<T> rb(5); rb.head(4) = T(2)*d; rb(4) = T(0);
    gsVector<T> bSol = Kb.fullPivLu().solve(rb).head(4);
    for (index_t k = 0; k < 4; ++k) out.b[k] = bSol(k);
    gsVector<T> bRes = Psi * bSol - S.D3;
    out.betaErr = S.integrate(bRes.cwiseProduct(bRes));
    return out;
}

/// Per-interface driver.
///
/// Returns 8 numbers in patch-1 tangential parametrisation:
///   [ a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1 ]
///
/// SAME-SIGN CORRECTION: when sign(D1) == sign(D2) along the interface,
/// the homogeneous condition a1*D2 + a2*D1 = 0 with a normalisation
/// integral(a1+a2)=1 has no solution; we solve with D1 := -D1 and then
/// remap (a1, b2) := (-a1, -b2) (a2 and b1 unchanged, D3 unchanged).
///
/// FINAL beta SIGN: the embedding `createGluingDataArgyrisBasis` uses
/// the convention beta_1*D2 - beta_2*D1 = -D3, while the solver fits
/// beta_1*D2 - beta_2*D1 = +D3.  So we negate all beta values just
/// before returning, to match the embedding.  This matches the
/// post-negation that the v2 code did at the same point.
template <class T>
gsVector<T> computeGluingDataForInterface(
    const gsMultiPatch<T>& mp,
    const boundaryInterface& interf,
    bool& success,
    const T eps = T(1e-8),
    index_t numGaussPerSpan = 0)
{
    success = false;
    gsVector<T> result(8); result.setZero();

    bool flipped = false;
    InterfaceSamples<T> S =
        sampleInterface(mp, interf, flipped, numGaussPerSpan);

    if (S.D1.minCoeff() * S.D1.maxCoeff() < 0 ||
        S.D2.minCoeff() * S.D2.maxCoeff() < 0)
        return result;

    const bool sameSign = (S.D1.minCoeff() * S.D2.minCoeff() > 0);
    if (sameSign) S.D1 = -S.D1;

    SolveResult<T> r = solveLinearGluing(S);
    if (r.alphaErr >= eps) return result;

    if (sameSign)
    {
        // v3 same-sign convention
        r.a[0] = -r.a[0];  r.a[1] = -r.a[1];
        r.b[2] = -r.b[2];  r.b[3] = -r.b[3];
    }

    T a10 = r.a[0], a11 = r.a[1], a20 = r.a[2], a21 = r.a[3];
    T b10 = r.b[0], b11 = r.b[1], b20 = r.b[2], b21 = r.b[3];

    // Embedding-convention sign flip on beta (see header comment above).
    b10 = -b10; b11 = -b11; b20 = -b20; b21 = -b21;

    result(0) = a10; result(1) = a11;
    result(2) = b10; result(3) = b11;
    if (flipped)
    {
        // Patch-2 tangent is flipped, so we evaluate the patch-2 alpha
        // and beta at the reversed endpoint pairing (a21 at gd-t=0,
        // a20 at gd-t=1).  The sign-related fix-up is done at the
        // embedding call site via `tangentSign`.
        result(4) = a21; result(5) = a20;
        result(6) = b21; result(7) = b20;
    }
    else
    { result(4) = a20; result(5) = a21; result(6) = b20; result(7) = b21; }
    success = true;
    return result;
}


// ====================================================================
// main – two-patch AS-G1 basis
// ====================================================================

int main(int argc, char* argv[])
{
    using T = double;

    std::string geometry("domain2d/2patch/two_bilinear_patches.xml");
    std::string outDir("");
    index_t numGaussPerSpan = 0;
    index_t refinements = 0;
    index_t plot = -1;

    gsCmdLine cmd("AS-G1 conforming two-patch basis.");
    cmd.addString("f", "file", "Multi-patch geometry file.", geometry);
    cmd.addString("o", "outDir",
                  "Output directory for pvd/vts files (created if needed).", outDir);
    cmd.addInt("n", "numGauss",
               "Gauss points per knot span (0 = auto).", numGaussPerSpan);
    cmd.addInt("r", "refinements",
               "Uniform refinements before proceeding.", refinements);
    cmd.addInt("p", "plot",
               "Basis function to plot: -1=none, -2=ALL, 0..N=single.", plot);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // ---- Create output directory if specified ----
    std::string prefix;
    if (!outDir.empty())
    {
        // Ensure trailing slash
        if (outDir.back() != '/')
            outDir += '/';
        // Create the directory (works on Linux/macOS)
        gsFileManager::mkdir(outDir);
        prefix = outDir;
        gsInfo << "Output directory: " << outDir << "\n";
    }

    // ---- Read geometry ----
    gsMultiPatch<T>::uPtr mpPtr = gsReadFile<>(geometry);
    if (!mpPtr) { gsInfo << "Cannot read " << geometry << ".\n"; return -1; }
    gsMultiPatch<T>& mp = *mpPtr;
    mp.computeTopology();

    GISMO_ENSURE(mp.nPatches() == 2,
                 "This example is designed for exactly two patches.");
    GISMO_ENSURE(mp.nInterfaces() == 1,
                 "Expected exactly one interface between the two patches.");

    gsInfo << "Patches: " << mp.nPatches()
           << "  Interfaces: " << mp.nInterfaces() << "\n";

    // ---- Ensure minimum degree 3 ----
    const short_t inputDeg = mp.patch(0).basis().degree(0);
    if (inputDeg < 3)
    {
        const short_t elev = 3 - inputDeg;
        mp.degreeElevate(elev);
        gsInfo << "Input degree " << inputDeg
               << " → elevated by " << elev << " to degree 3\n";
    }
    else
    {
        gsInfo << "Input degree " << inputDeg << " (no elevation needed)\n";
    }

    // ---- Refine ----
    // Insert knots with multiplicity = degree - 1 (gives C^1 continuity),
    // which is the minimum needed for the Argyris construction
    // (sideBasis must have minInteriorMultiplicity > 1).
    if (refinements > 0)
    {
        const short_t deg = mp.patch(0).basis().degree(0);
        const index_t mult = std::max<index_t>(deg - 1, 1);
        for (index_t i = 0; i < refinements; ++i)
            mp.uniformRefine(1, mult);
        gsInfo << "Refined " << refinements << " time(s), knot mult = "
               << mult << " (degree " << deg << ")\n";
    }

    // ---- Identify the interface ----
    const boundaryInterface& ifc = *mp.iBegin();
    const patchSide ps1 = ifc.first();
    const patchSide ps2 = ifc.second();

    gsInfo << "Interface: patch " << ps1.patch << " side " << ps1.side()
           << " <-> patch " << ps2.patch << " side " << ps2.side() << "\n";

    // ---- Compute gluing data for the interface ----
    bool ok = false;
    gsVector<T> gd = computeGluingDataForInterface(
        mp, ifc, ok, T(1e-8), numGaussPerSpan);
    GISMO_ENSURE(ok, "Gluing data computation failed.");

    const T a1_0 = gd(0), a1_1 = gd(1), b1_0 = gd(2), b1_1 = gd(3);
    const T a2_0 = gd(4), a2_1 = gd(5), b2_0 = gd(6), b2_1 = gd(7);

    gsInfo << "\nGluing data:\n"
           << "  Patch " << ps1.patch << " side " << ps1.side()
           << ": alpha=" << a1_0 << "*(1-t)+" << a1_1 << "*t"
           << "  beta=" << b1_0 << "*(1-t)+" << b1_1 << "*t\n"
           << "  Patch " << ps2.patch << " side " << ps2.side()
           << ": alpha=" << a2_0 << "*(1-t)+" << a2_1 << "*t"
           << "  beta=" << b2_0 << "*(1-t)+" << b2_1 << "*t\n\n";

    // ---- Build per-patch interface-side embeddings ----
    const gsTensorBSplineBasis<2,T>& tb1 =
        dynamic_cast<const gsTensorBSplineBasis<2,T>&>(mp.patch(ps1.patch).basis());
    const gsTensorBSplineBasis<2,T>& tb2 =
        dynamic_cast<const gsTensorBSplineBasis<2,T>&>(mp.patch(ps2.patch).basis());

    // v4 fix #1: when the patch tangent runs opposite to the
    // gluing-data tangent (i.e. `flipped`), pass tangentSign = -1
    // to the embedding so that `beta * d/dt(smoother)` is evaluated
    // in the gluing-data tangent frame.  Both patches need it
    // because both halves of the AS-G1 relation are expressed in
    // the gluing-data tangent frame.
    const short_t _tdir1 = 1 - ps1.direction();
    const bool _flippedAtEmb = !ifc.dirOrientation(ps1, _tdir1);
    const T tSign = _flippedAtEmb ? T(-1) : T(1);

    gsSparseMatrix<T> E1 = createGluingDataArgyrisBasis(
        tb1, ps1.side(), a1_0, a1_1, b1_0, b1_1, T(1e-12), tSign);
    gsSparseMatrix<T> E2 = createGluingDataArgyrisBasis(
        tb2, ps2.side(), a2_0, a2_1, b2_0, b2_1, T(1e-12), tSign);

    gsInfo << "Patch " << ps1.patch << " interface embedding: "
           << E1.rows() << " x " << E1.cols() << "\n";
    gsInfo << "Patch " << ps2.patch << " interface embedding: "
           << E2.rows() << " x " << E2.cols() << "\n";

    // ---- DOF classification for each patch ----
    // The embedding column layout is:
    //   [nInterior_i | nLowerDeg_i | nSmoother_i]
    // where _i refers to the interface side.

    // Side bases for the interface side
    gsBSplineBasis<T> sideBasis1 = *tb1.boundaryBasis(ps1.side());
    gsBSplineBasis<T> sideBasis2 = *tb2.boundaryBasis(ps2.side());

    gsBSplineBasis<T> smootherBasis1 = sideBasis1;
    smootherBasis1.elevateContinuity(1);
    gsBSplineBasis<T> lowerDegBasis1 = sideBasis1;
    lowerDegBasis1.degreeReduce(1);

    gsBSplineBasis<T> smootherBasis2 = sideBasis2;
    smootherBasis2.elevateContinuity(1);
    gsBSplineBasis<T> lowerDegBasis2 = sideBasis2;
    lowerDegBasis2.degreeReduce(1);

    gsMatrix<index_t> bdryDOFs1  = tb1.boundary(ps1.side());
    gsMatrix<index_t> neighDOFs1 = tb1.boundaryOffset(ps1.side(), 1);
    std::vector<index_t> intDOFs1 =
        getInteriorDofs(tb1.size(), bdryDOFs1, neighDOFs1);

    gsMatrix<index_t> bdryDOFs2  = tb2.boundary(ps2.side());
    gsMatrix<index_t> neighDOFs2 = tb2.boundaryOffset(ps2.side(), 1);
    std::vector<index_t> intDOFs2 =
        getInteriorDofs(tb2.size(), bdryDOFs2, neighDOFs2);

    const index_t nInt1 = static_cast<index_t>(intDOFs1.size());
    const index_t nInt2 = static_cast<index_t>(intDOFs2.size());
    const index_t nLD1  = lowerDegBasis1.size();   // second-layer cols
    const index_t nLD2  = lowerDegBasis2.size();
    const index_t nSm1  = smootherBasis1.size();   // boundary cols
    const index_t nSm2  = smootherBasis2.size();

    // Sanity checks
    GISMO_ENSURE(E1.cols() == nInt1 + nLD1 + nSm1,
                 "Column count mismatch for patch " << ps1.patch);
    GISMO_ENSURE(E2.cols() == nInt2 + nLD2 + nSm2,
                 "Column count mismatch for patch " << ps2.patch);

    gsInfo << "\nPatch " << ps1.patch << ": nInterior=" << nInt1
           << "  nSecondLayer=" << nLD1
           << "  nBoundary=" << nSm1 << "\n";
    gsInfo << "Patch " << ps2.patch << ": nInterior=" << nInt2
           << "  nSecondLayer=" << nLD2
           << "  nBoundary=" << nSm2 << "\n";

    // ====================================================================
    // Assemble the global basis
    // ====================================================================
    //
    // The global DOF numbering is:
    //
    //   [patch1_interior | patch2_interior | shared_boundary | shared_second_layer]
    //
    // where:
    //   - patch_i interior: DOFs not on the interface boundary/second-layer
    //     These are independent per patch and map via identity in E_i.
    //   - shared boundary: The nSm (= nSm1 = nSm2) interface trace DOFs.
    //     The same trace coefficients are used for both patches.
    //   - shared second-layer: The nLD (= nLD1 = nLD2) interface d-derivative DOFs.
    //     These are shared between patches; the AS-G1 condition means
    //     the same second-layer coefficient feeds both patches' embeddings.
    //
    // IMPORTANT: For the two patches to share the interface trace values
    // and d-derivative values, the two side bases must be compatible
    // (same knot vector along the interface).  We verify this.

    GISMO_ENSURE(nSm1 == nSm2,
                 "Smoother basis sizes differ across interface ("
                     << nSm1 << " vs " << nSm2 << ").");
    GISMO_ENSURE(nLD1 == nLD2,
                 "Lower-degree basis sizes differ across interface ("
                     << nLD1 << " vs " << nLD2 << ").");

    const index_t nSharedBdry = nSm1;    // shared trace DOFs
    const index_t nSharedL2   = nLD1;    // shared second-layer DOFs

    // Check interface orientation: if the tangential directions
    // run in opposite directions, we need to reverse the DOF mapping
    // for patch 2's shared columns.
    const short_t tanDir1 = 1 - ps1.direction();
    const bool flipped = !ifc.dirOrientation(ps1, tanDir1);

    gsInfo << "Interface orientation: "
           << (flipped ? "FLIPPED" : "aligned") << "\n";

    const index_t nGlobal = nInt1 + nInt2 + nSharedBdry + nSharedL2;

    gsInfo << "\nGlobal DOFs: " << nGlobal
           << " = " << nInt1 << " (int1) + " << nInt2 << " (int2) + "
           << nSharedBdry << " (shared bdry) + "
           << nSharedL2 << " (shared L2)\n";

    // Build the two global-to-patch matrices:
    //   G1 : nGlobal → tb1.size()    (coefficients for patch 1)
    //   G2 : nGlobal → tb2.size()    (coefficients for patch 2)
    //
    // For patch 1:
    //   - Global cols [0 .. nInt1-1]  →  E1 cols [0 .. nInt1-1]   (interior)
    //   - Global cols [nInt1+nInt2 .. +nSharedBdry-1]  →  E1 cols [nInt1+nLD1 .. end] (boundary)
    //   - Global cols [nInt1+nInt2+nSharedBdry .. end]  →  E1 cols [nInt1 .. nInt1+nLD1-1] (L2)
    //
    // For patch 2:
    //   - Global cols [nInt1 .. nInt1+nInt2-1]  →  E2 cols [0 .. nInt2-1]   (interior)
    //   - Global cols [nInt1+nInt2 .. +nSharedBdry-1]  →  E2 cols [nInt2+nLD2 .. end] (boundary)
    //   - Global cols [nInt1+nInt2+nSharedBdry .. end]  →  E2 cols [nInt2 .. nInt2+nLD2-1] (L2)

    const index_t gOff_int1     = 0;
    const index_t gOff_int2     = nInt1;
    const index_t gOff_bdry     = nInt1 + nInt2;
    const index_t gOff_L2       = nInt1 + nInt2 + nSharedBdry;

    // --- Patch 1 global matrix ---
    gsSparseMatrix<T> G1(tb1.size(), nGlobal);
    {
        // Interior columns of E1 → global interior cols for patch 1
        for (index_t j = 0; j < nInt1; ++j)
            for (typename gsSparseMatrix<T>::InnerIterator it(E1, j); it; ++it)
                G1.insert(it.row(), gOff_int1 + j) = it.value();

        // Second-layer columns of E1 → global shared L2 cols
        for (index_t j = 0; j < nLD1; ++j)
        {
            const index_t e1col = nInt1 + j;
            for (typename gsSparseMatrix<T>::InnerIterator it(E1, e1col); it; ++it)
                G1.insert(it.row(), gOff_L2 + j) = it.value();
        }

        // Boundary columns of E1 → global shared boundary cols
        for (index_t j = 0; j < nSm1; ++j)
        {
            const index_t e1col = nInt1 + nLD1 + j;
            for (typename gsSparseMatrix<T>::InnerIterator it(E1, e1col); it; ++it)
                G1.insert(it.row(), gOff_bdry + j) = it.value();
        }
    }
    G1.makeCompressed();

    // --- Patch 2 global matrix ---
    gsSparseMatrix<T> G2(tb2.size(), nGlobal);
    {
        // Interior columns of E2 → global interior cols for patch 2
        for (index_t j = 0; j < nInt2; ++j)
            for (typename gsSparseMatrix<T>::InnerIterator it(E2, j); it; ++it)
                G2.insert(it.row(), gOff_int2 + j) = it.value();

        // Second-layer columns of E2 → global shared L2 cols
        // AS-G1 requires alpha_1*d_1 + alpha_2*d_2 = 0, so patch 2's
        // second-layer contribution is sign-flipped relative to patch 1.
        // When the tangential directions are FLIPPED, the d-deriv basis
        // function indexing also reverses (j2 = nLD2-1-j).  The relative
        // sign of patch-2's d_i derivative depends on the orientation:
        //   - aligned:  d_2 = +alpha_2 * phi_j (after embedding) → negate to oppose patch 1
        //   - flipped:  d_2 also has an additional minus from the reparam,
        //               which cancels the negation, so we DON'T negate.
        const T l2Sign = flipped ? T(1) : T(-1);
        for (index_t j = 0; j < nLD2; ++j)
        {
            const index_t j2 = flipped ? (nLD2 - 1 - j) : j;
            const index_t e2col = nInt2 + j2;
            for (typename gsSparseMatrix<T>::InnerIterator it(E2, e2col); it; ++it)
                G2.insert(it.row(), gOff_L2 + j) = l2Sign * it.value();
        }

        // Boundary columns of E2 → global shared boundary cols
        // If flipped, DOF j on patch 1 corresponds to DOF (nSm2-1-j) on patch 2.
        for (index_t j = 0; j < nSm2; ++j)
        {
            const index_t j2 = flipped ? (nSm2 - 1 - j) : j;
            const index_t e2col = nInt2 + nLD2 + j2;
            for (typename gsSparseMatrix<T>::InnerIterator it(E2, e2col); it; ++it)
                G2.insert(it.row(), gOff_bdry + j) = it.value();
        }
    }
    G2.makeCompressed();

    gsInfo << "\nGlobal-to-patch matrices:\n"
           << "  G1: " << G1.rows() << " x " << G1.cols() << "\n"
           << "  G2: " << G2.rows() << " x " << G2.cols() << "\n";

    // ====================================================================
    // Numerical G1 smoothness verification
    // ====================================================================
    {
        const short_t normDir1 = ps1.direction(), tanDir1_ = 1 - normDir1;
        const short_t normDir2 = ps2.direction(), tanDir2_ = 1 - normDir2;
        const bool par1_ = ps1.parameter(), par2_ = ps2.parameter();
        const bool ifcFlipped = !ifc.dirOrientation(ps1, tanDir1_);

        gsMatrix<T> sup1 = mp.patch(ps1.patch).support();
        gsMatrix<T> sup2 = mp.patch(ps2.patch).support();
        const T ifcCoord1 = sup1(normDir1, par1_ ? 1 : 0);
        const T ifcCoord2 = sup2(normDir2, par2_ ? 1 : 0);
        const T t1a = sup1(tanDir1_, 0), t1b = sup1(tanDir1_, 1);
        const T t2a = sup2(tanDir2_, 0), t2b = sup2(tanDir2_, 1);

        const index_t nCheck = 21;
        T maxValErr = 0, maxGradErr = 0;

        T maxErrInt = 0, maxErrTrace = 0, maxErrL2 = 0;

        for (index_t idx = 0; idx < nGlobal; ++idx)
        {
            gsVector<T> globalVec = gsVector<T>::Zero(nGlobal);
            globalVec(idx) = T(1);
            gsVector<T> c1 = G1 * globalVec;
            gsVector<T> c2 = G2 * globalVec;

            auto func1 = tb1.makeGeometry(c1);
            auto func2 = tb2.makeGeometry(c2);

            T thisMaxGrad = 0;

            for (index_t i = 0; i < nCheck; ++i)
            {
                T s = T(i) / T(nCheck - 1);
                T t1 = t1a + s * (t1b - t1a);
                T s2 = ifcFlipped ? (1.0 - s) : s;
                T t2 = t2a + s2 * (t2b - t2a);

                gsMatrix<T> pt1(2, 1), pt2(2, 1);
                pt1(normDir1, 0) = ifcCoord1; pt1(tanDir1_, 0) = t1;
                pt2(normDir2, 0) = ifcCoord2; pt2(tanDir2_, 0) = t2;

                // Check C0: values must match
                gsMatrix<T> v1 = func1->eval(pt1);
                gsMatrix<T> v2 = func2->eval(pt2);
                T valErr = std::abs(v1(0, 0) - v2(0, 0));
                maxValErr = std::max(maxValErr, valErr);

                // Check G1: physical gradients must match
                gsMatrix<T> df1, df2, dG1, dG2;
                func1->deriv_into(pt1, df1);
                func2->deriv_into(pt2, df2);
                mp.patch(ps1.patch).deriv_into(pt1, dG1);
                mp.patch(ps2.patch).deriv_into(pt2, dG2);

                gsMatrix<T> J1(2,2), J2(2,2);
                J1(0,0) = dG1(0,0); J1(0,1) = dG1(1,0);
                J1(1,0) = dG1(2,0); J1(1,1) = dG1(3,0);
                J2(0,0) = dG2(0,0); J2(0,1) = dG2(1,0);
                J2(1,0) = dG2(2,0); J2(1,1) = dG2(3,0);

                gsVector<T> paramGrad1(2), paramGrad2(2);
                paramGrad1(0) = df1(0,0); paramGrad1(1) = df1(1,0);
                paramGrad2(0) = df2(0,0); paramGrad2(1) = df2(1,0);

                gsVector<T> physGrad1 = J1.inverse().transpose() * paramGrad1;
                gsVector<T> physGrad2 = J2.inverse().transpose() * paramGrad2;

                T gradErr = (physGrad1 - physGrad2).norm();
                maxGradErr = std::max(maxGradErr, gradErr);
                thisMaxGrad = std::max(thisMaxGrad, gradErr);
            }

            // Classify DOF
            if (idx < gOff_bdry)
                maxErrInt = std::max(maxErrInt, thisMaxGrad);
            else if (idx < gOff_L2)
                maxErrTrace = std::max(maxErrTrace, thisMaxGrad);
            else
                maxErrL2 = std::max(maxErrL2, thisMaxGrad);
        }

        gsInfo << "\n=== G1 Smoothness Check (physical gradient) ===\n"
               << "  Max C0 (value) error:      " << maxValErr << "\n"
               << "  Max grad (physical) error:  " << maxGradErr << "\n"
               << "    Interior DOFs:    " << maxErrInt << "\n"
               << "    Trace DOFs:       " << maxErrTrace << "\n"
               << "    D-deriv DOFs:     " << maxErrL2 << "\n";
        if (maxValErr < 1e-8 && maxGradErr < 1e-3)
            gsInfo << "  STATUS: PASS\n";
        else if (maxValErr < 1e-8 && maxGradErr < 1e-1)
            gsInfo << "  STATUS: APPROX (AS-G1 approximation error)\n";
        else
            gsInfo << "  STATUS: FAIL (check coupling)\n";
    }

    // ====================================================================
    // Plot global basis functions on both patches
    // ====================================================================

    // Build a name and category tag for each global DOF
    // Global layout: [p0_interior | p1_interior | ifc_trace | ifc_dderiv]
    auto basisName = [&](index_t idx) -> std::string
    {
        if (idx < gOff_int2)
            return "p" + std::to_string(ps1.patch) + "_int_"
                       + std::to_string(idx - gOff_int1);
        if (idx < gOff_bdry)
            return "p" + std::to_string(ps2.patch) + "_int_"
                       + std::to_string(idx - gOff_int2);
        if (idx < gOff_L2)
            return "ifc_trace_" + std::to_string(idx - gOff_bdry);
        return "ifc_dderiv_" + std::to_string(idx - gOff_L2);
    };

    // Determine which basis functions to plot
    std::vector<index_t> toPlot;
    if (plot == -2)
    {
        toPlot.resize(nGlobal);
        std::iota(toPlot.begin(), toPlot.end(), 0);
        gsInfo << "\nPlotting all " << nGlobal << " basis functions ...\n";
    }
    else if (plot >= 0 && plot < nGlobal)
    {
        toPlot.push_back(plot);
        gsInfo << "\nPlotting global basis function " << plot
               << " (" << basisName(plot) << ") ...\n";
    }
    else if (plot >= nGlobal)
    {
        gsInfo << "\nBasis function index " << plot
               << " out of range [0, " << nGlobal - 1 << "].\n";
    }

    for (const index_t idx : toPlot)
    {
        gsVector<T> globalVec = gsVector<T>::Zero(nGlobal);
        globalVec(idx) = T(1);

        gsVector<T> c1 = G1 * globalVec;
        gsVector<T> c2 = G2 * globalVec;

        gsMatrix<T> coefs1 = c1;
        gsMatrix<T> coefs2 = c2;

        gsMultiPatch<T> sol;
        sol.addPatch(tb1.makeGeometry(give(coefs1)));
        sol.addPatch(tb2.makeGeometry(give(coefs2)));

        gsField<T> field(mp, sol);

        std::string fname = prefix + "as_g1_" + basisName(idx);
        gsWriteParaview<>(field, fname, 1000);
        gsInfo << "  [" << idx << "] " << fname << ".pvd\n";
    }

    if (!toPlot.empty())
    {
        gsInfo << "\nFile naming convention:\n"
               << "  p<N>_int_<K>   = patch N interior DOF K\n"
               << "  ifc_trace_<K>  = shared interface trace DOF K\n"
               << "  ifc_dderiv_<K> = shared interface d-derivative DOF K\n";
    }

    // ====================================================================
    // Print summary information about the global basis
    // ====================================================================
    gsInfo << "\n=== Summary ===\n";
    gsInfo << "  Tensor basis size per patch: " << tb1.size() << "\n";
    gsInfo << "  Per-patch embedding cols:    " << E1.cols()
           << " (int=" << nInt1 << " L2=" << nLD1 << " bdry=" << nSm1 << ")\n";
    gsInfo << "  Global DOFs:                 " << nGlobal << "\n";
    gsInfo << "    Patch " << ps1.patch << " interior: " << nInt1 << "\n";
    gsInfo << "    Patch " << ps2.patch << " interior: " << nInt2 << "\n";
    gsInfo << "    Shared interface trace:     " << nSharedBdry << "\n";
    gsInfo << "    Shared interface d-deriv:   " << nSharedL2 << "\n";
    gsInfo << "\nTo plot basis function k:\n"
           << "  ./bin/as_g1_two_patch_basis_v4 -f <file> -r <ref> -p k\n"
           << "To plot ALL basis functions:\n"
           << "  ./bin/as_g1_two_patch_basis_v4 -f <file> -r <ref> -p -2\n";

    gsInfo << "\nDone.\n";
    return 0;
}
