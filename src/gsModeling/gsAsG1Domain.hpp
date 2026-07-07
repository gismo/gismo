/** @file gsAsG1Domain.hpp

    @brief Derive gluing data for an AS-G1 geometry.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

namespace gismo {

// ====================================================================
// Gluing-data helpers
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
    const gsBSplineBasis<T>& bb =
        dynamic_cast<const gsBSplineBasis<T>&>(geo.basis().component(dir));
    return bb.knots().breaks();
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
    const T eps = T(1e-8),
    index_t numGaussPerSpan = 0)
{
    gsVector<T> result(8); result.setZero();

    bool flipped = false;
    InterfaceSamples<T> S =
        sampleInterface(mp, interf, flipped, numGaussPerSpan);

    GISMO_ENSURE (S.D1.minCoeff() * S.D1.maxCoeff() > 0 &&
        S.D2.minCoeff() * S.D2.maxCoeff() > 0, "Not AS-G1");

    const bool sameSign = (S.D1.minCoeff() * S.D2.minCoeff() > 0);
    if (sameSign) S.D1 = -S.D1;

    SolveResult<T> r = solveLinearGluing(S);
    GISMO_ENSURE (r.alphaErr <= eps, "Not AS-G1");

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

    if (flipped)
    {
        // Patch-2 tangent is flipped, so we evaluate the patch-2 alpha
        // and beta at the reversed endpoint pairing (a21 at gd-t=0,
        // a20 at gd-t=1).  The sign-related fix-up is done at the
        // embedding call site via `tangentSign`.
        result(0) = a10; result(1) = a11;
        result(2) = b10; result(3) = b11;
        result(4) = a21; result(5) = a20;
        result(6) = b21; result(7) = b20;
    }
    else
    {
        result(0) = a10; result(1) = a11;
        result(2) = b10; result(3) = b11;
        result(4) = a20; result(5) = a21;
        result(6) = b20; result(7) = b21;
    }
    return result;
}


// ====================================================================
// Multi-patch driver: assemble the (nPatches x 16) matrix
//   columns per side s in {west=0, east=1, south=2, north=3}:
//     4*s+0 = alpha_0    4*s+1 = alpha_1
//     4*s+2 = beta_0     4*s+3 = beta_1
// Boundary sides keep their default (alpha=1, beta=0).
// ====================================================================
template <class T>
gsMatrix<T> computeGluingData(
    const gsMultiPatch<T>& mp, T eps = T(1e-8), index_t numGaussPerSpan = 0)
{
    const index_t nP = mp.nPatches();
    gsMatrix<T> M(nP, 16);
    for (index_t p = 0; p < nP; ++p)
        for (index_t s = 0; s < 4; ++s) {
            M(p, 4*s+0) = T(1); M(p, 4*s+1) = T(1);
            M(p, 4*s+2) = T(0); M(p, 4*s+3) = T(0);
        }

    auto sideCol = [](const patchSide& ps) -> index_t
        { return 4 * (static_cast<index_t>(ps.side()) - 1); };

    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const boundaryInterface& interf = *it;

        gsVector<T> gd = computeGluingDataForInterface(
            mp, interf, eps, numGaussPerSpan);

        const index_t c1 = sideCol(interf.first());
        const index_t c2 = sideCol(interf.second());
        M.row(interf.first().patch ).segment(c1, 4) = gd.head(4).transpose();
        M.row(interf.second().patch).segment(c2, 4) = gd.tail(4).transpose();
    }
    return M;
}


} // namespace gismo
