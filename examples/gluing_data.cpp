/** @file gluing_data_v4.cpp

    @brief Compute AS-G1 gluing data for multi-patch geometries (v4).

    ## v4 changes relative to v3

    The gluing-data **computation** is mathematically unchanged from v3:
    we solve exactly the same linear least-squares system for
    `(alpha_k, beta_k)`.  v4 only differs in:

      * Extended documentation explaining the tangent-frame convention
        and what consumers (Argyris embeddings) must do when the
        interface orientation flag `flipped` is true.
      * The per-interface output now prints the recommended
        `tangentSign` value that consumers should pass to their
        `createGluingDataArgyrisBasis(..., tangentSign)` call.

    ## What this computes

    For every interface of a 2D multi-patch geometry, finds linear
    functions alpha_k(t), beta_k(t) (k=1,2) along the shared edge that
    satisfy the AS-G1 conditions

        alpha_1(t) * D2(t)  +  alpha_2(t) * D1(t)  =  0,           (G0)
        beta_1(t)  * D2(t)  -  beta_2(t)  * D1(t)  =  D3(t),       (G1)

    with the normalisation  integral(alpha_1 + alpha_2) dt = 1,
    where

        D1 = s1 * det( dG1/dn , dG1/dt )
        D2 = s2 * det( dG2/dn , dG2/dt )
        D3 = s1*s2 * det( dG1/dn , dG2/dn )

    and s_k = (-1)^{par_k} corrects for the side parameter.

    ## Tangent-frame convention (important!)

    The parameter `t` runs along the **patch-1 natural tangent** mapped
    to [0,1].  Patch 2's samples are aligned to this parameter by
    reversing patch 2's natural tangent whenever
    `flipped = !ifc.dirOrientation(ps1, tDir1)` is true.  Therefore:

      * `(alpha_1, beta_1)` are linear coefficients in patch 1's own
        natural tangent (always the same as the gluing-data tangent).
      * `(alpha_2, beta_2)` are linear coefficients in the **gluing-data
        tangent**, not in patch 2's natural tangent.  When `flipped`
        is true, patch 2's natural tangent runs opposite to the
        gluing-data tangent.

    A consumer that builds the Argyris-type trace embedding
    `createGluingDataArgyrisBasis(tb, side, a0, a1, b0, b1,
                                  eps, tangentSign)` must therefore
    pass `tangentSign = flipped ? -1 : +1` for **both** patches when
    the gluing data was computed in the gluing-data frame as above.

    ## Return convention
    Per interface, 8 numbers in patch-1 tangential parametrisation:
        [ a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1 ]
    where _0/_1 are values at t=0 / t=1.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <iostream>
#include <iomanip>
#include <set>
#include <algorithm>
#include <gismo.h>

using namespace gismo;

// ====================================================================
// Tiny helpers
// ====================================================================

/// 2D determinant of two column vectors.
template <class T>
inline T det2(const gsVector<T>& a, const gsVector<T>& b)
{ return a(0)*b(1) - a(1)*b(0); }

/// Pick the partial derivative d_dir(f) at column `col` from the matrix
/// returned by gsGeometry::deriv_into (layout: row = dir + 2*comp).
template <class T>
inline gsVector<T> partial(const gsMatrix<T>& d, index_t dir, index_t col)
{
    gsVector<T> v(2);
    v(0) = d(dir,     col);
    v(1) = d(dir + 2, col);
    return v;
}

/// Knot-span breaks of the 1D component basis in direction `dir`.
template <class T>
std::vector<T> breaksOf(const gsGeometry<T>& g, short_t dir)
{
    const gsBSplineBasis<T>* bb =
        dynamic_cast<const gsBSplineBasis<T>*>(&g.basis().component(dir));
    GISMO_ENSURE(bb, "gluing_data_v4 currently supports B-spline bases only.");
    return bb->knots().breaks();
}

// ====================================================================
// Sampled interface data
// ====================================================================
//
// Everything needed to solve (G0)+(G1) at quadrature points along the
// interface, expressed in patch-1's tangential parameter mapped to [0,1].
//
template <class T>
struct InterfaceSamples
{
    gsVector<T> D1, D2, D3;   // sampled determinants
    gsVector<T> t;            // normalised parameter in [0,1]
    gsVector<T> w;            // Gauss weights (already on physical span)

    index_t size() const { return D1.size(); }

    /// integral( f )  over the interface
    T integrate(const gsVector<T>& f) const { return w.dot(f); }
};

// ====================================================================
// Single linear least-squares solve for (alpha, beta)
// ====================================================================
//
// Step A.  Solve for alpha = (a1_0, a1_1, a2_0, a2_1) by
//          minimising  || a1*D2 + a2*D1 ||^2  on the linear hat basis
//          {1-t, t}, with constraint  integral(a1+a2) = 1.
//
// Step B.  With alpha fixed, solve for beta = (b1_0, b1_1, b2_0, b2_1)
//          by minimising  || b1*D2 - b2*D1 - D3 ||^2  with gauge
//          constraint  integral(a1*b1 - a2*b2) = 0.
//
// Returns the residuals (alphaErr, betaErr).  Success is judged by the
// caller (typically alphaErr < eps).
//
template <class T>
struct SolveResult { T a[4]; T b[4]; T alphaErr; T betaErr; };

template <class T>
SolveResult<T> solveLinearGluing(const InterfaceSamples<T>& S)
{
    const index_t N = S.size();

    // ----- Build basis matrix Phi = [ (1-t)*D2 , t*D2 , (1-t)*D1 , t*D1 ]
    gsMatrix<T> Phi(N, 4);
    for (index_t i = 0; i < N; ++i)
    {
        const T p0 = T(1) - S.t(i), p1 = S.t(i);
        Phi(i,0) = p0 * S.D2(i);
        Phi(i,1) = p1 * S.D2(i);
        Phi(i,2) = p0 * S.D1(i);
        Phi(i,3) = p1 * S.D1(i);
    }

    // ----- Weighted Gram matrix G_{jk} = integral( Phi_j * Phi_k )
    gsMatrix<T> G = Phi.transpose() * S.w.asDiagonal() * Phi;

    // ----- Tikhonov bias toward constant alpha (resolves the null
    //       direction when D2/D1 is constant; otherwise negligible).
    //       Penalises  (a1_1 - a1_0)^2 + (a2_1 - a2_0)^2 .
    const T tikh = T(1e-10) * G.diagonal().cwiseAbs().maxCoeff();
    G(0,0) += tikh; G(1,1) += tikh; G(0,1) -= tikh; G(1,0) -= tikh;
    G(2,2) += tikh; G(3,3) += tikh; G(2,3) -= tikh; G(3,2) -= tikh;

    // ----- Constraint integral(a1+a2)=1 :  c^T x = 1, c = (1/2,1/2,1/2,1/2)
    gsVector<T> c(4); c.setConstant(T(0.5));

    // ----- KKT system for alpha
    gsMatrix<T> K(5,5); K.setZero();
    K.block(0,0,4,4) = T(2) * G;
    K.block(0,4,4,1) = c;
    K.block(4,0,1,4) = c.transpose();
    gsVector<T> r(5); r.setZero(); r(4) = T(1);

    gsVector<T> aSol = K.fullPivLu().solve(r).head(4);

    SolveResult<T> out;
    for (index_t k = 0; k < 4; ++k) out.a[k] = aSol(k);

    // alpha residual
    gsVector<T> aRes = Phi * aSol;
    out.alphaErr = S.integrate(aRes.cwiseProduct(aRes));

    // ----- Build basis Psi for beta: same as Phi but with a sign flip
    //       on the D1 columns, because the beta equation is  b1*D2 - b2*D1 = D3.
    gsMatrix<T> Psi = Phi;
    Psi.col(2) *= T(-1);
    Psi.col(3) *= T(-1);

    gsMatrix<T> H = Psi.transpose() * S.w.asDiagonal() * Psi;
    gsVector<T> d = Psi.transpose() * S.w.asDiagonal() * S.D3;

    // Tikhonov bias toward beta_1 = beta_2 (resolves the null direction
    // when alpha_1 = alpha_2 const; otherwise negligible).
    //   penalises (b1_0 - b2_0)^2 + (b1_1 - b2_1)^2
    const T tikhB = T(1e-10) * H.diagonal().cwiseAbs().maxCoeff();
    H(0,0) += tikhB; H(2,2) += tikhB; H(0,2) -= tikhB; H(2,0) -= tikhB;
    H(1,1) += tikhB; H(3,3) += tikhB; H(1,3) -= tikhB; H(3,1) -= tikhB;

    // gauge constraint  integral( a1*b1 - a2*b2 ) = 0
    //   e_k = integral( alpha_k(t) * basis_k(t) ),  with sign +1 for b1, -1 for b2
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

    // beta residual
    gsVector<T> bRes = Psi * bSol - S.D3;
    out.betaErr = S.integrate(bRes.cwiseProduct(bRes));

    return out;
}

// ====================================================================
// Sample D1,D2,D3 along an interface
// ====================================================================
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

    // Sign factors (parametrisation to canonical "west")
    const T s1 = par1 ? T(-1) : T(1);
    const T s2 = par2 ? T(-1) : T(1);

    const gsMatrix<T> sup1 = g1.support();
    const gsMatrix<T> sup2 = g2.support();
    const T n1 = sup1(nDir1, par1 ? 1 : 0);   // normal coord on patch 1
    const T n2 = sup2(nDir2, par2 ? 1 : 0);   // normal coord on patch 2

    tangentialFlipped = !interf.dirOrientation(ps1, tDir1);

    // ----- Merge knot-span breaks from both sides into patch-1 tangent range
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

    // ----- Gauss quadrature on the merged knot vector
    const short_t deg = std::max(g1.basis().degree(tDir1), g2.basis().degree(tDir2));
    const index_t nGauss = (numGaussPerSpan > 0) ? numGaussPerSpan : 2*deg + 1;
    gsGaussRule<T> rule(nGauss);
    gsMatrix<T> nodes; gsVector<T> wts;
    rule.mapToAll(brk, nodes, wts);
    const index_t N = nodes.cols();

    // ----- Evaluate derivatives at the quadrature points
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

    // ----- Assemble samples
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

// ====================================================================
// Per-interface driver: sample, solve, sign-correct, pack
// ====================================================================
//
// Returns 8 numbers in patch-1 tangential parametrisation:
//     [ a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1 ]
//
// SAME-SIGN CORRECTION
// --------------------
// If D1 and D2 have the same sign, the homogeneous condition (G0)
// cannot hold with both alpha_k > 0 and the integral normalisation
// integral(a1+a2)=1.  We instead set D1 := -D1 and solve; the result
// satisfies
//     a1' * D2 + a2' * (-D1) = 0          and     integral(a1'+a2') = 1.
// Translating back to the original D1:
//     a1 = -a1',  a2 = a2'   (matches a1*D2 + a2*D1 = 0)
//     b1 =  b1',  b2 = -b2'  (matches b1*D2 - b2*D1 = D3, with D3
//                             left untouched -- self-consistent choice)
// Note: after this remap, integral(a1+a2) is generally NOT 1 (e.g. 0 in
// the symmetric case).  The downstream beta gauge constraint
// integral(a1*b1 - a2*b2)=0 is HOMOGENEOUS in (alpha,beta), so any
// uniform rescale of both leaves it satisfied; we therefore keep the
// "unit-magnitude" values produced by the same-sign branch instead of
// blowing them up by 1/integral(a1+a2).  Downstream consumers of this
// data should be aware that the gauge integral is not pinned in this
// branch.
//
template <class T>
gsVector<T> computeGluingDataForInterface(
    const gsMultiPatch<T>& mp,
    const boundaryInterface& interf,
    bool& success,
    T eps = T(1e-8),
    index_t numGaussPerSpan = 0)
{
    gsVector<T> result(8); result.setZero();
    success = false;

    bool tFlip = false;
    InterfaceSamples<T> S = sampleInterface(mp, interf, tFlip, numGaussPerSpan);

    if (S.D1.minCoeff()*S.D1.maxCoeff() < 0 ||
        S.D2.minCoeff()*S.D2.maxCoeff() < 0)
    {
        gsInfo << "  WARNING: D1 or D2 changes sign on the interface.\n";
        return result;
    }

    const bool sameSign = (S.D1.minCoeff() * S.D2.minCoeff() > 0);
    if (sameSign) S.D1 = -S.D1;     // see comment block above

    SolveResult<T> r = solveLinearGluing(S);
    if (r.alphaErr >= eps)
    {
        gsInfo << "  FAILED  alphaErr=" << r.alphaErr << "\n";
        return result;
    }

    // Apply same-sign correction (no gauge rescale -- see header comment)
    if (sameSign)
    {
        r.a[0] = -r.a[0];  r.a[1] = -r.a[1];
        r.b[2] = -r.b[2];  r.b[3] = -r.b[3];
    }

    gsInfo << std::fixed << std::setprecision(4)
           << "  alpha_1(t) = " << r.a[0] << "*(1-t) + " << r.a[1] << "*t\n"
           << "  alpha_2(t) = " << r.a[2] << "*(1-t) + " << r.a[3] << "*t\n"
           << "  beta_1(t)  = " << r.b[0] << "*(1-t) + " << r.b[1] << "*t\n"
           << "  beta_2(t)  = " << r.b[2] << "*(1-t) + " << r.b[3] << "*t\n"
           << "  alphaErr="    << r.alphaErr
           << "  betaErr="     << r.betaErr
           << (sameSign ? "  [same-sign branch]\n" : "\n");

    // Pack into [a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1].
    // If the second patch traverses the interface in the opposite
    // direction (tFlip), swap the endpoints of its alpha/beta.
    result(0) = r.a[0]; result(1) = r.a[1];
    result(2) = r.b[0]; result(3) = r.b[1];
    if (tFlip) {
        result(4) = r.a[3]; result(5) = r.a[2];
        result(6) = r.b[3]; result(7) = r.b[2];
    } else {
        result(4) = r.a[2]; result(5) = r.a[3];
        result(6) = r.b[2]; result(7) = r.b[3];
    }
    success = true;
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
gsMatrix<T> computeGluingDataMatrix(
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
        gsInfo << "\nInterface " << interf << "\n";

        // Detect orientation for the tangentSign hint
        const patchSide _ps1 = interf.first();
        const short_t  _tDir1 = 1 - _ps1.direction();
        const bool     _flipped = !interf.dirOrientation(_ps1, _tDir1);
        gsInfo << "  Orientation: " << (_flipped ? "FLIPPED" : "aligned") << "\n"
               << "  Embedding consumers should use tangentSign = "
               << (_flipped ? "-1" : "+1") << " for BOTH patches.\n";

        bool ok = false;
        gsVector<T> gd = computeGluingDataForInterface(
            mp, interf, ok, eps, numGaussPerSpan);
        if (!ok) continue;

        const index_t c1 = sideCol(interf.first());
        const index_t c2 = sideCol(interf.second());
        M.row(interf.first().patch ).segment(c1, 4) = gd.head(4).transpose();
        M.row(interf.second().patch).segment(c2, 4) = gd.tail(4).transpose();
    }
    return M;
}

// ====================================================================
// Entry point
// ====================================================================
int main(int argc, char* argv[])
{
    using T = double;

    std::string geometry("domain2d/2patch/two_bicubic_patches.xml");
    index_t numGaussPerSpan = 0;
    real_t eps = 1e-8;

    gsCmdLine cmd("Compute AS-G1 gluing data (v4 — adds tangentSign hint).");
    cmd.addString("f", "file",     "Input multi-patch file.",        geometry);
    cmd.addInt   ("n", "numGauss", "Gauss points per span (0=auto).", numGaussPerSpan);
    cmd.addReal  ("e", "eps",      "Residual tolerance.",             eps);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
    if (!mpPtr) { gsInfo << "Cannot read " << geometry << ".\n"; return -1; }
    gsMultiPatch<>& mp = *mpPtr;
    mp.computeTopology();

    gsInfo << "Patches: "    << mp.nPatches()
           << "  Interfaces: " << mp.nInterfaces() << "\n";

    gsMatrix<T> M = computeGluingDataMatrix<T>(mp, T(eps), numGaussPerSpan);

    static const char* sideNames[4] = {"west", "east", "south", "north"};
    gsInfo << "\n=== Gluing data matrix (" << M.rows() << " x " << M.cols() << ") ===\n";
    for (index_t p = 0; p < M.rows(); ++p) {
        gsInfo << "Patch " << p << ":\n";
        for (index_t s = 0; s < 4; ++s)
            gsInfo << "  " << std::setw(5) << sideNames[s]
                   << ":  alpha=" << M(p,4*s+0) << "*(1-t)+" << M(p,4*s+1) << "*t"
                   << "   beta="  << M(p,4*s+2) << "*(1-t)+" << M(p,4*s+3) << "*t\n";
    }
    return 0;
}
