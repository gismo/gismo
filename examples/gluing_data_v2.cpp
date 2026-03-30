/** @file gluing_data_v2.cpp

    @brief Compute AS-G1 gluing data for multi-patch geometries.

    ## Overview
    This file computes the gluing data (alpha and beta functions) needed
    for AS (Approximate Smooth) G1 constructions on multi-patch spline
    geometries.  In a G1 construction, two patches sharing an interface
    must satisfy compatibility conditions that involve linear functions
    alpha and beta along each interface edge.

    ## Algorithm
    For each interface, we identify the normal and tangential directions
    from the patchSide info and evaluate derivatives directly on the
    original (unmodified) patches -- no reparametrisation is performed.
    Three determinant quantities are evaluated at Gauss quadrature
    points along the interface:

      D1(t) = s1 * det( d_norm G1, d_tan G1 )
      D2(t) = s2 * det( d_norm G2, d_tan G2 )
      D3(t) = s1*s2 * det( d_norm G1, d_norm G2 )

    where s_k = (-1)^{par_k} is a sign factor that ensures consistency
    with the canonical (reparametrised-to-west) convention.

    The alpha functions (a1, a2) are found by minimising the residual of
      a1*D2 + a2*D1 = 0   subject to   integral(a1 + a2) = 1.

    A two-stage approach is used:
      Step 1:  Try constant alpha.
      Step 1b: If Step 1 succeeds, solve for linear beta with b1 = b2.
      Step 2:  If Step 1 fails, solve for general linear alpha via KKT.
      Step 2b: Solve for general linear beta with gauge constraint
               integral(a1*b1 - a2*b2) = 0.

    The beta functions (b1, b2) satisfy   b1*D2 - b2*D1 = D3.

    ## Return convention
    The per-interface function returns 8 values:
      [a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1]
    where 1/2 = first/second patch side, _0/_1 = value at t=0 / t=1.
    The parameter t runs along the tangential direction of each patch
    in its original orientation.

    ## Matrix layout
    The matrix function collects all interfaces into an (nPatches x 16)
    matrix, with 4 columns per side of each patch:
      side index: west=0, east=1, south=2, north=3
      columns 4*s+0 = alpha_0,  4*s+1 = alpha_1
      columns 4*s+2 = beta_0,   4*s+3 = beta_1
    Initialised with alpha=1, beta=0 (identity for boundary sides).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#include <iostream>
#include <set>
#include <algorithm>
#include <gismo.h>

using namespace gismo;

// ====================================================================
// Small helpers
// ====================================================================

/// 2D cross-determinant  det([a | b])
template <class T>
inline T det2(const gsVector<T>& a, const gsVector<T>& b)
{ return a(0)*b(1) - a(1)*b(0); }

/// Extract partial-derivative vector from deriv_into output.
///   deriv_into layout (2D geometry, parDim=2):
///     row 0 = d_0 f_x,  row 1 = d_1 f_x,
///     row 2 = d_0 f_y,  row 3 = d_1 f_y.
/// Returns the 2D vector (d_dir f_x, d_dir f_y) at column col.
template <class T>
inline gsVector<T> getPartial(const gsMatrix<T>& d, index_t dir, index_t col)
{
    gsVector<T> v(2);
    v(0) = d(dir,     col);
    v(1) = d(dir + 2, col);
    return v;
}

/// Collect knot-span break points for the 1D component basis in
/// direction dir.  Falls back to support-interval endpoints when
/// the component basis is not a gsBSplineBasis.
template <class T>
std::vector<T> collectBreaks(const gsGeometry<T>& geo, short_t dir)
{
    const gsBSplineBasis<T>* bb =
        dynamic_cast<const gsBSplineBasis<T>*>(&geo.basis().component(dir));
    if (bb)
        return bb->knots().breaks();

    gsMatrix<T> sup = geo.support();
    std::vector<T> brk;
    brk.push_back(sup(dir, 0));
    brk.push_back(sup(dir, 1));
    return brk;
}

// ====================================================================
// Solver helpers (replace goto with functions that return bool)
// ====================================================================

/// Step 1 + 1b: Try constant alpha, then equal-beta.
/// On success writes a10..b21 and returns true.
template <class T>
bool solveConstantAlpha(
    const gsVector<T>& D1, const gsVector<T>& D2, const gsVector<T>& D3,
    const gsVector<T>& w,  const gsVector<T>& t_vals,
    T& a10, T& a11, T& a20, T& a21,
    T& b10, T& b11, T& b20, T& b21,
    T eps)
{
    const index_t N = D1.size();
    auto integrate = [&](const gsVector<T>& f) -> T { return w.dot(f); };

    T intD1D1 = integrate((D1.array()*D1.array()).matrix());
    T intD1D2 = integrate((D1.array()*D2.array()).matrix());
    T intD2D2 = integrate((D2.array()*D2.array()).matrix());

    T denom = intD1D1 - 2*intD1D2 + intD2D2;
    T a1c, a2c, err;
    if (std::abs(denom) < 1e-30)
    { a1c = 0.5;  a2c = 0.5;  err = intD1D1; }
    else
    {
        a1c = (intD1D1 - intD1D2) / denom;
        a2c = 1.0 - a1c;
        err = a1c*a1c*intD2D2 + 2*a1c*a2c*intD1D2 + a2c*a2c*intD1D1;
    }
    gsInfo << "  alpha_1=" << a1c << "  alpha_2=" << a2c
           << "  err=" << err << "\n";

    if (err >= eps)
        return false;

    a10 = a1c;  a11 = a1c;
    a20 = a2c;  a21 = a2c;

    // Step 1b: linear beta with b1 = b2
    gsInfo << "--- Step 1b: linear beta (b1=b2) ---\n";
    gsVector<T> S = D2 - D1;
    gsMatrix<T> Ab(2,2);  Ab.setZero();
    gsVector<T> rb(2);    rb.setZero();
    for (index_t i = 0; i < N; ++i)
    {
        T phi[2] = {1-t_vals(i), t_vals(i)};
        for (index_t j = 0; j < 2; ++j)
        {
            for (index_t k = 0; k < 2; ++k)
                Ab(j,k) += w(i) * phi[j]*phi[k] * S(i)*S(i);
            rb(j) += w(i) * phi[j] * S(i) * D3(i);
        }
    }
    gsVector<T> bc = Ab.fullPivLu().solve(rb);
    b10 = bc(0);  b11 = bc(1);
    b20 = b10;    b21 = b11;

    gsVector<T> bfit(N);
    for (index_t i = 0; i < N; ++i)
    {
        T bv = b10*(1-t_vals(i)) + b11*t_vals(i);
        T ri = bv*S(i) - D3(i);
        bfit(i) = ri*ri;
    }
    T berr = integrate(bfit);
    gsInfo << "  beta=" << b10 << "*(1-t)+" << b11 << "*t"
           << "  err=" << berr << "\n";

    return (berr < eps);
}

/// Step 2 + 2b: General linear alpha via KKT, then general linear beta.
/// On success writes a10..b21 and returns true.
template <class T>
bool solveLinearAlpha(
    const gsVector<T>& D1, const gsVector<T>& D2, const gsVector<T>& D3,
    const gsVector<T>& w,  const gsVector<T>& t_vals,
    T& a10, T& a11, T& a20, T& a21,
    T& b10, T& b11, T& b20, T& b21,
    T eps)
{
    const index_t N = D1.size();
    auto integrate = [&](const gsVector<T>& f) -> T { return w.dot(f); };

    // Basis vectors: f = [phi0*D2, phi1*D2, phi0*D1, phi1*D1]
    gsMatrix<T> fv(N, 4);
    for (index_t i = 0; i < N; ++i)
    {
        T p0 = 1-t_vals(i), p1 = t_vals(i);
        fv(i,0) = p0*D2(i);
        fv(i,1) = p1*D2(i);
        fv(i,2) = p0*D1(i);
        fv(i,3) = p1*D1(i);
    }

    // Gram matrix
    gsMatrix<T> G(4,4);  G.setZero();
    for (index_t j = 0; j < 4; ++j)
        for (index_t k = j; k < 4; ++k)
        {
            gsVector<T> pr(N);
            for (index_t i = 0; i < N; ++i)
                pr(i) = fv(i,j)*fv(i,k);
            G(j,k) = integrate(pr);
            G(k,j) = G(j,k);
        }

    // Constraint: integral(a1+a2)=1  =>  c^T x = 1
    gsVector<T> c(4);  c.setConstant(0.5);

    // KKT system
    gsMatrix<T> KKT(5,5);  KKT.setZero();
    KKT.block(0,0,4,4) = 2.0*G;
    KKT.block(0,4,4,1) = c;
    KKT.block(4,0,1,4) = c.transpose();
    gsVector<T> rhs(5);  rhs.setZero();  rhs(4) = 1.0;

    gsVector<T> sol = KKT.fullPivLu().solve(rhs);
    a10 = sol(0);  a11 = sol(1);  a20 = sol(2);  a21 = sol(3);

    // Alpha error
    gsVector<T> aerr_v(N);
    for (index_t i = 0; i < N; ++i)
    {
        T a1v = a10*(1-t_vals(i)) + a11*t_vals(i);
        T a2v = a20*(1-t_vals(i)) + a21*t_vals(i);
        T r   = a1v*D2(i) + a2v*D1(i);
        aerr_v(i) = r*r;
    }
    T aerr = integrate(aerr_v);
    gsInfo << "  alpha_1 = " << a10 << "*(1-t) + " << a11 << "*t\n";
    gsInfo << "  alpha_2 = " << a20 << "*(1-t) + " << a21 << "*t\n";
    gsInfo << "  err = " << aerr << "\n";

    if (aerr >= eps)
        return false;

    // Step 2b: general linear beta
    gsInfo << "--- Step 2b: linear beta ---\n";

    // Reuse fv with sign flips for beta: g = [phi0*D2, phi1*D2, -phi0*D1, -phi1*D1]
    gsMatrix<T> H(4,4);  H.setZero();
    gsVector<T> d_vec(4);  d_vec.setZero();
    for (index_t j = 0; j < 4; ++j)
    {
        T sj = (j < 2) ? 1.0 : -1.0;
        for (index_t k = j; k < 4; ++k)
        {
            T sk = (k < 2) ? 1.0 : -1.0;
            gsVector<T> pr(N);
            for (index_t i = 0; i < N; ++i)
                pr(i) = sj*fv(i,j) * sk*fv(i,k);
            H(j,k) = integrate(pr);
            H(k,j) = H(j,k);
        }
        gsVector<T> pr(N);
        for (index_t i = 0; i < N; ++i)
            pr(i) = sj*fv(i,j) * D3(i);
        d_vec(j) = integrate(pr);
    }

    // Gauge constraint: integral(a1*b1 - a2*b2) = 0
    gsVector<T> e(4);
    {
        gsVector<T> tmp(N);
        for (index_t i = 0; i < N; ++i)
            tmp(i) = (a10*(1-t_vals(i))+a11*t_vals(i)) * (1-t_vals(i));
        e(0) = integrate(tmp);
        for (index_t i = 0; i < N; ++i)
            tmp(i) = (a10*(1-t_vals(i))+a11*t_vals(i)) * t_vals(i);
        e(1) = integrate(tmp);
        for (index_t i = 0; i < N; ++i)
            tmp(i) = -(a20*(1-t_vals(i))+a21*t_vals(i)) * (1-t_vals(i));
        e(2) = integrate(tmp);
        for (index_t i = 0; i < N; ++i)
            tmp(i) = -(a20*(1-t_vals(i))+a21*t_vals(i)) * t_vals(i);
        e(3) = integrate(tmp);
    }

    // KKT: [2H e; e^T 0] [y;mu] = [2d;0]
    gsMatrix<T> KKT2(5,5);  KKT2.setZero();
    KKT2.block(0,0,4,4) = 2.0*H;
    KKT2.block(0,4,4,1) = e;
    KKT2.block(4,0,1,4) = e.transpose();
    gsVector<T> rhs2(5);
    rhs2.head(4) = 2.0*d_vec;
    rhs2(4) = 0.0;

    gsVector<T> sol2 = KKT2.fullPivLu().solve(rhs2);
    b10 = sol2(0);  b11 = sol2(1);
    b20 = sol2(2);  b21 = sol2(3);

    // Beta error
    gsVector<T> bfit(N);
    for (index_t i = 0; i < N; ++i)
    {
        T p0 = 1-t_vals(i), p1 = t_vals(i);
        T b1v = b10*p0 + b11*p1;
        T b2v = b20*p0 + b21*p1;
        T r   = b1v*D2(i) - b2v*D1(i) - D3(i);
        bfit(i) = r*r;
    }
    T berr = integrate(bfit);
    gsInfo << "  beta_1 = " << b10 << "*(1-t) + " << b11 << "*t\n";
    gsInfo << "  beta_2 = " << b20 << "*(1-t) + " << b21 << "*t\n";
    gsInfo << "  err = " << berr << "\n";

    return (berr < eps);
}

// ====================================================================
// Per-interface computation (generic -- works with any gsGeometry)
// ====================================================================

/// Compute gluing data for a single interface.
///
/// Works with any gsGeometry<T> (not just gsTensorBSpline).
/// No reparametrisation is performed; instead, sign factors
/// s_k = (-1)^{par_k} correct for the orientation.
///
/// Returns a vector of 8 values:
///   [a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1]
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

    // Normal direction = direction perpendicular to the interface side
    // Tangential direction = the other one
    const short_t normDir1 = ps1.direction();
    const short_t tanDir1  = 1 - normDir1;
    const short_t normDir2 = ps2.direction();
    const short_t tanDir2  = 1 - normDir2;

    // par = false(0) means interface at start of normal dir,
    //        true(1) means interface at end of normal dir.
    const bool par1 = ps1.parameter();
    const bool par2 = ps2.parameter();

    // Sign factors: after reparametrisation to west, the Jacobian
    // determinant picks up factor (-1)^par.  We apply this sign
    // directly to D1, D2, D3 instead of reparametrising.
    const T signD1 = par1 ? -1.0 : 1.0;
    const T signD2 = par2 ? -1.0 : 1.0;
    const T signD3 = signD1 * signD2;

    // Support (parameter domain) of each patch
    gsMatrix<T> sup1 = geo1.support();
    gsMatrix<T> sup2 = geo2.support();

    // Interface coordinate: value of normal coordinate at the interface
    const T ifcCoord1 = sup1(normDir1, par1 ? 1 : 0);
    const T ifcCoord2 = sup2(normDir2, par2 ? 1 : 0);

    // Tangential flip: does the tangential direction of patch 2 run
    // opposite to that of patch 1?
    const bool flipped = !interf.dirOrientation(ps1, tanDir1);

    gsInfo << "Interface: " << interf << "\n"
           << "  Patch " << ps1.patch << " side " << ps1
           << " (normDir=" << normDir1 << " par=" << par1
           << " signD=" << signD1 << ")\n"
           << "  Patch " << ps2.patch << " side " << ps2
           << " (normDir=" << normDir2 << " par=" << par2
           << " signD=" << signD2 << ")\n"
           << "  Tangential flipped: " << flipped << "\n";

    // ---- Gauss quadrature along tangential direction of patch 1 ----
    const short_t deg1 = geo1.basis().degree(tanDir1);
    const short_t deg2 = geo2.basis().degree(tanDir2);
    const index_t nGauss = numGaussPerSpan > 0
                               ? numGaussPerSpan
                               : 2*std::max(deg1, deg2) + 1;

    std::vector<T> breaks1 = collectBreaks(geo1, tanDir1);
    std::vector<T> breaks2 = collectBreaks(geo2, tanDir2);

    // Map patch2 breaks into patch1 tangential range
    const T t1a = sup1(tanDir1, 0), t1b = sup1(tanDir1, 1);
    const T t2a = sup2(tanDir2, 0), t2b = sup2(tanDir2, 1);
    for (T& br : breaks2)
    {
        T s = (br - t2a) / (t2b - t2a);   // normalise to [0,1]
        if (flipped) s = 1.0 - s;
        br = t1a + s * (t1b - t1a);        // map to patch1 range
    }

    std::set<T> breakSet(breaks1.begin(), breaks1.end());
    breakSet.insert(breaks2.begin(), breaks2.end());
    std::vector<T> mergedBreaks(breakSet.begin(), breakSet.end());

    gsGaussRule<T> gaussRule(nGauss);
    gsMatrix<T> gaussNodes;
    gsVector<T> w;
    gaussRule.mapToAll(mergedBreaks, gaussNodes, w);
    const index_t N = gaussNodes.cols();

    gsInfo << "Gauss quadrature: " << mergedBreaks.size()-1
           << " spans x " << nGauss << " pts = " << N << " total\n";

    // ---- Evaluate derivatives at quadrature points ----
    gsMatrix<T> pts1(2, N), pts2(2, N);
    for (index_t i = 0; i < N; ++i)
    {
        T t = gaussNodes(0, i);   // parameter along patch1 tangential

        // Patch 1: fix normal coord at interface, vary tangential
        pts1(normDir1, i) = ifcCoord1;
        pts1(tanDir1,  i) = t;

        // Patch 2: map t into patch2's tangential range
        T s = (t - t1a) / (t1b - t1a);       // normalise to [0,1]
        if (flipped) s = 1.0 - s;
        T t2 = t2a + s * (t2b - t2a);

        pts2(normDir2, i) = ifcCoord2;
        pts2(tanDir2,  i) = t2;
    }

    gsMatrix<T> derivs1, derivs2;
    geo1.deriv_into(pts1, derivs1);
    geo2.deriv_into(pts2, derivs2);

    // ---- Compute D1, D2, D3 with sign corrections ----
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

        // Normalise t to [0,1] for the linear basis {1-t, t}
        t_vals(i) = (gaussNodes(0,i) - t1a) / (t1b - t1a);
    }

    // Sanity check
    if (D1.minCoeff()*D1.maxCoeff() < 0 || D2.minCoeff()*D2.maxCoeff() < 0)
    {
        gsInfo << "WARNING: D1 or D2 changes sign.\n";
        return result;
    }

    gsInfo << "D1 range: [" << D1.minCoeff() << ", " << D1.maxCoeff() << "]\n"
           << "D2 range: [" << D2.minCoeff() << ", " << D2.maxCoeff() << "]\n"
           << "D3 range: [" << D3.minCoeff() << ", " << D3.maxCoeff() << "]\n";

    // ---- Solve for alpha and beta ----
    T a10, a11, a20, a21, b10, b11, b20, b21;

    gsInfo << "\n--- Step 1: constant alpha ---\n";
    if (solveConstantAlpha(D1, D2, D3, w, t_vals,
                           a10, a11, a20, a21,
                           b10, b11, b20, b21, eps))
    {
        gsInfo << "  SUCCESS (Step 1).\n";
        success = true;
    }
    else
    {
        gsInfo << "  -> Step 2: linear alpha\n";
        if (solveLinearAlpha(D1, D2, D3, w, t_vals,
                             a10, a11, a20, a21,
                             b10, b11, b20, b21, eps))
        {
            gsInfo << "  SUCCESS (Step 2).\n";
            success = true;
        }
        else
            gsInfo << "  FAILED: not AS-G1.\n";
    }

    if (!success)
        return result;

    gsInfo << "\n=== Gluing data (patch1 tangential parametrisation) ===\n"
           << "  alpha_1(t) = " << a10 << "*(1-t) + " << a11 << "*t\n"
           << "  alpha_2(t) = " << a20 << "*(1-t) + " << a21 << "*t\n"
           << "  beta_1(t)  = " << b10 << "*(1-t) + " << b11 << "*t\n"
           << "  beta_2(t)  = " << b20 << "*(1-t) + " << b21 << "*t\n";

    // ---- Pack result, handling tangential flip for patch 2 ----
    result(0) = a10;  result(1) = a11;
    result(2) = b10;  result(3) = b11;
    if (flipped)
    {
        result(4) = a21;  result(5) = a20;
        result(6) = b21;  result(7) = b20;
    }
    else
    {
        result(4) = a20;  result(5) = a21;
        result(6) = b20;  result(7) = b21;
    }
    return result;
}

// ====================================================================
// Collect all interfaces into a matrix
// ====================================================================

template <class T>
gsMatrix<T> computeGluingDataMatrix(
    const gsMultiPatch<T>& mp,
    const T eps = 1e-8,
    index_t numGaussPerSpan = 0)
{
    const index_t nP = mp.nPatches();
    gsMatrix<T> mat(nP, 16);

    // Default: alpha=1, beta=0 for every side
    for (index_t p = 0; p < nP; ++p)
        for (index_t s = 0; s < 4; ++s)
        {
            mat(p, 4*s+0) = 1.0;
            mat(p, 4*s+1) = 1.0;
            mat(p, 4*s+2) = 0.0;
            mat(p, 4*s+3) = 0.0;
        }

    auto sideToCol = [](const patchSide& ps) -> index_t
    { return static_cast<index_t>(ps.side()) - 1; };

    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const boundaryInterface& interf = *it;
        bool ok = false;
        gsVector<T> gd = computeGluingDataForInterface(
            mp, interf, ok, eps, numGaussPerSpan);
        if (!ok) { gsInfo << "Failed for interface " << interf << "\n"; continue; }

        const patchSide& ps1 = interf.first();
        const patchSide& ps2 = interf.second();
        const index_t col1 = 4 * sideToCol(ps1);
        const index_t col2 = 4 * sideToCol(ps2);

        mat(ps1.patch, col1+0) = gd(0);
        mat(ps1.patch, col1+1) = gd(1);
        mat(ps1.patch, col1+2) = gd(2);
        mat(ps1.patch, col1+3) = gd(3);

        mat(ps2.patch, col2+0) = gd(4);
        mat(ps2.patch, col2+1) = gd(5);
        mat(ps2.patch, col2+2) = gd(6);
        mat(ps2.patch, col2+3) = gd(7);
    }
    return mat;
}

// ====================================================================
// main
// ====================================================================

int main(int argc, char* argv[])
{
    using T = double;

    std::string geometry("domain2d/two_bicubic_patches.xml");
    index_t numGaussPerSpan = 0;

    gsCmdLine cmd("Compute gluing data (generic, no reparametrisation).");
    cmd.addString("f", "file", "G+Smo input multi-patch file.", geometry);
    cmd.addInt("n", "numGauss",
               "Number of Gauss points per knot span (0=auto).", numGaussPerSpan);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
    if (!mpPtr) { gsInfo << "Cannot read " << geometry << ".\n"; return -1; }
    gsMultiPatch<>& mp = *mpPtr;
    mp.computeTopology();

    gsInfo << "Patches: " << mp.nPatches()
           << "  Interfaces: " << mp.nInterfaces() << "\n\n";

    gsMatrix<T> mat = computeGluingDataMatrix(mp);

    const char* sideNames[4] = {"west", "east", "south", "north"};
    gsInfo << "=== Gluing data matrix (" << mat.rows()
           << " x " << mat.cols() << ") ===\n\n";
    for (index_t p = 0; p < mat.rows(); ++p)
    {
        gsInfo << "Patch " << p << ":\n";
        for (index_t s = 0; s < 4; ++s)
        {
            T a0 = mat(p, 4*s+0), a1 = mat(p, 4*s+1);
            T b0 = mat(p, 4*s+2), b1 = mat(p, 4*s+3);
            gsInfo << "  " << std::setw(5) << sideNames[s]
                   << ":  alpha=" << a0 << "*(1-t)+" << a1 << "*t"
                   << "   beta=" << b0 << "*(1-t)+" << b1 << "*t\n";
        }
        gsInfo << "\n";
    }
    gsInfo << "Raw matrix:\n" << mat << "\n";

    return 0;
}
