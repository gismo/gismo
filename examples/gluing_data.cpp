/** @file gluing_data.cpp

    @brief Compute AS-G1 gluing data for multi-patch geometries.

    ## Overview
    This file computes the gluing data (alpha and beta functions) needed
    for AS (Approximate Smooth) G1 constructions on multi-patch spline
    geometries.  In a G1 construction, two patches sharing an interface
    must satisfy compatibility conditions that involve linear functions
    alpha and beta along each interface edge.

    ## Algorithm
    For each interface, both adjacent patches are reparametrized so that
    the shared edge sits at u=0 (west side).  After reparametrization,
    the normal direction is u (index 0) and the tangential direction is
    v (index 1) for both patches.  Three determinant quantities are
    evaluated at Gauss quadrature points along the interface:

      D1(t) = det( d_u G1, d_v G1 )    -- Jacobian determinant of patch 1
      D2(t) = det( d_u G2, d_v G2 )    -- Jacobian determinant of patch 2
      D3(t) = det( d_u G1, d_u G2 )    -- cross of normal derivatives

    The alpha functions (a1, a2) are found by minimising the residual of
    the alpha condition  a1*D2 + a2*D1 = 0  under the normalisation
    constraint  integral(a1 + a2) = 1.  A two-stage approach is used:

      Step 1:  Try a constant alpha  a1 = a2 = const.
      Step 1b: If Step 1 succeeds, solve for linear beta with b1 = b2.
      Step 2:  If Step 1 fails, solve for general linear alpha via KKT.
      Step 2b: Solve for general linear beta with a gauge constraint
               integral(a1*b1 - a2*b2) = 0.

    The beta functions (b1, b2) satisfy the beta condition
      b1*D2 - b2*D1 = D3.

    ## Return convention
    The per-interface function returns 8 values:
      [a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1]
    where 1/2 = first/second patch side, _0/_1 = value at t=0 / t=1.
    The parameter t runs along the tangential direction in the
    **original** patch orientation (before reparametrization).

    ## Matrix layout
    The matrix function collects all interfaces into an (nPatches x 16)
    matrix, with 4 columns per side of each patch:
      side index: west=0, east=1, south=2, north=3
      columns 4*s+0 = alpha_0,  4*s+1 = alpha_1
      columns 4*s+2 = beta_0,   4*s+3 = beta_1
    Initialised with alpha=1, beta=0 (identity for boundary sides with
    no adjacent patch).

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
//  reparametriseToWest
// ====================================================================
//
// Purpose:
//   Reparametrize a 2D tensor B-spline patch so that the given
//   interface side becomes the west side (u = 0).  After the call,
//   the normal direction to the interface is u (index 0) and the
//   tangential direction along the interface is v (index 1).
//
// Why we do this:
//   By moving the interface to a canonical position (west = u=0),
//   the rest of the algorithm can uniformly treat every interface
//   the same way: the normal derivative is always d/du and the
//   tangential derivative is always d/dv.
//
// How the tangential parameter maps back to the original patch:
//   Case            | Operation       | Tangential v = original ...
//   ---------------------------------------------------------------
//   west   (u=0)    | nothing         | v  (same direction)
//   east   (u=1)    | reverse u       | v  (same direction)
//   south  (v=0)    | swap u <-> v    | v  (= original u, same dir)
//   north  (v=1)    | swap + reverse u| v  (= original u, same dir)
//
//   In every case the tangential direction is orientation-preserving
//   with respect to the original parametric direction along the
//   interface.  This is important: the returned coefficients _0/_1
//   will correspond to the correct ends of the original interface.
//
// Parameters:
//   geo  [in/out]  The tensor B-spline patch to reparametrize.
//   ps   [in]      The patchSide describing which side is the
//                   interface (direction 0/1 and parameter false/true).
// ====================================================================
template <class T>
void reparametriseToWest(gsTensorBSpline<2,T> & geo, const patchSide & ps)
{
    const short_t dir = ps.direction();
    const bool    par = ps.parameter();

    // ---------------------------------------------------------------
    // Helper lambda: swap the two parametric directions u <-> v.
    //
    // This is equivalent to transposing the control point grid.
    // If the original grid has nu control points in u and nv in v,
    // the new grid has nv points in the first direction and nu in
    // the second.  The knot vectors are also swapped.
    //
    // Control point indexing convention:
    //   original: coef index = i + j * nu   (i=u index, j=v index)
    //   after:    coef index = j + i * nv   (j becomes new "u", i new "v")
    // ---------------------------------------------------------------
    auto swapUV = [](gsTensorBSpline<2,T> & g)
    {
        const index_t nu = g.basis().size(0);
        const index_t nv = g.basis().size(1);
        const index_t d  = g.geoDim();
        gsMatrix<T> oldC = g.coefs();
        gsMatrix<T> newC(nv*nu, d);
        for (index_t j = 0; j < nv; ++j)
            for (index_t i = 0; i < nu; ++i)
                newC.row(j + i*nv) = oldC.row(i + j*nu);

        gsKnotVector<T> kvU = g.knots(0);
        gsKnotVector<T> kvV = g.knots(1);
        g = gsTensorBSpline<2,T>(kvV, kvU, newC);
    };

    // ---------------------------------------------------------------
    // Helper lambda: reverse parametric direction d (0 = u, 1 = v).
    //
    // Reverses the ordering of control points along direction d and
    // reverses the corresponding knot vector.  This maps  t -> 1-t
    // for the affected direction, flipping the orientation.
    //
    // For example, reversing u turns an east interface (u=1) into a
    // west interface (u=0) while keeping v in the same direction.
    // ---------------------------------------------------------------
    auto reverseDir = [](gsTensorBSpline<2,T> & g, short_t d)
    {
        const index_t nu = g.basis().size(0);
        const index_t nv = g.basis().size(1);
        const index_t gd = g.geoDim();
        gsMatrix<T> oldC = g.coefs();
        gsMatrix<T> newC(nu*nv, gd);
        for (index_t j = 0; j < nv; ++j)
            for (index_t i = 0; i < nu; ++i)
            {
                index_t iNew = (d == 0) ? (nu-1-i) : i;
                index_t jNew = (d == 1) ? (nv-1-j) : j;
                newC.row(iNew + jNew*nu) = oldC.row(i + j*nu);
            }

        gsKnotVector<T> kvU = g.knots(0);
        gsKnotVector<T> kvV = g.knots(1);
        if (d == 0) kvU.reverse();
        else        kvV.reverse();
        g = gsTensorBSpline<2,T>(kvU, kvV, newC);
    };

    // Dispatch based on which side is the interface:
    //   dir=0, par=false  ->  west  (u=0):  already at west, do nothing
    //   dir=0, par=true   ->  east  (u=1):  reverse u  =>  u=1 maps to u=0
    //   dir=1, par=false  ->  south (v=0):  swap u<->v =>  old v=0 is now u=0
    //   dir=1, par=true   ->  north (v=1):  swap + reverse u  =>  old v=1 is now u=0
    if (dir == 0 && !par)       {}                              // west — done
    else if (dir == 0 && par)   { reverseDir(geo, 0); }        // east
    else if (dir == 1 && !par)  { swapUV(geo); }               // south
    else                        { swapUV(geo); reverseDir(geo, 0); } // north
}

// ====================================================================
//  computeGluingDataForInterface
// ====================================================================
//
// Purpose:
//   Compute the alpha and beta gluing data for a single interface
//   between two patches.  Returns a vector of 8 coefficients that
//   fully define the linear alpha and beta functions for both patches
//   on this interface.
//
// Algorithm outline:
//   1. Copy the two patches and reparametrize both to west.
//   2. Check that the boundary curves match (detect orientation flip).
//   3. Set up composite Gauss quadrature along the interface.
//   4. Evaluate the Jacobian-derived quantities D1, D2, D3.
//   5. Solve for alpha and beta in two stages (constant -> linear).
//   6. If the tangential direction was flipped for patch 2, swap
//      the _0 and _1 coefficients of patch 2 so that they refer
//      to the original (pre-reparametrization) orientation.
//
// Parameters:
//   mp       [in]   The multipatch geometry.
//   interf   [in]   The boundaryInterface to process.
//   success  [out]  Whether the solve succeeded within tolerance.
//   eps      [in]   Tolerance for residual checks (default 1e-8).
//   numGaussPerSpan [in]  Gauss points per knot span (0 = auto).
//
// Returns:
//   gsVector of 8 values on success, or zeros on failure:
//     [a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1]
//   where 1/2 = first/second patch of the interface, and
//   _0/_1 = value at tangential parameter t=0 / t=1 in original
//   patch orientation.
// ====================================================================
template <class T>
gsVector<T> computeGluingDataForInterface(
    const gsMultiPatch<T> & mp,
    const boundaryInterface & interf,
    bool & success,
    const T eps = 1e-8,
    index_t numGaussPerSpan = 0)
{
    success = false;
    gsVector<T> result(8);
    result.setZero();

    const patchSide firstPS  = interf.first();
    const patchSide secondPS = interf.second();

    // --- Step 0: Cast and copy the two patches ---
    // We need gsTensorBSpline<2,T> specifically because we manipulate
    // control points and knot vectors directly during reparametrization.
    // We make copies so the original multipatch is not modified.
    const auto* p1 = dynamic_cast<const gsTensorBSpline<2,T>*>(
                         &mp.patch(firstPS.patch));
    const auto* p2 = dynamic_cast<const gsTensorBSpline<2,T>*>(
                         &mp.patch(secondPS.patch));
    if (!p1 || !p2)
    {
        gsInfo << "Patches must be gsTensorBSpline<2,T>.\n";
        return result;
    }
    gsTensorBSpline<2,T> geo1 = *p1;
    gsTensorBSpline<2,T> geo2 = *p2;

    // --- Step 1a: Reparametrize both patches so interface is at u=0 (west) ---
    // After this, for both geo1 and geo2:
    //   - The interface sits at u = 0  (the west boundary)
    //   - The normal direction to the interface is u  (index 0)
    //   - The tangential direction along the interface is v  (index 1)
    reparametriseToWest(geo1, firstPS);
    reparametriseToWest(geo2, secondPS);

    // After reparametrization: normal = 0 (u), tangential = 1 (v)
    const index_t normDir = 0;  // u-direction = normal to interface
    const index_t tanDir  = 1;  // v-direction = tangential along interface

    // support(dir, 0) = lower bound, support(dir, 1) = upper bound
    gsMatrix<T> support1 = geo1.support();
    gsMatrix<T> support2 = geo2.support();

    // --- Step 1b: Check boundary matching and detect orientation flip ---
    // After reparametrization, the boundary curves of both patches at
    // u=0 should be the same geometric curve.  However, they might be
    // traversed in opposite directions (the tangential parameter runs
    // in opposite directions on the two patches).
    //
    // We sample the boundaries at Nq points and compare:
    //   Case 1: qVals1 == qVals2         ->  same direction (flipped=false)
    //   Case 2: qVals1 == reverse(qVals2) ->  opposite direction (flipped=true)
    //           In this case, we reverse the v-direction of geo2 so that
    //           the boundary curves match directly.  Later, we swap the
    //           _0 and _1 coefficients of patch 2 to undo this flip.
    bool flipped = false;
    {
        const index_t Nq = 5;  // number of test points for matching check
        gsMatrix<T> qPts1(2, Nq), qPts2(2, Nq);
        for (index_t i = 0; i < Nq; ++i)
        {
            T tc = static_cast<T>(i) / (Nq - 1);
            qPts1(normDir, i) = support1(normDir, 0);
            qPts1(tanDir,  i) = support1(tanDir, 0)*(1-tc)
                              + support1(tanDir, 1)*tc;
            qPts2(normDir, i) = support2(normDir, 0);
            qPts2(tanDir,  i) = support2(tanDir, 0)*(1-tc)
                              + support2(tanDir, 1)*tc;
        }
        gsMatrix<T> qVals1, qVals2;
        geo1.eval_into(qPts1, qVals1);
        geo2.eval_into(qPts2, qVals2);

        const T tol = 1e-10;
        if ((qVals1 - qVals2).norm() < tol)
        {
            // Boundary curves match directly — same tangential direction
            flipped = false;
        }
        else if ((qVals1 - qVals2.rowwise().reverse()).norm() < tol)
        {
            // Boundary curves match only when one is reversed — flip detected.
            // We reverse the v-direction of geo2 so that computations below
            // can proceed with matching boundaries.  The flip will be undone
            // at the end by swapping _0 and _1 for patch 2's coefficients.
            flipped = true;
            // Reverse the v-direction of geo2 to make boundaries match
            const index_t nu = geo2.basis().size(0);
            const index_t nv = geo2.basis().size(1);
            gsMatrix<T> oldC = geo2.coefs();
            gsMatrix<T> newC(nu*nv, geo2.geoDim());
            for (index_t j = 0; j < nv; ++j)
                for (index_t i = 0; i < nu; ++i)
                    newC.row(i + (nv-1-j)*nu) = oldC.row(i + j*nu);
            gsKnotVector<T> kvU = geo2.knots(0);
            gsKnotVector<T> kvV = geo2.knots(1);
            kvV.reverse();
            geo2 = gsTensorBSpline<2,T>(kvU, kvV, newC);
            support2 = geo2.support();
        }
        else
        {
            gsInfo << "Boundary values do NOT match for interface "
                   << interf << "!\n";
            return result;
        }
    }

    // --- Step 2: Set up Gauss quadrature along the interface (v-direction) ---
    //
    // We use composite Gauss quadrature with knot spans as sub-intervals.
    // The knot break points from both patches are merged so that every
    // sub-interval is a single span in both bases.  The default number
    // of Gauss points per sub-interval is  2*max(deg1, deg2) + 1,
    // which is enough to integrate products of degree-p basis functions
    // exactly.
    const short_t deg1 = geo1.basis().degree(tanDir);
    const short_t deg2 = geo2.basis().degree(tanDir);
    const index_t nGauss = numGaussPerSpan > 0
                               ? numGaussPerSpan
                               : 2*std::max(deg1, deg2) + 1;

    // Collect all knot break points from both patches along the tangential
    // direction.  Using a set automatically removes duplicates and sorts.
    std::set<T> breakSet;
    {
        const auto* bb1 = dynamic_cast<const gsBSplineBasis<T>*>(
                              &geo1.basis().component(tanDir));
        if (bb1) { auto b = bb1->knots().breaks(); breakSet.insert(b.begin(), b.end()); }
        else { breakSet.insert(support1(tanDir,0)); breakSet.insert(support1(tanDir,1)); }

        const auto* bb2 = dynamic_cast<const gsBSplineBasis<T>*>(
                              &geo2.basis().component(tanDir));
        if (bb2) { auto b = bb2->knots().breaks(); breakSet.insert(b.begin(), b.end()); }
        else { breakSet.insert(support2(tanDir,0)); breakSet.insert(support2(tanDir,1)); }
    }
    std::vector<T> mergedBreaks(breakSet.begin(), breakSet.end());

    // Map Gauss nodes and weights onto all sub-intervals
    gsGaussRule<T> gaussRule(nGauss);
    gsMatrix<T> gaussNodes;  // 1 x N matrix of quadrature abscissae
    gsVector<T> w;            // N-vector of quadrature weights
    gaussRule.mapToAll(mergedBreaks, gaussNodes, w);
    const index_t N = gaussNodes.cols();  // total number of quadrature points

    // --- Step 3: Evaluate geometry derivatives at quadrature points ---
    //
    // We evaluate derivatives at points on the interface (u = 0) for
    // both patches.  The evaluation points are (u=0, v=gaussNode[i]).
    gsMatrix<T> pts1(2, N), pts2(2, N);
    for (index_t i = 0; i < N; ++i)
    {
        pts1(normDir, i) = support1(normDir, 0);  // u = 0 (interface)
        pts1(tanDir,  i) = gaussNodes(0, i);      // v = quadrature point
        pts2(normDir, i) = support2(normDir, 0);  // u = 0 (interface)
        pts2(tanDir,  i) = gaussNodes(0, i);      // v = quadrature point
    }
    gsMatrix<T> derivs1, derivs2;
    geo1.deriv_into(pts1, derivs1);  // 4 x N matrix
    geo2.deriv_into(pts2, derivs2);  // 4 x N matrix

    // deriv_into layout for 2D geometry (geoDim=2, parDim=2):
    //   row 0 = d_u f_x   (partial of x w.r.t. u)
    //   row 1 = d_v f_x   (partial of x w.r.t. v)
    //   row 2 = d_u f_y   (partial of y w.r.t. u)
    //   row 3 = d_v f_y   (partial of y w.r.t. v)
    //
    // getPartial(derivs, dir, i) extracts the 2D vector
    //   (d_dir f_x,  d_dir f_y)  at quadrature point i.
    auto getPartial = [](const gsMatrix<T>& d, index_t dir, index_t i) -> gsVector<T>
    {
        gsVector<T> v(2);
        v(0) = d(dir,     i);   // d_dir f_x
        v(1) = d(dir + 2, i);   // d_dir f_y
        return v;
    };

    // 2x2 determinant of two 2D column vectors: det([a | b])
    auto det2 = [](const gsVector<T>& a, const gsVector<T>& b) -> T
    { return a(0)*b(1) - a(1)*b(0); };

    // --- Step 4: Compute D1, D2, D3 at every quadrature point ---
    //
    //   D1(t) = det( d_u G1, d_v G1 )   -- Jacobian determinant of patch 1
    //           This is the signed area element of patch 1 at the interface.
    //
    //   D2(t) = det( d_u G2, d_v G2 )   -- Jacobian determinant of patch 2
    //           This is the signed area element of patch 2 at the interface.
    //
    //   D3(t) = det( d_u G1, d_u G2 )   -- determinant of the two normal
    //           derivatives.  This measures how the normals of the two
    //           patches differ across the interface.
    //
    //   t_vals stores the quadrature abscissa for use as the linear
    //   interpolation parameter in the alpha/beta basis functions:
    //     f(t) = f_0 * (1 - t) + f_1 * t
    gsVector<T> D1(N), D2(N), D3(N), t_vals(N);
    for (index_t i = 0; i < N; ++i)
    {
        gsVector<T> dG1dn = getPartial(derivs1, normDir, i);  // d_u G1
        gsVector<T> dG1dt = getPartial(derivs1, tanDir,  i);  // d_v G1
        gsVector<T> dG2dn = getPartial(derivs2, normDir, i);  // d_u G2
        gsVector<T> dG2dt = getPartial(derivs2, tanDir,  i);  // d_v G2

        D1(i) = det2(dG1dn, dG1dt);   // det(d_u G1, d_v G1)
        D2(i) = det2(dG2dn, dG2dt);   // det(d_u G2, d_v G2)
        D3(i) = det2(dG1dn, dG2dn);   // det(d_u G1, d_u G2)
        t_vals(i) = gaussNodes(0, i);  // tangential parameter value
    }

    // Sanity check: D1 and D2 should not change sign along the interface.
    // A sign change would indicate a degenerate or self-intersecting patch.
    if (D1.minCoeff()*D1.maxCoeff() < 0 || D2.minCoeff()*D2.maxCoeff() < 0)
    {
        gsInfo << "WARNING: D1 or D2 changes sign for interface "
               << interf << ".\n";
        return result;
    }

    // Helper: compute the integral of a function sampled at quadrature points.
    auto integrate = [&](const gsVector<T>& f) -> T { return w.dot(f); };

    // --- Step 5: Solve for alpha and beta ---
    //
    // We seek linear functions (in the tangential parameter t):
    //   a1(t) = a1_0 * (1-t) + a1_1 * t      (alpha for patch 1)
    //   a2(t) = a2_0 * (1-t) + a2_1 * t      (alpha for patch 2)
    //   b1(t) = b1_0 * (1-t) + b1_1 * t      (beta for patch 1)
    //   b2(t) = b2_0 * (1-t) + b2_1 * t      (beta for patch 2)
    //
    // These must satisfy:
    //   Alpha condition:  a1(t)*D2(t) + a2(t)*D1(t) = 0   for all t
    //   Beta condition:   b1(t)*D2(t) - b2(t)*D1(t) = D3(t)  for all t
    //   Normalisation:    integral(a1 + a2) = 1
    //   Gauge (Step 2b):  integral(a1*b1 - a2*b2) = 0
    //
    // The solution proceeds in two stages:
    //   Stage A (Steps 1 + 1b): Try constant alpha first (cheaper).
    //   Stage B (Steps 2 + 2b): If constant alpha fails, use general linear.
    T a10, a11, a20, a21;  // alpha coefficients
    T b10, b11, b20, b21;  // beta coefficients

    // ---- Step 1: Try constant alpha  a1(t) = a1c,  a2(t) = a2c ----
    //
    // For constant alpha, the alpha condition  a1c*D2 + a2c*D1 = 0
    // becomes an overdetermined system.  We minimise:
    //   min  integral( (a1c*D2 + a2c*D1)^2 )   s.t.  a1c + a2c = 1
    //
    // Expanding the integral and substituting a2c = 1 - a1c:
    //   a1c = (intD1D1 - intD1D2) / (intD1D1 - 2*intD1D2 + intD2D2)
    {
        T intD1D1 = integrate((D1.array()*D1.array()).matrix());
        T intD1D2 = integrate((D1.array()*D2.array()).matrix());
        T intD2D2 = integrate((D2.array()*D2.array()).matrix());

        T denom = intD1D1 - 2*intD1D2 + intD2D2;
        T a1c, a2c, err;
        if (std::abs(denom) < 1e-30)
        {
            // Degenerate case: D1 ≈ D2 everywhere, so a1=a2=0.5
            a1c = 0.5; a2c = 0.5; err = intD1D1;
        }
        else
        {
            a1c = (intD1D1 - intD1D2) / denom;
            a2c = 1.0 - a1c;
            // Compute the residual  integral( (a1c*D2 + a2c*D1)^2 )
            err = a1c*a1c*intD2D2 + 2*a1c*a2c*intD1D2 + a2c*a2c*intD1D1;
        }

        if (err < eps)
        {
            // Constant alpha succeeded!  Set linear alpha to the constant.
            a10 = a1c; a11 = a1c;
            a20 = a2c; a21 = a2c;

            // Step 1b: Solve for linear beta with b1 = b2
            //
            // When alpha is constant, the beta condition simplifies.
            // Since a1c*D2 + a2c*D1 ≈ 0  =>  D1 ≈ -(a1c/a2c)*D2,
            // the beta condition  b1*D2 - b2*D1 = D3  with b1=b2=b
            // becomes:  b*(D2 - D1) = D3,  i.e.  b*S = D3  where S = D2-D1.
            //
            // We do a linear least-squares fit for b(t) = b_0*(1-t) + b_1*t:
            //   min integral( (b(t)*S(t) - D3(t))^2 )
            gsVector<T> S = D2 - D1;
            gsMatrix<T> Ab(2,2); Ab.setZero();
            gsVector<T> rb(2);   rb.setZero();
            for (index_t i = 0; i < N; ++i)
            {
                T phi[2] = {1-t_vals(i), t_vals(i)};  // linear basis: {1-t, t}
                for (index_t j = 0; j < 2; ++j)
                {
                    for (index_t k = 0; k < 2; ++k)
                        Ab(j,k) += w(i) * phi[j]*phi[k] * S(i)*S(i);
                    rb(j) += w(i) * phi[j] * S(i) * D3(i);
                }
            }
            gsVector<T> bc = Ab.fullPivLu().solve(rb);
            b10 = bc(0); b11 = bc(1);
            b20 = b10;   b21 = b11;   // b1 = b2 in this step

            // Check the beta residual
            gsVector<T> bfit(N);
            for (index_t i = 0; i < N; ++i)
            {
                T bv = b10*(1-t_vals(i)) + b11*t_vals(i);
                T r  = bv*S(i) - D3(i);
                bfit(i) = r*r;
            }
            T berr = integrate(bfit);
            if (berr < eps) { success = true; goto done; }
            // If beta residual is too large, fall through to Step 2
        }
    }

    // ---- Step 2: General linear alpha via constrained least squares ----
    //
    // Constant alpha failed (residual too large).  Now we solve for
    // general linear alpha:
    //   a1(t) = a10*(1-t) + a11*t
    //   a2(t) = a20*(1-t) + a21*t
    //
    // We minimise integral( (a1*D2 + a2*D1)^2 )  subject to
    //   integral(a1 + a2) = 1   (normalisation constraint)
    //
    // The integrand is written as a linear combination of 4 basis
    // functions:
    //   fv(:,0) = (1-t)*D2    (coefficient a10)
    //   fv(:,1) =    t *D2    (coefficient a11)
    //   fv(:,2) = (1-t)*D1    (coefficient a20)
    //   fv(:,3) =    t *D1    (coefficient a21)
    //
    // so that  a1*D2 + a2*D1 = sum_j  x_j * fv(:,j)  where x = [a10,a11,a20,a21].
    {
        gsMatrix<T> fv(N, 4);
        for (index_t i = 0; i < N; ++i)
        {
            T p0 = 1-t_vals(i), p1 = t_vals(i);
            fv(i,0) = p0*D2(i);   // basis function for a10
            fv(i,1) = p1*D2(i);   // basis function for a11
            fv(i,2) = p0*D1(i);   // basis function for a20
            fv(i,3) = p1*D1(i);   // basis function for a21
        }

        // Gram matrix:  G(j,k) = integral( fv(:,j) * fv(:,k) )
        gsMatrix<T> G(4,4); G.setZero();
        for (index_t j = 0; j < 4; ++j)
            for (index_t k = j; k < 4; ++k)
            {
                gsVector<T> pr(N);
                for (index_t i = 0; i < N; ++i) pr(i) = fv(i,j)*fv(i,k);
                G(j,k) = integrate(pr);
                G(k,j) = G(j,k);  // symmetric
            }

        // Normalisation constraint:  integral(a1 + a2) = 1
        // Since a1+a2 = (a10+a20)*(1-t) + (a11+a21)*t,
        // and integral((1-t)) = integral(t) = 0.5 (for t in [0,1]),
        // the constraint is c^T x = 1 with c = [0.5, 0.5, 0.5, 0.5].
        gsVector<T> c(4); c.setConstant(0.5);

        // Solve the KKT system for equality-constrained least squares:
        //   [ 2G   c ] [ x ]   [ 0 ]
        //   [ c^T  0 ] [ λ ] = [ 1 ]
        gsMatrix<T> KKT(5,5); KKT.setZero();
        KKT.block(0,0,4,4) = 2.0*G;
        KKT.block(0,4,4,1) = c;
        KKT.block(4,0,1,4) = c.transpose();
        gsVector<T> rhs(5); rhs.setZero(); rhs(4) = 1.0;

        gsVector<T> sol = KKT.fullPivLu().solve(rhs);
        a10 = sol(0); a11 = sol(1); a20 = sol(2); a21 = sol(3);

        // Check the alpha residual: integral( (a1*D2 + a2*D1)^2 )
        gsVector<T> aerr_v(N);
        for (index_t i = 0; i < N; ++i)
        {
            T a1v = a10*(1-t_vals(i)) + a11*t_vals(i);
            T a2v = a20*(1-t_vals(i)) + a21*t_vals(i);
            T r   = a1v*D2(i) + a2v*D1(i);
            aerr_v(i) = r*r;
        }
        T aerr = integrate(aerr_v);

        if (aerr >= eps)
        {
            gsInfo << "Alpha error too large for interface "
                   << interf << " (err=" << aerr << ").\n";
            return result;
        }

        // Step 2b: General linear beta with gauge constraint
        //
        // Now solve for beta:
        //   b1(t) = b10*(1-t) + b11*t
        //   b2(t) = b20*(1-t) + b21*t
        //
        // Minimise  integral( (b1*D2 - b2*D1 - D3)^2 )
        // subject to the gauge constraint:
        //   integral( a1*b1 - a2*b2 ) = 0
        //
        // This is again written as  b1*D2 - b2*D1  using the same
        // basis functions fv, but with sign flips for the D1 terms:
        //   fv(:,0)*b10 + fv(:,1)*b11 - fv(:,2)*b20 - fv(:,3)*b21 = D3
        //
        // The signs sj = +1 for j<2 (b1 terms) and sj = -1 for j>=2
        // (b2 terms) encode the minus sign in front of b2*D1.
        gsMatrix<T> H(4,4); H.setZero();
        gsVector<T> d_vec(4); d_vec.setZero();
        for (index_t j = 0; j < 4; ++j)
        {
            T sj = (j < 2) ? 1.0 : -1.0;  // +1 for b1, -1 for b2
            for (index_t k = j; k < 4; ++k)
            {
                T sk = (k < 2) ? 1.0 : -1.0;
                // H(j,k) = integral( sj*fv_j * sk*fv_k )
                gsVector<T> pr(N);
                for (index_t i = 0; i < N; ++i)
                    pr(i) = sj*fv(i,j) * sk*fv(i,k);
                H(j,k) = integrate(pr);
                H(k,j) = H(j,k);
            }
            // d_vec(j) = integral( sj*fv_j * D3 )
            gsVector<T> pr(N);
            for (index_t i = 0; i < N; ++i)
                pr(i) = sj*fv(i,j) * D3(i);
            d_vec(j) = integrate(pr);
        }

        // Gauge constraint:  integral( a1*b1 - a2*b2 ) = 0
        //
        // This is a linear constraint in the unknowns [b10, b11, b20, b21]:
        //   e(0) = integral( a1(t) * (1-t) )       -- coefficient of b10
        //   e(1) = integral( a1(t) * t )            -- coefficient of b11
        //   e(2) = integral( -a2(t) * (1-t) )       -- coefficient of b20
        //   e(3) = integral( -a2(t) * t )            -- coefficient of b21
        gsVector<T> e(4);
        {
            gsVector<T> tmp(N);
            for (index_t i = 0; i < N; ++i)
                tmp(i) = (a10*(1-t_vals(i)) + a11*t_vals(i)) * (1-t_vals(i));
            e(0) = integrate(tmp);
            for (index_t i = 0; i < N; ++i)
                tmp(i) = (a10*(1-t_vals(i)) + a11*t_vals(i)) * t_vals(i);
            e(1) = integrate(tmp);
            for (index_t i = 0; i < N; ++i)
                tmp(i) = -(a20*(1-t_vals(i)) + a21*t_vals(i)) * (1-t_vals(i));
            e(2) = integrate(tmp);
            for (index_t i = 0; i < N; ++i)
                tmp(i) = -(a20*(1-t_vals(i)) + a21*t_vals(i)) * t_vals(i);
            e(3) = integrate(tmp);
        }

        // KKT system for the beta solve:
        //   [ 2H   e ] [ b ]   [ 2*d_vec ]
        //   [ e^T  0 ] [ μ ] = [    0    ]
        gsMatrix<T> KKT2(5,5); KKT2.setZero();
        KKT2.block(0,0,4,4) = 2.0*H;
        KKT2.block(0,4,4,1) = e;
        KKT2.block(4,0,1,4) = e.transpose();
        gsVector<T> rhs2(5);
        rhs2.head(4) = 2.0*d_vec;
        rhs2(4) = 0.0;

        gsVector<T> sol2 = KKT2.fullPivLu().solve(rhs2);
        b10 = sol2(0); b11 = sol2(1);
        b20 = sol2(2); b21 = sol2(3);

        // Check the beta residual: integral( (b1*D2 - b2*D1 - D3)^2 )
        gsVector<T> bfit(N);
        for (index_t i = 0; i < N; ++i)
        {
            T b1v = b10*(1-t_vals(i)) + b11*t_vals(i);
            T b2v = b20*(1-t_vals(i)) + b21*t_vals(i);
            T r   = b1v*D2(i) - b2v*D1(i) - D3(i);
            bfit(i) = r*r;
        }
        T berr = integrate(bfit);

        if (berr >= eps)
        {
            gsInfo << "Beta error too large for interface "
                   << interf << " (err=" << berr << ").\n";
            return result;
        }
        success = true;
    }

done:
    // --- Step 6: Pack the result and handle orientation flip ---
    //
    // Result layout: [a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1]
    //
    // The coefficients were computed in the reparametrized tangential
    // parameter v (after reparametrisation).  The reparametrisation
    // preserves the tangential direction in all cases:
    //   west/east  side: tangential = v  (same as original v)
    //   south/north side: tangential = v  (= original u, same direction)
    // Therefore _0 corresponds to original tangential=0 and _1 to
    // tangential=1 for patch 1 (no flip was applied to patch 1).
    //
    // For patch 2: if the boundary-matching step (Step 1b) detected a
    // flip and reversed the v-direction of geo2, the _0 and _1 indices
    // are swapped relative to the original orientation.  We swap them
    // back here so that the returned coefficients always refer to the
    // **original** (input) orientation of the patch.
    result(0) = a10;  result(1) = a11;   // patch 1 alpha (no flip)
    result(2) = b10;  result(3) = b11;   // patch 1 beta  (no flip)
    if (flipped)
    {
        // Undo the flip: swap _0 <-> _1 for patch 2
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
//  computeGluingDataMatrix
// ====================================================================
//
// Purpose:
//   Compute gluing data for ALL interfaces of a multipatch geometry
//   and collect the results into a single matrix for easy access.
//
// Returns:
//   An (nPatches x 16) matrix.  Each row corresponds to a patch.
//   Each row has 4 blocks of 4 columns, one block per side:
//
//     Block s (4 columns):   alpha_0, alpha_1, beta_0, beta_1
//
//   Side index s:  west=0, east=1, south=2, north=3
//   So the column mapping is:
//     col 4*s+0 = alpha at t=0 for side s
//     col 4*s+1 = alpha at t=1 for side s
//     col 4*s+2 = beta  at t=0 for side s
//     col 4*s+3 = beta  at t=1 for side s
//
//   Initialisation: alpha=1, beta=0 for all sides.  This is the
//   identity (trivial) gluing data, appropriate for boundary sides
//   that do not participate in any interface.  Only sides that are
//   part of an interface get overwritten with computed values.
//
// Parameters:
//   mp               [in]  The multipatch geometry.
//   eps              [in]  Tolerance for residual checks (default 1e-8).
//   numGaussPerSpan  [in]  Gauss points per knot span (0 = auto).
// ====================================================================
template <class T>
gsMatrix<T> computeGluingDataMatrix(
    const gsMultiPatch<T> & mp,
    const T eps = 1e-8,
    index_t numGaussPerSpan = 0)
{
    const index_t nP = mp.nPatches();
    gsMatrix<T> mat(nP, 16);

    // Initialise: alpha=1, beta=0 for all 4 sides of every patch.
    // This is the trivial/identity gluing data for boundary sides.
    for (index_t p = 0; p < nP; ++p)
        for (index_t s = 0; s < 4; ++s)
        {
            mat(p, 4*s + 0) = 1.0;  // alpha_0 = 1
            mat(p, 4*s + 1) = 1.0;  // alpha_1 = 1
            mat(p, 4*s + 2) = 0.0;  // beta_0  = 0
            mat(p, 4*s + 3) = 0.0;  // beta_1  = 0
        }

    // Map patchSide enum value to column-block index:
    //   boundary::side enum: west=1, east=2, south=3, north=4
    //   column block index:  0, 1, 2, 3
    // So column block = side_enum_value - 1.
    auto sideToCol = [](const patchSide & ps) -> index_t
    {
        return static_cast<index_t>(ps.side()) - 1;
    };

    // Loop over all interfaces and compute gluing data for each one.
    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const boundaryInterface & interf = *it;
        bool ok = false;
        gsVector<T> gd = computeGluingDataForInterface(mp, interf, ok, eps,
                                                       numGaussPerSpan);
        if (!ok)
        {
            gsInfo << "Failed for interface " << interf << "\n";
            continue;  // leave the default (identity) values for this interface
        }

        // Unpack: gd = [a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1]
        // Store into the matrix at the correct patch row and side columns.
        const patchSide & ps1 = interf.first();
        const patchSide & ps2 = interf.second();
        const index_t col1 = 4 * sideToCol(ps1);  // column offset for patch 1's side
        const index_t col2 = 4 * sideToCol(ps2);  // column offset for patch 2's side

        mat(ps1.patch, col1 + 0) = gd(0);  // a1_0  (alpha_0 for patch 1)
        mat(ps1.patch, col1 + 1) = gd(1);  // a1_1  (alpha_1 for patch 1)
        mat(ps1.patch, col1 + 2) = gd(2);  // b1_0  (beta_0 for patch 1)
        mat(ps1.patch, col1 + 3) = gd(3);  // b1_1  (beta_1 for patch 1)

        mat(ps2.patch, col2 + 0) = gd(4);  // a2_0  (alpha_0 for patch 2)
        mat(ps2.patch, col2 + 1) = gd(5);  // a2_1  (alpha_1 for patch 2)
        mat(ps2.patch, col2 + 2) = gd(6);  // b2_0  (beta_0 for patch 2)
        mat(ps2.patch, col2 + 3) = gd(7);  // b2_1  (beta_1 for patch 2)
    }

    return mat;
}

// ====================================================================
// ====================================================================
//  main
// ====================================================================
//
// Command-line interface:
//   -f <file>      Multipatch geometry XML file (default: two_bicubic_patches)
//   -n <int>       Number of Gauss points per knot span (0 = auto)
//
// Workflow:
//   1. Read the multipatch geometry from the XML file.
//   2. Compute topology (interface detection).
//   3. Call computeGluingDataMatrix to compute all alpha/beta values.
//   4. Print a human-readable summary and the raw matrix.
// ====================================================================
int main(int argc, char* argv[])
{
    using T = double;

    // Default geometry file (can be overridden with -f flag)
    std::string geometry("domain2d/two_bicubic_patches.xml");
    index_t numGaussPerSpan = 0;  // 0 = auto-detect from polynomial degree

    gsCmdLine cmd("Compute AS-G1 gluing data matrix.");
    cmd.addString("f", "file", "G+Smo input multipatch file.", geometry);
    cmd.addInt("n", "numGauss",
               "Number of Gauss points per knot span (0=auto).", numGaussPerSpan);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Read geometry and compute topology (discovers interfaces automatically)
    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
    if (!mpPtr) { gsInfo << "Cannot read " << geometry << ".\n"; return -1; }
    gsMultiPatch<>& mp = *mpPtr;
    mp.computeTopology();

    gsInfo << "Patches: " << mp.nPatches()
           << "  Interfaces: " << mp.nInterfaces() << "\n\n";

    // Compute gluing data matrix for all interfaces
    gsMatrix<T> mat = computeGluingDataMatrix(mp);

    // Print
    const char* sideNames[4] = {"west", "east", "south", "north"};
    // Print human-readable summary: for each patch and each side,
    // show the alpha and beta functions as  f(t) = f_0*(1-t) + f_1*t.
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

    // Also print the raw matrix for programmatic use
    gsInfo << "Raw matrix:\n" << mat << "\n";

    return 0;
}