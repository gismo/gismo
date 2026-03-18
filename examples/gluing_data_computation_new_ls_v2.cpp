/** @file gluing_data_computation_new_ls.cpp

    @brief Compute gluing data via least-squares fitting of linear
           alpha/beta functions to the AS-G1 condition.

    Orientation convention (Left/Right):
      Left patch:  interface at u=0 (west).
                   normal = du (dir 0, points into patch),
                   tangential = dv (dir 1, points "upward").
      Right patch: interface at v=0 (south).
                   normal = dv (dir 1, points into patch),
                   tangential = du (dir 0, points "upward").

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
// Helper functions for reparametrisation
// ====================================================================

/// Swap the two parametric directions (transpose the CP grid)
template <class T>
void swapUV(gsTensorBSpline<2,T> & g)
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
}

/// Reverse direction \a d (0=u, 1=v)
template <class T>
void reverseDir(gsTensorBSpline<2,T> & g, short_t d)
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
}

/// Reparametrize the LEFT patch so that the given \a side becomes
/// "west" (u=0).  After the call:
///   normal direction  = u (index 0) -- points into the patch
///   tangential direction = v (index 1) -- points "upward"
///
/// west  (u=0): nothing
/// east  (u=1): reverse u
/// south (v=0): swap u<->v
/// north (v=1): swap, then reverse u
template <class T>
void reparametriseToWest(gsTensorBSpline<2,T> & geo, const patchSide & ps)
{
    const short_t dir = ps.direction();
    const bool    par = ps.parameter();

    if (dir == 0 && !par)       // west (u=0) -- nothing
    { }
    else if (dir == 0 && par)   // east (u=1) -- reverse u
    { reverseDir(geo, 0); }
    else if (dir == 1 && !par)  // south (v=0) -- swap
    { swapUV(geo); }
    else                        // north (v=1) -- swap, reverse u
    { swapUV(geo); reverseDir(geo, 0); }
}

/// Reparametrize the RIGHT patch so that the given \a side becomes
/// "south" (v=0).  After the call:
///   normal direction  = v (index 1) -- points into the patch
///   tangential direction = u (index 0) -- points "upward"
///
/// south (v=0): nothing
/// north (v=1): reverse v
/// west  (u=0): swap u<->v
/// east  (u=1): swap, then reverse v
template <class T>
void reparametriseToSouth(gsTensorBSpline<2,T> & geo, const patchSide & ps)
{
    const short_t dir = ps.direction();
    const bool    par = ps.parameter();

    if (dir == 1 && !par)       // south (v=0) -- nothing
    { }
    else if (dir == 1 && par)   // north (v=1) -- reverse v
    { reverseDir(geo, 1); }
    else if (dir == 0 && !par)  // west (u=0) -- swap
    { swapUV(geo); }
    else                        // east (u=1) -- swap, reverse v
    { swapUV(geo); reverseDir(geo, 1); }
}


int main(int argc, char* argv[])
{
    using T = double;

    std::string geometry("domain2d/two_bicubic_patches.xml");
    index_t numGaussPerSpan = 0;

    gsCmdLine cmd("Compute gluing data via least-squares.");
    cmd.addString("f", "file", "G+Smo input multi patch file.", geometry);
    cmd.addInt("n", "numGauss",
               "Number of Gauss points per knot span (0=auto).", numGaussPerSpan);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ------------------------------------------------------------------
    // Read geometry
    // ------------------------------------------------------------------
    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
    if (!mpPtr) { gsInfo << "Cannot read " << geometry << ".\n"; return -1; }
    gsMultiPatch<>& mp = *mpPtr;
    mp.computeTopology();

    if (mp.nPatches() != 2)
    { gsInfo << "Need exactly 2 patches, got " << mp.nPatches() << ".\n"; return -1; }
    if (mp.nInterfaces() != 1)
    { gsInfo << "Need exactly 1 interface, got " << mp.nInterfaces() << ".\n"; return -1; }

    // ------------------------------------------------------------------
    // Process every interface
    // ------------------------------------------------------------------
    for (auto it = mp.iBegin(); it < mp.iEnd(); ++it)
    {
        const boundaryInterface & interf = *it;
        gsInfo << "Interface: " << interf << "\n";

        patchSide firstPS  = interf.first();
        patchSide secondPS = interf.second();

        // ==============================================================
        // 1.  Make working copies and reparametrize:
        //       Left patch (geo1):  interface -> west  (u=0)
        //       Right patch (geo2): interface -> south (v=0)
        //
        //     After reparametrisation:
        //       geo1: normDir1=0 (u), tanDir1=1 (v), interface at u=0
        //       geo2: normDir2=1 (v), tanDir2=0 (u), interface at v=0
        // ==============================================================
        gsTensorBSpline<2,T> geo1, geo2;
        {
            //@DR.S DO I NEED THIS CHECK? Maybe rewrite for more general geometries?
            
            const gsTensorBSpline<2,T>* p1 =
                dynamic_cast<const gsTensorBSpline<2,T>*>(&mp.patch(firstPS.patch));
            const gsTensorBSpline<2,T>* p2 =
                dynamic_cast<const gsTensorBSpline<2,T>*>(&mp.patch(secondPS.patch));
            if (!p1 || !p2)
            {
                gsInfo << "Patches must be TensorBSpline2.\n";
                return -1;
            }
            geo1 = *p1;
            geo2 = *p2;
        }

        gsInfo << "Before reparametrisation:\n";
        gsInfo << "  Patch " << firstPS.patch  << " side " << firstPS  << "\n";
        gsInfo << "  Patch " << secondPS.patch << " side " << secondPS << "\n";

        // LEFT patch -> west (u=0)
        reparametriseToWest(geo1, firstPS);
        // RIGHT patch -> south (v=0)
        reparametriseToSouth(geo2, secondPS);

        // Direction indices after reparametrisation
        const index_t normDir1 = 0;  // u for left patch
        const index_t tanDir1  = 1;  // v for left patch
        const index_t normDir2 = 1;  // v for right patch
        const index_t tanDir2  = 0;  // u for right patch

        gsMatrix<T> support1 = geo1.support();
        gsMatrix<T> support2 = geo2.support();

        // ==============================================================
        // 1b. Verify boundary values match
        //     Left patch interface: u=0, parametrized by v
        //     Right patch interface: v=0, parametrized by u
        //     Both tangential parameters should trace the same curve
        // ==============================================================
        const index_t Nq = 5;
        gsMatrix<T> qPts1(2, Nq), qPts2(2, Nq);
        for (index_t i = 0; i < Nq; ++i)
        {
            T tc = static_cast<T>(i) / (Nq - 1);
            // Left patch: u=0, v varies
            qPts1(normDir1, i) = support1(normDir1, 0);  // u=0
            qPts1(tanDir1,  i) = support1(tanDir1, 0)*(1-tc)
                               + support1(tanDir1, 1)*tc;
            // Right patch: v=0, u varies
            qPts2(normDir2, i) = support2(normDir2, 0);  // v=0
            qPts2(tanDir2,  i) = support2(tanDir2, 0)*(1-tc)
                               + support2(tanDir2, 1)*tc;
        }
        gsMatrix<T> qVals1, qVals2;
        geo1.eval_into(qPts1, qVals1);
        geo2.eval_into(qPts2, qVals2);

        const T tol = 1e-10;
        if ((qVals1 - qVals2).norm() < tol)
        {
            gsInfo << "Boundary values match directly.\n";
        }
        else if ((qVals1 - qVals2.rowwise().reverse()).norm() < tol)
        {//NOT SURE IF THIS IS NECESSARY, BUT IN CASE THE TANGENTIAL DIRECTION 
        //OF THE RIGHT PATCH IS REVERSED, WE CAN STILL MATCH THE BOUNDARY VALUES

            gsInfo << "Boundary values match after flipping tangential of patch 2.\n";
            // Reverse the u-direction of geo2 (tangential direction of right patch)
            reverseDir(geo2, tanDir2);
            // Recompute support
            support2 = geo2.support();
        }
        else
        {
            gsInfo << "Boundary values do NOT match!\n";
            gsInfo << "  qVals1:\n" << qVals1 << "\n";
            gsInfo << "  qVals2:\n" << qVals2 << "\n";
            return -1;
        }

        // ==============================================================
        // 2.  Set up Gauss quadrature along the interface
        //     The interface parameter runs along:
        //       v for geo1 (tanDir1=1), u for geo2 (tanDir2=0)
        //     We use the tangential parameter range of geo1 (v).
        // ==============================================================
        const short_t deg1 = geo1.basis().degree(tanDir1);
        const short_t deg2 = geo2.basis().degree(tanDir2);
        const short_t maxDeg = std::max(deg1, deg2);
        const index_t nGauss = numGaussPerSpan > 0
                                   ? numGaussPerSpan
                                   : 2*maxDeg + 1;

        // Merge knot breaks from tangential directions of both patches
        std::vector<T> breaks1, breaks2;
        {
            const gsBSplineBasis<T>* bb =
                dynamic_cast<const gsBSplineBasis<T>*>(&geo1.basis().component(tanDir1));
            if (bb) breaks1 = bb->knots().breaks();
            else  { breaks1.push_back(support1(tanDir1,0));
                    breaks1.push_back(support1(tanDir1,1)); }
        }
        {
            const gsBSplineBasis<T>* bb =
                dynamic_cast<const gsBSplineBasis<T>*>(&geo2.basis().component(tanDir2));
            if (bb) breaks2 = bb->knots().breaks();
            else  { breaks2.push_back(support2(tanDir2,0));
                    breaks2.push_back(support2(tanDir2,1)); }
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

        // ==============================================================
        // 3.  Evaluate geometry derivatives at quadrature points
        //     Left patch: pts at (u=0, v=t)
        //     Right patch: pts at (u=t, v=0)
        // ==============================================================
        gsMatrix<T> pts1(2, N), pts2(2, N);
        for (index_t i = 0; i < N; ++i)
        {
            // Left patch: u=0, v=gaussNodes
            pts1(normDir1, i) = support1(normDir1, 0);  // u=0
            pts1(tanDir1,  i) = gaussNodes(0, i);
            // Right patch: v=0, u=gaussNodes
            pts2(normDir2, i) = support2(normDir2, 0);  // v=0
            pts2(tanDir2,  i) = gaussNodes(0, i);
        }

        gsMatrix<T> derivs1, derivs2;
        geo1.deriv_into(pts1, derivs1);
        geo2.deriv_into(pts2, derivs2);

        // Extract 2D partial derivative vector from deriv_into output
        //   deriv_into layout: row 0 = d_u f_x,  row 1 = d_v f_x,
        //                      row 2 = d_u f_y,  row 3 = d_v f_y
        auto getPartial = [](const gsMatrix<T>& d, index_t dir, index_t i) -> gsVector<T>
        {
            gsVector<T> v(2);
            v(0) = d(dir,     i);   // d_dir f_x
            v(1) = d(dir + 2, i);   // d_dir f_y
            return v;
        };

        // 2D cross-determinant
        auto det2 = [](const gsVector<T>& a, const gsVector<T>& b) -> T
        { return a(0)*b(1) - a(1)*b(0); };

        // Parameter values (for linear basis phi0=1-t, phi1=t)
        gsVector<T> t_vals(N);
        for (index_t i = 0; i < N; ++i)
            t_vals(i) = gaussNodes(0, i);

        // ==============================================================
        // 4.  Compute D1, D2, D3
        //     D1 = det(d_n G1, d_t G1) = det(d_u G1, d_v G1)
        //     D2 = det(d_n G2, d_t G2) = det(d_v G2, d_u G2)
        //     D3 = det(d_n G1, d_n G2) = det(d_u G1, d_v G2)
        // ==============================================================
        gsVector<T> D1(N), D2(N), D3(N);
        for (index_t i = 0; i < N; ++i)
        {
            gsVector<T> dG1dn = getPartial(derivs1, normDir1, i);  // d_u G1
            gsVector<T> dG1dt = getPartial(derivs1, tanDir1,  i);  // d_v G1
            gsVector<T> dG2dn = getPartial(derivs2, normDir2, i);  // d_v G2
            gsVector<T> dG2dt = getPartial(derivs2, tanDir2,  i);  // d_u G2

            D1(i) = det2(dG1dn, dG1dt);   // det(d_u G1, d_v G1)
            D2(i) = det2(dG2dn, dG2dt);   // det(d_v G2, d_u G2)
            D3(i) = det2(dG1dn, dG2dn);   // det(d_u G1, d_v G2)
        }

        // Check that D1, D2 don't change sign
        if (D1.minCoeff() * D1.maxCoeff() < 0 ||
            D2.minCoeff() * D2.maxCoeff() < 0)
        {
            gsInfo << "WARNING: D1 or D2 changes sign -- geometry may be degenerate.\n";
            continue;
        }

        gsInfo << "D1 range: [" << D1.minCoeff() << ", " << D1.maxCoeff() << "]\n";
        gsInfo << "D2 range: [" << D2.minCoeff() << ", " << D2.maxCoeff() << "]\n";
        gsInfo << "D3 range: [" << D3.minCoeff() << ", " << D3.maxCoeff() << "]\n";

        // Weighted L2 integral via Gauss quadrature
        auto integrate = [&](const gsVector<T>& f) -> T { return w.dot(f); };

        // ==============================================================
        // 5.  Compute gluing data: alpha and beta
        //
        //   AS-G1 vector condition (new convention):
        //     a2 * d_u G1  +  a1 * d_v G2  +  beta * d_v G1 = 0
        //   where beta = a1*b2 + a2*b1.
        //
        //   Determinant conditions (same form as before):
        //     a1*D2 + a2*D1 = 0         (alpha condition)
        //     b1*D2 - b2*D1 = D3        (beta condition)
        //   with constraints:
        //     integral(a1 + a2) = 1
        //     integral(a1*b1 - a2*b2) = 0  (gauge, Step 2b only)
        // ==============================================================
        const T eps = 1e-8;
        T a10, a11, a20, a21;
        T b10, b11, b20, b21;
        bool success = false;

        // ------ Step 1: constant alpha ------
        gsInfo << "\n--- Step 1: constant alpha ---\n";
        {
            T intD1D1 = integrate((D1.array()*D1.array()).matrix());
            T intD1D2 = integrate((D1.array()*D2.array()).matrix());
            T intD2D2 = integrate((D2.array()*D2.array()).matrix());

            // Minimize integral(a1*D2 + a2*D1)^2 with a1+a2=1
            //   = integral(a1*(D2-D1) + D1)^2
            // d/da1 = 0 => a1 = (intD1D1 - intD1D2) / (intD1D1 - 2*intD1D2 + intD2D2)
            T denom = intD1D1 - 2*intD1D2 + intD2D2;
            T a1c, a2c;
            T err;
            if (std::abs(denom) < 1e-30)
            {
                a1c = 0.5; a2c = 0.5;
                err = intD1D1;
            }
            else
            {
                a1c = (intD1D1 - intD1D2) / denom;
                a2c = 1.0 - a1c;
                err = a1c*a1c*intD2D2 + 2*a1c*a2c*intD1D2 + a2c*a2c*intD1D1;
            }
            gsInfo << "  alpha_1=" << a1c << "  alpha_2=" << a2c
                   << "  err=" << err << "\n";

            if (err < eps)
            {
                a10 = a1c; a11 = a1c;
                a20 = a2c; a21 = a2c;

                // -- Step 1b: linear beta with b1=b2 --
                //   b1*D2 - b2*D1 = D3 with b1=b2=b: b*(D2-D1) = D3
                gsInfo << "--- Step 1b: linear beta (b1=b2) ---\n";
                gsVector<T> S = D2 - D1;
                gsMatrix<T> Ab(2,2); Ab.setZero();
                gsVector<T> rb(2);   rb.setZero();
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
                b10 = bc(0); b11 = bc(1);
                b20 = b10;   b21 = b11;

                // Fitting error
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

                if (berr < eps)
                { gsInfo << "  SUCCESS (Step 1).\n"; success = true; }
                else
                { gsInfo << "  Beta fit too large -> Step 2.\n"; }
            }
            else
            { gsInfo << "  Error too large -> Step 2.\n"; }
        }

        // ------ Step 2: linear alpha ------
        if (!success)
        {
            gsInfo << "\n--- Step 2: linear alpha ---\n";

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

            // Gram matrix G_jk = integral(f_j * f_k)
            gsMatrix<T> G(4,4); G.setZero();
            for (index_t j = 0; j < 4; ++j)
                for (index_t k = j; k < 4; ++k)
                {
                    gsVector<T> pr(N);
                    for (index_t i = 0; i < N; ++i) pr(i) = fv(i,j)*fv(i,k);
                    G(j,k) = integrate(pr);
                    G(k,j) = G(j,k);
                }

            // Constraint: integral(a1+a2) = 1  =>  c=[1/2,1/2,1/2,1/2], c^T x = 1
            gsVector<T> c(4); c.setConstant(0.5);

            // KKT system
            gsMatrix<T> KKT(5,5); KKT.setZero();
            KKT.block(0,0,4,4) = 2.0*G;
            KKT.block(0,4,4,1) = c;
            KKT.block(4,0,1,4) = c.transpose();
            gsVector<T> rhs(5); rhs.setZero(); rhs(4) = 1.0;

            gsInfo << "KKT matrix:\n" << KKT << "\n";

            gsVector<T> sol = KKT.fullPivLu().solve(rhs);
            a10 = sol(0); a11 = sol(1); a20 = sol(2); a21 = sol(3);

            // Compute error
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

            if (aerr < eps)
            {
                // -- Step 2b: linear b1, b2 --
                gsInfo << "--- Step 2b: linear beta ---\n";

                // Basis: g = [phi0*D2, phi1*D2, -phi0*D1, -phi1*D1]
                //   b1*D2 - b2*D1 = b10*(phi0*D2) + b11*(phi1*D2)
                //                 + b20*(-phi0*D1) + b21*(-phi1*D1)
                gsMatrix<T> gv(N, 4);
                for (index_t i = 0; i < N; ++i)
                {
                    T p0 = 1-t_vals(i), p1 = t_vals(i);
                    gv(i,0) =  p0*D2(i);
                    gv(i,1) =  p1*D2(i);
                    gv(i,2) = -p0*D1(i);
                    gv(i,3) = -p1*D1(i);
                }

                gsMatrix<T> H(4,4); H.setZero();
                gsVector<T> d_vec(4); d_vec.setZero();
                for (index_t j = 0; j < 4; ++j)
                {
                    for (index_t k = j; k < 4; ++k)
                    {
                        gsVector<T> pr(N);
                        for (index_t i = 0; i < N; ++i) pr(i) = gv(i,j)*gv(i,k);
                        H(j,k) = integrate(pr);
                        H(k,j) = H(j,k);
                    }
                    gsVector<T> pr(N);
                    for (index_t i = 0; i < N; ++i) pr(i) = gv(i,j)*D3(i);
                    d_vec(j) = integrate(pr);
                }

                // Constraint: integral(a1*b1 - a2*b2) = 0
                gsVector<T> e(4);
                {
                    gsVector<T> tmp(N);
                    for (index_t i = 0; i < N; ++i)
                    {
                        T p0 = 1-t_vals(i), p1 = t_vals(i);
                        T a1 = a10*p0 + a11*p1;
                        tmp(i) = a1*p0;
                    }
                    e(0) = integrate(tmp);
                    for (index_t i = 0; i < N; ++i)
                        tmp(i) = (a10*(1-t_vals(i))+a11*t_vals(i))*t_vals(i);
                    e(1) = integrate(tmp);
                    for (index_t i = 0; i < N; ++i)
                    {
                        T p0 = 1-t_vals(i), p1 = t_vals(i);
                        T a2 = a20*p0 + a21*p1;
                        tmp(i) = -a2*p0;
                    }
                    e(2) = integrate(tmp);
                    for (index_t i = 0; i < N; ++i)
                        tmp(i) = -(a20*(1-t_vals(i))+a21*t_vals(i))*t_vals(i);
                    e(3) = integrate(tmp);
                }

                // KKT: [2H e; e^T 0] [y;mu] = [2d;0]
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

                // Fitting error
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

                if (berr < eps) { gsInfo << "  SUCCESS (Step 2).\n"; success = true; }
                else            { gsInfo << "  Beta fit too large -- not AS-G1.\n"; }
            }
            else
            { gsInfo << "  Alpha error too large -- not AS-G1.\n"; }
        }

        if (!success)
        {
            gsInfo << "\nFailed: geometry is not AS-G1.\n";
            continue;
        }

        // ==============================================================
        // 6.  Output
        // ==============================================================
        gsInfo << "\n=== Gluing data ===\n";
        gsInfo << "  alpha_1(t) = " << a10 << "*(1-t) + " << a11 << "*t\n";
        gsInfo << "  alpha_2(t) = " << a20 << "*(1-t) + " << a21 << "*t\n";
        gsInfo << "  beta_1(t)  = " << b10 << "*(1-t) + " << b11 << "*t\n";
        gsInfo << "  beta_2(t)  = " << b20 << "*(1-t) + " << b21 << "*t\n";

        // ==============================================================
        // 7.  Verification at uniform check points
        //
        //   The AS-G1 vector condition (new convention):
        //     a2 * d_u G1  +  a1 * d_v G2  +  beta * d_v G1 = 0
        //   where beta = a1*b2 + a2*b1
        //
        //   Left patch evaluated at u=0, right patch at v=0.
        //   d_u G1 = normal of left,  d_v G1 = tangential of left
        //   d_v G2 = normal of right, d_u G2 = tangential of right
        // ==============================================================
        gsInfo << "\n=== Verification ===\n";
        const index_t Nchk = 21;
        T maxRes = 0;

        gsMatrix<T> cPts1(2, Nchk), cPts2(2, Nchk);
        for (index_t i = 0; i < Nchk; ++i)
        {
            T tc = static_cast<T>(i) / (Nchk - 1);
            // Left patch: u=0, v varies
            cPts1(normDir1, i) = support1(normDir1, 0);
            cPts1(tanDir1,  i) = support1(tanDir1, 0)*(1-tc)
                               + support1(tanDir1, 1)*tc;
            // Right patch: v=0, u varies
            cPts2(normDir2, i) = support2(normDir2, 0);
            cPts2(tanDir2,  i) = support2(tanDir2, 0)*(1-tc)
                               + support2(tanDir2, 1)*tc;
        }

        gsMatrix<T> cD1, cD2;
        geo1.deriv_into(cPts1, cD1);
        geo2.deriv_into(cPts2, cD2);

        for (index_t i = 0; i < Nchk; ++i)
        {
            T tc = static_cast<T>(i) / (Nchk - 1);
            T p0 = 1-tc, p1 = tc;
            T a1 = a10*p0 + a11*p1;
            T a2 = a20*p0 + a21*p1;
            T b1 = b10*p0 + b11*p1;
            T b2 = b20*p0 + b21*p1;
            T beta_comb = a1*b2 + a2*b1;

            gsVector<T> dG1dn = getPartial(cD1, normDir1, i);  // d_u G1
            gsVector<T> dG1dt = getPartial(cD1, tanDir1,  i);  // d_v G1
            gsVector<T> dG2dn = getPartial(cD2, normDir2, i);  // d_v G2

            // AS-G1: a2*d_u(G1) + a1*d_v(G2) + beta*d_v(G1) = 0
            gsVector<T> res = a2*dG1dn + a1*dG2dn + beta_comb*dG1dt;
            T rn = res.norm();
            maxRes = std::max(maxRes, rn);

            if (Nchk <= 21)
                gsInfo << "  t=" << tc << "  |res|=" << rn << "\n";
        }

        gsInfo << "\n  Max residual = " << maxRes << "\n";
        if (maxRes < 1e-10)
            gsInfo << "  => AS-G1 SATISFIED (exact).\n";
        else if (maxRes < 1e-4)
            gsInfo << "  => AS-G1 SATISFIED.\n";
        else
            gsInfo << "  => AS-G1 NOT satisfied.\n";
    }

    return 0;
}
