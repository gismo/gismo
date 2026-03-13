/** @file get_gluing_data_least_squares.cpp

    @brief Compute gluing data via least-squares fitting of linear
           alpha/beta functions to the AS-G1 condition.

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

/// Reparametrize a 2D tensor B-spline patch so that the given \a side
/// becomes the "west" side (u=0).  After the call the normal direction
/// is u (index 0) and the tangential direction is v (index 1).
///
/// For south (v=0): swap u<->v, then reverse the new u so that the
///   former south edge sits at u=0 with outward normal -u.
/// For east (u=1): reverse u so east becomes west.
/// For north (v=1): swap u<->v, then the former north edge sits at
///   u=1 in the new layout, so reverse u.
/// For west (u=0): nothing to do.
template <class T>
void reparametrizeToWest(gsTensorBSpline<2,T> & geo, const patchSide & ps)
{
    const short_t dir = ps.direction();
    const bool    par = ps.parameter();   // false=0, true=1

    // Helper: swap the two parametric directions (transpose the CP grid)
    auto swapUV = [](gsTensorBSpline<2,T> & g)
    {
        const index_t nu = g.basis().size(0);
        const index_t nv = g.basis().size(1);
        const index_t d  = g.geoDim();
        gsMatrix<T> oldC = g.coefs();          // nu*nv × d
        gsMatrix<T> newC(nv*nu, d);
        for (index_t j = 0; j < nv; ++j)
            for (index_t i = 0; i < nu; ++i)
                newC.row(j + i*nv) = oldC.row(i + j*nu);

        gsKnotVector<T> kvU = g.knots(0);
        gsKnotVector<T> kvV = g.knots(1);
        g = gsTensorBSpline<2,T>(kvV, kvU, newC);
    };

    // Helper: reverse direction \a d (0=u, 1=v)
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

    if (dir == 0 && !par)       // west  (u=0) — already correct
    {
        // nothing
    }
    else if (dir == 0 && par)   // east  (u=1) — reverse u
    {
        reverseDir(geo, 0);
    }
    else if (dir == 1 && !par)  // south (v=0) — swap u<->v
    {                           // old south (v=0) → new west (u=0)
        swapUV(geo);
    }
    else                        // north (v=1) — swap, then reverse u
    {                           // old north (v=1) → new east (u=1)
        swapUV(geo);            // → reverse u to get west (u=0)
        reverseDir(geo, 0);
    }
}


int main(int argc, char* argv[])
{
    using T = double;

    std::string geometry("domain2d/two_bicubic_patches.xml");
    index_t numGaussPerSpan = 0;  // 0 = automatic (deg+1 per span)

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
        // 1.  Make working copies and reparametrize so that for both
        //     patches the interface sits at u=0 (west).
        //     After reparametrisation:
        //       normal  direction = u  (index 0)
        //       tangent direction = v  (index 1)
        // ==============================================================
        gsTensorBSpline<2,T> geo1, geo2;
        {
            // Try to cast; works for TensorBSpline2 geometries
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

        reparametrizeToWest(geo1, firstPS);
        reparametrizeToWest(geo2, secondPS);

        // After reparametrisation both interfaces are at u=0 (west).
        // Normal direction = 0 (u), tangential direction = 1 (v).
        const index_t normDir = 0;
        const index_t tanDir  = 1;

        // Verify that boundary values match (possibly after flipping v)
        gsMatrix<T> support1 = geo1.support();
        gsMatrix<T> support2 = geo2.support();

        // Quick check at 5 uniform points
        const index_t Nq = 5;
        gsMatrix<T> qPts1(2, Nq), qPts2(2, Nq);
        for (index_t i = 0; i < Nq; ++i)
        {
            T tc = static_cast<T>(i) / (Nq - 1);
            qPts1(tanDir,  i) = support1(tanDir, 0)*(1-tc) + support1(tanDir, 1)*tc;
            qPts1(normDir, i) = support1(normDir, 0);  // u=0
            qPts2(tanDir,  i) = support2(tanDir, 0)*(1-tc) + support2(tanDir, 1)*tc;
            qPts2(normDir, i) = support2(normDir, 0);  // u=0
        }
        gsMatrix<T> qVals1, qVals2;
        geo1.eval_into(qPts1, qVals1);
        geo2.eval_into(qPts2, qVals2);

        bool flipped = false;
        (void)flipped; // may be used later
        const T tol = 1e-10;
        if ((qVals1 - qVals2).norm() < tol)
        {
            gsInfo << "Boundary values match directly.\n";
        }
        else if ((qVals1 - qVals2.rowwise().reverse()).norm() < tol)
        {
            gsInfo << "Boundary values match after flipping v of patch 2.\n";
            // Reverse the v-direction of geo2
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
            flipped = true;
        }
        else
        {
            gsInfo << "Boundary values do NOT match!\n";
            gsInfo << "  qVals1:\n" << qVals1 << "\n";
            gsInfo << "  qVals2:\n" << qVals2 << "\n";
            return -1;
        }

        // ==============================================================
        // 2.  Set up Gauss quadrature along the interface (v-direction)
        // ==============================================================
        const short_t deg1 = geo1.basis().degree(tanDir);
        const short_t deg2 = geo2.basis().degree(tanDir);
        const short_t maxDeg = std::max(deg1, deg2);
        // We need to integrate products like (α·D)² where α is linear
        // (deg 1) and D has degree 2*geoDeg-1 from Jacobian determinants.
        // The Gram matrix requires exact integration of degree
        // 2*(1 + 2*maxDeg - 1) = 4*maxDeg.  Gauss rule with n points
        // integrates degree 2n-1 exactly, so n ≥ (4*maxDeg+1)/2.
        // We use 2*maxDeg + 1 to be safe.
        const index_t nGauss = numGaussPerSpan > 0
                                   ? numGaussPerSpan
                                   : 2*maxDeg + 1;

        // Merge knot breaks from both patches' tangential direction
        std::vector<T> breaks1, breaks2;
        {
            const gsBSplineBasis<T>* bb =
                dynamic_cast<const gsBSplineBasis<T>*>(&geo1.basis().component(tanDir));
            if (bb) breaks1 = bb->knots().breaks();
            else  { breaks1.push_back(support1(tanDir,0));
                    breaks1.push_back(support1(tanDir,1)); }
        }
        {
            const gsBSplineBasis<T>* bb =
                dynamic_cast<const gsBSplineBasis<T>*>(&geo2.basis().component(tanDir));
            if (bb) breaks2 = bb->knots().breaks();
            else  { breaks2.push_back(support2(tanDir,0));
                    breaks2.push_back(support2(tanDir,1)); }
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
               << " spans × " << nGauss << " pts = " << N << " total\n";

        // ==============================================================
        // 3.  Evaluate geometry derivatives at quadrature points
        // ==============================================================
        gsMatrix<T> pts1(2, N), pts2(2, N);
        for (index_t i = 0; i < N; ++i)
        {
            pts1(tanDir,  i) = gaussNodes(0, i);
            pts1(normDir, i) = support1(normDir, 0);  // u=0
            pts2(tanDir,  i) = gaussNodes(0, i);
            pts2(normDir, i) = support2(normDir, 0);  // u=0
        }

        gsMatrix<T> derivs1, derivs2;
        geo1.deriv_into(pts1, derivs1);
        geo2.deriv_into(pts2, derivs2);

        // Extract 2D partial derivative vector from deriv_into output
        //   deriv_into layout: row 0 = ∂_u f_x,  row 1 = ∂_v f_x,
        //                      row 2 = ∂_u f_y,  row 3 = ∂_v f_y
        auto getPartial = [](const gsMatrix<T>& d, index_t dir, index_t i) -> gsVector<T>
        {
            gsVector<T> v(2);
            v(0) = d(dir,     i);   // ∂_dir f_x
            v(1) = d(dir + 2, i);   // ∂_dir f_y
            return v;
        };

        // 2D cross-determinant
        auto det2 = [](const gsVector<T>& a, const gsVector<T>& b) -> T
        { return a(0)*b(1) - a(1)*b(0); };

        // Parameter values (for linear basis φ₀=1-t, φ₁=t)
        gsVector<T> t_vals(N);
        for (index_t i = 0; i < N; ++i)
            t_vals(i) = gaussNodes(0, i);

        // ==============================================================
        // 4.  Compute D1, D2, D3
        //     D1 = |∂_n G1, ∂_t G1|   (Jacobian of patch 1)
        //     D2 = |∂_n G2, ∂_t G2|   (Jacobian of patch 2)
        //     D3 = |∂_n G1, ∂_n G2|   (cross of normal derivatives)
        // ==============================================================
        gsVector<T> D1(N), D2(N), D3(N);
        for (index_t i = 0; i < N; ++i)
        {
            gsVector<T> dG1dn = getPartial(derivs1, normDir, i);
            gsVector<T> dG1dt = getPartial(derivs1, tanDir,  i);
            gsVector<T> dG2dn = getPartial(derivs2, normDir, i);
            gsVector<T> dG2dt = getPartial(derivs2, tanDir,  i);

            D1(i) = det2(dG1dn, dG1dt);
            D2(i) = det2(dG2dn, dG2dt);
            D3(i) = det2(dG1dn, dG2dn);
        }

        // D1, D2 are signed Jacobian determinants.  They may be
        // negative (when the parametrisation is left-handed).
        // We check that they don't change sign (no degeneracy).
        if (D1.minCoeff() * D1.maxCoeff() < 0 ||
            D2.minCoeff() * D2.maxCoeff() < 0)
        {
            gsInfo << "WARNING: D1 or D2 changes sign — geometry may be degenerate.\n";
            continue;
        }

        gsInfo << "D1 range: [" << D1.minCoeff() << ", " << D1.maxCoeff() << "]\n";
        gsInfo << "D2 range: [" << D2.minCoeff() << ", " << D2.maxCoeff() << "]\n";
        gsInfo << "D3 range: [" << D3.minCoeff() << ", " << D3.maxCoeff() << "]\n";

        // Weighted L² integral via Gauss quadrature
        auto integrate = [&](const gsVector<T>& f) -> T { return w.dot(f); };

        // ==============================================================
        // 5.  Compute gluing data: alpha and beta
        //
        //   AS-G1 determinant conditions (with raw signed D values):
        //     α₁ D₂ + α₂ D₁ = 0        (alpha condition)
        //     β₁ D₂ − β₂ D₁ = D₃        (beta  condition)
        //   with constraints:
        //     ∫(α₁ + α₂) = 1
        //     ∫(α₁ β₁ − α₂ β₂) = 0     (gauge, Step 2b only)
        // ==============================================================
        const T eps = 1e-8;
        T a10, a11, a20, a21;
        T b10, b11, b20, b21;
        bool success = false;

        // ------ Step 1: constant alpha ------
        // Condition: α₁D₂ + α₂D₁ = 0 with α₁+α₂ = 1
        // → α₁D₂ + (1-α₁)D₁ = 0  → α₁(D₂+D₁) = -D₁ ... wait
        // Actually minimize ∫(α₁D₂ + α₂D₁)² with α₂=1-α₁
        //   = ∫(α₁(D₂-D₁) + D₁)²
        // ∂/∂α₁ = 0 → α₁ = (∫D₁² - ∫D₁D₂) / ∫(D₂-D₁)²
        gsInfo << "\n--- Step 1: constant alpha ---\n";
        {
            T intD1D1 = integrate((D1.array()*D1.array()).matrix());
            T intD1D2 = integrate((D1.array()*D2.array()).matrix());
            T intD2D2 = integrate((D2.array()*D2.array()).matrix());

            T denom = intD1D1 - 2*intD1D2 + intD2D2;   // ∫(D₁-D₂)²
            T a1c, a2c;
            T err;
            if (std::abs(denom) < 1e-30)
            {
                // D₁ = D₂ everywhere, α₁D₂+α₂D₁ = (α₁+α₂)D₁ = D₁ ≠ 0
                // Cannot satisfy unless D₁=0
                a1c = 0.5; a2c = 0.5;
                err = intD1D1;   // won't be < eps
            }
            else
            {
                a1c = (intD1D1 - intD1D2) / denom;
                a2c = 1.0 - a1c;
                // error = ∫(α₁D₂ + α₂D₁)²
                err = a1c*a1c*intD2D2 + 2*a1c*a2c*intD1D2 + a2c*a2c*intD1D1;
            }
            gsInfo << "  alpha_1=" << a1c << "  alpha_2=" << a2c
                   << "  err=" << err << "\n";

            if (err < eps)
            {
                a10 = a1c; a11 = a1c;
                a20 = a2c; a21 = a2c;

                // -- Step 1b: linear beta with β₁=β₂ --
                //   With β₁=β₂=β: β₁D₂−β₂D₁ = β(D₂−D₁)
                //   Fit: ‖β(D₂−D₁) − D₃‖² → min
                gsInfo << "--- Step 1b: linear beta (β₁=β₂) ---\n";
                gsVector<T> S = D2 - D1;    // D₂ − D₁
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

            // Basis vectors: f = [φ₀D₂, φ₁D₂, φ₀D₁, φ₁D₁]
            //   α₁D₂ + α₂D₁ = a10(φ₀D₂) + a11(φ₁D₂) + a20(φ₀D₁) + a21(φ₁D₁)
            gsMatrix<T> fv(N, 4);
            for (index_t i = 0; i < N; ++i)
            {
                T p0 = 1-t_vals(i), p1 = t_vals(i);
                fv(i,0) = p0*D2(i);
                fv(i,1) = p1*D2(i);
                fv(i,2) = p0*D1(i);
                fv(i,3) = p1*D1(i);
            }

            // Gram matrix G_jk = ∫ f_j f_k
            gsMatrix<T> G(4,4); G.setZero();
            for (index_t j = 0; j < 4; ++j)
                for (index_t k = j; k < 4; ++k)
                {
                    gsVector<T> pr(N);
                    for (index_t i = 0; i < N; ++i) pr(i) = fv(i,j)*fv(i,k);
                    G(j,k) = integrate(pr);
                    G(k,j) = G(j,k);
                }

            // Constraint: ∫(α₁+α₂) = 1  ⟹  c=[½,½,½,½], c^T x = 1
            gsVector<T> c(4); c.setConstant(0.5);

            // KKT system: [2G c; c^T 0] [x;λ] = [0;1]
            gsMatrix<T> KKT(5,5); KKT.setZero();
            KKT.block(0,0,4,4) = 2.0*G;
            KKT.block(0,4,4,1) = c;
            KKT.block(4,0,1,4) = c.transpose();
            gsVector<T> rhs(5); rhs.setZero(); rhs(4) = 1.0;

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
                // -- Step 2b: linear β₁, β₂ --
                //   ‖β₁D₂ − β₂D₁ − D₃‖² → min
                //   with ∫(α₁β₁ − α₂β₂) = 0
                gsInfo << "--- Step 2b: linear beta ---\n";

                // Basis: g = [φ₀D₂, φ₁D₂, −φ₀D₁, −φ₁D₁]
                //   β₁D₂−β₂D₁ = b10(φ₀D₂) + b11(φ₁D₂) + b20(−φ₀D₁) + b21(−φ₁D₁)
                //              = b10(φ₀D₂) + b11(φ₁D₂) - b20(φ₀D₁) - b21(φ₁D₁)
                //   So b20,b21 enter with a minus sign for D₁.
                gsMatrix<T> gv(N, 4);
                for (index_t i = 0; i < N; ++i)
                {
                    T p0 = 1-t_vals(i), p1 = t_vals(i);
                    gv(i,0) =  p0*D2(i);   // β₁ coeff (1-t) × D₂
                    gv(i,1) =  p1*D2(i);   // β₁ coeff t     × D₂
                    gv(i,2) = -p0*D1(i);   // β₂ coeff (1-t) × (−D₁)
                    gv(i,3) = -p1*D1(i);   // β₂ coeff t     × (−D₁)
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

                // Constraint: ∫(α₁β₁ − α₂β₂) = 0
                //   = b10∫(α₁φ₀) + b11∫(α₁φ₁) − b20∫(α₂φ₀) − b21∫(α₂φ₁)
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

                // KKT: [2H e; e^T 0] [y;μ] = [2d;0]
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
                else            { gsInfo << "  Beta fit too large — not AS-G1.\n"; }
            }
            else
            { gsInfo << "  Alpha error too large — not AS-G1.\n"; }
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
        //   The AS-G1 vector condition (using the reparametrised patches
        //   with interface at u=0 for both):
        //
        //     α₂·∂_u G₁ + α₁·∂_u G₂ + (α₁β₂+α₂β₁)·∂_v G₁ = 0
        //
        //   The determinant form (with raw signed D values) gives:
        //     α₁D₂ + α₂D₁ = 0   and   β₁D₂ − β₂D₁ = D₃
        // ==============================================================
        gsInfo << "\n=== Verification ===\n";
        const index_t Nchk = 21;
        T maxRes = 0;

        gsMatrix<T> cPts1(2, Nchk), cPts2(2, Nchk);
        for (index_t i = 0; i < Nchk; ++i)
        {
            T tc = static_cast<T>(i) / (Nchk - 1);
            cPts1(tanDir,  i) = support1(tanDir, 0)*(1-tc) + support1(tanDir, 1)*tc;
            cPts1(normDir, i) = support1(normDir, 0);
            cPts2(tanDir,  i) = support2(tanDir, 0)*(1-tc) + support2(tanDir, 1)*tc;
            cPts2(normDir, i) = support2(normDir, 0);
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

            gsVector<T> dG1dn = getPartial(cD1, normDir, i);
            gsVector<T> dG1dt = getPartial(cD1, tanDir,  i);
            gsVector<T> dG2dn = getPartial(cD2, normDir, i);

            // AS-G1: α₂·∂_u G₁ + α₁·∂_u G₂ + β·∂_v G₁ = 0  (raw normals)
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

