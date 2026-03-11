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
#include <map>
#include <algorithm>
#include <numeric>
#include <gismo.h>

using namespace gismo;


int main(int argc, char* argv[])
{
    using T = double;
    

    //std::string geometry("domain2d/simple_two_patch_geometry.xml");
    //std::string geometry("domain2d/two_bicubic_patches.xml");
    std::string geometry("domain2d/two_patches_from_gluing_data.xml");
    //std::string geometry("domain2d/two_bilinear_patches.xml");
    index_t numGaussPerSpan = 0;  // 0 = automatic (deg+1 per span)
   
    gsCmdLine cmd("Compute gluing data via least-squares.");
    cmd.addString("f", "file", "G+Smo input multi patch file.", geometry);
    cmd.addInt("n", "numGauss", "Number of Gauss points per knot span (0=auto).", numGaussPerSpan);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    
    // Read geometry
    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
    if (!mpPtr)
    {
        gsInfo << "No geometry found in file " << geometry << ".\n";
        return -1;
    }
    gsMultiPatch<>& mp = *mpPtr;
    
    mp.computeTopology(); // if the xml file does not have a topology, we need to compute it first
    
    if (mp.nPatches() != 2)
    {
        gsInfo << "The file contains " << mp.nPatches() << " patches. This is not ok!\n";
        return -1;
    }
    
    if (mp.nInterfaces() != 1)
    {
        gsInfo << "The mp contains " << mp.nInterfaces() << " interfaces.  This is not ok!\n";
        return -1;
    }
    
    for (auto it = mp.iBegin(); it<mp.iEnd(); ++it)
    {
        const boundaryInterface& interf = *it;
        gsInfo << "The one and only interfaces is:\n" << interf << "\n";
        
        patchSide firstPS = interf.first();
        patchSide secondPS = interf.second();
        
        
        // A. First patch
        const gsGeometry<>& firstGeo = mp.patch(firstPS.patch);
        gsInfo << "firstGeo is " << firstGeo << "\n";

        // B. Second patch
        const gsGeometry<>& secondGeo = mp.patch(secondPS.patch);
        gsInfo << "secondGeo is " << secondGeo << "\n";
        
        // ================================================================
        // Set up Gauss quadrature along the interface
        //
        // We merge the knot break points from both patches' tangential
        // directions so the quadrature rule respects all knot spans.
        // With enough Gauss points per span, integrals of polynomial
        // products become exact.
        // ================================================================
        const index_t normDir1 = firstPS.direction();
        const index_t ambDir1  = 1 - normDir1;
        const index_t normDir2 = secondPS.direction();
        const index_t ambDir2  = 1 - normDir2;
        
        // Get degree along tangential direction for each patch
        const short_t deg1 = firstGeo.basis().degree(ambDir1);
        const short_t deg2 = secondGeo.basis().degree(ambDir2);
        const short_t maxDeg = std::max(deg1, deg2);
        
        // Determine Gauss nodes per span
        // Integrand degree: derivatives are deg-1, products of two are 2*(deg-1),
        // times linear basis φ is 2*(deg-1)+1. Gauss rule with n nodes integrates
        // degree 2n-1 exactly, so n = deg suffices. Use deg+1 for safety.
        index_t nGauss = numGaussPerSpan > 0 ? numGaussPerSpan : maxDeg + 1;
        
        // Merge knot breaks from both patches
        std::vector<T> breaks1, breaks2;
        {
            gsMatrix<T> sup1 = firstGeo.support();
            const gsBSplineBasis<T>* bb1 =
                dynamic_cast<const gsBSplineBasis<T>*>(&firstGeo.basis().component(ambDir1));
            if (bb1)
                breaks1 = bb1->knots().breaks();
            else
            {
                breaks1.push_back(sup1(ambDir1, 0));
                breaks1.push_back(sup1(ambDir1, 1));
            }
        }
        {
            gsMatrix<T> sup2 = secondGeo.support();
            const gsBSplineBasis<T>* bb2 =
                dynamic_cast<const gsBSplineBasis<T>*>(&secondGeo.basis().component(ambDir2));
            if (bb2)
                breaks2 = bb2->knots().breaks();
            else
            {
                breaks2.push_back(sup2(ambDir2, 0));
                breaks2.push_back(sup2(ambDir2, 1));
            }
        }
        
        // Merge and sort unique breaks (both must be on [0,1] for this to make sense)
        std::set<T> breakSet(breaks1.begin(), breaks1.end());
        breakSet.insert(breaks2.begin(), breaks2.end());
        std::vector<T> mergedBreaks(breakSet.begin(), breakSet.end());
        
        gsInfo << "Knot breaks patch 1 (tangential dir " << ambDir1 << "):";
        for (auto b : breaks1) gsInfo << " " << b;
        gsInfo << "\n";
        gsInfo << "Knot breaks patch 2 (tangential dir " << ambDir2 << "):";
        for (auto b : breaks2) gsInfo << " " << b;
        gsInfo << "\n";
        gsInfo << "Merged breaks:";
        for (auto b : mergedBreaks) gsInfo << " " << b;
        gsInfo << "\n";
        gsInfo << "Gauss points per span: " << nGauss << "\n";
        
        // Create Gauss quadrature mapped to all knot spans
        gsGaussRule<T> gaussRule(nGauss);
        gsMatrix<T> gaussNodes1D;
        gsVector<T> w;  // Gauss weights
        gaussRule.mapToAll(mergedBreaks, gaussNodes1D, w);
        
        const index_t N = gaussNodes1D.cols();
        gsInfo << "Total quadrature points: " << N << "\n";
        
        // Build 2D evaluation points for each patch
        gsMatrix<T> support1 = firstGeo.support();
        gsMatrix<T> support2 = secondGeo.support();
        
        gsMatrix<T> pts1(2, N), pts2(2, N);
        for (index_t i = 0; i < N; ++i)
        {
            pts1(ambDir1,  i) = gaussNodes1D(0, i);
            pts1(normDir1, i) = support1(normDir1, firstPS.parameter());
            pts2(ambDir2,  i) = gaussNodes1D(0, i);
            pts2(normDir2, i) = support2(normDir2, secondPS.parameter());
        }
        
        // Evaluate geometry values and derivatives at Gauss points
        gsMatrix<T> values1, derivatives1, values2, derivatives2;
        firstGeo.eval_into(pts1, values1);
        firstGeo.deriv_into(pts1, derivatives1);
        secondGeo.eval_into(pts2, values2);
        secondGeo.deriv_into(pts2, derivatives2);
        
        gsInfo << "values1:\n" << values1 << "\n";
        gsInfo << "derivatives1:\n" << derivatives1 << "\n";
        gsInfo << "values2:\n" << values2 << "\n";
        gsInfo << "derivatives2:\n" << derivatives2 << "\n";

        // Check if boundary values match (possibly after flipping)
        const T tol = 1e-10;
        bool flipped = false;
        if ((values1 - values2).norm() < tol)
        {
            gsInfo << "Boundary values match directly.\n";
        }
        else if ((values1 - values2.rowwise().reverse()).norm() < tol)
        {
            gsInfo << "Boundary values match after flipping. Reversing values2 and derivatives2.\n";
            values2 = values2.rowwise().reverse().eval();
            derivatives2 = derivatives2.rowwise().reverse().eval();
            // When the interface is traversed in opposite directions, the tangential
            // derivative changes sign (dG/dt -> -dG/dt), while the normal derivative
            // keeps its sign.
            for (index_t i = 0; i < N; ++i)
            {
                derivatives2(ambDir2, i)     *= -1.0;  // flip ∂_t f_x
                derivatives2(ambDir2 + 2, i) *= -1.0;  // flip ∂_t f_y
            }
            flipped = true;
        }
        else
        {
            gsInfo << "Boundary values do NOT match. Something is wrong!\n";
            return -1;
        }
        
     

        // ================================================================
        // Compute gluing data via the algorithm in gluing_data_computation_method.txt
        //
        // Using Gauss quadrature for exact integration of polynomial products.
        //
        // The AS-G1 condition (vector form):
        //   alpha_2(t)*sign1*dG1/dn + alpha_1(t)*sign2*dG2/dn + beta(t)*dG1/dt = 0
        // ================================================================
        
        // Outward normal sign: +1 if parameter=1 (east/north), -1 if parameter=0 (west/south)
        const T sign1 = firstPS.parameter()  ? 1.0 : -1.0;
        const T sign2 = secondPS.parameter() ? 1.0 : -1.0;
        gsInfo << "Normal signs: sign1=" << sign1 << " (patch " << firstPS.patch 
               << ", " << (firstPS.parameter()?"east/north":"west/south") << "), sign2=" 
               << sign2 << " (patch " << secondPS.patch 
               << ", " << (secondPS.parameter()?"east/north":"west/south") << ")\n";
        
        auto getPartial = [](const gsMatrix<T>& derivs, index_t dir, index_t idx) -> gsVector<T> {
            // deriv_into stores: row 0 = ∂_u f_x, row 1 = ∂_v f_x, row 2 = ∂_u f_y, row 3 = ∂_v f_y
            // We want the 2D vector (∂_dir f_x, ∂_dir f_y)
            gsVector<T> v(2);
            v(0) = derivs(dir,     idx);  // ∂_dir f_x
            v(1) = derivs(dir + 2, idx);  // ∂_dir f_y
            return v;
        };
        
        // 2D determinant |a, b| = a(0)*b(1) - a(1)*b(0)
        auto det2D = [](const gsVector<T>& a, const gsVector<T>& b) -> T {
            return a(0)*b(1) - a(1)*b(0);
        };

        // --- Step 0: Compute D1, D2, D3 at Gauss quadrature points ---
        gsVector<T> D1(N), D2(N), D3(N);
        // t_vals: the parameter values at quadrature points (needed for basis functions φ₀=1-t, φ₁=t)
        gsVector<T> t_vals(N);
        
        for (index_t i = 0; i < N; ++i)
        {
            t_vals(i) = gaussNodes1D(0, i);
            
            gsVector<T> dG1dn = getPartial(derivatives1, normDir1, i);
            gsVector<T> dG1dt = getPartial(derivatives1, ambDir1,  i);
            gsVector<T> dG2dn = getPartial(derivatives2, normDir2, i);
            gsVector<T> dG2dt = getPartial(derivatives2, ambDir2,  i);
            
            D1(i) = det2D(dG1dn, dG1dt);  // |partial_1 G^(1), partial_2 G^(1)|
            D2(i) = det2D(dG2dn, dG2dt);  // |partial_1 G^(2), partial_2 G^(2)|
            D3(i) = det2D(dG1dt, dG2dt);  // |partial_2 G^(1), partial_2 G^(2)|
        }
        
        gsInfo << "\nStep 0: D1, D2, D3 at quadrature points (first/last 3):\n";
        for (index_t i = 0; i < N; ++i)
        {
            if (i < 3 || i >= N-3)
                gsInfo << "  t=" << t_vals(i) 
                       << ": D1=" << D1(i) << ", D2=" << D2(i) << ", D3=" << D3(i) << "\n";
            else if (i == 3)
                gsInfo << "  ...\n";
        }
        
        // Ensure D1 > 0 and D2 > 0 (flip orientation if needed)
        T signD1 = (D1.sum() > 0) ? 1.0 : -1.0;
        T signD2 = (D2.sum() > 0) ? 1.0 : -1.0;
        if (signD1 < 0) { D1 *= -1.0; D3 *= -1.0; }
        if (signD2 < 0) { D2 *= -1.0; D3 *= -1.0; }
        
        if (D1.minCoeff() <= 0 || D2.minCoeff() <= 0)
            gsInfo << "WARNING: D1 or D2 is non-positive -- geometry may be degenerate!\n";
        
        // --- Gauss quadrature integration: int f(t) dt ≈ sum_i w_i * f(t_i) ---
        auto gaussIntegral = [&](const gsVector<T>& f) -> T {
            return w.dot(f);
        };
        
        const T eps = 1e-8;
        
        // Gluing data coefficients (linear functions in Bernstein form on [0,1])
        // alpha_k(t) = ak0*(1-t) + ak1*t,  beta_k(t) = bk0*(1-t) + bk1*t
        T a10, a11, a20, a21;
        T b10, b11, b20, b21;
        bool success = false;
        
        // ================================================================
        // Step 1: Constant alpha_1, alpha_2
        //   Constraint: alpha_1 + alpha_2 = 1
        //   Minimize: ||alpha_1*D2 - alpha_2*D1||^2
        //   Closed form: alpha_1 = int(D1*S) / int(S^2),  S = D1+D2
        // ================================================================
        gsInfo << "\n--- Step 1: Constant alpha ---\n";
        {
            gsVector<T> S = D1 + D2;
            T intD1S = gaussIntegral((D1.array()*S.array()).matrix());
            T intSS  = gaussIntegral((S.array()*S.array()).matrix());
            T alpha1c = intD1S / intSS;
            T alpha2c = 1.0 - alpha1c;

            gsVector<T> res = alpha1c*D2 - alpha2c*D1;
            T err = gaussIntegral((res.array()*res.array()).matrix());
            
            gsInfo << "  alpha_1 = " << alpha1c << ",  alpha_2 = " << alpha2c << "\n";
            gsInfo << "  ||alpha_1*D2 - alpha_2*D1||^2 = " << err << "\n";
            
            if (err < eps)
            {
                gsInfo << "  Constant alphas accepted.\n";
                a10 = alpha1c; a11 = alpha1c;
                a20 = alpha2c; a21 = alpha2c;
                
                // ============================================================
                // Step 1b: Linear beta with beta_1 = beta_2 = beta
                //   Since beta_1 = beta_2 and alpha_1+alpha_2 = 1, the combined
                //   beta is: alpha_1*beta_2 + alpha_2*beta_1 = beta*(alpha_1+alpha_2) = beta.
                //
                //   Compute beta_target(t) from the vector AS-G1 equation:
                //     alpha_2*sign1*dG1dn + alpha_1*sign2*dG2dn + beta*dG1dt = 0
                //     => beta_target = -(alpha_2*sign1*dG1dn + alpha_1*sign2*dG2dn)·dG1dt / |dG1dt|²
                //
                //   Then fit linear beta(t) = b0*(1-t) + b1*t to beta_target via L².
                // ============================================================
                gsInfo << "\n--- Step 1b: Linear beta (beta_1 = beta_2) ---\n";
                {
                    // Compute beta_target at each sample point
                    gsVector<T> beta_target(N);
                    T maxResPerp = 0;
                    for (index_t i = 0; i < N; ++i)
                    {
                        gsVector<T> dG1dn = getPartial(derivatives1, normDir1, i);
                        gsVector<T> dG1dt = getPartial(derivatives1, ambDir1,  i);
                        gsVector<T> dG2dn = getPartial(derivatives2, normDir2, i);

                        gsVector<T> v = sign1 * alpha2c * dG1dn + sign2 * alpha1c * dG2dn;
                        T dt_sq = dG1dt.squaredNorm();
                        beta_target(i) = (dt_sq > 1e-14) ? -v.dot(dG1dt) / dt_sq : 0.0;
                        T rp = (v + beta_target(i) * dG1dt).norm();
                        maxResPerp = std::max(maxResPerp, rp);
                    }

                    gsInfo << "  Max perpendicular residual = " << maxResPerp << "\n";

                    if (maxResPerp > 1e-4)
                    {
                        gsInfo << "  Perpendicular residual too large. NOT AS-G1!\n";
                    }
                    else
                    {
                        // Fit linear beta(t) = b0*(1-t) + b1*t to beta_target via L²
                        // using Gauss quadrature weights
                        gsMatrix<T> Ab(2,2); Ab.setZero();
                        gsVector<T> rhsb(2); rhsb.setZero();
                        for (index_t i = 0; i < N; ++i)
                        {
                            T phi[2] = {1-t_vals(i), t_vals(i)};
                            T wi = w(i);
                            for (index_t j = 0; j < 2; ++j)
                            {
                                for (index_t k = 0; k < 2; ++k)
                                    Ab(j,k) += wi * phi[j]*phi[k];
                                rhsb(j) += wi * phi[j] * beta_target(i);
                            }
                        }
                        gsVector<T> bc = Ab.fullPivLu().solve(rhsb);
                        b10 = bc(0); b11 = bc(1);
                        b20 = b10;   b21 = b11;

                        // Compute fitting error using Gauss quadrature
                        gsVector<T> berr_vec(N);
                        for (index_t i = 0; i < N; ++i)
                        {
                            T bv = b10*(1-t_vals(i)) + b11*t_vals(i);
                            T ri = bv - beta_target(i);
                            berr_vec(i) = ri*ri;
                        }
                        T berr = gaussIntegral(berr_vec);
                        gsInfo << "  beta: b0=" << b10 << ", b1=" << b11 << "\n";
                        gsInfo << "  ||beta_linear - beta_target||^2 = " << berr << "\n";

                        if (berr < eps) { gsInfo << "  SUCCESS (Step 1b).\n"; success = true; }
                        else            { gsInfo << "  Beta fit error too large -> Step 2.\n"; }
                    }
                }
            }
            else
            {
                gsInfo << "  Error too large -> Step 2.\n";
            }
        }
        
        // ================================================================
        // Step 2: Linear alpha_1, alpha_2
        //   alpha_k(t) = ak0*(1-t) + ak1*t
        //   Constraint: int_0^1 (alpha_1+alpha_2) dt = 1
        //     => (a10+a20)/2 + (a11+a21)/2 = 1
        //   Minimize: ||alpha_1*D2 - alpha_2*D1||^2
        //   x = [a10, a11, a20, a21]
        //   f_j: phi_0*D2, phi_1*D2, -phi_0*D1, -phi_1*D1
        //   G_jk = int f_j*f_k,  c = [1/2,1/2,1/2,1/2]
        //   KKT system: [2G c; c^T 0] [x; lam] = [0; 1]
        // ================================================================
        if (!success)
        {
            gsInfo << "\n--- Step 2: Linear alpha ---\n";
            
            gsMatrix<T> fv(N, 4);
            for (index_t i = 0; i < N; ++i)
            {
                T p0 = 1-t_vals(i), p1 = t_vals(i);
                fv(i,0) =  p0*D2(i);
                fv(i,1) =  p1*D2(i);
                fv(i,2) = -p0*D1(i);
                fv(i,3) = -p1*D1(i);
            }
            //build G (Gram matrix) for the quadratic form ||alpha_1*D2 - alpha_2*D1||^2 = x^T G x
            gsMatrix<T> G(4,4); G.setZero();
            for (index_t j = 0; j < 4; ++j)
                for (index_t k = j; k < 4; ++k)
                {
                    gsVector<T> pr(N);
                    for (index_t i = 0; i < N; ++i) pr(i) = fv(i,j)*fv(i,k);
                    G(j,k) = gaussIntegral(pr);
                    G(k,j) = G(j,k);
                }

            // Build constraint vector c = [1/2, 1/2, 1/2, 1/2]
            gsVector<T> c(4); c << 0.5, 0.5, 0.5, 0.5;
            
            // Build KKT system
            gsMatrix<T> KKT(5,5); KKT.setZero();
            KKT.block(0,0,4,4) = 2.0*G;
            KKT.block(0,4,4,1) = c;
            KKT.block(4,0,1,4) = c.transpose();

            gsInfo << "KKT matrix:\n" << KKT << "\n";
            
            gsVector<T> rhs_kkt(5); rhs_kkt.setZero(); rhs_kkt(4) = 1.0;
            gsVector<T> sol = KKT.fullPivLu().solve(rhs_kkt); // Solve for [x; lambda]
            a10 = sol(0); a11 = sol(1); a20 = sol(2); a21 = sol(3);
            
            gsVector<T> ar(N), ar2(N);
            for (index_t i = 0; i < N; ++i)
            {
                T a1v = a10*(1-t_vals(i)) + a11*t_vals(i);
                T a2v = a20*(1-t_vals(i)) + a21*t_vals(i);
                ar(i) = a1v*D2(i) - a2v*D1(i);
                ar2(i) = ar(i)*ar(i);
            }
            T aerr = gaussIntegral(ar2);
            
            gsInfo << "  alpha_1: a10=" << a10 << ", a11=" << a11 << "\n";
            gsInfo << "  alpha_2: a20=" << a20 << ", a21=" << a21 << "\n";
            gsInfo << "  ||alpha_1*D2 - alpha_2*D1||^2 = " << aerr << "\n";
            
            if (aerr < eps)
            {
                gsInfo << "  Linear alphas accepted.\n";
                
                // ============================================================
                // Step 2b: Linear beta_1, beta_2 (vector AS-G1)
                //
                //   Compute beta_comb(t) from the vector AS-G1 equation:
                //     alpha_2(t)*sign1*dG1dn + alpha_1(t)*sign2*dG2dn + beta_comb(t)*dG1dt = 0
                //     => beta_comb(t) = -(alpha_2*sign1*dG1dn + alpha_1*sign2*dG2dn)·dG1dt / |dG1dt|²
                //
                //   Then fit linear beta_1(t), beta_2(t) such that
                //     alpha_1(t)*beta_2(t) + alpha_2(t)*beta_1(t) ≈ beta_comb(t)
                //
                //   with constraint: int_0^1 (alpha_1*beta_1 - alpha_2*beta_2) dt = 0
                //
                //   y = [b10, b11, b20, b21]
                //   basis functions phi_j(t) for combined beta:
                //     g_j(t) = alpha_1(t)*phi_j(t)   for j=2,3  (contributes to beta_2)
                //            + alpha_2(t)*phi_j(t)   for j=0,1  (contributes to beta_1)
                //
                //   So: alpha_1*beta_2 + alpha_2*beta_1
                //       = b20*alpha_1*phi_0 + b21*alpha_1*phi_1 + b10*alpha_2*phi_0 + b11*alpha_2*phi_1
                //
                //   H_jk = int g_j*g_k,  d_j = int g_j*beta_comb
                //   e_j from constraint
                //   KKT: [2H e; e^T 0] [y; mu] = [2d; 0]
                // ============================================================
                gsInfo << "\n--- Step 2b: Linear beta with constraint (vector AS-G1) ---\n";
                
                // Compute beta_comb at each sample point from vector AS-G1
                gsVector<T> beta_comb(N);
                gsVector<T> res_perp_2b(N);
                for (index_t i = 0; i < N; ++i)
                {
                    T p0 = 1-t_vals(i), p1 = t_vals(i);
                    T a1v = a10*p0 + a11*p1;
                    T a2v = a20*p0 + a21*p1;
                    
                    gsVector<T> dG1dn = getPartial(derivatives1, normDir1, i);
                    gsVector<T> dG1dt = getPartial(derivatives1, ambDir1,  i);
                    gsVector<T> dG2dn = getPartial(derivatives2, normDir2, i);
                    
                    gsVector<T> v = sign1 * a2v * dG1dn + sign2 * a1v * dG2dn;
                    T dt_sq = dG1dt.squaredNorm();
                    beta_comb(i) = (dt_sq > 1e-14) ? -v.dot(dG1dt) / dt_sq : 0.0;
                    gsVector<T> res = v + beta_comb(i) * dG1dt;
                    res_perp_2b(i) = res.norm();
                }
                
                T maxResPerp2b = res_perp_2b.maxCoeff();
                gsInfo << "  Max perpendicular residual from projection = " << maxResPerp2b << "\n";
                
                if (maxResPerp2b > 1e-4)
                {
                    gsInfo << "  Perpendicular residual too large. NOT AS-G1!\n";
                }
                else
                {
                    // Build basis functions g_j(t) for the combined beta:
                    //   combined_beta = alpha_1*beta_2 + alpha_2*beta_1
                    //                 = b10*(alpha_2*phi_0) + b11*(alpha_2*phi_1) 
                    //                 + b20*(alpha_1*phi_0) + b21*(alpha_1*phi_1)
                    // So g_0 = alpha_2*phi_0, g_1 = alpha_2*phi_1, g_2 = alpha_1*phi_0, g_3 = alpha_1*phi_1
                    gsMatrix<T> gv(N, 4);
                    for (index_t i = 0; i < N; ++i)
                    {
                        T p0 = 1-t_vals(i), p1 = t_vals(i);
                        T a1v = a10*p0 + a11*p1;
                        T a2v = a20*p0 + a21*p1;
                        gv(i,0) = a2v*p0;   // beta_1 contrib: alpha_2*(1-t)
                        gv(i,1) = a2v*p1;   // beta_1 contrib: alpha_2*t
                        gv(i,2) = a1v*p0;   // beta_2 contrib: alpha_1*(1-t)
                        gv(i,3) = a1v*p1;   // beta_2 contrib: alpha_1*t
                    }
                    
                    gsMatrix<T> H(4,4); H.setZero();
                    gsVector<T> d_vec(4); d_vec.setZero();
                    for (index_t j = 0; j < 4; ++j)
                    {
                        for (index_t k = j; k < 4; ++k)
                        {
                            gsVector<T> pr(N);
                            for (index_t i = 0; i < N; ++i) pr(i) = gv(i,j)*gv(i,k);
                            H(j,k) = gaussIntegral(pr);
                            H(k,j) = H(j,k);
                        }
                        gsVector<T> pr_r(N);
                        for (index_t i = 0; i < N; ++i) pr_r(i) = gv(i,j)*beta_comb(i);
                        d_vec(j) = gaussIntegral(pr_r);
                    }
                    
                    // Constraint: int(alpha_1*beta_1 - alpha_2*beta_2) = 0
                    //   = b10*int(a1*phi_0) + b11*int(a1*phi_1) - b20*int(a2*phi_0) - b21*int(a2*phi_1)
                    gsVector<T> e_vec(4);
                    {
                        gsVector<T> a1p0(N), a1p1(N), a2p0(N), a2p1(N);
                        for (index_t i = 0; i < N; ++i)
                        {
                            T p0 = 1-t_vals(i), p1 = t_vals(i);
                            T a1v = a10*p0 + a11*p1;
                            T a2v = a20*p0 + a21*p1;
                            a1p0(i) = a1v*p0; a1p1(i) = a1v*p1;
                            a2p0(i) = a2v*p0; a2p1(i) = a2v*p1;
                        }
                        e_vec(0) =  gaussIntegral(a1p0);
                        e_vec(1) =  gaussIntegral(a1p1);
                        e_vec(2) = -gaussIntegral(a2p0);
                        e_vec(3) = -gaussIntegral(a2p1);
                    }
                    
                    gsMatrix<T> KKT2(5,5); KKT2.setZero();
                    KKT2.block(0,0,4,4) = 2.0*H;
                    KKT2.block(0,4,4,1) = e_vec;
                    KKT2.block(4,0,1,4) = e_vec.transpose();
                    
                    gsVector<T> rhs2(5);
                    rhs2.head(4) = 2.0*d_vec;
                    rhs2(4) = 0.0;
                    
                    gsVector<T> sol2 = KKT2.fullPivLu().solve(rhs2);
                    b10 = sol2(0); b11 = sol2(1);
                    b20 = sol2(2); b21 = sol2(3);
                    
                    // Check how well the combined beta matches beta_comb
                    gsVector<T> bres(N), bres2(N);
                    for (index_t i = 0; i < N; ++i)
                    {
                        T p0 = 1-t_vals(i), p1 = t_vals(i);
                        T a1v = a10*p0 + a11*p1;
                        T a2v = a20*p0 + a21*p1;
                        T b1v = b10*p0 + b11*p1;
                        T b2v = b20*p0 + b21*p1;
                        T combined = a1v*b2v + a2v*b1v;
                        bres(i) = combined - beta_comb(i);
                        bres2(i) = bres(i)*bres(i);
                    }
                    T berr = gaussIntegral(bres2);
                    
                    gsInfo << "  beta_1: b10=" << b10 << ", b11=" << b11 << "\n";
                    gsInfo << "  beta_2: b20=" << b20 << ", b21=" << b21 << "\n";
                    gsInfo << "  ||alpha_1*beta_2 + alpha_2*beta_1 - beta_comb||^2 = " << berr << "\n";
                    
                    if (berr < eps) { gsInfo << "  SUCCESS (Step 2b).\n"; success = true; }
                    else            { gsInfo << "  Beta fit error too large. NOT AS-G1!\n"; }
                }
            }
            else
            {
                gsInfo << "  Alpha error too large. NOT AS-G1!\n";
            }
        }
        
        if (!success)
        {
            gsInfo << "\nFailed to compute gluing data. The geometry is not AS-G1.\n";
            continue;
        }
        
        // ================================================================
        // Output results
        // ================================================================
        gsInfo << "\n=== Gluing data results ===\n";
        gsInfo << "  alpha_1(t) = " << a10 << "*(1-t) + " << a11 << "*t\n";
        gsInfo << "  alpha_2(t) = " << a20 << "*(1-t) + " << a21 << "*t\n";
        gsInfo << "  beta_1(t)  = " << b10 << "*(1-t) + " << b11 << "*t\n";
        gsInfo << "  beta_2(t)  = " << b20 << "*(1-t) + " << b21 << "*t\n";
        
        // ================================================================
        // Verify the AS-G1 condition at a fine set of uniform check points
        // (independent of the Gauss quadrature used for integration).
        //
        //   alpha_2(t)*sign1*dG1/dn + alpha_1(t)*sign2*dG2/dn + beta_comb(t)*dG1/dt = 0
        //
        // where beta_comb = alpha_1*beta_2 + alpha_2*beta_1
        // ================================================================
        gsInfo << "\n=== Verification (at 21 uniform check points) ===\n";
        
        const index_t Ncheck = 21;
        T maxResPerp = 0, maxResFull = 0;
        
        // Build 2D check points for both patches
        gsMatrix<T> chkPts1(2, Ncheck), chkPts2(2, Ncheck);
        for (index_t i = 0; i < Ncheck; ++i)
        {
            T tc = static_cast<T>(i) / (Ncheck - 1.);
            chkPts1(ambDir1,  i) = support1(ambDir1, 0) * (1. - tc) + support1(ambDir1, 1) * tc;
            chkPts1(normDir1, i) = support1(normDir1, firstPS.parameter());
            // If flipped, reverse the parameter direction for patch 2
            T tc2 = flipped ? (1.0 - tc) : tc;
            chkPts2(ambDir2,  i) = support2(ambDir2, 0) * (1. - tc2) + support2(ambDir2, 1) * tc2;
            chkPts2(normDir2, i) = support2(normDir2, secondPS.parameter());
        }
        
        gsMatrix<T> chkDerivs1, chkDerivs2;
        firstGeo.deriv_into(chkPts1, chkDerivs1);
        secondGeo.deriv_into(chkPts2, chkDerivs2);
        
        // Apply the same tangential derivative sign flip as during integration
        if (flipped)
        {
            for (index_t i = 0; i < Ncheck; ++i)
            {
                chkDerivs2(ambDir2, i)     *= -1.0;  // flip ∂_t f_x
                chkDerivs2(ambDir2 + 2, i) *= -1.0;  // flip ∂_t f_y
            }
        }
        
        for (index_t i = 0; i < Ncheck; ++i)
        {
            T t = static_cast<T>(i) / (Ncheck - 1.);
            T p0 = 1-t, p1 = t;
            T a1v = a10*p0 + a11*p1;
            T a2v = a20*p0 + a21*p1;
            T b1v = b10*p0 + b11*p1;
            T b2v = b20*p0 + b21*p1;
            T beta_comb_v = a1v*b2v + a2v*b1v;
            
            gsVector<T> dG1dn = getPartial(chkDerivs1, normDir1, i);
            gsVector<T> dG1dt = getPartial(chkDerivs1, ambDir1,  i);
            gsVector<T> dG2dn = getPartial(chkDerivs2, normDir2, i);
            
            // LHS of AS-G1: alpha_2*sign1*dG1dn + alpha_1*sign2*dG2dn
            gsVector<T> v = sign1 * a2v * dG1dn + sign2 * a1v * dG2dn;
            
            // Perpendicular residual (residual after projection)
            T dt_sq = dG1dt.squaredNorm();
            T beta_proj = (dt_sq > 1e-14) ? -v.dot(dG1dt) / dt_sq : 0.0;
            gsVector<T> res_perp = v + beta_proj * dG1dt;
            T resPerp = res_perp.norm();
            maxResPerp = std::max(maxResPerp, resPerp);
            
            // Full residual using computed combined beta
            gsVector<T> res_full = v + beta_comb_v * dG1dt;
            T resFull = res_full.norm();
            maxResFull = std::max(maxResFull, resFull);
            
            gsInfo << "  t=" << t << ": a1=" << a1v << ", a2=" << a2v 
                   << ", b1=" << b1v << ", b2=" << b2v
                   << ", beta_comb=" << beta_comb_v
                   << ", beta_proj=" << beta_proj
                   << ", |res_perp|=" << resPerp
                   << ", |res_full|=" << resFull << "\n";
        }
        gsInfo << "\n  Max perpendicular residual = " << maxResPerp << "\n";
        gsInfo << "  Max full residual (with fitted beta) = " << maxResFull << "\n";
        
        if (maxResFull < 1e-10)
            gsInfo << "  => AS-G1 condition SATISFIED (exact).\n";
        else if (maxResFull < 1e-4)
            gsInfo << "  => AS-G1 condition SATISFIED.\n";
        else if (maxResPerp < 1e-4)
            gsInfo << "  => AS-G1 perpendicular OK, but beta linear fit has error = " << maxResFull << ".\n";
        else
            gsInfo << "  => AS-G1 condition NOT satisfied (residual too large).\n";
    }
    
    return 0;
}

