/** @file gsProjection_test.cpp

    @brief Tests for gsProjection (L2, H1, H2 norms).

    Three test steps:
      Step 1 - Polynomial exactness: if f ∈ Vₕ the projection error is ≈ ε_machine.
      Step 2 - Galerkin orthogonality: the linear-system residual is ≈ ε_machine.
      Step 3 - Convergence rates: h-refinement gives the expected EOC for each norm.

    Geometry: unit square with identity map (ensures polynomial exactness in Step 1).
    Test function: f = sin(πx)sin(πy) for Steps 2-3.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#include "gismo_unittest.h"
#include <gsAssembler/gsProjection.h>
#include <gsSolver/gsSolverUtils.h>

using namespace gismo;

// Compute L2 error, H1 seminorm, and Laplacian seminorm between f_ex and a
// projected solution given by coefs on mb.basis(0) over domain mp.
static void computeErrors(const gsMultiBasis<> & mb,
                          const gsMultiPatch<> & mp,
                          const gsFunctionSet<> & f_ex,
                          const gsMatrix<>      & coefs,
                          real_t & l2err,
                          real_t & h1semi,
                          real_t & lapl_semi)
{
    gsGeometry<>::uPtr sol = mb.basis(0).makeGeometry(coefs);
    gsExprEvaluator<> ev;
    ev.setIntegrationDomain(mb.domain());
    auto G    = ev.getMap(mp);
    auto u_ex = ev.getVariable(f_ex);
    auto u_h  = ev.getVariable(*sol);
    l2err     = math::sqrt(ev.integral((u_ex - u_h).sqNorm() * meas(G)));
    h1semi    = math::sqrt(ev.integral((igrad(u_ex, G) - igrad(u_h, G)).sqNorm() * meas(G)));
    lapl_semi = math::sqrt(ev.integral((ilapl(u_ex, G) - ilapl(u_h, G)).sqNorm() * meas(G)));
}

// Unit square with affine identity geometry map.
static gsMultiPatch<> makeSquare()
{
    gsMultiPatch<> mp;
    mp.addPatch(*gsNurbsCreator<>::BSplineSquare());
    return mp;
}

SUITE(gsProjection_test)
{

    // ===================================================================
    // Step 1: Polynomial Exactness
    // If f ∈ Vₕ then u_h = f exactly; error should be ≈ ε_machine.
    // ===================================================================
    TEST(polynomial_exactness)
    {
        const index_t degree = 3;
        gsMultiPatch<> mp = makeSquare();
        gsMultiBasis<> mb(mp);
        mb.setDegree(degree);

        gsMatrix<> coefs;
        real_t err;

        // f = x² + xy lives in Vₕ for degree ≥ 2
        gsFunctionExpr<> f_poly_0("x^2+x*y", 2);
        gsFunctionExpr<> f_poly_1("2*x+y", "x", 2);   // df/dx = 2x+y, df/dy = x
        gsFunctionExpr<> f_poly_2("2", "0", "0", 2); // d^2f/dx^2 = 2, d^2f/dx
        gsFunctionWithDerivatives<> f_poly(f_poly_0, f_poly_1, f_poly_2);

        err = gsProjection<ProjectionNorm::L2, real_t>::project(mb, mp, f_poly, coefs);
        CHECK_CLOSE(math::sqrt(err), 0.0, 1e-7);

        err = gsProjection<ProjectionNorm::H1, real_t>::project(mb, mp, f_poly, coefs);
        CHECK_CLOSE(math::sqrt(err), 0.0, 1e-7);

        // f = x³ + y³  →  Δf = 6x + 6y ∈ Vₕ for degree ≥ 3
        gsFunctionExpr<> f_poly_h2_0("x^3+y^3", 2);
        gsFunctionExpr<> f_poly_h2_1("3*x^2", "3*y^2", 2); // df/dx = 3x², df/dy = 3y²
        gsFunctionExpr<> f_poly_h2_2("6*x", "6*y", "0", 2); // d²f/dx² = 6x, d²f/dxdy = 0, d²f/dy² = 0
        gsFunctionWithDerivatives<> f_poly_h2(f_poly_h2_0, f_poly_h2_1, f_poly_h2_2);
        err = gsProjection<ProjectionNorm::H2, real_t>::project(mb, mp, f_poly_h2, coefs);
    }

    TEST(lumped_system_produces_diagonal_matrix)
    {
        const index_t degree = 2;
        gsMultiPatch<> mp = makeSquare();
        gsMultiBasis<> mb(mp);
        mb.setDegree(degree);
        mb.uniformRefine(2);

        gsFunctionExpr<> f_smooth_0("sin(pi*x)*sin(pi*y)", 2);
        gsFunctionExpr<> f_smooth_1("pi*cos(pi*x)*sin(pi*y)", "pi*sin(pi*x)*cos(pi*y)", 2); // df/dx = πcos(πx)sin(πy), df/dy = πsin(πx)cos(πy)
        gsFunctionExpr<> f_smooth_2("-pi^2*sin(pi*x)*sin(pi*y)", "-pi^2*sin(pi*x)*sin(pi*y)", "-pi^2*sin(pi*x)*sin(pi*y)", 2); // d²f/dx² = -π²sin(πx)sin(πy), d²f/dxdy = -π²sin(πx)sin(πy), d²f/dy² = -π²sin(πx)sin(πy)
        gsFunctionWithDerivatives<> f_smooth(f_smooth_0, f_smooth_1, f_smooth_2);

        gsOptionList options;
        options.addSwitch("Lumped", "Use a lumped mass matrix for the projection system", true);

        gsSparseMatrix<> M;
        gsMatrix<> b;
        gsProjection<ProjectionNorm::L2, real_t>::system(mb, mp, f_smooth, M, b, gsBoundaryConditions<>(), options);

        CHECK(M.rows() > 0);
        CHECK(M.rows() == M.cols());
        CHECK(b.rows() == M.rows());

        // Verify the lumped matrix is diagonal
        bool isDiag = true;
        for (int k = 0; k < M.outerSize(); ++k)
            for (gsSparseMatrix<>::InnerIterator it(M, k); it; ++it)
                if (it.row() != it.col() && math::abs(it.value()) > 1e-14)
                    isDiag = false;
        CHECK(isDiag);

        gsMatrix<> u_h = gsSparseSolver<real_t>::SimplicialLDLT().compute(M).solve(b);
        CHECK((M * u_h - b).norm() / b.norm() < 1e-10);
    }

    // ===================================================================
    // Step 2: Galerkin Orthogonality
    // The assembled linear system should be solved to machine precision:
    // ||M u_h - b|| / ||b|| ≈ ε_machine.
    // ===================================================================
    TEST(galerkin_orthogonality)
    {
        const index_t degree = 3;
        gsMultiPatch<> mp = makeSquare();
        gsMultiBasis<> mb(mp);
        mb.setDegree(degree);
        mb.uniformRefine(2);

        gsFunctionExpr<> f_smooth_0("sin(pi*x)*sin(pi*y)", 2);
        gsFunctionExpr<> f_smooth_1("pi*cos(pi*x)*sin(pi*y)", "pi*sin(pi*x)*cos(pi*y)", 2); // df/dx = πcos(πx)sin(πy), df/dy = πsin(πx)cos(πy)
        gsFunctionExpr<> f_smooth_2("-pi^2*sin(pi*x)*sin(pi*y)", "-pi^2*sin(pi*x)*sin(pi*y)", "-pi^2*sin(pi*x)*sin(pi*y)", 2); // d²f/dx² = -π²sin(πx)sin(πy), d²f/dxdy = -π²sin(πx)sin(πy), d²f/dy² = -π²sin(πx)sin(πy)
        gsFunctionWithDerivatives<> f_smooth(f_smooth_0, f_smooth_1, f_smooth_2);

        gsSparseMatrix<> M;
        gsMatrix<>       b, u_h;

        // L2
        gsProjection<ProjectionNorm::L2, real_t>::system(mb, mp, f_smooth, M, b);
        u_h = gsSparseSolver<real_t>::SimplicialLDLT().compute(M).solve(b);
        CHECK((M * u_h - b).norm() / b.norm() < 1e-10);

        // H1
        gsProjection<ProjectionNorm::H1, real_t>::system(mb, mp, f_smooth, M, b);
        u_h = gsSparseSolver<real_t>::SimplicialLDLT().compute(M).solve(b);
        CHECK((M * u_h - b).norm() / b.norm() < 1e-10);

        // H2
        gsProjection<ProjectionNorm::H2, real_t>::system(mb, mp, f_smooth, M, b);
        u_h = gsSparseSolver<real_t>::SimplicialLDLT().compute(M).solve(b);
        CHECK((M * u_h - b).norm() / b.norm() < 1e-10);
    }

    // ===================================================================
    // Step 3: Optimal Convergence Rates
    //
    // Degree p=3, f = sin(πx)sin(πy).  Expected EOC (asymptotic):
    //   L2-projection:  ||e||_L2 → O(h^{p+1}=4),  ||e||_H1 → O(h^p=3)
    //   H1-projection:  ||e||_H1 → O(h^p=3)
    //   H2-projection:  ||Δe||   → O(h^{p-1}=2)
    //
    // One pre-refinement skips the pre-asymptotic regime; 3 loop steps give
    // 3 data points for gsSolverUtils::convergenceRateLS.
    // ===================================================================
    TEST(convergence_L2_projection)
    {
        const index_t degree    = 3;
        const index_t maxIter   = 3;
        const real_t  tol_l2   = 3.5;  // expected p+1=4
        const real_t  tol_h1   = 2.5;  // expected p=3

        gsMultiPatch<> mp = makeSquare();
        gsFunctionExpr<> f_smooth_0("sin(pi*x)*sin(pi*y)", 2);
        gsFunctionExpr<> f_smooth_1("pi*cos(pi*x)*sin(pi*y)", "pi*sin(pi*x)*cos(pi*y)", 2);
        gsFunctionExpr<> f_smooth_2("-pi^2*sin(pi*x)*sin(pi*y)", "-pi^2*sin(pi*x)*sin(pi*y)", "-pi^2*sin(pi*x)*sin(pi*y)", 2);
        gsFunctionWithDerivatives<> f_smooth(f_smooth_0, f_smooth_1, f_smooth_2);

        std::vector<real_t> h_list, l2_list, h1_list;

        gsMultiBasis<> mb(mp);
        mb.setDegree(degree);
        mb.uniformRefine(); // skip pre-asymptotic first level

        for (index_t i = 0; i < maxIter; ++i)
        {
            mb.uniformRefine();

            gsMatrix<> coefs;
            gsProjection<ProjectionNorm::L2, real_t>::project(mb, mp, f_smooth, coefs);

            real_t l2err, h1semi, dummy;
            computeErrors(mb, mp, f_smooth, coefs, l2err, h1semi, dummy);

            h_list.push_back(math::pow((real_t)mb.basis(0).size(),
                                       -1.0 / mb.dim()));
            l2_list.push_back(l2err);
            h1_list.push_back(h1semi);
        }

        real_t eoc_l2 = gsSolverUtils<>::convergenceRateLS(l2_list, h_list);
        real_t eoc_h1 = gsSolverUtils<>::convergenceRateLS(h1_list, h_list);

        CHECK(eoc_l2 > tol_l2);
        CHECK(eoc_h1 > tol_h1);
    }

    TEST(convergence_H1_projection)
    {
        const index_t degree  = 3;
        const index_t maxIter = 3;
        const real_t  tol_h1  = 2.5;  // expected p=3

        gsMultiPatch<> mp = makeSquare();
        gsFunctionExpr<> f_smooth_0("sin(pi*x)*sin(pi*y)", 2);
        gsFunctionExpr<> f_smooth_1("pi*cos(pi*x)*sin(pi*y)", "pi*sin(pi*x)*cos(pi*y)", 2);
        gsFunctionExpr<> f_smooth_2("-pi^2*sin(pi*x)*sin(pi*y)", "-pi^2*sin(pi*x)*sin(pi*y)", "-pi^2*sin(pi*x)*sin(pi*y)", 2);
        gsFunctionWithDerivatives<> f_smooth(f_smooth_0, f_smooth_1, f_smooth_2);

        std::vector<real_t> h_list, h1_list;

        gsMultiBasis<> mb(mp);
        mb.setDegree(degree);
        mb.uniformRefine();

        for (index_t i = 0; i < maxIter; ++i)
        {
            mb.uniformRefine();

            gsMatrix<> coefs;
            gsProjection<ProjectionNorm::H1, real_t>::project(mb, mp, f_smooth, coefs);

            real_t dummy, h1semi, dummy2;
            computeErrors(mb, mp, f_smooth, coefs, dummy, h1semi, dummy2);

            h_list.push_back(math::pow((real_t)mb.basis(0).size(),
                                       -1.0 / mb.dim()));
            h1_list.push_back(h1semi);
        }

        real_t eoc_h1 = gsSolverUtils<>::convergenceRateLS(h1_list, h_list);
        CHECK(eoc_h1 > tol_h1);
    }

    TEST(convergence_H2_projection)
    {
        const index_t degree   = 3;
        const index_t maxIter  = 3;
        const real_t  tol_lapl = 1.7;  // expected p-1=2

        gsMultiPatch<> mp = makeSquare();
        gsFunctionExpr<> f_smooth_0("sin(pi*x)*sin(pi*y)", 2);
        gsFunctionExpr<> f_smooth_1("pi*cos(pi*x)*sin(pi*y)", "pi*sin(pi*x)*cos(pi*y)", 2);
        gsFunctionExpr<> f_smooth_2("-pi^2*sin(pi*x)*sin(pi*y)", "-pi^2*sin(pi*x)*sin(pi*y)", "-pi^2*sin(pi*x)*sin(pi*y)", 2);
        gsFunctionWithDerivatives<> f_smooth(f_smooth_0, f_smooth_1, f_smooth_2);

        std::vector<real_t> h_list, lapl_list;

        gsMultiBasis<> mb(mp);
        mb.setDegree(degree);
        mb.uniformRefine();

        for (index_t i = 0; i < maxIter; ++i)
        {
            mb.uniformRefine();

            gsMatrix<> coefs;
            gsProjection<ProjectionNorm::H2, real_t>::project(mb, mp, f_smooth, coefs);

            real_t dummy, dummy2, lapl_semi;
            computeErrors(mb, mp, f_smooth, coefs, dummy, dummy2, lapl_semi);

            h_list.push_back(math::pow((real_t)mb.basis(0).size(),
                                       -1.0 / mb.dim()));
            lapl_list.push_back(lapl_semi);
        }

        real_t eoc_lapl = gsSolverUtils<>::convergenceRateLS(lapl_list, h_list);
        CHECK(eoc_lapl > tol_lapl);
    }

} // SUITE
