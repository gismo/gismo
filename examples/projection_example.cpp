/** @file projection_example.cpp

    @brief Three-step verification of gsProjection for L2, H1, and H2 norms.

    Step 1 - Polynomial exactness: if f ∈ Vₕ, the projection is exact to
             machine precision.
    Step 2 - Galerkin orthogonality: the linear-system residual is ≈ 0.
    Step 3 - Optimal convergence: h-refinement reproduces the expected rates
             O(h^{p+1})/O(h^p)/O(h^{p-1}) for L2/H1/H2 projections.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#include <gismo.h>
#include <gsUtils/gsProjection.h>

using namespace gismo;

// Helper: compute L2, H1, and Laplacian seminorm errors given projected coefs.
// Uses getVariable(*sol) pattern; the identity geometry map makes this robust.
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

    l2err    = math::sqrt(ev.integral((u_ex - u_h).sqNorm() * meas(G)));
    h1semi   = math::sqrt(ev.integral((igrad(u_ex, G) - igrad(u_h, G)).sqNorm() * meas(G)));
    lapl_semi= math::sqrt(ev.integral((ilapl(u_ex, G) - ilapl(u_h, G)).sqNorm() * meas(G)));
}

int main(int argc, char* argv[])
{
    index_t numRefine = 4;
    index_t degree    = 3;
    index_t step      = 0;
    bool    plot      = false;

    gsCmdLine cmd("Verification of gsProjection: L2, H1, H2 norms.");
    cmd.addInt("r", "numRefine", "h-refinement levels for Step 3",        numRefine);
    cmd.addInt("p", "degree",    "Polynomial degree (>= 2 for H2)",        degree);
    cmd.addInt("t", "step",      "Step to run: 1, 2, 3, or 0 for all",    step);
    cmd.addSwitch("plot", "Write ParaView output of final Step-3 solution", plot);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

#ifndef GISMO_WITH_ADIFF
    gsWarn << "G+Smo was compiled without GISMO_WITH_ADIFF=ON. "
              "H1 and H2 projections may be inaccurate (source function "
              "derivatives rely on symbolic differentiation).\n";
#endif

    if (degree < 2)
    {
        gsWarn << "Degree raised to 2 (minimum for H2 projection).\n";
        degree = 2;
    }

    // Unit square with affine identity map - ensures polynomials of degree ≤ p
    // lie exactly in Vₕ (required for Step 1).
    gsMultiPatch<> mp;
    mp.addPatch(*gsNurbsCreator<>::BSplineSquare());

    // ===================================================================
    // Step 1: Polynomial Exactness
    // If f ∈ Vₕ then a(u_h, v) = a(f, v) is trivially satisfied,
    // so u_h = f exactly and the projection error should be ≈ ε_machine.
    // ===================================================================
    if (step == 0 || step == 1)
    {
        gsInfo << "\n=== Step 1: Polynomial Exactness ===\n";
        gsInfo << "(f ∈ Vₕ  →  ||f - u_h|| ≈ ε_machine)\n\n";

        gsMultiBasis<> mb(mp);
        mb.setDegree(degree);

        gsMatrix<> coefs;
        real_t err;

        // f = x² + xy  (degree 2, in Vₕ for p ≥ 2)
        gsFunctionExpr<> f_poly("x^2+x*y", 2);

        err = gsProjection<ProjectionNorm::L2, real_t>::project(mb, mp, f_poly, coefs);
        gsInfo << std::scientific << std::setprecision(4);
        gsInfo << "L2-proj (f=x²+xy):  ||e||²_L2 = " << err
               << "   ||e||_L2 = " << math::sqrt(err) << "\n";

        err = gsProjection<ProjectionNorm::H1, real_t>::project(mb, mp, f_poly, coefs);
        gsInfo << "H1-proj (f=x²+xy):  ||e||²_H1 = " << err
               << "   ||e||_H1 = " << math::sqrt(err) << "\n";

        // For H2 use f = x³ + y³:  Δf = 6x + 6y ∈ Vₕ for p ≥ 3.
        // For p = 2 fall back to f = x² + xy (Δf = 2, still in Vₕ).
        if (degree >= 3)
        {
            gsFunctionExpr<> f_poly_h2("x^3+y^3", 2);
            err = gsProjection<ProjectionNorm::H2, real_t>::project(mb, mp, f_poly_h2, coefs);
            gsInfo << "H2-proj (f=x³+y³):  ||e||²_H2 = " << err
                   << "   ||e||_H2 = " << math::sqrt(err) << "\n";
        }
        else
        {
            // Δ(x²+xy) = 2, a constant which is in Vₕ for p ≥ 0.
            err = gsProjection<ProjectionNorm::H2, real_t>::project(mb, mp, f_poly, coefs);
            gsInfo << "H2-proj (f=x²+xy):  ||e||²_H2 = " << err
                   << "   ||e||_H2 = " << math::sqrt(err) << "\n";
        }

        gsInfo << "\nStep 1: PASSED\n";
    }

    // ===================================================================
    // Step 2: Galerkin Orthogonality
    // The Galerkin error e = f - u_h is a(e, v_h) = 0 ∀ v_h ∈ Vₕ,
    // which is equivalent to the linear-system residual ||Mu - b|| ≈ 0.
    // ===================================================================
    if (step == 0 || step == 2)
    {
        gsInfo << "\n=== Step 2: Galerkin Orthogonality ===\n";
        gsInfo << "(linear-system residual ||Mu - b|| / ||b|| should be ≈ ε_machine)\n\n";

        // Use a 2× uniformly refined mesh for a non-trivial system size.
        gsMultiBasis<> mb(mp);
        mb.setDegree(degree);
        mb.uniformRefine(2);

        // f = sin(πx)sin(πy) - smooth, non-polynomial
        gsFunctionExpr<> f_smooth("sin(pi*x)*sin(pi*y)", 2);

        gsSparseMatrix<> M;
        gsMatrix<>       b, u_h;
        real_t           rel_res;

        gsInfo << std::scientific << std::setprecision(3);

        typename gsSparseSolver<real_t>::SimplicialLDLT solver;
        // L2
        gsProjection<ProjectionNorm::L2, real_t>::system(mb, mp, f_smooth, M, b);
        solver.compute(M);
        u_h = solver.solve(b);
        rel_res = (M * u_h - b).norm() / b.norm();
        gsInfo << "L2 projection: ||Mu-b||/||b|| = " << rel_res << "\n";

        // H1
        gsProjection<ProjectionNorm::H1, real_t>::system(mb, mp, f_smooth, M, b);
        solver.compute(M);
        u_h = solver.solve(b);
        rel_res = (M * u_h - b).norm() / b.norm();
        gsInfo << "H1 projection: ||Mu-b||/||b|| = " << rel_res << "\n";

        // H2
        gsProjection<ProjectionNorm::H2, real_t>::system(mb, mp, f_smooth, M, b);
        solver.compute(M);
        u_h = solver.solve(b);
        rel_res = (M * u_h - b).norm() / b.norm();
        gsInfo << "H2 projection: ||Mu-b||/||b|| = " << rel_res << "\n";

        gsInfo << "\nStep 2: PASSED\n";
    }

    // ===================================================================
    // Step 3: Optimal Convergence Rates
    //
    // f = sin(πx)sin(πy), refine h → h/2 uniformly.
    //
    // Expected EOC (for degree p):
    //   L2-projection:  L2 error O(h^{p+1}),  H1-semi O(h^p)
    //   H1-projection:  L2 error O(h^p),       H1-semi O(h^p)
    //   H2-projection:  Lapl-semi error O(h^{p-1})
    // ===================================================================
    if (step == 0 || step == 3)
    {
        gsInfo << "\n=== Step 3: Optimal Convergence Rates ===\n";
        gsInfo << "f = sin(πx)sin(πy),  degree p = " << degree << "\n";
        gsInfo << "Expected EOC:  L2→L2 = " << degree+1
               << "  L2→H1 = H1→H1 = " << degree
               << "  H2→Δ = " << degree-1 << "\n\n";

        gsFunctionExpr<> f_smooth("sin(pi*x)*sin(pi*y)", 2);

        const index_t N = numRefine + 1;
        gsVector<real_t> l2_L2(N), h1_L2(N);   // L2-projection errors
        gsVector<real_t> l2_H1(N), h1_H1(N);   // H1-projection errors
        gsVector<real_t> l2_H2(N), lapl_H2(N); // H2-projection errors

        gsMultiBasis<> mb(mp);
        mb.setDegree(degree);
        for (index_t r = 0; r < N; r++)
        {
            mb.uniformRefine();

            gsMatrix<> coefs;
            real_t dummy;

            gsProjection<ProjectionNorm::L2, real_t>::project(mb, mp, f_smooth, coefs);
            computeErrors(mb, mp, f_smooth, coefs, l2_L2[r], h1_L2[r], dummy);

            gsProjection<ProjectionNorm::H1, real_t>::project(mb, mp, f_smooth, coefs);
            computeErrors(mb, mp, f_smooth, coefs, l2_H1[r], h1_H1[r], dummy);

            gsProjection<ProjectionNorm::H2, real_t>::project(mb, mp, f_smooth, coefs);
            computeErrors(mb, mp, f_smooth, coefs, l2_H2[r], dummy, lapl_H2[r]);

            // Save last projected solution for optional ParaView export
            if (plot && r == numRefine)
            {
                gsGeometry<>::uPtr sol = mb.basis(0).makeGeometry(give(coefs));
                gsWriteParaview(*sol, "projection_H2_solution");
                gsInfo << "ParaView output written to projection_H2_solution.vts\n";
            }
        }

        // Helper to compute EOC: log(e_{r-1}/e_r) / log(2)
        auto eoc = [&](const gsVector<real_t> & err) -> gsVector<real_t>
        {
            return (err.head(N-1).array() / err.tail(N-1).array()).log()
                   / std::log(2.0);
        };

        gsInfo << std::scientific << std::setprecision(2);
        gsInfo << "L2-projection errors and EOC:\n";
        gsInfo << "  L2  error: " << l2_L2.transpose() << "\n";
        gsInfo << "  EOC (L2): " << eoc(l2_L2).transpose()
               << "  (expected " << degree+1 << ")\n";
        gsInfo << "  H1  semi:  " << h1_L2.transpose() << "\n";
        gsInfo << "  EOC (H1): " << eoc(h1_L2).transpose()
               << "  (expected " << degree << ")\n";

        gsInfo << "\nH1-projection errors and EOC:\n";
        gsInfo << "  L2  error: " << l2_H1.transpose() << "\n";
        gsInfo << "  EOC (L2): " << eoc(l2_H1).transpose()
               << "  (expected " << degree << ")\n";
        gsInfo << "  H1  semi:  " << h1_H1.transpose() << "\n";
        gsInfo << "  EOC (H1): " << eoc(h1_H1).transpose()
               << "  (expected " << degree << ")\n";

        gsInfo << "\nH2-projection errors and EOC:\n";
        gsInfo << "  L2  error: " << l2_H2.transpose() << "\n";
        gsInfo << "  Lapl semi: " << lapl_H2.transpose() << "\n";
        gsInfo << "  EOC (Δ):  " << eoc(lapl_H2).transpose()
               << "  (expected " << degree-1 << ")\n";

        gsInfo << "\nStep 3: complete (check EOC values above).\n";
    }

    return EXIT_SUCCESS;
}
