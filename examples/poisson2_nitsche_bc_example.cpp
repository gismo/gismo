/** @file poisson2_nitsche_bc_example.cpp

    @brief Poisson on [0,1]^2 with weak Dirichlet imposition by symmetric Nitsche.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    ./build/bin/poisson2_nitsche_bc_example -r 7 -e 1 --plot
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    real_t  penalty    = -1; // negative => automatic

    gsCmdLine cmd("Poisson 2D with Nitsche weak Dirichlet BCs.");
    cmd.addInt("e", "degreeElevation",
               "Number of degree elevation steps", numElevate);
    cmd.addInt("r", "uniformRefine",
               "Number of uniform refinement loops", numRefine);
    cmd.addReal("p", "penalty", "Nitsche penalty parameter", penalty);
    cmd.addSwitch("plot", "Create a ParaView file", plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineSquare());

    gsFunctionExpr<> u_exact("sin(pi*x)*sin(pi*y)", 2);
    gsFunctionExpr<> f_rhs("2*pi^2*sin(pi*x)*sin(pi*y)", 2);
    gsFunctionExpr<> g_Dir("sin(pi*x)*sin(pi*y)", 2);

    gsBoundaryConditions<> bc;
    for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
        bc.addCondition(*bit, condition_type::weak_dirichlet, &g_Dir);
    bc.setGeoMap(mp);

    gsMultiBasis<> dbasis(mp, true);
    dbasis.setDegree(dbasis.maxCwiseDegree() + numElevate);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    typedef gsExprAssembler<>::element     element;

    gsExprAssembler<> A(1, 1);
    A.setIntegrationDomain(dbasis.domain());
    gsExprEvaluator<> ev(A);

    geometryMap G = A.getMap(mp);
    space u = A.getSpace(dbasis);
    element el = A.getElement();

    auto ff   = A.getCoeff(f_rhs, G);
    auto u_ex = ev.getVariable(u_exact, G);

    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver;
#else
    gsSparseSolver<>::CGDiagonal solver;
#endif

    gsVector<> l2err(numRefine + 1), h1err(numRefine + 1);
    gsInfo << "(dot1=assembled, dot2=solved, dot3=error)\n\nDoFs: ";

    for (int r = 0; r <= numRefine; ++r)
    {
        dbasis.uniformRefine();

        u.setup(bc, dirichlet::none, 0);
        A.initSystem();
        A.computePattern(igrad(u) * igrad(u).tr());

        const index_t deg = dbasis.maxCwiseDegree();
        const real_t alpha = (penalty < 0)
            ? real_t(2.5) * (deg + 2) * (deg + 1)
            : penalty;

        A.assemble(
            igrad(u, G) * igrad(u, G).tr() * meas(G),
            u * ff * meas(G)
        );

        auto g_D = A.getBdrFunction(G);

        // Symmetric Nitsche terms (matrix part)
        A.assembleBdr(
            bc.get("Weak Dirichlet"),
            - (igrad(u, G) * nv(G)) * u.tr()
            - u * (igrad(u, G) * nv(G)).tr()
            + alpha / el.area(G) * u * u.tr() * tv(G).norm()
        );

        // Symmetric Nitsche terms (rhs part)
        A.assembleBdr(
            bc.get("Weak Dirichlet"),
            - (igrad(u, G) * nv(G).normalized()) * g_D * nv(G).norm()
            + alpha / el.area(G) * u * g_D * tv(G).norm()
        );

        gsInfo << A.numDofs() << "." << std::flush;

        solver.compute(A.matrix());
        solVector = solver.solve(A.rhs());
        gsInfo << "." << std::flush;

        l2err[r] = math::sqrt(ev.integral((u_ex - u_sol).sqNorm() * meas(G)));
        h1err[r] = l2err[r] + math::sqrt(ev.integral((igrad(u_ex) - igrad(u_sol, G)).sqNorm() * meas(G)));
        gsInfo << ". " << std::flush;
    }

    gsInfo << "\n\nL2 error: " << std::scientific << std::setprecision(3)
           << l2err.transpose() << "\n";
    gsInfo << "H1 error: " << std::scientific << h1err.transpose() << "\n";

    if (numRefine > 0)
    {
        gsInfo << "\nEoC (L2): " << std::fixed << std::setprecision(2)
               << (l2err.head(numRefine).array() / l2err.tail(numRefine).array()).log().transpose() / std::log(2.0)
               << "\n";
        gsInfo << "EoC (H1): " << std::fixed << std::setprecision(2)
               << (h1err.head(numRefine).array() / h1err.tail(numRefine).array()).log().transpose() / std::log(2.0)
               << "\n";
    }

    if (plot)
    {
        ev.options().setSwitch("plot.elements", true);
        ev.options().setInt("plot.npts", 1000);
        ev.writeParaview(u_sol, G, "poisson2_nitsche_bc");
    }

    return EXIT_SUCCESS;
}
