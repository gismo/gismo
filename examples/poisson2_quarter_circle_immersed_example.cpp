/** @file poisson2_quarter_circle_immersed_example.cpp

    @brief Poisson on an immersed quarter-circle domain inside [0,1]^2.

    Embedding box: [0,1]x[0,1]
    Physical domain: Omega = {(x,y): x>0, y>0, x^2+y^2<1}
    ./build/bin/poisson2_quarter_circle_immersed_example -r 7 -e 1 --plot 

    This file is part of the G+Smo library.
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    real_t penalty = 1e3;

    gsCmdLine cmd("Poisson immersed quarter-circle in [0,1]^2.");
    cmd.addInt("e", "degreeElevation", "Number of degree elevation steps", numElevate);
    cmd.addInt("r", "uniformRefine", "Number of uniform refinement steps (before solve)", numRefine);
    cmd.addReal("p", "penalty", "Weak Dirichlet penalty", penalty);
    cmd.addSwitch("plot", "Create a ParaView file", plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // Embedding box [0,1]^2
    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineSquare());

    gsFunctionExpr<> u_exact("sin(pi*x)*sin(pi*y)", 2);
    gsFunctionExpr<> f_rhs("2*pi^2*sin(pi*x)*sin(pi*y)", 2);

    // Quarter-circle level-set: inside when max(-x,-y,x^2+y^2-1) < 0
    gsFunctionExpr<> impl_fun("max(max(-x,-y),x^2+y^2-1)", 2);
    // gsFunctionExpr<> impl_fun("1 - x^2 - y^2", 2);


    gsBoundaryConditions<> bc;
    bc.addCondition(0, boundary::west, condition_type::weak_dirichlet, &u_exact);
    bc.setGeoMap(mp);

    gsMultiBasis<> dbasis(mp, true);
    dbasis.setDegree(dbasis.maxCwiseDegree() + numElevate);
    for (index_t i = 0; i < numRefine; ++i)
        dbasis.uniformRefine();

    gsTensorBSplineBasis<2,real_t> * tbsPtr =
        dynamic_cast<gsTensorBSplineBasis<2,real_t>*>(&dbasis.basis(0));
    GISMO_ENSURE(tbsPtr, "Basis is not a tensor B-spline basis");

    gsImplicitTrimmedDomain<2,real_t> tr_domain(impl_fun, *tbsPtr);

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1, 1);
    A.options().setInt("quRule", gsQuadrature::AlgoimRule);
    A.setIntegrationDomain(memory::make_shared_not_owned(&tr_domain));
    gsExprEvaluator<> ev(A);

    geometryMap G = A.getMap(mp);
    space u = A.getSpace(dbasis);

    auto ff = A.getCoeff(f_rhs, G);
    auto u_ex = ev.getVariable(u_exact, G);

    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver;
#else
    gsSparseSolver<>::CGDiagonal solver;
#endif

    gsInfo << "(dot1=assembled, dot2=solved, dot3=error)\n\nDoFs: ";

    u.setup(bc, dirichlet::none, 0);
    A.initSystem();
    A.computePattern(igrad(u) * igrad(u).tr());

    A.assemble(
        igrad(u, G) * igrad(u, G).tr() * meas(G),
        u * ff * meas(G)
    );

    auto g_D = A.getBdrFunction(G);
    A.assembleBdr(
        bc.get("Weak Dirichlet"),
        penalty * u * u.tr(),
        penalty * u * g_D
    );

    gsInfo << A.numDofs() << "." << std::flush;

    solver.compute(A.matrix());
    solVector = solver.solve(A.rhs());
    gsInfo << "." << std::flush;

    const real_t l2err = math::sqrt(ev.integral((u_ex - u_sol).sqNorm() * meas(G)));
    const real_t h1err = l2err +
        math::sqrt(ev.integral((igrad(u_ex) - igrad(u_sol, G)).sqNorm() * meas(G)));

    gsInfo << ".\n\nL2 error: " << std::scientific << std::setprecision(3) << l2err << "\n";
    gsInfo << "H1 error: " << std::scientific << h1err << "\n";

    if (plot)
    {
        A.setIntegrationDomain(dbasis.domain());
        ev.options().setSwitch("plot.elements", true);
        ev.options().setInt("plot.npts", 1000);
        ev.writeParaview(u_sol, G, "poisson2_quarter_circle_immersed");
    }

    return EXIT_SUCCESS;
}
