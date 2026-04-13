/** @file poisson2_quarter_circle_immersed_example.cpp

    @brief Poisson on an immersed quarter-circle domain inside [0,1]^2.

    Embedding box: [0,1]x[0,1]
    Physical domain: Omega = {(x,y): x>0, y>0, x^2+y^2<1}
    ./build/bin/poisson2_nitsche_immersed_example -r 7 -e 1 --plot 

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
    // gsFunctionExpr<> impl_fun("max(max(-x,-y),x^2+y^2-1)", 2);
    // gsFunctionExpr<> impl_fun("1 - x^2 - y^2", 2);
    gsFunctionExpr<> impl_fun("-0.16 + (x-0.5)^2 + (y-0.5)^2", 2);
    // Square
    gsFunctionExpr<> impl_square("x*(x-1)*y*(y-1)",2);
    // gsFunctionExpr<> normal_immersed("2*(x-0.5)","2*(y-0.5)",2);
    gsFunctionExpr<> normal_immersed_unit(
        "2*(x-0.5)/sqrt(4*(x-0.5)^2 + 4*(y-0.5)^2)",
        "2*(y-0.5)/sqrt(4*(x-0.5)^2 + 4*(y-0.5)^2)",
        2
    );


    gsBoundaryConditions<> bc;
    // bc.addCondition(0, boundary::west, condition_type::weak_dirichlet, &u_exact);
    // bc.addCondition(0, boundary::south, condition_type::weak_dirichlet, &u_exact);
    // bc.addCondition(0, boundary::none,  condition_type::weak_dirichlet, &u_exact); // immersed interface
    bc.setGeoMap(mp);

    gsMultiBasis<> dbasis(mp, true);
    dbasis.setDegree(dbasis.maxCwiseDegree() + numElevate);
    
    const index_t deg = dbasis.maxCwiseDegree();
    // const real_t gamma = (penalty < 0)? real_t(2.5) * (deg + 2) * (deg + 1): penalty;
    const real_t gamma = 1e+8;


    #ifdef GISMO_WITH_PARDISO
        gsSparseSolver<>::PardisoLDLT solver;
    #else
        gsSparseSolver<>::CGDiagonal solver;
    #endif

    gsVector<> l2err(numRefine + 1), h1err(numRefine + 1);
    gsInfo<< "Degree of the basis: " << dbasis.maxCwiseDegree() << "\n";
    gsInfo << "(dot1=assembled, dot2=solved, dot3=error)\n\nDoFs: ";

    for(int r = 0; r <= numRefine; ++r)
    {        
        dbasis.uniformRefine();

        // Determine maximum mesh size
        real_t hmax = 0;
        for (size_t p=0; p!=dbasis.nBases(); p++)
        {
            hmax = math::max(hmax, dbasis.basis(p).getMaxCellLength());
        }

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

        std::vector<patchSide> bdr_immersed(1); // vector of size 1
        bdr_immersed[0]= patchSide(0,boundary::none); // assign the first and only element

        geometryMap G = A.getMap(mp);
        space u = A.getSpace(dbasis);

        auto ff    = A.getCoeff(f_rhs, G);
        auto n_imm = A.getCoeff(normal_immersed_unit, G); // grad(phi) — normalized?
        auto u_ex = ev.getVariable(u_exact, G);
        auto impl = A.getCoeff(impl_fun,G);
        gsMatrix<> solVector;
        solution u_sol = A.getSolution(u, solVector);

        u.setup(bc, dirichlet::none, 0);
        A.initSystem();
        A.computePattern(igrad(u) * igrad(u).tr());
        // A.getCoeff()

        A.assemble(
            igrad(u, G) * igrad(u, G).tr() * meas(G),
            u * ff * meas(G)
        );

        auto g_D = A.getCoeff(u_exact, G); // prescribe solution at the immersed boundary!

        // Symmetric Nitsche terms (matrix part)
        // A.assembleBdr(
        //     bdr_immersed,
        //     - (igrad(u, G) * igrad(impl,G)) * u.tr()
        //     - u * (igrad(u, G) * igrad(impl,G)).tr()
        //     + gamma / hmax * u * u.tr() * meas(G)
        // );

        A.assembleBdr(
            bdr_immersed,
            - (igrad(u, G) * n_imm) * u.tr()
            - u * (igrad(u, G) * n_imm).tr()
            + gamma / hmax * u * u.tr() * meas(G)
        );

        // Symmetric Nitsche terms (rhs part)
        // unit outward normal: n_imm.normalized()
        // A.assembleBdr(
        //     bdr_immersed,
        //     - (igrad(u, G) * igrad(impl,G)) * g_D * meas(G)
        //     + gamma / hmax * u * g_D * meas(G)
        // );

        // if i use igrad(impl,G), I get some issues in the computation! 

        A.assembleBdr(
            bdr_immersed,
            - (igrad(u, G) * n_imm) * g_D * meas(G)
            + gamma / hmax * u * g_D * meas(G)
        );

        gsInfo << A.numDofs() << "." << std::flush;

        solver.compute(A.matrix());
        solVector = solver.solve(A.rhs());
        gsInfo << "." << std::flush;

        l2err[r] = math::sqrt(ev.integral((u_ex - u_sol).sqNorm() * meas(G)));
        h1err[r] = l2err[r] + math::sqrt(ev.integral((igrad(u_ex) - igrad(u_sol, G)).sqNorm() * meas(G)));
        gsInfo << ". " << std::flush;
    
        if (plot && r == numRefine)
        {
            ev.options().setSwitch("plot.elements", true);
            ev.options().setInt("plot.npts", 100000);
            ev.writeParaview(u_sol, G, "poisson2_immersed_circle");
        }

        gsInfo << "Penalty: " << gamma/hmax << "\n";

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
    
    return EXIT_SUCCESS;
}
