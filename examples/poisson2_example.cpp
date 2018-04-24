/** @file poisson2_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, S. Takacs
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    std::string fn("pde/poisson2d_bvp.xml");
    index_t numRefine  = 5;
    index_t numElevate = 0;
    bool rates = false;
    std::string preconder("none");
    real_t tol = 1.e-8;
    index_t maxIter = 200;
    bool plot = false;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addString("f", "file",            "XML file providing the geometry.",                               fn         );
    cmd.addInt   ("e", "degreeElevation", "Number of degree elevation steps to perform before solving "
                                          "(0: equalize degree in all directions)",                         numElevate );
    cmd.addInt   ("r", "uniformRefine",   "Number of uniform h-refinement steps to perform before solving", numRefine  );
    cmd.addSwitch(     "rates",           "Iterate over all h-refinements and provide rates",               rates      );
    cmd.addString("p", "preconder",       "Preconditioner to be used by iterative solver",                  preconder  );
    cmd.addReal  ("t", "tol",             "Tolerance of iterative solver",                                  tol        );
    cmd.addInt   ("n", "maxIter",         "Maximum number of iterations in iterative solver",               maxIter    );
    cmd.addSwitch(     "plot",            "Create a ParaView visualization file with the solution",         plot       );

    cmd.getValues(argc,argv);
    //! [Parse command line]

    if ( preconder != "none" && preconder != "j" && preconder != "gs" )
    {
        gsInfo << "Unknwon preconditioner chosen. Known are only:\n"
                    "  \"none\" ... No preconditioner, i.e., identity matrix.\n"
                    "  \"j\"    ... Jacobi preconditioner.\n"
                    "  \"gs\"   ... Symmetric Gauss Seidel preconditioner.\n"
                    "\n";
        return EXIT_FAILURE;
    }

    //! [Read input file]

    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> mp;
    fd.getId(0, mp); // id=0: Multipatch domain

    gsFunctionExpr<> f;
    fd.getId(1, f); // id=1: source function
    gsInfo<<"Source function "<< f << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // id=2: boundary conditions
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    gsOptionList Aopt;
    fd.getId(4, Aopt); // id=4: assembler options

    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

    // h-refine each basis
    if (!rates)
    {
        for (index_t r = 0; r < numRefine; ++r)
            dbasis.uniformRefine();
        numRefine = 0;
    }

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);
    A.setOptions(Aopt);

    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space u = A.getSpace(dbasis);
    u.setInterfaceCont(0);
    u.addBc( bc.get("Dirichlet") );

    // Set the source term
    variable ff = A.getCoeff(f, G);

    // Recover manufactured solution
    gsFunctionExpr<> ms;
    fd.getId(3, ms); // id=3: reference solution
    //gsInfo<<"Exact solution: "<< ms << "\n";
    variable u_ex = ev.getVariable(ms, G);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);
    //! [Problem setup]

    //! [Solver loop]
    gsVector<> l2err(numRefine+1), h1err(numRefine+1);

    gsInfo<< "(dot1=assembled, dot2=solved, dot3=got_error)\n"
        "\nDoFs: ";

    for (index_t r=0; r<=numRefine; ++r)
    {
        if (r!=0)
            dbasis.uniformRefine();

        // Initialize the system
        A.initSystem();

        gsInfo<< A.numDofs() <<std::flush;

        // Compute the system matrix and right-hand side
        A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );

        // Enforce Neumann conditions to right-hand side
        variable g_N = A.getBdrFunction();
        A.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );
        //gsInfo<<"Sparse Matrix:\n"<< A.matrix().toDense() <<"\n";
        //gsInfo<<"Rhs vector:\n"<< A.rhs().transpose() <<"\n";

        gsInfo<< "." <<std::flush;// Assemblying done

        gsLinearOperator<>::Ptr preconditioner;

        if (preconder=="j")
            preconditioner = makeJacobiOp( A.matrix() );
        else if (preconder=="gs")
            preconditioner = makeSymmetricGaussSeidelOp( A.matrix() );

        gsConjugateGradient<> solver( A.matrix(), preconditioner );
        solver.setTolerance(tol);
        solver.setMaxIterations(maxIter);
        gsMatrix<> errorHistory;
        solVector.clear();
        solver.solveDetailed( A.rhs(), solVector, errorHistory );

        const bool success = solver.error() <= solver.tolerance();

        if (success)
            gsInfo<< "." <<std::flush; // Linear solving done
        else
            gsInfo<< "!" <<std::flush; // Linear solving failed

        l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
        h1err[r]= l2err[r] +
        math::sqrt(ev.integral( ( igrad(u_ex) - grad(u_sol)*jac(G).inv() ).sqNorm() * meas(G) ));

        gsInfo<< ". " <<std::flush; // Error computations done

        if (!rates)
        {
            if (success)
                gsInfo << "\n\nCG solved the problem with ";
            else
                gsInfo << "\n\nCG did not reach the desired error goal within ";

            gsInfo << ( errorHistory.rows() - 1 ) << " iterations:\n";
            if (errorHistory.rows() < 20)
                gsInfo << errorHistory.transpose() << "\n\n";
            else
                gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

        }

    } //for loop

    //! [Solver loop]

    //! [Error and convergence rates]
    gsInfo<< "\n\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";

    if (numRefine > 0)
    {
        gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
            << ( l2err.head(numRefine).array() /
                l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

        gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
            <<( h1err.head(numRefine).array() /
                h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    }
    //! [Error and convergence rates]

    // if (save)

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", true);
        ev.writeParaview( u_sol   , G, "solution");
        //ev.writeParaview( u_ex    , G, "solution_ex");
        //ev.writeParaview( u, G, "aa");

        gsFileManager::open("solution.pvd");
    }
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main
