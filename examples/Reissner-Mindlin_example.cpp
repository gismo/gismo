/** @file Reissner-Mindlin_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Yang Xia & Hugo Verhelst
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 0;
    bool last = true;
    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    gsMultiPatch<> mp;
    // Unit square
    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.addAutoBoundaries();
    mp.embed(3);
    real_t E_modulus = 1e0;
    real_t thickness = 1e0;
    real_t PoissonRatio = 0.0;

    gsVector<> tmp(3);
    tmp << 0, 0, 0;
    gsConstantFunction<> f(tmp,3);

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);

    for (index_t i=0; i!=3; ++i) // i=0: x-direction, i=1: y-direction, i=2: z-direction
    {
        // side, condition_type, function, unknown, parametric?, component
        // unknown: displacement or rotation
        // component: component of n-dimensional space, e.g. the i-th displacement

        // Displacements
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i );
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i );
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i );

        // Rotations
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 1, false, i );
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 1, false, i );
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 1, false, i );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 1, false, i );
    }
    tmp << 0,0,-1;

    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    gsMultiPatch<> mp_def = mp;

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine-1; ++r)
            dbasis.uniformRefine();
        numRefine = 0;
    }

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(2,2);
    // A.setOptions(Aopt);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis); // is compatible with new B-bar?? --> test with Poisson first
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap Gori = A.getMap(mp);
    geometryMap Gdef = A.getMap(mp_def);

    // Set the discretization space
    space u = A.getSpace(dbasis,3,0);
    space theta = A.getSpace(dbasis,3,1);

    // Set the source term
    variable ff = A.getCoeff(f, Gori);

    // Solution vector and solution variable
    gsMatrix<> solVectorU,solVectorT;
    solution u_sol = A.getSolution(u, solVectorU);
    solution theta_sol = A.getSolution(theta, solVectorT);

    //! [Problem setup]

    //! [Solver loop]
    gsSparseSolver<>::CGDiagonal solver;

    gsVector<> l2err(numRefine+1), h1err(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=got_error)\n"
        "\nDoFs: ";
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();

        u.setup(bc, dirichlet::interpolation, 0);
        theta.setup(bc, dirichlet::interpolation, 0); // will setup take the unknowns into account????

        // Initialize the system
        A.initSystem(false);

        gsInfo<< A.numDofs() <<std::flush;

        // Compute the system matrix and right-hand side
        A.assemble( (           igrad(u, Gori) * igrad(u, Gori).tr()                    ) * meas(Gori),
                    (           u * ff                    ) * meas(Gori) ); // assembly depends on u, theta, u_sol, theta_sol

        // Enforce Neumann conditions to right-hand side
        variable g_N = A.getBdrFunction();
        A.assembleRhsBc(
                        u * g_N.val() * nv(Gori).norm(),
                            bc.neumannSides() );

        gsInfo<< "." <<std::flush;// Assemblying done

        solver.compute( A.matrix() );
        gsVector<> solVector = solver.solve(A.rhs());

        gsInfo<< "." <<std::flush; // Linear solving done

        // l2err[r]= math::sqrt( ev.integral( (u_sol).sqNorm() * meas(G) ) );
        // h1err[r]= l2err[r] +
        // math::sqrt(ev.integral( ( igrad(u_ex) - grad(u_sol)*jac(G).inv() ).sqNorm() * meas(G) ));

        gsInfo<< ". " <<std::flush; // Error computations done

    } //for loop

    //! [Solver loop]

    //! [Error and convergence rates]
    // gsInfo<< "\n\nL2 error: "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    // gsInfo<< "H1 error: "<<std::scientific<<h1err.transpose()<<"\n";

    // if (!last && numRefine>0)
    // {
    //     gsInfo<< "\nEoC (L2): " << std::fixed<<std::setprecision(2)
    //           << ( l2err.head(numRefine).array() /
    //                l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

    //     gsInfo<<   "EoC (H1): "<< std::fixed<<std::setprecision(2)
    //           <<( h1err.head(numRefine).array() /
    //               h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    // }
    //! [Error and convergence rates]

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", true);
        ev.writeParaview( u_sol   , Gori, "solution");
        //ev.writeParaview( u_ex    , G, "solution_ex");
        //ev.writeParaview( u, G, "aa");

        gsFileManager::open("solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main
