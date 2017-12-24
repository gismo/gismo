/** @file poisson_example.cpp

    @brief Tutorial on how to use G+Smo to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[]) 
{
    //! [Parse command line]
    bool plot = false;
    // Number for h-refinement of the computational (trial/test) basis.
    int numRefine  = 5;
    // Number for p-refinement of the computational (trial/test) basis.
    int numElevate = 0;

    bool last = false;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    cmd.getValues(argc,argv);
    //! [Parse command line]

    //! [Read input file]

    gsFileData<> fd("pde/poisson2d_bvp.xml");
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";
    
    gsMultiPatch<> mp;
    fd.getId(0, mp);
    //gsInfo<<"Computational domain: "<< mp << "\n";
    //gsInfo<<"Topology:\n"<< mp.topology() <<"\n";

    gsFunctionExpr<> f;
    fd.getId(1, f);
    //gsInfo<<"Source function "<< f << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(200, bc);
    //gsInfo<<"Boundary conditions:\n"<< bc <<"\n";

    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis;
    for (unsigned i = 0; i < mp.nPatches(); ++i)
        dbasis.addBasis( mp.patch(i).basis().source().clone() );
        //dbasis.addBasis( mp.patch(i).basis().clone() );
    
    dbasis.degreeElevate(1,0);
    dbasis.setTopology(mp);
    //gsInfo<<"Topology:\n"<< mp.topology() <<"\n";

    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( numElevate > -1 )
    {
        // Find maximum degree with respect to all the variables
        int max_tmp = dbasis.maxCwiseDegree();
        
        // Elevate all degrees uniformly
        max_tmp += numElevate;
        dbasis.setDegree(max_tmp);
    }

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
    gsExprAssembler<real_t> A(1,1);
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    typedef gsExprAssembler<real_t>::variable    variable;
    typedef gsExprAssembler<real_t>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<real_t> ev(A);
    
    // Set the geometry map
    geometryMap G = A.setMap(mp);

    // Set the discretization space
    variable    u  = A.setSpace(dbasis, bc);

    // Set the source term
    variable    ff = A.setCoeff(f, G);

    // Recover manufactured solution
    gsFunctionExpr<> ms;
    fd.getId(100, ms);
    //gsInfo<<"Exact solution: "<< ms << "\n";
    variable u_ex = ev.setVariable(ms, G);    

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

    gsSparseSolver<>::CGDiagonal solver;
    
    //! [Problem setup]

    //! [Solver loop]
    gsVector<> l2err(numRefine+1), h1err(numRefine+1);
    gsInfo<< "\nDoFs: ";
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();
            
        // Initialize the system
        A.initSystem();
        
        gsInfo<< A.numDofs() <<std::flush;
            
        // Compute the system matrix and right-hand side
        A.assembleLhsRhs( igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );
        
        // Enforce Neumann conditions to right-hand side
        variable g_N = A.getBdrFunction();
        A.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );
        
        gsInfo<< "." <<std::flush;
        
        //gsInfo<<"Sparse Matrix:\n"<< A.matrix().toDense() <<"\n";
        //gsInfo<<"Rhs vector:\n"<< A.rhs().transpose() <<"\n";
        //gsInfo<<"Number of degrees of freedom: "<< A.numDofs() <<"\n";
        
        solver.compute( A.matrix() );
        solVector = solver.solve(A.rhs());

        gsInfo<< "." <<std::flush;
        
        //solution u_sol = A.getSolution(u, solVector); // solVector is not copied
        
        l2err[r]= math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
        //gsInfo<< "* The L2 error: "<< l2err[r] <<"\n";
        
        h1err[r]= //l2err[r] +
            math::sqrt( ev.integral( ( grad(u_ex) - grad(u_sol)*jac(G).inv() ).sqNorm() * meas(G) ) );
        //gsInfo<< "* The H1 error: "<< h1err[r] <<"\n";

        gsInfo<< ". " <<std::flush;        

    } //for loop

    //! [Solver loop]

    //! [Error and convergence rates]
    gsInfo<< "\n\nL2 error: " <<std::scientific<< l2err.transpose() <<"\n";
    gsInfo<< "H1 error: " <<std::scientific<< h1err.transpose() <<"\n";

    if (!last)
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
    {
        
    }

    //! [Export visualiuation in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.writeParaview( u_sol   , G, "solution"   , 3000, true);
        ev.writeParaview( u_ex    , G, "solution_ex", 3000, true);
        
        //ev.writeParaview( u, G, "aa", 3000, true); ???
        return system("paraview solution.pvd &");
    }
    //! [Export visualiuation in ParaView]
    
    return EXIT_SUCCESS;

}// end main

