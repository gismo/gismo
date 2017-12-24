/** @file tutorialPoisson.cpp

    @brief Tutorial on how to use G+Smo to solve the Poisson equation,
    see the \ref PoissonTutorial

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
# include <gismo.h>
# include <gsAssembler/gsExprEvaluator.h>
# include <gsAssembler/gsExprAssembler.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[]) 
{
    //! [Parse command line]
    bool plot = false;
    // Number for h-refinement of the computational (trial/test) basis.
    int numRefine  = 1;
    // Number for p-refinement of the computational (trial/test) basis.
    int numElevate = 0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    cmd.getValues(argc,argv);
    //! [Parse command line]

    // Load input file
    gsFileData<> fd(GISMO_DATA_DIR "pde/poisson2d_bvp.xml");

    gsMultiPatch<> mp;
    // topology is OK
    //gsReadFile<>(GISMO_DATA_DIR "planar/quarter_annulus_2p.xml", mp);
    fd.getId(0, mp);
    mp.computeTopology(); // topology was wrong!
//    gsInfo<<"Computational domain: "<< mp << "\n";
    gsWrite(mp, "q2patch");

    //! [Function data]
    gsFunctionExpr<> f;
    fd.getId(1, f);
    gsInfo<<"Source function "<< f << "\n";
    //! [Function data]
  
    //! [Boundary conditions]
    gsBoundaryConditions<> bc;
    fd.getId(200, bc);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";
    //! [Boundary conditions]

    //! [Refinement]
    
    gsMultiBasis<> dbasis;
    for (unsigned i = 0; i < mp.nPatches(); ++i)
        dbasis.addBasis( mp.patch(i).basis().source().clone() );
    dbasis.setTopology(mp);
    gsInfo<<"B0:\n"<< dbasis.basis(0) <<"\n";
    gsInfo<<"B1:\n"<< dbasis.basis(1) <<"\n";
    gsInfo<<"Topology:\n"<< mp.topology() <<"\n";

    // h-refine each basis
    for (int i = 0; i < numRefine; ++i)
        dbasis.uniformRefine();

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

    gsInfo<<"Discretization basis: "<< dbasis.basis(0) <<"\n";
    
    //! [Assemble]
    gsExprAssembler<real_t> A(1,1);
    typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    typedef gsExprAssembler<real_t>::variable    variable;
    //typedef gsExprAssembler<real_t>::space       space;
    typedef gsExprAssembler<real_t>::element     element;
    typedef gsExprAssembler<real_t>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);

    // Set the geometry map
    geometryMap G = A.setMap(mp);

    // Set the discretization basis
    variable    u  = A.setSpace(dbasis, bc);

    // Set the source term
    variable    ff = A.setCoeff(f, G);
    
    // Initialize the system
    A.initSystem();
    //gsInfo<<"Active options:\n"<< A.options() <<"\n";

    // Compute the system matrix and right-hand side
    A.assembleLhsRhs( igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );

    // Enforce Neumann conditions
    variable g_N = A.getBdrFunction();
    A.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );
    
    if ( A.numDofs() < 20 )
    {
        gsInfo<<"Sparse Matrix:\n"<< A.matrix().toDense() <<"\n";
        gsInfo<<"Rhs vector:\n"<< A.rhs().transpose() <<"\n";
    }
    else
        gsInfo<<"Number of degrees of freedom: "<< A.numDofs() <<"\n";

    gsInfo<<"Solving...\n";
    gsSparseSolver<>::CGDiagonal solver( A.matrix() );
    gsMatrix<> solVector = solver.solve( A.rhs()    );
    solution u_sol = A.getSolution(u, solVector); // solVector is not copied
    
    gsExprEvaluator<real_t> ev(A);
    ev.options().setInt("quB", 2);

    // Recover manufactured solution
    gsFunctionExpr<> ms;
    fd.getId(100, ms);
    gsInfo<<"Exact solution: "<< ms << "\n";
    variable u_ex = ev.setVariable(ms, G);    
        
    gsMultiPatch<> ss;
    u_sol.extract(ss);
    gsField<> sf(mp, ss, true);
    variable u_sol1 = ev.setVariable(ss); //again

    real_t l2norm = math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
    gsInfo<< "* The L2 distance: "<<l2norm << " vs "
          << gsNormL<2>(sf,ms).compute() <<"\n";

    ev.options().setInt("quB", 1);

    
//    ev.integral( u_sol ); // (( ? bug) sol
        
//    ev.integral( ( u_sol1) ); // ((!!bug) sol1

    gsInfo<< "*----------------------------------- \n";
        
    gsVector<> pt(2); pt<<0.5,0.5;
//    gsDebugVar( u_sol1.rows() );

    //gsDebugVar( u_sol1.data().flags );
    gsInfo<< "*--- a0: "<<ev.eval( (u_sol1), pt )<<"\n";
    gsInfo<< "*--- b0: "<<ev.eval( (u_sol ), pt )<<"\n";

    // !!! BUG:...................
    gsInfo<< "*--- a: "<<ev.eval( grad(u_sol1), pt )<<"\n";
    
    gsInfo<< "*--- b: "<<ev.eval( grad(u_sol ), pt )<<"\n";

    gsInfo<< "*----------------------------------- \n";

    // !!! Depends on grad_expr::rows() ?????????????
    
    gsInfo<< "* sol1: "<<ev.integral( ( u_sol1 ).sqNorm() )<<"\n";
    gsInfo<< "* grad sol1: "<<ev.integral( grad( u_sol1 ).sqNorm() )<<"\n";
    gsInfo<< "* sol: "<<ev.integral( ( u_sol ).sqNorm() )<<"\n";
    gsInfo<< "* grad sol: "<<ev.integral( grad( u_sol ).sqNorm() )<<"\n";
    
    gsInfo<< "* ZERO0: "<<ev.integral( ( u_sol1 - u_sol ).sqNorm() )<<"\n";
    
    gsInfo<< "* ZERO1: "<<ev.integral( ( grad(u_sol1) - grad(u_sol) ).sqNorm() )<<"\n";
    
    real_t h1norm = math::sqrt(
        ev.integral( ( grad(u_ex) - grad(u_sol)*jac(G).ginv() ).sqNorm() * meas(G) ) ) ;
//        ev.integral( ( grad(u_ex) - grad(u_sol1)*jac(G).ginv() ).sqNorm() * meas(G) ) );
    gsInfo<< "* The H1 distance: "<<h1norm << " vs "
          << gsSeminormH1<real_t>(sf,ms).compute() <<"\n";

    ss.patch(0).coefs().setZero();
    ss.patch(1).coefs().setZero();
    h1norm = math::sqrt(
        ev.integral( ( grad(u_ex) ).sqNorm() * meas(G) ) );
    gsInfo<< "* The H1 norm: "<<h1norm << " vs "
          << gsSeminormH1<real_t>(sf,ms).compute() <<"\n";

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.writeParaview( u_sol   , G, "solution"   , 3000, true);
        ev.writeParaview( u_ex    , G, "solution_ex", 3000, true);

        //ev.writeParaview( u, G, "aa", 3000, true); ???
        //return system("paraview solution.pvd &");
    }
    
    return EXIT_SUCCESS;

}// end main

