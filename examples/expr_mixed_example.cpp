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
    int numElevate = 2;

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
    A.assembleRhsBc(u * unv(G) * g_N.val() * nv(G).norm(), bc.neumannSides() );
    
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
    gsInfo<<"Mf: "<< ms << "\n";
    variable u_ex = ev.setVariable(ms, G);
    
    gsMultiPatch<> ss;
    u_sol.extract(ss);
    variable u_sol1 = ev.setVariable(ss); // again as variable
    gsField<> sf(mp, ss, true);
    real_t l2norm1 = math::sqrt( ev.integral( (u_sol * u_sol ) * meas(G) ) );
    real_t l2norm11 = math::sqrt( ev.integral( (u_sol1 * u_sol1 ) * meas(G) ) );
    gsInfo<< "* The norm of approx. sol.: "<<l2norm1 <<"="<< l2norm11
          <<" vs " << gsNormL<2>(sf).compute() <<"\n";
    l2norm1 = math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) );
    gsInfo<< "* The L2 distance: "<<l2norm1 << " vs "
          << gsNormL<2>(sf,ms).compute() <<"\n";
    gsInfo<< "* Diff should be zero.: "<<
        ev.integral( u_sol-u_sol1 ) <<"\n";
    ss.patch(0).coefs().setZero();
    ss.patch(1).coefs().setZero();
    l2norm1 = math::sqrt( ev.integral( u_ex * u_ex  * meas(G) ) );
    gsInfo<< "* The norm of exact sol.: "<<l2norm1 << " vs "
          << gsNormL<2>(sf,ms).compute() <<"\n";

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.writeParaview( u_sol   , G, "solution"   , 3000, true);
        ev.writeParaview( u_ex    , G, "solution_ex", 3000, true);

        //ev.writeParaview( u, G, "aa", 3000, true); ???
        //return system("paraview solution.pvd &");
    }

// --------------------------------------------------------

    {
        
    gsExprAssembler<real_t> mf(2,2);

    // Elements used for numerical integration
    mf.setIntegrationElements(dbasis);

    // Set the geometry map
    geometryMap G = mf.setMap(mp);

    // Set the discretization bases
    gsMultiBasis<> dbasis1 = dbasis;
    dbasis1.degreeDecrease();
    variable    u  = mf.setSpace(dbasis1, 1       , 0); // deg-1, cont-1
    variable    v  = mf.setSpace(dbasis , mp.dim(), 1); // deg, cont    

    // Set the source term
    variable    ff = mf.setCoeff(f, G);

    // Initialize the system
    mf.initSystem();

    // Compute the system matrix and right-hand side
    mf.assembleLhs   ( v*v.tr() * meas(G)         );
    mf.assembleLhs   ( idiv(v,G) * u.tr() * meas(G) );
    mf.assembleLhsRhs( u * idiv(v,G).tr() * meas(G) , - u * ff * meas(G));

    // Apply natural conditions
    variable gg = mf.getBdrFunction();
    mf.assembleRhsBc(       // BC PROBLEM ????
        //v * unv(G) * gg.val() * nv(G).norm()
        v * nv(G) * gg.val()
        , bc.container("Dirichlet") );
    
    if ( mf.numDofs() < 20 )
    {
        gsInfo<<"Sparse Matrix:\n"<< mf.matrix().toDense() <<"\n";
        gsInfo<<"Rhs vector:\n"<< mf.rhs().transpose() <<"\n";
    }
    else
        gsInfo<<"Number of degrees of freedom: "<< mf.numDofs() <<"\n";
    
    gsInfo<<"Solving...\n";
    //gsDebugVar( mf.matrix().toDense().determinant() );
    gsSparseSolver<>::QR qr( mf.matrix() );
    gsMatrix<> solVector = qr.solve( mf.rhs() );
    solution u_sol = mf.getSolution(u, solVector); // solVector is not copied

    gsExprEvaluator<real_t> ev(mf);
    ev.options().setInt("quB", 2);
    variable u_ex = ev.setVariable(ms, G);
    gsInfo<< "* The L2 distance: "<<
       math::sqrt( ev.integral( (u_ex - u_sol).sqNorm() * meas(G) ) ) <<"\n";
 
    
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.writeParaview( u_sol   , G, "solution_mf", 3000, true);
        return system("paraview solution_mf.pvd &");
    }
    
    }
    
    return EXIT_SUCCESS;

}// end main

