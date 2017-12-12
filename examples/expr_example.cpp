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

    //! [Parse command line]
    cmd.getValues(argc,argv);

    // Load input file
    gsFileData<> fd(GISMO_DATA_DIR "pde/poisson2d_bvp.xml");

    gsMultiPatch<> mp;

    // topology is OK
    //gsReadFile<>(GISMO_DATA_DIR "planar/quarter_annulus_2p.xml", mp);

    fd.getId(0, mp);
    mp.computeTopology();

    gsInfo<<"Computational domain: "<< mp << "\n";
    

    //! [Function data]
    gsFunctionExpr<> f, g_N, g_D;

    // Read source function
    fd.getId(1, f);
    gsInfo<<"Source function "<< f << "\n";

    // // Dirichlet data
    // fd.getId(0, g_D);
    // gsInfo<<"Dirichlet data:"<< g_D <<"\n\n";

    // // Neumann data
    // fd.getId(0, g_N);
    // gsInfo<<"Neumann data"<< g_N <<"\n\n";

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


    // h-refine each basis (4, one for each patch)
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
    gsExprAssembler<real_t> assembler(1,1);
    typedef gsExprAssembler<real_t>::geometryMap geometryMap;
    typedef gsExprAssembler<real_t>::variable    variable;
    typedef gsExprAssembler<real_t>::space       space;
    typedef gsExprAssembler<real_t>::element     element;
    typedef gsExprAssembler<real_t>::solution    solution;

    // Elements used for numerical integration
    assembler.setIntegrationElements(dbasis);

    // Set the geometry map
    geometryMap G = assembler.setMap(mp);

    // Set the discretization basis
    variable    u  = assembler.setSpace(dbasis, bc);

    // Set the source term
    variable    ff = assembler.setCoeff(f, G);
    
    // Initialize the system
    assembler.initSystem();
    //gsInfo<<"Active options:\n"<< assembler.options() <<"\n";

    // Compute the system matrix and right-hand side
    assembler.assembleLhsRhs( igrad(u, G) * igrad(u, G).tr() * meas(G),
                              u * ff * meas(G) );

    if ( assembler.numDofs() < 20 )
    {
        gsInfo<<"Sparse Matrix:\n"<< assembler.matrix().toDense() <<"\n";
        gsInfo<<"Rhs vector:\n"<< assembler.rhs().transpose() <<"\n";
    }
    else
        gsInfo<<"Number of degrees of freedom: "<< assembler.numDofs() <<"\n";

    /*
    // Accumulate the Neumann BC contributions to the right-hand side
    variable gg = assembler.getBdrFunction();
    assembler.assembleRhsBc( u * gg, bc.neumannSides() );
    */

    gsInfo<<"Solving...\n";
    gsSparseSolver<>::CGDiagonal solver( assembler.matrix() );
    gsMatrix<> solVector = solver.solve( assembler.rhs()    );

    solution u_sol = assembler.getSolution(u, solVector);
    
    gsExprEvaluator<real_t> ev(assembler);

    gsFunctionExpr<> ms;
    fd.getId(100, ms);
    gsInfo<<"Mf: "<< ms << "\n";
    variable u_ex = ev.setVariable(ms, G);
    
    // real_t l2norm = ev.integral( (u_ex).sqr() * meas(G) );
    // gsInfo<< "* The squared L2 norm [ u.sqr() * meas(G) ]: "<<l2norm<<"\n"; 
    l2norm = math::sqrt( ev.integral( (u_sol).sqr() * meas(G) ) );
    gsInfo<< "* The squared L2 norm [ u.sqr() * meas(G) ]: "<<l2norm<<"\n";    
    l2norm = ev.integral( (u_sol-u_ex).sqNorm() * meas(G) );
    gsInfo<< "* The squared L2 norm [ u.sqr() * meas(G) ]: "<<l2norm<<"\n";    
    gsMultiPatch<> ss;
    u_sol.extract(ss);
    gsField<> sf(mp, ss, true);

    gsNormL<2> norm2(sf);
    gsInfo<< "* The squared norm of disc. solution [ gsNorm ]: "<<norm2.compute()<<"\n";    

    gsNormL<2> norm(sf,ms, false);
    gsInfo<< "* The squared L2 norm [ gsNorm ]: "<<norm.compute()<<"\n";    


    if (plot)
    {
        // Read exact solution
        fd.getId(100, f);

        gsInfo<<"Plotting in Paraview...\n";
        gsMultiPatch<> ss;
        u_sol.extract(ss);
        gsField<> sf(mp, ss, true);
        gsWriteParaview<>(sf, "solField");
        gsField<> sfe(mp, ms, false );
        gsWriteParaview<>(sfe, "solField_ex");

/*
        variable    u_ex = assembler.setCoeff(f, G);
        ev.writeParaview( u   , "solution"   , 3000, true);
        ev.writeParaview( u_ex, "solution_ex", 3000, true);
        return system("paraview solution.pvd &");
*/
    }

    return EXIT_SUCCESS;

}// end main

