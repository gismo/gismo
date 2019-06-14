/** @file assembly_example.cpp

    @brief Tutorial on how to use the gsAssembler class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

# include <gismo.h>

#include <gsAssembler/gsAssembler.h>      // included here for demonstration
#include <gsAssembler/gsVisitorPoisson.h>
#include <gsAssembler/gsVisitorNitsche.h>
#include <gsAssembler/gsVisitorNeumann.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;

    gsCmdLine cmd("Tutorial on assemblying a Poisson problem.");
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]


    // Grab a pre-defined grid of patches
    gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(2, 2, 0.5);
    gsInfo << "The domain is: "<< patches <<"\n";

    //! [Boundary conditions]
    gsFunctionExpr<> g("sin(pi*x*1)*sin(pi*y*2)+pi/10",
                       //"sin(pi*x*3)*sin(pi*y*4)-pi/10",
                       2);

    gsBoundaryConditions<> bcInfo;

    for (gsMultiPatch<>::const_biterator
             bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &g );
    }
    //! [Boundary conditions]


    //! [Refinement]
    // Copy basis from the multi-patch ( one per patch)
    gsMultiBasis<> splinebasis( patches );

    // Number for h-refinement of the computational (trial/test) basis.
    int numRefine  = 2;

    // Number for p-refinement of the computational (trial/test) basis.
    int numElevate = 2;

    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( numElevate > -1 )
    {
        // Find maximum degree with respect to all the variables
        int max_tmp = splinebasis.minCwiseDegree();
        // Elevate all degrees uniformly
        max_tmp += numElevate;
        splinebasis.setDegree(max_tmp);
    }

    // h-refine each basis (4, one for each patch)
    for (int i = 0; i < numRefine; ++i)
        splinebasis.uniformRefine();
    //! [Refinement]


    //! [Poisson Pde]
    gsFunctionExpr<> f("((pi*1)^2 + (pi*2)^2)*sin(pi*x*1)*sin(pi*y*2)",
                       //"((pi*3)^2 + (pi*4)^2)*sin(pi*x*3)*sin(pi*y*4)",
                       2);

    gsPoissonPde<> ppde(patches, bcInfo, f);
    //! [Poisson Pde]


    //! [Assembler]
    gsOptionList opt = gsAssembler<>::defaultOptions();
    opt.setInt("DirichletValues"  , dirichlet::l2Projection);
    opt.setInt("DirichletStrategy", dirichlet::elimination );
    opt.setInt("InterfaceStrategy", iFace    ::conforming  );
    opt.setReal("quA", 1.0);
    opt.setInt ("quB", 1  );
    opt.setReal("bdA", 2.0);
    opt.setInt ("bdB", 1  );
    gsInfo << "Assembler "<< opt;

    // Does 3 things:
    // 1. computes dirichlet dof values (if needed)
    // 2. Provides some routines that can be called for system assembly
    // 3. reconstructs the solution as a function, given a solution vector
    gsAssembler<real_t> PA;
    PA.initialize(ppde, splinebasis, opt);
    //! [Assembler]

    //! [Dof mapper]
    gsDofMapper mapper; // Gets the indices mapped from Basis --> Matrix

    splinebasis.getMapper((dirichlet::strategy)opt.getInt("DirichletStrategy"),
                          (iFace::    strategy)opt.getInt("InterfaceStrategy"),
                          bcInfo, mapper, 0);

    mapper.print();
    //! [Dof mapper]


    //! [Fill matrix manually]
    gsSparseSystem<> sys(mapper);
    sys.reserve(10, ppde.numRhs() ); // reserving enough space is crutial for performance!
    gsInfo << sys;
    PA.setSparseSystem(sys);

    // After the system is set the boundary dofs can be computed
    PA.computeDirichletDofs();


    PA.push<gsVisitorPoisson<real_t> >(); //bases: implied in sparsesystem + m_bases

    PA.push<gsVisitorNeumann<real_t> >(ppde.bc().neumannSides() );

    PA.finalize();

    // try this for small matrices
    //gsInfo<< PA.matrix().toDense() <<"\n"

    //! [Fill matrix manually]

    gsInfo << "Assembled a system (matrix and load vector) with "
           << PA.numDofs() << " dofs.\n";

    //! [Solve]
    // Initialize the conjugate gradient solver
    gsInfo << "Solving...\n";

    gsSparseSolver<>::CGDiagonal solver( PA.matrix() );
    gsMatrix<> solVector = solver.solve( PA.rhs() );

    gsInfo << "Solved the system with CG solver.\n";
    //! [Solve]

    //! [Construct solution]
    // Construct the solution as a scalar field
    gsMultiPatch<> mpsol;
    PA.constructSolution(solVector, mpsol);
    gsField<> sol( PA.patches(), mpsol);
    //! [Construct solution]

    if (plot)
    {
        //! [Plot in Paraview]
        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>(sol, "poisson2d", 1000);
        const gsField<> exact( PA.patches(), g, false );
        gsWriteParaview<>( exact, "poisson2d_exact", 1000);

        // Run paraview
        gsFileManager::open("poisson2d.pvd");
        //! [Plot in Paraview]
    }
    else
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    }

    return 0;
}
