/** @file heatEquation_example.cpp

    @brief Solves the heat equation using time-stepping

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Moore, A. Mantzaflaris
*/

#include <gismo.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false;
    gsCmdLine cmd("Testing the heat equation.");
    cmd.addSwitch("plot", "Plot the result in ParaView.", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Source function
    gsConstantFunction<> f(0,2);
    gsInfo<<"Source function is: "<< f << "\n";

    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches(*gsNurbsCreator<>::BSplineSquareDeg(2));
    patches.computeTopology();

    // Boundary conditions
    gsBoundaryConditions<> bcInfo;
    gsConstantFunction<> g_N(1,2); // Neumann
    gsConstantFunction<> g_D(0,2); // (Dirichlet
    bcInfo.addCondition(0, boundary::west,  condition_type::neumann  , &g_N);
    bcInfo.addCondition(0, boundary::east,  condition_type::dirichlet, &g_D);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &g_D);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g_D);

    gsMultiBasis<> refine_bases( patches );
    // Number for h-refinement of the computational (trial/test) basis.
    int numRefine  = 2;

    // Number for p-refinement of the computational (trial/test) basis.
    int numElevate = 0;

    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( numElevate > -1 )
    {
        // Find maximum degree with respect to all the variables
        int tmp = refine_bases.maxDegree(0);
        for (short_t j = 1; j < patches.parDim(); ++j )
            if ( tmp < refine_bases.maxDegree(j) )
                tmp = refine_bases.maxDegree(j);

        // Elevate all degrees uniformly
        tmp += numElevate;
        refine_bases.setDegree(tmp);
    }

    // h-refine the basis
    for (int i = 0; i < numRefine; ++i)
        refine_bases.uniformRefine();

    // Determines the theta-scheme used for time integration
    // (eg. Forward/backward Euler or Crank Nicolson(theta=0.5)
    real_t theta = 0.5;

    gsPoissonPde<> pde(patches, bcInfo, f);
    // Assembler (constructs stationary matrix and right-hand side vector)
    gsPoissonAssembler<> stationary(pde, refine_bases);
    stationary.options().setInt("DirichletStrategy", dirichlet::elimination);
    stationary.options().setInt("InterfaceStrategy", iFace::glue);
    gsHeatEquation<real_t> assembler(stationary);
    assembler.setTheta(theta);
    gsInfo<<assembler.options()<<"\n";

    // A Conjugate Gradient linear solver with a diagonal (Jacobi) preconditionner
    gsSparseSolver<>::CGDiagonal solver;

    // Generate system matrix and load vector
    gsInfo<<"Assembling mass and stiffness...\n";
    assembler.assemble();

    gsMatrix<> Sol, Rhs;
    int ndof = assembler.numDofs();
    real_t endTime = 0.1;
    int numSteps = 40;
    Sol.setZero(ndof, 1); // Initial solution

    real_t Dt = endTime / numSteps ;

    const std::string baseName("heat_eq_solution");
    gsParaviewCollection collection(baseName);

    std::string fileName;

    if ( plot )
    {
        //sol = assembler.constructSolution(Sol); // same as next line
        gsField<> sol = stationary.constructSolution(Sol);
        fileName = baseName + "0";
        gsWriteParaview<>(sol, fileName, 1000, true);
        collection.addTimestep(fileName,0,"0.vts");
    }

    for ( int i = 1; i<=numSteps; ++i) // for all timesteps
    {
        // Compute the system for the timestep i (rhs is assumed constant wrt time)
        assembler.nextTimeStep(Sol, Dt);
        gsInfo<<"Solving timestep "<< i*Dt<<".\n";

        // Solve for current timestep, overwrite previous solution
        Sol = solver.compute( assembler.matrix() ).solve( assembler.rhs() );

        // Obtain current solution as an isogeometric field
        //sol = assembler.constructSolution(Sol); // same as next line
        gsField<> sol = stationary.constructSolution(Sol);

        if ( plot )
        {
            // Plot the snapshot to paraview
            fileName = baseName + util::to_string(i);
            gsWriteParaview<>(sol, fileName, 1000, true);
            collection.addTimestep(fileName,i,"0.vts");
        }
    }

    //gsInfo<< " time = "<<endTime<<"\n";

    if ( plot )
    {
        collection.save();
        gsFileManager::open("heat_eq_solution.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return  EXIT_SUCCESS;
}
