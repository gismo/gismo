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
    gsCmdLine cmd("Testing the heat equation.");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  
    // Source function
    gsConstantFunction<> f(1,2);
    gsInfo << "Source function is: " << f << "\n";

    // Define Geometry, must be a gsMultiPatch object
    gsMultiPatch<> patches(*gsNurbsCreator<>::BSplineSquareDeg(2));
    patches.computeTopology();

    // Boundary conditions
    gsBoundaryConditions<> bcInfo;
    gsConstantFunction<> g_N(1,2); // Neumann
    gsConstantFunction<> g_D(0,2); // Dirichlet
    bcInfo.setGeoMap(patches);
    bcInfo.addCondition(0, boundary::west,  condition_type::neumann  , &g_N);
    bcInfo.addCondition(0, boundary::east,  condition_type::dirichlet, &g_D);
    bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &g_D);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g_D);
    gsInfo<<"Boundary conditions:\n"<< bcInfo <<"\n";
    
    gsMultiBasis<> bases( patches );

    // Number for h-refinement of the computational (trial/test) basis.
    int numRefine  = 2;

    // Number for p-refinement of the computational (trial/test) basis.
    int numElevate = 0;

    // Elevate and p-refine the basis to order k + numElevate
    // where k is the highest degree in the bases
    if ( numElevate > -1 )
    {
        // Find maximum degree with respect to all the variables
        int tmp = bases.maxDegree(0);
        for (short_t j = 1; j < patches.parDim(); ++j )
            if ( tmp < bases.maxDegree(j) )
                tmp = bases.maxDegree(j);

        // Elevate all degrees uniformly
        tmp += numElevate;
        bases.setDegree(tmp);
    }

    // h-refine the basis
    for (int i = 0; i < numRefine; ++i)
        bases.uniformRefine();

    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n";
    
    real_t theta = 0.5;
    real_t endTime = 0.1;
    int numSteps = 40;
    
    real_t Dt = endTime / numSteps ;

    const std::string baseName("heat_eq_solution");
    gsParaviewCollection collection(baseName);

    // Generate system matrix and load vector
    gsInfo << "Assembling mass and stiffness...\n";
   
    gsExprAssembler<> K(1,1);
    gsExprAssembler<> M(1,1);

    gsInfo<<"Active options:\n"<< K.options() <<"\n";
    gsInfo<<"Active options:\n"<< M.options() <<"\n";
    
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    K.setIntegrationElements(bases);
    M.setIntegrationElements(bases);   

    gsExprEvaluator<> evK(K);
    gsExprEvaluator<> evM(M);
    
    // Set the geometry map
    geometryMap G_K = K.getMap(patches);
    geometryMap G_M = M.getMap(patches);

    // Set the discretization space
    space u_K = K.getSpace(bases);
    space u_M = M.getSpace(bases);
    //    u_K.setInterfaceCont(0);
    //    u_M.setInterfaceCont(0);
    //    u_K.addBc( bcInfo.get("Dirichlet") );
    //    u_M.addBc( bcInfo.get("Dirichlet") );

    u_K.setup(bcInfo, dirichlet::interpolation, 0);
    u_M.setup(bcInfo, dirichlet::interpolation, 0);

    // Set the source term
    variable ff_K = K.getCoeff(f, G_K);
    variable ff_M = M.getCoeff(f, G_M);

    K.initSystem(false);
    M.initSystem(false);
    K.assemble( igrad(u_K, G_K) * igrad(u_K, G_K).tr() * meas(G_K), u_K * ff_K * meas(G_K) );
    M.assemble( u_M * u_M.tr() * meas(G_M), u_M * ff_M * meas(G_M) );

    // Enforce Neumann conditions to right-hand side
    variable g_Neumann = K.getBdrFunction();
    K.assembleRhsBc(u_K * g_Neumann.val() * nv(G_K).norm(), bcInfo.neumannSides() );

    // A Conjugate Gradient linear solver with a diagonal (Jacobi) preconditionner
    gsSparseSolver<>::CGDiagonal solver;
    gsMatrix<> Sol(M.numDofs(), 1);
    
    for ( int i = 1; i<=numSteps; ++i) // for all timesteps
    {
        // Compute the system for the timestep i (rhs is assumed constant wrt time)
        gsInfo << "Solving timestep " << i*Dt << ".\n";
        Sol = solver.compute(M.matrix() +
                             Dt*theta*K.matrix()
                             ).solve(Dt*K.rhs() +
                                     (M.matrix()-Dt*(1.0-theta)*K.matrix())*Sol);
    }

    gsInfo << "Norm of the solution" << std::endl;
    gsInfo << Sol.norm() << std::endl;

    return  EXIT_SUCCESS;
}
