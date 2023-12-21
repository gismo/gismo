/** @file flow-over-heated-plate.cpp

    @brief Heat equation participant for the PreCICE example "flow over heated plate"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
// #include <gsPreCICE/gsPreCICEFunction.h>
#include <gsPreCICE/gsPreCICEVectorFunction.h>

#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsThinShellAssembler.h>

#ifdef gsStructuralAnalysis_ENABLED
#include <gsStructuralAnalysis/gsTimeIntegrator.h>
#endif

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t plotmod = 1;
    index_t numRefine  = 0;
    index_t numElevate = 0;
    std::string precice_config;

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read input file]
    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");

    // Generate domain
    gsMultiPatch<> patches;
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(-0.05,0.0,0.05,1.0));

    // p-refine
    if (numElevate!=0)
        patches.degreeElevate(numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        patches.uniformRefine();
    numRefine = 0;

    // Create bases
    gsMultiBasis<> bases(patches);//true: poly-splines (not NURBS)

    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n";

    real_t rho = 3000;
    real_t E = 4e6;
    real_t nu = 0.3;
    real_t mu = E / (2.0 * (1.0 + nu));
    real_t lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

    // Set the interface for the precice coupling
    std::vector<patchSide> couplingInterfaces(1);
    couplingInterfaces[0] = patchSide(0,boundary::east);
    // couplingInterfaces[1] = patchSide(0,boundary::north);
    // couplingInterfaces[2] = patchSide(0,boundary::west);

    // TEMPORARY: get assembler options. Should come from the gsElasticityAssembler, but the problem is that that one needs to be defined using the BCs, which need the interface...

    gsExprAssembler<> A(1,1);
    gsOptionList quadOptions = A.options();

    // Get the quadrature points
    gsMatrix<> uv = gsQuadrature::getAllNodes(bases.basis(0),quadOptions,couplingInterfaces);
    gsMatrix<> xy = patches.patch(0).eval(uv);

    // Define precice interface
    gsPreCICE<real_t> interface("Solid", precice_config);
    interface.addMesh("Solid-Mesh",xy);
    real_t precice_dt = interface.initialize();

    index_t meshID = interface.getMeshID("Solid-Mesh");
    index_t dispID = interface.getDataID("Displacement",meshID);
    index_t forceID = interface.getDataID("Stress",meshID);

// ----------------------------------------------------------------------------------------------

    // Define boundary conditions
    gsBoundaryConditions<> bcInfo;
    // Dirichlet side
    gsConstantFunction<> g_D(0,patches.geoDim());
    // Coupling side
    // gsFunctionExpr<> g_C("1","0",patches.geoDim());
    gsPreCICEVectorFunction<real_t> g_C(&interface,meshID,forceID,patches,patches.geoDim(),false);
    // Add all BCs
    // Coupling interface
    // bcInfo.addCondition(0, boundary::north,  condition_type::neumann , &g_C);
    bcInfo.addCondition(0, boundary::east,  condition_type::neumann  , &g_C);
    // bcInfo.addCondition(0, boundary::west,  condition_type::neumann  , &g_C);
    // Bottom side (prescribed temp)
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g_D, 0, false, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g_D, 0, false, 1);
    // Assign geometry map
    bcInfo.setGeoMap(patches);

// ----------------------------------------------------------------------------------------------

    // source function, rhs
    gsConstantFunction<> Stress(0.,0.,2);

    gsConstantFunction<> Thickness(1.,2);
    gsConstantFunction<> YoungsModulus(E,2);
    gsConstantFunction<> PoissonsRatio(nu,2);
    gsConstantFunction<> Density(rho,2);

    gsMaterialMatrixLinear<2,real_t> materialMatrix(patches,
                                                    Thickness,
                                                    YoungsModulus,
                                                    PoissonsRatio,
                                                    Density);

    gsThinShellAssembler<2,real_t,false>  assembler(patches,
                                                    bases,
                                                    bcInfo,
                                                    Stress,
                                                    &materialMatrix);

    // assemble mass
    assembler.assembleMass();
    // assemble stiffness
    assembler.assemble();

    gsMatrix<real_t> Minv;
    gsSparseMatrix<> M = assembler.massMatrix();
    gsSparseMatrix<> K = assembler.matrix();
    gsSparseMatrix<> K_T;

    // Time step
    real_t dt = precice_dt;

    // Project u_wall as initial condition (violates Dirichlet side on precice interface)
    // RHS of the projection
    gsMatrix<> solVector;
    solVector.setZero(assembler.numDofs(),1);

    // Assemble the RHS
    gsMatrix<> F = assembler.rhs();
    gsMatrix<> F0;
    gsMatrix<> F_checkpoint, U_checkpoint, V_checkpoint, A_checkpoint;

    F_checkpoint = F0 = F;
    U_checkpoint = gsVector<real_t>::Zero(assembler.numDofs(),1);
    V_checkpoint = gsVector<real_t>::Zero(assembler.numDofs(),1);
    A_checkpoint = gsVector<real_t>::Zero(assembler.numDofs(),1);

    // Define the solution collection for Paraview
    gsParaviewCollection collection("solution");

    index_t timestep = 0;
    index_t timestep_checkpoint = 0;
    real_t time = 0;

    // Plot initial solution
    if (plot)
    {
        gsMultiPatch<> solution;
        assembler.constructSolution(solVector,solution);

        // solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
        gsField<> solField(patches,solution);
        std::string fileName = "./solution" + util::to_string(timestep);
        gsWriteParaview<>(solField, fileName, 500);
        fileName = "solution" + util::to_string(timestep) + "0";
        collection.addTimestep(fileName,time,".vts");
    }


    // Function for the Residual
    std::function<gsMatrix<real_t> (real_t) > Forcing;
    Forcing = [&F,&g_C,&xy](real_t time)
    {
        return F;
    };

    // gsTimeIntegrator<real_t> timeIntegrator;
    // if (!nonlinear)
    // timeIntegrator = gsTimeIntegrator<real_t>(M,K,Forcing,dt);
    // else
    //     timeIntegrator = gsTimeIntegrator<real_t>(M,Jacobian,Residual,dt);

    gsTimeIntegrator<real_t> timeIntegrator(M,K,Forcing,dt);
    timeIntegrator.setMethod("Newmark");
    timeIntegrator.setDisplacement(U_checkpoint);
    timeIntegrator.setVelocity(V_checkpoint);
    timeIntegrator.setAcceleration(A_checkpoint);
    timeIntegrator.constructSolution();

    // Time integration loop
    while (interface.isCouplingOngoing())
    {
        if (interface.isActionRequired(interface.actionWriteIterationCheckpoint()))
        {
            U_checkpoint = timeIntegrator.displacements();
            V_checkpoint = timeIntegrator.velocities();
            A_checkpoint = timeIntegrator.accelerations();

            gsInfo<<"Checkpoint written:\n";
            gsInfo<<"\t ||U|| = "<<timeIntegrator.displacements().norm()<<"\n";
            gsInfo<<"\t ||V|| = "<<timeIntegrator.velocities().norm()<<"\n";
            gsInfo<<"\t ||A|| = "<<timeIntegrator.accelerations().norm()<<"\n";


            timestep_checkpoint = timestep;
            interface.markActionFulfilled(interface.actionWriteIterationCheckpoint());
        }

        assembler.assemble();
        F = assembler.rhs();

        // solve gismo timestep
        gsInfo << "Solving timestep " << time << "...";
        timeIntegrator.step();
        timeIntegrator.constructSolution();
        solVector = timeIntegrator.displacements();
        gsInfo<<"Finished\n";

        // potentially adjust non-matching timestep sizes
        dt = std::min(dt,precice_dt);
        timeIntegrator.setTimeStep(dt);

        gsMultiPatch<> solution;
        assembler.constructSolution(solVector,solution);
        // write heat fluxes to interface
        gsMatrix<> result(patches.geoDim(),uv.cols());
        for (index_t k=0; k!=uv.cols(); k++)
            result.col(k) = solution.patch(0).eval(uv.col(k));
        interface.writeBlockVectorData(meshID,dispID,xy,result);

        // do the coupling
        precice_dt = interface.advance(dt);

        if (interface.isActionRequired(interface.actionReadIterationCheckpoint()))
        {
            /// Not converged. gsTimeIntegrator should NOT advance
            timeIntegrator.setTime(time);
            timeIntegrator.setTimeStep(dt);
            timeIntegrator.setDisplacement(U_checkpoint);
            timeIntegrator.setVelocity(V_checkpoint);
            timeIntegrator.setAcceleration(A_checkpoint);
            timeIntegrator.constructSolution();

            gsInfo<<"Checkpoint loaded:\n";
            gsInfo<<"\t ||U|| = "<<timeIntegrator.displacements().norm()<<"\n";
            gsInfo<<"\t ||V|| = "<<timeIntegrator.velocities().norm()<<"\n";
            gsInfo<<"\t ||A|| = "<<timeIntegrator.accelerations().norm()<<"\n";

            timestep = timestep_checkpoint;
            interface.markActionFulfilled(interface.actionReadIterationCheckpoint());
        }
        else
        {
            // gsTimeIntegrator advances
            // advance variables
            time += dt;
            timestep++;

            if (timestep % plotmod==0 && plot)
            {
                // solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
                gsField<> solField(patches,solution);
                std::string fileName = "./solution" + util::to_string(timestep);
                gsWriteParaview<>(solField, fileName, 500);
                fileName = "solution" + util::to_string(timestep) + "0";
                collection.addTimestep(fileName,time,".vts");
            }
        }


    }

    if (plot)
    {
        collection.save();
    }


    return  EXIT_SUCCESS;
}
