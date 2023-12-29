//
// Created by Jingya Li on 21/12/2023.
//

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
// #include <gsPreCICE/gsPreCICEFunction.h>
#include <gsPreCICE/gsPreCICEVectorFunction.h>

#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElasticityAssembler.h>

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

    // Create bases
    gsMultiBasis<> bases(patches);//true: poly-splines (not NURBS)

    // Set degree
    bases.setDegree( bases.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        bases.uniformRefine();
    numRefine = 0;

    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n";

    real_t rho = 3000;
    real_t E = 4e6;
    real_t nu = 0.3;
    real_t mu = E / (2.0 * (1.0 + nu));
    real_t lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));

    // Set the interface for the precice coupling
    std::vector<patchSide> couplingInterfaces(3);
    couplingInterfaces[0] = patchSide(0,boundary::east);
    couplingInterfaces[1] = patchSide(0,boundary::north);
    couplingInterfaces[2] = patchSide(0,boundary::west);


    // Set the dimension of the points
    gsMatrix<> nodes;
    // Start iteration over elements
    gsVector<> tmp;

    // TEMPORARY: get assembler options. Should come from the gsElasticityAssembler, but the problem is that that one needs to be defined using the BCs, which need the interface...

    gsExprAssembler<> A(1,1);
    gsOptionList quadOptions = A.options();

    index_t quadSize = 0;

    for (std::vector<patchSide>::const_iterator it = couplingInterfaces.begin(); it!=couplingInterfaces.end(); it++)
    {
        // Get a domain iterator on the coupling interface
        typename gsBasis<real_t>::domainIter domIt = bases.basis(it->patch).makeDomainIterator(it->side());

        // First obtain the size of all quadrature points
        typename gsQuadRule<real_t>::uPtr QuRule; // Quadrature rule  ---->OUT
        for (; domIt->good(); domIt->next() )
        {
            QuRule = gsQuadrature::getPtr(bases.basis(it->patch), quadOptions,it->side().direction());
            quadSize+=QuRule->numNodes();
        }
    }

    // Initialize parametric coordinates
    gsMatrix<> uv1(patches.domainDim(),quadSize);
    gsMatrix<> uv2(patches.domainDim(),quadSize);
    // Initialize physical coordinates
    gsMatrix<> xy1(patches.targetDim(),quadSize);
    gsMatrix<> xy2(patches.targetDim(),quadSize);

    // Grab all quadrature points
    index_t offset = 0;

    for (std::vector<patchSide>::const_iterator it = couplingInterfaces.begin(); it!=couplingInterfaces.end(); it++)
    {
        // Get a domain iterator on the coupling interface
        typename gsBasis<real_t>::domainIter domIt = bases.basis(it->patch).makeDomainIterator(it->side());
        typename gsQuadRule<real_t>::uPtr QuRule; // Quadrature rule  ---->OUT
        for (domIt->reset(); domIt->good(); domIt->next())
        {
            QuRule = gsQuadrature::getPtr(bases.basis(it->patch), quadOptions,it->side().direction());
            // Map the Quadrature rule to the element
            QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                           nodes, tmp);
            uv.block(0,offset,patches.domainDim(),QuRule->numNodes()) = nodes;

            gsMatrix<> tmp2;
            patches.patch(it->patch).eval_into(nodes,tmp2);
            xy.block(0,offset,patches.targetDim(),QuRule->numNodes()) = patches.patch(it->patch).eval(nodes);
            offset += QuRule->numNodes();
        }
    }

    // Define precice interface
    gsPreCICE<real_t> interface1("Solid1", precice_config);
    interface1.addMesh("Solid-Mesh-Downstream",xy1);
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
    gsPreCICEVectorFunction<real_t> g_C(&interface,meshID,forceID,patches,patches.geoDim());
    // Add all BCs
    // Coupling interface
    bcInfo.addCondition(0, boundary::north,  condition_type::neumann , &g_C);
    bcInfo.addCondition(0, boundary::east,  condition_type::neumann  , &g_C);
    bcInfo.addCondition(0, boundary::west,  condition_type::neumann  , &g_C);
    // Bottom side (prescribed temp)
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g_D, 0);
    bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, &g_D, 1);
    // Assign geometry map
    bcInfo.setGeoMap(patches);

// ----------------------------------------------------------------------------------------------
    // source function, rhs
    gsConstantFunction<> g(0.,0.,2);

    // source function, rhs
    gsConstantFunction<> gravity(0.,0.,2);

    // creating mass assembler
    gsMassAssembler<real_t> massAssembler(patches,bases,bcInfo,gravity);
    massAssembler.options().setReal("Density",rho);
    massAssembler.assemble();

    // creating stiffness assembler.
    gsElasticityAssembler<real_t> assembler(patches,bases,bcInfo,g);
    assembler.options().setReal("YoungsModulus",E);
    assembler.options().setReal("PoissonsRatio",nu);
    assembler.assemble();

    gsMatrix<real_t> Minv;
    gsSparseMatrix<> M = massAssembler.matrix();
    gsSparseMatrix<> K = assembler.matrix();
    gsSparseMatrix<> K_T;

    // Time step
    real_t dt = precice_dt;

    // Project u_wall as initial condition (violates Dirichlet side on precice interface)
    // RHS of the projection
    gsMatrix<> solVector;
    solVector.setZero(assembler.numDofs(),1);

    std::vector<gsMatrix<> > fixedDofs = assembler.allFixedDofs();
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
        assembler.constructSolution(solVector,fixedDofs,solution);

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
        gsDebugVar(g_C.eval(xy));
        gsDebugVar(F.transpose());
        return F;
    };

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
        assembler.constructSolution(solVector,fixedDofs,solution);
        // write heat fluxes to interface
        gsMatrix<> result(patches.geoDim(),uv.cols());
        for (index_t k=0; k!=uv.cols(); k++)
        {
            // gsDebugVar(ev.eval(nv(G),uv.col(k)));
            result.col(k) = solution.patch(0).eval(uv.col(k));
        }
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