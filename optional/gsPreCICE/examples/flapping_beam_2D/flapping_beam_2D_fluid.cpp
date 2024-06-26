/// This is the fluid-structure interaction benchmark FSI2 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006.
///
///
///
///TODO: FIX THE MEMORY LEAK 
/**

GISMO_DEBUG: flapping_beam_2D_solid.cpp:147, forceControlPoints.cols(): 
608
GISMO_DEBUG: flapping_beam_2D_solid.cpp:148, basis->size(): 
612

**/
/// Author: A.Shamanskiy (2016 - ...., TU Kaiserslautern)
#include <gismo.h>
#include <gsElasticity/gsElasticityAssembler.h>
#include <gsElasticity/gsElTimeIntegrator.h>
#include <gsElasticity/gsNsAssembler.h>
#include <gsElasticity/gsNsTimeIntegrator.h>
#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsALE.h>
#include <gsElasticity/gsPartitionedFSI.h>
#include <gsElasticity/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/gsGeoUtils.h>


#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEUtils.h>
#include <gsPreCICE/gsPreCICEFunction.h>
// #include <gsPreCICE/gsPreCICEVectorFunction.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    gsInfo << "Testing the two-way fluid-structure interaction solver in 2D.\n";

    //=====================================//
                // Input //
    //=====================================//

    // flow parameters
    std::string filenameFlow = "flappingBeam_flow.xml";
    real_t viscosity = 0.001;
    real_t meanVelocity = 1.;
    real_t densityFluid = 1.0e3;
    // ALE parameters
    real_t meshPR = 0.4; // poisson ratio for ALE
    real_t meshStiff = 2.5;
    index_t ALEmethod = ale_method::TINE;
    bool oneWay = false;
    // space discretization
    index_t numUniRef = 3;
    // time integration
    real_t timeStep = 0.01;
    real_t timeSpan = 15.;
    real_t thetaFluid = 0.5;
    bool imexOrNewton = true;
    bool warmUp = false;
    // output parameters
    index_t numPlotPoints = 0.;
    index_t verbosity = 0;

    std::string precice_config("../precice_config.xml");

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the two-way fluid-structure interaction solver in 2D.");
    cmd.addReal("m","mesh","Poisson's ratio for ALE",meshPR);
    cmd.addReal("x","chi","Local stiffening degree for ALE",meshStiff);
    cmd.addSwitch("o","oneway","Run as a oneway coupled simulation: beam-to-flow",oneWay);
    cmd.addInt("a","ale","ALE mesh method: 0 - HE, 1 - IHE, 2 - LE, 3 - ILE, 4 - TINE, 5 - BHE",ALEmethod);
    cmd.addInt("r","refine","Number of uniform refinement applications",numUniRef);
    cmd.addReal("t","time","Time span, sec",timeSpan);
    cmd.addReal("s","step","Time step",timeStep);
    cmd.addSwitch("w","warmup","Use large time steps during the first 2 seconds",warmUp);
    cmd.addInt("p","points","Number of points to plot to Paraview",numPlotPoints);
    cmd.addInt("v","verbosity","Amount of info printed to the prompt: 0 - none, 1 - crucial, 2 - all",verbosity);
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    // scanning geometry
    gsMultiPatch<> geoFlow;
    gsReadFile<>(filenameFlow, geoFlow);
    // I deform only three flow patches adjacent to the FSI interface
    gsMultiPatch<> geoALE;
    for (index_t p = 0; p < 3; ++p)
        geoALE.addPatch(geoFlow.patch(p+3).clone());
    geoALE.computeTopology();

    // creating bases
    for (index_t i = 0; i < numUniRef; ++i)
    {
        geoFlow.uniformRefine();
        geoALE.uniformRefine();
    }
    gsMultiBasis<> basisPressure(geoFlow);
    // I use subgrid elements (because degree elevation is not implemented for geometries, so I cant use Taylor-Hood)
    geoALE.uniformRefine();
    geoFlow.uniformRefine();
    gsMultiBasis<> basisVelocity(geoFlow);
    gsMultiBasis<> basisALE(geoALE);

    //=============================================//
        // Define preCICE setup //
    //=============================================//

    std::string participantName = "Fluid";
    gsPreCICE<real_t> participant(participantName, precice_config);

    // Mesh provided by solid
    std::string GeometryKnotMesh                = "GeometryKnotMesh";
    std::string GeometryControlPointMesh        = "GeometryControlPointMesh";

    // Mesh provided by fluid
    std::string ForceKnotMesh                   = "ForceKnotMesh";
    std::string ForceControlPointMesh           = "ForceControlPointMesh";

    // Data provided by solid
    std::string GeometryControlPointData        = "GeometryControlPointData";

    // Data provided by fluid
    std::string ForceControlPointData           = "ForceControlPointData";

    // Setup bounding box onto the force mesh
    gsMatrix<> bbox(2,2);
    bbox << -1e300, 1e300, // X dimension limits
            -1e300, 1e300; // Y dimension limits
    bbox.transposeInPlace();
    bbox.resize(1,bbox.rows()*bbox.cols());
    participant.setMeshAccessRegion(GeometryControlPointMesh, bbox);

    // Setup force control point mesh (and force knot vectors)

    GISMO_ASSERT(((dynamic_cast<gsTensorNurbsBasis<2, real_t>*>(&basisVelocity.basis(4))->knots(0) )
                == (dynamic_cast<gsTensorNurbsBasis<2, real_t>*>(&basisVelocity.basis(3))->knots(0))), "Knot vectors don't match!!!");


    gsTensorBSplineBasis<2, real_t> forceBasis(dynamic_cast<gsTensorNurbsBasis<2, real_t>*>(&basisVelocity.basis(4))->knots(0),
                                               dynamic_cast<gsTensorNurbsBasis<2, real_t>*>(&basisVelocity.basis(5))->knots(1));

    gsMatrix<> forceKnotMatrix = knotsToMatrix(forceBasis);


    participant.addMesh(ForceKnotMesh, forceKnotMatrix);

    gsMatrix<> forceControlPoints(2,forceBasis.size());
    gsDebugVar(forceControlPoints.rows());

    gsVector<index_t> forceControlPointsIDs;


    // TODO: This communicates all of the control points but we only need to communicate the boundaries
    participant.addMesh(ForceControlPointMesh,forceControlPoints, forceControlPointsIDs);

    // Initialize preCICE (send fluid mesh to API)
    real_t precice_dt = participant.initialize();

    //Get the geometry mesh from the API
    gsVector<index_t> geometryKnotIDs;
    gsMatrix<> geometryKnots;
    participant.getMeshVertexIDsAndCoordinates(GeometryKnotMesh, geometryKnotIDs, geometryKnots);

    // Gives a full tensor product basis
    gsBasis<> * geometryKnotBasis = knotMatrixToBasis<real_t>(geometryKnots).get();

    // Receive geometry control points
    gsVector<index_t> geometryControlPointIDs;
    gsMatrix<> geometryControlPoints;
    participant.getMeshVertexIDsAndCoordinates(GeometryControlPointMesh, geometryControlPointIDs, geometryControlPoints);

    // Reconstruct the beam geometry in its own mesh
    gsMultiPatch<> beam;

    beam.addPatch(give(geometryKnotBasis->makeGeometry(geometryControlPoints.transpose())));

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // source function, rhs
    gsConstantFunction<> gZero(0.,0.,2);
    // inflow velocity profile U(y) = 1.5*U_mean*y*(H-y)/(H/2)^2; channel height H = 0.41
    gsFunctionExpr<> inflow(util::to_string(meanVelocity) + "*6*y*(0.41-y)/0.41^2",2);

    // containers for solution as IGA functions
    gsMultiPatch<> velFlow, presFlow, dispBeam, dispALE, velALE;

    dispBeam.addPatch(forceBasis.makeGeometry(forceControlPoints.transpose()));
   
    // boundary conditions: flow
    gsBoundaryConditions<> bcInfoFlow;
    bcInfoFlow.addCondition(0,boundary::west,condition_type::dirichlet,&inflow,0);
    bcInfoFlow.addCondition(0,boundary::west,condition_type::dirichlet,0,1);
    for (index_t d = 0; d < 2; ++d)
    {   // no slip conditions
        bcInfoFlow.addCondition(0,boundary::east,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(1,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(1,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(2,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(2,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(3,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(3,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(4,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(4,boundary::north,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(5,boundary::west,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(6,boundary::south,condition_type::dirichlet,0,d);
        bcInfoFlow.addCondition(6,boundary::north,condition_type::dirichlet,0,d);
    }
    // flow to beam interface: these Neumann boundary condtions contain references to flow and ALE solutions;
    // by updating them, we update the boundary load as well
    gsFsiBoundaryLoad<real_t> fSouth(geoALE,patchSide(1,boundary::north),
                                     dispALE,patchSide(1,boundary::north),
                                     velFlow,patchSide(4,boundary::north),
                                     presFlow,patchSide(4,boundary::north),
                                     viscosity,densityFluid,true);
    
    gsFsiBoundaryLoad<real_t> fEast(geoALE,patchSide(2,boundary::west),
                                     dispALE,patchSide(2,boundary::west),
                                     velFlow,patchSide(5,boundary::west),
                                     presFlow,patchSide(5,boundary::west),
                                     viscosity,densityFluid,true);

    gsFsiBoundaryLoad<real_t> fNorth(geoALE,patchSide(0,boundary::south),
                                     dispALE,patchSide(0,boundary::south),
                                     velFlow,patchSide(3,boundary::south),
                                     presFlow,patchSide(3,boundary::south),
                                     viscosity,densityFluid,true);

    // gsFsiLoad<real_t> fSouth(geoALE,dispALE,1,boundary::north,
    //                          velFlow,presFlow,4,viscosity,densityFluid,true);
    // gsFsiLoad<real_t> fEast(geoALE,dispALE,2,boundary::west,
    //                         velFlow,presFlow,5,viscosity,densityFluid,true);
    // gsFsiLoad<real_t> fNorth(geoALE,dispALE,0,boundary::south,
    //                          velFlow,presFlow,3,viscosity,densityFluid,true);

    // beam to ALE interface: ALE module contains a reference to the beam displacement field;
    // by updating the displacement field, we update the displacement of the FSI interface in ALE computations
    gsBoundaryInterface interfaceBeam2ALE;
    interfaceBeam2ALE.addInterfaceSide(0,boundary::north,0,boundary::south);
    interfaceBeam2ALE.addInterfaceSide(0,boundary::south,1,boundary::north);
    interfaceBeam2ALE.addInterfaceSide(0,boundary::east,2,boundary::west);

    // ALE to flow bdry interface: Navier-Stokes solver contatains a reference to the ALE velocity field,
    // and the FSI module use the ALE displacement to deform the flow geometry
    gsBoundaryInterface interfaceALE2Flow;
    interfaceALE2Flow.addInterfaceSide(0,boundary::south,3,boundary::south);
    interfaceALE2Flow.addInterfaceSide(1,boundary::north,4,boundary::north);
    interfaceALE2Flow.addInterfaceSide(2,boundary::west,5,boundary::west);
    interfaceALE2Flow.addPatches(0,3);
    interfaceALE2Flow.addPatches(1,4);
    interfaceALE2Flow.addPatches(2,5);

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // navier stokes solver in the current configuration
    gsNsAssembler<real_t> nsAssembler(geoFlow,basisVelocity,basisPressure,bcInfoFlow,gZero);
    nsAssembler.options().setReal("Viscosity",viscosity);
    nsAssembler.options().setReal("Density",densityFluid);
    gsMassAssembler<real_t> nsMassAssembler(geoFlow,basisVelocity,bcInfoFlow,gZero);
    nsMassAssembler.options().setReal("Density",densityFluid);
    gsNsTimeIntegrator<real_t> nsTimeSolver(nsAssembler,nsMassAssembler,&velALE,&interfaceALE2Flow);
    nsTimeSolver.options().setInt("Scheme",imexOrNewton ? time_integration::implicit_nonlinear : time_integration::implicit_linear);
    nsTimeSolver.options().setReal("Theta",thetaFluid);
    nsTimeSolver.options().setSwitch("ALE",true);
    gsInfo << "Initialized Navier-Stokes system with " << nsAssembler.numDofs() << " dofs.\n";
    // mesh deformation module
    gsALE<real_t> moduleALE(geoALE,dispBeam,interfaceBeam2ALE,ale_method::method(ALEmethod));
    moduleALE.options().setReal("LocalStiff",meshStiff);
    moduleALE.options().setReal("PoissonsRatio",meshPR);
    moduleALE.options().setSwitch("Check",false);
    gsInfo << "Initialized mesh deformation system with " << moduleALE.numDofs() << " dofs.\n";

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // isogeometric fields (geometry + solution)
    gsField<> velocityField(nsAssembler.patches(),velFlow);
    gsField<> pressureField(nsAssembler.patches(),presFlow);
    gsField<> aleField(geoALE,dispALE);

    // creating containers to plot several fields corresponding to the same geometry to one Paraview file
    std::map<std::string,const gsField<> *> fieldsFlow;
    fieldsFlow["Velocity"] = &velocityField;
    fieldsFlow["Pressure"] = &pressureField;
    std::map<std::string,const gsField<> *> fieldsALE;
    fieldsALE["ALE"] = &aleField;
    // paraview collection of time steps
    gsParaviewCollection collectionFlow("flapping_beam_2D_flow");
    gsParaviewCollection collectionALE("flapping_beam_2D_ALE");

    // gsProgressBar bar;
    // gsStopwatch totalClock, iterClock;

    //=============================================//
                   // Initial condtions //
    //=============================================//

    // we will change Dirichlet DoFs for warming up, so we save them here for later
    gsMatrix<> inflowDDoFs;
    nsAssembler.getFixedDofs(0,boundary::west,inflowDDoFs);
    nsAssembler.homogenizeFixedDofs(-1);

    // set initial velocity: zero free and fixed DoFs
    nsTimeSolver.setSolutionVector(gsMatrix<>::Zero(nsAssembler.numDofs(),1));
    nsTimeSolver.setFixedDofs(nsAssembler.allFixedDofs());

    // plotting initial condition
    nsTimeSolver.constructSolution(velFlow,presFlow);
    moduleALE.constructSolution(dispALE);

    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flapping_beam_2D_flow",collectionFlow,0,numPlotPoints);
        plotDeformation(geoALE,dispALE,"flapping_beam_2D_ALE",collectionALE,0);
    }

    //=============================================//
                   // Coupled simulation //
    //=============================================//

    real_t simTime = 0.;
    real_t numTimeStep = 0;
    real_t timeALE = 0.;
    real_t timeFlow = 0.;

    // totalClock.restart();

    gsInfo << "Running the simulation...\n";

    // Time integration loop
    // gsMultiPatch<> beamNew, beamOld, dbeam;
    gsMultiPatch<> dbeam = beam;
    // beamNew = beamOld = beam;
    // dbeam = beamNew;
    dbeam.patch(0).coefs().setZero();


    index_t timestep_checkpoint = 0;
    while (participant.isCouplingOngoing())
    {

        if (participant.requiresWritingCheckpoint())
        {
            nsTimeSolver.saveState();
            moduleALE.saveState();

            // IMPORTANT!!
            //////// do something with beamNew, beamOld
            //

            gsInfo << "Writing checkpoint \n";
        }

        if (timeFlow < 2.)
            nsAssembler.setFixedDofs(0,boundary::west,inflowDDoFs*(1-cos(EIGEN_PI*(timeFlow+precice_dt)/2.))/2);

        participant.readData(GeometryControlPointMesh,GeometryControlPointData,geometryControlPointIDs,geometryControlPoints);
        gsDebugVar(geometryControlPoints.transpose());

        // beamNew.patch(0).coefs() = beam.patch(0).coefs() + geometryControlPoints.transpose();

        // gsDebugVar(beamNew.patch(0).coefs().rows());
        // gsDebugVar(beamOld.patch(0).coefs().rows());
        dbeam.patch(0).coefs() = geometryControlPoints.transpose();
        dbeam.patch(0).coefs().setZero();
        // gsDebugVar(dbeam.patch(0).coefs().rows());
        // gsDebugVar(beamNew.patch(0).coefs().transpose());
        // gsDebugVar(beamOld.patch(0).coefs().transpose());


        // project dbeam onto dispBeam
        // Represent dbeam onto disp beam's basis
        gsQuasiInterpolate<real_t>::localIntpl(dispBeam.basis(0), dbeam.patch(0), dispBeam.patch(0).coefs());

        gsDebugVar(dispBeam.patch(0).coefs().transpose());


        /*
         * MESH DEFORMATION STEP
         */
        // Given dispBeam, move the mesh and compute its velocities
        // Store the old fluid mesh in the velocity mesh
        moduleALE.constructSolution(velALE);
        gsDebugVar(velALE.patch(0).coefs().transpose());
        // Compute the fluid mesh displacements
        moduleALE.updateMesh();

        // Update the fluid mesh displacements
        moduleALE.constructSolution(dispALE);

        gsDebugVar(dispALE.patch(0).coefs().transpose());
        gsDebugVar(dispALE.patch(1).coefs().transpose());
        gsDebugVar(dispALE.patch(2).coefs().transpose());

        // Update the fluid mesh velocities (v = (u_t - u_t-1) / dt)
        for (size_t p = 0; p < velALE.nPatches(); ++p)
            velALE.patch(p).coefs() = (dispALE.patch(p).coefs() - velALE.patch(p).coefs()) / timeStep;

        gsDebugVar(velALE.patch(0).coefs().transpose());
        gsDebugVar(velALE.patch(0).coefs().transpose());
        gsDebugVar(velALE.patch(0).coefs().transpose());


        // apply new ALE deformation to the flow domain
        for (size_t p = 0; p < nsTimeSolver.aleInterface().patches.size(); ++p)
        {
            index_t pFlow = nsTimeSolver.aleInterface().patches[p].second;
            index_t pALE = nsTimeSolver.aleInterface().patches[p].first;
            nsTimeSolver.assembler().patches().patch(pFlow).coefs() += dispALE.patch(pALE).coefs();
            nsTimeSolver.mAssembler().patches().patch(pFlow).coefs() += dispALE.patch(pALE).coefs();
        }

        /*
         * FLUID SOLVE
         */
        // set velocity boundary condition on the FSI interface; velocity comes from the ALE velocity;
        // FSI inteface info is contained in the Navier-Stokes solver
        for (size_t p = 0; p < nsTimeSolver.aleInterface().sidesA.size(); ++p)
        {
            index_t pFlow = nsTimeSolver.aleInterface().sidesB[p].patch;
            boxSide sFlow = nsTimeSolver.aleInterface().sidesB[p].side();
            index_t pALE = nsTimeSolver.aleInterface().sidesA[p].patch;
            boxSide sALE = nsTimeSolver.aleInterface().sidesA[p].side();
            nsTimeSolver.assembler().setFixedDofs(pFlow,sFlow,velALE.patch(pALE).boundary(sALE)->coefs());
        }

        nsTimeSolver.makeTimeStep(timeStep);

        /*
         *
         */

        // potentially adjust non-matching timestep sizes
        timeStep = std::min(timeStep,precice_dt);

        // Construct the velocity and pressure fields
        nsTimeSolver.constructSolution(velFlow,presFlow);
        gsDebugVar(velFlow.patch(3).coefs().transpose());
        gsDebugVar(velFlow.patch(4).coefs().transpose());
        gsDebugVar(velFlow.patch(5).coefs().transpose());
        gsWriteParaview(velFlow,"flapping_beam_2D_flow");

        gsMatrix<> coefsSouth;
        gsQuasiInterpolate<real_t>::localIntpl(*basisVelocity.basis(4).boundaryBasis(boundary::north), fSouth, coefsSouth);
        gsMatrix<> coefsEast;
        gsQuasiInterpolate<real_t>::localIntpl(*basisVelocity.basis(5).boundaryBasis(boundary::west), fEast, coefsEast);
        gsMatrix<> coefsNorth;
        gsQuasiInterpolate<real_t>::localIntpl(*basisVelocity.basis(3).boundaryBasis(boundary::south), fNorth, coefsNorth);

        // gsMatrix<> coefsSouth;
        // gsQuasiInterpolate<real_t>::localIntpl(basisVelocity.basis(4).boundaryBasis(), fSouth, coefsSouth);
        // gsMatrix<> coefsEast;
        // gsQuasiInterpolate<real_t>::localIntpl(*basisVelocity.basis(5), fEast, coefsEast);
        // gsMatrix<> coefsNorth;
        // gsQuasiInterpolate<real_t>::localIntpl(*basisVelocity.basis(3), fNorth, coefsNorth);


        // gsDebugVar(coefsSouth);

        forceControlPoints.setZero();

        gsMatrix<index_t> indexSouth = basisVelocity.basis(4).boundary(boundary::north);
        gsMatrix<index_t> indexEast = basisVelocity.basis(5).boundary(boundary::west);
        gsMatrix<index_t> indexNorth = basisVelocity.basis(3).boundary(boundary::south);

        for (index_t k =0 ; k != indexSouth.size(); ++k)
            forceControlPoints.col(indexSouth(k,0)) = coefsSouth.row(k).transpose();

        for (index_t k =0 ; k != indexEast.size(); ++k)
            forceControlPoints.col(indexEast(k,0)) = coefsEast.row(k).transpose();

        for (index_t k =0 ; k != indexNorth.size(); ++k)
            forceControlPoints.col(indexNorth(k,0)) = coefsNorth.row(k).transpose();

        gsDebugVar(forceControlPoints);


        // Write the force to the solid solver

        /*
         *
         * WRITE THE STRESS TO THE SOLID SOLVER
         *
         */
        // // Write the beam displacements to the fluid solver
        participant.writeData(ForceControlPointMesh,ForceControlPointData,forceControlPointsIDs,forceControlPoints);

        // do the coupling
        precice_dt =participant.advance(timeStep);



        if (participant.requiresReadingCheckpoint())
        {
            nsTimeSolver.recoverState();
            timeStep = timestep_checkpoint;

            gsInfo << "Read Checkpoint\n";

            // IMPORTANT!!
            //////// do something with beamNew, beamOld
        }
        else
        {
            // gsTimeIntegrator advances the time step
            // advance variables
            // beamOld = beamNew;
            timeALE += timeStep;
            timeFlow += timeStep;
            numTimeStep++;

            gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flapping_beam_2D_flow",collectionFlow,numTimeStep,1000);
            //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flapping_beam_2D_ALE",collectionALE,numTimeStep,1000);
            plotDeformation(geoALE,dispALE,"flapping_beam_2D_ALE",collectionALE,numTimeStep);
        }

    }



    // while (simTime < timeSpan)
    // {
    //     // bar.display(simTime/timeSpan);
    //     // change time step for the initial warm-up phase
    //     real_t tStep = (warmUp && simTime < 2.) ? 0.1 : timeStep;
    //     // smoothly change the inflow boundary condition
    //     if (simTime < 2.)
    //         nsAssembler.setFixedDofs(0,boundary::west,inflowDDoFs*(1-cos(EIGEN_PI*(simTime+tStep)/2.))/2);
    //     if (simTime > 7.)
    //         moduleALE.options().setSwitch("Check",true);

    //     if (!moduleFSI.makeTimeStep(tStep))
    //     {
    //         gsInfo << "Invalid ALE mapping. Terminated.\n";
    //         break;
    //     }

    //     // Iteration end
    //     simTime += tStep;
    //     timeALE += moduleFSI.timeALE();
    //     timeBeam += moduleFSI.timeEL();
    //     timeFlow += moduleFSI.timeNS();
    //     numTimeStep++;

    //     if (numPlotPoints > 0)
    //     {
    //         gsWriteParaviewMultiPhysicsTimeStep(fieldsFlow,"flapping_beam_2D_flow",collectionFlow,numTimeStep,numPlotPoints);
    //         gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flapping_beam_2D_beam",collectionBeam,numTimeStep,numPlotPoints);
    //         //gsWriteParaviewMultiPhysicsTimeStep(fieldsALE,"flapping_beam_2D_ALE",collectionALE,numTimeStep,numPlotPoints);
    //         plotDeformation(geoALE,dispALE,"flapping_beam_2D_ALE",collectionALE,numTimeStep);
    //     }
    //     writeLog(logFile,nsAssembler,velFlow,presFlow,dispBeam,geoALE,dispALE,
    //              simTime,timeALE,timeFlow,timeBeam, moduleFSI.numberIterations(),
    //              nsTimeSolver.numberIterations(),elTimeSolver.numberIterations(),
    //              moduleFSI.aitkenOmega(),moduleFSI.residualNormAbs(),moduleFSI.residualNormRel());
    // }

    // //=============================================//
    //                // Final touches //
    // //=============================================//

    // gsInfo << "Complete in: " << secToHMS(totalClock.stop())
    //        << ", ALE time: " << secToHMS(timeALE)
    //        << ", flow time: " << secToHMS(timeFlow)
    //        << ", beam time: " << secToHMS(timeBeam) << std::endl;

    // if (numPlotPoints > 0)
    // {
    //     collectionFlow.save();
    //     collectionBeam.save();
    //     collectionALE.save();
    //     gsInfo << "Open \"flapping_beam_2D_*.pvd\" in Paraview for visualization.\n";
    // }
    // logFile.close();
    // gsInfo << "Log file created in \"flapping_beam_2D.txt\".\n";

    return 0;
}
