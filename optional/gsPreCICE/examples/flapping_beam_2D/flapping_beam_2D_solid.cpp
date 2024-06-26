/// This is the fluid-structure interaction benchmark FSI2 from this paper:
/// "Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow"
/// Stefan Turek and Jaroslav Hron, <Fluid-Structure Interaction>, 2006.
///
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

    // beam parameters
    std::string filenameBeam = "flappingBeam_beam.xml";
    real_t youngsModulus = 1.4e6;
    real_t poissonsRatio = 0.4;
    real_t densitySolid = 1.0e4;
    real_t beamLoad = 0.0;
    // space discretization
    index_t numUniRef = 3;
    // time integration
    real_t timeStep = 0.01;
    real_t timeSpan = 15.;
    real_t thetaSolid = 1.;
    bool imexOrNewton = true;
    bool warmUp = false;
    // output parameters
    index_t numPlotPoints = 0.;
    index_t verbosity = 0;

    std::string precice_config("../precice_config.xml");

    // minimalistic user interface for terminal
    gsCmdLine cmd("Testing the two-way fluid-structure interaction solver in 2D.");
    cmd.addReal("l","load","Gravitation acceleration acting on the beam",beamLoad);
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
    gsMultiPatch<> geoBeam;
    gsReadFile<>(filenameBeam, geoBeam);

    // creating bases
    for (index_t i = 0; i < numUniRef; ++i)
        geoBeam.uniformRefine();
    gsMultiBasis<> basisDisplacement(geoBeam);


    //=============================================//
        // Define preCICE setup //
    //=============================================//

    std::string participantName = "Solid";
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
    participant.setMeshAccessRegion(ForceControlPointMesh, bbox);

    // Setup geometry control point mesh
    gsVector<index_t> geometryControlPointIDs;
    gsMatrix<> geometryControlPoints = geoBeam.patch(0).coefs().transpose();
    participant.addMesh(GeometryControlPointMesh, geometryControlPoints, geometryControlPointIDs);

    // Setup the geometry knot mesh
    gsMatrix<> geometryKnots = knotsToMatrix(basisDisplacement.basis(0));
    participant.addMesh(GeometryKnotMesh,geometryKnots);

    // Initialize preCICE (send solid mesh to API)
    real_t precice_dt = participant.initialize();


    //Get the force mesh from the API
    gsVector<index_t> forceKnotIDs;
    gsMatrix<> forceKnots;
    participant.getMeshVertexIDsAndCoordinates(ForceKnotMesh,forceKnotIDs,forceKnots);

    // Gives a full tensor product basis
    gsBasis<> * basis = knotMatrixToBasis<real_t>(forceKnots).get();

    gsMatrix<> forceControlPoints(basis->size(),2);
    forceControlPoints.setZero();


    // Gives the coefficients of the control points ONLY ON THE BOUNDARIES

    // IMPORTANT: THE control points now are all control points!!!
    gsVector<index_t> forceControlPointIDs;
    participant.getMeshVertexIDsAndCoordinates(ForceControlPointMesh, forceControlPointIDs,forceControlPoints);

    /*
     * TODO: ADD forceControlPointsBdr on the right rows of forceControlPoints
     */



    // // Step 2: Regenerate the geometry
    gsDebugVar(forceControlPoints.cols());
    gsDebugVar(basis->size());
    gsMultiPatch<> forceMesh; //Geometry object belongs to gsFunctionSet
    forceMesh.addPatch(give(basis->makeGeometry(forceControlPoints.transpose())));


    // NOTE: forceMesh should have domainDim 2!!

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // source function, rhs
    gsConstantFunction<> g(0.,beamLoad,2);

    gsConstantFunction<> gZero(0.,0.,2);
    // inflow velocity profile U(y) = 1.5*U_mean*y*(H-y)/(H/2)^2; channel height H = 0.41
    // gsFunctionExpr<> inflow(util::to_string(meanVelocity) + "*6*y*(0.41-y)/0.41^2",2);

    // containers for solution as IGA functions
    gsMultiPatch<> dispBeam;
    // boundary conditions: beam
    gsBoundaryConditions<> bcInfoBeam;
    for (index_t d = 0; d < 2; ++d)
        bcInfoBeam.addCondition(0,boundary::west,condition_type::dirichlet,0,d);
    // flow to beam interface: these Neumann boundary condtions contain references to flow and ALE solutions;
    // by updating them, we update the boundary load as well

    bcInfoBeam.addCondition(0,boundary::south,condition_type::neumann,&forceMesh.patch(0));
    bcInfoBeam.addCondition(0,boundary::east,condition_type::neumann,&forceMesh.patch(0));
    bcInfoBeam.addCondition(0,boundary::north,condition_type::neumann,&forceMesh.patch(0));

    //=============================================//
          // Setting assemblers and solvers //
    //=============================================//

    // elasticity solver: beam
    gsElasticityAssembler<real_t> elAssembler(geoBeam,basisDisplacement,bcInfoBeam,g);
    elAssembler.options().setReal("YoungsModulus",youngsModulus);
    elAssembler.options().setReal("PoissonsRatio",poissonsRatio);
    elAssembler.options().setInt("MaterialLaw",material_law::saint_venant_kirchhoff);
    gsMassAssembler<real_t> elMassAssembler(geoBeam,basisDisplacement,bcInfoBeam,gZero);
    elMassAssembler.options().setReal("Density",densitySolid);
    gsElTimeIntegrator<real_t> elTimeSolver(elAssembler,elMassAssembler);
    elTimeSolver.options().setInt("Scheme",time_integration::implicit_nonlinear);
    elTimeSolver.options().setReal("Beta",thetaSolid/2);
    elTimeSolver.options().setReal("Gamma",thetaSolid);
    gsInfo << "Initialized elasticity system with " << elAssembler.numDofs() << " dofs.\n";

    //=============================================//
             // Setting output and auxilary //
    //=============================================//

    // beam stress field
    gsPiecewiseFunction<> stresses;
    // isogeometric fields (geometry + solution)
    gsField<> displacementField(geoBeam,dispBeam);
    gsField<> stressField(geoBeam,stresses,true);

    // creating containers to plot several fields corresponding to the same geometry to one Paraview file
    std::map<std::string,const gsField<> *> fieldsBeam;
    fieldsBeam["Displacement"] = &displacementField;
    fieldsBeam["von Mises"] = &stressField;
    // paraview collection of time steps
    gsParaviewCollection collectionBeam("flapping_beam_2D_beam");

    // gsProgressBar bar;
    // gsStopwatch totalClock, iterClock;

    //=============================================//
                   // Initial condtions //
    //=============================================//

    // set initial velocity: zero free and fixed DoFs
    elTimeSolver.setDisplacementVector(gsMatrix<>::Zero(elAssembler.numDofs(),1));
    elTimeSolver.setVelocityVector(gsMatrix<>::Zero(elAssembler.numDofs(),1));

    // plotting initial condition
    elTimeSolver.constructSolution(dispBeam);
    elAssembler.constructCauchyStresses(dispBeam,stresses,stress_components::von_mises);

    if (numPlotPoints > 0)
    {
        gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flapping_beam_2D_beam",collectionBeam,0,1000);
    }

    //=============================================//
                   // Coupled simulation //
    //=============================================//

    real_t simTime = 0.;
    real_t numTimeStep = 0;
    real_t timeBeam = 0.;

    index_t timestep_checkpoint = 0;

    // totalClock.restart();

    gsInfo << "Running the simulation...\n";
    // Time integration loop
    while (participant.isCouplingOngoing())
    {
        if (participant.requiresWritingCheckpoint())
        {
            elTimeSolver.saveState();
        }

        participant.readData(ForceControlPointMesh,ForceControlPointData,forceControlPointIDs,forceControlPoints);

        /*
         * TODO: ADD forceControlPointsBdr on the right rows of forceControlPoints
         */

        forceMesh.patch(0).coefs() = forceControlPoints.transpose();

        // Perform a time integration step of the solid solver
        elTimeSolver.makeTimeStep(timeStep);

        // potentially adjust non-matching timestep sizes
        timeStep = std::min(timeStep,precice_dt);

        // Construct the displacment field
        elTimeSolver.constructSolution(dispBeam);

        gsDebugVar(geometryControlPoints.transpose());
        geometryControlPoints = dispBeam.patch(0).coefs().transpose();
        gsDebugVar(geometryControlPoints.transpose());
        gsDebugVar(dispBeam.patch(0).coefs().transpose());

        // Write the beam displacements to the fluid solver
        participant.writeData(GeometryControlPointMesh,GeometryControlPointData,geometryControlPointIDs,geometryControlPoints);

        // do the coupling
        precice_dt =participant.advance(timeStep);

        if (participant.requiresReadingCheckpoint())
        {
            elTimeSolver.recoverState();
            timeStep = timestep_checkpoint;
        }
        else
        {
            // gsTimeIntegrator advances the time step
            // advance variables
            timeBeam += timeStep;
            numTimeStep++;

            gsWriteParaviewMultiPhysicsTimeStep(fieldsBeam,"flapping_beam_2D_beam",collectionBeam,numTimeStep,1000);
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
