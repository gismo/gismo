/** @file test-communication-IGA-IGA-fluid.cpp
 * 
 * 	@brief test communicate spline through preCICE. To see if the spline force function can be reconstructed
 * 	
 * 	This file is a part of G+Smo library.
 * 
 * 	This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 * 
 * 	Author(s): J. Li
 * 
 **/

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEUtils.h>
#include <gsPreCICE/gsPreCICEFunction.h>
#include <gsPreCICE/gsLookupFunction.h>
// #include <gsPreCICE/gsPreCICEVectorFunction.h>

// #include <gsElasticity/gsMassAssembler.h>
// #include <gsElasticity/gsElasticityAssembler.h>


using namespace gismo;

// template<typename T>
// T ProjectForceMesh(      const gsMatrix<>       & forceControlPoints, 
//                          const gsMultiBasis<>   & forceBasis,
//                          const gsMultiBasis<>   & geometryBasis,
//                          const index_t          & targetDim,
//                          gsMatrix<T>            & results)
// {
//     gsExprAssembler<> A(1,1);

//     gsMatrix<T> solVector;

//     // Set the discretization space
//     space u = A.getSpace(forceBasis, targetDim);

//     solution sol = A.getSolution(u, solVector);
//     auto f = A.getCoordinates(u);

//     u.setup(0);
//     A.initSystem();



//     gsMatrix<> matrixForce, matrixGeometry;
// }


int main(int argc, char *argv[])
{
    bool plot = false;
    bool get_readTime = false;
    bool get_writeTime = false;
    index_t plotmod = 1;
    index_t numRefine  = 1;
    index_t numElevate = 0;
    int method = 3; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson, 6 RK4

    real_t rho = 3000;
    real_t E = 4e6;
    real_t nu = 0.3;

    //! [Parse command line]
    std::string precice_config("../precice_config.xml");

    gsCmdLine cmd("Coupled heat equation using PreCICE.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m", "method","1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson",method);
    cmd.addSwitch("readTime", "Get the read time", get_readTime);
    cmd.addSwitch("writeTime", "Get the write time", get_writeTime);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");
    //! [Read input file]
    std::string participantName = "Solid";
    gsPreCICE<real_t> participant(participantName, precice_config);

    //Generate domain
    gsMultiPatch<> patches;
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0,0.0,0.5,1.0));

    //Embed the 2D geometry to 3D
    gsMultiPatch<> solutions;
    patches.addAutoBoundaries();
    patches.embed(3);

    // Create bases
        // p-refine
    if (numElevate!=0)
        patches.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        patches.uniformRefine();

    // Create bases
    gsMultiBasis<> bases(patches);//true: poly-splines (not NURBS)

    gsInfo << "Patches: "<< patches.nPatches() <<", degree: "<< bases.minCwiseDegree() <<"\n";



    /*
     * Data initialization
     *
     * This participant manages the geometry. The follow meshes and data are made available:
     *
     * - Meshes:
     *   + KnotMesh             This mesh contains the knots as mesh vertices
     *   + ControlPointMesh:    This mesh contains the control points as mesh vertices
     *   + ForceMesh:           This mesh contains the integration points as mesh vertices
     *
     * - Data:
     *   + ControlPointData:    This data is defined on the ControlPointMesh and stores the displacement of the control points
     *   + ForceData:           This data is defined on the ForceMesh and stores pressure/forces
     */
    std::string KnotMesh        = "KnotMesh";
    std::string ControlPointMesh= "ControlPointMesh";
    std::string ForceMesh       = "ForceMesh";

    std::string ControlPointData= "ControlPointData";
    std::string ForceData       = "ForceData";

    // Step 1: write the meshes to PreCICE
    // Step 1a: KnotMesh
    // get the knots in a matrix, using the utility function knotsToMatrix
    gsMatrix<> knots = knotsToMatrix(bases.basis(0)); //
    participant.addMesh(KnotMesh,knots);

    // Step 1b: ControlPointMesh
    // get the control points, in the format where every column is a control point
    gsVector<index_t> controlPointIDs; // needed for writing
    gsMatrix<> controlPoints = patches.patch(0).coefs().transpose();
    participant.addMesh(ControlPointMesh,controlPoints,controlPointIDs);

    // Step 1c: ForceMesh
    // Get the quadrature nodes on the coupling interface
    gsOptionList quadOptions = gsExprAssembler<>::defaultOptions();

    // Get the quadrature points
    gsMatrix<> quadPoints = gsQuadrature::getAllNodes(bases.basis(0),quadOptions);
    gsVector<index_t> quadPointsIDs;
    participant.addMesh(ForceMesh,quadPoints,quadPointsIDs);
    gsMatrix<> quadPointsData(3,quadPoints.cols());
    quadPointsData.setZero();


    // Step 2 (not needed???)
    // Needed for direct mesh coupling
    // gsMatrix<> bbox = bases.basis(0).support();
    // bbox.transposeInPlace();
    // bbox.resize(1,bbox.rows()*bbox.cols());
    // participant.setMeshAccessRegion(ForceMesh,bbox);

    // gsMatrix<> bbox(3,2);
    // bbox << 0,1,
    //         0,1,
    //         0,1;
    // bbox.transposeInPlace();
    // bbox.resize(1,bbox.rows()*bbox.cols());
    // participant.setMeshAccessRegion(ForceMesh,bbox);
    // Step 3: initialize the participant
    real_t precice_dt = participant.initialize();



    

    // // Step 2: Regenerate the geometry


    // // Define boundary condition for solid mesh
    // gsBoundaryConditions<> bcInfo;
    // //Bottom Side
    // bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, -1);
    // bcInfo.addCondition(0, boundary::south, condition_type::clamped, nullptr, 2);

    // // Assign geometry map
    // bcInfo.setGeoMap(patches);


     // Time integration loop

    real_t t=0;
    real_t dt = precice_dt;
    real_t t_read = 0;
    real_t t_write = 0;

    index_t timestep = 0;

    real_t read_time = 0;
    real_t write_time = 0;
    index_t read_calls = 0;
    index_t write_calls = 0;
    gsStopwatch timer;
    while (participant.isCouplingOngoing())
    {
        if (participant.requiresWritingCheckpoint())
        {
            gsInfo<<"Checkpoint written:\n";

        }

        timer.restart();
        participant.readData(ForceMesh,ForceData,quadPointsIDs,quadPointsData);
        read_time += timer.stop();
        read_calls++;

        if (get_readTime)
            t_read += participant.readTime();

        // potentially adjust non-matching timestep sizes
        dt = std::min(dt,precice_dt);


        controlPoints.setRandom();

        timer.restart();
        participant.writeData(ControlPointMesh,ControlPointData,controlPointIDs,controlPoints);
        write_time += timer.stop();
        write_calls++;

        if (get_writeTime)
            t_write +=participant.writeTime();


        // do the coupling
        precice_dt =participant.advance(dt);

        if (participant.requiresReadingCheckpoint())
        {
            gsInfo <<"Called Read Checkpoint \n";
        }
        else
        {
            // gsTimeIntegrator advances the time step                   
            // advance variables
            t += dt;
            timestep++;

        }
    }
    if (get_readTime)
    {
        gsInfo << "Solid Read time: "<<read_time<<" in "<<read_calls<<"calls.\n";
    }

    if (get_writeTime)
    {
        gsInfo << "Solid Write time: " <<write_time<<" in "<<write_calls<<"calls.\n";
    }

    std::ofstream file("./IGA-vertex-solid_times.csv");
    file<<"read time, read calls, write time, write calls\n";
    file<<read_time<<","<<read_calls<<","<<write_time<<","<<write_calls<<"\n";
    file.close();

    return  EXIT_SUCCESS;

 }   
