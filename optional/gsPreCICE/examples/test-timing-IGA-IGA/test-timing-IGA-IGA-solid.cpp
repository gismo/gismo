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
     * This participant receives mesh information (knot vector, control points) from the solid,
     * and it creates its own mesh based on the information. 
     * And writes pressure, reads displacements from the solid participant.
     * The follow meshes and data are made available:
     *
     * - Meshes:
     *   + GeometryControlPointMesh:        This mesh contains the control points as mesh vertices.
     *   + GeometryKnotMesh:                This mesh contains the knots as mesh vertices.
     *   + ForceKnotMesh:           		This mesh contains the knots as mesh vertices for force mesh.
     *   + ForceControlPointMesh:   		This mesh contains the control points as mesh vertices for force mesh.
     *
     *
     * - Data:
     *   + GeometryControlPointData:        This data is defined on the GeometryControlPointMesh and stores the displacement of the control points
     *	 + GeometryKnotData
     *   + ForceControlPointData:   		This data is defined on the ForceControlPointMesh and stores the change of pressure/stress
     *	 + ForceKnotData
     */

    std::string GeometryKnotMesh               	= "GeometryKnotMesh";
    std::string GeometryControlPointMesh		= "GeometryControlPointMesh";
    std::string ForceKnotMesh           		= "ForceKnotMesh";
    std::string ForceControlPointMesh   		= "ForceControlPointMesh";

    // std::string StressData       = "StressData";
    std::string GeometryControlPointData        = "GeometryControlPointData";
    // std::string GeometryKnotData				= "GeometryKnotData";
    // std::string ForceKnotData           		= "ForceKnotMeshData";
    std::string ForceControlPointData   		= "ForceControlPointData";

    gsMatrix<> bbox(3,2);
    bbox << -1e300, 1e300, // X dimension limits
            -1e300, 1e300, // Y dimension limits
            -1e300, 1e300; // Z dimension limits
    bbox.transposeInPlace();
    bbox.resize(1,bbox.rows()*bbox.cols());
    participant.setMeshAccessRegion(ForceControlPointMesh, bbox);

    // Setup geometry control point mesh
    gsVector<index_t> geometryControlPointIDs;
    gsMatrix<> geometryControlPoints = patches.patch(0).coefs().transpose();
    participant.addMesh(GeometryControlPointMesh, geometryControlPoints, geometryControlPointIDs);


    // Setup the geometry knot mesh
    gsMatrix<> geometryKnots = knotsToMatrix(bases.basis(0)); 


    participant.addMesh(GeometryKnotMesh,geometryKnots);

    real_t precice_dt = participant.initialize();

    //G[et the force mesh from the coupling interface
    gsVector<index_t> forceKnotIDs;
    gsMatrix<> forceKnots;
    participant.getMeshVertexIDsAndCoordinates(ForceKnotMesh,forceKnotIDs,forceKnots);


    gsVector<index_t> forceControlPointIDs;
    gsMatrix<> forceControlPoints;
    participant.getMeshVertexIDsAndCoordinates(ForceControlPointMesh, forceControlPointIDs,forceControlPoints);
    



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

    while (participant.isCouplingOngoing())
    {
        if (participant.requiresWritingCheckpoint())
        {
            gsInfo<<"Checkpoint written:\n";

        }
        participant.readData(ForceControlPointMesh,ForceControlPointData,forceControlPointIDs,forceControlPoints);

        if (get_readTime)
            t_read += participant.readTime();

        // potentially adjust non-matching timestep sizes
        dt = std::min(dt,precice_dt);


        geometryControlPoints.setRandom();
        participant.writeData(GeometryControlPointMesh,GeometryControlPointData,geometryControlPointIDs,geometryControlPoints);
        
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
        gsInfo << "Solid Read time: " << t_read << "\n";
    }

    if (get_writeTime)
    {
        gsInfo << "Solid Write time: " << t_write << "\n";
    }

    return  EXIT_SUCCESS;
 }   
