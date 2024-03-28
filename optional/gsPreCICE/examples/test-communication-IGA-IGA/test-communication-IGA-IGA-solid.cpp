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
// #include <gsPreCICE/gsPreCICEVectorFunction.h>

#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElasticityAssembler.h>


#ifdef gsStructuralAnalysis_ENABLED
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicExplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicImplicitEuler.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicNewmark.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicBathe.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicWilson.h>
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsDynamicRK4.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>
#endif

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false;
    index_t plotmod = 1;
    index_t numRefine  = 1;
    index_t numElevate = 0;
    int method = 3; // 1: Explicit Euler, 2: Implicit Euler, 3: Newmark, 4: Bathe, 5: Wilson, 6 RK4

    //! [Parse command line]
    std::string precice_config("../precice_config.xml");

    gsCmdLine cmd("Coupled heat equation using PreCICE.");
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
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
    std::string GeometryKnotData				= "GeometryKnotData";
    std::string ForceKnotData           		= "ForceKnotMeshData";
    std::string ForceControlPointData   		= "ForceControlPointData";

    gsMatrix<> bbox(2,2);
    bbox << -1e300, 1e300, // X dimension limits
            -1e300, 1e300; // Y dimension limits
    bbox.transposeInPlace();
    bbox.resize(1,bbox.rows()*bbox.cols());
    participant.setMeshAccessRegion(ForceControlPointMesh, bbox);

    // Setup geometry control point mesh
    gsVector<index_t> geometryControlPointIDs;
    participant.addMesh(GeometryControlPointMesh, patches.patch(0).coefs().transpose(), geometryControlPointIDs);
    
    // Setup the geometry knot mesh
    gsMatrix<> geometryKnots = knotsToMatrix(bases.basis(0)); 
    gsDebugVar(geometryKnots);
    participant.addMesh(GeometryKnotMesh,geometryKnots);

    real_t precice_dt = participant.initialize();

    //Get the force mesh from the coupling interface
    gsVector<index_t> forceKnotIDs;
    gsMatrix<> forceKnots;
    participant.getMeshVertexIDsAndCoordinates(ForceKnotMesh,forceKnotIDs,forceKnots);
    gsDebugVar(forceKnots);

    gsBasis<> * basis = knotMatrixToBasis<real_t>(forceKnots).get();

    gsVector<index_t> forceControlPointIDs;
    gsMatrix<> forceControlPoints;
    participant.getMeshVertexIDsAndCoordinates(ForceControlPointMesh, forceControlPointIDs,forceControlPoints);

    gsDebugVar(forceControlPoints.dim());

    // // Step 2: Regenerate the geometry
    // gsMultiPatch<> forceMesh;
    // forceMesh.addPatch(give(basis->makeGeometry(forceControlPoints.transpose())));

	return EXIT_SUCCESS;
}
