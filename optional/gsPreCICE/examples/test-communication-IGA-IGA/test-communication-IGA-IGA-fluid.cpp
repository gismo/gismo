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

    //! [Parse command line]
    std::string precice_config("../precice_config.xml");

    gsCmdLine cmd("Coupled heat equation using PreCICE.");
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");
    //! [Read input file]
    std::string participantName = "Fluid";
    gsPreCICE<real_t> participant(participantName, precice_config);


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
     *   + GeometryControlPointData:        This data is defined on the ControlPointMesh and stores the displacement of the control points
     *   + ForceControlPointData:   		This data is defined on the ForceControlPointMesh and stores the change of pressure/stress
     */


    // Meshes
    std::string GeometryControlPointMesh		= "GeometryControlPointMesh";
    std::string GeometryKnotMesh				= "GeometryKnotMesh";


    std::string ForceKnotMesh					= "ForceKnotMesh";
    std::string ForceControlPointMesh			= "ForceControlPointMesh";


    // Data
    std::string GeometryControlPointData		= "GeometryControlPointData";
    std::string ForceControlPointData			= "ForceControlPointData";

    gsMatrix<> bbox(3,2);
    bbox << -1e300, 1e300, // X dimension limits
            -1e300, 1e300, // Y dimension limits
            -1e300, 1e300; // Z dimension limits
    bbox.transposeInPlace();
    bbox.resize(1,bbox.rows()*bbox.cols());

    // Step 2: Add force mesh

    // Step 2.1: Define force knot mesh
    gsVector<index_t> forceKnotIDs;
    gsMultiPatch<> forceMesh;
    gsKnotVector<> kv(0,1,10,3,1); //begin, end, # interiors, mult, by default =1
    gsTensorBSplineBasis<2, real_t> forceBasis(kv,kv);
    gsMatrix<> forceControlPoints(forceBasis.size(), 2);
    forceControlPoints.setConstant(-1e4);
    gsDebugVar(forceControlPoints.dim());
    forceMesh.addPatch(give(forceBasis.makeGeometry(forceControlPoints)));


    gsMultiBasis<> bases(forceMesh);
    gsMatrix<> forceKnots = knotsToMatrix(bases.basis(0));
    gsDebugVar(forceKnots.dim());


    // Regenerate the geometry in another solver
    // forceMesh.addPatch(give(forceBasis.makeGeometry(forceControlPoints)));

    // Step 2.1: Define force control point mesh
    gsVector<index_t> forceControlPointIDs;
    participant.addMesh(ForceKnotMesh,forceKnots,forceKnotIDs);
    participant.addMesh(ForceControlPointMesh,forceControlPoints.transpose(),forceControlPointIDs);
    participant.setMeshAccessRegion(GeometryControlPointMesh, bbox);


    real_t precice_dt = participant.initialize();

    // Step 1: Collect geometry from PreCICE
    gsVector<index_t> geometryKnotIDs;
    gsMatrix<> geometryKnots;
    participant.getMeshVertexIDsAndCoordinates(GeometryKnotMesh, geometryKnotIDs, geometryKnots);
    gsDebugVar(geometryKnots);


    gsBasis<> * geometryKnotBasis = knotMatrixToBasis<real_t>(geometryKnots).get();

    gsVector<index_t> geometryControlPointIDs;
    gsMatrix<> geometryControlPoints;
    participant.getMeshVertexIDsAndCoordinates(GeometryControlPointMesh, geometryControlPointIDs, geometryControlPoints);
    gsDebugVar(geometryControlPoints.dim());




	return EXIT_SUCCESS;
}	