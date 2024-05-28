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
    bool plot = false;
    bool write = false;
    bool get_readTime = false;
    bool get_writeTime = false;
    index_t plotmod = 1;
    index_t loadCase = 0;
    std::string precice_config("../precice_config.xml");

    gsCmdLine cmd("Coupled heat equation using PreCICE.");
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addInt("l","loadCase", "Load case: 0=constant load, 1='spring' load", loadCase);
    cmd.addSwitch("write", "Create a file with point data", write);
    cmd.addSwitch("readTime", "Get the read time", get_readTime);
    cmd.addSwitch("writeTime", "Get the write time", get_writeTime);
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
    // std::string GeometryKnotData                = "GeometryKnotData";
    // std::string ForceKnotData                   = "ForceKnotMeshData";

    gsMatrix<> bbox(3,2);
    bbox << -1e300, 1e300, // X dimension limits
            -1e300, 1e300, // Y dimension limits
            -1e300, 1e300; // Z dimension limits
    bbox.transposeInPlace();
    bbox.resize(1,bbox.rows()*bbox.cols());

    // Step 2: Add force mesh

    // Step 2.1: Define force knot mesh
    gsVector<index_t> forceKnotIDs;
    gsKnotVector<> kv(0,1,10,3,1); //begin, end, # interiors, mult, by default =1
    gsTensorBSplineBasis<2, real_t> forceBasis(kv,kv);
    gsMatrix<> forceControlPoints(forceBasis.size(), 3);
    forceControlPoints.setZero();
    forceControlPoints.col(2).setConstant(-5e3);
    forceControlPoints.transposeInPlace();

    gsMatrix<> forceKnots = knotsToMatrix(forceBasis);




    // gsMultiBasis<> bases(forceMesh);
    // gsMatrix<> forceKnots = knotsToMatrix(bases.basis(0));
    // gsDebugVar(forceKnots.dim());
    // 
    


    // Regenerate the geometry in another solver
    // forceMesh.addPatch(give(forceBasis.makeGeometry(forceControlPoints)));

    // Step 2.1: Define force control point mesh
    gsVector<index_t> forceControlPointIDs;
    participant.addMesh(ForceKnotMesh,forceKnots,forceKnotIDs);
    participant.addMesh(ForceControlPointMesh,forceControlPoints,forceControlPointIDs);
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

    gsMultiPatch<> mp, deformation;
    gsDebugVar(*geometryKnotBasis);
    gsDebugVar(geometryKnotBasis->size());
    gsDebugVar(geometryControlPoints.rows());
    mp.addPatch(give(geometryKnotBasis->makeGeometry(geometryControlPoints.transpose())));

    deformation = mp;


    real_t t=0;
    real_t dt = precice_dt;
    real_t t_read = 0;
    real_t t_write = 0;

    index_t timestep = 0;

    // Define the solution collection for Paraview
    gsParaviewCollection collection("./output/solution");

    gsMatrix<> points(2,1);
    points.col(0)<<0.5,1;

    gsStructuralAnalysisOutput<real_t> writer("./output/pointData.csv",points);
    writer.init({"x","y","z"},{"time"});

    // Time integration loop
    while(participant.isCouplingOngoing())
    {

        // Read control point displacements
        participant.readData(GeometryControlPointMesh, GeometryControlPointData, geometryControlPointIDs, geometryControlPoints);
        deformation.patch(0).coefs() = geometryControlPoints.transpose();
        if (get_readTime)
            t_read += participant.readTime();

        /*
         * Two projection approaches:
         * - gsQuasiInterpolate<real_t>::localIntpl(forceBasis, deformation.patch(0), coefs);
         * - gsL2Projection<real_t>::projectFunction(forceBasis, deformation, mp, coefs)
         *
         * Check if forceBasis + coefs gives the same as deformation.patch(0)
         */
        if (loadCase == 0)
        {
            gsInfo << "Load case 0: Constant load\n";
        }    
        else if (loadCase == 1)
        {
            gsInfo << "Load case 1: Spring load\n";
            if(timestep > 0)
            {
                gsMatrix<> coefs;
                gsL2Projection<real_t>::projectFunction(forceBasis, deformation, mp, coefs);
                coefs.resize(forceControlPoints.rows(),forceControlPoints.cols());
                forceControlPoints.row(2) = -coefs.row(2);
            }    
        }
        else
        {
            GISMO_ERROR("Load case "<<loadCase<<" unknown.");
        }


        participant.writeData(ForceControlPointMesh,ForceControlPointData,forceControlPointIDs,forceControlPoints);

        if (get_writeTime)
            t_write +=participant.writeTime();
        // bool requiresWritingCheckpoint and requiresReadingCheckpoint are **REQUIRED** in PreCICE. 
        // And the requiresWritingCheckpoint() should be called before advance()
        gsMatrix<> pointDataMatrix;
        gsMatrix<> otherDataMatrix(1,1);
        if (participant.requiresWritingCheckpoint())
        {
            gsInfo<<"Reading Checkpoint\n";
        }

        // do the coupling
        precice_dt = participant.advance(dt);

        dt = std::min(precice_dt, dt);

        if (participant.requiresReadingCheckpoint())
        {    
            gsInfo<<"Writing Checkpoint \n";
        }
        else
        {
            t += dt;
            timestep++;

            gsField<> solution(mp,deformation);
            if (timestep % plotmod==0 && plot)
            {
                // solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
                std::string fileName = "./output/solution" + util::to_string(timestep);
                gsWriteParaview<>(solution, fileName, 500);
                fileName = "solution" + util::to_string(timestep) + "0";
                collection.addTimestep(fileName,t,".vts");
            }
            if (write)
            {
                solution.patch(0).eval_into(points,pointDataMatrix);
                otherDataMatrix<<t;
                writer.add(pointDataMatrix,otherDataMatrix);
            }    

            // solution.patch(0).eval_into(points,pointDataMatrix);
            // otherDataMatrix<<time;
            // writer.add(pointDataMatrix,otherDataMatrix);
        }
    }
    if (get_readTime)
    {
        gsInfo << "Fluid Read time: " << t_read << "\n";
    }

    if (get_writeTime)
    {
        gsInfo << "Fluid Write time: " << t_write << "\n";
    }

    if (plot)
    {
        // Refer to gsParaviewCollection
        collection.save();
    }



	return EXIT_SUCCESS;
}	
