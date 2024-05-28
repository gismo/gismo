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
    std::string KnotMesh        = "KnotMesh";
    std::string ControlPointMesh= "ControlPointMesh";
    std::string ForceMesh       = "ForceMesh";

    std::string ControlPointData= "ControlPointData";
    std::string ForceData       = "ForceData";

    gsMatrix<> bbox(3,2);
    bbox << -1e300, 1e300, // X dimension limits
            -1e300, 1e300, // Y dimension limits
            -1e300, 1e300; // Z dimension limits
    bbox.transposeInPlace();
    bbox.resize(1,bbox.rows()*bbox.cols());

    participant.setMeshAccessRegion(ControlPointMesh,bbox);

    // No mesh is defined from this participant, so we initialize immediately
    real_t precice_dt = participant.initialize();

    // Step 1: Collect geometry from PreCICE



    gsVector<index_t> controlPointIDs;
    gsMatrix<> controlPoints;
    participant.getMeshVertexIDsAndCoordinates(ControlPointMesh,controlPointIDs,controlPoints);



    gsVector<index_t> quadPointIDs;
    gsMatrix<> quadPoints;
    participant.getMeshVertexIDsAndCoordinates(ForceMesh,quadPointIDs,quadPoints);
    gsMatrix<> quadPointData(controlPoints.rows(),quadPoints.cols());



    real_t t=0;
    real_t dt = precice_dt;
    real_t t_read = 0;
    real_t t_write = 0;

    index_t timestep = 0;


    // Time integration loop
    while(participant.isCouplingOngoing())
    {

        
    	if (participant.requiresWritingCheckpoint())
        {
            gsInfo<<"Checkpoint written:\n";

        }

        // Read control point displacements
        participant.readData(ControlPointMesh,ControlPointData,controlPointIDs,controlPoints);

        if (get_readTime)
            t_read += participant.readTime();

    
        quadPointData.setRandom();
        participant.writeData(ForceMesh,ForceData,quadPointIDs,quadPointData);

        if (get_writeTime)
            t_write +=participant.writeTime();
        // bool requiresWritingCheckpoint and requiresReadingCheckpoint are **REQUIRED** in PreCICE. 
        // And the requiresWritingCheckpoint() should be called before advance()

        if (participant.requiresWritingCheckpoint())
        {
            gsInfo<<"Reading Checkpoint\n";
        }

        // do the coupling
        precice_dt = participant.advance(dt);

        dt = std::min(precice_dt, dt);

       
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
        gsInfo << "Fluid Read time: " << t_read << "\n";
    }

    if (get_writeTime)
    {
        gsInfo << "Fluid Write time: " << t_write << "\n";
    }



	return EXIT_SUCCESS;
}	
