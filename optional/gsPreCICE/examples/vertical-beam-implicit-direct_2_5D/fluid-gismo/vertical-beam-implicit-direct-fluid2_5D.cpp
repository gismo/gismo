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
#include <gsPreCICE/gsPreCICEUtils.h>
#include <gsPreCICE/gsPreCICEFunction.h>
// #include <gsPreCICE/gsPreCICEVectorFunction.h>

#include <gsElasticity/gsMassAssembler.h>
#include <gsElasticity/gsElasticityAssembler.h>

#ifdef gsStructuralAnalysis_ENABLED
#include <gsStructuralAnalysis/src/gsDynamicSolvers/gsTimeIntegrator.h>
#include <gsStructuralAnalysis/src/gsStructuralAnalysisTools/gsStructuralAnalysisUtils.h>
#endif

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t plotmod = 1;
    std::string precice_config;

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  
    //! [Read input file]
    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");

    /*
     * Initialize the preCICE participant
     *
     *
     */
    std::string participantName = "Fluid";
    gsPreCICE<real_t> participant(participantName, precice_config);

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

    gsMatrix<> bbox(3,2); //We need to specify the dimension before initialization
    // bbox.col(0).setConstant(-1e300);
    // bbox.col(1).setConstant(1e300);
    bbox << -1e300, 1e300, // X dimension limits
            -1e300, 1e300, // Y dimension limits
            -1e300, 1e300; // Z dimension limits
    bbox.transposeInPlace();
    bbox.resize(1,bbox.rows()*bbox.cols());
    participant.setMeshAccessRegion(ControlPointMesh,bbox);

    // No mesh is defined from this participant, so we initialize immediately
    real_t precice_dt = participant.initialize();

    /*
     * Collect the geometry
     *
     */

    gsVector<index_t> knotIDs;
    gsMatrix<> knots;
    participant.getMeshVertexIDsAndCoordinates(KnotMesh,knotIDs,knots);

    gsBasis<> * basis = knotMatrixToBasis<real_t>(knots).get();

    gsVector<index_t> controlPointIDs;
    gsMatrix<> controlPoints;
    participant.getMeshVertexIDsAndCoordinates(ControlPointMesh,controlPointIDs,controlPoints);

    gsMultiPatch<> mp, deformation;
    // mp.addPatch(give(basis->makeGeometry(CPs.transpose())));

    mp.addPatch(give(basis->makeGeometry(controlPoints.transpose())));
    deformation = mp;
    deformation.patch(0).coefs().setZero();

    gsVector<index_t> quadPointIDs;
    gsMatrix<> quadPoints;
    participant.getMeshVertexIDsAndCoordinates(ForceMesh,quadPointIDs,quadPoints);
    gsMatrix<> quadPointData(mp.geoDim(),quadPoints.cols());

    real_t t = 0, dt = precice_dt;
    index_t timestep = 0;
    // Define the solution collection for Paraview
    gsParaviewCollection collection("./output/solution");
    // Time integration loop
    while (participant.isCouplingOngoing())
    {
        if (participant.requiresWritingCheckpoint())
            gsInfo<<"Writing Checkpoint\n";

        // Read control point displacements
        participant.readData(ControlPointMesh,ControlPointData,controlPointIDs,controlPoints);
        deformation.patch(0).coefs() = controlPoints.transpose();

        // Write data at the quadrature points
        quadPointData.setZero();
        quadPointData.row(1).setConstant(-1e3);
        gsDebugVar(quadPointData);
        participant.writeData(ForceMesh,ForceData,quadPointIDs,quadPointData);

        // do the coupling
        precice_dt =participant.advance(dt);

        dt = std::min(precice_dt, dt);

        if (participant.requiresReadingCheckpoint())
            gsInfo<<"Reading Checkpoint\n";
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

            // solution.patch(0).eval_into(points,pointDataMatrix);
            // otherDataMatrix<<time;
            // writer.add(pointDataMatrix,otherDataMatrix);
        }
    }

    if (plot)
    {
        collection.save();
    }


    return  EXIT_SUCCESS;
}
