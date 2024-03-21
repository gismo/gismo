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
    index_t plotmod = 1;
    std::string precice_config;

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  
    //! [Read input file]
    GISMO_ASSERT(gsFileManager::fileExists(precice_config),"No precice config file has been defined");

    // Generate domain
    gsMultiPatch<> patches;
    patches.addPatch(gsNurbsCreator<>::BSplineRectangle(0.0,0.0,0.5,1.0));

    patches.addAutoBoundaries();
    patches.embed(3);

    // Generate discrete fluid mesh
    gsMatrix<> bbox = patches.patch(0).support();



    int numPoints = 100;
    gsMatrix<> pointGrid = gsPointGrid<>(bbox, numPoints);
 
    gsMatrix<> mesh = patches.patch(0).eval(pointGrid);



    // Embed the 2D geometry to 3D

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
     * This participant creates its own mesh, and it writes and reads displacements and stresses on its mesh.
     * The follow meshes and data are made available:
     *
     * - Meshes:
     *   + FluidMesh:               This mesh contains the integration points as mesh vertices
     *
     * - Data:
     *   + DisplacementData:    This data is defined on the FluidMesh and stores the displacement at the integration points
     *   + StressData:           This data is defined on the FluidMesh and stores pressure/forces at the integration points
     */
    std::string FluidMesh        = "FluidMesh";
    std::string StressData       = "StressData";
    std::string DisplacementData = "DisplacementData";


    // Step 1: FluidMesh
    // Get the quadrature nodes on the coupling interface
    gsOptionList quadOptions = gsExprAssembler<>::defaultOptions();

    // Get the quadrature points
    gsVector<index_t> FluidMeshIDs;
    // gsMatrix<> quadPoints = gsQuadrature::getAllNodes(bases.basis(0),quadOptions);
    gsDebugVar(mesh);
    participant.addMesh(FluidMesh,mesh,FluidMeshIDs);

    // Step 2: initialize the participant
    real_t precice_dt = participant.initialize();
    gsDebugVar("Initialize Fluid");

    /*
     * Collect the geometry
     *
     */

    real_t t = 0, dt = precice_dt;
    index_t timestep = 0;
    // Define the solution collection for Paraview
    gsParaviewCollection collection("./output/solution");

    gsMatrix<> StressPointData(patches.geoDim(),mesh.cols());
    gsDebugVar("Got here 2");

    // Time integration loop
    while (participant.isCouplingOngoing())
    {
        if (participant.requiresWritingCheckpoint())
            gsInfo<<"Writing Checkpoint\n";

        // Read control point displacements
        gsMatrix<> meshPointDisplacements;
        gsDebugVar("Got here 2");
        participant.readData(FluidMesh,DisplacementData,FluidMeshIDs,meshPointDisplacements);

        // Write data at the quadrature points
        StressPointData.setZero();
        StressPointData.row(3).setConstant(-1e4);
        gsDebugVar("Got here 2");
        participant.writeData(FluidMesh,StressData,FluidMeshIDs,StressPointData);
        gsDebugVar("Got here 2");

        // do the coupling
        precice_dt =participant.advance(dt);

        dt = std::min(precice_dt, dt);

        if (participant.requiresReadingCheckpoint())
            gsInfo<<"Reading Checkpoint\n";
        else
        {
            t += dt;
            timestep++;

            // gsField<> solution(mp,deformation);
            if (timestep % plotmod==0 && plot)
            {
                // // solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
                // std::string fileName = "./output/solution" + util::to_string(timestep);
                // gsWriteParaview<>(solution, fileName, 500);

                // fileName = "solution" + util::to_string(timestep) + "0";
                // collection.addTimestep(fileName,t,".vts");
            }

            // solution.patch(0).eval_into(points,pointDataMatrix);
            // otherDataMatrix<<time;
            // writer.add(pointDataMatrix,otherDataMatrix);
        }
    }

    if (plot)
    {
        // collection.save();
    }


    return  EXIT_SUCCESS;
}
