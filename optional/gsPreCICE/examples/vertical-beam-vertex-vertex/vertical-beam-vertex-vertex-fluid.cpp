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
    index_t loadCase = 0;
    bool get_readTime = false;
    bool get_writeTime = false;
    std::string precice_config;

    gsCmdLine cmd("Flow over heated plate for PreCICE.");
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt("m","plotmod", "Modulo for plotting, i.e. if plotmod==1, plots every timestep", plotmod);
    cmd.addInt("l","loadCase", "Load case: 0=constant load, 1='spring' load", loadCase);
    cmd.addSwitch("readTime", "Get the read time", get_readTime);
    cmd.addSwitch("writeTime", "Get the write time", get_writeTime);
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
    // Get the quadrature points
    gsVector<index_t> FluidMeshIDs;
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
    real_t t_read = 0;
    real_t t_write = 0;
    index_t timestep = 0;
    // Define the solution collection for Paraview
    gsParaviewCollection collection("./output/solution");

    gsMatrix<> StressPointData(patches.geoDim(),mesh.cols());

    // Time integration loop
    while (participant.isCouplingOngoing())
    {
        if (participant.requiresWritingCheckpoint())
            gsInfo<<"Writing Checkpoint\n";

        // Read control point displacements
        gsMatrix<> meshPointDisplacements;
        participant.readData(FluidMesh,DisplacementData,FluidMeshIDs,meshPointDisplacements);

        if (get_readTime)
            t_read += participant.readTime();


        if (loadCase==0)
        {
            // Write data at the quadrature points
            StressPointData.setZero();
            StressPointData.row(2).setConstant(-5e3);
        }
        else if (loadCase==1)
        {
            StressPointData.setZero();
            for (index_t k = 0; k!=meshPointDisplacements.cols(); k++)
            {
                StressPointData(2,k) = -meshPointDisplacements(2,k);
            }
            // impact loading at t==0
            if (t==0)
                StressPointData.row(2).array() += -5e3;
        }
        else
            GISMO_ERROR("Load case "<<loadCase<<" unknown.");


        participant.writeData(FluidMesh,StressData,FluidMeshIDs,StressPointData);

        if (get_writeTime)
            t_write +=participant.writeTime();

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

    if (get_readTime)
    {
        gsInfo << "Read time: " << t_read << "\n";
    }

    if (get_writeTime)
    {
        gsInfo << "Write time: " << t_write << "\n";
    }

    if (plot)
    {
        // collection.save();
    }


    return  EXIT_SUCCESS;
}
