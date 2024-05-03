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
     * This participant receives mesh information (knot vector, control points) from the solid,
     * and it creates its own mesh based on the information. 
     * And writes pressure, reads displacements from the solid participant.
     * The follow meshes and data are made available:
     *
     * - Meshes:
     *   + ControlPointMesh:        This mesh contains the control points as mesh vertices.
     *   + KnotMesh:                This mesh contains the knots as mesh vertices.
     *   + ForceKnotMesh:           This mesh contains the knots as mesh vertices for force mesh.
     *   + ForceControlPointMesh:   This mesh contains the control points as mesh vertices for force mesh.
     *
     *
     * - Data:
     *   + ControlPointData:        This data is defined on the ControlPointMesh and stores the displacement of the control points
     *   + ForceCOntrolPointData:   This data is defined on the ForceControlPointMesh and stores the change of pressure/stress
     */

    // Meshes
    std::string GeometryControlPointMesh        = "GeometryControlPointMesh";
    std::string GeometryKnotMesh                = "GeometryKnotMesh";
    std::string ForceKnotMesh                   = "ForceKnotMesh";
    std::string ForceControlPointMesh           = "ForceControlPointMesh";


    // Data
    std::string GeometryControlPointData        = "GeometryControlPointData"; 
    std::string ForceControlPointData           = "ForceControlPointData";

    // std::string ForceData        = "ForceData"; // There is no force data in this context, the force information will be written as a spline geometry


    gsMatrix<> bbox(3,2); //We need to specify the dimension before initialization
    // bbox.col(0).setConstant(-1e300);
    // bbox.col(1).setConstant(1e300);
    bbox << -1e300, 1e300, // X dimension limits
            -1e300, 1e300, // Y dimension limits
            -1e300, 1e300; // Z dimension limits
    bbox.transposeInPlace();
    bbox.resize(1,bbox.rows()*bbox.cols());
    participant.setMeshAccessRegion(GeometryControlPointMesh,bbox);

    // Step 1: Collect geometry from PreCICE
    gsVector<index_t> knotIDs;
    gsMatrix<> knots;
    participant.getMeshVertexIDsAndCoordinates(GeometryKnotMesh,knotIDs,knots);

    gsBasis<> * basis = knotMatrixToBasis<real_t>(knots).get();

    gsVector<index_t> controlPointIDs;
    gsMatrix<> controlPoints;
    participant.getMeshVertexIDsAndCoordinates(GeometryControlPointMesh,controlPointIDs,controlPoints);


    // Step 2: Initialize the force mesh
    gsVector<index_t> forceKnotIDs;
    gsMultiPatch<> forceMesh;
    gsKnotVector<> forceKnots(0,1,10,3,1); //begin, end, # interiors, mult, by default =1
    gsTensorBSplineBasis<2, real_t> forceBasis(forceKnots,forceKnots);

    gsMatrix<> forceControlPoints(forceBasis.size(), 1);
    forceControlPoints.setConstant(-1e4);

    forceMesh.addPatch(give(forceBasis.makeGeometry(forceControlPoints)));



    // gsBoundaryConditions<> bcInfo;


    //Bottom side (fixed)
    // bcInfo.addCondition(0, boundary::south, condition_type::dirichlet, nullptr, -1);
    // bcInfo.addCondition(0, boundary::south, condition_type::clamped, nullptr, 2);

    // gsFunctionExpr<> surfForce("0", "0", "-1e4",3);

    // bcInfo.addCondition(0, boundary::north, condition_type::dirichlet, &surfForce); 

    // bcInfo.setGeoMap(forceMesh);

    // Step 4a: write force knot mesh to PreCICE 
    gsMultiBasis<> bases(forceMesh);
    gsMatrix<> forceKnots = knotsToMatrix(bases.basis(0));
    participant.addMesh(ForceKnotMesh, forceKnots);

    // Step 4b: write force control point mesh to PreCICE
    gsVector<index_t> ForceControlPointIDs;
    gsMatrix<> forceControlPoints = forceMesh.patch(0).coefs().transpose();
    participant.addMesh(ForceControlPointMesh,forceControlPoints, ForceControlPointIDs);




    // The force mesh has been added to the participant
    real_t precice_dt = participant.initialize();



    // gsVector<index_t> knotIDs;
    // gsMatrix<> knots;
    // participant.getMeshVertexIDsAndCoordinates(KnotMesh,knotIDs,knots);

    // gsBasis<> * basis = knotMatrixToBasis<real_t>(knots).get();

    // gsVector<index_t> controlPointIDs;
    // gsMatrix<> controlPoints;
    // participant.getMeshVertexIDsAndCoordinates(ControlPointMesh,controlPointIDs,controlPoints);

    // gsVector<index_t> forceControlPointIDs;
    // gsMatrix<> ForceControlPoints;
    // partipant.getMeshVertexIDsAndCoordinates(KnotMesh,forceControlPointIDs,ForceControlPoints);


    gsMultiPatch<> mp, deformation;


    // mp.addPatch(give(basis->makeGeometry(CPs.transpose())));

    mp.addPatch(give(basis->makeGeometry(controlPoints.transpose())));
    deformation = mp;
    deformation.patch(0).coefs().setZero();




    // gsVector<index_t> quadPointIDs;
    // gsMatrix<> quadPoints;
    // participant.getMeshVertexIDsAndCoordinates(ForceMesh,quadPointIDs,quadPoints);
    // gsMatrix<> quadPointData(mp.geoDim(),quadPoints.cols());

    // real_t t = 0, dt = precice_dt;
    // index_t timestep = 0;
    // // Define the solution collection for Paraview
    // gsParaviewCollection collection("./output/solution");
    // // Time integration loop
    // while (participant.isCouplingOngoing())
    // {
    //     if (participant.requiresWritingCheckpoint())
    //         gsInfo<<"Writing Checkpoint\n";

    //     // Read control point displacements
    //     participant.readData(ControlPointMesh,ControlPointData,controlPointIDs,controlPoints);
    //     deformation.patch(0).coefs() = controlPoints.transpose();

    //     // Write data at the quadrature points
    //     quadPointData.setZero();
    //     quadPointData.row(2).setConstant(-1e4);
    //     participant.writeData(ForceMesh,ForceData,quadPointIDs,quadPointData);

    //     // do the coupling
    //     precice_dt =participant.advance(dt);

    //     dt = std::min(precice_dt, dt);

    //     if (participant.requiresReadingCheckpoint())
    //         gsInfo<<"Reading Checkpoint\n";
    //     else
    //     {
    //         t += dt;
    //         timestep++;

    //         gsField<> solution(mp,deformation);
    //         if (timestep % plotmod==0 && plot)
    //         {
    //             // solution.patch(0).coefs() -= patches.patch(0).coefs();// assuming 1 patch here
    //             std::string fileName = "./output/solution" + util::to_string(timestep);
    //             gsWriteParaview<>(solution, fileName, 500);
    //             fileName = "solution" + util::to_string(timestep) + "0";
    //             collection.addTimestep(fileName,t,".vts");
    //         }

    //         // solution.patch(0).eval_into(points,pointDataMatrix);
    //         // otherDataMatrix<<time;
    //         // writer.add(pointDataMatrix,otherDataMatrix);
    //     }
    // }

    if (plot)
    {
        // collection.save();
    }


    return  EXIT_SUCCESS;
}
