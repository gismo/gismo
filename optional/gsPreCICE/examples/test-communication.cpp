/** @file heat-equation-coupling.cpp

    @brief Heat equation participant for a double coupled heat equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsPreCICE/gsPreCICE.h>
#include <gsPreCICE/gsPreCICEFunction.h>
#include <random>


using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    std::string participantName = "SolverOne";
    std::string precice_config("../precice_config.xml");

    gsCmdLine cmd("Coupled heat equation using PreCICE.");
    cmd.addString("n","name","name of this participant",participantName);
    cmd.addString( "c", "config", "PreCICE config file", precice_config );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read input file]

    std::string otherName;
    if (participantName == "SolverOne")
        otherName = "SolverTwo";
    else if (participantName == "SolverTwo")
        otherName = "SolverOne";
    else
        GISMO_ERROR("Invalid names");

    std::string participantDataName = participantName + "Data";
    std::string otherDataName = otherName + "Data";

    std::string participantMeshName = participantName + "Mesh";
    std::string otherMeshName = otherName + "Mesh";

    // Set random seed based on name
    std::srand(std::hash<std::string>{}(participantName));
    gsMatrix<> points;
    points.setRandom(2,10);

    gsPreCICE<real_t> interface(participantName, precice_config);
    interface.addMesh(participantMeshName,points);

    if (interface.requiresInitialData())
        gsDebugVar("Hi");

    gsMatrix<> bbox(2,2);
    // bbox.row(0)<<points.row(0).minCoeff(),points.row(0).maxCoeff();
    // bbox.row(1)<<points.row(1).minCoeff(),points.row(1).maxCoeff();
    bbox.row(0)<<-1,1;
    bbox.row(1)<<-1,1;
    bbox.transposeInPlace();
    bbox.resize(1,bbox.rows()*bbox.cols());
    gsDebugVar(bbox);

    interface.setMeshAccessRegion(otherMeshName,bbox);
    real_t precice_dt = interface.initialize();


    std::vector<index_t> writeIDs;
    gsMatrix<> writePoints;
    interface.getMeshVertexIDsAndCoordinates(otherMeshName,writeIDs,writePoints);

    gsInfo<<"My mesh is:\n"
          <<points
          <<"\n";

    gsInfo<<"Their mesh is:\n"
          <<writePoints
          <<"\n";


    gsMatrix<> writeData;
    writeData.setRandom(1,writeIDs.size());
    interface.writeData(otherMeshName,otherDataName,writeIDs,writeData);
    gsInfo<<"Wrote the following data ("<<participantName<<" -> "<<otherName<<"):\n"
          <<writeData
          <<"\n";

    interface.advance(precice_dt);

    gsMatrix<> readData;
    interface.readData(participantMeshName,participantDataName,points,readData);



    gsInfo<<"Read the following data ("<<otherName<<" -> "<<participantName<<"):\n"
          <<readData
          <<"\n";

        // interface.writeData(writeMeshName,writeName=="Dirichlet" ? fluxName : tempName,writeIDs,result);
    return  EXIT_SUCCESS;
}
