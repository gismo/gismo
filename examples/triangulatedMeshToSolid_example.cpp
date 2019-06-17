/** @file triangulatedMeshToSolid_example.cpp

    @brief Recover the features of a triangulated surface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Pauley
*/


#include <gismo.h>
#include <iostream>
#include <fstream>

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot  = false; // If user gives --plot as argument, paraview file is generated and launched on exit
    bool toxml = false;
    bool noSmooth = false;
    bool writePatchNumbers = false;
    std::string filename = "off/mushroom_triangulated.off";
    real_t cutoffAngle = 40.0;
    real_t innerAngle = 15.0;
    real_t patchAreaWeight = 0.2;
    real_t mergeSmallPatches = 2;
    index_t degree = 4;
    index_t interiorPts = 5;
    index_t wEdge = 10;
    index_t wInterior = 1;

    gsCmdLine cmd("Recover the features of a triangulated surface.");
    cmd.addPlainString("filename", "File containing the input mesh", filename);
    cmd.addSwitch("plot", "Output mesh in ParaView format", plot);
    cmd.addSwitch("xml", "Output solid to xml file", toxml);
    cmd.addSwitch("nosmooth", "Do not smooth corners of faces", noSmooth);
    cmd.addSwitch("writenums", "Write patch numbers to patchnumbers.txt", writePatchNumbers);
    cmd.addInt("f", "interiorweight","Surface fitting: weight for interior points.",
               wInterior);
    cmd.addInt("e", "edgeweight","Surface fitting: weight for edges.",
               wEdge);
    cmd.addInt("i", "interior","Interior knot points of each dimension of trimmed surface.",
               interiorPts);
    cmd.addInt("d", "degree","Degree of each dimension of trimmed surface.",
               degree);
    cmd.addReal("p", "paw","Patch area weight.", patchAreaWeight);
    cmd.addReal("m", "msp","Merge small patches.", mergeSmallPatches);
    cmd.addReal("y", "innerangle","Cutoff angle (degrees) for second pass.",
                innerAngle);
    cmd.addReal("c", "cutoff","Cutoff angle (degrees).", cutoffAngle);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    filename = gsFileManager::find(filename);
    if ( filename.empty() )
    {
        gsInfo << "Could not find file.\n";
        return EXIT_SUCCESS;
    }

    // decide on the base filename
    size_t nameStartIdx = filename.rfind('/');
    if ( nameStartIdx == std::string::npos)
        nameStartIdx = 0;
    else
        nameStartIdx+= 1;
    size_t nameEndIdx = filename.rfind('.');
    std::string baseName = filename.substr(nameStartIdx, nameEndIdx - nameStartIdx);
    gsInfo << "--- Settings ---\n";
    gsInfo << "Processing " << baseName << "\n";
    gsInfo << "Cutoff angle is " << cutoffAngle << "\n";
    gsInfo << "Master surfaces have degree " << degree << " and " << interiorPts << " interior knot points.\n";
    gsInfo << "Surface fit: edge weighting is " << wEdge << " and interior point weighting is " << wInterior << ".\n";
    if(noSmooth) gsInfo << "Will NOT smooth corners of faces\n";
    else gsInfo << "WILL smooth corners of faces\n";
    if(plot) gsInfo << "WILL create paraview plot\n";
    else gsInfo << "Will NOT create paraview plot\n";
    gsInfo << "----------------\n";

    gsMesh<>::uPtr m = gsReadFile<>(filename);
    if (m)
        gsInfo<< "Got "<< *m <<"\n";
    else
    {
        gsInfo<< "Problem encountered in file "<<filename<<", quitting." <<"\n";
        return EXIT_SUCCESS;
    }
    std::vector< gsMesh<> * > paraMeshes;// filled inside toSolid
    std::vector< gsMesh<> * > fitMeshes; // filled inside toSolid
    std::vector< gsMesh<> * > patchMeshes; // filled inside toSolid

    std::vector< std::vector< gsVertex<>* > > iPoints;
    std::vector< std::vector< gsVertex<>* > > oPoints;
    std::vector< std::vector< std::vector<gsVertex<>* > > >  innerBdrys;
    std::vector< std::vector< gsVertex<> > > innerBdrysMassP;
    std::vector< std::vector< bool > > oPointsConvexFlag;

    // set up the gsTriMeshToSolid class to perform the feature detection
    gsTriMeshToSolid<> tmts(m.get());
    bool non_manifold, warning_borders;

    // compute the features
    tmts.getFeatures(cutoffAngle, non_manifold, warning_borders);
    m->cleanMesh();

    // give every face a patch number
    tmts.calcPatchNumbers();

    // improve quality by further dividing and merging patches according to more complex rules
    tmts.storeNeighboringFaces();
    tmts.divideAndMergePatches(innerAngle, patchAreaWeight, mergeSmallPatches);

    // recompute patch numbers
    tmts.calcPatchNumbers();

    if(writePatchNumbers)
    {
        // write patch numbers
        gsInfo << "Writing patch numbers...\n";
        std::ofstream pn("patchnumbers.txt");

        for(std::vector<gsMeshElement<>::gsFaceHandle >::const_iterator it(m->faces().begin());it!=m->faces().end();++it)
        {
            pn << (**it).faceIdentity << "\n";
        }
    }

    // get the patches
    tmts.getFaces(iPoints, oPoints, innerBdrys, innerBdrysMassP, oPointsConvexFlag);

    gsSolid<> sl;
    //if you need a higher number of interior points when fitting the surface increase the 8th input variable(now 5)
    tmts.toSolid(sl,iPoints,oPoints,innerBdrys,innerBdrysMassP,oPointsConvexFlag,paraMeshes,fitMeshes,patchMeshes,degree,interiorPts,true,300,true,wEdge,wInterior,1,noSmooth);
    gsInfo<<sl<<'\n';

    if (toxml)
    {
        gsInfo << "Writing xml file..." << "\n";

        gsFileData<> newdata;
        newdata << sl;
        newdata.dump(baseName);

        //write sharpness information to a file:
        // first line: #patches
        // next #patches lines: #outer points of patch i, sharpness of first point ... sharpness of last point of patch i
        // next #patches lines: #holes of patch i
        // next #holes of patch i lines: #points of hole j of patch i, sharpness of first point ... sharpness of last point of hole j of patch i
        std::string outputFilename;
        outputFilename = "corner_data.txt";
        std::ofstream myfile (outputFilename.c_str());
        if (myfile.is_open())
        {
            myfile << oPoints.size()<<"\n";

            for(size_t i=0;i<oPoints.size();i++)
            {
                myfile << oPoints[i].size()<<" ";
                for(size_t j=0;j<oPoints[i].size();j++)
                    myfile<<(oPoints[i][j]->numEdges>2)<<" ";
                myfile << "\n";

            }
            myfile << "\n";

            for(size_t i=0;i<innerBdrys.size();i++)
            {
                myfile << innerBdrys[i].size()<<" ";
                for(size_t j=0;j<innerBdrys[i].size();j++)
                {
                    myfile << innerBdrys[i][j].size()<<" ";
                    for(size_t k=0;k<innerBdrys[i][j].size();k++)
                        myfile<<(innerBdrys[i][j][k]->numEdges>2)<<" ";
                }
                myfile << "\n";
            }

            myfile.close();
        }

    }

    if (plot)
    {
        // Write a paraview file
        gsInfo<<"Writing paraview file..." << "\n";

        gsWriteParaview( *m, "output");
    }
    else
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    }

    // free meshes
    freeAll(fitMeshes);
    freeAll(paraMeshes);
    freeAll(patchMeshes);

    return EXIT_SUCCESS;
}
