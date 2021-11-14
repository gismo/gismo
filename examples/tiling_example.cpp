/** @file tiling_example.cpp

    @brief Tutorial using tiles.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <iostream>
#include <gismo.h>
#include <gsCore/gsGeometry.h>
#include <gsNurbs/gsNurbsCreator.h>

using namespace gismo;

// template <typename T>
// gsMultiPatch<T> makeGrid(gsGeometry<T> & g, const index_t m = 0.0, const index_t n = 0.0)
// {
//     return makeGrid(g.clone(),m,n);
// }


int main(int argc, char* argv[])
{

    std::string input("surfaces/simple.xml");
    std::string plotname("tiles");
    std::string writename("tiles");

    bool alternate = false;
    bool plot = false;
    bool write = false;

    index_t N = 6;
    index_t M = 4;

    index_t numElevate = 0;
    index_t numHRef = 0;

    gsCmdLine cmd("Tutorial on gsGeometry class.");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addInt("N","numX","Number in horizontal direction",N);
    cmd.addInt("M","numY","Number in vertical direction",M);
    cmd.addInt("r","numRef","Number of refinements of original tile",numHRef);
    cmd.addInt("e","numEl","Number of elevations of original tile",numElevate);
    cmd.addSwitch("plot", "Plot the geometry to Paraview", plot);
    cmd.addSwitch("write", "Write the geometry to XML", write);
    cmd.addSwitch("alternate", "Alternate tiles (mirroring)", alternate);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======================================================================
    // reading the geometry
    // ======================================================================

    gsFileData<> fileData(input);



    // ======================================================================
    // writing to paraview
    // ======================================================================

    gsMultiPatch<> mp, mp_copy;
    std::vector<gsMultiPatch<real_t>> container(2);

    if (input.empty())
    {
        mp.clear();
        mp.addPatch(gsNurbsCreator<>::BSplineTrapezium(1,0.5,1.0,0.0,0.0));
        mp.addPatch(gsNurbsCreator<>::BSplineTrapezium(1,0.5,1.0,0.0,90.));
        mp.addPatch(gsNurbsCreator<>::BSplineTrapezium(1,0.5,1.0,0.0,180));
        mp.addPatch(gsNurbsCreator<>::BSplineTrapezium(1,0.5,1.0,0.0,-90));

        // mp.clear();
        // mp.addPatch(gsNurbsCreator<>::BSplineTrapezium(1,0.5,1.0,0.0,0.0));
        // mp.addPatch(gsNurbsCreator<>::BSplineTrapezium(1,0.5,1.0,0.0,90.));
        // mp.addPatch(gsNurbsCreator<>::BSplineTrapezium(1,0.5,1.0,0.0,180));
        // mp.addPatch(gsNurbsCreator<>::BSplineTrapezium(1,0.5,1.0,0.0,-90));
        // mp_copy = mp;

        // mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,1,1,45));
        // gsWriteParaview(mp,output + "_a",1000,true,true);

        // mp.clear();
        // mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,1,1,0));
        // rotate2D(mp.patch(0),45.0,0.0,0.0);
        // gsWriteParaview(mp,output + "_b",1000,true,true);


        // mp.clear();
        // mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,1,1,0));
        // mp.addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,1,1,0));
        // shift2D(mp.patch(1),1.0,0.);
        // shift2D(mp,1.0,1.0);
        // makeGrid(mp,3,3);
        // gsWriteParaview(mp,output + "_c",1000,false);

        alternate = true;
    }
    else
        gsReadFile<>(input,mp);

    for (index_t k = 0; k!=numHRef; k++)
        mp.uniformRefine();

    for (index_t k = 0; k!=numElevate; k++)
        mp.degreeElevate();

    gsInfo<<"MultiPatch:\n"<<mp;
    for (size_t p=0; p!=mp.nPatches(); p++)
        gsInfo<<"\t - Patch "<<p<<":\n"<<mp.patch(p)<<"\n";

    if (alternate)
    {
        mp_copy = mp;
        mp.clear();
        mp = mp_copy;
        gsNurbsCreator<>::mirror2D(mp,1);
        container[0] = mp_copy;
        container[1] = mp;
        mp = gsNurbsCreator<>::makeGrid(container,N,M);
    }
    else
        gsNurbsCreator<>::makeGrid(mp,N,M);

    mp.computeTopology();

    gsInfo<<"Tiling finished";

    gsMatrix<> bbox, pbbox;
    mp.boundingBox(bbox);
    // std::vector<index_t> boundaryPatches(mp.nPatches());
    for (size_t p = 0; p!= mp.nPatches(); p++)
    {
        gsMultiPatch<> mp_tmp(mp.patch(p));
        mp_tmp.boundingBox(pbbox);
        if ( (bbox(0,0) - pbbox(0,0)) ==0 )
            gsInfo<<"Patch "<<p<<" is a boundary patch on west side!\n";
        if ( (bbox(0,1) - pbbox(0,1)) ==0 )
            gsInfo<<"Patch "<<p<<" is a boundary patch on east side!\n";
        if ( (bbox(1,0) - pbbox(1,0)) ==0 )
            gsInfo<<"Patch "<<p<<" is a boundary patch on south side!\n";
        if ( (bbox(1,1) - pbbox(1,1)) ==0 )
            gsInfo<<"Patch "<<p<<" is a boundary patch on north side!\n";
    }

    if (plot)
        gsWriteParaview(mp,plotname,10,true);

    if (write)
        gsWrite(mp,writename);





    return 0;
}


