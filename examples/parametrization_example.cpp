/** @file parametrization_example.cpp

    @brief Tutorial on how to use gsParametrization

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    std::string parametrizationMethod("shape"); // shape, uniform, distance
    // shape: best method, shape of the mesh is preserved, smooth surface fitting
    // uniform: the lambdas according to floater's algorithm are set to 1/d, where d is the number of neighbours
    // distance: the lambdas according to floater's algorithm are set to the relative distances between the point and its neighbours
    std::string boundaryMethod("restrict"); // chords, corners, smallest, restrict, opposite, distributed
    //chords: choose boundary points distributed on the unit square boundary wrt the chord lengths, no input needed
    //corners: choose 4 boundary corners for the corner points of the unit square and the rest of the boundary corners is distributed on the four edges of the unit square wrt the chord lengths, input is the corner numbers, e.g. 1,2,3,4 for the first, second, third and fourth boundary point
    //smallest: choose 4 boundary corners as the vertices with the smallest inner angles, no input needed
    //restrict: choose first boundary corner as the one with the smallest inner angle, then restrict an environment of this point (range) and search for the second smallest boundary corner point in all the others except the restricted ones, etc., input range e.g. r=0.1 for 1/10 of whole boundary length is restricted around the already chosen corner point
    // opposite: choose the boundary corner points such that they are nearly opposite of each other, input range e.g. r=0.1 for 1/10 of whole boundary length around exact opposite point on boundary is possible
    // distributed: choose the smallest inner angle corners (number for how much to choose) and choose four corners s.t. they are as evenly distributed as possible, input number n=6 for choosing 6 boundary vertices with smallest inner angles and then find 4 of them s.t. evenly distributed
    std::string filenameIn("/home/adamanta/vogl/SVN/gismoParam/examples/R60_01_Segel_sheet_closed.stl");
    std::string filenameOut("flatMesh");
    real_t range = 0.1; // in case of restrict or opposite
    int number = 4; // number of corners, in case of distributed
    std::vector<int> corners; // in case of corners

    gsCmdLine cmd("parametrization_example Command line");
    cmd.addString("m", "parametrizationMethod", "parametrization methods: {shape, uniform, distance}\n"
                                                "shape: best method, shape of the mesh is preserved, smooth surface fitting\n"
                                                "uniform: the lambdas according to floater's algorithm are set to 1/d, where d is the number of neighbours\n"
                                                "distance: the lambdas according to floater's algorithm are set to the relative distances between the point and its neighbours",
                  parametrizationMethod);
    cmd.addString("b", "boundaryMethod", "boundary methodes: {chords, corners, smallest, restrict, opposite, distributed}\n"
                                         "chords: choose boundary points distributed on the unit square boundary wrt the chord lengths, no input needed\n"
                                         "corners: choose 4 boundary corners for the corner points of the unit square and the rest of the boundary corners is distributed on the four edges of the unit square wrt the chord lengths, input is the corner numbers, e.g. 1,2,3,4 for the first, second, third and fourth boundary point"
                                         "smallest: choose 4 boundary corners as the vertices with the smallest inner angles, no input needed"
                                         "restrict: choose first boundary corner as the one with the smallest inner angle, then restrict an environment of this point (range) and search for the second smallest boundary corner point in all the others except the restricted ones, etc., input range e.g. r=0.1 for 1/10 of whole boundary length is restricted around the already chosen corner point"
                                         "opposite: choose the boundary corner points such that they are nearly opposite of each other, input range e.g. r=0.1 for 1/10 of whole boundary length around exact opposite point on boundary is possible"
                                         "distributed: choose the smallest inner angle corners (number for how much to choose) and choose four corners s.t. they are as evenly distributed as possible, input number n=6 for choosing 6 boundary vertices with smallest inner angles and then find 4 of them s.t. evenly distributed",
                  boundaryMethod);
    cmd.addString("f", "filenameIn", "input file name", filenameIn);
    cmd.addString("o", "filenameOut", "output file name", filenameOut);
    cmd.addReal("r", "range", "in case of restrict or opposite", range);
    cmd.addInt("n", "number", "number of corners, in case of corners", number);
    cmd.addMultiInt("c", "corners", "vector for corners, call it every time for an entry (-c 3 -c 1 -c 2 => {3,1,2})", corners);
    cmd.getValues(argc, argv);

    gsOptionList ol = cmd.getOptionList();

    gsFileData<> fd(ol.getString("filenameIn"));
    gsMesh<real_t>::uPtr mm = fd.getFirst<gsMesh<real_t> >();

    gsParametrization<real_t> pm(*mm, ol);

    gsMesh<> mesh = pm.createFlatMesh();
    gsMatrix<> uv = pm.createUVmatrix();
    gsMatrix<> xyz = pm.createXYZmatrix();

    gsInfo << mesh << "\n";
    gsInfo << uv << "\n";
    gsInfo << xyz << "\n";

    gsWriteParaview(mesh, ol.getString("filenameOut"));
    gsFileManager::open(ol.getString("filenameOut") + ".pvd");

    return 0;
}