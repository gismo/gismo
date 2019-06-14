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
    bool paraview = false;
    index_t parametrizationMethod(1); // 1:shape, 2:uniform, 3:distance
    // shape: best method, shape of the mesh is preserved, smooth surface fitting
    // uniform: the lambdas according to Floater's algorithm are set to 1/d, where d is the number of neighbours
    // distance: the lambdas according to Floater's algorithm are set to the relative distances between the point and its neighbours
    index_t boundaryMethod(4); // 1:chords, 2:corners, 3:smallest, 4:restrict, 5:opposite, 6:distributed
    //chords: choose boundary points distributed on the unit square boundary wrt the chord lengths, no input needed
    //corners: choose 4 boundary corners for the corner points of the unit square and the rest of the boundary corners is distributed on the four edges of the unit square wrt the chord lengths, input is the corner numbers, e.g. 1,2,3,4 for the first, second, third and fourth boundary point
    //smallest: choose 4 boundary corners as the vertices with the smallest inner angles, no input needed
    //restrict: choose first boundary corner as the one with the smallest inner angle, then restrict an environment of this point (range) and search for the second smallest boundary corner point in all the others except the restricted ones, etc., input range e.g. r=0.1 for 1/10 of whole boundary length is restricted around the already chosen corner point
    // opposite: choose the boundary corner points such that they are nearly opposite of each other, input range e.g. r=0.1 for 1/10 of whole boundary length around exact opposite point on boundary is possible
    // distributed: choose the smallest inner angle corners (number for how much to choose) and choose four corners s.t. they are as evenly distributed as possible, input number n=6 for choosing 6 boundary vertices with smallest inner angles and then find 4 of them s.t. evenly distributed
    std::string filenameIn("stl/norm.stl");
    std::string filenameOut("flatMesh");
    real_t range = 0.1; // in case of restrict or opposite
    index_t number = 4; // number of corners, in case of distributed
    std::vector<index_t> corners; // in case of corners

    gsCmdLine cmd("parametrization_example Command line");
    cmd.addInt("m", "parametrizationMethod", "parametrization methods: {1: shape, 2: uniform, 3: distance}\n"
                                                "shape: best method, shape of the mesh is preserved, smooth surface fitting\n"
                                                "uniform: the lambdas according to Floater's algorithm are set to 1/d, where d is the number of neighbours\n"
                                                "distance: the lambdas according to Floater's algorithm are set to the relative distances between the point and its neighbours",
                  parametrizationMethod);
    cmd.addInt("b", "boundaryMethod", "boundary methodes: {1: chords, 2: corners, 3: smallest, 4: restrict, 5: opposite, 6: distributed}\n"
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
    cmd.addSwitch("plot","Plot with Paraview",paraview);
    cmd.getValues(argc, argv);

    gsOptionList ol = cmd.getOptionList();

    gsFileData<> fd(ol.getString("filenameIn"));

    gsInfo << "Reading input into gsMesh<real_t>::uPtr: ";
    gsStopwatch stopwatch;
    gsMesh<real_t>::uPtr mm = fd.getFirst<gsMesh<real_t> >();
    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    gsInfo << "creating gsParametrization<real_t>       ";
    stopwatch.restart();
    gsParametrization<real_t> pm(*mm);
    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    pm.setOptions(ol);

    gsInfo << "gsParametrization::compute()             ";
    stopwatch.restart();
    pm.compute();
    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    gsInfo << "gsParametrization::createFlatMesh()      ";
    stopwatch.restart();
    gsMesh<> mesh = pm.createFlatMesh();
    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    gsInfo << "gsParametrization::createUVmatrix()      ";
    stopwatch.restart();
    gsMatrix<> uv = pm.createUVmatrix();
    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    gsInfo << "gsParametrization::createXYZmatrix()     ";
    stopwatch.restart();
    gsMatrix<> xyz = pm.createXYZmatrix();
    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    gsInfo << mesh << "\n";
    gsInfo << uv << "\n";
    gsInfo << xyz << "\n";

    if(paraview)
    {
        gsWriteParaview(mesh, ol.getString("filenameOut"));
        gsFileManager::open(ol.getString("filenameOut") + ".pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return 0;
}
