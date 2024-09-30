/** @file parametrization_example.cpp

    @brief Tutorial on how to use gsParametrization

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl, D. Mokris
*/

#include <gismo.h>

#include <gsModeling/gsParametrization.h>
#include <gsModeling/gsPeriodicOverlap.h>
#include <gsModeling/gsPeriodicStitch.h>

using namespace gismo;

template <class T>
void readParsAndPts(const std::string& filename,
                    gsMatrix<T>& pars,
                    gsMatrix<T>& pts)
{
    gsFileData<T> fd(filename);
    // Cf. https://stackoverflow.com/questions/3505713/c-template-compilation-error-expected-primary-expression-before-token
    fd.template getId<gsMatrix<T> >(0, pars);
    fd.template getId<gsMatrix<T> >(1, pts);

    GISMO_ASSERT(pars.cols() == pts.cols(), "The numbers of parameters and points differ.");
}

template <class T>
void readPts(const std::string& filename,
             gsMatrix<T>& pts)
{
    gsFileData<T> fd(filename);
    fd.template getId<gsMatrix<T> >(0, pts);
}

template <class T>
typename gsParametrization<T>::uPtr newPeriodicOverlap(const gsMesh<T>& mm, const std::string& filenameV0, const std::string& filenameV1,
						       const std::string& filenameOverlap, const gsOptionList& ol)
{
    gsMatrix<real_t> verticesV0, paramsV0, verticesV1, paramsV1;
    readParsAndPts(filenameV0, paramsV0, verticesV0);
    readParsAndPts(filenameV1, paramsV1, verticesV1);

    gsFileData<real_t> fd_overlap(filenameOverlap);
    gsMesh<real_t> overlap = *(fd_overlap.getFirst<gsMesh<real_t> >());

    return typename gsPeriodicOverlap<T>::uPtr(new gsPeriodicOverlap<T>(mm, verticesV0, paramsV0, verticesV1, paramsV1, overlap, ol));

    // Note: paramsV0 and paramsV1 go out of scope here; that's why gs::PeriodicParametrization<T>::m_paramsV0 and [...]V1 cannot be const-references.
}

template <class T>
typename gsParametrization<T>::uPtr newPeriodicStitch(const gsMesh<T>& mm, const std::string& filenameV0, const std::string& filenameV1,
						      const std::string& filenameStitch, const gsOptionList& ol)
{
    gsMatrix<real_t> verticesV0, paramsV0, verticesV1, paramsV1;
    readParsAndPts(filenameV0, paramsV0, verticesV0);
    readParsAndPts(filenameV1, paramsV1, verticesV1);

    gsMatrix<real_t> stitchVertices;
    readPts(filenameStitch, stitchVertices);

    return typename gsPeriodicStitch<T>::uPtr(new gsPeriodicStitch<T>(mm, verticesV0, paramsV0, verticesV1, paramsV1, stitchVertices, ol));
}


int main(int argc, char *argv[])
{
    bool paraview = false;
    bool fitting = false;
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
    std::string filenameV0, filenameV1, filenameOverlap, filenameStitch;
    real_t range = 0.1; // in case of restrict or opposite
    index_t number = 4; // number of corners, in case of distributed
    std::vector<index_t> corners; // in case of corners

    gsCmdLine cmd("Hi, this is the parametrization_example.\n"
		  "Usage examples:\n"
		  "0. The standard parametrization method:\n"
		  "> ./bin/parametrization_example\n"
		  "1. Periodic parametrization (overlap method):\n"
		  "> ./bin/parametrization_example -f parametrization/powerplant-mesh.stl -m 1 -o powerplant -d parametrization/powerplant-bottom.xml -t parametrization/powerplant-top.xml -l parametrization/powerplant-overlap.stl --fit --plot\n"
		  "2. Periodic parametrization (stitch method):\n"
		  "> ./bin/parametrization_example -f parametrization/powerplant-mesh.stl -m 1 -o powerplant -d parametrization/powerplant-bottom.xml -t parametrization/powerplant-top.xml -s parametrization/powerplant-stitch.xml --fit --plot\n");
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
    cmd.addString("d", "v0", ".xml file containing points and u-parameters of points with v=0", filenameV0);
    cmd.addString("t", "v1", ".xml file containing points and u-parameters of points with v=1", filenameV1);
    
    cmd.addString("l", "overlap", ".stl file of the overlap for periodicity; must not be combined with -s", filenameOverlap);
    cmd.addString("s", "stitch", ".xml file with the vertices on the stitch for periodicity; must not be combined with -l.", filenameStitch);
    cmd.addReal("r", "range", "in case of restrict or opposite", range);
    cmd.addInt("n", "number", "number of corners, in case of corners", number);
    cmd.addMultiInt("c", "corners", "vector for corners, call it every time for an entry (-c 3 -c 1 -c 2 => {3,1,2})", corners);
    cmd.addSwitch("plot","Plot with Paraview",paraview);
    cmd.addSwitch("fit", "Create an .xml file suitable for surface fitting with G+Smo.", fitting);
    cmd.getValues(argc, argv);

    gsFileData<> fd(cmd.getString("filenameIn"));

    gsInfo << "Reading input into gsMesh<real_t>::uPtr: ";
    gsStopwatch stopwatch;
    gsMesh<real_t>::uPtr mm = fd.getFirst<gsMesh<real_t> >();
    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    gsInfo << "creating gsParametrization<real_t>       ";
    stopwatch.restart();

    gsParametrization<real_t>::uPtr pm;

    if(cmd.askString("overlap", "").compare("") > 0)
        pm = newPeriodicOverlap(*mm, filenameV0, filenameV1, filenameOverlap, cmd);
    else if(cmd.askString("stitch", "").compare("") > 0)
        pm = newPeriodicStitch(*mm, filenameV0, filenameV1, filenameStitch, cmd);
    else
        pm = gsParametrization<real_t>::uPtr(new gsParametrization<real_t>(*mm, cmd));

    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    pm->setOptions(cmd);

    gsInfo << "gsParametrization::compute()             ";
    stopwatch.restart();

    pm->compute();

    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    gsInfo << "gsParametrization::createFlatMesh()      ";
    gsMesh<> flatMesh;

    stopwatch.restart();
    flatMesh = pm->createFlatMesh();

    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    gsInfo << "gsParametrization::createUVmatrix()      ";
    stopwatch.restart();
    gsMatrix<> uv;
    uv=pm->createUVmatrix();
    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    gsInfo << "gsParametrization::createXYZmatrix()     ";
    stopwatch.restart();
    gsMatrix<> xyz;
    xyz = pm->createXYZmatrix();
    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    if((cmd.askString("overlap", "").compare("") > 0) ||
       (cmd.askString("stitch", "").compare("") > 0))
        pm->restrictMatrices(uv, xyz);

    if(paraview)
    {
        gsInfo << "Writing to Paraview.\n";

        // .pvd with the flat mesh
        gsWriteParaview(flatMesh, cmd.getString("filenameOut"));
        gsFileManager::open(cmd.getString("filenameOut") + ".pvd");

        // .vtk with the vertices coloured according to the parameters
        // Note: calling gsWriteParaview directly with the uv matrix
        // would not do, as the vertices are in different order than
        pm->writeTexturedMesh(cmd.getString("filenameOut"));
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    if(fitting)
    {
        gsInfo << "Writing to G+Smo XML." << std::endl;
        gsFileData<> output;
        output << uv;
        output << xyz;
        output.save(cmd.getString("filenameOut"));
    }

    return 0;
}
