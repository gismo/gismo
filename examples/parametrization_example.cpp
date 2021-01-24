/** @file parametrization_example.cpp

    @brief Tutorial on how to use gsParametrization

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl, D. Mokris
*/

#include <gismo.h>

#include "gsModeling/gsParametrization.h"
#include "gsModeling/gsPeriodicParametrizationOverlap.h"
#include "gsModeling/gsPeriodicParametrizationStitch.h"
//#include "gsModeling/gsFreeBoundaryParametrization.h"
//#include "gsModeling/gsIterativeParametrization.h"

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
    std::string filenameV0, filenameV1, fileCorners, filenameOverlap, filenameStitch;
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
    cmd.addString("d", "v0", ".xml file containing points and u-parameters of points with v=0", filenameV0);
    cmd.addString("t", "v1", ".xml file containing points and u-parameters of points with v=1", filenameV1);
    cmd.addString("x", "fileCorners", "file with 3D coordinates of the corners", fileCorners);
    
    cmd.addString("l", "overlap", ".stl file of the overlap for periodicity; must not be combined with -s", filenameOverlap);
    cmd.addString("s", "stitch", ".xml file with the vertices on the stitch for periodicity; must not be combined with -l.", filenameStitch);
    cmd.addReal("r", "range", "in case of restrict or opposite", range);
    cmd.addInt("n", "number", "number of corners, in case of corners", number);
    cmd.addMultiInt("c", "corners", "vector for corners, call it every time for an entry (-c 3 -c 1 -c 2 => {3,1,2})", corners);
    cmd.addSwitch("plot","Plot with Paraview",paraview);
    cmd.addSwitch("fit", "Create an .xml file suitable for surface fitting with G+Smo.", fitting);
    cmd.getValues(argc, argv);

    gsOptionList ol = cmd.getOptionList();

    gsFileData<> fd(ol.getString("filenameIn"));

    gsInfo << "Reading input into gsMesh<real_t>::uPtr: ";
    gsStopwatch stopwatch;
    gsMesh<real_t>::uPtr mm = fd.getFirst<gsMesh<real_t> >();
    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    enum domainOptions {standard, overlap, stitch, iterative, free} domainMethod;

    if( ol.askString("overlap", "").compare("") > 0 )
	domainMethod = overlap;
    else if( ol.askString("stitch", "").compare("") > 0 )
	domainMethod = stitch;
    // else if( ol.askString("fileCorners", "").compare("") > 0)
    // 	domainMethod = free;
    // else
    // 	domainMethod = iterative;
    else
	domainMethod = standard;

    gsInfo << "DomainMethod set to " << domainMethod << "." << std::endl;

    gsInfo << "creating gsParametrization<real_t>       ";
    stopwatch.restart();

    gsPeriodicParametrizationOverlap<real_t>* pm_over;
    gsPeriodicParametrizationStitch<real_t>* pm_stitch;
    //gsFreeBoundaryParametrization<real_t>* pm_free;
    //gsIterativeParametrization<real_t>*      pm_iter;
    gsParametrization<real_t>* pm;

    if(domainMethod == overlap)
    {
	pm_over = new gsPeriodicParametrizationOverlap<real_t>(*mm, ol);
	pm = pm_over;
    }
    else if(domainMethod == stitch)
    {
	pm_stitch = new gsPeriodicParametrizationStitch<real_t>(*mm, ol);
	pm = pm_stitch;
    }
    // else if(domainMethod == free)
    // {
    // 	pm_free = new gsFreeBoundaryParametrization<real_t>(*mm, ol);
    // 	pm = pm_free;
    // }
    // else if(domainMethod == iterative)
    // {
    // 	pm_iter = new gsIterativeParametrization<real_t>(*mm, ol);
    // 	pm = pm_iter;
    // }
    else
    {
	pm = new gsParametrization<real_t>(*mm, ol);
    }

    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    pm->setOptions(ol);

    gsInfo << "gsParametrization::compute()             ";
    stopwatch.restart();

    if( domainMethod == overlap )
	pm_over->compute_periodic_overlap(filenameV0, filenameV1, filenameOverlap);
    else if( domainMethod == stitch )
    {
	gsMatrix<real_t> verticesV0, paramsV0, verticesV1, paramsV1, stitchVertices;
	readParsAndPts(filenameV0, paramsV0, verticesV0);
	readParsAndPts(filenameV1, paramsV1, verticesV1);
	readPts(filenameStitch,    stitchVertices);

	pm_stitch->compute(verticesV0, paramsV0, verticesV1, paramsV1, stitchVertices);
    }
    // else if(domainMethod == free)
    // 	pm_free->compute_free_boundary();
    // else if(domainMethod == iterative)
    // 	pm_iter->compute();
    else
	pm->compute();

    stopwatch.stop();
    gsInfo << stopwatch << "\n";

    gsInfo << "gsParametrization::createFlatMesh()      ";
    gsMesh<> flatMesh;

    stopwatch.restart();
    if( domainMethod == overlap )
	flatMesh = pm_over->createFlatMesh(true);
    else if( domainMethod == stitch )
	flatMesh = pm_stitch->createFlatMesh(true);
    else
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

    if( domainMethod == overlap)
	pm_over->restrictMatrices(uv, xyz);
    else if( domainMethod == stitch)
	pm_stitch->restrictMatrices(uv, xyz);

    if(paraview)
    {
	gsInfo << "Writing to Paraview." << std::endl;

	// .pvd with the flat mesh
	gsWriteParaview(flatMesh, ol.getString("filenameOut"));

	// .vtk with the vertices coloured according to the parameters
	// Note: calling gsWriteParaview directly with the uv matrix
	// would not do, as the vertices are in different order than
	pm->writeTexturedMesh(ol.getString("filenameOut"));
	pm->writeSTL(*mm, ol.getString("filenameOut"));
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
	output.save(ol.getString("filenameOut"));
    }

    delete pm;

    return 0;
}
