/** @file embeddedCurve_example.cpp

    @brief Example for embedded curve

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Chianese, H.M. Verhelst
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;

template <class T>
T curveCoordinate(const gsGeometry<T> & geometry, const T & value, const short_t & dir);

// template <class T>
// T findRoot();

int main(int argc, char *argv[])
{
    // bool plot = false; // If set to true, paraview file is generated and launched on exit
    // bool trim = false; // If set to true, trim/merge operations are displayed
    // bool intersect = false; // If set to true, intersection example is displayed

    gsCmdLine cmd("TODO");
    // cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    // cmd.addSwitch("trim", "Basic trim/merge operations", trim);
    // cmd.addSwitch("intersect", "Intersection operations", intersect);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    // Make a BSpline curve
    gsKnotVector<> kv(0, 1, 1, 3);//start,end,interior knots, start/end multiplicites of knots1
    gsMatrix<> coefs(4, 3);
    coefs << 0, 0, 0,
             1, 2, 3,
             2, 1, 4,
             4, 4, 4;

    gsBSpline<> curve( kv, give(coefs));

    // Print the Bspline curve
    gsInfo << "I am a " << curve << "\n";


    // gsInfo<<curveCoordinate(coefs);

    return EXIT_SUCCESS;
}

template <class T>
T curveCoordinate(const gsGeometry<T> & geometry, const T & value, const short_t & dir)
{
    GISMO_ASSERT(geometry.parDim() == 1, "The geometry must be a curve with parameter dimension 1");
    GISMO_ASSERT(dir < geometry.targetDim(), "The direction must be less than the target dimension of the geometry");
    gsMatrix<T> u(1,1);
    u(0,0) = value;
    return geometry.eval(u)(dir,0);
    // return geometry.deriv(u)(dir,0);
};