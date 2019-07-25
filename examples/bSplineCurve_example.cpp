/** @file bSplineCurve_example.cpp

    @brief Tutorial on gsBSpline class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false; // If set to true, paraview file is generated and launched on exit


    gsCmdLine cmd("Tutorial 01 shows the use of BSpline curves.");
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);

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

    if (plot)
    {
        // Output a paraview file
        gsWriteParaview( curve, "bsplinecurve", 100);
        gsFileManager::open("bsplinecurve.pvd");
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return 0;
}
