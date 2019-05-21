/** @file gsPatchFromBoundary.cpp

    @brief Constructs a patch (surface, volume,..) from a set of
    boundaries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

using namespace gismo;


int main(int argc, char* argv[])
{
    bool save      = false;
    index_t method = 0;
    real_t tol = 1e-4;
    std::string fn = "curves3d/curve_boundary.xml";

    // Read input from command line arguments
    gsCmdLine cmd("Constructs a patch given a domain boundary.");
    cmd.addPlainString("filename", "File containing boundary data", fn);
    cmd.addInt("m","method" ,"Method: 0 Coons' patch (default), 1 Spring patch, 2: Cross-Ap. patch", method);
    cmd.addReal  ("t","tolerance","Tolerance for identifing patch interfaces", tol);
    cmd.addSwitch("save", "Save result in XML format", save);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Load XML file
    gsMultiPatch<> boundary;
    gsReadFile<>(fn, boundary);
    GISMO_ENSURE(!boundary.empty(), "The gsMultiPatch is empty - maybe file is missing or corrupt.");
    gsInfo<<"Got "<< boundary <<"\n";
    boundary.computeTopology(tol);
    GISMO_ENSURE( boundary.isClosed(), "The boundary is not closed, adjust tolerance.");
    boundary.closeGaps(tol);

    switch (method)
    {
    case 1:
    {
        gsInfo<<"Using spring patch construction.\n";
        gsSpringPatch<real_t> spring(boundary);
        gsInfo<<"Created a " << spring.compute() <<"\n";
        if (save) gsWrite(spring.result(), "result_patch");
        break;
    }
    case 2:
    {
        gsInfo<<"Using cross approximation construction.\n";
        gsCrossApPatch<real_t> cross(boundary);
        gsInfo<<"Created a " << cross.compute() <<"\n";
        if (save) gsWrite(cross.result(), "result_patch");
        break;
    }
    case 0:
    default:
        gsInfo<<"Using Coons' patch construction.\n";
        gsCoonsPatch<real_t> coons(boundary);
        gsInfo<<"Created a " << coons.compute() <<"\n";
        if (save) gsWrite(coons.result(), "result_patch");
        break;
    }

    if (save)
        gsInfo << "Result saved to result_patch.xml\n";
    else
        gsInfo << "Done. No output created, re-run with --save to get xml "
                  "file containing the data.\n";
    return 0;
}
