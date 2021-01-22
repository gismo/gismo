/** @file poisson2_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    std::string save;
    real_t tol = 1e-4;
    bool verbose = false;
    int startIndex = 0;

    gsCmdLine cmd("Adds topology to multipatch file.");
    cmd.addReal("T","tol", "tolerance", tol);
    cmd.addString("f","file", "File from which the geometries must be read", save);
    cmd.addSwitch("verbose","Verbose output", verbose);
    cmd.addInt("S", "startIndex", "Index From which the geometry Indeces must start", startIndex);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    gsInfo << "\n";

    if (verbose) {
        if (save.empty()) {
            gsInfo << "gsError: \t" << "No file path specified";
        }
    }
    GISMO_ENSURE(!save.empty(), "No file path specified.");

    gsFileData<> fd(save);

    if (verbose) {
        if (!fd.template has<gsMultiPatch<real_t>>()) {
            gsInfo << "gsError: \t" << "File does not have a Multipatch object!";
        }
    }
    GISMO_ENSURE(fd.template has<gsMultiPatch<real_t>>(),"File does not have a Multipatch object!");

    if (verbose)
        gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    std::vector< memory::unique_ptr<gsMultiPatch<real_t>> > vec = fd.template getAll<gsMultiPatch<real_t>>();

    index_t num = fd.template count<gsMultiPatch<real_t>>(); // id=0: Multipatch domain
    GISMO_ENSURE(num==1,"Number of multipatch objects in XML should be 1, but is "<<num);

    gsMultiPatch<> mp;
    fd.template getFirst<gsMultiPatch<real_t>>(mp); // Multipatch domain
    mp.clearTopology();
    mp.computeTopology(tol);


    std::string TopologyData;

    mp.topology().printLikeXML(std::cout, startIndex);


    if (verbose)
    {
        gsInfo<<"Writing geometry to XML...\n";
        gsInfo<<"Path = "<<save<<"\n";
    }

    //gsWrite(mp,save);
    return EXIT_SUCCESS;

}// end main
