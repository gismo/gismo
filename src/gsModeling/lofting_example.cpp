/** @file lofting_example.cpp

    @brief Demonstrates adaptive fitting of data samples

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Options with default values
    
    std::string fn = "curves3d/loft_exm0.xml";
    real_t lambda = 1e-07;

    // Reading options from the command line
    gsCmdLine cmd("Loft parametrized B-spline cuves. Expected input file is an XML "
            "file containing the curves to be lofted.");
    
    cmd.addString("f", "filedata", "Input sample data", fn);
    cmd.addReal("s", "lambda", "smoothing coefficient", lambda);
    

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    
    std::vector<gsBSpline<real_t>> container;
    //! [Read data]
    gsFileData<> fd_in(fn);
    
    // This is for outputing an XML file, if requested
    gsFileData<> fd;


    // Create knot-vectors without interior knots
    // gsKnotVector<> u_knots (u_min, u_max, 0, deg_x+1 ) ;
    // gsKnotVector<> v_knots (v_min, v_max, 0, deg_y+1 ) ;


    // Create hierarchical refinement object
    // gsLofting<real_t> loft( container );


    return 0;
}
