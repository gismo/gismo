/** @file gsParasolidBCURVE.cpp

    @brief Test for the PK_BCURVE from the Siemens's Parasolid geometric kernel.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <iostream>
#include <string>

#include <gismo.h>


#include <gsParasolid/gsWriteParasolid.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    std::string input("curves3d/bspline3d_curve_01.xml");
    std::string output("out");

    gsCmdLine cmd("Ploting a G+Smo B-spline curve to parasolid PK_BCURVE.");
    cmd.addString("i", "input", "Input file", input);
    cmd.addString("o", "output", "Output file", output);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    std::cout << " \n\nInput arguments: \n\n"
	      << "input: " << input << "\n\n"
	      << "output: " << output << "\n\n"
	      << "--------------------------------------------------\n" << std::endl;

    memory::unique_ptr< gsBSpline<> > curve = gsReadFile<>(input);

    extensions::gsWriteParasolid(*curve, output);

    return 0;
}


