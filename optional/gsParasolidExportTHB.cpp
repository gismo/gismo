/** @file gsParasolidExportMesh.cpp

    @brief Test for exporting a THB-spline surface to the Parasolid geometric kernel via Gabor's algorithm.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <iostream>
#include <string>

#include <gismo.h>
#include <gsIO/gsIOUtils.h>

#include <gsParasolid/gsWriteParasolid.h>

using namespace gismo;

int main(int argc, char *argv[])
{

    std::string input("surfaces/thbs_face_3levels.xml");
    std::string output("out");

    gsCmdLine cmd("Exporting a THB-spline surface to parasolid.");
    cmd.addString("i", "input", "Input file", input);
    cmd.addString("o", "output", "Output file", output);
    
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    std::cout << "\n\nInput arguments: \n\n"
	      << "input: " << input << "\n\n"
	      << "output: " << output << "\n\n"
	      << "--------------------------------------------------\n" << std::endl;


    memory::unique_ptr< gsTHBSpline<2> > thb = gsReadFile<>(input);
    
    extensions::gsWriteParasolid(*thb, output);
    
    return 0;
}
    
    







