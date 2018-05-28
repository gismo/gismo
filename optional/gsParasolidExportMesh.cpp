/** @file gsParasolidExportMesh.cpp

    @brief Test for exporting gsMesh to the Siemens's Parasolid geometric kernel.

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

    gsCmdLine cmd("Ploting a geometry to paradolid wire geometry.");
    cmd.addString("i", "input", "Input file", input);
    cmd.addString("o", "output", "Output file", output);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    std::cout << " \n\nInput arguments: \n\n"
	      << "input: " << input << "\n\n"
	      << "output: " << output << "\n\n"
	      << "--------------------------------------------------\n" << std::endl;
    
    gsGeometry<>::uPtr geom = gsReadFile<>(input);
    
    gsMesh<> mesh;
    
    makeMesh<>( geom->basis(), mesh, 5);
    
    geom->evaluateMesh(mesh);
    
    extensions::gsWriteParasolid(mesh, output);
    
    return 0;
}

    
