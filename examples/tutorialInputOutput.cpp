/** @file tutorialInputOutput.cpp

    @brief Tutorial on how to use input output facilities of G+Smo library.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

// For detailed options about visualization look at gsView.cpp example. 



#include <iostream>
#include <gismo.h>

using namespace gismo;


int main(int argc, char* argv[])
{

    std::string input(GISMO_DATA_DIR "/curves3d/bspline3d_curve_01.xml");
    std::string output("out");
    
    // for information about command line arguments 
    // look at tutorialCommandLineArg.cpp
        
    gsCmdLine cmd("Tutorial Input Output");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addString("o", "output", "Name of the output file", output);
    bool ok = cmd.getValues(argc,argv);

    if (!ok)
    {
        gsWarn << "Error during parsing the command line!";
        return 1;
    }
    
    // reading the geometry
    
    gsFileData<> fileData(input);
    
    gsGeometry<>* pGeom = NULL;
    if (fileData.has< gsGeometry<> >())
    {
        pGeom = fileData.getFirst< gsGeometry<> >();
    }
    else
    {
        gsWarn << "Input file doesn't have a geometry inside.\n";
        return -1;
    }

        
    if (pGeom == NULL)
    {
        gsWarn << "Didn't find any geometry.\n";
        return -1;
    }


    
    // printing the geometry

    gsInfo << "The file contains: \n" << *pGeom << "\n";
    

    // writing a paraview file
            
    const std::string out = output + "Paraview";
    gsWriteParaview(*pGeom, out);
    gsInfo << "Wrote a paraview file: " << out << "\n";
    
    
    // writing a G+Smo .xml file
            
    gsFileData<> fd;
    fd << *pGeom;
    fd.dump(output);
    gsInfo << "Wrote G+Smo file: " << output << "\n";
    

    delete pGeom;
    return 0;
}



