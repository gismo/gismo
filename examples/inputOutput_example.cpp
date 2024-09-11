/** @file inputOutput_example.cpp

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

    //! [Parse command line]
    std::string input("curves3d/bspline3d_curve_01.xml");
    std::string output("");
    bool bezier=false; 

    gsCmdLine cmd("Tutorial Input Output");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addString("o", "output", "Name of the output file", output);
    cmd.addSwitch("b","bezier", "Output using Bezier elements", bezier);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]


    std::string ext;
    ext = (bezier) ? ".vtu" : ".vtp";

    //! [Read geometry]
    if (!gsFileManager::fileExists(input))
    {
        gsWarn << "The file cannot be found!\n";
        return EXIT_FAILURE;
    }

    gsInfo << "Read file \"" << input << "\"\n";

    gsFileData<> fileData(input);

    gsGeometry<>::uPtr pGeom;
    if (fileData.has< gsGeometry<> >())
    {
        pGeom = fileData.getFirst< gsGeometry<> >();
    }
    else
    {
        gsWarn << "Input file doesn't have a geometry.\n";
        return EXIT_FAILURE;
    }
    //! [Read geometry]

    //! [Print geometry]
    gsInfo << "The file contains: \n" << *pGeom << "\n";
    //! [Print geometry]

    if ( output.empty() )
    {
        gsInfo << "Call program with option -o <basename> to write data to files\n";
        gsInfo << "<basename>Paraview"+ext+", <basename>Paraview.pvd, <basename>.xml, <basename>ControlPoints.csv\n";
        return EXIT_SUCCESS;
    }

    // writing a ParaView file
    const std::string out = output + "Paraview";
    if (!fileData.has< gsMultiPatch<> >()) 
    {
        if (!bezier)
        {
            gsWriteParaview(*pGeom, out);
        }
        else
        {
            gsMultiPatch<> mPatch;
            mPatch.addPatch(*pGeom);
            gsWriteParaviewBezier(mPatch, out);
        }
        gsInfo << "Wrote ParaView files: " << out << ext << ", " << out << ".pvd\n";
    }
    else
    {
        gsMultiPatch<> mPatch;
        fileData.getFirst(mPatch);
        if (!bezier)
        {
            gsWriteParaview(mPatch, out);
        }
        else
        {
            bool singleFile = (mPatch.nPatches() > 10 );
            gsWriteParaviewBezier(mPatch, out, singleFile);
        }
        gsInfo << "Wrote ParaView file: " << out << ".pvd\n";
    }

    //! [Write geometry]
    // writing a G+Smo .xml file
    gsFileData<> fd;
    fd << *pGeom;
    // output is a string. The extention .xml is added automatically
    fd.save(output);
    gsInfo << "Wrote G+Smo file:     " << output << ".xml \n";
    //! [Write geometry]

    //! [Write control points to .csv]
    // writing a .csv file
    if (pGeom->targetDim()==2)
        gsWriteCsv(output+".csv", pGeom->coefs(), {"X","Y"});
    else if (pGeom->targetDim()==3)
        gsWriteCsv(output+".csv", pGeom->coefs(), {"X","Y","Z"});
    gsInfo << "Wrote CSV file:       " << output+".csv \n";
    //! [Write control points to .csv]

    return EXIT_SUCCESS;
}
