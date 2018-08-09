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

static const real_t pi = (3.141592653589793238462643383279502);

using namespace gismo;


int main(int argc, char* argv[])
{

    //! [Parse command line]
    std::string input("motor_conforming.xml");
    std::string output("");

    gsCmdLine cmd("Tutorial Input Output");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addString("o", "output", "Name of the output file", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read geometry]
    if (!gsFileManager::fileExists(input))
    {
        gsWarn << "The file cannot be found!\n";
        return EXIT_FAILURE;
    }

    gsInfo << "Read file \"" << input << "\"\n";

    gsFileData<> fileData(input);

    gsMultiPatch<>::uPtr pGeom;
    if (fileData.has< gsMultiPatch<> >())
    {
        pGeom = fileData.getFirst< gsMultiPatch<> >();
    }
    else
    {
        gsWarn << "Input file doesn't have a geometry.\n";
        return EXIT_FAILURE;
    }
    //! [Read geometry]

    /*
     * create the full cross section of the electric motor
    gsMultiPatch<> mpFull;

    gsMultiPatch<> mpDub(*pGeom);
    mpDub.clearTopology();

    for(unsigned rot = 0; rot < 4; rot++)
    {
        for (unsigned i = 0; i < pGeom->nPatches(); i++)
        {
            gsGeometry <real_t> & geo = (mpDub.patch(i));
            geo.rotate(rot * pi / 2.);
            mpFull.addPatch(geo);
        }

    }

    mpFull.computeTopology();
     */

    gsPiecewiseFunction<real_t> one;
    gsFunctionExpr<real_t> patchOne("1.0", 2);
    for(index_t n = 0; n < pGeom->nPatches(); n++)
        one.addPiece(patchOne);


    //! [Print geometry]
    //gsInfo << "The file contains: \n" << mpFull << "\n";
    //! [Print geometry]

    if ( output.empty() )
    {
        gsInfo << "Call program with option -o <basename> to write data to files\n";
        gsInfo << "<basename>Paraview.vtp, <basename>Paraview.pvd, <basename>.xml\n";
        return EXIT_SUCCESS;
    }

    // writing a paraview file
    const std::string out = output + "Paraview";
    gsField<> field(*pGeom, one);
    gsWriteParaview(field, out, 200);
    gsInfo << "Wrote paraview files: " << out << ".vtp, " << out << ".pvd\n";

    //! [Write geometry]
    // writing a G+Smo .xml file
    //gsFileData<> fd;
    //fd << mpFull;
    // output is a string. The extention .xml is added automatically
    //fd.save(output);
    //gsInfo << "Wrote G+Smo file:     " << output << ".xml \n";
    //! [Write geometry]

    return EXIT_SUCCESS;
}
