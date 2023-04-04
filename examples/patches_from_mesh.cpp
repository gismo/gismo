/** @file refit_patches.cpp

    @brief Computes patches from structured (tensor-product) data samples by fitting.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    std::string filename("off/neon_side.off");
    std::string out;
    bool plot = false;

    gsCmdLine cmd("Computes patches from structured (tensor-product) data samples by fitting. Give a file path to an XML or 3dm file to refit the patches!");
    cmd.addPlainString("filename", "File containing multipatch input (.xml).", filename);
    cmd.addSwitch("plot", "plot results", plot);
    cmd.addString("w","write", "write results", out);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(!filename.empty(),"Filename is empty");
    gsStopwatch time;
    gsFileData<> fd(filename);

    gsSurfMesh HEmesh;
    fd.getFirst<gsSurfMesh>(HEmesh);

    gsMultiPatch<> mp = HEmesh.linear_patches();

    gsInfo<<"Reading mesh:\t"<<time.stop()<<" seconds\n";
    if (plot) gsWriteParaview(HEmesh,"HEmesh");

    time.restart();


    gsInfo<<"Making multipatch:\t"<<time.stop()<<" seconds\n";
    time.restart();
    if (plot)
    {
        gsWriteParaview(mp,"mp",1000,true);
        gsInfo<<"Plotting multipatch:\t"<<time.stop()<<" seconds\n";
    }
    time.restart();
    if (!out.empty())
    {
        gsWrite(mp,out);
        gsInfo<<"Writing multipatch:\t"<<time.stop()<<" seconds\n";
    }

    return EXIT_SUCCESS;
}
