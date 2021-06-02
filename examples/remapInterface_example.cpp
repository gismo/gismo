/** @file remapInterface_example.cpp

    @brief Some tests for gsRemapInterface

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <gismo.h>

#include <gsAssembler/gsRemapInterface.h>

using namespace gismo;


void showCorners( const gsGeometry<>& geo )
{
    gsMatrix<> gr = geo.parameterRange();
    const real_t &xmin = gr(0,0), &xmax = gr(0,1);
    const real_t &ymin = gr(1,0), &ymax = gr(1,1);

    gsMatrix<> in(4,2);
    in << xmin,ymin,    xmax,ymin,    xmin,ymax,    xmax,ymax;
    //gsMatrix<> out;
    //geo.eval_into(in.transpose(),out);
    gsInfo << in.transpose() << "\n";

}

int main(int argc, char* argv[])
{

    index_t checkAffine = 3;
    std::string fn("domain2d/yeti_mp2.xml");

    gsCmdLine cmd("remapInterface_example");
    cmd.addInt("c","checkAffine","Number of inner grid points (per direction) to check for affine.\n"
                                 "Iff 0, always assumed affine. Iff -1, never assumed affine. Default: 3",
                                 checkAffine);
    cmd.addString("g","geo","Name of geometry file", fn);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    if ( ! gsFileManager::fileExists(fn) )
    {
        gsInfo << "Geometry file could not be found.\n";
        gsInfo << "I was searching in the current directory and in: " << gsFileManager::getSearchPaths() << "\n";
        return EXIT_FAILURE;
    }

    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(fn);
    gsMultiPatch<>& mp = *mpPtr;

    for (size_t i=0; i!=mp.nPatches(); ++i)
    {
        gsInfo << "Patch " << i << ":\n"; showCorners(mp[i]);
    }

    gsMultiBasis<> mb(mp); // extract basis

    gsInfo << "Number of interfaces is: " << mp.nInterfaces() << "\n";
    index_t jj=0;
    for ( auto mptr = mp.iBegin(); mptr != mp.iEnd(); ++mptr)
    {
        gsInfo << "******************************************\n   Interface number      "  << jj << "\n\n";
        ++jj;
        const boundaryInterface &bi = *mptr;
        gsInfo << "First: " << bi.first() << "\n";
        showCorners(mp[bi.first().patch]);
        gsInfo << "\n\nSecond: " << bi.second() << "\n";
        showCorners(mp[bi.second().patch]);
        gsInfo << "\n\n";
        gsRemapInterface<real_t> ri(mp,mb,bi,checkAffine);
        gsInfo << ri << "\n";

    }
    return 0;
}
