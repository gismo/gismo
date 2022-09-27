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
    std::string filename("domain2d/yeti_mp2.xml");
    bool plot = false;
    real_t tol = 1e-5;
    // real_t gtol = 1e-6;
    // bool reparam = false, gaps = true;

    gsCmdLine cmd("Computes patches from structured (tensor-product) data samples by fitting. Give a file path to an XML or 3dm file to refit the patches!");
    cmd.addPlainString("filename", "File containing multipatch input (.xml).", filename);
    cmd.addReal  ("t","tolerance","Tolerance for identifing patch interfaces", tol);
    cmd.addSwitch("plot", "plot results", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<>::uPtr mp = gsReadFile<>(filename);
    if (mp->nInterfaces()==0 && mp->nBoundary()==0)
        mp->computeTopology(tol,true);

    gsInfo <<" Got "<< *mp <<" \n" ;
    // STEP 1: Get curve network with merged linear interfaces
    gsInfo<<"Loading curve network...";
    mp->constructInterfaceRep();
    mp->constructBoundaryRep();
    auto & irep = mp->interfaceRep();
    auto & brep = mp->boundaryRep();

    // outputing...
    gsMultiPatch<> crv_net, iface_net, bnd_net;
    for (auto it = irep.begin(); it!=irep.end(); ++it)
    {
        iface_net.addPatch((*it->second));
        crv_net.addPatch((*it->second));
    }
    for (auto it = brep.begin(); it!=brep.end(); ++it)
    {
        bnd_net.addPatch((*it->second));
        crv_net.addPatch((*it->second));
    }

    if (plot) gsWriteParaview(iface_net,"iface_net",200);
    gsWrite(iface_net,"iface_net");
    if (plot) gsWriteParaview(bnd_net,"bnd_net",200);
    gsWrite(bnd_net,"bnd_net");
    if (plot) gsWriteParaview(crv_net,"crv_net",200);
    gsWrite(crv_net,"crv_net");

    gsInfo<<"Finished\n";
    if (!plot)
        gsInfo<<"Result not plotted but written to iface_net.xml, bnd_net.xml, crv_net.xml. Call with --plot to plot the results directly in ParaView\n";

    return EXIT_SUCCESS;
}
