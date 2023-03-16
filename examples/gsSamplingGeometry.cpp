/** @file gsSamplingGeometry.cpp

    @brief Take a gsGeometry and make a point-cloud out of it.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <gismo.h>
#include <string>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Options with default values
    bool save     = false;
    std::string in_name("in_geometry");
    index_t numSamples(10000);
    int numPts = 64;
    int numPatch = 0;
    real_t umin = 0.;
    real_t umax = 0.;
    real_t vmin = 0.;
    real_t vmax = 0.;
    bool plot_mesh = true;
    bool plot_net = false;
    std::string fn = "fitting/deepdrawingC.xml";

    // Reading options from the command line
    gsCmdLine cmd("Fit parametrized sample data with a surface patch. Expected input file is an XML "
            "file containing two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
            "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
            "3 x N matrix. Every column represents a point (x,y,z) in space.");
    cmd.addSwitch("save", "Save result in XML format", save);
    cmd.addString("t", "data", "Input sample data", fn);
    cmd.addInt("m", "psize", "number of gridded samples", numPts);
    cmd.addInt("p", "patch", "number of the patch to be gridded", numPatch);
    cmd.addReal("a", "u-min", "Select manually the parametric domain [a,b]x[c,d]: a value.", umin);
    cmd.addReal("b", "u-max", "Select manually the parametric domain [a,b]x[c,d]: b value.", umax);
    cmd.addReal("c", "v-min", "Select manually the parametric domain [a,b]x[c,d]: c value.", vmin);
    cmd.addReal("d", "v-max", "Select manually the parametric domain [a,b]x[c,d]: d value.", vmax);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read data]
    gsFileData<> fd_in(fn);

    std::vector<gsGeometry<>::uPtr> geo = fd_in.getAll< gsGeometry<> >();
    if ( ! geo.empty() ){
        gsInfo<< "Got "<< geo.size() <<" patch"<<(geo.size() == 1 ? "." : "es.") <<"\n";
    gsWriteParaview(memory::get_raw(geo), in_name, numSamples, plot_mesh, plot_net);
    }
    else{
        gsInfo<< "----------------------------------------------" <<"\n";
        gsInfo<< "No geometry has been found in input, quitting." <<"\n";
        gsInfo<< "----------------------------------------------" <<"\n";
        return 0;
    }

    if ( numPatch > geo.size()-1 ){
        gsInfo<< "--------------------------------------------------------------------" <<"\n";
        gsInfo<< "Patch number "<< numPatch <<" has not been found in input, quitting." <<"\n";
        gsInfo<< "--------------------------------------------------------------------" <<"\n";
        return 0;
    }

    gsMatrix<> ab(2,2);
    if (umin + umax + vmin + vmax > 0.)
    {
        ab << umin, umax, vmin, vmax;
        gsInfo << "Selected grid diagonal corners:\n" << ab << "\n";
    }
    else
    {
        ab = geo[numPatch]->support();
        gsInfo << "Grid support diagonal corners:\n" << ab << "\n";
    }

    gsVector<unsigned> numPtsVec(2);
    numPtsVec<<numPts,numPts;
    gsVector<> a = ab.col(0);
    gsVector<> b = ab.col(1);
    gsMatrix<> uv = gsPointGrid(a,b, numPtsVec);
    gsInfo << "Uniform grid:\n" << uv.rows() << " x " << uv.cols() << "\n";
    gsWriteParaviewPoints(uv, "parameters");

    gsMatrix<> xyz;
    geo[numPatch]->eval_into(uv, xyz);
    gsWriteParaviewPoints(xyz, "points");
    gsInfo << "Spatial points:\n" << xyz.rows() << " x " << xyz.cols() << "\n";

    gsFileData<> fd;
    fd << xyz;
    fd << uv;

    fd.dump("dump_out");

    for(index_t p=0; p != geo.size(); p++){
    gsMatrix<> pp;
    geo[p]->eval_into(uv, pp);
    gsWriteParaviewPoints(pp, "datapatch" + internal::to_string(p));
    }

    return 0;
}
