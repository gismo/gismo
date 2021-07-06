/** @file dg_example.cpp

    @brief Some tests for gsVisitorDg

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <gismo.h>

using namespace gismo;


int main(int argc, char* argv[])
{

    real_t y = 1, s = 1, t = 1;
    real_t alpha = 1, beta  = 1, penalty = -1;
    bool oneSided = false;

    gsCmdLine cmd("dg_example");
    cmd.addReal("y","Domain.Height","",y);
    cmd.addReal("s","Domain.Width1","",s);
    cmd.addReal("t","Domain.Width2","",t);

    cmd.addReal("a","DG.Alpha","",alpha);
    cmd.addReal("b","DG.Beta","",beta);
    cmd.addReal("d","DG.Penalty","",penalty);
    cmd.addSwitch("o","DG.OneSided","",oneSided);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << cmd.getOptionList() << "\n";

    // xlow, ylow, xup, yup, rotate
    std::vector< gsGeometry<>* > pc;
    pc.push_back( gsNurbsCreator<>::BSplineRectangle(0,-s,y,0).release() );
    pc.push_back( gsNurbsCreator<>::BSplineRectangle(0, 0,y,t).release() );
    gsMultiPatch<> mp(pc); // consumes ptrs
    mp.computeTopology();

    gsMultiBasis<> mb(mp); // extract basis

    {
        const boundaryInterface &bi = *(mp.iBegin());

        gsInfo << "First: " << bi.first() << "\n";
        gsInfo << "Second: " << bi.second() << "\n";

        gsOptionList opt = gsGenericAssembler<>::defaultOptions();
        opt.setReal( "DG.Alpha", alpha );
        opt.setReal( "DG.Beta", beta );
        opt.setReal( "DG.Penalty", penalty );
        opt.setSwitch( "DG.OneSided", oneSided );
        opt.setInt( "InterfaceStrategy", iFace::dg );
        gsGenericAssembler<> ass(mp,mb,opt);
        gsSparseMatrix<> dgmat = ass.assembleDG(bi);

        gsInfo << std::fixed << std::setprecision(1) << dgmat.toDense() << "\n\n";
    }

    {
        const boundaryInterface &bi = mp.iBegin()->getInverse();

        gsInfo << "First: " << bi.first() << "\n";
        gsInfo << "Second: " << bi.second() << "\n";

        gsOptionList opt = gsGenericAssembler<>::defaultOptions();
        opt.setReal( "DG.Alpha", alpha );
        opt.setReal( "DG.Beta", beta );
        opt.setReal( "DG.Penalty", penalty );
        opt.setSwitch( "DG.OneSided", oneSided );
        opt.setInt( "InterfaceStrategy", iFace::dg );
        gsGenericAssembler<> ass(mp,mb,opt);
        gsSparseMatrix<> dgmat = ass.assembleDG(bi);

        gsInfo << std::fixed << std::setprecision(1) << dgmat.toDense() << "\n\n";
    }



    return 0;
}
