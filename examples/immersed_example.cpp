/** @file basis_example.cpp

    @brief Tutorial on gsBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <iostream>
#include <gismo.h>

using namespace gismo;

int main(int argc, char* argv[])
{

    gsCmdLine cmd("Tutorial on gsBasis class.");
    // cmd.addPlainString("input", "G+Smo input basis file.", input);
    // cmd.addString("o", "output", "Name of the output file.", output);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::NurbsCircle());

    // defines the inside or the outside of the domain
    gsFunctionExpr<> inOut("if (x^2+y^2-1 < 0, 1, 0)",2);
    // gsWriteParaview(inOut,mp,"inOut");

    gsKnotVector<> kv(0,1,0,2);
    gsTensorBSplineBasis<2,real_t> tbs(kv,kv);



    gsVector<unsigned,2> upp;
    upp.setConstant(kv.numElements());

    // make the tree
    gsHDomain<2,unsigned> tree;
    tree.init(upp,10);

    gsDebugVar(upp);


    return 0;
}



