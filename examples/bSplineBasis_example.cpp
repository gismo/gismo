/** @file bSplineBasis_example.cpp

    @brief Tutorial on the gsBSplineBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

// Look also in tuturialBasis for more functionality of gsBSplineBasis.

#include <iostream>
#include <string>
#include <gismo.h>

using namespace gismo;

// forward declaration of some utility functions
void print(const gsBSplineBasis<>& bsb, const std::string& name);
void printToParaview(const gsBSplineBasis<>& bsb, const std::string& name);

int main(int argc, char* argv[])
{
    // ======================================================================
    // different construction of a knot vector
    // ======================================================================

    /*
    std::vector<real_t> uu {-1, -1,  0,  0, 0,  1,  1, 1,  2,  2};
    gsKnotVector<> kv2(uu,4);
    gsDebugVar(kv2);
    gsBSplineBasis<> b2(kv2);
    gsMatrix<> trans;
    trans.setIdentity(b2.size(), b2.size());
    auto g2 = b2.makeGeometry(give(trans));
    b2.uniformRefine(1,3);
    auto g3 = b2.interpolateAtAnchors( g2->eval(b2.anchors()) );
    gsDebugVar(g3->basis());
    gsDebugVar(g3->coefs());
    //gsDebugVar(g3->coefs().pruned(1e-9));
    */

        gsKnotVector<> kv(0,1,0,1,1,2);
    gsTensorBSplineBasis<2> b(kv,kv);
    gsDebugVar(b);
    gsMatrix<real_t> trans;
    trans.setIdentity(b.size(), b.size());
    auto g = b.makeGeometry(give(trans));
    g->degreeElevate(2);
    gsDebugVar(g->basis());
    gsDebugVar(g->coefs().pruned(1e-9));
    
    return 0;
}

void printToParaview(const gsBSplineBasis<>& bsb,
                     const std::string& name)
{
    gsInfo << "Writing bsb0 to paraview in a file: " << name << "\n\n";
    gsWriteParaview(bsb, name);
}

void print(const gsBSplineBasis<>& bsb,
           const std::string& name)
{
    gsInfo << name << ": \n";
    bsb.print(gsInfo);
    gsInfo << "\n\n";
}


