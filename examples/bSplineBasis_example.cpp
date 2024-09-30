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


    real_t a = 0; // starting knot
    real_t b = 1; // ending knot
    index_t interior = 4; // number of interior knots
    index_t multEnd = 3; // multiplicity at the two end knots
    bool paraview = false;

    gsCmdLine cmd("This is a tutorial on the gsBSplineBasis class.");
    cmd.addReal("","starting","Starting knot",a);
    cmd.addReal("","ending","Ending knot",b);
    cmd.addInt("n","interior","Number of interior knots",interior);
    cmd.addInt("m","mult","Multiplicity at the two end knots",multEnd);
    cmd.addSwitch("plot","Plot with paraview",paraview);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "------------- Constructions -----------------------------\n";

    int degree = multEnd - 1;
    gsKnotVector<> kv(a, b, interior, multEnd);

    gsBSplineBasis<> bsb0(kv);
    print(bsb0, "bsb0");

    gsBSplineBasis<> bsb1(a, b, interior, degree);
    print(bsb1, "bsb1");


    // ======================================================================
    // some properties
    // ======================================================================


    gsInfo << "------------- Some properties    -----------------------\n\n";

    gsInfo << "bsb0.size(): " << bsb0.size() << "\n\n"
              << "bsb0.numElements(): " << bsb0.numElements() << "\n\n"
              << "bsb0.degree(): " << bsb0.degree() << "\n\n";


    // ======================================================================
    // some operations
    // ======================================================================

    gsInfo << "------------- Some operations    -----------------------\n\n";

    const gsKnotVector<>& knots = bsb0.knots();
    gsInfo << "Knots: \n";
    knots.print(gsInfo);
    gsInfo << "\n\n";

    if (paraview)
        printToParaview(bsb0, "basis");

    gsInfo << "bsb0.uniformRefine()\n";
    bsb0.uniformRefine();
    if (paraview)
        printToParaview(bsb0, "basisRefined");

    gsInfo << "bsb0.degreeElevate()\n";
    bsb0.degreeElevate();
    if (paraview)
        printToParaview(bsb0, "basisElevated");
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "files containing the solution.\n";

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


