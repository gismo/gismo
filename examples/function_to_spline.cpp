/** @file gsInterpolateMap.cpp

    @brief Creates a B-spline patch from an input mathematical function



    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Programmed during the workshop https://github.com/gismo/gismo/wiki/G-Smo-developer-days-2015
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char* argv[])
{
    index_t p = 3; // Degree
    index_t k = 3; // Number of interor knots (more knots --> better approximation)
    index_t s = 1000; // samples for plotting
    std::string func_name("sin(pi*x)*cos(pi*y)");
    index_t d = 3; //dimension
    gsCmdLine cmd("Input is coordinate functions X(x,y,z), Y(x,y,z). Z(x,y,z), or"
                  "X(x,y), Y(x,y). Z(x,y), or X(x), Y(x). Z(x). This is controlled by"
                  "the dimension parameter d. The parameters take values in the interval [0,1]");

    cmd.addInt("p", "degree", "this is the degree", p);
    cmd.addInt("k", "knots", "This is the number of interior knots", k);
    cmd.addInt("d", "dim", "this is the parametric dimension", d);
    cmd.addInt("s", "samples", "this is samples for paraview", s);
    cmd.addString("F", "func", "The function to interpolate", func_name);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Define a function R^d --> R^2
    gsFunctionExpr<> func(func_name, d);

    // Define a B-spline space of dimension 2
    // start, end, numInt, multiplicities at end-points
    gsKnotVector<> KV(0.0, 1.0, k, p+1);
    gsBasis<>::uPtr tBasis;

    switch (d) // static dispatch..
    {
    case 1:
        tBasis = gsBSplineBasis<>::make(KV);
        break;
    case 2:
        tBasis = memory::make_unique(new gsTensorBSplineBasis<2>(KV,KV));
        break;
    case 3:
        tBasis = memory::make_unique(new gsTensorBSplineBasis<3>(KV,KV,KV));
        break;
    default:
    {
        gsWarn<<"Dimension must be 1, 2 or 3.";
        return 0;
    }
    };

    // gsFileData<> fd_in(basis_name);
    // gsBasis<>::uPtr tBasis = fd_in.getFirst< gsBasis<> >();

    gsInfo <<"We are going to project "<< func <<"\n";

    writeSingleGeometry(func,tBasis->support(),"Func",100000);

    gsGeometry<>::uPtr geom;
    switch (d) // static dispatch..
    {
    case 1:
        geom = gsNurbsCreator<>::BSplineUnitInterval(1);
        break;
    case 2:
        geom = gsNurbsCreator<>::BSplineSquare();
        break;
    case 3:
        geom = gsNurbsCreator<>::BSplineCube();
        break;
    default:
    {
        gsWarn<<"Dimension must be 1, 2 or 3.";
        return 0;
    }
    };

    gsMatrix<> coefs;
    gsL2Projection<real_t>::project(*tBasis, *geom, func, coefs);
    gsGeometry<>::uPtr projected = tBasis->makeGeometry(coefs);

    gsInfo << "Result :"<< *projected <<"\n";

    // Save the result as an XML file
    gsFileData<> fd;
    gsMultiPatch<> mp;
    mp.addPatch(give(projected));
    fd << mp;
    fd.save("projected_spline");
    gsInfo <<"Result saved as projected_spline.xml \n";

    // Produce Paraview file for the projected
    //gsWriteParaview( *projected, support, "projected", s );

    return 0;
}

