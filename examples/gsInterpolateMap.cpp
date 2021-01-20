/** @file gsInterpolateMap.cpp

    @brief Creates a B-spline patch from an input mathematical function

    Input is coordinate functions X(x,y,z), Y(x,y,z). Z(x,y,z), or
    X(x,y), Y(x,y). Z(x,y), or X(x), Y(x). Z(x). This is controlled by
    the dimension parameter d. The parameters take values in the interval [0,1]

    For instance, we may interpolate:

    - a circle
    ./bin/gsInterpolateMap -d 1 -X "cos(2*pi*x)" -Y "sin(2*pi*x)" -Z "0"

    - a helical curve
    ./bin/gsInterpolateMap -d 1 -X "cos(2*pi*x)" -Y "sin(2*pi*x)" -Z "x"

    - a helical surface
    ./bin/gsInterpolateMap -d 2 -X "x*cos(2*pi*y)" -Y "x*sin(2*pi*y)" -Z "2*pi*y"

    - Hyperboloid
    ./bin/gsInterpolateMap -d 2 -X "sqrt(1+4*(y-1/2)^2)*cos(2*pi*x)" -Y "4*(y-1/2)" -Z "sqrt(1+4*(y-1/2)^2)*sin(2*pi*x)"  -k 8

    - A bottle
    ./bin/gsInterpolateMap -d 2 -X "(2+sin(4*pi*y))*cos(2*pi*x)" -Y "10*y" -Z "(2+sin(4*pi*y))*sin(2*pi*x)"

    - Mobius band
    ./bin/gsInterpolateMap -d 2 -X "(4+2*(y-0.5)*cos(pi*x))*cos(2*pi*x)" -Y "(4+2*(y-0.5)*cos(pi*x))*sin(2*pi*x)" -Z "2*(y-0.5)*sin(pi*x)"

    - An part of an annulus
    ./bin/gsInterpolateMap -d 2 -X "cos(x)*(y+1)" -Y "sin(x)*(y+1)" -Z "0"

    - A quarter of an annulus
    ./bin/gsInterpolateMap -d 2 -X "cos(0.5*pi*x)*(y+1.5)" -Y "sin(0.5*pi*x)*(y+1.5)"

    - a cylinder
    ./bin/gsInterpolateMap -d 2 -X "5*cos(2*pi*y)" -Y "5*sin(2*pi*y)" -Z "15*x"

    - Sphere
    ./bin/gsInterpolateMap -d 2 -X "cos(pi*x)*cos(2*pi*y)" -Y "sin(pi*x)*cos(2*pi*y)" -Z "sin(2*pi*y)" -k 10

    - Torus
    ./bin/gsInterpolateMap -d 2 -X "(3+cos(2*pi*y))*cos(2*pi*x)" -Y "(3+cos(2*pi*y))*sin(2*pi*x)" -Z "sin(2*pi*y)" -k 8

    - Cone
    ./bin/gsInterpolateMap -d 2 -X "cosh(2*pi*x)*cos(2*pi*y)" -Y "cosh(2*pi*x)*sin(2*pi*y)" -Z "sinh(2*pi*x)"

    - 3D pipe
    ./bin/gsInterpolateMap -d 3 -X "cos(2*pi*x)*(y+2)" -Y "sin(2*pi*x)*(y+2)" -Z "4*z"

    - a volume
    ./bin/gsInterpolateMap -d 3 -X "cos(x)*(y+1)" -Y "sin(x)*(y+1)" -Z "z*2"

    - a cylindrical volume
    ./bin/gsInterpolateMap -d 3 -X "x*cos(2*pi*y)" -Y "x*sin(2*pi*y)" -Z "3*z"

    - a degenerate volume
    ./bin/gsInterpolateMap -d 3 -X "5*cos(2*pi*x)" -Y "5*sin(2*pi*y)" -Z "15*z"

    ./bin/gsInterpolateMap -d 3 -X "(z+1)*sqrt(1+4*(y-1/2)^2)*cos(2*pi*x)" -Y "(z+1)*4*(y-1/2)" -Z "(z+1)*sqrt(1+4*(y-1/2)^2)*sin(2*pi*x)"  -k 4

    Try more from http://virtualmathmuseum.org/Surface/gallery_o.html
    (take care of the parameter range ;) )

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
    std::string func_name_x("x");
    std::string func_name_y("y");
    std::string func_name_z("z");
    index_t d = 3; //dimension

    gsCmdLine cmd("Input is coordinate functions X(x,y,z), Y(x,y,z). Z(x,y,z), or"
                  "X(x,y), Y(x,y). Z(x,y), or X(x), Y(x). Z(x). This is controlled by"
                  "the dimension parameter d. The parameters take values in the interval [0,1]");

    cmd.addInt("p", "degree", "this is the degree", p);
    cmd.addInt("k", "knots", "This is the number of interior knots", k);
    cmd.addInt("d", "dim", "this is the parametric dimension", d);
    cmd.addInt("s", "samples", "this is samples for paraview", s);
    cmd.addString("X", "f1", "The X-coordinate of the function", func_name_x);
    cmd.addString("Y", "f2", "The Y-coordinate of the function", func_name_y);
    cmd.addString("Z", "f3", "The Z-coordinate of the function", func_name_z);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Define a function R^d --> R^2
    gsFunctionExpr<> func(func_name_x, func_name_y, func_name_z, d);

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

    gsInfo <<"We are going to interpolate "<< func <<"\n";

    // Support: 2x2 matrix, first column is the lower left
    // point of the parametric domain (eg. 0,0)
    // second column is the upper right corner of the support (eg. 1,1)
    // gsMatrix<> support = tBasis->support();

    // Points to interpolate at (Greville points):
    gsMatrix<> intGrid = tBasis->anchors();
    gsInfo <<"Int. grid dim: "<< intGrid.dim() <<"\n";

    // Evaluate f at the Greville points
    gsMatrix<> fValues = func.eval(intGrid);
    gsInfo <<"Function values dim: "<< fValues.dim() <<"\n";

    // Returns a geometry with basis = tBasis
    // and coefficients being
    // computed as the interpolant of \a funct
    gsGeometry<>::uPtr interpolant = tBasis->interpolateAtAnchors(fValues);

    gsInfo << "Result :"<< *interpolant <<"\n";

    // Save the result as an XML file
    gsFileData<> fd;
    fd << *interpolant;
    std::string remark("Made by interpolation from function: ( ");
    for (index_t i = 0; i<3; ++i)
    {
        remark += func.expression(i);
        remark += ( i==2 ? " )" : " , ");
    }
    fd.addComment( remark );
    fd.save("interpolant_spline");
    gsInfo <<"Result saved as interpolant_spline.xml \n";

    // Produce Paraview file for the interpolant
    //gsWriteParaview( *interpolant, support, "interpolant", s );

    return 0;
}

