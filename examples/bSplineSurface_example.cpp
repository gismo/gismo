/** @file bSplineSurface_example.cpp

    @brief Tutorial on gsTensorBSpline class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <cmath>
#include <iostream>

#include <gismo.h>


using namespace gismo;

const double PI = 3.14159265;

int main(int argc, char* argv[])
{
    index_t n = 5;
    index_t m = 5;
    index_t degree = 3;
    std::string output("");

    gsCmdLine cmd("Tutorial on gsTensorBSpline class.");
    cmd.addInt   ("n", "dof1", "Number of basis function in one direction"  , n);
    cmd.addInt   ("m", "dof2", "Number of basis function in other direction", m);
    cmd.addInt   ("d", "degree", "Degree of a surface", degree);
    cmd.addString("o", "output", "Name of the output file.", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Adjust values to the minimum required
    degree = math::max( (index_t)(0), degree    );
    n      = math::max(n, degree + 1);
    m      = math::max(m, degree + 1);

    gsInfo << "----------------------\n\n"
              << "n: " << n << "\n\n"
              << "m: " << m << "\n\n"
              << "degree: " << degree << "\n\n"
              << "output: " << output << "\n\n"
              << "----------------------\n\n";

    // 1. construction of a knot vector for each direction
    gsKnotVector<> kv1(0, 1, n - degree - 1, degree + 1);
    gsKnotVector<> kv2(0, 1, m - degree - 1, degree + 1);

    // 2. construction of a basis
    gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);

    // 3. construction of a coefficients
    gsMatrix<> greville = basis.anchors();
    gsMatrix<> coefs (greville.cols(), 3);

    for (index_t col = 0; col != greville.cols(); col++)
    {
        real_t x = greville(0, col);
        real_t y = greville(1, col);

        coefs(col, 0) = x;
        coefs(col, 1) = y;
        coefs(col, 2) = math::sin(x * 2 * PI) * math::sin(y * 2 * PI);
    }

    // 4. putting basis and coefficients toghether
    gsTensorBSpline<2, real_t>  surface(basis, coefs);


    // 5. saving surface, basis and control net to a file
    if (output != "")
    {
        std::string out = output + "Geometry";
        gsInfo << "Writing the surface to a paraview file: " << out
                  << "\n\n";

        gsWriteParaview(surface, out);

        out = output + "Basis";
        gsInfo << "Writing the basis to a paraview file: " << out
                  << "\n\n";

        gsWriteParaview(basis, out);


        out = output + "ContolNet";
        gsInfo << "Writing the control net to a paraview file: " << out
                  << "\n" << "\n";

        gsMesh<> mesh;
        surface.controlNet(mesh);
        gsWriteParaview(mesh, out);
    }
    else
    {
        gsInfo << "Done. No output created, re-run with --output <filename> to get a ParaView "
                  "file containing the solution.\n";
    }

    return 0;
}
