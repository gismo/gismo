/** @file geom_integrals.cpp

    @brief Tensor product BSpline surface fitting with HLBFGS optimization.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <iostream>

#include <gismo.h>

using namespace gismo;


int main(int argc, char *argv[])
{


    std::string fn = "../filedata/surfaces/simple.xml"; // f
    std::string output("");

    gsCmdLine cmd("Computing integrals.");
    cmd.addString("f", "filename", "G+Smo input geometry file.", fn);
    cmd.addString("o", "output", "Name of the output file", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsFileData<> fileData(fn);

    gsGeometry<>::uPtr pGeom;
    pGeom = fileData.getFirst< gsGeometry<> >();

    // G+Smo geometries contains basis and coefficients
    gsTensorBSplineBasis<2>& basis = static_cast<gsTensorBSplineBasis<2> &>(pGeom->basis());
    //gsTensorBSplineBasis<2, real_t>& basis = pGeom->basis();

    gsInfo << "\nBasis: \n" << basis << "\n";

    const gsMatrix<>& coefs = pGeom->coefs();
    gsInfo << "\nCoefficients: \n" << coefs << "\n" << "\n";
    //! [printing the geometry]

    gsTensorBSpline<2, real_t> original(basis, coefs);

    gsInfo << "Input geometry:\n " << original << "\n";

    gsInfo << "------------------------------------------------------------\n";



    gsExprAssembler<> A(1,1);
    gsMultiBasis<> dbasis(basis);
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    auto G = A.getMap(original);

    // Set the discretization space
    auto u = A.getSpace(dbasis);
    gsInfo << "basis: " << basis << "\n";
    //u.setup(bc, dirichlet::interpolation, 0);
    u.setup();

    gsInfo << "space setup successful.\n";
    A.initSystem();
    gsInfo<< A.numDofs() << "\n";


    A.assemble(u*meas(G));
    gsVector<> Rhs = A.rhs();
    gsInfo << "Integral of each basis function:\n" << Rhs << "\n";


    A.assemble(u * u.tr() * meas(G));
    gsInfo << "Integral of the square of each basis function:\n" << A.matrix().diagonal() << "\n";



    return EXIT_SUCCESS;

}
