/** @file geometry_example.cpp

    @brief Tutorial on gsGeometry abstract class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <iostream>
#include <gismo.h>

using namespace gismo;


// Returns the string with the size of a matrix.
template <typename T>
std::string size(const gsMatrix<T>& matrix)
{
    std::string result = "(" + util::to_string(matrix.rows()) + " x " +
        util::to_string(matrix.cols()) + ")";

    return result;
}


int main(int argc, char* argv[])
{

    std::string input("surfaces/simple.xml");
    std::string output("");

    gsCmdLine cmd("Tutorial on gsGeometry class.");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addString("o", "output", "Name of the output file", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======================================================================
    // reading the geometry
    // ======================================================================

    gsFileData<> fileData(input);

    gsGeometry<>::uPtr pGeom;
    if (fileData.has< gsGeometry<> >())
    {
        pGeom = fileData.getFirst< gsGeometry<> >();
    }
    else
    {
        gsInfo << "Input file doesn't have a geometry inside." << "\n";
        return -1;
    }


    if (!pGeom)
    {
        gsInfo << "Didn't find any geometry." << "\n";
        return -1;
    }

    // ======================================================================
    // printing some information about the basis
    // ======================================================================


    // ----------------------------------------------------------------------
    // printing the geometry
    // ----------------------------------------------------------------------

    //! [printing the geometry]
    gsInfo << "The file contains: \n" << *pGeom << "\n";

    // G+Smo geometries contains basis and coefficients
    const gsBasis<>& basis = pGeom->basis();
    gsInfo << "\nBasis: \n" << basis << "\n";

    const gsMatrix<>& coefs = pGeom->coefs();
    gsInfo << "\nCoefficients: \n" << coefs << "\n" << "\n";
    //! [printing the geometry]


    // ----------------------------------------------------------------------
    // printing some properties about the basis
    // ----------------------------------------------------------------------

    //! [printing properties]
    gsInfo << "Dimension of the parameter space: " << pGeom->parDim() << "\n"
              << "Dimension of the geometry: " << pGeom->geoDim() << "\n";

    // support of the geometry, this is the same as gsBasis::support
    // (dim x 2 matrix, the parametric domain)
    gsMatrix<> support = pGeom->support();
    gsInfo << "Support: \n"
              << support << "\n" << "\n";
    //! [printing properties]


    // ======================================================================
    // evaluation
    // ======================================================================


    // ----------------------------------------------------------------------
    // values, 1st derivatives, and 2nd derivatives
    // ----------------------------------------------------------------------

    //! [values and derivatives]
    gsMatrix<> u = 0.3 * support.col(0) + 0.7 * support.col(1);
    gsInfo << "u " << size(u) << ": \n" << u << "\n" << "\n";

    // geoDim x 1 matrix
    gsMatrix<> value = pGeom->eval(u);
    // geoDim x parDim matrix (columns represent gradients)
    gsMatrix<> der1 = pGeom->deriv(u);
    // [geoDim * (parDim + parDim * (parDim - 1) / 2)] x 1 matrix
    gsMatrix<> der2 = pGeom->deriv2(u);

    gsInfo << "Value at u " << size(value) << ": \n"
              << value
              << "\n\nDerivative at u " << size(der1) << ": \n"
              << der1
              << "\n\nSecond derivative at u " << size(der2) << ": \n"
              << der2
              << "\n" << "\n";
    //! [values and derivatives]

    gsInfo << "\nFor more information about evaluation "
              << "(and order of derivatives) look at doxygen documentation."
              << "\n" << "\n";


    // ======================================================================
    // contol net
    // ======================================================================

    //! [control net]
    // mesh holds the control net of a geometry
    // mesh is a set of vertices and lines (connections between vertices)
    gsMesh<> mesh;
    pGeom->controlNet(mesh);
    //! [control net]

    // ======================================================================
    // writing to paraview
    // ======================================================================

    if (!output.empty())
    {
        //! [write to paraview]
        std::string out = output + "Geometry";
        gsInfo << "Writing the geometry to a paraview file: " << out
                  << "\n\n";

        gsWriteParaview(*pGeom, out);

        out = output + "Basis";
        gsInfo << "Writing the basis to a paraview file: " << out
                  << "\n\n";

        gsWriteParaview(basis, out);

        out = output + "ContolNet";
        gsInfo << "Writing the control net to a paraview file: " << out
                  << "\n" << "\n";

        gsWriteParaview(mesh, out);
        //! [write to paraview]
    }
    else
        gsInfo << "Done. No output created, re-run with --output <fn> to get a ParaView "
                  "file containing the solution.\n";

    return 0;
}


