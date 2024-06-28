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

    std::string input("bspbasis/tpBSpline2_02.xml");
    std::string output("");

    gsCmdLine cmd("Tutorial on gsBasis class.");
    cmd.addPlainString("input", "G+Smo input basis file.", input);
    cmd.addString("o", "output", "Name of the output file.", output);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======================================================================
    // reading the basis
    // ======================================================================

    gsFileData<> fileData(input);

    gsBasis<>::uPtr pBasis;
    if (fileData.has< gsBasis<> >())
    {
        pBasis = fileData.getFirst< gsBasis<> >();
    }
    else
    {
        gsInfo << "Input file doesn't have a basis inside.\n";
        return -1;
    }


    if (!pBasis)
    {
        gsInfo << "Didn't find any basis.\n";
        return -1;
    }


    // ======================================================================
    // printing some information about the basis
    // ======================================================================


    // printing the basis
    gsInfo << "The file contains: \n" << *pBasis << "\n";


    // printing some properties of the basis
    gsInfo << "Dimension of the parameter space: " << pBasis->dim() << "\n"
           << "Number of basis functions: " << pBasis->size() << "\n"
           << "Number of elements: " << pBasis->numElements() << "\n"
           << "Max degree of the basis: " << pBasis->maxDegree() << "\n"
           << "Min degree of the basis: " << pBasis->minDegree() << "\n"
           << "\n";


    // support of the basis
    // (dim x 2 matrix, the parametric domain)
    gsMatrix<> support = pBasis->support();
    gsInfo << "Support: \n"
           << support << "\n\n";



    // ======================================================================
    // evaluation
    // ======================================================================


    // ----------------------------------------------------------------------
    // values
    // ----------------------------------------------------------------------

    gsMatrix<> u = 0.3 * support.col(0) + 0.7 * support.col(1);
    gsInfo << "u " << size(u) << ": \n" << u << "\n\n";

    // indices of active (nonzero) functions at parameter u
    gsMatrix<index_t> active = pBasis->active(u);
    gsInfo << "Active functions at u " << size(active) << ": \n"
           << active << "\n\n";

    gsInfo << "Is number 2 active at the point ? " <<pBasis->isActive(0,u.col(0)) << ": \n"
           << active << "\n\n";


    // values of all active functions at u
    gsMatrix<> values = pBasis->eval(u);
    gsInfo << "Values at u " << size(values) << ": \n"
           << values << "\n\n";

    // values of single basis functions
    for (index_t i = 0; i != active.rows(); i++)
    {
        gsMatrix<> val = pBasis->evalSingle(active(i), u);

        gsInfo << "basis fun. index:  " << active(i)
               << "   value: " << val(0, 0) << "\n";
    }
    gsInfo << "\n";


    // ----------------------------------------------------------------------
    // derivatives
    // ----------------------------------------------------------------------


    gsMatrix<> derivs = pBasis->deriv(u);
    gsInfo << "Derivatives at u " << size(derivs) << ": \n"
           << derivs << "\n\n";


    // derivatives of single basis function
    for (index_t i = 0; i != active.rows(); i++)
    {
        gsMatrix<> der = pBasis->derivSingle(active(i), u);

        gsInfo << "basis fun. index:  " << active(i)
               << "   value: " << std::setw(15) <<  der(0, 0) << "\n";

        for (index_t row = 1; row != der.rows(); row++)
        {
            gsInfo << std::setw(46) << der(row, 0) << "\n";
        }
    }
    gsInfo << "\n";


    // ----------------------------------------------------------------------
    // second derivatives
    // ----------------------------------------------------------------------


    gsMatrix<> derivs2 = pBasis->deriv2(u);
    gsInfo << "Second derivatives at u " << size(derivs2) << ": \n"
           << derivs2 << "\n\n";

    for (index_t i = 0; i != active.rows(); i++)
    {
        gsMatrix<> der2 = pBasis->deriv2Single(active(i), u);

        gsInfo << "basis fun. index:  " << active(i)
               << "   value: " << std::setw(15) << der2(0, 0) << "\n";

        for (index_t row = 1; row != der2.rows(); row++)
        {
            gsInfo << std::setw(46) << der2(row, 0) << "\n";
        }
    }

    gsInfo << "\nFor more information about evaluation "
           << "(and order of derivatives) look at doxygen documentation."
           << "\n\n";



    // ======================================================================
    // writing to a paraview file
    // ======================================================================

    if (output != "")
    {
        gsInfo << "Writing the basis to a paraview file: " << output
               << "\n\n";
        gsWriteParaview(*pBasis, output);
    }
    else
    {
        gsInfo << "Done. No output created, re-run with --output <filename> to get a ParaView "
                  "file containing the solution.\n";
    }

    return 0;
}



