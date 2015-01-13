/** @file tutorialGeometry

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
    std::string result = "(" + internal::toString(matrix.rows()) + " x " + 
        internal::toString(matrix.cols()) + ")";

    return result;
}


int main(int argc, char* argv[])
{

    std::string input("");
    std::string output("");

    try
    {
        // for information about command line arguments 
        // look at tutorialCommandLineArg.cpp
        
        gsCmdLine cmd("Tutorial on gsGeometry class.");
        
        gsArgValPlain<std::string> inArg("input", "G+Smo input geometry file.", 
                                         false, GISMO_DATA_DIR
                                         "surfaces/simple.xml", 
                                         "file", cmd);
        
        gsArgVal<std::string> outArg("o", "output", "Name of the output file.",
                                     false, "", "string", cmd);

        cmd.parse(argc, argv);
        input = inArg.getValue();
        output = outArg.getValue();
    }
    catch (gsArgException& e)
    {
        std::cout << "Error: " << e.error() << " " << e.argId() << std::endl;
        return -1;
    }
    
    // ======================================================================
    // reading the geometry
    // ======================================================================

    gsFileData<> fileData(input);
    
    gsGeometry<>* pGeom = NULL;
    if (fileData.has< gsGeometry<> >())
    {
        pGeom = fileData.getFirst< gsGeometry<> >();
    }
    else
    {
        std::cout << "Input file doesn't have a geometry inside." << std::endl;
        return -1;
    }

        
    if (pGeom == NULL)
    {
        std::cout << "Didn't find any geometry." << std::endl;
        return -1;
    }
    
    // ======================================================================
    // printing some information about the basis
    // ======================================================================

    
    // ----------------------------------------------------------------------
    // printing the geometry
    // ----------------------------------------------------------------------

    std::cout << "The file contains: \n" << *pGeom << std::endl;
    
    // G+Smo geometries contains basis and coefficients
    const gsBasis<>& basis = pGeom->basis();
    std::cout << "\nBasis: \n" << basis << std::endl;
    
    const gsMatrix<>& coefs = pGeom->coefs();
    std::cout << "\nCoefficients: \n" << coefs << "\n" << std::endl;
    
    
    // ----------------------------------------------------------------------
    // printing some properties about the basis
    // ----------------------------------------------------------------------
    
    std::cout << "Dimension of the parameter space: " << pGeom->parDim() << "\n"
              << "Dimension of the geometry: " << pGeom->geoDim() << "\n";

    // support of the geometry, this is the same as gsBasis::support
    // (dim x 2 matrix, the parametric domain)
    gsMatrix<> support = pGeom->support();
    std::cout << "Support: \n"
              << support << "\n" << std::endl;


    
    // ======================================================================
    // evaluation 
    // ======================================================================
    

    // ----------------------------------------------------------------------
    // values, 1st derivatives, and 2nd derivatives
    // ----------------------------------------------------------------------
    
    gsMatrix<> u = 0.3 * support.col(0) + 0.7 * support.col(1);
    std::cout << "u " << size(u) << ": \n" << u << "\n" << std::endl;

    // geoDim x 1 matrix
    gsMatrix<> value = pGeom->eval(u); 
    // geoDim x parDim matrix (columns represent gradients)
    gsMatrix<> der1 = pGeom->deriv(u);
    // [geoDim * (parDim + parDim * (parDim - 1) / 2)] x 1 matrix
    gsMatrix<> der2 = pGeom->deriv2(u);
    
    std::cout << "Value at u " << size(value) << ": \n"
              << value
              << "\n\nDerivative at u " << size(der1) << ": \n"
              << der1 
              << "\n\nSecond derivative at u " << size(der2) << ": \n"
              << der2
              << "\n" << std::endl;
    
    std::cout << "\nFor more information about evaluation "
              << "(and order of derivatives) look at doxygen documentation." 
              << "\n" << std::endl;

    
    // ======================================================================
    // contol net
    // ======================================================================
    
    // mesh holds the control net of a geometry
    // mesh is a set of vertices and lines (connections between vertices)
    gsMesh<> mesh; 
    pGeom->controlNet(mesh);
    
    
    // ======================================================================
    // writing to paraview
    // ======================================================================
    
    if (output != "")
    {
        std::string out = output + "Geometry";
        std::cout << "Writing the geometry to a paraview file: " << out 
                  << "\n\n";
        
        gsWriteParaview(*pGeom, out);
        
        out = output + "Basis";
        std::cout << "Writing the basis to a paraview file: " << out 
                  << "\n\n";
        
        gsWriteParaview(basis, out);
        
        out = output + "ContolNet";
        std::cout << "Writing the control net to a paraview file: " << out 
                  << "\n" << std::endl;
        
        gsWriteParaview(mesh, out);

    }

    delete pGeom;
    return 0;
}
