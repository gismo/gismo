/** @file tutorialBasis

    @brief Tutorial on gsBasis abstract class.

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
    
    // for information about command line arguments 
    // look at tutorialCommandLineArg.cpp
        
    gsCmdLine cmd("Tutorial about gsBasis class.");
        
    gsArgValPlain<std::string> inArg("input", "G+Smo input basis file.", 
                                         false, GISMO_DATA_DIR
                                         "bspbasis/tpBSpline2_02.xml", 
                                         "file", cmd);

    cmd.addString("o", "output", "Name of the output file.", output);
    
    bool ok = cmd.getValues(argc,argv);
    input = inArg.getValue();
    if (!ok) {
        std::cout << "Error during parsing command line!";
        return 1;
    }
    
    // ======================================================================
    // reading the basis
    // ======================================================================
    
    gsFileData<> fileData(input);
    
    gsBasis<>* pBasis = NULL;
    if (fileData.has< gsBasis<> >())
    {
        pBasis = fileData.getFirst< gsBasis<> >();
    }
    else
    {
        std::cout << "Input file doesn't have a basis inside." << std::endl;
        return -1;
    }

        
    if (pBasis == NULL)
    {
        std::cout << "Didn't find any basis." << std::endl;
        return -1;
    }


    // ======================================================================
    // printing some information about the basis
    // ======================================================================
    
    
    // printing the basis
    std::cout << "The file contains: \n" << *pBasis << std::endl;
            
    
    // printing some properties of the basis
    std::cout << "Dimension of the parameter space: " << pBasis->dim() << "\n"
              << "Number of basis functions: " << pBasis->size() << "\n"
              << "Number of elements: " << pBasis->numElements() << "\n"
              << "Max degree of the basis: " << pBasis->maxDegree() << "\n"
              << "Min degree of the basis: " << pBasis->minDegree() << "\n"
              << std::endl;


    // support of the basis 
    // (dim x 2 matrix, the parametric domain)
    gsMatrix<> support = pBasis->support();
    std::cout << "Support: \n"
              << support << "\n" << std::endl;
    

    
    // ======================================================================
    // evaluation 
    // ======================================================================
    

    // ----------------------------------------------------------------------
    // values
    // ----------------------------------------------------------------------

    gsMatrix<> u = 0.3 * support.col(0) + 0.7 * support.col(1);
    std::cout << "u " << size(u) << ": \n" << u << "\n" << std::endl;
    
    // indices of active (nonzero) functions at parameter u
    gsMatrix<unsigned> active = pBasis->active(u);
    std::cout << "Active functions at u " << size(active) << ": \n" 
              << active << "\n" << std::endl;
    

    // values of all active functions at u
    gsMatrix<> values = pBasis->eval(u);
    std::cout << "Values at u " << size(values) << ": \n"
              << values << "\n" << std::endl;
    
    // values of single basis functions
    for (index_t i = 0; i != active.rows(); i++)
    {
        gsMatrix<> val = pBasis->evalSingle(active(i), u);
        
        std::cout << "basis fun. index:  " << active(i) 
                  << "   value: " << val(0, 0) << "\n";
    }
    std::cout << std::endl;


    // ----------------------------------------------------------------------
    // derivatives
    // ----------------------------------------------------------------------
    

    gsMatrix<> derivs = pBasis->deriv(u);
    std::cout << "Derivatives at u " << size(derivs) << ": \n"
              << derivs << "\n" << std::endl;
    

    // derivatives of single basis function
    for (index_t i = 0; i != active.rows(); i++)
    {
        gsMatrix<> der = pBasis->derivSingle(active(i), u);
        
        std::cout << "basis fun. index:  " << active(i)
                  << "   value: " << std::setw(15) <<  der(0, 0) << "\n";
        
        for (index_t row = 1; row != der.rows(); row++)
        {
            std::cout << std::setw(46) << der(row, 0) << "\n";
        }
    }
    std::cout << std::endl;


    // ----------------------------------------------------------------------
    // second derivatives
    // ----------------------------------------------------------------------
    

    gsMatrix<> derivs2 = pBasis->deriv2(u);
    std::cout << "Second derivatives at u " << size(derivs2) << ": \n"
              << derivs2 << "\n" << std::endl;
    
    for (index_t i = 0; i != active.rows(); i++)
    {
        gsMatrix<> der2 = pBasis->deriv2Single(active(i), u);

        std::cout << "basis fun. index:  " << active(i)
                  << "   value: " << std::setw(15) << der2(0, 0) << "\n";
        
        for (index_t row = 1; row != der2.rows(); row++)
        {
            std::cout << std::setw(46) << der2(row, 0) << "\n";
        }
    }

    std::cout << "\nFor more information about evaluation "
              << "(and order of derivatives) look at doxygen documentation." 
              << "\n" << std::endl;
    


    // ======================================================================
    // writing to a paraview file
    // ======================================================================
    
    if (output != "")
    {
        std::cout << "Writing the basis to a paraview file: " << output 
                  << "\n" << std::endl;
        gsWriteParaview(*pBasis, output);
    }
    
    delete pBasis;
    return 0;
}



