/** @file gsParaviewDataSet.h

    @brief Provides a helper class to write Paraview (.vts) files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Karampatzakis, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsMSplines/gsMappedBasis.h>   // Only to make linker happy
#include <gsCore/gsDofMapper.h>         // Only to make linker happy
#include <gsAssembler/gsExprHelper.h>  
#include <gsAssembler/gsExprEvaluator.h>

#include<fstream>

namespace gismo 
{

class GISMO_EXPORT gsParaviewDataSet // a collection of .vts files 
{
private:
    std::string m_basename;
    gsMultiPatch<real_t> * m_geometry;
    gsExprEvaluator<real_t> * m_evaltr;
    gsOptionList m_options;
    
public:
    gsParaviewDataSet(std::string basename, 
                      gsMultiPatch<real_t> * geometry,
                      gsExprEvaluator<real_t> * eval, gsOptionList options);

                       


    template<class E>
    void addField(const expr::_expr<E> & expr,
                  std::string label)
    {
        // evaluates the expression and appends it to the vts files
        //for every patch
        unsigned nPts = m_options.askInt("numPoints",1000);
        unsigned precision = m_options.askInt("precision",5);

        std::vector<std::string> tags = m_evaltr->expr2vtk(expr, label,nPts,precision);
        std::vector<std::string> fnames = filenames();

        for ( index_t k=0; k!=m_geometry->nPieces(); k++) // For every patch.
        {
            std::ofstream file;
            file.open( fnames[k].c_str(), std::ios_base::app); // Append to file
            file << tags[k];
            file.close(); 
        }
    }

    // Just here to stop the recursion
    void addFields(std::vector<std::string> labels){} 


    // The recursive case: we take a number, alongside
    // some other numbers, and produce their sum.
    template <class E, typename... Rest>
    void addFields(std::vector<std::string> labels, const expr::_expr<E> & expr, Rest... rest) {
        // keep all but first label 
        std::vector<std::string> newlabels(labels.cbegin()+1, labels.cend());
        

        addField(   expr, labels[0]);       // Add the expression 'expr' with it's corresponding label ( first one )
        addFields(   newlabels, rest...);   // Recursion
    }

    std::vector<std::string> filenames();

    void save();

    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addInt("numPoints", "Number of points per-patch.", 1000);
        opt.addInt("precision", "Number of decimal digits.", 5);
        return opt;
    }

    gsOptionList & options() {return m_options;}


    /// @brief  Evaluates gsFunctionSet over all pieces( patches ) and returns all <DataArray> xml tags as a vector of strings
    /// @tparam T 
    /// @param funSet gsFunctionSet to be evaluated
    /// @param nPts   Number of evaluation points, per patch.
    /// @param precision Number of decimal points in xml output
    /// @return Vector of strings of all <DataArrays>
    template< class T>
    static std::vector<std::string> toVTK(const gsFunctionSet<T> & funSet, unsigned nPts=1000, unsigned precision=5)
    {   
        std::vector<std::string> out;
        std::stringstream dataArray;
        dataArray.setf( std::ios::fixed ); // write floating point values in fixed-point notation.
        dataArray.precision(precision);

        gsMatrix<T> evalPoint;

        // Loop over all patches
        for ( index_t i=0; i != funSet.nPieces(); ++i )
        {
            gsGridIterator<T,CUBE> grid(funSet.basis(i).support(), nPts);

            // Evaluate the MultiPatch for every parametric point of the grid iterator
            gsMatrix<T> xyzPoints( funSet.targetDim(), grid.numPoints());
            index_t col = 0;
            for( grid.reset(); grid; ++grid )
            {
                evalPoint = *grid; // ..
                xyzPoints.col(col) =  funSet.piece(i).eval(evalPoint);
                col++;
            }

            // Format as vtk xml string
            dataArray <<"<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n";
            // For every point
            for ( index_t j=0; j<xyzPoints.cols(); ++j)
            {
                for ( index_t i=0; i!=xyzPoints.rows(); ++i)
                    dataArray<< xyzPoints(i,j) <<" ";
                for ( index_t i=xyzPoints.rows(); i<3; ++i)
                    dataArray<<"0 ";
            }
            dataArray <<"\n</DataArray>\n";

            out.push_back( dataArray.str() );
            dataArray.str(std::string()); // Clear the dataArray stringstream
        }
        return out; 
    }
    

};
} // End namespace gismo
