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
    const gsMultiPatch<real_t> * m_geometry;
    gsExprEvaluator<real_t> * m_evaltr;
    gsOptionList m_options;
    
public:
    gsParaviewDataSet(std::string basename, 
                      gsMultiPatch<real_t> * const geometry,
                      gsExprEvaluator<real_t> * eval,
                      gsOptionList options);

    gsParaviewDataSet():m_basename(""),
                        m_geometry(nullptr),
                        m_evaltr(nullptr),
                        m_options(defaultOptions())
                        {}

                       
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

    template<class T>
    void addField(const gsField<T> field, std::string label)
    {
        GISMO_ENSURE((field.parDim()  == m_geometry->domainDim() && 
                      field.geoDim()  == m_geometry->targetDim() &&
                      field.nPieces() == m_geometry->nPieces() && 
                      field.patches().coefsSize() == m_geometry->coefsSize()),
                    "Provided gsField and stored geometry are not compatible!" );
        // evaluates the field  and appends it to the vts files
        //for every patch
        unsigned nPts = m_options.askInt("numPoints",1000);
        unsigned precision = m_options.askInt("precision",5);

        std::vector<std::string> tags = toVTK( field, nPts, precision, label);
        std::vector<std::string> fnames = filenames();

        for ( index_t k=0; k!=m_geometry->nPieces(); k++) // For every patch.
        {
            std::ofstream file;
            file.open( fnames[k].c_str(), std::ios_base::app); // Append to file
            file << tags[k];
            file.close(); 
        }
    }

    template <class T, typename... Rest>
    void addFields(std::vector<std::string> labels, const gsField<T> field, Rest... rest) {
        // keep all but first label 
        std::vector<std::string> newlabels(labels.cbegin()+1, labels.cend());
        
        addField(   field, labels[0]);       // Add the expression 'expr' with it's corresponding label ( first one )
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

private:
    /// @brief  Evaluates gsFunctionSet over all pieces( patches ) and returns all <DataArray> xml tags as a vector of strings
    /// @tparam T 
    /// @param funSet gsFunctionSet to be evaluated
    /// @param nPts   Number of evaluation points, per patch.
    /// @param precision Number of decimal points in xml output
    /// @return Vector of strings of all <DataArrays>
    template< class T>
    static std::vector<std::string> toVTK(const gsFunctionSet<T> & funSet, unsigned nPts=1000, unsigned precision=5, std::string label="")
    {   
        std::vector<std::string> out;
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

            out.push_back( toDataArray(xyzPoints, label, precision)  );
        }
        return out; 
    }

    template< class T>
    static std::vector<std::string> toVTK(const gsField<T> & field, unsigned nPts=1000, unsigned precision=5, std::string label="")
    {   
        std::vector<std::string> out;
        gsMatrix<T> evalPoint;

        // Loop over all patches
        for ( index_t i=0; i != field.nPieces(); ++i )
        {
            gsGridIterator<T,CUBE> grid(field.fields().basis(i).support(), nPts);

            // Evaluate the MultiPatch for every parametric point of the grid iterator
            gsMatrix<T> xyzPoints( field.dim(), grid.numPoints());
            index_t col = 0;
            for( grid.reset(); grid; ++grid )
            {
                evalPoint = *grid; // ..
                xyzPoints.col(col) =  field.value(evalPoint, i);
                col++;
            }

            out.push_back( toDataArray(xyzPoints, label, precision)  );
        }
        return out; 
    }

    /// @brief Formats the coordinates of points as a <DataArray> xml tag for ParaView export.
    /// @tparam T Arithmetic type
    /// @param points A gsMatrix<T> with the coordinates of the points, stored column-wise. Its size is (numDims, numPoints)
    /// @param label A string with the label of the data
    /// @param precision Number of decimal points in xml output
    /// @return The raw xml string 
    template<class T>
    static std::string toDataArray(const gsMatrix<T> & points, const std::string label, unsigned precision)
    {
        std::stringstream stream;
        stream.setf( std::ios::fixed ); // write floating point values in fixed-point notation.
        stream.precision(precision); 
        // Format as vtk xml string
        stream <<"<DataArray type=\"Float32\" format=\"ascii\" ";
        if ( "" != label )
            stream << "Name=\"" << label <<"\" ";
        stream << "NumberOfComponents=\"3\">\n";
        // For every point
        for ( index_t j=0; j<points.cols(); ++j)
        {
            for ( index_t i=0; i!=points.rows(); ++i)
                stream<< points(i,j) <<" ";
            for ( index_t i=points.rows(); i<3; ++i)
                stream<<"0 ";
        }
        stream <<"\n</DataArray>\n";

        return stream.str();
    }

};
} // End namespace gismo
