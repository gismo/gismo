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
/**
    \brief This class represents a group  of vtk (Paraview) files
    that refer to one multiPatch, for one timestep.

    This class is used by gsParaviewCollection to manage said files, 
    but can be used by the user explicitly as well.

    \ingroup IO
*/
class GISMO_EXPORT gsParaviewDataSet // a collection of .vts files 
{
private:
    std::string m_basename;
    std::vector<std::string> m_filenames;
    const gsMultiPatch<real_t> * m_geometry;
    gsExprEvaluator<real_t> * m_evaltr;
    gsOptionList m_options;
    bool m_isSaved;
    
public:
    /// @brief Basic constructor
    /// @param basename The basename that will be used to create all the individual filenames
    /// @param geometry A gsMultiPatch of the geometry that will be exported and where the fields are defined
    /// @param eval Optional. A gsExprEvaluator, necessary when working with gsExpressions for evaluation purposes
    /// @param options A set of options, if unspecified, defaultOptions() is called.
    gsParaviewDataSet(std::string basename, 
                      gsMultiPatch<real_t> * const geometry,
                      gsExprEvaluator<real_t> * eval=nullptr,
                      gsOptionList options=defaultOptions());

    gsParaviewDataSet():m_basename(""),
                        m_geometry(nullptr),
                        m_evaltr(nullptr),
                        m_options(defaultOptions()),
                        m_isSaved(false)
                        {}
                   
    /// @brief Evaluates an expression, and writes that data to the vtk files.
    /// @tparam E 
    /// @param expr The gsExpression to be evaluated
    /// @param label The name that will be displayed in Paraview for this field.   
    template<class E>
    void addField(const expr::_expr<E> & expr,
                  std::string label)
    {
        GISMO_ENSURE( !m_isSaved, "You cannot add more fields if the gsParaviewDataSet has been saved.");
        // evaluates the expression and appends it to the vts files
        //for every patch
        unsigned nPts = m_options.askInt("numPoints",1000);
        unsigned precision = m_options.askInt("precision",5);

        //gsExprEvaluator<real_t> ev;
        //gsMultiBasis<real_t> mb(*m_geometry);
        //ev.setIntegrationElements(mb);
        std::vector<std::string> tags = toVTK(expr,nPts,precision,label);
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


    /// @brief Recursive form of addField()
    /// @tparam E 
    /// @tparam ...Rest 
    /// @param labels Vector of strings, containing the names of the fields as they will be shown in ParaView.
    /// @param expr The expressions to be evaluated ( arbitrary number of them )
    /// @param ...rest 
    template <class E, typename... Rest>
    void addFields(std::vector<std::string> labels, const expr::_expr<E> & expr, Rest... rest) {
        // keep all but first label 
        GISMO_ENSURE( sizeof...(Rest) == labels.size() - 1, "The length of labels must match the number of expressions provided" );
        std::vector<std::string> newlabels(labels.cbegin()+1, labels.cend());
        

        addField(   expr, labels[0]);       // Add the expression 'expr' with it's corresponding label ( first one )
        addFields(   newlabels, rest...);   // Recursion
    }

    /// @brief Evaluates a gsField ( the function part ), and writes that data to the vtk files.
    /// @tparam T
    /// @param field The gsField to be evaluated
    /// @param label The name that will be displayed in Paraview for this field.   
    template<class T>
    void addField(const gsField<T> field, std::string label)
    {
        GISMO_ENSURE( !m_isSaved, "You cannot add more fields if the gsParaviewDataSet has been saved.");
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

    /// @brief Recursive form of addField()
    /// @tparam T 
    /// @tparam ...Rest 
    /// @param labels Vector of strings, containing the names of the fields as they will be shown in ParaView.
    /// @param field The gsFields to be evaluated ( arbitrary number of them )
    /// @param ...rest 
    template <class T, typename... Rest>
    void addFields(std::vector<std::string> labels, const gsField<T> field, Rest... rest) {
        // keep all but first label 
        std::vector<std::string> newlabels(labels.cbegin()+1, labels.cend());
        
        addField(   field, labels[0]);       // Add the expression 'expr' with it's corresponding label ( first one )
        addFields(   newlabels, rest...);   // Recursion
    }

    /// @brief Returns the names of the files created by this gsParaviewDataSet.
    /// @return A vector of strings
    const std::vector<std::string> filenames();

    void save();

    bool isEmpty();

    bool isSaved();

    /// @brief Accessor to the current options.
    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addInt("numPoints", "Number of points per-patch.", 1000);
        opt.addInt("precision", "Number of decimal digits.", 5);
        opt.addInt("plotElements.resolution", "Drawing resolution for element mesh.", -1);
        opt.addSwitch("makeSubfolder", "Export vtk files to subfolder ( below the .pvd file ).", true);
        opt.addString("subfolder","Name of subfolder where the vtk files will be stored.", "");
        opt.addSwitch("plotElements", "Controls plotting of element mesh.", false);
        opt.addSwitch("plotControlNet", "Controls plotting of control point grid.", false);
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
        gsMatrix<T> evalPoint, xyzPoints;

        // Loop over all patches
        for ( index_t i=0; i != funSet.nPieces(); ++i )
        {
            gsGridIterator<T,CUBE> grid(funSet.piece(i).support(), nPts);

            // Evaluate the MultiPatch for every parametric point of the grid iterator
            xyzPoints.resize( funSet.targetDim(), grid.numPoints() );
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
        gsMatrix<T> evalPoint, xyzPoints;

        // Loop over all patches
        for ( index_t i=0; i != field.nPieces(); ++i )
        {
            gsGridIterator<T,CUBE> grid(field.fields().piece(i).support(), nPts);

            // Evaluate the MultiPatch for every parametric point of the grid iterator
            xyzPoints.resize( field.dim(), grid.numPoints());
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

    /// @brief  Evaluates one expression over all patches and returns all <DataArray> xml tags as a vector of strings
    /// @tparam T 
    /// @param expr Expression to be evaluated
    /// @param label The label with which the expression will appear in Paraview
    /// @return Vector of strings of all <DataArrays>
    template<class E>
    std::vector<std::string> toVTK(const expr::_expr<E> & expr,
                                            unsigned nPts=1000,
                                            unsigned precision=5,
                                            std::string label="SolutionField")
    {   
        std::vector<std::string> out;
        std::stringstream dataArray;
        dataArray.setf( std::ios::fixed ); // write floating point values in fixed-point notation.
        dataArray.precision(precision);
        // m_exprdata->parse(expr);

        //if false, embed topology ?
        const index_t n = m_evaltr->exprData()->multiBasis().nBases();

        gsMatrix<real_t> pts, vals, ab;

        for ( index_t i=0; i != n; ++i )
        {
            ab = m_evaltr->exprData()->multiBasis().piece(i).support();
            gsGridIterator<real_t,CUBE> pt(ab, nPts);
            m_evaltr->eval(expr, pt, i);
            nPts = pt.numPoints();
            
            vals = m_evaltr->allValues(m_evaltr->elementwise().size()/nPts, nPts);

            dataArray <<"<DataArray type=\"Float32\" Name=\""<< label <<"\" format=\"ascii\" NumberOfComponents=\""<< ( vals.rows()==1 ? 1 : 3) <<"\">\n";
            if ( vals.rows()==1 )
                for ( index_t j=0; j<vals.cols(); ++j)
                    dataArray<< vals.at(j) <<" ";
            else
            {
                for ( index_t j=0; j<vals.cols(); ++j)
                {
                    for ( index_t i=0; i!=vals.rows(); ++i)
                        dataArray<< vals(i,j) <<" ";
                    for ( index_t i=vals.rows(); i<3; ++i)
                        dataArray<<"0 ";
                }
            }
            dataArray <<"\n</DataArray>\n";
            out.push_back( dataArray.str() );
            dataArray.str(std::string()); // Clear the dataArray stringstream
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

    void initFilenames();

};
} // End namespace gismo
