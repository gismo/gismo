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
#include <gsIO/gsIOUtils.h>

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
    template <class E>
    void addField(const expr::_expr<E>& expr, std::string label) {
        GISMO_ENSURE(!m_isSaved,
                     "You cannot add more fields if the gsParaviewDataSet has "
                     "been saved.");
        // evaluates the expression and appends it to the vts files
        // for every patch
        const unsigned nPts = m_options.askInt("numPoints", 1000);
        const unsigned precision = m_options.askInt("precision", 5);
        const bool export_base64 = m_options.askSwitch("base64", false);

        // gsExprEvaluator<real_t> ev;
        // gsMultiBasis<real_t> mb(*m_geometry);
        // ev.setIntegrationElements(mb);
        const std::vector<std::string> tags =
            toVTK(expr, nPts, precision, label, export_base64);
        std::vector<std::string> fnames = filenames();

        for (index_t k = 0; k != m_geometry->nPieces();
             k++)  // For every patch.
        {
            std::ofstream file;
            file.open(fnames[k].c_str(), std::ios_base::app);  // Append to file
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

    /// @brief Evaluates a gsField ( the function part ), and writes that data
    /// to the vtk files.
    /// @tparam T
    /// @param field The gsField to be evaluated
    /// @param label The name that will be displayed in Paraview for this field.
    template <class T>
    void addField(const gsField<T> field, std::string label) {
        GISMO_ENSURE(!m_isSaved,
                     "You cannot add more fields if the gsParaviewDataSet has "
                     "been saved.");
        GISMO_ENSURE(
            (field.parDim() == m_geometry->domainDim() &&
             field.geoDim() == m_geometry->targetDim() &&
             field.nPieces() == m_geometry->nPieces() &&
             field.patches().coefsSize() == m_geometry->coefsSize()),
            "Provided gsField and stored geometry are not compatible!");
        // evaluates the field  and appends it to the vts files
        // for every patch
        const unsigned nPts = m_options.askInt("numPoints", 1000);
        const unsigned precision = m_options.askInt("precision", 5);
        const bool export_base64 = m_options.askSwitch("base64", false);

        const std::vector<std::string> tags =
            toVTK(field, nPts, precision, label, export_base64);
        const std::vector<std::string> fnames = filenames();

        for (index_t k = 0; k != m_geometry->nPieces();
             k++)  // For every patch.
        {
            std::ofstream file;
            file.open(fnames[k].c_str(), std::ios_base::app);  // Append to file
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
        opt.addSwitch("base64", "Export in base64 binary format", false);
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
    /// @param label
    /// @param export_base64 export as base64 encoded string
    /// @return Vector of strings of all <DataArrays>
 template <class T>
 static std::vector<std::string> toVTK(const gsFunctionSet<T>& funSet,
                                       unsigned nPts = 1000,
                                       unsigned precision = 5,
                                       std::string label = "",
                                       const bool& export_base64 = false) {
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

            out.push_back( toDataArray(xyzPoints, label, precision, export_base64)  );
        }
        return out;
 }

 template <class T>
 static std::vector<std::string> toVTK(const gsField<T>& field,
                                       unsigned nPts = 1000,
                                       unsigned precision = 5,
                                       std::string label = "",
                                       const bool& export_base64 = false) {
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

            out.push_back(
                toDataArray(xyzPoints, label, precision, export_base64));
        }
        return out;
 }

 /// @brief  Evaluates one expression over all patches and returns all
 /// <DataArray> xml tags as a vector of strings, as no points are exported, no
 /// need to enforce 1:3D fields
 /// @tparam T
 /// @param expr Expression to be evaluated
 /// @param label The label with which the expression will appear in Paraview
 /// @param export_base64 export as base64 encoded string
 /// @return Vector of strings of all <DataArrays>
 template <class E>
 std::vector<std::string> toVTK(const expr::_expr<E>& expr,
                                unsigned nPts = 1000, unsigned precision = 5,
                                std::string label = "SolutionField",
                                const bool& export_base64 = false) {
     std::vector<std::string> out;
     std::stringstream data_array_stream;
     // Write floating points in fixed value notation (@todo: necessary?),
     // will be ignored in binary export
     data_array_stream.setf(std::ios::fixed);
     data_array_stream.precision(precision);

     // if false, embed topology ?
     const index_t n = m_evaltr->exprData()->multiBasis().nBases();

     gsMatrix<real_t> evaluated_values, bounding_box_dimensions;

     for (index_t i = 0; i != n; ++i) {
         // Get bounding box and sample evaluation points on current patch
         bounding_box_dimensions =
             m_evaltr->exprData()->multiBasis().piece(i).support();
         gsGridIterator<real_t, CUBE> grid_iterator(bounding_box_dimensions,
                                                    nPts);
         m_evaltr->eval(expr, grid_iterator, i);

         // Evaluate Expression on grid_points
         evaluated_values = m_evaltr->allValues(
             m_evaltr->elementwise().size() / grid_iterator.numPoints(),
             grid_iterator.numPoints());

         // Data-Header
         if (export_base64) {
             // Only floating types are supported in this function
             const std::string vtk_typename = []() {
                 if (std::is_same<real_t, float>::value) {
                     return std::string("Float32");
                 } else if (std::is_same<real_t, double>::value) {
                     return std::string("Float64");
                 } else {
                     // Will only be triggered in debug mode, but export will
                     // still work, keyword needs to be added manually
                     GISMO_ASSERT(false,
                                  "Unspported floating point type requested");
                 }
             }();

             // Header
             data_array_stream << "<DataArray type=\"" << vtk_typename
                               << "\" format=\"binary\" ";
             if ("" != label) {
                 data_array_stream << "Name=\"" << label << "\" ";
             };
             // N-dimensional exports are supported (reshape to col might
             // be required)
             data_array_stream << "NumberOfComponents=\""
                               << evaluated_values.rows() << "\">\n";

             // Prepend the number of bytes to be expected (using a
             // single-item array of unsigned 64 integers)
             data_array_stream
                 << Base64::Encode(std::vector<uint64_t>{
                        evaluated_values.cols() * evaluated_values.rows() *
                        sizeof(real_t{})})
                        // Write the actual data
                        // Vector-valued data is stored column-wise
                        + Base64::Encode(evaluated_values, false);
         } else {
             data_array_stream << "<DataArray type=\"Float32\" Name=\"" << label
                               << "\" format=\"ascii\" NumberOfComponents=\""
                               << (evaluated_values.rows() == 1 ? 1 : 3)
                               << "\">\n";
             if (evaluated_values.rows() == 1)
                 for (index_t j = 0; j < evaluated_values.cols(); ++j)
                     data_array_stream << evaluated_values.at(j) << " ";
             else {
                 for (index_t j = 0; j < evaluated_values.cols(); ++j) {
                     for (index_t k = 0; k != evaluated_values.rows(); ++k)
                         data_array_stream << evaluated_values(k, j) << " ";
                     for (index_t k = evaluated_values.rows(); k < 3; ++k)
                         data_array_stream << "0 ";
                 }
             }
         }
         data_array_stream << "\n</DataArray>\n";
         out.push_back(data_array_stream.str());
         data_array_stream.str(
             std::string());  // Clear the data_array_stream stringstream
     }
     return out;
 }

 /// @brief Formats the coordinates of points as a <DataArray> xml tag for
 /// ParaView export.
 /// @tparam T Arithmetic type
 /// @param points A gsMatrix<T> with the coordinates of the points, stored
 /// column-wise. Its size is (numDims, numPoints)
 /// @param label A string with the label of the data
 /// @param precision Number of decimal points in xml output
 /// @param export_base64 (defaults true) export as base63 encoded string -
 /// ignore precision
 /// @return The raw xml string
 template <class T>
 static std::string toDataArray(const gsMatrix<T>& points,
                                const std::string label, unsigned precision,
                                const bool& export_base64) {
     std::stringstream stream;

     if (export_base64) {
         // Only floating types are supported in this function
         const std::string vtk_typename = []() {
             if (std::is_same<T, float>::value) {
                 return std::string("Float32");
             } else if (std::is_same<T, double>::value) {
                 return std::string("Float64");
             }
         }();

         // Write data to vector to ensure it is 3D :/
         std::vector<T> copy_of_matrix;
         copy_of_matrix.reserve(points.cols() * 3);
         // Note : Matrix is transposed
         for (index_t j = 0; j < points.cols(); ++j) {
             for (index_t i = 0; i < points.rows(); ++i) {
                 copy_of_matrix.push_back(points(i, j));
             }
             for (index_t i = points.rows(); i < 3; ++i)
                 copy_of_matrix.push_back(0);
         }

         // Header
         stream << "<DataArray type=\"" << vtk_typename
                << "\" format=\"binary\" ";
         if ("" != label) {
             stream << "Name=\"" << label << "\" ";
         };
         // Only 3D exports are supported
         stream << "NumberOfComponents=\"3\">\n";
         // Prepend the number of bytes to be expected (using a single-item
         // array of unsigned 64 integers)
         stream << Base64::Encode(std::vector<uint64_t>{copy_of_matrix.size() *
                                                        sizeof(T{})})
                       // Write the actual data
                       + Base64::Encode(copy_of_matrix);
     } else {
         stream.setf(std::ios::fixed);  // write floating point values in
                                        // fixed-point notation.
         stream.precision(precision);
         // Format as vtk xml string
         stream << "<DataArray type=\"Float32\" format=\"ascii\" ";
         if ("" != label) stream << "Name=\"" << label << "\" ";
         stream << "NumberOfComponents=\"3\">\n";
         // For every point
         for (index_t j = 0; j < points.cols(); ++j) {
             for (index_t i = 0; i != points.rows(); ++i)
                 stream << points(i, j) << " ";
             for (index_t i = points.rows(); i < 3; ++i) stream << "0 ";
         }
     }
     stream << "\n</DataArray>\n";

     return stream.str();
 }

 void initFilenames();
};
} // End namespace gismo
