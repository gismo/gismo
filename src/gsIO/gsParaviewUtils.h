/** @file gsParaviewUtils.h

    @brief ParaView output Utilities

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Zwar, C. Karampatzakis
*/
#pragma once

// #include <gismo.h>
// #include <gsCore/gsForwardDeclarations.h>
// #include <gsMSplines/gsMappedBasis.h>   // Only to make linker happy
// #include <gsCore/gsDofMapper.h>         // Only to make linker happy
// #include <gsAssembler/gsExprHelper.h>  
// #include <gsAssembler/gsExprEvaluator.h>
// #include <gsIO/gsIOUtils.h>
#include <gsLsdyna/gsBextFormat.h>
#include <gsIO/gsBase64.h>

#define VTK_BEZIER_QUADRILATERAL 77

namespace gismo
{
    /// @brief  Evaluates gsFunctionSet over all pieces( patches ) and returns all <DataArray> xml tags as a vector of strings
    /// @tparam T 
    /// @param funSet gsFunctionSet to be evaluated
    /// @param nPts   Number of evaluation points, per patch.
    /// @param precision Number of decimal points in xml output
    /// @param label
    /// @param export_base64 export as base64 encoded string
    /// @return Vector of strings of all <DataArrays>
    template <class T>
    std::vector<std::string> toVTK(const gsFunctionSet<T>& funSet,
                                   unsigned nPts = 1000,
                                   unsigned precision = 5,
                                   std::string label = "",
                                   const bool& export_base64 = false);
/*     {
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

            if (""!=label)
                out.push_back( toDataArray(xyzPoints, {{"Name",label}}, precision, export_base64)  );
            else
                out.push_back( toDataArray(xyzPoints, {{"",""}}, precision, export_base64)  );
        }
        return out;
    } */


    template <class T>
    std::vector<std::string> toVTK(const gsField<T>& field,
                                   unsigned nPts = 1000,
                                   unsigned precision = 5,
                                   std::string label = "",
                                   const bool& export_base64 = false);
/*     {
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
    } */

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
                                   const gsExprEvaluator<> & evaltr,
                                   unsigned nPts = 1000, unsigned precision = 5,
                                   std::string label = "SolutionField",
                                   const bool& export_base64 = false);
/*     {
        std::vector<std::string> out;
        std::stringstream data_array_stream;
        // Write floating points in fixed value notation (@todo: necessary?),
        // will be ignored in binary export
        data_array_stream.setf(std::ios::fixed);
        data_array_stream.precision(precision);

        // if false, embed topology ?
        const index_t n = evaltr.exprData()->multiBasis().nBases();

        gsMatrix<real_t> evaluated_values, bounding_box_dimensions;

        for (index_t i = 0; i != n; ++i) {
            // Get bounding box and sample evaluation points on current patch
            bounding_box_dimensions =
                evaltr.exprData()->multiBasis().piece(i).support();
            gsGridIterator<real_t, CUBE> grid_iterator(bounding_box_dimensions,
                                                    nPts);
            evaltr.eval(expr, grid_iterator, i);

            // Evaluate Expression on grid_points
            evaluated_values = evaltr.allValues(
                evaltr.elementwise().size() / grid_iterator.numPoints(),
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
                        GISMO_ERROR("Unspported floating point type requested");
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
                    << Base64::Encode(std::vector<uint64_t>(1,
                    evaluated_values.cols() * evaluated_values.rows() *
                                                            sizeof(real_t)))
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
    } */

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
    std::string toDataArray(const gsMatrix<T>& points,
                            std::map<std::string, std::string> attributes={{"",""}},
                            unsigned precision = 5,
                            const bool& export_base64=false);
/*     {
        std::stringstream stream;

        // Determing 'type' attribute based on input
        const std::string vtk_typename = []() {
            if (std::is_same<T, float>::value) {
                return std::string("Float32");
            } else if (std::is_same<T, double>::value) {
                return std::string("Float64");
            } else if (std::is_same<T, short int>::value) {
                return std::string("Int8");
            } else if (std::is_same<T, unsigned short int>::value) {
                return std::string("UInt8");
            } else if (std::is_same<T, int>::value) {
                return std::string("Int16");
            } else if (std::is_same<T, unsigned int>::value) {
                return std::string("UInt16");
            } else if (std::is_same<T, long int>::value) {
                return std::string("Int32");
            } else if (std::is_same<T, unsigned long int>::value) {
                return std::string("UInt32");
            } else if (std::is_same<T, long long int>::value) {
                return std::string("Int64");
            } else if (std::is_same<T, unsigned long long int>::value) {
                return std::string("UInt64");
            }
        }();


        // Header
        stream << "<DataArray type=\"" << vtk_typename
            << "\" format=\"";
        if (export_base64) stream << "binary\" ";
        else stream << "ascii\" ";

        // Write attributes
        for (auto const& block : attributes)
        {
            if (block.first!="")
            stream << block.first <<"=\""<< block.second <<"\" ";
        }
        stream << "NumberOfComponents=\"" << points.rows() << "\">\n";


        if (export_base64) {

            // Write data to vector to ensure it is 3D :/
            std::vector<T> copy_of_matrix;
            copy_of_matrix.reserve(points.cols() * points.rows());
            // Note : Matrix is transposed
            for (index_t j = 0; j < points.cols(); ++j) {
                for (index_t i = 0; i < points.rows(); ++i) {
                    copy_of_matrix.push_back(points(i, j));
                }
            }

            // Prepend the number of bytes to be expected (using a single-item
            // array of unsigned 64 integers)
            stream << Base64::Encode(std::vector<uint64_t>(1,copy_of_matrix.size() *
                                                        sizeof(T)))
                        // Write the actual data
                        + Base64::Encode(copy_of_matrix);
        } else {
            stream.setf(std::ios::fixed);  // write floating point values in
                                        // fixed-point notation.
            stream.precision(precision);
            // Format as vtk xml string
            stream << "<DataArray type=\"Float32\" format=\"ascii\" ";
            stream << "NumberOfComponents=\"3\">\n";
            // For every point
            for (index_t j = 0; j < points.cols(); ++j) {
                for (index_t i = 0; i != points.rows(); ++i)
                    stream << points(i, j) << " ";
            }
        }
        stream << "\n</DataArray>\n";

        return stream.str();
    } */



    /// @brief ID transformation between G+Smo and vtk  control point IDs 
    /// @param nU Number of control points in u parametric direction
    /// @param nV Number of control points in u parametric direction
    /// @return 
    gsMatrix<index_t> vtkIDTransform(index_t nU, index_t nV);
/*     {
        // T converts coefs from G+Smo's convetnion to ParaView's convention
        gsMatrix<index_t> T(nU*nV, nU*nV);
        T.setZero();

        // T( Paraview , gismo  )
        // Corners ( always 0-3 )
        T(0,0) = 1;
        T(1,nU-1) = 1;
        T(2,nU*nV-1) = 1;
        T(3,nU*nV-nU) = 1;

        // Edges
        for (index_t i=1;i<nU-1;++i) // Parallel to u
        {
            T(3+i,i) = 1;
            T(3+ (nU-2) + (nV-2) + i, nU*(nV-1) + i  ) = 1;


        }
        for (index_t j=1;j<nV-1;++j) // Parallel to v
        {
            T( 3 + (nU-2) + j, (j+1)*nU-1) = 1 ;
            T( 3 + 2*(nU-2) + (nV-2) + j, nU*j) = 1;
        }
        // Internal
        for (index_t i=0;i<nU-2;++i)
        {
            for (index_t j=0;j<nV-2;++j)
            {
                T(2*( nU + nV) - 4 + j*(nU-2) + i , nU*(j+1)+(i+1)) = 1;
            }
        }
        return T;
    } */

    // std::string toDataArray(gsMatrix<index_t> mat, std::map<std::string, std::string> attributes={{"",""}}, index_t ind=0)
    // {
    //     std::stringstream stream;
    //     stream.setf(std::ios::fixed);  // write floating point values in fixed-point notation.
    //     // stream.precision(5);
    //     // stream << std::setfill('0') << std::setw(16);
    //     // Format as vtk xml string
    //     stream << INDENT(ind)  << "<DataArray type=\"Float32\" format=\"ascii\" ";
    //     for (auto const& block : attributes)
    //     {
    //         if (block.first!="")
    //         stream << block.first <<"=\""<< block.second <<"\" ";
    //     }
    //     stream <<">\n";
    //     ++ind;
    //     // For every point
    //     for (index_t i = 0; i < mat.rows(); ++i)
    //     {
    //         stream << INDENT(ind);
    //         for (index_t j = 0; j != mat.cols(); ++j)
    //             stream << mat(i, j) << " ";
    //         for (index_t j = mat.cols(); j < 3; ++j) stream << "0 ";
    //         if (mat.rows()-1==i) break;
    //         stream << "\n";
    //     }
    //     --ind;
    //     stream << "\n"<<INDENT(ind)<<"</DataArray>\n";

    //     return stream.str();
    // }

    // std::string toDataArray(gsMatrix<real_t> mat, std::map<std::string, std::string> attributes={{"",""}}, index_t ind=0)
    // {
    //     std::stringstream stream;

    //     stream.setf(std::ios::fixed);  // write floating point values in fixed-point notation.
    //     // stream.precision(5);
    //     // stream << std::setfill('0') << std::setw(16);
    //     // Format as vtk xml string
    //     stream << INDENT(ind) <<"<DataArray type=\"Float32\" format=\"ascii\" ";
    //     stream << "NumberOfComponents=\""<< /* TODO mat.cols() */ 3 << "\" ";
    //     for (auto const& block : attributes)
    //     {
    //         if (block.first!="")
    //         stream << block.first <<"=\""<< block.second <<"\" ";
    //     }
    //     stream <<">\n";
    //     ++ ind; 
    //     // For every point
    //     for (index_t i = 0; i < mat.rows(); ++i)
    //     {
    //         stream << INDENT(ind);
    //         for (index_t j = 0; j != mat.cols(); ++j)
    //             stream << mat(i, j) << " ";
    //         for (index_t j = mat.cols(); j < 3; ++j) stream << "0 ";
    //         if (mat.rows()-1==i) break;
    //         stream << "\n";
    //     }
    //     --ind;
    //     stream << "\n"<<INDENT(ind)<<"</DataArray>\n";

    //     return stream.str();
    // }


    /// @brief Converts an integer to a 'DataArray' xml tag, which is returned as a string.
    /// @param num The integer to be formatted
    /// @param attributes Optional, map of strings, with attribute name mapping to attribute value.
    /// @param ind Optional, indentation level for the resulting string.
    /// @return 
    std::string toDataArray(index_t num, std::map<std::string, std::string> attributes={{"",""}});
/*     {
        std::stringstream stream;

        stream.setf(std::ios::fixed);  // write floating point values in fixed-point notation.
        // stream.precision(5);
        // stream << std::setfill('0') << std::setw(16);
        // Format as vtk xml string
        stream  <<"<DataArray type=\"Float32\" format=\"ascii\" ";
        for (auto const& block : attributes)
        {
            if (block.first!="")
            stream << block.first <<"=\""<< block.second <<"\" ";
        }
        stream <<">\n" << num << "\n</DataArray>\n";

        return stream.str();
    } */

    std::string ΕlBlock2vtk( ElementBlock ElBlock, const gsMatrix<> geomCoefs);
/*     {
        // Number of all control points of resulting Bezier elements of this block
        index_t totalPoints = ((ElBlock.PR+1) * (ElBlock.PS+1)) * ElBlock.numElements;
        std::stringstream stream;

        // Control point ID transformation matrix
        // from G+Smo notation to Paraview notation
        gsMatrix<> T = vtkIDTransform(ElBlock.PR+1, ElBlock.PS+1);

        // Setup matrices with Cell data
        gsMatrix<index_t> degrees(ElBlock.numElements,3);
        degrees.rowwise() = (gsVector3d<index_t>()<< ElBlock.PR, ElBlock.PS, 0).finished().transpose();

        gsMatrix<index_t> connectivity(1,totalPoints);
        connectivity.asVector().setLinSpaced(0,totalPoints-1);

        gsMatrix<index_t> offsets(1,ElBlock.numElements);
        offsets.asVector().setLinSpaced((ElBlock.PR+1) * (ElBlock.PS+1),totalPoints);

        gsMatrix<index_t> types(1,ElBlock.numElements);
        types.setOnes();
        types*= VTK_BEZIER_QUADRILATERAL;

        // Loop over all elements of the Element Block
        auto Ait = ElBlock.actives.begin();        // Actives Iterator
        auto Cit = ElBlock.coefVectors.begin();    // Coefficients Iteratos

        index_t i = 0;
        gsMatrix<real_t> newCoefs(totalPoints, geomCoefs.cols());
        for(; Ait != ElBlock.actives.end() && Cit != ElBlock.coefVectors.end(); ++Ait, ++Cit)
        {
            newCoefs.middleRows(i*(ElBlock.PR+1) * (ElBlock.PS+1),(ElBlock.PR+1) * (ElBlock.PS+1)) =
                T * *Cit * geomCoefs(Ait->asVector(),gsEigen::all);
            // ID transform * bezier extraction operator * original control points
            ++i;
        } 

        // Format to xml
        stream <<"<Piece NumberOfPoints=\""<< totalPoints<<"\" NumberOfCells=\""<<ElBlock.numElements<<"\">\n";
        // TODO remove
        stream << "<PointData>\n";
        stream << toDataArray( newCoefs, {{"Name","Position"}});
        stream << "</PointData>\n";

        stream << "<CellData HigherOrderDegrees=\"HigherOrderDegrees\">\n";
        stream << "" << toDataArray(degrees,{{"Name","HigherOrderDegrees"},{"NumberOfComponents","3"}});
        stream << "</CellData>\n";

        stream << "<Points>\n";
        stream << ""<< toDataArray( newCoefs, {{"Name","Points"}});
        stream << "</Points>\n";

        stream << "<Cells>\n";
        stream << toDataArray(connectivity, {{"Name","connectivity"}});
        stream << toDataArray(offsets, {{"Name","offsets"}});
        stream << toDataArray(types, {{"Name","types"}});
        stream << "</Cells>\n";
        stream << "</Piece>\n";

        return stream.str(); 
    } */


    template<class T>
    std::string geom2BezVtk(const gsGeometry<T> & geom);
/*     {
        std::map<index_t, ElementBlock> ElBlocks = BezierOperator<2>(geom.basis() );

        std::stringstream stream;
        stream << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n"
        << "\t<UnstructuredGrid>\n";


        for (auto const& pair : ElBlocks)
        {
            stream << ΕlBlock2vtk(pair.second, geom.coefs());
        }
        stream << "\t</UnstructuredGrid>\n";
        stream << "</VTKFile>\n";
        return stream.str(); 
    } */
} // namespace gismo

#undef VTK_BEZIER_QUADRILATERAL

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsParaviewUtils.hpp)
#endif
