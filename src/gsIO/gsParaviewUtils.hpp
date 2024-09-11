/** @file gsParaviewUtils.h

    @brief ParaView output Utilities

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Zwar, C. Karampatzakis
*/
#pragma once

// #include <gsCore/gsForwardDeclarations.h>
#include<gsCore/gsField.h>
#include <gsMSplines/gsMappedBasis.h>   // Only to make linker happy
#include <gsCore/gsDofMapper.h>         // Only to make linker happy
#include <gsCore/gsLinearAlgebra.h>         // Only to make linker happy
#include <gsAssembler/gsExprHelper.h>  
#include <gsAssembler/gsExprEvaluator.h>
#include <gsIO/gsIOUtils.h>


#include <gsLsdyna/gsBextFormat.h>

#include<fstream>
#include<iostream>


namespace gismo
{
    template <class T>
    std::vector<std::string> toVTK(const gsFunctionSet<T>& funSet,
                                   unsigned nPts,
                                   unsigned precision,
                                   std::string label,
                                   const bool& export_base64)
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

            if (xyzPoints.rows() < 3)
                // Pad matrix with zeros
                xyzPoints.conservativeResizeLike(gsEigen::MatrixXd::Zero(3,xyzPoints.cols()));

            if (""!=label)
                out.push_back( toDataArray(xyzPoints, {{"Name",label}}, precision, export_base64)  );
            else
                out.push_back( toDataArray(xyzPoints, {{"",""}}, precision, export_base64)  );
        }
        return out;
    }

    template <class T>
    std::vector<std::string> toVTK(const gsField<T>& field,
                                   unsigned nPts,
                                   unsigned precision,
                                   std::string label,
                                   const bool& export_base64)
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

            if (""!=label)
                out.push_back(
                    toDataArray(xyzPoints, {{"Name", label}}, precision, export_base64));
            else
                out.push_back(
                    toDataArray(xyzPoints, {{"", ""}}, precision, export_base64));

        }
        return out;
    }

    // template <class E>
    // std::vector<std::string> toVTK(const expr::_expr<E>& expr,
    //                                const gsExprEvaluator<> * evaltr,
    //                                unsigned nPts,
    //                                unsigned precision,
    //                                std::string label,
    //                                const bool& export_base64)
    // {
    //     std::vector<std::string> out;
    //     std::stringstream data_array_stream;
    //     // Write floating points in fixed value notation (@todo: necessary?),
    //     // will be ignored in binary export
    //     data_array_stream.setf(std::ios::fixed);
    //     data_array_stream.precision(precision);

    //     // if false, embed topology ?
    //     const index_t n = evaltr->exprData()->multiBasis().nBases();

    //     gsMatrix<real_t> evaluated_values, bounding_box_dimensions;

    //     for (index_t i = 0; i != n; ++i) {
    //         // Get bounding box and sample evaluation points on current patch
    //         bounding_box_dimensions =
    //             evaltr->exprData()->multiBasis().piece(i).support();
    //         gsGridIterator<real_t, CUBE> grid_iterator(bounding_box_dimensions,
    //                                                 nPts);
    //         evaltr->eval(expr, grid_iterator, i);

    //         // Evaluate Expression on grid_points
    //         evaluated_values = evaltr->allValues(
    //             evaltr->elementwise().size() / grid_iterator.numPoints(),
    //             grid_iterator.numPoints());

    //         // Data-Header
    //         if (export_base64) {
    //             // Only floating types are supported in this function
    //             const std::string vtk_typename = []() {
    //                 if (std::is_same<real_t, float>::value) {
    //                     return std::string("Float32");
    //                 } else if (std::is_same<real_t, double>::value) {
    //                     return std::string("Float64");
    //                 } else {
    //                     // Will only be triggered in debug mode, but export will
    //                     // still work, keyword needs to be added manually
    //                     GISMO_ERROR("Unspported floating point type requested");
    //                 }
    //             }();

    //             // Header
    //             data_array_stream << "<DataArray type=\"" << vtk_typename
    //                             << "\" format=\"binary\" ";
    //             if ("" != label) {
    //                 data_array_stream << "Name=\"" << label << "\" ";
    //             };
    //             // N-dimensional exports are supported (reshape to col might
    //             // be required)
    //             data_array_stream << "NumberOfComponents=\""
    //                             << evaluated_values.rows() << "\">\n";

    //             // Prepend the number of bytes to be expected (using a
    //             // single-item array of unsigned 64 integers)
    //             data_array_stream
    //                 << Base64::Encode(std::vector<uint64_t>(1,
    //                 evaluated_values.cols() * evaluated_values.rows() *
    //                                                         sizeof(real_t)))
    //                     // Write the actual data
    //                     // Vector-valued data is stored column-wise
    //                     + Base64::Encode(evaluated_values, false);
    //         } else {
    //             data_array_stream << "<DataArray type=\"Float32\" Name=\"" << label
    //                             << "\" format=\"ascii\" NumberOfComponents=\""
    //                             << (evaluated_values.rows() == 1 ? 1 : 3)
    //                             << "\">\n";
    //             if (evaluated_values.rows() == 1)
    //                 for (index_t j = 0; j < evaluated_values.cols(); ++j)
    //                     data_array_stream << evaluated_values.at(j) << " ";
    //             else {
    //                 for (index_t j = 0; j < evaluated_values.cols(); ++j) {
    //                     for (index_t k = 0; k != evaluated_values.rows(); ++k)
    //                         data_array_stream << evaluated_values(k, j) << " ";
    //                     for (index_t k = evaluated_values.rows(); k < 3; ++k)
    //                         data_array_stream << "0 ";
    //                 }
    //             }
    //         }
    //         data_array_stream << "\n</DataArray>\n";
    //         out.push_back(data_array_stream.str());
    //         data_array_stream.str(
    //             std::string());  // Clear the data_array_stream stringstream
    //     }
    //     return out;
    // }

/*    template <class T>
    std::string toDataArray(const gsMatrix<T>& matrix,
                            std::map<std::string,std::string> attributes,
                            unsigned precision,
                            const bool& export_base64)
    {
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
        if (matrix.rows()>1)
            stream << "NumberOfComponents=\"" << matrix.rows() << "\" ";
        stream << ">\n";


        if (export_base64) {

            // Write data to vector to ensure it is 3D :/
            std::vector<T> copy_of_matrix;
            copy_of_matrix.reserve(matrix.cols() * matrix.rows());
            // Note : Matrix is transposed
            for (index_t j = 0; j < matrix.cols(); ++j) {
                for (index_t i = 0; i < matrix.rows(); ++i) {
                    copy_of_matrix.push_back(matrix(i, j));
                }
            }

            // Prepend the number of bytes to be expected (using a single-item
            // array of unsigned 64 integers)
            stream << Base64::Encode(std::vector<uint64_t>(1,copy_of_matrix.size() *ß
                                                        sizeof(T)))
                        // Write the actual data
                        + Base64::Encode(copy_of_matrix);
        } else {
            stream.setf(std::ios::fixed);  // write floating point values in
                                        // fixed-point notation.
            stream.precision(precision);
            // Format as vtk xml string
            // stream << "<DataArray type=\"Float32\" format=\"ascii\" ";
            // stream << "NumberOfComponents=\"3\">\n";
            // For every point
            for (index_t j = 0; j < matrix.cols(); ++j) {
                for (index_t i = 0; i != matrix.rows(); ++i)
                    stream << matrix(i, j) << " ";
            }
        }
        stream << "\n</DataArray>\n";

        return stream.str();
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


    template<class T>
    // std::string BezierVTK(const gsGeometry<T> & geom)
    std::vector<std::string> BezierVTK(const gsMultiPatch<T> & mPatch)
    {
        std::vector<std::string> out;
        for (index_t p=0;p<mPatch.nPatches();++p)
        {
            std::stringstream stream;
            stream << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n"
            << "<UnstructuredGrid>\n";

            std::map<index_t, ElementBlock> ElBlocks = BezierOperator(mPatch.patch(p).basis() );
            for (auto const& pair : ElBlocks)
            {
                stream << ΕlBlock2vtk(pair.second, mPatch.patch(p).coefs());
            }
            stream << "</UnstructuredGrid>\n";
            stream << "</VTKFile>\n";
            out.push_back(stream.str()); 
        }
        return out;
    }

} // namespace gismo

