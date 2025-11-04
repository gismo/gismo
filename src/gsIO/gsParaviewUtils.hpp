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
#include <gsExpressions/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsIO/gsBase64.h>


#include<fstream>
#include<iostream>

#define VTK_BEZIER_QUADRILATERAL 77

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

    template <class MatrixType>
    std::string toDataArray(const MatrixType & matrix,
                            std::map<std::string, std::string> attributes,
                            unsigned precision,
                            const bool& export_base64)
    {
        std::stringstream stream;

        typedef typename MatrixType::Scalar T;

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
            stream << Base64::Encode(std::vector<uint64_t>(1,copy_of_matrix.size() *
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
    }


    template<class T>
    // std::string BezierVTK(const gsGeometry<T> & geom)
    std::string BezierVTK(const gsMultiPatch<T> & mPatch)
    {
        std::stringstream stream;

        stream << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n"
        << "<UnstructuredGrid>\n";


        const gsMultiPatch<T> bezierExt = mPatch.extractBezier();
        index_t totalPoints = bezierExt.coefsSize();

        // Set up matrices with cell data
        gsMatrix<index_t> degrees(bezierExt.nPatches(),3);

        gsMatrix<index_t> connectivity(1,totalPoints);
        connectivity.asVector().setLinSpaced(0,totalPoints-1);

        gsMatrix<index_t> offsets(1,bezierExt.nPatches());
        index_t offset=0;

        gsMatrix<index_t> types(1, bezierExt.nPatches());
        types.setOnes();
        types*= VTK_BEZIER_QUADRILATERAL;

        gsMatrix<T> coefs(bezierExt.coefsSize(), bezierExt.targetDim());
        gsMatrix<T> weights(bezierExt.coefsSize(),1);
        weights.setOnes();

        // For every patch / bezier element
        for (size_t p=0;p<bezierExt.nPatches();++p)
        {
            // Control point ID transformation matrix
            // from G+Smo notation to Paraview notation
            gsMatrix<> IdTransform = vtkIDTransform(bezierExt.patch(p).degree(0)+1, bezierExt.patch(p).degree(1)+1);

            // Fill matrices with Cell data
            degrees.row(p) << bezierExt.patch(p).degree(0), bezierExt.patch(p).degree(1), 0;
            offsets(p) = offset + bezierExt.patch(p).coefsSize();

            // Tranform control points to vtk ordering and concatenate in coefs
            coefs.middleRows(offset,bezierExt.patch(p).coefsSize()) = IdTransform * bezierExt.patch(p).coefs();
            if (bezierExt.basis(p).isRational())
                weights.middleRows(offset,bezierExt.patch(p).coefsSize()) = IdTransform * bezierExt.basis(p).weights();
            offset += bezierExt.patch(p).coefsSize();
        }

        if (coefs.cols() == 2)
            // Pad matrix with zeros
            coefs.conservativeResizeLike(gsEigen::MatrixXd::Zero(coefs.rows(),3));

        // Format to xml
        stream <<"<Piece NumberOfPoints=\""<< totalPoints<<"\" NumberOfCells=\""<< bezierExt.nPatches()<<"\">\n";

        stream << "<PointData RationalWeights=\"RationalWeights\">\n";
        stream << "" << toDataArray( weights.transpose() ,{{"Name","RationalWeights"}});
        stream << "</PointData>\n";


        stream << "<CellData HigherOrderDegrees=\"HigherOrderDegrees\">\n";
        stream << "" << toDataArray(degrees.transpose(),{{"Name","HigherOrderDegrees"}});
        stream << "</CellData>\n";

        stream << "<Points>\n";
        stream << ""<< toDataArray( coefs.transpose(), {{"Name","Points"}});
        stream << "</Points>\n";

        stream << "<Cells>\n";
        stream << toDataArray(connectivity, {{"Name","connectivity"}});
        stream << toDataArray(offsets, {{"Name","offsets"}});
        stream << toDataArray(types, {{"Name","types"}});
        stream << "</Cells>\n";
        stream << "</Piece>\n";
        stream << "</UnstructuredGrid>\n";
        stream << "</VTKFile>\n";

        return stream.str();
    }

} // namespace gismo
#undef VTK_BEZIER_QUADRILATERAL

