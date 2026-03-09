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
#include <gsCore/gsForwardDeclarations.h>
#include <gsMSplines/gsMappedBasis.h>   // Only to make linker happy
#include <gsCore/gsDofMapper.h>         // Only to make linker happy
#include <gsExpressions/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>
// #include <gsIO/gsIOUtils.h>



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


    template <class T>
    std::vector<std::string> toVTK(const gsField<T>& field,
                                   unsigned nPts = 1000,
                                   unsigned precision = 5,
                                   std::string label = "",
                                   const bool& export_base64 = false);


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
    template <class MatrixType>
    std::string toDataArray(const MatrixType & matrix,
                            std::map<std::string, std::string> attributes={{"",""}},
                            unsigned precision = 5,
                            const bool& export_base64=false);

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
                                   gsExprEvaluator<> * evaltr,
                                   unsigned nPts = 1000,
                                   unsigned precision = 5,
                                   std::string label = "SolutionField",
                                   const bool& export_base64 = false)
    {
        std::vector<std::string> out;

        // if false, embed topology ?
        const index_t n = evaltr->exprData()->domain().nPieces();

        gsMatrix<real_t> evaluated_values, bounding_box_dimensions;

        for (index_t i = 0; i != n; ++i) {
            // Get bounding box and sample evaluation points on current patch
            bounding_box_dimensions =
                evaltr->exprData()->domain().subdomain(i)->boundingBox();
            gsGridIterator<real_t, CUBE> grid_iterator(bounding_box_dimensions,
                                                    nPts);
            evaltr->eval(expr, grid_iterator, i);

            // Evaluate Expression on grid_points
            evaluated_values = evaltr->allValues(
                evaltr->elementwise().size() / grid_iterator.numPoints(),
                grid_iterator.numPoints());


            GISMO_ASSERT(evaluated_values.rows() <= 3, "The expression can be scalar or have at most 3 components.");
            if (evaluated_values.rows() == 2)
                // Pad matrix with zeros
                evaluated_values.conservativeResizeLike(gsEigen::MatrixXd::Zero(3,evaluated_values.cols()));

            out.push_back(toDataArray(evaluated_values,{{"Name",label}}, precision,export_base64));

        }
        return out;
    }





    /// @brief ID transformation between G+Smo and vtk  control point IDs
    /// @param nU Number of control points in u parametric direction
    /// @param nV Number of control points in u parametric direction
    /// @return
    GISMO_EXPORT gsMatrix<real_t> vtkIDTransform(index_t nU, index_t nV);



    /// @brief Converts an integer to a 'DataArray' xml tag, which is returned as a string.
    /// @param num The integer to be formatted
    /// @param attributes Optional, map of strings, with attribute name mapping to attribute value.
    /// @param ind Optional, indentation level for the resulting string.
    /// @return
    GISMO_EXPORT std::string toDataArray(index_t num, std::map<std::string, std::string> attributes={{"",""}});



    template<class T>
    std::string BezierVTK(const gsMultiPatch<T> & mPatch);
} // namespace gismo

#undef VTK_BEZIER_QUADRILATERAL

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsParaviewUtils.hpp)
#endif
