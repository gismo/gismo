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
    std::vector<std::string> toParaview(const gsFunctionSet<T>& funSet,
                                   unsigned nPts = 1000,
                                   unsigned precision = 5,
                                   std::string label = "",
                                   const bool& export_base64 = false);


    template <class T>
    std::vector<std::string> toParaview(const gsField<T>& field,
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
