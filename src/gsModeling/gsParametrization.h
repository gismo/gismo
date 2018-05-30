/** @file gsParametrization.h

    @brief Class that maintains parametrization

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{
/**
 @brief Class that maintains parametrization
  This class Parametrization stores the mesh information and the two-dimensional parameter points.
  The parameter points are stored in a vector, where the i-th vector element is the parameter point for the vertex with index i.
  This means that the first n elements in the vector are the inner parameter points, the rest of them are the boundary parameter points.

  The parametrization gets constructed with the name of the STL-file containing the mesh, the boundary method and the parametrization method.
  For boundary methods one can choose between
   chords
   corners
   smallest
   opposite
   restrict
   distributed
  and for parametrization method one can choose between
   uniform
   shape
   distance

  There are functions for returning the number of vertices and the number of inner vertices.
  Also every parameter point can be returned.
  The parametrization can be printed by printing all parameter points.

  \ingroup Modeling
 */
template<class T>
class gsParametrization
{
private:
    class MeshInfo
    {
    public:
        MeshInfo() {}
        /**
         * @brief Constructor
         * This constructor uses the STL-file named 'filename' for construction.
         *
         * @param[in] filename const std::string& filename of the STL-file containing the triangle mesh
         **/
        MeshInfo(const std::string &filename) {

        }
    };
public:
    /**
     * @brief Default constructor
     */
    gsParametrization() : m_mesh()
    {}

    /**
     * @brief Constructor
     * This constructor takes as input the filename of the STL-file, the boundary method and the parametrization method.
     * For these properties it calculates the parameter points concerning to Floater's algorithm.
     *
     * @param[in] optionList A gsOptionList with values:
     *
     *  - String filenameIn: name of STL-file
     *  - String filenameOut: output filename (XML)
     *  - String boundaryMethod: name of the boundary calculating method
     *      - chords: choose boundary points distributed on the unit square boundary wrt the chord lengths, no input
     *      needed
     *      - corners: choose 4 boundary corners for the corner points of the unit square and the rest of the boundary
     *      corners is distributed on the four edges of the unit square wrt the chord lengths, input is the corner
     *      numbers, e.g. 1,2,3,4 for the first, second, third and fourth boundary point
     *      - smallest: choose 4 boundary corners as the vertices with the smallest inner angles, no input needed
     *      - restrict: choose first boundary corner as the one with the smallest inner angle, then restrict an
     *      environment of this point (range) and search for the second smallest boundary corner point in all the
     *      others except the restricted ones, etc., input range e.g. r=0.1 for 1/10 of whole boundary length is
     *      restricted around the already chosen corner point
     *      - opposite: choose the boundary corner points such that they are nearly opposite of each other, input range
     *      e.g. r=0.1 for 1/10 of whole boundary length around exact opposite point on boundary is possible
     *      - distributed: choose the smallest inner angle corners (number for how much to choose) and choose four
     *      corners s.t. they are as evenly distributed as possible, input number n=6 for choosing 6 boundary vertices
     *      with smallest inner angles and then find 4 of them s.t. evenly distributed
     *   - String parametrizationMethod:
     *      - shape: best method, shape of the mesh is preserved, smooth surface fitting
     *      - uniform: the lambdas according to floater's algorithm are set to 1/d, where d is the number of neighbours
     *      - distance: the lambdas according to floater's algorithm are set to the relative distances between the point
     *      and its neighbours
     *   - vector<int> corners: vector for corners, call it every time for an entry
     *   - real_t range: in case of restrict or opposite
     *   - index_t number: number of corners, in case of corners
     * @param[in] filename const std::string& - name of STL-file
     * @param[in] boundaryMethod const std::string& - name of the boundary calculating method
     * @param[in] parametrizationMethod const std::string& - name of the parametrization calculating method
     */
    gsParametrization(gsOptionList optionList);

private:
    MeshInfo m_mesh;
    std::vector<gsPoint2D> m_parameterPoints;
}; // class gsParametrization

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsParametrization.hpp)
#endif
