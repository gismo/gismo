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
#include <gsUtils/gsMesh/gsHalfEdgeMesh.h>

namespace gismo
{
/**
* @brief Class that maintains parametrization
* This class Parametrization stores the mesh information and the two-dimensional parameter points.
* The parameter points are stored in a vector, where the i-th vector element is the parameter point for the vertex with index i.
* This means that the first n elements in the vector are the inner parameter points, the rest of them are the boundary parameter points.
*
* The parametrization gets constructed from a gsHalfEdgeMesh object, the boundary method and the parametrization method.
* For boundary methods one can choose between
*  chords
*  corners
*  smallest
*  opposite
*  restrict
*  distributed
* and for parametrization method one can choose between
*  uniform
*  shape
*  distance
*
* There are functions for returning the number of vertices and the number of inner vertices.
* Also every parameter point can be returned.
* The parametrization can be printed by printing all parameter points.
*/
template<class T>
class gsParametrization
{
public:
    typedef gsPoint<2, T> Point2D;

    // if we use std::vector with static Eigen classes, the second template parameter is needed
	typedef std::vector<Point2D, typename Point2D::aalloc> VectorType;

private:
    gsHalfEdgeMesh<T> m_mesh;     ///< mesh information
	VectorType m_parameterPoints; ///< parameter points
    gsOptionList m_options;

public:

    /// Constructor using the input mesh and (possibly) options
    explicit gsParametrization(gsMesh<T> &mesh, const gsOptionList & list = defaultOptions());

    /// @brief Returns the list of default options for gsParametrization
    static gsOptionList defaultOptions();

    /// Main function which performs the computation
    gsParametrization<T>& compute();

    /**
     * Parametric Coordinates u,v from 0..1
     * @return
     */
    gsMatrix<T> createUVmatrix();

    /**
     * Corresponding mapped values in R3 to parametric coordinates.
     * @return
     */
    gsMatrix<T> createXYZmatrix();

    /**
     * Creates a flat mesh
     * @return
     */
    gsMesh<T> createFlatMesh();

    gsOptionList& options() { return m_options; }

    gsParametrization<T>& setOptions(const gsOptionList& list);

private:

    /**
     * @brief Class that maintains the local neighbourhood properties.
     *
     * As the Floater Algorithm needs some information concerning neighbourhood properties, the LocalNeighbourhood class extracts these information from the triangle mesh.
     * The idea is to have a class object for every inner vertex of the mesh, such that the lambdas appearing in LocalParametrization can be calculated with the information stored in that object.
     * The needed informations are the vertex index, the neighbours given by a chain, the angles between the neighbours and the lengths of the halfedges connecting the vertex index with its neighbours.
     *
     * A local neighbourhood can be constructed be default and by using a MeshInfo object and the vertex index of the local neighbourhood.
     *
     * There are comparison operators == and !=.
     * With the getter functions the vertex index, number of neighbours, vertex indices of neighbours, angles and neighbour distances can be obtained.
     * In case the vertex is a boundary vertex, which is indicated by setting the bool value to FALSE in the constructor, the inner angle can be calculated, too.
     * By inner angle, the sum of angles between neighbours is meant.
     *
     * A local neighbourhood can be printed.
     *
     */
    class LocalNeighbourhood
    {

    public:
        /**
         * @brief Constructor
         *
         * This constructor needs a MeshInfo object and a vertex index as an input.
         * Optional a bool value can be set to 0, which indicates that the vertex is a boundary vertex. With this it should be ensured, that a local neighbourhood is not accidentally constructed for a boundary vertex.
         *
         * It is tested whether vertexIndex > 1 and, vertexIndex < n (inner vertex) or the input innerVertex == 0. Otherwise an error message is printed.
         *
         * For construction of m_neighbours all opposite halfedges are found, chained and angles between origin of the current halfedge, vertex and end of the halfedge are calculated. The neighbour distances are given by the lengths of the halfedges connecting vertex index with its neighbours.
         *
         * A LocalNeighbourhood object can be printed.
         *
         * @param[in] meshInfo const MeshInfo& - mesh information
         * @param[in] vertexIndex const int - vertex index
         * @param[in] innerVertex const bool - optional bool value that indicates if vertex is inner vertex or not
         */
        LocalNeighbourhood(const gsHalfEdgeMesh<T> &meshInfo,
                           const size_t vertexIndex,
                           const bool innerVertex = 1);

        /**
         * @brief Get vertex index
         *
         * @return vertex index
         */
        size_t getVertexIndex() const;

        /**
         * @brief Get number of neighbours
         *
         * This method returns the number of neighbours, which is given by the number of vertices of the chain m_neighbours. For that getNumberOfVertices() of the Chain class is used.
         *
         * @return number of neighbours
         */
        size_t getNumberOfNeighbours() const;

        /**
         * @brief Get vertex indices of neighbours
         *
         * This method returns a list of the indices of the neighbours, which is given by the indices of the chain m_neighbours. For that getVertexIndices() of the Chain class is used.
         *
         * @return list of integer values where each integer is a neighbour's vertex index
         */
        const std::list<size_t> getVertexIndicesOfNeighbours() const;

        /**
         * @brief Get angles
         *
         * This method returns a list of the angles between the vector from the vertex to a neighbour and the vector from the vertex to the next neighbour.
         * The angles are stored in m_angles and were calculated when constructing the object.
         *
         * @return angles stored in a list of T values
         */
        const std::list<T> &getAngles() const;

        /**
         * @brief Get inner angle
         *
         * This method is used for boundary values, meaning the neighbours chain is not closed.
         * It returns the sum of all angles stored in m_angles.
         *
         * @return inner angle, which is the sum of all angles
         */
        T getInnerAngle() const;

        /**
         * @brief Get neighbour distances
         *
         * This method returns a list of T values representing the lengths of the halfedges connecting the vertex with its neighbours stored in m_neighbourDistances.
         *
         * @return list of neighbour distances
         */
        std::list<T> getNeighbourDistances() const;

    private:
        size_t m_vertexIndex; ///< vertex index
        typename gsHalfEdgeMesh<T>::Chain m_neighbours; ///< chain of neighbours
        std::list<T> m_angles; ///< list of angles between neighbours
        std::list<T> m_neighbourDistances; ///< list of distances to neighbours
    };

    /**
     * @brief Class maintains local parametrization
     * This class represents a local parametrization for a point in the triangle mesh, which is identified by the vertex index.
     * The parametrization is given by the weights lambda(i,j) which is the weight of vertex x(i) regarding x(j) according to Floater's algorithm.
     *
     * An object gets constructed by a MeshInfo object, a local neighbourhood and a parametrization method.
     * There is a function for returning the lambdas.
     * A local parametrization is outputted by printing it's positiv lambdas.
     * */
    class LocalParametrization
    {

    public:
        /**
         * @brief Constructor
         * Using this constructor one needs to input mesh information, a local neighbourhood and a parametrization method.
         *
         * @param[in] meshInfo gsHalfEdgeMesh object
         * @param[in] localNeighbourhood local neighbourhood stores the needed information about the neighbours
         * @param[in] parametrizationMethod method used for parametrization, one can choose between
         * 1. shape,
         * 2. uniform,
         * 3. distance
         * */
        LocalParametrization(const gsHalfEdgeMesh<T> &meshInfo,
                             const LocalNeighbourhood &localNeighbourhood,
                             const size_t parametrizationMethod = 2);

        /**
         * @brief Get lambdas
         * The lambdas are returned.
         *
         * @return lambdas
         */
        const std::vector<T> &getLambdas() const;

    private:
        /**
         * @brief Calculate lambdas
         * The lambdas according to Floater's algorithm are calculated.
         *
         * @param[in] N const int - number of vertices of triangle mesh
         * @param[in] points std::vector<Point2D>& - two-dimensional points that have same angles ratio as mesh neighbours
         */
        void calculateLambdas(const size_t N, VectorType& points);

        size_t m_vertexIndex; ///< vertex index
        std::vector<T> m_lambdas; ///< lambdas

    };

    /**
     * @brief Class that maintains neighbourhood information of triangle mesh.
     * Represents the neighbourhood properties of a vertex in the triangle mesh.
     *
     * For the parametrization according to Floater's Algorithm two linear equation systems have to be solved in order to calculate the coordinates of the parameter points.
     * All the information needed to construct these linear equation systems is stored in a Neighbourhood object.
     *
     * Every object of Neighbourhood therefore stores mesh information, the local parametrization for every inner vertex and the local neighbourhood for every boundary vertex.
     * The linear equations are given by
     *  Au=b1 and Av=b2,
     * where A is the matrix with a(i,i) = 1 and a(i,j) = -lambda(i,j) for j != i and
     * the right-hand sides are given by the sum of all lambda(i,j)*u(j) or lambda(i,j)*v(j) for all boundary points, which are calculated beforehand.
     * The boudary points can be calculated in different ways like choosing some particular boundary corner points and distribute the rest of the points evenly or according to their distance.
     *
     * An object is constructed from a gsHalfEdgeMesh object and the desired parametrization method, that can be chosen from 'uniform', 'shape' and 'distance'.
     * Basically the construction is just about constructing the MeshInfo object, the vector of localParametrization objects and the vector of LocalNeighbourhood objects.
     *
     * There are getter functions for all information that is needed to formulate the equation system. E. g. one can get the number of vertices, number of inner vertices, boundary length, number of boundary halfedges, halfedge lengths and lambdas.
     * Furthermore there are functions to find the boundary corner points according to the method.
     *
     * A neighbourhood can be printed.
     * */
    class Neighbourhood
    {

    public:
        /**
         * @brief Default constructor
         */
        //Neighbourhood() { }

        /**
         * @brief Constructor
         *
         * This constructor takes as an input the filename and a parametrization method.
         * The LocalParametrization object then is constructed according to this method.
         *
         * @param[in] meshInfo const gsHalfEdgeMesh<T> object
         * @param[in] parametrizationMethod const size_t - {1:shape, 2:uniform, 3:distance}
         */
        explicit Neighbourhood(const gsHalfEdgeMesh<T> &meshInfo,
                               const size_t parametrizationMethod = 2);

        /**
         * @brief Get vector of lambdas
         *
         * This method returns a vector that stores the lambdas.
         *
         * @return vector of lambdas
         */
        const std::vector<T> &getLambdas(const size_t i) const;

        /**
         * @brief Get boundary corners depending on the method
         *
         * This method returns a vector with the indices of the boundary corners.
         * ...........................................................................................................
         * @return vector of boundary corners
         */
        const std::vector<index_t>
        getBoundaryCorners(const size_t method, const T range = 0.1, const size_t number = 4) const;

        /**
         * @brief
         */
        static const Point2D findPointOnBoundary(T w, size_t index);

    private:
        std::vector<T> midpoints(const size_t numberOfCorners, const T length) const;
        void searchAreas(const T range,
                         std::vector<std::pair<T, size_t> > &sortedAngles,
                         std::vector<index_t> &corners) const;
        void takeCornersWithSmallestAngles(size_t number,
                                           std::vector<std::pair<T, size_t> > &sortedAngles,
                                           std::vector<index_t> &corners) const;

        const gsHalfEdgeMesh<T> & m_basicInfos;
        std::vector<LocalParametrization> m_localParametrizations;
        std::vector<LocalNeighbourhood> m_localBoundaryNeighbourhoods;
    };


private:
    /**
    * @brief Get parameter point
    * Returns the parameter point with given vertex index.
    *
    * @param[in] vertexIndex int - vertex index
    * @return two-dimensional parameter point
    */
    const Point2D &getParameterPoint(size_t vertexIndex) const;

    /**
    * @brief Constructs linear equation system and solves it
    *
    * The last step in Floater's algorithm is solving the linear equation system to obtain the parameter values for the inner vertices of the triangle mesh.
    * This method constructs the above-mentioned system using information from neighbourhood. The matrix is given by
    *  a(i,i) = 1
    *  a(i,j) = -lambda(i,j) for j!=i
    * and the right hand side is calculated using the boundary parameters found beforehand. The parameter values are multiplied with corresponding lambda values and summed up.
    * In the last step the system is solved and the parameter points are stored in m_parameterPoints.
    *
    * @param[in] neighbourhood const Neighbourhood& - neighbourhood information of the mesh
    * @param[in] n const int - number of inner vertices and therefore size of the square matrix
    * @param[in] N const int - number of the vertices and therefore N-n is the size of the right-hand-side vector
    */
    void constructAndSolveEquationSystem(const Neighbourhood &neighbourhood, const size_t n, const size_t N);

    void calculate(const size_t boundaryMethod,
                   const size_t paraMethod,
                   const std::vector<index_t> &cornersInput,
                   const T rangeInput,
                   const size_t numberInput);

    T findLengthOfPositionPart(const size_t position,
                                    const size_t numberOfPositions,
                                    const std::vector<index_t> &bounds,
                                    const std::vector<T> &lengths);

    bool rangeCheck(const std::vector<index_t> &corners, const size_t minimum, const size_t maximum);

}; // class gsParametrization

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsParametrization.hpp)
#endif
