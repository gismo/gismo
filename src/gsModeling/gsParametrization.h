/** @file gsParametrization.h

    @brief Class that maintains parametrization

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl, D. Mokris
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
    explicit gsParametrization(gsMesh<T> &mesh,
			       const gsOptionList & list = defaultOptions(),
			       bool periodic = false);

    /// @brief Returns the list of default options for gsParametrization
    static gsOptionList defaultOptions();

    /// Main function which performs the computation
    gsParametrization<T>& compute();

    /// Analogous main function for periodic parametrizations.
    gsParametrization<T>& compute_periodic_overlap(std::string bottomFile,
						   std::string topFile,
						   std::string overlapFile,
						   std::vector<size_t>& left,
						   std::vector<size_t>& right);

    /// Periodic parametrization using Pierre's trick.
    gsParametrization<T>& compute_periodic_stitch(std::string bottomFile,
						  std::string topFile,
						  std::string stitchFile,
						  std::vector<std::vector<size_t> >& posCorrections);

    gsParametrization<T>& compute_free_boundary();

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

    void restrictMatrices(gsMatrix<T>& uv, gsMatrix<T>& xyz,
			  real_t uMin = 0, real_t uMax = 1, real_t vMin = 0, real_t vMax = 1)
    {
	std::vector<index_t> goodCols;
	for(index_t j=0; j<uv.cols(); j++)
	{
	    if(uv(0, j) >= uMin && uv(0, j) <= uMax && uv(1, j) >= vMin && uv(1, j) <= vMax)
		goodCols.push_back(j);
	}

	// The following is possible, because j <= goodCols[j].
	size_t newSize = goodCols.size();
	for(size_t j=0; j<newSize; j++)
	{
	    uv.col(j)  =  uv.col(goodCols[j]);
	    xyz.col(j) = xyz.col(goodCols[j]);
	}

	gsInfo << newSize << " points remain." << std::endl;

	uv.conservativeResize(2, newSize);
	xyz.conservativeResize(3, newSize);
    }

    /// For Pierre's way.
    void restrictMatrices_2(gsMatrix<T>& uv, gsMatrix<T>& xyz,
			    real_t uMin = 0, real_t uMax = 1)
    {
	real_t uLength = uMax - uMin;
	for(index_t j=0; j<uv.cols(); j++)
	{
	    real_t u = uv(0, j);

	    if(u < uMin)
		uv(0, j) += uLength;
	    else if(u > uMax)
		uv(0 ,j) -= uLength;
	}
    }

    /**
     * Creates a flat mesh
     * @return
     */
    gsMesh<T> createFlatMesh() const;

    /**
     * Creates a flat mesh out of a periodic parametrization created by the overlap method.
     * @param left Indices of the vertices on the left boundary of the parameter domain.
     * @param right Indices of the vertices on the right boundary of the parameter domain.
     */ // TODO: Clarify, which indexing!
    gsMesh<T> createFlatMesh(const std::vector<size_t>& left,
			     const std::vector<size_t>& right,
			     bool restrict = false) const;

    /**
     * Creates a flat mesh out of a periodic parametrization created by a the stitch method.
     * @param posCorrections Positive corrections from the stitch algorithm.
     * @param restrict If set to true, the mesh is restricted to [0, 1]^2.
     */
    gsMesh<T> createFlatMesh(const std::vector<std::vector<size_t> >& posCorrections,
			     bool restrict = false) const;

    /**
     * Creates a copy of the mesh @a original with the u-coordinates shifted by @a uShift.
     * Useful for demonstration purposes and debugging.
     */
    gsMesh<T> createShiftedCopy(const gsMesh<T>& original, real_t uShift) const;

    void writeTexturedMesh(std::string filename) const;

    void writeSTL (const gsMesh<T>& mesh, std::string filename) const;

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

	/// Can be probably integrated into the standard constructor.
        explicit Neighbourhood(const gsHalfEdgeMesh<T> &meshInfo,
			       const std::vector<size_t>& stitchIndices,
			       std::vector<std::vector<size_t> >& posCorrections,
			       std::vector<std::vector<size_t> >& negCorrections,
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

	std::vector<size_t> computeCorrections(const std::vector<size_t>& stitchIndices,
					       const LocalNeighbourhood& localNeighbourhood) const;

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
    void constructAndSolveEquationSystem(const Neighbourhood &neighbourhood,
					 const size_t n,
					 const size_t N);

    /** Similar to @a constructAndSolveEquationSystem but using the NxN system. */
    void constructAndSolveEquationSystem_2(const Neighbourhood &neighbourhood,
					   const size_t n,
					   const size_t N);

    /** Similar to @a _2 but does the periodic thing through the twin trick.*/     // TODO: Remove
    void constructAndSolveEquationSystem(const Neighbourhood &neighbourhood,
					 const size_t n,
					 const size_t N,
					 const std::vector<std::pair<size_t, size_t> >& twins);

    /** Similar to @a constructAndSolveEquationSystem but works for periodic meshes using
     * the corrections.
     */ // TODO: Explain the parameters.
    void constructAndSolveEquationSystem(const Neighbourhood &neighbourhood,
					 const size_t n,
					 const size_t N,
					 const std::vector<std::vector<size_t> >& posCorrections,
					 const std::vector<std::vector<size_t> >& negCorrections);

    std::vector<size_t> readIndices(const std::string& filename) const;
    
    void constructAndSolveEquationSystem(const Neighbourhood &neighbourhood,
					 const size_t n,
					 const size_t N,
					 const std::vector<size_t>& corners,
					 const std::vector<size_t>& botBoundary,
					 const std::vector<size_t>& rgtBoundary,
					 const std::vector<size_t>& topBoundary,
					 const std::vector<size_t>& lftBoundary);

    void constructTwins(std::vector<std::pair<size_t, size_t> >& twins,
			const gsMesh<T>& overlapMesh,
			typename gsMesh<T>::gsVertexHandle u0vMin,
			typename gsMesh<T>::gsVertexHandle u0vMax,
			typename gsMesh<T>::gsVertexHandle u1vMin,
			typename gsMesh<T>::gsVertexHandle u1vMax);

    void calculate(const size_t boundaryMethod,
                   const size_t paraMethod,
                   const std::vector<index_t> &cornersInput,
                   const T rangeInput,
                   const size_t numberInput);

    void calculate_periodic_overlap(const size_t paraMethod,
				    const std::vector<size_t>& indicesU0,
				    const std::vector<T>& valuesU0,
				    const std::vector<size_t>& indicesU1,
				    const std::vector<T>& valuesU1,
				    const gsMesh<T>& overlapMesh,
				    std::vector<size_t>& left,
				    std::vector<size_t>& right);

    void calculate_periodic_stitch(const size_t paraMethod,
				   const std::vector<size_t>& indicesV0,
				   const std::vector<T>& valuesV0,
				   const std::vector<size_t>& indicesV1,
				   const std::vector<T>& valuesV1,
				   const std::vector<size_t>& stitchIndices,
				   std::vector<std::vector<size_t> >& posCorrections);

    std::vector<size_t> getSide(const std::list<size_t>& boundary, size_t beg, size_t end) const;

    void calculate_free_boundary(const size_t paraMethod,
				 const std::string& fileCorners);

    T findLengthOfPositionPart(const size_t position,
                                    const size_t numberOfPositions,
                                    const std::vector<index_t> &bounds,
                                    const std::vector<T> &lengths);

    bool rangeCheck(const std::vector<index_t> &corners, const size_t minimum, const size_t maximum);

    /// Helper function to constructAndSolveEquationSystem with twins.
    void updateLambdasWithTwins(std::vector<T>& lambdas,
				const std::vector<std::pair<size_t, size_t> >& twins,
				size_t vertexId) const;

    real_t correspondingV(const typename gsMesh<T>::VertexHandle& v0,
			  const typename gsMesh<T>::VertexHandle& v1,
			  real_t u) const;

    void addThreeFlatTrianglesOneOut(gsMesh<T>& mesh,
				     const typename gsMesh<T>::VertexHandle& v0,
				     const typename gsMesh<T>::VertexHandle& v1,
				     const typename gsMesh<T>::VertexHandle& v2) const;

    void addThreeFlatTrianglesTwoOut(gsMesh<T>& mesh,
				     const typename gsMesh<T>::VertexHandle& v0,
				     const typename gsMesh<T>::VertexHandle& v1,
				     const typename gsMesh<T>::VertexHandle& v2) const;

    void addOneFlatTriangleNotIntersectingBoundary(gsMesh<T>& mesh,
						   typename gsMesh<T>::VertexHandle& v0,
						   typename gsMesh<T>::VertexHandle& v1,
						   typename gsMesh<T>::VertexHandle& v2) const;

    // TODO: Get rid of this function.
    gsMesh<T> createMidMesh(const std::vector<size_t>& right,
			    const std::vector<size_t>& left) const;

    gsMesh<T> createRestrictedFlatMesh(const gsHalfEdgeMesh<T>& unfolded) const;

}; // class gsParametrization

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsParametrization.hpp)
#endif
