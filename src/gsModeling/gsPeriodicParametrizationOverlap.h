/** @file gsPeriodicParametrizationOverlap.h

    @brief Implementation of periodic Floater parametrization using overlaps.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
*/

#pragma once

#include <gsModeling/gsParametrization.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

/**
 * A class for computing periodic parametrizations of closed
 * (cylinder-like) surface meshes. The result will be periodic in the
 * u-direction and the parameter domain will be [0, 1]^2.
 *
 * Main idea
 * =========
 *
 * The user prescribes the u-values of the points on the lower and
 * upper boundaries (i.e., v = 0 and v = 1, respectively). In
 * addition, they are required to provide an *overlap mesh*: an
 * edge-connected chain of triangles connecting the periodic interface
 * on the bottom and top side of the mesh, which serves as an initial
 * guess for the periodic interface. The triangles corresponding to
 * the overlap are repeated, which allows ``unwrapping" the mesh as.
 *
 * The equations for the inner vertices (depicted with empty circles
 * in the figure below) are assembled by the methods of the parent
 * class. The vertices from the overlap are set to have their
 * u-coordinate smaller or greater by one than that of their twin,
 * depending on whether they are on the left or the right hand side of
 * the overlap (empty and full squares, respectively).
 * 
 * @image html  gsPeriodicParametrizationOverlap-scheme.png
 *
 * @image latex gsPeriodicParametrizationOverlap-scheme.pdf
 *
 * TODO: Mention the paper once it passes the review.
 */
template <class T>
class GISMO_EXPORT gsPeriodicParametrizationOverlap : public gsParametrization<T>
{
    typedef typename gsParametrization<T>::Neighbourhood Neighbourhood;

public:
    /** Constructor
     * @param mesh The surface mesh to be parametrized.
     * @param list List of the method options.
     */
    explicit gsPeriodicParametrizationOverlap(gsMesh<T> &mesh,
					      const gsOptionList &list = gsParametrization<T>::defaultOptions())
	: gsParametrization<T>(mesh, list)
    {
	// Note: m_twins and m_overlapHEM are initiated later on.
    }
    
    /// Analogous main function for periodic parametrizations.
    // TODO: Move the file inputs outside.
    gsPeriodicParametrizationOverlap<T>& compute_periodic_overlap(std::string bottomFile,
								  std::string topFile,
								  std::string overlapFile);

protected:
    /**
     * Analogous to the @a calculate method of the parent class.
     *
     * @param paraMethod which parametrization method to use
     * @param indicesV0 indices of the vertices on the bottom (i.e., v = 0) boundary
     * @param valuesV0 their prescribed parameters
     * @param indicesV1 indices of the vertices on the upper (i.e., v = 1) boundary
     * @param valuesV1 their prescribed parameters
     */
    void calculate_periodic_overlap(const size_t paraMethod,
				    const std::vector<size_t>& indicesV0,
				    const std::vector<T>& valuesV0,
				    const std::vector<size_t>& indicesV1,
				    const std::vector<T>& valuesV1);

    /// Finds the twin of the vertex nr. @a vertexId in the mesh.
    size_t findTwin(size_t vertexId) const
    {
	return this->m_mesh.getVertexIndex(m_overlapHEM.getVertexUnsorted(vertexId));
    }

    /**
     * Helper function to @a constructTwins. Constructs twins on the
     * left or right side of the overlap.
     * @param currentNrAllVertices number of input vertices plus the
     * number of twins constructed so far
     * @param vertexIndices indices of the boundary vertices of the
     * overlap mesh
     * @param from the index (in the overlap mesh) of the starting
     * vertex of the side
     * @param to the index of the end vertex of the side
     * @param rightHandSide flag to decide, whether we are working on
     * the right hand side (true) or left hand side (false) of the
     * overlap.
     */
    void constructTwinsBetween(size_t& currentNrAllVertices,
			       std::list<size_t> vertexIndices,
			       size_t from, size_t to, bool rightHandSide);

    /**
     * Helper function to @a constructTwins, preprocesses the inputs for
     * the overloaded method by finding the indices of the vertices @a from and @a to.
     */
    void constructTwinsBetween(size_t& currentNrAllVertices,
			       std::list<size_t> vertexIndices,
			       typename gsMesh<T>::gsVertexHandle from,
			       typename gsMesh<T>::gsVertexHandle to,
			       bool rightHandSide)
    {
	return constructTwinsBetween(currentNrAllVertices,
				     vertexIndices,
				     m_overlapHEM.findVertex(from),
				     m_overlapHEM.findVertex(to),
				     rightHandSide);
    }

    /**
     * Construct the twins.
     * @param uMinv0 the vertex (of the input mesh) with parameter v = 0 and lowest u value
     * @param uMaxv0 the vertex with parameter v = 0 and highest u value
     * @param uMinv1 the vertex with parameter v = 1 and lowest u value
     * @param uMaxv1 the vertex with parameter v = 1 and highest u value
     */
    void constructTwins(typename gsMesh<T>::gsVertexHandle uMinv0,
			typename gsMesh<T>::gsVertexHandle uMaxv0,
			typename gsMesh<T>::gsVertexHandle uMinv1,
			typename gsMesh<T>::gsVertexHandle uMaxv1);

    /** Analogous to the overloaded function from the parent class but
     * does the periodic thing through the twin trick.*/
    // TODO: Rename (currently, this shadows the parent class implementation).
    void constructAndSolveEquationSystem(const Neighbourhood &neighbourhood,
					 const size_t n,
					 const size_t N);

    /** Helper function to constructAndSolveEquationSystem with twins.
     * Takes the coefficients computed by the standard method and
     * shifts them to the twins.
     * @param lambdas the coefficients of the linear combination for
     * the current vertex
     * @param vertex the id of the current vertex
     */
    void updateLambdasWithTwins(std::vector<T>& lambdas,
				size_t vertexId) const;

    // From here on the visualisation functions.
public:
    /** Creates a flat (i.e., 2D) mesh out of a periodic
     * parametrization created by the overlap method.
     *
     * @param restrict If set to true, the mesh is restricted to [0, 1]^2.
     */
    gsMesh<T> createFlatMesh(bool restrict) const;

    // Cf. https://stackoverflow.com/questions/18100865/compiler-cant-find-base-class-method-when-called-from-derived-and-the-derived
    //using gsParametrization<T>::createFlatMesh;

protected:
    /** Creates a flat (i.e., 2D) mesh with the overlap triangles on
     * both sides of the parametric domain.
     *
     * @param left Indices of the vertices on the left boundary of the
     * overlap mesh (using the same indexing as the twins).
     *
     * @param right Indices of the vertices on the right boundary of
     * the overlap mesh (using the same indexing as the twins).
     */
    gsMesh<T> createExtendedFlatMesh(const std::vector<size_t>& right,
				     const std::vector<size_t>& left) const;


protected: // members

    std::vector<std::pair<size_t, size_t> > m_twins; /**< Every twin pair is a pair of vertices
						      * with (almost) the same coordinates.
						      * The first in the pair has the u-coordinate
						      * smaller by one than the second.
						      */

    gsHalfEdgeMesh<T> m_overlapHEM; ///< The mesh of the overlap, cf. the introduction to the class.
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPeriodicParametrizationOverlap.hpp)
#endif
