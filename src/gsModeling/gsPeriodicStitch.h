/** @file gsPeriodicStitch.h

    @brief Implementation of periodic Floater parametrization using a
    stitch. This class is an alternative to gsPeriodicOverlap.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
*/

#pragma once

#include <gsIO/gsOptionList.h>
#include <gsModeling/gsPeriodicParametrization.h>

namespace gismo
{

/**
 * A class for computing periodic parametrizations of closed
 * (cylinder-like) surface meshes. The result will be periodic in the
 * u-direction and the parameter domain will be [0, 1]^2. This class is an alternative to
 * gsPeriodicOverlap. The idea is adapted from
 *
 * Tong, Y., Alliez, P., Cohen-Steiner, D., Desbrun, M.: Designing
 * quadrangulations with discrete harmonic forms, in: Sheffer, A.,
 * Polthier, K. (Eds.), Symposium on Geometry Processing,
 * Eurographics. pp. 201–210, 2006,
 * 
 * where it was used constructing discrete harmonic mappings on
 * arbitrary topology.
 *
 * Main idea
 * =========
 *
 * The user prescribes the u-values of the points on the lower and
 * upper boundaries (i.e., v = 0 and v = 1, respectively). In
 * addition, they are required to provide a *stitch*: a list of
 * vertices forming a simple polyline connecting the vertices with the
 * highest u-value on v=0 and v=1, which serves as a first guess for
 * the location of the periodic interface.
 *
 * The equations for the inner vertices (depicted with empty circles
 * in the figure below) are assembled by the methods of the parent
 * class. For vertices on the stitch (full squares) the u-coordinate
 * of their neighbours to the right (empty squares) is temporarily
 * increased by one; similarly, when computing the empty-squared
 * neighbours, the u-coordinate from the full-square neighbours is
 * decreased by one.
 * 
 * @image html  gsPeriodicParametrizationStitch.png
 *
 * @image latex gsPeriodicParametrizationStitch.pdf
 
 */
template <class T>
class GISMO_EXPORT gsPeriodicStitch : public gsPeriodicParametrization<T>
{

public:

    typedef memory::shared_ptr<gsPeriodicStitch<T> > uPtr;

protected:

    /**
     * Modification of the corresponding class from
     * gsParametrization<T>. Given the indices of the stitch vertices,
     * a matrix of corrections is produced according to vertices being
     * neighbours across the interface.
     */
    class Neighbourhood : public gsParametrization<T>::Neighbourhood
    {
    public:
        typedef typename gsParametrization<T>::LocalNeighbourhood LocalNeighbourhood;

        /** Constructor.
         * @param meshInfo: surface mesh (as in the parent class)
         * @param stitchIndices: indices of the vertices forming the stitch
         * @param[out] corrections: a reference to @a m_corrections of
         * gsPeriodicStitch that gets filled here
         * @param parametrizationMethod: parametrization method (as in the parent class)
         */
        explicit Neighbourhood(const gsHalfEdgeMesh<T> &meshInfo,
                               const std::vector<size_t>& stitchIndices,
                               gsSparseMatrix<int>& corrections,
                               const size_t parametrizationMethod = 2);

    protected:
        /**
         * For a stitch vertex finds its neighbours across the interface.
         * @param stitchIndices indices of the stitch vertices
         * @param localNeighbourhood local neighbourhood of the vertex in question
         */
        std::vector<size_t> computeCorrections(const std::vector<size_t>& stitchIndices,
                                               const LocalNeighbourhood& localNeighbourhood) const;
    };

public:

    /** Constructor
     * @param mesh the surface mesh to be parametrized
     * @param verticesV0 matrix with three rows containing the vertices with v=0
     * @param paramsV0 matrix with one row containing the u-parameters of vertices with v=0
     * @param verticesV1 matrix with three rows containing the vertices with v=1
     * @param paramsV1 matrix with one row containing the u-parameters of vertices with v=1
     * @param stitchVertices matrix with three rows containing the vertices on the stitch
     * @param list list of the method options
     */
    explicit gsPeriodicStitch(const gsMesh<T> &mesh,
                              const gsMatrix<T>& verticesV0,
                              const gsMatrix<T>& paramsV0,
                              const gsMatrix<T>& verticesV1,
                              const gsMatrix<T>& paramsV1,
                              const gsMatrix<T>& stitchVertices,
                              const gsOptionList &list = gsParametrization<T>::defaultOptions())
        : gsPeriodicParametrization<T>(mesh, verticesV0, paramsV0, verticesV1, paramsV1, list),
        m_stitchIndices(this->indices(stitchVertices))
    {
        // Note: m_corrections gets filled later.
    }

    /// Computes the periodic parametrization.
    void compute();

protected:
    /**
     * Calculation itself
     * @param paraMethod parametrization method (cf. gsParametrization<T>)
     */
    void calculate(const size_t paraMethod);

    /** Similar to @a constructAndSolveEquationSystem but works for
     * periodic meshes using the corrections.
     */
    void constructAndSolveEquationSystem(const Neighbourhood &neighbourhood,
                                         const size_t n,
                                         const size_t N);

public:

    /**
     * Creates a flat mesh out of a periodic parametrization created by a the stitch method.
     * The mesh is then restricted to [0, 1]^2.
     */
    gsMesh<T> createFlatMesh() const
    {
        gsMesh<T> unfolded = createUnfoldedFlatMesh();
        typename gsPeriodicParametrization<T>::FlatMesh display(unfolded);
        return display.createRestrictedFlatMesh();
    }

protected:
    /**
     * Creates an unfolded mesh in the sense that the vertices on the
     * stitch are present twice, once in their original parameters and
     * once with u decreased by one.
     */
    gsMesh<T> createUnfoldedFlatMesh() const;

    /**
     * Query, whether the vertex with index @a vertexIndex (in the
     * numbering of Floater's paper, i.e., starting from 1) is on the
     * stitch or not.
     */
    bool isOnStitch(size_t vertexIndex) const
    {
        return (std::find(m_stitchIndices.begin(), m_stitchIndices.end(), vertexIndex)
                != m_stitchIndices.end());
    }

    /**
     * Query, whether an edge between vertices with indices @a beg and
     * @a end corrects a stitch vertex with its neighbour to the right
     * (i.e., whether it connects and empty square and full square
     * vertex in the figure documenting this class).
     */
    bool edgeIsInCorrections(index_t beg, index_t end) const
    {
        return ((m_corrections(beg, end) ==  1) ||
                (m_corrections(beg, end) == -1) ||
                (m_corrections(end, beg) ==  1) ||
                (m_corrections(end, beg) == -1));
        // Actually, the first two conditions should be enough.         
    }

protected:

    /**
     * Slot (i, j) is equal to 1 iff vertex with id i-1 is on a stitch
     * and vertex with id j-1 is among its neighbours across the
     * interface. Slot (j, i) is then set to -1.
     */
    gsSparseMatrix<int> m_corrections;

    const std::vector<size_t> m_stitchIndices; ///< indices of the vertices on the stitch
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPeriodicStitch.hpp)
#endif
