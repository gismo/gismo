/** @file gsAlmostC1.h

    @brief Creates the D-Patch smoothing matrix.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsBoxTopology.h>
#include <gsCore/gsMultiPatch.h>

#include <gsMSplines/gsMappedBasis.h>

namespace gismo
{


/**
 * @brief      Constructs the D-Patch, from which the transformation matrix can be called
 *
 * @tparam     d     parametric dimension
 */
template<short_t d,class T>
class gsAlmostC1  //: public gsMappedGeom<d,T>
{
    protected:
        const gsMultiPatch<T> & m_patches;
        gsMultiPatch<T> m_RefPatches;
        gsMultiBasis<T> m_bases, m_Bbases;
        gsMultiBasis<T> m_globalBasis;
        gsMappedBasis<d,T> m_MBasis;
        std::vector<gsBasis<T> *> m_basisContainer;
        mutable gsSparseMatrix<T> m_tMatrix;
        mutable std::vector<boxCorner> m_bcorners;
        mutable std::vector<patchCorner> m_pcorners;
        mutable boxCorner m_bcorner;
        mutable patchCorner m_pcorner;
        mutable std::vector< std::vector<gsMatrix<index_t> > > m_sides; // structure: m_sides[s][0] is the 0th line of indices of interface s, m_sides[s][k] is the kth line of indices
        mutable std::vector<bool> m_sideCheck;
        mutable std::vector<bool> m_vertCheck;
        mutable std::vector<bool> m_basisCheck;

        mutable gsDofMapper m_mapModified,m_mapOriginal;

        mutable gsSparseMatrix<T> m_matrix;

        mutable size_t m_size;

        mutable std::vector<std::pair<patchCorner,index_t>> m_bVertices;
        mutable std::vector<std::pair<patchCorner,index_t>> m_iVertices;
        // mutable std::vector<patchSide> m_boundaries;
        // mutable std::vector<patchSide> m_interfaces;

        mutable gsMatrix<T> m_coefs;

        mutable bool m_verbose;

        std::vector<patchCorner> m_C0s;

        #define PI 3.141592653589793

        // This will store the triangles as soon as they are computed
        std::map<patchCorner,gsMatrix<T,3,3>> m_triangles;

    public:

        /// Shared pointer for gsAlmostC1
        typedef memory::shared_ptr< gsAlmostC1 > Ptr;

        /// Unique pointer for gsAlmostC1
        typedef memory::unique_ptr< gsAlmostC1 > uPtr;

        /// Empty constructor
        gsAlmostC1() : m_patches(gsMultiPatch<T>())
        { }


        /**
         * @brief      Default constructor
         *
         * @param      mp    Multipatch of the geometry
         */
        gsAlmostC1(gsMultiPatch<T> const & mp) ;

        /// Copy constructor
        gsAlmostC1( const gsAlmostC1& other );

        GISMO_CLONE_FUNCTION(gsAlmostC1)

        virtual ~gsAlmostC1();

        /**
         * @brief      Allow verbosity
         */
        void verbose() { m_verbose = true;}

        /**
         * @brief       Returns the basis that is used for the D-Patch. Could be THB refined.
         *
         */
        gsMultiBasis<T> localBasis() const {return m_bases;}

        /**
         * @brief       Returns the basis on which the D-Patch is applied
         *
         */
        gsMultiBasis<T> makeGlobalBasis();

        /**
         * @brief       Returns the multipatch that is used for the D-Patch
         *
         */
        gsMultiPatch<T> getGeometry() const {return m_patches;}

        /**
         * @brief       Returns for each basis function if it is free or eliminated
         *
         * Returns for each basis function if it is free or eliminated and checks if the internal mapper is defined correctly
         */
        void mapperInfo() const;

        /**
         * @brief      Returns information about a vertex
         *
         * @param[in]  corner  The \ref patchCorner
         *
         * @return     Prints the patch number, the corner index, the valence and if the corner is an interior or a boundary vertex
         */
        const void vertexInfo(patchCorner corner) const;

        /**
         * @brief      Returns information about a vertex
         *
         * @param[in]  patch  The \ref patchSide
         *
         * @return     Prints the patch number, the side index, the valence and if the side is a boundary or an interface
         */
        const void sideInfo(patchSide side) const;

        /**
         * @brief       Returns information for all the sides in the topology.
         *
         * Returns for all the patches and for all sides (1 to 4) if it is a boundary or an interface.
         */
        const void sideInfo() const;

        /**
         * @brief       Returns information for all the corners in the topology.
         *
         * Returns for all the patches and for all corners (1 to 4) the valence and if it is an interior vertex or a boundary vertex.
         */
        const void cornerInfo() const;

        /**
         * @brief       Computes the C1 coefficients for pre-multiplication to make the multipatch
         *
         * Takes the coefficients which are tagged as "free" in the modified DoFMapper (m_mapModified) and when a boundary vertex with valence=3 is present, this one is shifted.
         *
         */
        gsMatrix<T> preCoefficients();

        /**
         * @brief       Computes the local coefficients and puts them in one big matrix
         */
        gsMatrix<T> allCoefficients() const;

        /**
         * @brief       Exports a single modified patch with index \a patch
         *
         * The patch is obtained by transforming the coefficients of the D-Patch to the original basis, such that the original basis functions can be used to plot the geometry (and the patch structure will remain intact).
         * To construct the geometry, the coefficients for the C1 basis are multiplied with the transpose of the transformation matrix. The C1 coefficients are obtained with \ref preCoefficients().
         *
         */
        gsGeometry<T>* exportPatch(index_t patch, bool computeCoefs=true);

        /**
         * @brief      Exports the modified geometry to a @a gsMultiPatch object
         *
         * @return     A multipatch with the geometry
         */
        gsMultiPatch<T> exportToPatches();

        /**
         * @brief      Returns the smoothing matrix into \a matrix
         *
         * @param      matrix  The matrix
         */
        const void matrix_into(gsSparseMatrix<T> & matrix) const
        { matrix = m_matrix; }

        /**
         * @brief      Returns the smoothing matrix
         *
         * The number of columns of the matrix corresponds to the number of basis functions in the local basis; this is the sum of all the basis functions over all the patches.
         * The number of rows of the matrix corresponds to the number of global basis functions, i.e. the number of basis functions corresponding to the D-Patch.
         * Multiplying the basis with the local basis function values gives the values of the basis functions in the global basis.
         *
         * @return     A matrix \a result to transfer local to global coefficients
         */
        const gsSparseMatrix<T> matrix() const
        {
            gsSparseMatrix<T> result; matrix_into(result);
            return result;
        }


    protected:

        gsMatrix<T> _getNormals(const std::vector<patchCorner> & corners) const;

        void _toBarycentricCoordinates(const gsMatrix<T> & Cs, gsMatrix<T> & u) const;

        std::tuple<gsMatrix<T>,gsMatrix<T>,gsMatrix<index_t>> _makeTriangle(const patchCorner & corner) const;

        gsMatrix<T,3,3> _getRotationMatrix(const gsVector<T,3> & a, const gsVector<T,3> & b) const;

        /**
         * @brief      Computes the index of a basis function using sides as reference
         *
         * @param[in]  index1  The index of the basis function parallel to the first side
         * @param[in]  side1   The first side
         * @param[in]  index2  The index of the basis function parallel to the second side
         * @param[in]  side2   The second side
         *
         * @return     Index that is \a index1 in direction of \a side1 and \a index2 in direction of \a side2
         */
        const index_t _indexFromSides(index_t index1, const patchSide side1, index_t index2, const patchSide side2);


        /**
         * @brief      Computes the index of a basis function taking one corner and one side as reference
         *
         * @param[in]  index   Offset of the basis function parallel to the side \a side, measured from \a corner
         * @param[in]  corner  The corner to be measured from
         * @param[in]  side    The side which contains \a corner
         * @param[in]  offset  The offset from the side (orthogonal to the side)
         *
         * @return     Index of \a index places from \a corner along \a side, with offset \a offset
         */
        const gsVector<index_t> _indicesFromVert(index_t index, const patchCorner corner, const patchSide side, index_t offset = 0);


        /**
         * @brief      Computes the index of a basis function taking one corner and one side as reference
         *
         * @param[in]  bases   (optional) Multibasis to evaluate the index on
         * @param[in]  index   Offset of the basis function parallel to the side \a side, measured from \a corner
         * @param[in]  corner  The corner to be measured from
         * @param[in]  side    The side which contains \a corner
         * @param[in]  offset  The offset from the side (orthogonal to the side)
         * @param[in]  levelOffset  The level to be computed from. \a levelOffset = 0 returns the deepest THB level, and any positive number will give one level coarser
         *
         * @return     Index of \a index places from \a corner along \a side, with offset \a offset and with offset \a levelOffset from the deepest level
         */
        const index_t _indexFromVert(gsMultiBasis<T> bases, index_t index, const patchCorner corner, const patchSide side, index_t offset = 0, index_t levelOffset = 0) const;
        const index_t _indexFromVert(index_t index, const patchCorner corner, const patchSide side, index_t offset = 0, index_t levelOffset = 0) const;


        /**
         * @brief      Computes the index of a basis function taking one corner and one side as reference (multiple indices)
         *
         * @param[in]  bases   (optional) Multibasis to evaluate the index on
         * @param[in]  index        Vector with offsets of the basis function parallel to the side \a side, measured from \a corner
         * @param[in]  corner       The corner
         * @param[in]  side         The side
         * @param[in]  offset       The offset
         * @param[in]  levelOffset  The level offset
         *
         * @return     { description_of_the_return_value }
         */
        const std::vector<index_t> _indexFromVert(gsMultiBasis<T> bases, std::vector<index_t> index, const patchCorner corner, const patchSide side, index_t offset = 0, index_t levelOffset = 0) const;
        const std::vector<index_t> _indexFromVert(std::vector<index_t> index, const patchCorner corner, const patchSide side, index_t offset = 0, index_t levelOffset = 0) const;


        /**
         * @brief      Computes the index of a basis function in a \ref patchCorner
         *
         * @param[in]  corner  The corner
         *
         * @return     Local index of the basis function
         */
        const index_t _indexInCorner(const patchCorner corner)
        {return _indexFromVert(0,corner,patchSide(corner.patch,1)); }

        /**
         * @brief      Returns the valence and whether a corner is interior or boundary
         *
         * @param[in]  corner  The \ref patchCorner
         *
         * @return     A pair with .first giving the valence and .second being true if the vertex is interior and false if the vertex is on a boundary
         */
        const std::pair<index_t,bool> _vertexData(const patchCorner corner) const;

        /**
         * @brief      Gets the valence.
         *
         * @param[in]  corner  The corner
         *
         * @return     The valence.
         */
        const index_t _getValence( patchCorner corner) const
        { return this->_vertexData(corner).first; }

        /**
         * @brief      Determines whether the specified corner is interior vertex.
         *
         * @param[in]  corner  The corner
         *
         * @return     True if the specified corner is interior vertex, False otherwise.
         */
        const bool _isInteriorVertex( patchCorner corner) const
        { return this->_vertexData(corner).second; }

        /**
         * @brief      Computes global index of the side
         *
         * @param[in]  patch    The patch number
         * @param[in]  bside    The \ref boxSide
         *
         * @return     Returns a global index of the side
         */
        const index_t _sideIndex( index_t patch,  boxSide bside)     const
        { return 4*patch + bside - 1; }
        /**
         * @brief      Computes global index of the side
         *
         * @param[in]  pside    The \ref patchSide
         *
         * @return     Returns a global index of the side
         */
        const index_t _sideIndex( patchSide pside)     const
        { return _sideIndex( pside.patch , pside.side() ); }

        /**
         * @brief      Computes global index of the corner
         *
         * @param[in]  patch    The patch number
         * @param[in]  corner   The \ref boxCorner
         *
         * @return     Returns a global index of the corner
         */
        const index_t _vertIndex( index_t patch,  boxCorner corner)  const
        { return 4*patch + corner -1; }

        /**
         * @brief      Computes global index of the corner
         *
         * @param[in]  pcorner   The \ref patchCorner
         *
         * @return     Returns a global index of the side
         */
        const index_t _vertIndex( patchCorner pcorner)     const
        { return _vertIndex( pcorner.patch , pcorner.corner() ); }

        void _getLowestCorners(std::vector<patchCorner> & pcorners, index_t n = 3) const;
        void _removeLowestCorners(std::vector<patchCorner> & pcorners, index_t n = 3) const;

    protected:
        /**
         * @brief      Prepares the THB basis if needed.
         *
         * This function constructs THB refinements on the places where they are needed, i.e. around EVs. It also constructs the transfer matrix (m_tMatrix) forms the transformation between the original B-spline basis and the THB-Spline basis.
         */
        void _makeTHB();

        /**
         * @brief      Computes D-Patch smoothing
         *
         * Given a basis with THB refinement around the EVs, this function computes the D-Patch smoothing
         */
        void _computeDPatch();

        /**
         * @brief      Makes the Pi matrix
         *
         * This matrix is used to transform the coefficients of the D-Patch smoothing matrix
         *
         * @param[in]  valence  The valence
         *
         * @return     Matrix for smoothing around an EV}
         */
        gsMatrix<T> _makePi(index_t valence);

        /**
         * @brief      Initializes the matrix, the basis and the mappers
         */
        void _initialize();

        /**
         * @brief      Computes the modified mapper
         *
         * The modified mapper is computed based on the elimination of different functions with different conditions.
         * 1) For interfaces, it eliminates all the nodes except the first two and the last two
         * 2) For boundaries, there is no elimination
         * 3) For vertices, there are few options
         *  a) Boundary vertices
         *      i)  Valence 1: No eliminations
         *      ii) Valence 2: the two outer-most basis functions on the interface (the one at the vertex and the one next to it on the interface) are both eliminated
         *      iii)Valence 3: On all the patches, the basis functions corresponding to the vertex are coupled to eachother. The basis functions next to this one (on an interface OR on a boundary) are eliminated
         *  b) Interior vertices: all basis functions along the interface are eliminated if not done so
         */
        void _computeMapper(); // also initialize the mappers!

        /**
         * @brief      Handles a vertex in the global matrix
         *
         * We use the following notation convention (per patch!):
         * b00 is the basis function at the vertex
         * b10 is the basis function next to the vertex along the first interface that connects to the vertex
         * b20 is the basis function next to b10 along the first interface that connects to the vertex
         * etc.
         *
         * b01 is the basis function next to the vertex along the second interface that connects to the vertex
         * b02 is the basis function next to b01 along the second interface that connects to the vertex
         * etc.
         *
         * b11 is the basis function with offset 1 from both interfaces and from the vertex itself
         * b22 is the basis function with offset 2 from both interfaces and from the vertex itself
         * etc.
         *
         * There are different options.
         * a) Boundary vertices
         *      i)  Valence 1: b00, b10, b01 and b00 all get weight 1.0 w.r.t the same basis function in the local basis
         *      ii) Valence 2: This case contains an interface between two patches. We use index k to denote the row basis functions along the interface. So k=0 corresponds to the basis functions on the boundary and k=1 corresponds to the basis functions with offset 1 from the boundaries. Using this convention, the functions bk1 in the local basis, are coupled to bk1 in the global basis with weight 1. The functions bk0 in the local basis (on the considered patch) are coupled to bk1 in the global basis with weight 0.5. The functions bk0 in the local basis (on the other patch) are coupled to bk1 (on the considered patch) in the global basis with weight 0.5.
         *      iii)Valence 3: In this case, the matched vertices on all the adjacent patches are treated in one go! Note that all the basis functions corresponding to the vertex (b00) on all patches are matched! We couple the b00 functions of all patches (in the local basis) with weight 1/4 to the b00 of the adjacent patch with the lowest number in the global basis. Then, the b11 on the considered patch is coupled with weight 1 to itself and with weight 0.25 to the b00s of the other patches. Then, we will handle the vertices where an interface and a boundary meet (there are two of these). For the patch corners that are on an interface, we find the b11 and b10 vertices (orthogonal to the interface) and we give all b10s weight 0.5 w.r.t. the b11s in the global basis (on both patches). Lastly, we add weight 0.5 for the b10s along the boundaries (so only for two patches) to the (matched) b00 basis function (all b00s refer to the same dof in the global basis).
         * b) Interior vertices (all valences):
         *      i)  b11 gets weight 1.0 w.r.t the same basis function in the local basis
         *      ii) all associated b00s in the local basis get weight 1/valence w.r.t. b11 in the global basis
         *      iii)the b10 and b01 in the local basis get weight 1/2 w.r.t. b11 in the global basis
         *
         * @param[in]  pcorner  The patchcorner
         */
        void _handleVertex(patchCorner pcorner);
        /**
         * @brief      Handles an interface in the global matrix
         *
         * Gives all the DoFs that have offset 1 (orthogonal) from the interface weight 1.0 w.r.t itself. All the DoFs ON the interface (on both patches) will have weight 0.5 to the DoF with offset 1.
         * This interface handling excludes the indices that are in the 0 and 1 ring around vertices.
         *
         * @param[in]  iface  The interface
         */
        void _handleInterface(boundaryInterface iface);
        /**
         * @brief      Handles a boundary in the global matrix
         *
         * Handles all DoFs on the boundary with unit-weight, except the ones in the 0 and 1 rings around the vertices.
         *
         * @param[in]  side  The boundary side
         */
        void _handleBoundary(patchSide side);
        /**
         * @brief      Handles the interior in the global matrix
         *
         * Gives all left-over DoFs, which are in the interior, weight 1 w.r.t. itself
         */
        void _handleInterior();
        /**
         * @brief      Prints which DoFs have been handled and which have been eliminated
         */
        void _whichHandled();
        /**
         * @brief      Calculates the smoothing matrix.
         */
        void _computeSmoothMatrix();

        // LEGACY
        /*
            void _computeSmoothMatrix2();
        */


};



}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAlmostC1.hpp)
#endif
