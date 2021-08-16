/** @file gsDPatchBasis.h

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

template<short_t d,class T>
class gsDPatchBasis : public gsMappedBasis<d,T>
{
    protected:
        using gsMappedBasis<d,T>::m_mapper;
        using gsMappedBasis<d,T>::m_topol;
        using gsMappedBasis<d,T>::m_bases;

    public:
        using gsMappedBasis<d,T>::size;
        using gsMappedBasis<d,T>::degree;
        using gsMappedBasis<d,T>::nPatches;

        typedef std::vector<gsBasis<T> *> BasisContainer;


    public:

        /// Shared pointer for gsDPatchBasis
        typedef memory::shared_ptr< gsDPatchBasis > Ptr;

        /// Unique pointer for gsDPatchBasis
        typedef memory::unique_ptr< gsDPatchBasis > uPtr;

        /// Empty constructor
        gsDPatchBasis() : m_patches(gsMultiPatch<T>())
        { }

        /// Default constructor
        gsDPatchBasis(gsMultiPatch<T> const & mp) ;

        /// Copy constructor
        gsDPatchBasis( const gsDPatchBasis& other );

        GISMO_CLONE_FUNCTION(gsDPatchBasis)

        virtual ~gsDPatchBasis();

        void verbose() { m_verbose = true;}

        /**
         * @brief       Returns the basis that is used for the D-Patch. Could be THB refined.
         *
         */
        gsMultiBasis<T> localBasis() const {return m_multiBasis;}

        /**
         * @brief       Returns the basis that is used for the D-Patch. Could be THB refined.
         *
         */
        gsMultiBasis<T> makeGlobalBasis();

        std::vector<gsBasis<T>*> container()
        {
            return m_basisContainer;
        }

        const gsMultiBasis<T> & multiBasis()
        {
            return m_globalBasis;
        }

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
        const index_t indexFromSides(index_t index1, const patchSide side1, index_t index2, const patchSide side2);


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
        const gsVector<index_t> indicesFromVert(index_t index, const patchCorner corner, const patchSide side, index_t offset = 0);


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
        const index_t indexFromVert(gsMultiBasis<T> bases, index_t index, const patchCorner corner, const patchSide side, index_t offset = 0, index_t levelOffset = 0);
        const index_t indexFromVert(index_t index, const patchCorner corner, const patchSide side, index_t offset = 0, index_t levelOffset = 0);


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
        const std::vector<index_t> indexFromVert(gsMultiBasis<T> bases, std::vector<index_t> index, const patchCorner corner, const patchSide side, index_t offset = 0, index_t levelOffset = 0);
        const std::vector<index_t> indexFromVert(std::vector<index_t> index, const patchCorner corner, const patchSide side, index_t offset = 0, index_t levelOffset = 0);


        /**
         * @brief      Computes the index of a basis function in a \ref patchCorner
         *
         * @param[in]  corner  The corner
         *
         * @return     Local index of the basis function
         */
        const index_t indexInCorner(const patchCorner corner)
        {return indexFromVert(0,corner,patchSide(corner.patch,1)); }

        /**
         * @brief      Returns the valence and whether a corner is interior or boundary
         *
         * @param[in]  corner  The \ref patchCorner
         *
         * @return     A pair with .first giving the valence and .second being true if the vertex is interior and false if the vertex is on a boundary
         */
        const std::pair<index_t,bool> vertexData(const patchCorner corner) const;

        /**
         * @brief      Gets the valence.
         *
         * @param[in]  corner  The corner
         *
         * @return     The valence.
         */
        const index_t getValence( patchCorner corner) const
        { return this->vertexData(corner).first; }

        /**
         * @brief      Determines whether the specified corner is interior vertex.
         *
         * @param[in]  corner  The corner
         *
         * @return     True if the specified corner is interior vertex, False otherwise.
         */
        const bool isInteriorVertex( patchCorner corner) const
        { return this->vertexData(corner).second; }

        /**
         * @brief      Computes global index of the side
         *
         * @param[in]  patch    The patch number
         * @param[in]  bside    The \ref boxSide
         *
         * @return     Returns a global index of the side
         */
        const index_t sideIndex( index_t patch,  boxSide bside)     const
        { return 4*patch + bside - 1; }
        /**
         * @brief      Computes global index of the side
         *
         * @param[in]  pside    The \ref patchSide
         *
         * @return     Returns a global index of the side
         */
        const index_t sideIndex( patchSide pside)     const
        { return sideIndex( pside.patch , pside.side() ); }

        /**
         * @brief      Computes global index of the corner
         *
         * @param[in]  patch    The patch number
         * @param[in]  corner   The \ref boxCorner
         *
         * @return     Returns a global index of the corner
         */
        const index_t vertIndex( index_t patch,  boxCorner corner)  const
        { return 4*patch + corner -1; }

        /**
         * @brief      Computes global index of the corner
         *
         * @param[in]  pcorner   The \ref patchCorner
         *
         * @return     Returns a global index of the side
         */
        const index_t sideIndex( patchCorner pcorner)     const
        { return vertIndex( pcorner.patch , pcorner.corner() ); }


    protected:
        /**
         * @brief      Makes a thb.
         */
        void makeTHB();

        /**
         * @brief      Calculates the d patch.
         */
        void computeDPatch();

        /**
         * @brief      Makes a pi.
         *
         * @param[in]  valence  The valence
         *
         * @return     { description_of_the_return_value }
         */
        gsMatrix<T> makePi(index_t valence);

        /**
         * @brief      Initializes the object.
         */
        void initialize(); // also initialize the mappers!

        /**
         * @brief      Initializes the mapper.
         */
        void initializeMapper(); // also initialize the mappers!

        /**
         * @brief      Calculates the smooth matrix.
         */
        void computeSmoothMatrix();

        /**
         * @brief      { function_description }
         *
         * @param[in]  pcorner  The pcorner
         */
        void handleVertex(patchCorner pcorner);
        /**
         * @brief      { function_description }
         *
         * @param[in]  iface  The interface
         */
        void handleInterface(boundaryInterface iface);
        /**
         * @brief      { function_description }
         *
         * @param[in]  side  The side
         */
        void handleBoundary(patchSide side);
        /**
         * @brief      { function_description }
         */
        void handleInterior();
        /**
         * @brief      { function_description }
         */
        void whichHandled();
        /**
         * @brief      Calculates the smooth matrix 2.
         */
        void computeSmoothMatrix2();

    protected:
        const gsMultiPatch<T> & m_patches;
        gsMultiPatch<T> m_RefPatches;
        gsMultiBasis<T> m_multiBasis, m_multiBasis0;
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

        #define PI 3.141592653589793


};



}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsDPatchBasis.hpp)
#endif
