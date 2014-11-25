/** @file gsDofMapper.h
    
    @brief Provides the gsDofMapper class for re-indexing DoFs.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Bressan, C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsBoundary.h>
#include <gsCore/gsExport.h>

namespace gismo
{

/** @brief
    Maintains a mapping from patch-local dofs to global dof indices
    and allows the elimination of individual dofs.

    A \em dof (degree of freedom) is, roughly speaking, an unknown in a
    discretization of a PDE. However, some dofs may be eliminated before
    solving the system and won't actually translate into unknowns.
    An example are dofs on Dirichlet boundaries.

    This class creates a mapping between an arbitrary number of per-patch
    local dofs to an enumeration of global dofs.
    Furthermore, dofs can also be marked as eliminated.

    Every global dof gets a unique number, forming
    a contiguous range starting from 0. This range has length gsDofMapper::size().
    The dofs are numbered in the following order:
      - first the \em free (non-eliminated) dofs
        (number: gsDofMapper::freeSize()),
      - then the eliminated dofs
        (number: gsDofMapper::boundarySize()).

    The eliminated dofs have their own 0-based numbering. The index of an eliminated
    global dof in this numbering can be queried with gsDofMapper::bindex().

    The object must be finalized before it is used,
    i.e. gsDofMapper::finalize() has to be called once before use.

*/

#define MAPPER_PATCH_DOF(a,b) m_dofs[m_offset[b]+a]

class GISMO_EXPORT gsDofMapper
{
public:

    /// Default empty constructor
    gsDofMapper();

    /**
     * @brief construct a dof mapper that identifies the degrees
     * of freedom in the multibasis
     * and eliminates the degrees of freedom on the dirichlet boundary
     *
     * @param bases
     * @param dirichlet
     */
    template<class T>
    gsDofMapper(
        const gsMultiBasis<T>         &bases,
        const gsBoundaryConditions<T> &dirichlet
         ) : m_shift(0)
    {
        init(bases,dirichlet);
    }

    /**
     * @brief construct a dof mapper that identifies the degrees
     * of freedom in the multibasis
     *
     * @param bases
     * @param dirichlet
     */
    template<class T>
    gsDofMapper(
        const gsMultiBasis<T>         &bases
        ) : m_shift(0)
    {
        init(bases);
    }


    /**
     * @brief construct a dof mapper that manages the the degrees
     * of freedom in a gsBasis
     *
     * @param basis
     */
    template<class T>
    gsDofMapper(
        const gsBasis<T>               &basis
         ) : m_shift(0)
    {
        initSingle(basis);
    }

    /// Initialize by a gsMultiBasis
    template <typename T>
    void init( const gsMultiBasis<T> & bases);

    /// Initialize by gsMultiBasis andboundary conditions
    template<class T>
    void init(
        const gsMultiBasis<T>         &basis,
        const gsBoundaryConditions<T> &dirichlet );

private:

    /// Initialize by a single basis patch
    template <typename T>
    void initSingle( const gsBasis<T> & basis);

    /// Initialize by vector of DoF indices
    void initPatchDofs(const gsVector<index_t> & patchDofSizes);

public:

    /// Called to initialize the gsDofMapper with matching interfaces
    /// after m_bases have already been set
    void setMatchingInterfaces(const gsBoxTopology & mp);


    /** \brief Calls matchDof() for all pairs of dofs on the given interface \a i between two patches.
     *
     * \warning Orientation is currently not implemented for 3D.
     */
    // 
    void matchInterface( boundaryInterface const & i );

    void matchInterface( index_t k1, index_t k2,
                         const gsMatrix<unsigned> & b1, 
                         const gsMatrix<unsigned> & b2, 
                         const gsVector<bool> &orient );

    /** \brief Calls matchDof() for all dofs on the given patch side
     * \a i ps. Thus, the whole set of dofs colapses to a single global dof
     *
     */
    void colapseDofs(index_t k, const gsMatrix<unsigned> & b );

    /** \brief Calls matchDof() for all dofs on the given patch sides
     * \a i psList. Thus, the whole side colapses to a single global dof
     *
     */
    void colapseSides( std::vector<patch_side> const & psList );

  

    /// Couples dof \a i of patch \a u with dof \a j of patch \a v such that they refer to the same global dof.
    void matchDof( index_t u, index_t i, index_t v, index_t j );


    // to do: put k at the end
    void markBoundary( index_t k, const gsMatrix<unsigned> & boundaryDofs );
    
    /// Mark the local dof \a i of patch \a k as eliminated.
    void eliminateDof( index_t i, index_t k );

    /// Must be called after all boundaries and interfaces have been marked to set up the dof numbering.
    void finalize();

    /// Print summary to cout
    void print() const;

    /// Set to identity
    void setIdentity(index_t nPatches, size_t nDofs);

    /// Add shifts
    void setShift(index_t shift);

    /** \brief Computes the global indices of the input local indices
     *
     * \param[in] locals a column matrix with the local indices
     * \param[in] patchIndex the index of the patch where the local indices belong to
     * \param[out] globals the global indices of the patch
     */
    void localToGlobal(const gsMatrix<unsigned>& locals,
                       index_t patchIndex,
                       gsMatrix<unsigned>& globals) const;


    /** \brief Returns the index associated to local dof \a i of patch \a k without shifts.
     *
     * \note This method only works after all interfaces and boundaries have been
     * marked and finalize() has been called.
     */

    inline index_t freeIndex( index_t i, index_t k = 0 ) const
    {
        GISMO_ASSERT(m_curElimId==0, "finalize() was not called on gsDofMapper");
        return MAPPER_PATCH_DOF(i,k);
    }

    /** \brief Returns the global dof index associated to local dof \a i of patch \a k.
     *
     * \note This method only works after all interfaces and boundaries have been
     * marked and finalize() has been called.
     */

    inline index_t index( index_t i, index_t k = 0 ) const
    {
        GISMO_ASSERT(m_curElimId==0, "finalize() was not called on gsDofMapper");
        return MAPPER_PATCH_DOF(i,k)+m_shift;
    }

    /// @brief Returns the boundary index of local dof \a i of patch \a k.
    ///
    /// Produces undefined results if local dof (i,k) does not lie on the boundary.
    inline index_t bindex( index_t i, index_t k = 0 ) const 
    { return global_to_bindex( index(i, k) );}

    /// @brief Returns the boundary index of global dof \a gl.
    ///
    /// Produces undefined results if dof \a gl does not lie on the boundary.
    inline index_t global_to_bindex( index_t gl ) const
    {
        GISMO_ASSERT( is_boundary_index( gl ), 
                      "global_to_bindex(): specified dof is not on the boundary");
        return gl - m_shift - freeSize();
    }

    /// Returns true if global dof \a gl is not eliminated.
    inline bool is_free_index( index_t gl ) const
    { return gl < freeSize() + m_shift; }

    /// Returns true if local dof \a i of patch \a k is not eliminated.
    inline bool is_free( index_t i, index_t k = 0 ) const 
    { return is_free_index( index(i, k) ); }

    /// Returns true if global dof \a gl is eliminated
    inline bool is_boundary_index( index_t gl ) const
    { return gl >= freeSize() + m_shift; }

    /// Returns true if local dof \a i of patch \a k is eliminated.
    inline bool is_boundary( index_t i, index_t k = 0 ) const 
    {return is_boundary_index( index(i, k) );}

    /// Returns the total number of dofs (free and eliminated).
    inline index_t size() const 
    { 
        GISMO_ENSURE(m_curElimId==0, "finalize() was not called on gsDofMapper");
        return m_numFreeDofs + m_numElimDofs; 
    }

    /// Returns the number of free (not eliminated) dofs.
    inline index_t freeSize() const 
    { 
        GISMO_ENSURE(m_curElimId==0, "finalize() was not called on gsDofMapper");
        return m_numFreeDofs; 
    }

    /// Returns the number of eliminated dofs.
    inline index_t boundarySize() const 
    { 
        GISMO_ENSURE(m_curElimId==0, "finalize() was not called on gsDofMapper");
        return m_numElimDofs; 
    }

private:

    // replace all references to oldIdx by newIdx
    inline void replaceDofGlobally(index_t oldIdx, index_t newIdx);

    void mergeDofsGlobally(index_t dof1, index_t dof2);

// Data members
private:

    // m_dofs/m_patchDofs stores for each patch the mapping from local to global dofs.
    //
    // During setup, the entries have a different meaning:
    //   0        -- regular free dof
    //   negative -- an eliminated dof
    //   positive -- a coupling dof
    // For nonzero entries, the value is an id which identifies the eliminated/coupling
    // group of the dof. Dofs with the same id will get the same dof index in the
    // final numbering stage in finalize().

    // Representation as vector of vectors, m_patchDofs[k][i]
    // corresponds to patch k patch-local basis index i 
    // std::vector< std::vector<index_t> > m_patchDofs;

    // Representation as a single vector plus offsets for patch-local
    // indices
    std::vector<index_t>     m_dofs;
    std::vector<std::size_t> m_offset;

    index_t m_shift;

    index_t m_numFreeDofs;
    index_t m_numElimDofs;

    // used during setup: running id for current eliminated/coupling dof
    // After finalize() is called m_curElimId takes the value zero.
    index_t m_curElimId;
    index_t m_curCouplingId;

}; // class gsDofMapper

} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsDofMapper.hpp)
#include GISMO_HPP_HEADER(gsDofMapper.cpp)
#endif

