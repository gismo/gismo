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

#define MAPPER_PATCH_DOF(a,b,c) m_dofs[c][m_offset[b]+a]

/** @brief Maintains a mapping from patch-local dofs to global dof indices
    and allows the elimination of individual dofs.

    A \em dof (degree of freedom) is, roughly speaking, an unknown in a
    discretization of a PDE. However, some dofs may be eliminated before
    solving the system and won't actually translate into unknowns.
    An example are dofs on Dirichlet boundaries.

    This class creates a mapping between an arbitrary number of
    per-patch local dofs to an enumeration of global dofs.
    Furthermore, dofs can also be marked as eliminated.

    This is a a many-to-one mapping: many patch-local dofs are mapped
    to a single global dfo (index).
    Every global dof gets a unique number, forming a continuous range
    starting from 0. This range has length gsDofMapper::size().
    The dofs are numbered in the following order:

    - first the standard \em free (non-eliminated) dofs, ie. dofs that
      are not coupled on the boundary. For the standard dofs there is
      unique pre-image pair (patch,localdof).

    - then the \em free dofs which are coupled with other dofs For the
      coupled dofs there is a list of pre-image pairs of the form
      (patch,localdof).

      Upto here we get all the dofs which are \em free (number:
      gsDofMapper::freeSize() ). Then a final group follows:

    - then the dofs that are on Dirichlet boundaries (number:
      gsDofMapper::boundarySize()). These dofs might have a unique
      pre-image or not.

    The boundary (eg. eliminated) dofs have their own 0-based
    numbering. The index of an boundary global dof in this numbering
    can be queried with gsDofMapper::bindex().

    The object must be finalized before it is used,
    i.e. gsDofMapper::finalize() has to be called once before use.

    \ingroup Core

*/
class GISMO_EXPORT gsDofMapper
{
public:

    /// Default empty constructor
    gsDofMapper();

    /**
     * @brief construct a dof mapper that identifies the degrees
     * of freedom in the multibasis
     * and eliminates the degrees of freedom on the Dirichlet boundary
     *
     * @param bases
     * @param dirichlet
     * @param unk
     */
    template<class T>
    gsDofMapper(
        const gsMultiBasis<T>         &bases,
        const gsBoundaryConditions<T> &dirichlet,
        int unk = 0
        ) : m_shift(0), m_bshift(0)
    {
      init(bases, dirichlet, unk); //obsolete, one component
    }

    /**
     * @brief construct a dof mapper that identifies the degrees
     * of freedom in the multibasis
     *
     * @param bases
     */
    template<class T>
    gsDofMapper(const gsMultiBasis<T> & bases, index_t nComp = 1) 
      : m_shift(0), m_bshift(0)
    {
      init(bases, nComp);
    }
    
    /**
     * @brief construct a dof mapper that identifies the degrees
     * of freedom for a vector of multibasis
     *
     * @param bases
     */
    template<class T>
    gsDofMapper(
        std::vector<const gsMultiBasis<T> *> const & bases
        ) : m_shift(0), m_bshift(0)
    {
        init(bases);
    }


    /**
     * @brief construct a dof mapper that manages the degrees
     * of freedom in a gsBasis
     *
     * @param basis
     */
    template<class T>
    gsDofMapper(const gsBasis<T> & basis, index_t nComp = 1) 
      : m_shift(0), m_bshift(0)
    {
      initSingle(basis, nComp);
    }

    /**
     * @brief construct a dof mapper with a given number of dofs per
     * patch
     *
     * @param patchDofSizes
     */
    gsDofMapper(const gsVector<index_t> &patchDofSizes, index_t nComp = 1) 
      : m_shift(0), m_bshift(0)
    {
        initPatchDofs(patchDofSizes, nComp);
    }

    /// Initialize by a gsMultiBasis
    template <typename T>
    void init(const gsMultiBasis<T> & bases, index_t nComp = 1);

    /// Initialize by a vector of gsMultiBasis.
    template <typename T>
    void init( std::vector<const gsMultiBasis<T> *> const & bases);

    /// Initialize by gsMultiBasis, boundary conditions and the index
    /// of the unknown to be eliminated
    template<class T>
    void init(const gsMultiBasis<T>         &basis,
	      const gsBoundaryConditions<T> &dirichlet, int unk = 0);

    void swap(gsDofMapper & other)
    {
        m_dofs  .swap(other.m_dofs);
        m_offset.swap(other.m_offset);

        std::swap(m_shift      , other.m_shift);
        std::swap(m_bshift     , other.m_bshift);
        std::swap(m_numFreeDofs, other.m_numFreeDofs);
        std::swap(m_numElimDofs, other.m_numElimDofs);
        std::swap(m_numCpldDofs, other.m_numCpldDofs);
        std::swap(m_curElimId  , other.m_curElimId);
        std::swap(m_tagged     , other.m_tagged);
    }

private:

    /// Initialize by a single basis patch
    template <typename T>
    void initSingle( const gsBasis<T> & basis, index_t nComp = 1);

    /// Initialize by vector of DoF indices and dimension
    void initPatchDofs(const gsVector<index_t> & patchDofSizes,
		       index_t nComp = 1);

public:

    /// Returns a vector taking flat local indices to global
    gsVector<index_t> asVector(index_t comp = 0) const;

    /** \brief Returns a vector taking global indices to flat local

        Assumes that the mapper is a permutation
    */
    gsVector<index_t> inverseAsVector(index_t comp = 0) const;

    /// Called to initialize the gsDofMapper with matching interfaces
    /// after m_bases have already been set
    void setMatchingInterfaces(const gsBoxTopology & mp);

    /** \brief Calls matchDof() for all dofs on the given patch side
     * \a i ps. Thus, the whole set of dofs collapses to a single
     * global dof
     *
     */
    void colapseDofs(index_t k, const gsMatrix<unsigned> & b, index_t comp = 0);

    /// \brief Couples dof \a i of patch \a u with dof \a j of patch
    /// \a v such that they refer to the same global dof at component
    /// \a comp.
    void matchDof( index_t u, index_t i, index_t v, index_t j, index_t comp = 0);

    /// \brief Couples dofs \a b1 of patch \a u with dofs \a b2 of patch
    /// \a v one by one such that they refer to the same global dof.
    void matchDofs(index_t u, const gsMatrix<index_t> & b1,
                   index_t v, const gsMatrix<index_t> & b2,
		           index_t comp = 0);

    /// Mark the local dof \a i of patch \a k as coupled.
    void markCoupled(index_t i, index_t k, index_t comp = 0);

    /// Mark a local dof \a i of patch \a k as tagged
    void markTagged(index_t i, index_t k, index_t comp = 0);

    /// Mark all coupled dofs as tagged
    void markCoupledAsTagged();

    /// Mark the local dofs \a boundaryDofs of patch \a k as eliminated.
    // to do: put k at the end
    void markBoundary(index_t k, const gsMatrix<index_t> & boundaryDofs, index_t comp = 0);

    /// Mark the local dof \a i of patch \a k as eliminated.
    void eliminateDof(index_t i, index_t k, index_t comp = 0);

    /// \brief Must be called after all boundaries and interfaces have
    /// been marked to set up the dof numbering.
    void finalize();

    /// \brief Checks whether finalize() has been called.
    bool isFinalized() const { return m_curElimId>=0; }

    /// \brief Returns true iff the mapper is a permuatation
    bool isPermutation() const { return static_cast<size_t>(size())==mapSize(); }

    /// \brief Print summary
    std::ostream& print( std::ostream& os = gsInfo ) const;

    ///\brief Set this mapping to be the identity
    void setIdentity(index_t nPatches, size_t nDofs, size_t nComp = 1);

    ///\brief Set the shift amount for the global numbering
    void setShift(index_t shift);

    /// \brief Permutes the mapped free indices according to permutation, i.e.,  dofs_perm[idx] = dofs_old[permutation[idx]]
    ///
    /// \warning Applying a permutation makes the functions regarding coupled dofs (cindex, is_coupled_index,.. ) invalid.
    /// The dofs are still coupled, but you have no way of extracting them. If you need this functions, first call
    /// markCoupledAsTagged() and then use the corresponding functions for tagged dofs.
    void permuteFreeDofs(const gsVector<index_t>& permutation, index_t comp = 0);

    ///\brief Returns the smallest value of the indices for \a comp
    index_t firstIndex(index_t comp = 0) const
    { return m_numFreeDofs[comp] + m_numElimDofs[comp] + m_shift; }

    ///\brief Returns one past the biggest value of the free indices
    index_t lastIndex() const { return m_shift + freeSize(); }

    ///\brief Set the shift amount for the boundary numbering
    void setBoundaryShift(index_t shift);

    /** \brief Computes the global indices of the input local indices
     *
     * \param[in] locals a column matrix with the local indices
     * \param[in] patchIndex the index of the patch where the local indices belong to
     * \param[out] globals the global indices of the patch
     */
    void localToGlobal(const gsMatrix<index_t>& locals,
                       index_t patchIndex,
                       gsMatrix<index_t>& globals, 
		               index_t comp = 0) const;

    /** \brief Computes the global indices of the input local indices
     *
     * \param[in] locals a column matrix with the local indices
     * \param[in] patchIndex the index of the patch where the local indices belong to
     * \param[out] globals the local-global correspondance
     * \param[out] numFree the number of free indices in \a local
     */
    void localToGlobal2(const gsMatrix<index_t>& locals,
                        index_t patchIndex,
                        gsMatrix<index_t>& globals,
                        index_t & numFree, 
		                index_t comp = 0) const;

    /** \brief Returns the index associated to local dof \a i of patch \a k without shifts.
     *
     * \note This method only works after all interfaces and boundaries have been
     * marked and finalize() has been called.
     */

    inline index_t freeIndex(index_t i, index_t k = 0, index_t c = 0) const
    {
        GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return MAPPER_PATCH_DOF(i,k,c);
    }

 index_t componentOf(index_t gl) const
 {
     GISMO_ASSERT(m_curElimId>=0,"finalize() was not called on gsDofMapper");
     //index_t c = 1; // could do some binary search
     //in place
     //while (gl >= m_numFreeDofs[c]+m_numElimDofs[c]) { ++c; }
     //elim
     return (gl<m_numFreeDofs.back() ? 
       std::distance(m_numFreeDofs.begin(), std::upper_bound(m_numFreeDofs.begin(), m_numFreeDofs.end(), gl))
	     : std::distance(m_numElimDofs.begin(),std::upper_bound(m_numElimDofs.begin(), m_numElimDofs.end(), gl-m_numFreeDofs.back())) ) - 1;
     //while (gl >= m_numFreeDofs[c] + m_shift) { ++c; } return c-1;
 }

    /** \brief Returns the global dof index associated to local dof \a i of patch \a k.
     *
     * \note This method only works after all interfaces and boundaries have been
     * marked and finalize() has been called.
     */

    inline index_t index(index_t i, index_t k = 0, index_t c = 0) const
    {
        GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return MAPPER_PATCH_DOF(i,k,c)+m_shift;
    }

    /// @brief Returns the boundary index of local dof \a i of patch \a k.
    ///
    /// Produces undefined results if local dof (i,k) does not lie on the boundary.
    inline index_t bindex(index_t i, index_t k = 0, index_t c = 0) const
    {
        GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return MAPPER_PATCH_DOF(i,k,c) - m_numFreeDofs.back()
	  - m_numElimDofs[c] + m_bshift;
    }

    /// Returns true iff all DoFs are considered as free
    bool allFree() const
    { return m_numFreeDofs.back()+m_numElimDofs.back()==m_curElimId; }

    /// @brief Returns the coupled dof index
    inline index_t cindex(index_t i, index_t k = 0, index_t c = 0) const
    {
        GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return MAPPER_PATCH_DOF(i,k,c) - m_numFreeDofs[c+1] 
	  + m_numCpldDofs[c+1];
    }

    /// @brief Returns the tagged dof index
    inline index_t tindex(index_t i, index_t k = 0, index_t c = 0) const
    {
        GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return std::distance(m_tagged.begin(),std::lower_bound(m_tagged.begin(),m_tagged.end(),MAPPER_PATCH_DOF(i,k,c)));
    }

    /// @brief Returns the boundary index of global dof \a gl.
    ///
    /// Produces undefined results if dof \a gl does not lie on the boundary.
    inline index_t global_to_bindex(index_t gl) const
    {
        GISMO_ASSERT( is_boundary_index( gl ),
                      "global_to_bindex(): dof "<<gl<<" is not on the boundary");

	gl -= m_numFreeDofs.back() + m_shift;
	const index_t c = std::distance(m_numElimDofs.begin(),
        std::upper_bound(m_numElimDofs.begin(), m_numElimDofs.end(),
	gl)) -1;
	return gl - m_numElimDofs[c] + m_bshift;
    }

    /// Returns true if global dof \a gl is not eliminated.
    inline bool is_free_index(index_t gl) const
    {
      return gl < m_curElimId + m_shift;
    }

    /// Returns true if local dof \a i of patch \a k is not eliminated.
    inline bool is_free( index_t i, index_t k = 0, index_t c = 0) const
    { return is_free_index( index(i, k, c) ); }

    /// Returns true if global dof \a gl is eliminated
    inline bool is_boundary_index( index_t gl ) const
    { return gl >= m_numFreeDofs.back() + m_shift; }

    /// Returns true if local dof \a i of patch \a k is eliminated.
    inline bool is_boundary(index_t i, index_t k = 0, index_t c = 0) const
    {return is_boundary_index( index(i, k, c) );}

    /// Returns true if local dof \a i of patch \a k is coupled.
    inline bool is_coupled( index_t i, index_t k = 0, index_t c = 0) const
    { return  is_coupled_index( index(i, k, c) ); }

    /// Returns true if \a gl is a coupled dof.
    inline bool is_coupled_index(index_t gl) const
    {
      const index_t gc = componentOf(gl);
      const index_t vv = m_numFreeDofs[gc+1] + m_shift;
      return  (gl < vv && // is a free dof and
	      (gl + m_numCpldDofs[gc+1] + 1 > vv) );  // is not standard dof
    }

    /// Returns true if local dof \a i of patch \a k is tagged.
    inline bool is_tagged(index_t i, index_t k = 0, index_t c = 0) const
    { return  is_tagged_index( index(i, k, c) ); }

    /// Returns true if \a gl is a tagged dof.
    inline bool is_tagged_index(index_t gl) const
    {
        return std::binary_search(m_tagged.begin(),m_tagged.end(),gl);
    }

    /// Returns the total number of dofs (free and eliminated).
    inline index_t size() const
    {
        GISMO_ENSURE(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return freeSize() + boundarySize();
    }

    /// Returns the total number of dofs (free and eliminated).
    inline index_t size(index_t comp) const
    {
        GISMO_ENSURE(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return m_numFreeDofs[comp+1]-m_numFreeDofs[comp] 
	  + m_numElimDofs[comp+1]-m_numElimDofs[comp];
    }

    /// Returns the number of free (not eliminated) dofs.
    inline index_t freeSize() const
    {
      return m_curElimId;
    }

    inline index_t freeSize(index_t comp) const
    {
      return m_numFreeDofs[comp+1]-m_numFreeDofs[comp] + 
	(allFree() ? m_numElimDofs[comp+1]-m_numElimDofs[comp] : 0 );
    }

    /// Returns the number of coupled (not eliminated) dofs.
    index_t coupledSize() const;

    /// Returns the number of tagged (not eliminated) dofs.
    index_t taggedSize() const;

    /// Returns the number of eliminated dofs.
    inline index_t boundarySize() const
    {
      GISMO_ENSURE(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return m_numElimDofs.back();
    }

    index_t boundarySizeWithDuplicates() const;

    /// Returns the offset corresponding to patch \a k
    size_t offset(int k) const {return m_offset[k];}

    /// Returns the number of patches present underneath the mapper
    size_t numPatches() const {return m_offset.size();}

    /// \brief Returns the total number of patch-local degrees of
    /// freedom that are being mapped
    size_t mapSize() const
    {return m_dofs.size() * m_dofs.front().size();}

    size_t componentsSize() const {return m_dofs.size();}

    /// \brief Returns the total number of patch-local DoFs
    /// that live on patch \a k for component \a c
    size_t patchSize(const index_t k, const index_t c = 0) const
    {
        const size_t k1(k+1);
        GISMO_ASSERT(k1<=numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
        if ( 1==m_offset.size() ) return  m_dofs[c].size();
        else if ( k1==m_offset.size() ) return (m_dofs[c].size() - m_offset.back());
        else return (m_offset[k1]-m_offset[k]);
    }

    /// \brief For \a gl being a global index, this function returns a
    /// vector of pairs (patch,dof) that contains all the pairs which
    /// map to \a gl
    void preImage(index_t gl, std::vector<std::pair<index_t,index_t> > & result) const;

    /// \brief Produces the inverse of the mapping on patch \a k
    /// assuming that the map is invertible on that patch
    std::map<index_t,index_t> inverseOnPatch(const index_t k) const;

    /// \brief For \a gl being a global index, this function returns
    /// true whenever \a gl corresponds to patch \a k
    bool indexOnPatch(const index_t gl, const index_t k) const;

    /// \brief For \a n being an index which is already offsetted, it
    /// returns the global index where it is mapped to by the dof
    /// mapper.
    index_t mapIndex(index_t n) const
    {
        return m_dofs[n/m_dofs.front().size()]
            [n%m_dofs.front().size()] + m_shift;
    }
private:

    void finalizeComp(const index_t comp);

    // replace all references to oldIdx by newIdx
    inline void replaceDofGlobally(index_t oldIdx, index_t newIdx);
    inline void replaceDofGlobally(index_t oldIdx, index_t newIdx, index_t comp);

    void mergeDofsGlobally(index_t dof1, index_t dof2);
    void mergeDofsGlobally(index_t dof1, index_t dof2, index_t comp);

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

    // Representation of each component as a single vector plus
    // offsets for patch-local indices
    std::vector<std::vector<index_t> >  m_dofs;

    /// Offsets
    std::vector<size_t> m_offset;

    /// Shifting of the global index (zero by default)
    index_t m_shift;

    /// Shifting of the boundary index (zero by default)
    index_t m_bshift;

    /// Offsets of free dofs, nComp+1
    std::vector<index_t> m_numFreeDofs;
    /// Offsets of eliminated dofs, nComp+1
    std::vector<index_t> m_numElimDofs;
    /// Offsets of coupled dofs, nComp+1
    std::vector<index_t> m_numCpldDofs;

    // used during setup: running id for current eliminated dof
    // After finalize() is called m_curElimId takes positive value.
    index_t m_curElimId;

    /// Stores the tagged indices
    std::vector<index_t> m_tagged;

}; // class gsDofMapper

/// Print (as string) a dofmapper structure
inline std::ostream& operator<<( std::ostream& os, const gsDofMapper& b )
{
    return b.print( os );
}


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsDofMapper.hpp)
#endif
