/** @file gsDofMapper.h

    @brief Provides the gsDofMapper class for re-indexing DoFs.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <type_traits>

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsBoundary.h>
#include <gsCore/gsExport.h>

namespace gismo
{

#define MAPPER_PATCH_DOF(a,b,c) m_dofs[c][off(b,c)+a]

/** @brief Maintains a mapping from patch-local dofs to global dof indices
    and allows the elimination of individual dofs.

    A \em dof (degree of freedom) is, roughly speaking, an unknown in a
    discretization of a PDE. However, some dofs may be eliminated before
    solving the system and won't actually translate into unknowns.
    An example are dofs on Dirichlet boundaries.

    This class creates a mapping between an arbitrary number of
    per-patch local dofs to an enumeration of global dofs.
    Furthermore, dofs can also be marked as eliminated.

    This is a many-to-one mapping: many patch-local dofs are mapped
    to a single global dof (index).
    Every global dof gets a unique number, forming a continuous range
    of length gsDofMapper::size() that starts at firstIndex(0), i.e. at
    0 unless the numbering is shifted (see below).
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

    The global numbering can be offset with setShift() and addShift().
    A global index is a token valid only for the mapper that produced
    it: shifted numberings of different mappers overlap (in a block
    system, the eliminated block of one space overlaps the free block
    of the next). For a free dof, the shifted index is also its
    row/column in the block system of the assembler that set the
    shift; the index of an eliminated dof has no meaning outside its
    mapper.

    Every function returning a global index (index(), asVector(),
    mapIndex(), firstIndex(), lastIndex(), the keys of
    inverseOnPatch()) returns a shifted one, and every function taking
    a global index \a gl (componentOf(), preImage(), anyPreImage(),
    indexOnPatch(), is_free_index(), is_boundary_index(),
    is_coupled_index(), is_tagged_index(), global_to_bindex()) expects
    a shifted one. Exceptions:
    - anyPreImages() and inverseAsVector() return dense vectors whose
      position <tt>gl - firstIndex(0)</tt> holds the entry of \a gl;
    - bindex(), global_to_bindex(), cindex() and tindex() return ranks
      among the eliminated, coupled and tagged dofs respectively, which
      do not depend on the shift.

    \ingroup Core

*/
class GISMO_EXPORT gsDofMapper
{
public:

    /// Default empty constructor
    gsDofMapper();

    /// \brief Constructs a mapper from the per-component patch sizes.
    /// \param patchDofSizes nPatches x nComp matrix; entry (k,c) is the number of DoFs of component c on patch k.
    ///        Components may have different sizes. Requires nPatches>0, nComp>0 and no negative entry.
    ///
    /// Templated such that it is only selected when giving an argument of type
    /// gsMatrix<index_t>. Anything else, including things that could be converted to
    /// a gsMatrix, goes to the "standard" constructor taking a vector.
    template<class T, typename std::enable_if<std::is_same<T, index_t>::value, int>::type = 0>
    explicit gsDofMapper(const gsMatrix<T> & patchDofSizes)
    { initPatchDofs(patchDofSizes); }

    /// \brief Uniform case: every one of the \a nComp components has the patch sizes \a patchDofSizes.
    /// Equivalent to gsDofMapper(patchDofSizes.replicate(1,nComp)).
    gsDofMapper(const gsVector<index_t> & patchDofSizes, index_t nComp = 1);

    void swap(gsDofMapper & other)
    {
        m_dofs  .swap(other.m_dofs);
        m_offset.swap(other.m_offset);

        std::swap(m_nPatches   , other.m_nPatches);
        std::swap(m_shift      , other.m_shift);
        std::swap(m_numFreeDofs, other.m_numFreeDofs);
        std::swap(m_numElimDofs, other.m_numElimDofs);
        std::swap(m_numCpldDofs, other.m_numCpldDofs);
        std::swap(m_curElimId  , other.m_curElimId);
        std::swap(m_tagged     , other.m_tagged);
    }

private:

    /// Initialize from the nPatches x nComp size matrix
    void initPatchDofs(const gsMatrix<index_t> & patchDofSizes);

public:

    /// Returns a vector taking flat local indices to global
    gsVector<index_t> asVector(index_t comp = 0) const;

    /** \brief Returns a vector taking global indices to flat local

        Position \a i of the result corresponds to the global index
        <tt>i + firstIndex(0)</tt>, i.e., global indices minus the shift.
        The entries are the flat local indices of component \a comp, or -1
        at positions belonging to other components.

        Requires isPermutation().
    */
    gsVector<index_t> inverseAsVector(index_t comp = 0) const;

    /** \brief Calls matchDof() for all dofs on the given patch side
     * \a i ps. Thus, the whole set of dofs collapses to a single
     * global dof
     *
     */
    void colapseDofs(index_t k, const gsMatrix<unsigned> & b, index_t comp = 0);

    /// \brief Couples dof \a i of patch \a u with dof \a j of patch
    /// \a v such that they refer to the same global dof at component
    /// \a comp.
    ///
    /// With <tt>comp == -1</tt> the same local indices are matched in every
    /// component; this is only valid when the components share one basis.
    /// A <tt>comp == -1</tt> call throws (GISMO_ENSURE, all builds) if a
    /// patch involved has a different size in some component; call once per
    /// component then.
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

    /// Mark the local dofs \a boundaryDofs of patch \a k as eliminated (comp = -1: see eliminateDof()).
    // to do: put k at the end
    void markBoundary(index_t k, const gsMatrix<index_t> & boundaryDofs, index_t comp = 0);

    /// Mark the local dof \a i of patch \a k as eliminated.
    ///
    /// With <tt>comp == -1</tt> the same local index is eliminated in every
    /// component; this is only valid when the components share one basis.
    /// A <tt>comp == -1</tt> call throws (GISMO_ENSURE, all builds) if patch
    /// \a k has a different size in some component; call once per component
    /// then.
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

    ///\brief Add a shift amount to the global numbering
    void addShift(index_t shift);

    ///\brief Returns the smallest (shifted) global index of the free
    ///dofs of component \a comp. Valid after finalize().
    index_t firstIndex(index_t comp = 0) const
    { return m_numFreeDofs[comp] + m_shift; }

    ///\brief Returns one past the biggest value of the free indices
    index_t lastIndex() const { return m_shift + freeSize(); }

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

    index_t componentOf(index_t gl) const
    {
         GISMO_ASSERT(m_curElimId>=0,"finalize() was not called on gsDofMapper");
         GISMO_ASSERT(gl >= m_shift && gl < m_shift + m_numFreeDofs.back() + m_numElimDofs.back(),
                      "componentOf(): global index "<<gl<<" outside ["<<m_shift<<","
                      <<m_shift + m_numFreeDofs.back() + m_numElimDofs.back()<<")");
         //index_t c = 1; // could do some binary search
         //in place
         //while (gl >= m_numFreeDofs[c]+m_numElimDofs[c]) { ++c; }
         //elim
         const index_t g = gl - m_shift;
         return (g<m_numFreeDofs.back() ?
                std::distance(m_numFreeDofs.begin(), std::upper_bound(m_numFreeDofs.begin(), m_numFreeDofs.end(), g))
              : std::distance(m_numElimDofs.begin(),std::upper_bound(m_numElimDofs.begin(), m_numElimDofs.end(), g-m_numFreeDofs.back())) ) - 1;
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
    /// The boundary index is the rank of the dof among the eliminated
    /// dofs, 0..boundarySize()-1, independent of the shift.
    /// Produces undefined results if local dof (i,k) does not lie on the boundary.
    inline index_t bindex(index_t i, index_t k = 0, index_t c = 0) const
    {
        GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return MAPPER_PATCH_DOF(i,k,c) - m_numFreeDofs.back();
    }

    /// Returns true iff all DoFs are considered as free
    bool allFree() const
    { return m_numFreeDofs.back()+m_numElimDofs.back()==m_curElimId; }

    /// @brief Returns the rank of local dof (i,k,c) among the coupled dofs,
    /// 0..coupledSize()-1 counted over all components, independent of the shift.
    /// Produces undefined results if the dof is not coupled.
    inline index_t cindex(index_t i, index_t k = 0, index_t c = 0) const
    {
        GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return MAPPER_PATCH_DOF(i,k,c) - m_numFreeDofs[c+1]
	  + m_numCpldDofs[c+1];
    }

    /// @brief Returns the rank of local dof (i,k,c) among the tagged dofs,
    /// 0..taggedSize()-1, independent of the shift. If the dof is not
    /// tagged, returns the position at which it would be inserted.
    inline index_t tindex(index_t i, index_t k = 0, index_t c = 0) const
    {
        GISMO_ASSERT(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return std::distance(m_tagged.begin(),std::lower_bound(m_tagged.begin(),m_tagged.end(),MAPPER_PATCH_DOF(i,k,c)));
    }

    /// @brief Returns the boundary index of global dof \a gl.
    ///
    /// Produces undefined results if dof \a gl does not lie on the boundary.
    /// Same numbering as bindex(): bindex(i,k,c) == global_to_bindex(index(i,k,c)).
    inline index_t global_to_bindex(index_t gl) const
    {
        GISMO_ASSERT( is_boundary_index( gl ),
                      "global_to_bindex(): dof "<<gl<<" is not on the boundary");
        return gl - m_numFreeDofs.back() - m_shift;
    }

    /// Returns true if global dof \a gl is not eliminated.
    /// \a gl must be a global index of this mapper (checked in debug mode).
    inline bool is_free_index(index_t gl) const
    {
        GISMO_ASSERT(gl >= m_shift && gl < m_shift + m_numFreeDofs.back() + m_numElimDofs.back(),
                     "is_free_index(): global index "<<gl<<" outside ["<<m_shift<<","
                     <<m_shift + m_numFreeDofs.back() + m_numElimDofs.back()<<")");
        return gl < m_curElimId + m_shift;
    }

    /// Returns true if local dof \a i of patch \a k is not eliminated.
    inline bool is_free( index_t i, index_t k = 0, index_t c = 0) const
    { return is_free_index( index(i, k, c) ); }

    /// Returns true if global dof \a gl is eliminated.
    /// \a gl must be a global index of this mapper (checked in debug mode).
    inline bool is_boundary_index( index_t gl ) const
    {
        GISMO_ASSERT(gl >= m_shift && gl < m_shift + m_numFreeDofs.back() + m_numElimDofs.back(),
                     "is_boundary_index(): global index "<<gl<<" outside ["<<m_shift<<","
                     <<m_shift + m_numFreeDofs.back() + m_numElimDofs.back()<<")");
        return gl >= m_numFreeDofs.back() + m_shift;
    }

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
      const index_t cp = m_numCpldDofs[gc+1] - m_numCpldDofs[gc];
      return  (gl < vv && // is a free dof and
	      (gl + cp + 1 > vv) );  // is not standard dof
    }

    /// Returns true if local dof \a i of patch \a k is tagged.
    inline bool is_tagged(index_t i, index_t k = 0, index_t c = 0) const
    { return  is_tagged_index( index(i, k, c) ); }

    /// Returns true if \a gl is a tagged dof.
    inline bool is_tagged_index(index_t gl) const
    {
        return std::binary_search(m_tagged.begin(),m_tagged.end(),gl-m_shift);
    }

    /// Returns the number of components present in the mapper
    inline index_t numComponents() const
    { return static_cast<index_t>(m_dofs.size()); }

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

    /// Returns the number of tagged dofs.
    index_t taggedSize() const;

    /// Returns the number of eliminated dofs.
    inline index_t boundarySize() const
    {
      GISMO_ENSURE(m_curElimId>=0, "finalize() was not called on gsDofMapper");
        return m_numElimDofs.back();
    }

    index_t boundarySizeWithDuplicates() const;

    /// \brief Returns the position of local DoF 0 of patch \a k in the
    /// patch-concatenated local numbering of component \a c (a local
    /// position, not a global index; no shift).
    size_t offset(int k, index_t c = 0) const { return off(k,c); }

    /// Returns the number of patches present underneath the mapper
    size_t numPatches() const {return static_cast<size_t>(m_nPatches);}

    /// \brief Returns the total number of patch-local degrees of
    /// freedom that are being mapped
    size_t mapSize() const
    {
        size_t n = 0;
        for (size_t c = 0; c != m_dofs.size(); ++c)
            n += m_dofs[c].size();
        return n;
    }

    size_t componentsSize() const {return m_dofs.size();}

    /// \brief Returns the total number of patch-local DoFs
    /// that live on patch \a k for component \a c
    size_t patchSize(const index_t k, const index_t c = 0) const
    {
        GISMO_ASSERT(static_cast<size_t>(k)+1<=numPatches(), "Invalid patch index "<< k <<" >= "<< numPatches() );
        return off(k+1,c) - off(k,c);
    }

    size_t totalSize(const index_t c = 0) const
    {
        return m_dofs[c].size();
    }

    /// \brief For \a gl being a global index, this function returns a
    /// vector of pairs (patch,dof) that contains all the pairs which
    /// map to \a gl
    void preImage(index_t gl, std::vector<std::pair<index_t,index_t> > & result) const;

    /// \brief For \a gl being a global index, this function returns a
    /// pair (patch,dof) that maps to \a gl
    std::pair<index_t,index_t> anyPreImage(index_t gl) const;

    /// \brief Returns, for every global index of the mapper, a pair
    /// (patch, local dof) of component \a comp that maps to it.
    ///
    /// The vector has size() entries; the entry of global index \a gl
    /// is at position <tt>gl - firstIndex(0)</tt>. Entries without a
    /// pre-image in component \a comp are (-1,0).
    std::vector<std::pair<index_t,index_t> > anyPreImages(index_t comp = 0) const;

    /// \brief Produces the inverse of the mapping on patch \a k
    /// assuming that the map is invertible on that patch.
    ///
    /// The keys are (shifted) global indices, the values are the local
    /// dofs of patch \a k, collected over all components. If several
    /// local dofs of patch \a k map to the same global index, which of
    /// them is stored is unspecified.
    std::map<index_t,index_t> inverseOnPatch(const index_t k) const;

    /// \brief For \a gl being a global index, this function returns
    /// true whenever \a gl corresponds to patch \a k
    bool indexOnPatch(const index_t gl, const index_t k, index_t & local) const;

    inline bool indexOnPatch(const index_t gl, const index_t k) const
    {
        index_t local;
        return indexOnPatch(gl, k, local);
    }

    /// \brief For \a n being an index which is already offsetted, it
    /// returns the global index where it is mapped to by the dof
    /// mapper. The index \a n runs over the concatenation of the
    /// components' patch-concatenated local numberings.
    inline index_t mapIndex(index_t n) const
    {
        size_t c = 0;
        while (static_cast<size_t>(n) >= m_dofs[c].size())
            n -= static_cast<index_t>(m_dofs[c++].size());
        return m_dofs[c][n] + m_shift;
    }

    /// \brief Returns all boundary dofs of component \a c on patch k (local dof indices)
    gsVector<index_t> findBoundary(const index_t k, const index_t c = 0) const;

    /// \brief Returns all free dofs of component \a c on patch k (local dof indices)
    gsVector<index_t> findFree(const index_t k, const index_t c = 0) const;

    /// \brief Returns all coupled dofs of component \a c on patch k (local dof
    /// indices); with \a j >= 0 only those coupled with a dof of patch \a j
    gsVector<index_t> findCoupled(const index_t k, const index_t j = -1, const index_t c = 0) const;

    /// \brief Returns all free, not coupled dofs of component \a c on patch k (local dof indices)
    gsVector<index_t> findFreeUncoupled(const index_t k, const index_t c = 0) const;

    /// \brief Returns all tagged dofs of component \a c on patch k (local dof indices)
    gsVector<index_t> findTagged(const index_t k, const index_t c = 0) const;

private:

    template<class Predicate, class Iterator>
      static gsVector<index_t> find_impl(Iterator istart, Iterator iend, Predicate pred);

    void finalizeComp(const index_t comp);

    /// Position of local dof 0 of patch \a k in m_dofs[\a c]
    size_t off(index_t k, index_t c) const
    {
        GISMO_ASSERT(k >= 0 && k <= m_nPatches && c >= 0
                     && static_cast<size_t>(c+1)*(m_nPatches+1) <= m_offset.size(),
                     "Invalid patch "<< k <<" or component "<< c <<" in gsDofMapper");
        return m_offset[c*(m_nPatches+1)+k];
    }

    /// Number of valid local indices on patch \a k of component \a c: the patch
    /// size, or every dof of \a c in the layout of setIdentity(), where all dofs
    /// are on the last patch and are addressed through any patch
    size_t localBound(index_t k, index_t c) const
    { return 0 == off(m_nPatches-1, c) ? off(m_nPatches, c) : off(k+1, c) - off(k, c); }

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

    /// Per-component patch offsets into m_dofs.
    ///
    /// Has nComp*(m_nPatches+1) entries; entry c*(m_nPatches+1)+k (read via
    /// off(k,c)) is the position of local dof 0 of patch k in m_dofs[c].
    /// off(0,c)==0, off(.,c) is non-decreasing in k, and the sentinel
    /// off(m_nPatches,c)==m_dofs[c].size().
    /// A default-constructed mapper has no component, m_nPatches==1 and
    /// m_offset=={0,0}, so that offset(0) and patchSize(0) are 0.
    std::vector<size_t> m_offset;

    /// Number of patches
    index_t m_nPatches;

    /// Shifting of the global index (zero by default)
    index_t m_shift;

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

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsDofMapper
   */
  void pybind11_init_gsDofMapper(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace gismo
