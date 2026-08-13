/** @file gsPartitionedDofMapper.h

    @brief DOF ownership and global permutation for a METIS-partitioned mesh.

    Every rank replicates the full METIS partition labelling (broadcast once,
    see gsMetisPartitioner::setPartLabels()), so ownership and the global
    permutation below are computable identically, without communication, on
    every rank. Given the partition-label assignment and a cyclic
    part -> rank map, this class answers, for every free global DOF:

      - which partition "owns" it (the lowest-numbered partition with an
        element active on it),
      - which rank that resolves to,
      - and its row in a PETSc layout where rows are grouped contiguously by
        owning rank (the permutation gismo's global-indexed matrix output
        must go through on insertion; see gsPetscLocalToGlobal.h).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsLinearAlgebra.h>

#include <limits>
#include <vector>

namespace gismo
{

/**
   @brief DOF ownership + global permutation derived from a METIS partition.

   Construction takes the element -> free-DOF incidence (as a
   vector-of-vectors, e.g. built locally from
   gsElementGraph::elementFreeDofsCSR()), the per-element partition labels, the
   total number of free DOFs, and the number of MPI ranks. It does not
   depend on gsElementGraph itself, nor on METIS's idx_t — everything here
   is plain index_t, computed once and replicated identically on every
   rank (see gsMetisPartitioner::makeDofMapper()).

   Ownership rule: \c ownerPart(g) is the smallest partition label among the
   partitions that have at least one element active on free DOF \a g.
   \c ownerRank(g) = \c rankOfPart(ownerPart(g), nranks), with the cyclic
   assignment \c rankOfPart(p,n) = p % n.

   Permutation: rows are grouped contiguously by owning rank, in ascending
   global-DOF order within a rank, i.e.
   \code
   perm(g) = rankOffset(ownerRank(g)) + (rank of g among the DOFs owned by ownerRank(g))
   \endcode
   \c Σ_r numOwnedDofs(r) == freeSize() always holds — every free DOF has
   exactly one owning rank — which is what makes \c permutation() valid
   input to \c MatSetSizes / a PETSc layout.

   \ingroup Domain
*/
class gsPartitionedDofMapper
{
public:

    gsPartitionedDofMapper() : m_nDofs(0), m_nranks(0) { }

    /**
       @param elemFreeDofs  Per-element list of active free global DOFs
                             (generic vector-of-vectors; nothing in this
                             library builds one from graph data anymore --
                             prefer the CSR-span constructor below, fed from
                             gsElementGraph::elementFreeDofsCSR()).
       @param partLabels    Per-element partition label (converted to
                             index_t at the METIS boundary — see
                             gsMetisPartitioner::makeDofMapper()).
       @param nDofs         Total number of free DOFs (gsDofMapper::freeSize()).
       @param nranks        Number of MPI ranks partitions are cyclically
                             assigned to.
    */
    gsPartitionedDofMapper(const std::vector<std::vector<index_t>>& elemFreeDofs,
                            const std::vector<index_t>&              partLabels,
                            index_t                                  nDofs,
                            index_t                                  nranks)
    {
        build(elemFreeDofs, partLabels, nDofs, nranks);
    }

    /**
       @brief CSR-span overload: same construction as above, but from a flat
       CSR pair (gsElementGraph::elementFreeDofsCSR()) instead of a
       vector-of-vectors, so the caller never has to materialize one
       std::vector<index_t> per element just to hand it to this constructor.

       @param dofsOffsets  CSR offsets into \a dofsValues; size
                             partLabels.size()+1 (one entry per element, plus
                             the end sentinel).
       @param dofsValues   Flat per-element free global DOF ids: element e's
                             DOFs are dofsValues[dofsOffsets[e] .. dofsOffsets[e+1]).
       @param partLabels    Per-element partition label.
       @param nDofs         Total number of free DOFs (gsDofMapper::freeSize()).
       @param nranks        Number of MPI ranks partitions are cyclically
                             assigned to.
    */
    gsPartitionedDofMapper(const std::vector<index_t>& dofsOffsets,
                            const std::vector<index_t>& dofsValues,
                            const std::vector<index_t>& partLabels,
                            index_t                      nDofs,
                            index_t                      nranks)
    {
        buildCSR(dofsOffsets, dofsValues, partLabels, nDofs, nranks);
    }

    /**
       @brief Build directly from per-DOF ownership ranges, without any
       element -> DOF incidence.

       For callers that stream the element mesh and fold ownership on the fly
       (see gsGeometricPartitioner) instead of storing a per-element DOF CSR.
       \a minPart[g] / \a maxPart[g] must be the min / max partition label over
       all elements active on free DOF \a g -- exactly the reduction the two
       constructors above perform internally, so the result is identical to
       constructing from the incidence those labels came from.

       \b Index convention: \a minPart and \a maxPart are indexed by
       <em>compact free index</em> in [0, nDofs), i.e. a gsDofMapper global
       index already shifted by gsDofMapper::firstIndex(). They are NOT indexed
       by raw global index.

       @param minPart  size nDofs; smallest partition label active on DOF g.
       @param maxPart  size nDofs; largest  partition label active on DOF g,
                       or -1 if no element touches g (rejected, see below).
       @param nDofs    Total number of free DOFs (gsDofMapper::freeSize()).
       @param nranks   Number of MPI ranks (cyclic part -> rank assignment).
    */
    static gsPartitionedDofMapper fromOwnershipRange(const std::vector<index_t>& minPart,
                                                     const std::vector<index_t>& maxPart,
                                                     index_t                     nDofs,
                                                     index_t                     nranks)
    {
        GISMO_ENSURE(nDofs >= 0, "nDofs must be non-negative.");
        GISMO_ENSURE(nranks > 0, "nranks must be positive.");
        GISMO_ENSURE(static_cast<index_t>(minPart.size()) == nDofs &&
                     static_cast<index_t>(maxPart.size()) == nDofs,
                     "minPart/maxPart must have exactly nDofs entries "
                     "(compact free indices in [0,nDofs)).");

        // GISMO_ENSURE, not GISMO_ASSERT: a DOF that no element touches still
        // yields a bijective permutation with consistent owned counts, so no
        // other invariant check can catch it -- and GISMO_ASSERT vanishes under
        // NDEBUG, i.e. in exactly the release builds that run at scale.
        //
        // The minPart[g] >= 0 half is memory safety, not just hygiene:
        // finishBuild computes m_ownerRank[g] = minPart[g] % nranks, which is
        // NEGATIVE for a negative label, and then does ++m_numOwned[ownerRank]
        // -- an out-of-bounds write that corrupts the heap (observed:
        // "double free or corruption (out)", exit 134, in BOTH debug and
        // release). minPart[g] <= maxPart[g] is the matching ordering
        // invariant; the two together subsume the old maxPart[g] >= 0 check.
        for (index_t g = 0; g < nDofs; ++g)
        {
            GISMO_ENSURE(maxPart[g] >= 0,
                         "Free DOF "<<g<<" is not active on any element -- "
                         "cannot determine ownership. (Did the producer skip "
                         "mapper components 1..numComponents()-1?)");
            GISMO_ENSURE(minPart[g] >= 0 && minPart[g] <= maxPart[g],
                         "Free DOF "<<g<<" has an invalid ownership range ["
                         <<minPart[g]<<","<<maxPart[g]<<"]: minPart must be "
                         "non-negative and not exceed maxPart.");
        }

        gsPartitionedDofMapper result;
        result.finishBuild(nDofs, nranks, minPart, maxPart);
        return result;
    }

    /// Cyclic part -> rank assignment; the single source of truth used both
    /// here and by callers that need to know which rank owns which partitions.
    static index_t rankOfPart(index_t part, index_t nranks)
    {
        GISMO_ASSERT(nranks > 0, "nranks must be positive.");
        return part % nranks;
    }

    /// Smallest partition label with an element active on free DOF \a g.
    index_t ownerPart(index_t g) const
    {
        GISMO_ASSERT(g >= 0 && g < m_nDofs, "DOF index out of range.");
        return m_ownerPart[g];
    }

    /// Rank owning free DOF \a g (== rankOfPart(ownerPart(g), nranks())).
    index_t ownerRank(index_t g) const
    {
        GISMO_ASSERT(g >= 0 && g < m_nDofs, "DOF index out of range.");
        return m_ownerRank[g];
    }

    /// True if \a g is touched by elements of more than one partition
    /// (a diagnostic, generalizing the old per-partition interior/interface
    /// classification to "does this DOF need cross-partition agreement").
    bool isInterfaceDof(index_t g) const
    {
        GISMO_ASSERT(g >= 0 && g < m_nDofs, "DOF index out of range.");
        return m_isInterface[g];
    }

    /// Number of free DOFs owned by \a rank. Σ_r numOwnedDofs(r) == freeSize().
    index_t numOwnedDofs(index_t rank) const
    {
        GISMO_ASSERT(rank >= 0 && rank < m_nranks, "Rank index out of range.");
        return m_numOwned[rank];
    }

    /// First permuted row index owned by \a rank (prefix sum of numOwnedDofs).
    index_t rankOffset(index_t rank) const
    {
        GISMO_ASSERT(rank >= 0 && rank < m_nranks, "Rank index out of range.");
        return m_rankOffset[rank];
    }

    /// Global permutation: perm(g) is the row free DOF \a g occupies in a
    /// PETSc layout where rows are grouped contiguously by owning rank.
    const gsVector<index_t>& permutation() const { return m_perm; }

    index_t freeSize() const { return m_nDofs; }
    index_t nranks()   const { return m_nranks; }

    /// \brief Returns the number of bytes actually allocated (capacity, not
    /// logical size) by this mapper's containers.
    size_t nBytes() const
    {
        size_t bytes = sizeof(*this);

        bytes += m_ownerPart.capacity()  * sizeof(index_t);
        bytes += m_ownerRank.capacity()  * sizeof(index_t);
        // std::vector<bool> is bit-packed: capacity() counts bits, not bytes;
        // dividing by 8 under-counts by at most one word.
        bytes += m_isInterface.capacity() / 8;
        bytes += m_numOwned.capacity()   * sizeof(index_t);
        bytes += m_rankOffset.capacity() * sizeof(index_t);
        // gsVector (Eigen) has no capacity() concept for a plain vector, so
        // size() is exact here, not a lower bound.
        bytes += static_cast<size_t>(m_perm.size()) * sizeof(index_t);

        return bytes;
    }

private:

    void build(const std::vector<std::vector<index_t>>& elemFreeDofs,
               const std::vector<index_t>&               partLabels,
               index_t                                   nDofs,
               index_t                                   nranks)
    {
        GISMO_ASSERT(elemFreeDofs.size() == partLabels.size(),
                     "elemFreeDofs and partLabels must have the same length "
                     "(one entry per element).");
        GISMO_ASSERT(nDofs >= 0, "nDofs must be non-negative.");
        GISMO_ASSERT(nranks > 0, "nranks must be positive.");

        const index_t minInit = std::numeric_limits<index_t>::max();
        std::vector<index_t> minPart(nDofs, minInit);
        std::vector<index_t> maxPart(nDofs, -1);

        const index_t nElem = static_cast<index_t>(elemFreeDofs.size());
        for (index_t e = 0; e < nElem; ++e)
        {
            const index_t p = partLabels[e];
            for (index_t g : elemFreeDofs[e])
            {
                GISMO_ASSERT(g >= 0 && g < nDofs, "elemFreeDofs entry out of range.");
                if (p < minPart[g]) minPart[g] = p;
                if (p > maxPart[g]) maxPart[g] = p;
            }
        }

        finishBuild(nDofs, nranks, minPart, maxPart);
    }

    // Same ownership computation as build(), but reads element->DOF
    // incidence from a flat CSR pair instead of a vector-of-vectors --
    // avoids materializing one std::vector<index_t> per element just to
    // feed it to this class (see gsElementGraph::elementFreeDofsCSR()).
    void buildCSR(const std::vector<index_t>& dofsOffsets,
                  const std::vector<index_t>& dofsValues,
                  const std::vector<index_t>& partLabels,
                  index_t                     nDofs,
                  index_t                     nranks)
    {
        GISMO_ASSERT(dofsOffsets.size() == partLabels.size() + 1,
                     "dofsOffsets must have one entry per element plus the "
                     "end sentinel (size partLabels.size()+1).");
        GISMO_ASSERT(dofsOffsets.empty() ||
                     static_cast<size_t>(dofsOffsets.back()) == dofsValues.size(),
                     "dofsOffsets.back() must equal dofsValues.size().");
        GISMO_ASSERT(nDofs >= 0, "nDofs must be non-negative.");
        GISMO_ASSERT(nranks > 0, "nranks must be positive.");

        const index_t minInit = std::numeric_limits<index_t>::max();
        std::vector<index_t> minPart(nDofs, minInit);
        std::vector<index_t> maxPart(nDofs, -1);

        const index_t nElem = static_cast<index_t>(partLabels.size());
        for (index_t e = 0; e < nElem; ++e)
        {
            const index_t p = partLabels[e];
            for (index_t k = dofsOffsets[e]; k < dofsOffsets[e + 1]; ++k)
            {
                const index_t g = dofsValues[k];
                GISMO_ASSERT(g >= 0 && g < nDofs, "dofsValues entry out of range.");
                if (p < minPart[g]) minPart[g] = p;
                if (p > maxPart[g]) maxPart[g] = p;
            }
        }

        finishBuild(nDofs, nranks, minPart, maxPart);
    }

    // Ownership/permutation computation shared by build() and buildCSR():
    // both reduce their input to the same (minPart, maxPart) per DOF before
    // reaching here.
    void finishBuild(index_t nDofs, index_t nranks,
                      const std::vector<index_t>& minPart,
                      const std::vector<index_t>& maxPart)
    {
        m_nDofs  = nDofs;
        m_nranks = nranks;

        m_ownerPart.resize(nDofs);
        m_ownerRank.resize(nDofs);
        m_isInterface.resize(nDofs);
        for (index_t g = 0; g < nDofs; ++g)
        {
            GISMO_ASSERT(maxPart[g] >= 0,
                         "Free DOF "<<g<<" is not active on any element — "
                         "cannot determine ownership.");
            m_ownerPart[g]   = minPart[g];
            m_ownerRank[g]   = rankOfPart(minPart[g], nranks);
            m_isInterface[g] = (minPart[g] != maxPart[g]);
        }

        // Owned counts and prefix offsets: Σ_r numOwnedDofs(r) == nDofs by
        // construction (every DOF has exactly one ownerRank).
        m_numOwned.assign(nranks, 0);
        for (index_t g = 0; g < nDofs; ++g)
            ++m_numOwned[m_ownerRank[g]];

        m_rankOffset.assign(nranks, 0);
        for (index_t r = 1; r < nranks; ++r)
            m_rankOffset[r] = m_rankOffset[r-1] + m_numOwned[r-1];

        // perm(g) = rankOffset(ownerRank(g)) + ownedLocalIndex(g), with
        // ownedLocalIndex assigned in ascending g order per rank — a pure
        // function of (elemFreeDofs, partLabels, nDofs, nranks), so it is
        // bit-identical across independent constructions (i.e. across ranks).
        std::vector<index_t> cursor(m_rankOffset);
        m_perm.resize(nDofs);
        for (index_t g = 0; g < nDofs; ++g)
        {
            const index_t r = m_ownerRank[g];
            m_perm[g] = cursor[r]++;
        }
    }

    index_t              m_nDofs;
    index_t              m_nranks;
    std::vector<index_t> m_ownerPart;
    std::vector<index_t> m_ownerRank;
    std::vector<bool>    m_isInterface;
    std::vector<index_t> m_numOwned;
    std::vector<index_t> m_rankOffset;
    gsVector<index_t>    m_perm;

}; // class gsPartitionedDofMapper

} // namespace gismo
