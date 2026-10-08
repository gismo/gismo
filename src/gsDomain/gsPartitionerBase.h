/** @file gsPartitionerBase.h

    @brief Base class for element partitioners: everything that depends on
    nothing but the per-element partition labels.

    A partitioner assigns one partition label in [0, nparts) to every element
    of a gsMultiBasis' domain. How those labels are produced (graph
    partitioning, space-filling curve, ...) is the only thing a derived class
    has to supply: partition() below is a template method that validates the
    problem, calls the pure-virtual computeLabels(), and validates the result.
    Everything downstream of the labels -- the per-partition subdomains, the
    element ids owned by an MPI rank, the combined per-rank subdomain -- is
    implemented here, once.

    A derived class may instead run distributed: each MPI rank then keeps only
    the sorted ids of the elements of its own parts (setOwnElements()) and no
    full label vector. Only ownedElements()/subdomainForRank() for that rank
    and gatherLabels() are available then.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDofMapper.h>
#include <gsDomain/gsDomain.h>
#include <gsDomain/gsIndexSubDomain.h>
#include <gsDomain/gsPartitionedDofMapper.h>

#include <vector>

namespace gismo
{

/**
   @brief Common base for element partitioners.

   Derived classes implement computeLabels(), which must fill the label
   vector (via setLabels()) with exactly one partition label per element of
   multiBasis().domain(), and makeDofMapper(), which turns the labelling into
   a gsPartitionedDofMapper (the way the element -> DOF incidence is obtained
   is partitioner-specific: stored CSR, streamed, ...).

   In distributed mode (distributed() == true) the object holds only the
   element ids of the parts owned by one rank, ownedElements(rank, nranks) and
   subdomainForRank(rank, nranks) answer for that (rank, nranks) pair alone,
   and labels() and subdomains() fail with GISMO_ENSURE; gatherLabels()
   (collective, derived-class specific) recovers the full labels.

   \code
   gsMetisPartitioner<real_t> part(mb, mapper, 8);  // derives from this class
   part.partition();
   auto subdomains = part.subdomains();             // vector<shared_ptr<gsDomain<real_t>>>
   A.setIntegrationDomain(subdomains[k]);
   \endcode

   \ingroup Domain
*/
template<class T>
class gsPartitionerBase
{
public:

    virtual ~gsPartitionerBase() { }

    /**
       @brief Partition the domain (template method).

       Validates the problem size, delegates the actual labelling to the
       derived class' computeLabels(), and checks that it produced exactly one
       label per element. After this call labels(), subdomains(),
       ownedElements(), subdomainForRank() and makeDofMapper() are valid.
       In distributed mode there is no full label vector to check; only
       ownedElements()/subdomainForRank() for the own (rank, nranks) and
       gatherLabels() are valid afterwards.

       The checks are GISMO_ENSURE (not GISMO_ASSERT) on purpose: a bad
       nparts must be caught in Release builds too, since that is the only
       build that runs at scale.
    */
    void partition()
    {
        const index_t N = static_cast<index_t>(m_mb.domain()->numElements());
        GISMO_ENSURE(N > 0, "Empty domain - nothing to partition.");
        GISMO_ENSURE(m_nparts > 0 && m_nparts <= N,
                     "nparts must be in [1, numElements] (nparts=" << m_nparts
                     << ", numElements=" << N << ").");

        computeLabels();

        if (!m_distributed)
            GISMO_ENSURE(static_cast<index_t>(m_labels.size()) == N,
                         "computeLabels() must produce exactly one label per element ("
                         << m_labels.size() << " labels, " << N << " elements).");

        // NOTE: no per-label range check here on purpose. It lives in
        // setLabels(), the one choke point every labelling and external-label
        // path passes through (rationale there); a check at this level would
        // run after m_labelsReady is set and would never see the external-label
        // path. A distributed result is validated in setOwnElements().
        m_labelsReady = true;
    }

    /**
       @brief Per-element partition labels (length = numElements).

       Element with global id e belongs to partition labels()[e]. Not
       available on a distributed partitioner (GISMO_ENSURE): use the
       collective gatherLabels() there.
    */
    const std::vector<index_t>& labels() const
    {
        requireLabels();
        GISMO_ENSURE(!m_distributed, "gsPartitionerBase::labels(): this partitioner is distributed (parallel constructor on more than one rank) and keeps only the calling rank's own elements; use gatherLabels() (collective) for the full labels, or the serial constructor.");
        return m_labels;
    }

    /**
       @brief The full per-element labels (length numElements,
       labels[e] = part of element e) on every rank.

       Collective in a distributed partitioner: every rank of its communicator
       must call it. The base default returns a copy of labels(), so
       partition() must have been called first; O(N) time and memory.
       Intended for diagnostics, since it materialises all N labels.
    */
    virtual std::vector<index_t> gatherLabels() const { return labels(); }

    /// @brief Number of partitions requested at construction.
    index_t nparts() const { return m_nparts; }

    /// @brief True when the object keeps only the calling rank's own elements
    /// and no full label vector. Serial partitioners are never distributed.
    bool distributed() const { return m_distributed; }

    /**
       @brief One gsIndexSubDomain per partition.

       The returned shared_ptrs can be passed directly to
       gsExprAssembler::setIntegrationDomain().

       Not available on a distributed partitioner (GISMO_ENSURE); use
       subdomainForRank() there. Complexity: O(N).
    */
    std::vector<typename gsDomain<T>::Ptr> subdomains() const
    {
        requireLabels();
        GISMO_ENSURE(!m_distributed, "gsPartitionerBase::subdomains(): this partitioner is distributed (parallel constructor on more than one rank) and has no full labelling; use subdomainForRank(comm.rank(), comm.size()), gatherLabels() (collective), or the serial constructor.");

        // Collect element ids per partition
        std::vector<std::vector<index_t>> partElems(m_nparts);
        for (index_t e = 0; e < static_cast<index_t>(m_labels.size()); ++e)
            partElems[m_labels[e]].push_back(e);

        // Build subdomain objects.  A single shared domain instance is co-owned
        // by all subdomains so the parent domain outlives every subdomain.
        typename gsDomain<T>::Ptr dom = m_mb.domain();
        std::vector<typename gsDomain<T>::Ptr> result;
        result.reserve(m_nparts);
        for (index_t k = 0; k < m_nparts; ++k)
            result.push_back(
                memory::make_shared( new gsIndexSubDomain<T>(dom, give(partElems[k])) ));
        return result;
    }

    /**
       @brief Sorted global element ids of all partitions owned by \a rank,
       under the cyclic assignment gsPartitionedDofMapper::rankOfPart(part,
       nranks) -- the same convention makeDofMapper() (and hence DOF
       ownership) uses, kept in exactly one place
       (gsPartitionedDofMapper::rankOfPart) rather than being re-derived here
       and in DOF ownership separately.

       Complexity: O(N) for N elements, one rankOfPart() modulo per element;
       O(n_own) copy of the stored list on a distributed partitioner, where
       only the (rank, nranks) pair the object was distributed over is
       available (GISMO_ENSURE otherwise). The result is strictly increasing.
    */
    std::vector<index_t> ownedElements(index_t rank, index_t nranks) const
    {
        requireLabels();
        if (m_distributed)
        {
            GISMO_ENSURE(rank == m_ownRank && nranks == m_ownNranks,
                         "gsPartitionerBase::ownedElements(" << rank << ", " << nranks
                         << "): only (rank, nranks) = (" << m_ownRank << ", " << m_ownNranks
                         << ") is available on a distributed partitioner; use gatherLabels() or the serial constructor.");
            return m_ownElements;
        }
        std::vector<index_t> result;
        for (index_t e = 0; e < static_cast<index_t>(m_labels.size()); ++e)
            if (gsPartitionedDofMapper::rankOfPart(m_labels[e], nranks) == rank)
                result.push_back(e); // e ascending -> result already sorted
        return result;
    }

    /**
       @brief One combined gsIndexSubDomain for \a rank -- the union of all
       partitions gsPartitionedDofMapper::rankOfPart assigns to it. Can be
       passed directly to gsExprAssembler::setIntegrationDomain(), replacing
       the hand-rolled ownedElements loop + gsSubDomain downcast a caller
       would otherwise need (see gsMetisPetscAssembly_example.cpp).

       Complexity: O(N) for the ownedElements() scan (O(n_own) on a
       distributed partitioner, for its own (rank, nranks) only), plus the
       gsIndexSubDomain construction cost. The index list is strictly
       increasing, so the constructor skips its sort. Allocates a fresh
       composite domain through gsMultiBasis::domain(), O(nPatches).
    */
    typename gsDomain<T>::Ptr subdomainForRank(index_t rank, index_t nranks) const
    {
        typename gsDomain<T>::Ptr dom = m_mb.domain();
        return memory::make_shared(
            new gsIndexSubDomain<T>(dom, ownedElements(rank, nranks)));
    }

    /**
       @brief Build the DOF-ownership map and global permutation for
       distributing this partitioning across \a nranks MPI ranks.

       Implemented by the derived class, since the element -> free-DOF
       incidence it needs is partitioner-specific (stored per-element CSR,
       streamed, ...).
    */
    virtual gsPartitionedDofMapper makeDofMapper(index_t nranks) const = 0;

protected:

    /**
       @brief Construct (does not partition yet -- call partition()).

       @param mb      Multi-patch basis.
       @param mapper  Finalized DOF mapper.
       @param nparts  Number of partitions.
    */
    gsPartitionerBase(const gsMultiBasis<T>& mb,
                      const gsDofMapper&     mapper,
                      index_t                nparts)
    : m_mb(mb), m_mapper(mapper), m_nparts(nparts), m_labelsReady(false),
      m_distributed(false), m_ownRank(-1), m_ownNranks(0)
    { }

    /// @brief Fill m_labels with one partition label per element. Called by
    /// partition(); implementations may also be fed externally-supplied
    /// labels through setLabels().
    virtual void computeLabels() = 0;

    /// @brief Set the labels (and mark them ready). Used by computeLabels()
    /// implementations and by derived setters accepting external labels.
    ///
    /// The range check lives HERE, not in partition(), because this is the one
    /// choke point through which labels ever become usable: every
    /// computeLabels() implementation routes through it, and so does every
    /// external-label setter (gsMetisPartitioner::setPartLabels). A check in
    /// partition() alone is bypassed twice over -- once by the external-label
    /// path, which never calls partition(), and once by a caller that catches
    /// partition()'s throw and then calls subdomains() anyway, since
    /// m_labelsReady was already true by then. Either route reaches the
    /// unguarded partElems[m_labels[e]] write in subdomains().
    ///
    /// Validate BEFORE mutating: on a bad label the object is left untouched
    /// and still usable, rather than half-updated with m_labelsReady set.
    void setLabels(std::vector<index_t> labels)
    {
        for (size_t e = 0; e != labels.size(); ++e)
            GISMO_ENSURE(labels[e] >= 0 && labels[e] < m_nparts,
                         "Partition label out of range for element " << e
                         << ": " << labels[e] << " not in [0," << m_nparts
                         << ").");
        m_labels      = give(labels);
        m_labelsReady = true;
        std::vector<index_t>().swap(m_ownElements);
        m_ownRank     = -1;
        m_ownNranks   = 0;
        m_distributed = false;
    }

    /// @brief Switch to distributed mode: keep only \a ids, the sorted global
    /// ids of the elements of the parts owned by \a rank out of \a nranks
    /// (gsPartitionedDofMapper::rankOfPart(part, nranks) == rank), and drop the
    /// full label vector.
    ///
    /// Layout: \a ids strictly increasing (hence unique), each in
    /// [0, numElements). The check is a local GISMO_ENSURE on a post-condition
    /// of the collective algorithm that produced the ids, not a data-dependent
    /// failure to be agreed: such an algorithm hands identical decisions to all
    /// ranks, so a failure here is a bug. Validates BEFORE mutating, like
    /// setLabels().
    ///
    /// Complexity: O(n_loc + nPatches) (the latter for the element count of
    /// the composite domain).
    void setOwnElements(std::vector<index_t> ids, index_t rank, index_t nranks)
    {
        const index_t N = static_cast<index_t>(m_mb.domain()->numElements());
        GISMO_ENSURE(nranks >= 1 && rank >= 0 && rank < nranks,
                     "setOwnElements: invalid (rank, nranks) = (" << rank << ", " << nranks << ").");
        for (size_t i = 0; i != ids.size(); ++i)
        {
            GISMO_ENSURE(ids[i] >= 0 && ids[i] < N,
                         "setOwnElements: element id " << ids[i] << " not in [0," << N << ").");
            GISMO_ENSURE(i == 0 || ids[i - 1] < ids[i],
                         "setOwnElements: element ids must be strictly increasing (position " << i << ").");
        }
        m_ownElements = give(ids);
        m_ownRank     = rank;
        m_ownNranks   = nranks;
        std::vector<index_t>().swap(m_labels);
        m_distributed = true;
        m_labelsReady = true;
    }

    /// @brief Precondition check shared by every label-dependent accessor.
    void requireLabels() const
    { GISMO_ASSERT(m_labelsReady, "Call partition() first."); }

    const gsMultiBasis<T>& multiBasis() const { return m_mb; }
    const gsDofMapper&     mapper()     const { return m_mapper; }

    const gsMultiBasis<T>& m_mb;
    const gsDofMapper&     m_mapper;
    const index_t          m_nparts;
    std::vector<index_t>   m_labels;
    bool                   m_labelsReady;
    bool                   m_distributed; ///< true when only the calling rank's own elements are kept, with no full label vector
    std::vector<index_t>   m_ownElements; ///< distributed mode: sorted, unique global ids of the elements of the parts owned by m_ownRank (rankOfPart(part, m_ownNranks) == m_ownRank); empty otherwise
    index_t                m_ownRank;     ///< distributed mode: the rank m_ownElements belongs to; -1 otherwise
    index_t                m_ownNranks;   ///< distributed mode: the number of ranks of that ownership; 0 otherwise

}; // class gsPartitionerBase

} // namespace gismo
