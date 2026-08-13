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

        GISMO_ENSURE(static_cast<index_t>(m_labels.size()) == N,
                     "computeLabels() must produce exactly one label per element ("
                     << m_labels.size() << " labels, " << N << " elements).");

        // NOTE: the per-label range check that used to live here has moved into
        // setLabels() -- see the rationale there. It was incomplete at this
        // level: computeLabels() reaches m_labelsReady through setLabels()
        // before this line runs, so catching the throw and calling subdomains()
        // still performed the out-of-bounds write, and the external-label path
        // never passes through partition() at all.
        m_labelsReady = true;
    }

    /**
       @brief Per-element partition labels (length = numElements).

       Element with global id e belongs to partition labels()[e].
    */
    const std::vector<index_t>& labels() const
    {
        requireLabels();
        return m_labels;
    }

    /// @brief Number of partitions requested at construction.
    index_t nparts() const { return m_nparts; }

    /**
       @brief One gsIndexSubDomain per partition.

       The returned shared_ptrs can be passed directly to
       gsExprAssembler::setIntegrationDomain().
    */
    std::vector<typename gsDomain<T>::Ptr> subdomains() const
    {
        requireLabels();

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
    */
    std::vector<index_t> ownedElements(index_t rank, index_t nranks) const
    {
        requireLabels();
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
    : m_mb(mb), m_mapper(mapper), m_nparts(nparts), m_labelsReady(false)
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

}; // class gsPartitionerBase

} // namespace gismo
