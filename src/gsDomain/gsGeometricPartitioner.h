/** @file gsGeometricPartitioner.h

    @brief Graph-free geometric partitioner for multi-patch IGA element meshes.

    Partitions the elements of a gsMultiBasis' domain into nparts balanced
    parts without ever building an element adjacency graph. Two streaming
    passes over the domain are enough:

      - Pass A computes one physical centroid and one integer weight per
        element (the geometry evaluation is batched in chunks, since a
        per-element evaluation would dominate the whole pass).
      - The label strategy (recursive coordinate bisection, or a
        Morton/Hilbert space-filling curve, see gsSpaceFillingCurve) turns
        the centroids + weights into one partition label per element.
      - Pass B streams the mesh a second time and folds per-DOF ownership
        directly into minPart[]/maxPart[], which
        gsPartitionedDofMapper::fromOwnershipRange() consumes. No per-element
        DOF list (CSR) is ever stored -- that storage is precisely what this
        class exists to avoid.

    This is a drop-in peer of gsMetisPartitioner: both derive from
    gsPartitionerBase<T>, so subdomains(), ownedElements() and
    subdomainForRank() are shared, and a caller can switch between them. A
    partitioner built with the communicator constructor on more than one rank
    is distributed: it exposes only the calling rank's own part (see
    gsPartitionerBase), and gatherLabels() for diagnostics.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsDofMapper.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsMultiPatch.h>
#include <gsDomain/gsDomain.h>
#include <gsDomain/gsDistributedPartition.h>
#include <gsDomain/gsIndexSubDomain.h>
#include <gsDomain/gsPartitionedDofMapper.h>
#include <gsDomain/gsPartitionerBase.h>
#include <gsParallel/gsMpi.h>
#include <gsUtils/gsSpaceFillingCurve.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace gismo
{

/**
   @brief Graph-free geometric partitioner for IGA element meshes.

   Assigns one partition label in [0, nparts) to every element of
   mb.domain(), using only element centroids and, by default, unit element
   weights; per-element free-DOF counts are available as an opt-in weighting
   via Options::weightByDofs. Three label strategies are available:

     - \c rcb     : recursive coordinate bisection (default),
     - \c hilbert : sort along a Hilbert curve, cut by weighted prefix sum,
     - \c morton  : the same, along a Morton (Z-order) curve.

   \code
   gsGeometricPartitioner<real_t> part(mp, mb, mapper, 8);
   part.partition();
   auto subdomains = part.subdomains();       // from gsPartitionerBase
   A.setIntegrationDomain(subdomains[k]);
   auto pdm = part.makeDofMapper(nranks);     // graph-free DOF ownership
   \endcode

   \note Every ordering used internally is a strict \em total order: the
   primary key is compared first and every tie is broken on the element id,
   so the labels do not depend on the sort algorithm or on the rank count.
   With the gsMpiComm constructor on P > 1 ranks the partitioner is
   distributed: each rank evaluates the centroids and weights of its own
   element slice only, the labelling runs as a collective state machine
   (gsDistributedRcb, gsDistributedCurve) and each rank ends up with the sorted
   ids of its own elements only. Use gatherLabels() when the full labels are
   needed (diagnostics).

   \ingroup Domain
*/
template<class T>
class gsGeometricPartitioner : public gsPartitionerBase<T>
{
    typedef gsPartitionerBase<T> Base;

public:

    /// Label strategy. Default is rcb.
    enum Strategy { rcb = 0, hilbert = 1, morton = 2 };

    struct Options
    {
        Strategy strategy      = rcb;   ///< label strategy
        bool     weightByDofs  = false; ///< opt in to weighting elements by their free-DOF count
        index_t  chunkSize     = 4096;  ///< elements per batched geometry evaluation
        bool     keepCentroids = false; ///< serial path only: retain centroids() after partition() (geoDim*8 B/element)
    };

    /**
       @brief Construct (does not partition yet -- call partition()).

       @param mp      Multi-patch geometry (supplies the physical centroids).
       @param mb      Multi-patch basis; its domain() is the element mesh.
       @param mapper  Finalized DOF mapper.
       @param nparts  Number of partitions.
       @param opts    Optional tuning parameters.

       \note mb.domain() builds a fresh gsCompositeDomain on every call, so
       the domain is obtained exactly once, here, and reused by the
       tensor-patch guard, pass A and pass B. Those three must agree on the
       element numbering.
    */
    gsGeometricPartitioner(const gsMultiPatch<T>& mp,
                           const gsMultiBasis<T>& mb,
                           const gsDofMapper&     mapper,
                           index_t                nparts,
                           Options                opts = Options{})
    : Base(mb, mapper, nparts), m_mp(mp), m_opts(opts), m_dom(mb.domain()),
      m_numElements(static_cast<index_t>(m_dom->numElements()))
    { }

    /**
       @brief Construct a partitioner whose pass A and labelling are
       distributed over the ranks of \a comm (does not partition yet -- call
       partition()).

       @param mp      Multi-patch geometry (supplies the physical centroids).
       @param mb      Multi-patch basis; its domain() is the element mesh.
       @param mapper  Finalized DOF mapper.
       @param nparts  Number of partitions; must equal comm.size() when
       comm.size() > 1, any value otherwise.
       @param comm    Communicator over which pass A and the labelling are
       split (comm.size() > 1: each rank keeps only its own element ids).
       @param opts    Optional tuning parameters.

       Parallel part: with P = comm.size() > 1 the object is always
       distributed. Rank r evaluates the geometry (and, with
       Options::weightByDofs, counts free DOFs) only on the contiguous element
       slice [N r/P, N (r+1)/P) of the N elements, which is O(N/P) element
       visits and evaluations per rank plus one O(dim + log nPatches) iterator
       jump to the slice start, and stores only that slice.

       Communication: one scalar max-reduction in pass A (it makes the
       iterator-id, slice-coverage and finite-centroid checks collective),
       then the requests of the label machine: the RCB machine issues at most
       4 + ceil(log2 nparts) (3 + 17 + maxPasses) reductions, the curve
       machine at most 1 + maxPasses + 1 + 3, with the request sizes and
       bounds stated in gsDistributedPartition.h (gsDistributedRcb,
       gsDistributedCurve). The last step is an all-to-all of about N/P
       element ids per rank.

       Memory per rank: O(N/P) -- the slice centroids (geoDim*sizeof(T) B per
       element), weights and curve keys -- plus O(P*256) select rows and the
       own ids; no allocation is sized by N. Replicated on every rank are the
       control-net bounding box and the multi-patch geometry. Complexity per
       rank: O(n_loc * passes) per bisection level for RCB, O(n_loc log P) per
       select pass for the curves, n_loc = N/P.

       Available after partition() on a distributed object: only
       subdomainForRank(comm.rank(), comm.size()), ownedElements() for the same
       pair, gatherLabels() (collective, O(N) memory), numElements(), nparts()
       and distributed(). labels(), subdomains(), elementWeights(),
       centroids() and makeDofMapper() fail with GISMO_ENSURE: dof numbering
       with this constructor is rendezvous-only.

       Distributed RCB needs IEEE float or double coordinates (its keys are
       order-preserving bit images); the curve strategies work for any T.

       Collective: every rank of \a comm must construct with identical
       arguments (mp, mb, mapper, nparts, opts) and call partition(). The
       communicator is copied (a gsMpiComm is a non-owning handle); the
       underlying MPI communicator must stay alive until partition() returns
       and during every later gatherLabels() call.
       A comm with size() <= 1 behaves exactly like the serial constructor and
       performs no communication.

       For the same slice centroids and weights, the RCB part sets equal those
       of the replicated recursion rcbSplit(), and the curve labels equal
       those of curveLabels(). They are not guaranteed to be bit-identical to
       those of the serial constructor, since the evaluation chunk boundaries
       differ and a batched eval_into() could in principle round differently.
    */
    gsGeometricPartitioner(const gsMultiPatch<T>& mp,
                           const gsMultiBasis<T>& mb,
                           const gsDofMapper&     mapper,
                           index_t                nparts,
                           const gsMpiComm&       comm,
                           Options                opts = Options{})
    : Base(mb, mapper, nparts), m_mp(mp), m_opts(opts), m_dom(mb.domain()),
      m_numElements(static_cast<index_t>(m_dom->numElements())), m_comm(comm)
    {
        GISMO_ENSURE(comm.size() >= 1 && comm.rank() >= 0 && comm.rank() < comm.size(),
                     "gsGeometricPartitioner: invalid communicator (size "
                     << comm.size() << ", rank " << comm.rank() << ").");
        if (comm.size() > 1)
            GISMO_ENSURE(nparts == comm.size(), "gsGeometricPartitioner: with a communicator of size "
                         << comm.size() << " nparts must equal comm.size() (got " << nparts << ").");
    }

    /// @brief "rcb" | "hilbert" | "morton"; fails on anything else.
    static Strategy strategyFromString(const std::string& s)
    {
        if ("rcb"     == s) return rcb;
        if ("hilbert" == s) return hilbert;
        if ("morton"  == s) return morton;
        GISMO_ERROR("gsGeometricPartitioner: unknown strategy \""<<s<<"\" "
                    "(expected \"rcb\", \"hilbert\" or \"morton\").");
    }

    /**
       @brief The full per-element labels on every rank. Collective on a
       distributed partitioner (every rank of the communicator must call it);
       the serial/size-1 case returns a copy of labels().

       Layout: the result has length N, labels[e] = part of element e. Each
       rank contributes its sorted own ids; part q is owned by rank q
       (nparts == comm.size()), so the gathered ids of source rank q, at
       displ[q] in the receive buffer, all get label q.

       Communication: one allgather of P counts and one allgatherv of N ids.
       Memory O(N): this is the only place of the P > 1 path with N-sized
       buffers, hence diagnostic use only. Requires N <= INT_MAX.
    */
    std::vector<index_t> gatherLabels() const override
    {
        if (!this->distributed()) return Base::gatherLabels();
        this->requireLabels();
        const index_t N = m_numElements;
        const int     P = m_comm.size();
        // N is replicated, so this fires on all ranks or on none
        GISMO_ENSURE(static_cast<int64_t>(N) <= static_cast<int64_t>(std::numeric_limits<int>::max()),
                     "gsGeometricPartitioner::gatherLabels(): the element count " << N
                     << " exceeds the MPI int count range.");
        int cnt = static_cast<int>(this->m_ownElements.size());
        std::vector<int> counts(P), displ(P);
        int rc = m_comm.allgather(&cnt, 1, counts.data());
        GISMO_ENSURE(0 == rc, "gsGeometricPartitioner::gatherLabels(): allgather of the own element counts failed.");
        int64_t total = 0;
        for (int q = 0; q != P; ++q) { displ[q] = static_cast<int>(total); total += counts[q]; }
        // counts is identical on every rank, so this check fires on all or none
        GISMO_ENSURE(total == static_cast<int64_t>(N),
                     "gsGeometricPartitioner::gatherLabels(): the own element lists do not cover the elements ("
                     << total << " of " << N << ").");
        std::vector<index_t> send(this->m_ownElements); // non-const buffer; MPI forbids aliasing
        std::vector<index_t> ids(N), labels(N);         // ids grouped by source rank q at displ[q]
        rc = m_comm.allgatherv(send.data(), cnt, ids.data(), counts.data(), displ.data());
        GISMO_ENSURE(0 == rc, "gsGeometricPartitioner::gatherLabels(): allgatherv of the own element lists failed.");
        for (int q = 0; q != P; ++q)
            for (int i = displ[q]; i != displ[q] + counts[q]; ++i)
                labels[ids[i]] = static_cast<index_t>(q);
        return labels;
    }

    /// @brief Per-element weight from pass A: all 1 by default, and the
    /// element's free-DOF count when Options::weightByDofs is set.
    /// Call partition() first. Serial path only: fails (GISMO_ENSURE) on a
    /// distributed partitioner, which keeps no per-element arrays.
    const std::vector<index_t>& elementWeights() const
    {
        GISMO_ENSURE(!this->distributed(),
                     "gsGeometricPartitioner::elementWeights(): not available on a distributed "
                     "partitioner (parallel constructor on more than one rank); use the serial constructor.");
        this->requireLabels();
        return m_weights;
    }

    /// @brief geoDim x numElements() physical centroids from pass A.
    /// Only available if Options::keepCentroids was set: the centroids are
    /// released at the end of partition() otherwise (they are geoDim*8
    /// bytes per element, which is the memory this class is built to save).
    /// Serial path only: fails (GISMO_ENSURE) on a distributed partitioner.
    const gsMatrix<T>& centroids() const
    {
        GISMO_ENSURE(!this->distributed(),
                     "gsGeometricPartitioner::centroids(): not available on a distributed "
                     "partitioner (parallel constructor on more than one rank); use the serial constructor.");
        GISMO_ENSURE(m_opts.keepCentroids,
                     "gsGeometricPartitioner::centroids(): the centroids are "
                     "released at the end of partition(); construct with "
                     "Options::keepCentroids = true to retain them.");
        return m_centroids;
    }

    /// @brief Number of elements of the partitioned domain.
    index_t numElements() const { return m_numElements; }

    /**
       @brief Pass B: graph-free DOF ownership -> gsPartitionedDofMapper.

       Streams the element mesh a second time and folds, for every free DOF,
       the minimum and maximum partition label of the elements active on it.
       That reduction is order-independent and idempotent under duplicates,
       so it reproduces gsPartitionedDofMapper's own (element -> DOF
       incidence) reduction exactly -- without materializing the incidence.

       Serial path only (it needs the full labels): fails (GISMO_ENSURE) on a
       distributed partitioner, where dof numbering is rendezvous-only.
    */
    gsPartitionedDofMapper makeDofMapper(index_t nranks) const override
    {
        GISMO_ENSURE(!this->distributed(),
                     "gsGeometricPartitioner::makeDofMapper(): not available on a distributed "
                     "partitioner (parallel constructor on more than one rank), whose dof "
                     "numbering is rendezvous-only; use the serial constructor.");
        this->requireLabels();

        const gsDofMapper&          mapper = this->mapper();
        const std::vector<index_t>& labels = this->labels();

        // Free global indices live in [base, base+nFree); minPart/maxPart are
        // indexed by the COMPACT free index g - base, which is what
        // gsPartitionedDofMapper::fromOwnershipRange() expects. base is
        // component-independent: all components share the one flat free range.
        const index_t base  = mapper.firstIndex();
        const index_t nFree = mapper.freeSize();
        const index_t nComp = mapper.numComponents();

        std::vector<index_t> minPart(nFree, std::numeric_limits<index_t>::max());
        std::vector<index_t> maxPart(nFree, -1);

        gsMatrix<index_t> locals, globals;
        gsMatrix<T>       centre;

        auto       it  = m_dom->beginAll();
        const auto end = m_dom->endAll();
        index_t expected = 0;
        for (; it != end; ++it, ++expected)
        {
            const index_t e = static_cast<index_t>(it.id());
            // GISMO_ENSURE, not GISMO_ASSERT: labels[e] and every array below
            // are indexed by id(), and this runs in Release at scale.
            GISMO_ENSURE(e == expected,
                "gsGeometricPartitioner: domain iterator id() is not contiguous "
                "in iteration order (got "<<e<<", expected "<<expected<<").");

            const index_t patch = static_cast<index_t>(it.patchIndex());
            const index_t p     = labels[e];

            centre = (it.lowerCorner() + it.upperCorner()) * T(0.5);
            this->multiBasis().piece(patch).active_into(centre, locals);

            // ALL components are walked, not just component 0. gsDofMapper
            // numbers a vector-valued space component-major over one flat free
            // range, and each component carries its own Dirichlet elimination,
            // so component c is NOT redundant with component 0. A
            // component-0-only pass would leave every DOF of components
            // 1..nComp-1 owned by no element -- which is invisible to
            // bijection/owned-count/cover checks and only shows up in
            // fromOwnershipRange()'s maxPart >= 0 test.
            for (index_t c = 0; c < nComp; ++c)
            {
                mapper.localToGlobal(locals, patch, globals, c);
                for (index_t i = 0; i < globals.size(); ++i)
                {
                    const index_t g = globals(i, 0);
                    if (mapper.is_free_index(g))
                    {
                        const index_t gc = g - base; // compact free index
                        if (p < minPart[gc]) minPart[gc] = p;
                        if (p > maxPart[gc]) maxPart[gc] = p;
                    }
                }
            }
        }

        return gsPartitionedDofMapper::fromOwnershipRange(minPart, maxPart,
                                                          nFree, nranks);
    }

    /// @brief Graph-free: there is no edge cut. Always -1 (drivers write it
    /// to the CSV's edgecut column so the schema stays fixed).
    index_t edgeCut() const { return -1; }

protected:

    /**
       @brief The base class' label hook: tensor-patch guard, pass A, then
       the selected label strategy.

       Called by gsPartitionerBase::partition(), which checks nparts before
       and, if not distributed, the label count afterwards; the label range
       is checked in setLabels(), the own ids in setOwnElements().
       Idempotent: every buffer used here is (re)sized from scratch.

       Serial and size-1 communicator: full labels over all N elements
       (setLabels()). Communicator with size() > 1: pass A on the rank's slice,
       then computeOwnElements(); the distributed RCB needs IEEE float/double
       coordinates, which is checked on every rank before any collective.
    */
    void computeLabels() override
    {
        const bool dist = m_comm.size() > 1;
        if (dist && hilbert != m_opts.strategy && morton != m_opts.strategy)
            GISMO_ENSURE(gsOrderedKey<T>::available, "gsGeometricPartitioner: the distributed RCB (parallel constructor, "
                         "comm.size() > 1) needs IEEE float or double coordinates; use the serial constructor or a "
                         "hilbert/morton strategy.");

        checkTensorPatches();
        passA();

        if (dist) { computeOwnElements(); return; }

        const index_t N = m_numElements;
        std::vector<index_t> labels(N, 0);

        // nparts == 1: every element is in part 0 (both strategies below
        // would produce exactly this, just less directly).
        if (this->nparts() > 1)
        {
            switch (m_opts.strategy)
            {
            case hilbert: labelsByCurve(labels, gsSpaceFillingCurve::Hilbert); break;
            case morton : labelsByCurve(labels, gsSpaceFillingCurve::Morton ); break;
            case rcb    :
            default     : labelsByRcb(labels);                                 break;
            }
        }

        this->setLabels(give(labels));

        // Release the centroids (geoDim*8 B/element: 400 MB at 16M elements
        // in 3D) unless the caller asked to keep them. m_weights (8 B/element)
        // stays: it is part of the public surface.
        if (!m_opts.keepCentroids)
        {
            gsMatrix<T> tmp;
            m_centroids.swap(tmp);
        }
    }

    /**
       @brief Distributed labelling of the rank's slice: runs the distributed
       RCB or curve machine through runCollectives() and stores the sorted ids
       of the elements of part m_comm.rank() (setOwnElements()).

       Requires pass A to have filled the slice centroids/weights, and
       nparts == m_comm.size() > 1 (part q is owned by rank q). Collectives:
       those of the machine, see gsDistributedRcb and gsDistributedCurve.
       Complexity and memory: O(n_loc) per rank for the keys and ids (plus
       O(P 256) select rows); nothing is sized by N. The slice centroids and
       weights are released afterwards, whatever Options::keepCentroids says.
    */
    void computeOwnElements()
    {
        const index_t N = m_numElements;
        const int64_t P = m_comm.size(), r = m_comm.rank();
        const index_t lo = sliceBegin(N, r, P);
        const index_t nLoc = static_cast<index_t>(m_centroids.cols());
        std::vector<index_t> own;
        switch (m_opts.strategy)
        {
        case hilbert:
        case morton:
        {
            std::vector<uint64_t> keys;
            curveKeys(hilbert == m_opts.strategy ? gsSpaceFillingCurve::Hilbert : gsSpaceFillingCurve::Morton, keys);
            gsDistributedCurve<T> mc(nLoc, lo, keys.data(), m_weights.data(), this->nparts(),
                                     static_cast<index_t>(r), static_cast<int64_t>(N));
            runCollectives(mc, m_comm);
            own = mc.releaseOwnElements();
            break;
        }
        case rcb:
        default:
        {
            gsDistributedRcb<T> mr(m_centroids, static_cast<short_t>(m_mp.geoDim()), lo, NULL, m_weights.data(),
                                   m_opts.weightByDofs, this->nparts(), static_cast<index_t>(r), static_cast<int64_t>(N));
            runCollectives(mr, m_comm);
            own = mr.releaseOwnElements();
            break;
        }
        }
        // The RCB machine reads the centroids and weights through pointers
        // until runCollectives returns, so they are released only here.
        { gsMatrix<T> tmp; m_centroids.swap(tmp); }
        std::vector<index_t>().swap(m_weights);
        this->setOwnElements(give(own), static_cast<index_t>(r), static_cast<index_t>(P));
    }

    /**
       @brief One bisection step of order[first,last) into \a k parts.

       Chooses the longest axis of the range's bounding box (ties to the
       lowest axis), reorders order[first,last) in place -- and nothing
       outside it -- and returns \c mid such that the left child is
       [first,mid) with kL = (k+1)/2 parts and the right child is [mid,last).
       With \a weighted the cut is found by a coordinate histogram on the
       weights; otherwise by an exact count-based split. The weighted case
       also falls back to the exact count split in two situations:
       -# a zero weight target, totalW*kL/k == 0 (i.e. totalW*kL < k,
          including totalW == 0); the histogram and the stable_partition are
          skipped;
       -# a degenerate histogram cut (mid == first or mid == last, n >= 2).

       Requires k >= 2 and first < last.

       Complexity: O(n) expected, n = last - first (std::nth_element). Weighted
       case: one O(n) weight sum; for a zero target that is followed directly
       by the nth_element split, otherwise by 16 O(n) histogram passes plus a
       std::stable_partition (O(n) with its temporary buffer, O(n log n) if
       the allocation fails), plus the nth_element split on a degenerate cut.
    */
    static index_t rcbBisect(const gsMatrix<T>& C, const std::vector<index_t>& w, bool weighted,
                             std::vector<index_t>& order, index_t first, index_t last, index_t k)
    {
        GISMO_ASSERT(k >= 2, "gsGeometricPartitioner::rcbBisect: k must be >= 2.");
        GISMO_ASSERT(first < last, "gsGeometricPartitioner::rcbBisect: empty range.");

        // ceil/floor, so nparts need not be a power of two
        const index_t kL = (k + 1) / 2;
        const index_t n = last - first;

        // Longest axis of the bounding box of this range; ties go to the
        // lowest axis index (strict >).
        const short_t geoDim = static_cast<short_t>(C.rows());
        short_t axis  = 0;
        T       bestE = -1, axLo = 0, axHi = 0;
        for (short_t d = 0; d != geoDim; ++d)
        {
            T lo = C(d, order[first]), hi = lo;
            for (index_t i = first + 1; i != last; ++i)
            {
                const T v = C(d, order[i]);
                if (v < lo) lo = v;
                if (v > hi) hi = v;
            }
            if (hi - lo > bestE) { bestE = hi - lo; axis = d; axLo = lo; axHi = hi; }
        }

        index_t mid = first;
        bool    exactSplit = !weighted;

        if (!exactSplit)
        {
            // Coordinate-histogram bisection on the weights: a bounded number
            // of O(n) passes, instead of an O(n log n) sort per level.
            int64_t totalW = 0;
            for (index_t i = first; i != last; ++i)
                totalW += static_cast<int64_t>(w[order[i]]);
            const int64_t target = (totalW * static_cast<int64_t>(kL))
                                 / static_cast<int64_t>(k);

            // A zero target is unreachable for the bisection (no midc satisfies
            // wL < 0): it would collapse to axLo + eps and isolate the elements
            // at the axis minimum. kL < k makes target == totalW possible only
            // for totalW == 0, so this test covers that case as well.
            if (0 == target)
                exactSplit = true;
            else
            {
                T lo = axLo, hi = axHi;
                for (int iter = 0; iter != 16; ++iter)
                {
                    const T midc = T(0.5) * (lo + hi);
                    int64_t wL = 0;
                    for (index_t i = first; i != last; ++i)
                        if (C(axis, order[i]) < midc)
                            wL += static_cast<int64_t>(w[order[i]]);
                    if (wL < target) lo = midc; else hi = midc;
                }
                const T cut = T(0.5) * (lo + hi);

                // Order-preserving partition: the outcome then depends on the data
                // only, never on the incoming order of `order`.
                typename std::vector<index_t>::iterator pivot =
                    std::stable_partition(order.begin() + first, order.begin() + last,
                        [&](index_t a) { return C(axis, a) < cut; });
                mid = static_cast<index_t>(pivot - order.begin());

                // Degenerate split (e.g. every centroid shares this coordinate):
                // fall back to the exact count-based split, otherwise one part
                // would swallow the whole range.
                if ((mid == first || mid == last) && n >= 2) exactSplit = true;
            }
        }

        if (exactSplit)
        {
            int64_t nL64 = (static_cast<int64_t>(n) * static_cast<int64_t>(kL)
                            + static_cast<int64_t>(k / 2))
                         / static_cast<int64_t>(k);
            if (nL64 < 0) nL64 = 0;
            if (nL64 > static_cast<int64_t>(n)) nL64 = n;
            const index_t nL = static_cast<index_t>(nL64);

            // Strict TOTAL order: coordinate along `axis`, ties broken on the
            // element id. std::nth_element is not stable, and the own-branch
            // descents of different ranks share the upper levels of the
            // recursion: every rank must agree bit for bit on every split it
            // evaluates, so a merely weak ordering would be a silent
            // correctness bug, not a quality nit.
            const short_t ax = axis;
            std::nth_element(order.begin() + first, order.begin() + first + nL,
                             order.begin() + last,
                [&C, ax](index_t a, index_t b)
                {
                    const T ca = C(ax, a), cb = C(ax, b);
                    return (ca < cb) || (ca == cb && a < b);
                });
            mid = first + nL;
        }
        return mid;
    }

    /**
       @brief Recursive coordinate bisection of order[first,last) into \a k
       parts labelled label0 .. label0+k-1, written to \a labels.

       Replicated recursion: visits every leaf. Empty ranges leave their
       parts empty, which is legal (gsIndexSubDomain and PETSc both accept an
       empty part).

       Complexity: O(n) expected per level (std::nth_element, or ~16 O(n)
       histogram passes in the weighted case) and O(log k) levels, i.e.
       O(N log k) expected overall. No sort per level.
    */
    static void rcbSplit(const gsMatrix<T>& C, const std::vector<index_t>& w, bool weighted,
                         std::vector<index_t>& labels, std::vector<index_t>& order,
                         index_t first, index_t last, index_t k, index_t label0)
    {
        if (1 == k)
        {
            for (index_t i = first; i != last; ++i) labels[order[i]] = label0;
            return;
        }

        const index_t kL = (k + 1) / 2;
        const index_t kR = k - kL;

        if (first == last) return;

        const index_t mid = rcbBisect(C, w, weighted, order, first, last, k);
        rcbSplit(C, w, weighted, labels, order, first, mid,  kL, label0);
        rcbSplit(C, w, weighted, labels, order, mid,   last, kR, label0 + kL);
    }

    /**
       @brief The leaf of the RCB recursion that carries label \a part:
       follows only the child whose label range contains \a part and returns
       the leaf range [first,last) of \a order. An empty range is legal.

       Precondition: label0 <= part < label0 + k.

       The leaf equals the one rcbSplit() produces for the same inputs:
       rcbBisect() permutes only order[first,last), so processing the sibling
       subtree in rcbSplit() never changes the contents or the order of the
       range descended into here; both recursions call rcbBisect() on
       identical inputs at every node on the path to \a part.

       Complexity: the sum of the range sizes along one root-to-leaf path,
       about N + N/2 + ... = 2N element visits for balanced (unweighted)
       splits, i.e. O(N) per call; worst case O(N ceil(log2 k)) when weighted
       splits are very unbalanced.
    */
    static std::pair<index_t,index_t> rcbOwnLeaf(const gsMatrix<T>& C, const std::vector<index_t>& w,
                         bool weighted, std::vector<index_t>& order, index_t first, index_t last,
                         index_t k, index_t label0, index_t part)
    {
        GISMO_ASSERT(label0 <= part && part < label0 + k,
                     "gsGeometricPartitioner::rcbOwnLeaf: part outside the label range.");
        while (k > 1 && first < last)
        {
            const index_t mid = rcbBisect(C, w, weighted, order, first, last, k);
            const index_t kL  = (k + 1) / 2;
            if (part < label0 + kL) { last = mid;                k  = kL; }
            else                    { first = mid; label0 += kL; k -= kL; }
        }
        return std::make_pair(first, last);
    }

    /**
       @brief Labels the elements by cutting the space-filling curve with
       keys \a keys into \a nparts parts of (nearly) equal weight.

       \a keys[e] and \a w[e] belong to element e (both of length N, indexed
       by element id); on return \a labels has size N and labels[e] is in
       [0, nparts). Precondition: w[e] >= 0 for every e.

       The elements are sorted by (key, id), a strict total order. With A(e)
       the exclusive prefix weight of e in that order, W = sum of w and
       P = nparts, the label is
       \f[ \mathrm{label}(e) = \#\{ j \in [1,P-1] : P\,A(e) \ge W\,j \}, \f]
       which is valid because A is non-decreasing along the curve. Labels are
       non-decreasing along the curve; a heavy element can make them skip, so
       some parts may be empty. If W == 0 (all weights zero) every element
       weighs 1 instead (W = N), the count split, as in rcbBisect().

       Exact int64 arithmetic (acc * P stays far below 2^63 for index_t
       weights at realistic element counts and weights).

       Complexity: O(N log N) for the sort plus O(N + P) for the cut, O(N)
       extra memory.
    */
    static void curveLabels(const std::vector<uint64_t>& keys, const std::vector<index_t>& w,
                            index_t nparts, std::vector<index_t>& labels)
    {
        GISMO_ASSERT(w.size() == keys.size(),
                     "gsGeometricPartitioner::curveLabels: keys and weights differ in size.");
        GISMO_ASSERT(nparts >= 1, "gsGeometricPartitioner::curveLabels: nparts must be >= 1.");

        const index_t N = static_cast<index_t>(keys.size());
        labels.resize(keys.size());

        std::vector<index_t> order(N);
        for (index_t i = 0; i != N; ++i) order[i] = i;

        // Strict TOTAL order: curve key, ties broken on the element id.
        // std::sort is not stable, and the labels are never broadcast --
        // every rank recomputes them independently and must agree bit for
        // bit, so a weak ordering would silently diverge DOF ownership.
        std::sort(order.begin(), order.end(),
            [&keys](index_t a, index_t b)
            { return (keys[a] < keys[b]) || (keys[a] == keys[b] && a < b); });

        // Cut the curve by weighted prefix sum, in exact integer arithmetic
        // (acc * nparts reaches ~2.2e12 at realistic sizes, hence int64_t).
        int64_t totalW = 0;
        for (index_t e = 0; e != N; ++e) totalW += static_cast<int64_t>(w[e]);
        // All-zero weights: unit weights, i.e. the count split. With W == 0
        // every test P*A >= W*j would hold and all elements would get P-1.
        const bool unit = (totalW == 0);
        if (unit) totalW = static_cast<int64_t>(N);

        const int64_t nparts64 = static_cast<int64_t>(nparts);
        int64_t acc = 0;
        index_t j   = 0;
        for (index_t pos = 0; pos != N; ++pos)
        {
            const index_t e = order[pos];
            while (j + 1 < nparts &&
                   acc * nparts64 >= totalW * static_cast<int64_t>(j + 1))
                ++j;
            labels[e] = j;                       // non-decreasing along the curve
            acc += unit ? 1 : static_cast<int64_t>(w[e]);
        }
    }

private:

    // ------------------------------------------------------------------
    // Preconditions
    // ------------------------------------------------------------------

    /// @brief Behavioural check that each patch domain's volume iterator
    /// reports its own patch index -- exactly the property pass A and pass B
    /// depend on (a type check against gsTensorDomain would be both narrower
    /// and less future-proof).
    void checkTensorPatches() const
    {
        const index_t nP = static_cast<index_t>(m_dom->nPieces());
        for (index_t p = 0; p < nP; ++p)
            GISMO_ENSURE(m_dom->subdomain(p)->beginAll().patchIndex() == p,
                "gsGeometricPartitioner: patch "<<p<<"'s domain iterator reports "
                "patchIndex() == "<<m_dom->subdomain(p)->beginAll().patchIndex()<<
                ", not "<<p<<". Only tensor-product patch domains propagate a patch "
                "index to the volume iterator (gsTensorDomainIterator.h); "
                "gsHDomainIterator (THB) and gsKnotDomainIterator (1-D) do not, so "
                "element->patch attribution and hence active_into()/localToGlobal() "
                "would silently run on patch 0 for every patch. Graph-free "
                "partitioning requires tensor-product patch domains.");
        // Note: for nPieces() == 1 this passes for any domain type (a
        // non-propagating iterator reports 0, and patch 0 is the only patch),
        // i.e. single-patch THB / 1-D geometries are correctly permitted.
    }

    /// @brief Agree collectively on a failure that may have occurred on a
    /// subset of the ranks (\a failure non-null on the failing ranks), so no
    /// rank is left waiting in a later collective call. Throws on all ranks.
    void ensureOnAllRanks(const char * failure) const
    {
        int bad = failure ? 1 : 0;
        bad = m_comm.max(bad);
        if (failure)
            gsWarn << "gsGeometricPartitioner, rank " << m_comm.rank() << ": " << failure << "\n";
        GISMO_ENSURE(0 == bad, "gsGeometricPartitioner: a check failed on at least "
                     "one rank, see the warning of that rank.");
    }

    /// First element id of slice q of P over N elements: floor(N q / P), in int64 (N*q overflows index_t).
    static index_t sliceBegin(index_t N, int64_t q, int64_t P)
    { return static_cast<index_t>(static_cast<int64_t>(N) * q / P); }

    /// Curve key of every column of m_centroids, keys[i] = key of column i
    /// (keys.size() == m_centroids.cols()). O(cols * geoDim) plus one curve
    /// encoding per column.
    void curveKeys(gsSpaceFillingCurve::Curve curve, std::vector<uint64_t>& keys) const
    {
        const short_t geoDim = static_cast<short_t>(m_centroids.rows());

        // gsSpaceFillingCurve is not templated: it works in real_t, so the
        // box and the points cross a T -> real_t boundary here.
        gsMatrix<real_t> box = m_bbox.template cast<real_t>();
        // The curve owns the bit budget: bits() == bitsPerAxis(curveDim()),
        // keyed on the number of NON-degenerate axes (a flat plate in 3D gets
        // 31 bits per axis, not 21). Do not recompute it from geoDim.
        const gsSpaceFillingCurve sfc(box, curve);

        keys.assign(m_centroids.cols(), 0);
        gsVector<real_t> pt(geoDim);
        for (index_t i = 0; i != static_cast<index_t>(m_centroids.cols()); ++i)
        {
            for (short_t d = 0; d != geoDim; ++d)
                pt[d] = static_cast<real_t>(m_centroids(d, i));
            keys[i] = sfc.encode(pt);
        }
    }

    // ------------------------------------------------------------------
    // Pass A: centroids + weights
    // ------------------------------------------------------------------

    /**
       @brief One physical centroid and one integer weight per element.

       The geometry evaluation is batched: gsDomainIterator::centerPoint()
       costs 3+2*parDim heap allocations per element already, so a
       per-element geometry evaluation on top of that would dominate the
       pass. Parametric centres are buffered and evaluated once per
       (patch, chunk) instead.

       Iteration is over beginAll()..endAll(), which -- unlike
       gsDomain::allElements() -- carries no OpenMP chunking, so this plain
       serial loop is correct and complete.

       Serial (and size <= 1 communicator): N centroids, no communication.
       With a communicator of size P > 1 rank r visits only its contiguous
       slice [N r/P, N (r+1)/P) (O(N/P) evaluations after an
       O(dim + log nPatches) iterator jump) and stores only that slice,
       O(n_loc (geoDim sizeof(T) + sizeof(index_t))) memory. One scalar
       max-reduction (ensureOnAllRanks) is the only collective: it agrees on
       the iterator-id, slice-coverage and finite-centroid checks.
    */
    void passA()
    {
        const index_t N      = m_numElements;
        const short_t parDim = m_mp.parDim();
        const short_t geoDim = m_mp.geoDim();

        // Element slice [lo, hi) handled by this rank; the whole range when serial.
        const bool    par = m_comm.size() > 1;
        const int64_t P   = par ? m_comm.size() : 1;
        const int64_t r   = par ? m_comm.rank() : 0;
        const index_t lo  = sliceBegin(N, r, P);
        const index_t hi  = sliceBegin(N, r + 1, P);
        const index_t chunk = math::min( math::max((index_t)1, m_opts.chunkSize),
                                         math::max((index_t)1, hi - lo) );

        // Bounding box of the geometry, computed once here and reused by the
        // space-filling-curve strategies.
        //
        // The finiteness check is mandatory: gsSpaceFillingCurve's
        // degenerate-axis test is extent[i] > relTol * maxExtent, which with
        // maxExtent == inf is false for EVERY axis, so curveDim() collapses to
        // 0 and every element is keyed 0 -- with no diagnostic at any level.
        // A non-finite coordinate breaks the RCB median split just as badly,
        // hence the check sits on the shared pass-A path.
        m_mp.boundingBox(m_bbox);
        GISMO_ENSURE(m_bbox.allFinite(),
                     "gsGeometricPartitioner: the multipatch bounding box is not finite ("
                     << m_bbox.transpose() << "). A non-finite extent silently collapses the "
                     "space-filling curve to a single cell (every element keyed 0), so this "
                     "is refused rather than partitioned.");

        m_centroids.resize(geoDim, hi - lo);   // column e - lo = centroid of element e
        m_weights.assign(hi - lo, 1);          // entry e - lo = weight of element e

        const gsDofMapper& mapper = this->mapper();
        const index_t      nComp  = mapper.numComponents();

        gsMatrix<T>          params(parDim, chunk), phys, u, centre;
        std::vector<index_t> bufElem(chunk);
        gsMatrix<index_t>    locals, globals;

        index_t curPatch = -1, nBuf = 0;

        // One geometry evaluation per (patch, chunk).
        auto flush = [&]()
        {
            if (0 == nBuf) return;
            u = params.leftCols(nBuf);           // materialize the block
            m_mp.patch(curPatch).eval_into(u, phys);
            // bufElem (global ids) is contiguous by construction: `expected`
            // increments by one per element and the id check below rejects
            // any gap, while a chunk is flushed before it can span two patches.
            // So bufElem[k] == bufElem[0] + k and the whole copy is one block
            // assignment. Measured 2026-08-13: the per-element form put ~38% of
            // pass A into Eigen Block construction and dense-assignment loops
            // (perf, self time, pass A isolated).
            m_centroids.middleCols(bufElem[0] - lo, nBuf) = phys;
            nBuf = 0;
        };

        auto       it  = m_dom->beginAll();
        const auto end = m_dom->endAll();
        if (lo > 0 && lo < N) it += lo;
        index_t expected = lo;
        std::string failure;
        for (; it != end && expected < hi; ++it, ++expected)
        {
            const index_t e = static_cast<index_t>(it.id());
            if (par && e != expected)
            {
                // Throwing here would leave the other ranks waiting in the
                // collective below; the failure is agreed on collectively.
                std::ostringstream os;
                os << "domain iterator id() is not contiguous in iteration order (got "
                   << e << ", expected " << expected << ")";
                failure = os.str();
                break;
            }
            // GISMO_ENSURE, not GISMO_ASSERT: every array here is indexed by
            // id(), and this runs in Release at scale, where a GISMO_ASSERT is
            // compiled out and a non-contiguous id corrupts silently.
            GISMO_ENSURE(e == expected,
                "gsGeometricPartitioner: domain iterator id() is not contiguous "
                "in iteration order (got "<<e<<", expected "<<expected<<").");

            const index_t patch = static_cast<index_t>(it.patchIndex());

            // A chunk may never span two patches (the evaluation is per-patch).
            if (nBuf == chunk || (nBuf > 0 && patch != curPatch)) flush();
            curPatch = patch;

            params.col(nBuf) = (it.lowerCorner() + it.upperCorner()) * T(0.5);
            bufElem[nBuf]    = e;

            if (m_opts.weightByDofs)
            {
                // active_into() takes a gsMatrix, not an expression, so this
                // branch -- and only this branch -- still materializes the
                // centre. The default path writes straight into params.
                centre = params.col(nBuf);

                // Same machinery as pass B -- all components, free indices
                // only -- but the DOF ids are counted, never stored.
                this->multiBasis().piece(patch).active_into(centre, locals);
                index_t cnt = 0;
                for (index_t c = 0; c < nComp; ++c)
                {
                    mapper.localToGlobal(locals, patch, globals, c);
                    for (index_t i = 0; i < globals.size(); ++i)
                        if (mapper.is_free_index(globals(i, 0))) ++cnt;
                }
                m_weights[e - lo] = cnt;
            }
            ++nBuf;
        }
        flush();

        if (par)
        {
            if (failure.empty() && expected != hi)
                failure = "domain iterator ended before the element slice was covered";
            if (failure.empty() && !m_centroids.allFinite())
                failure = "the geometry evaluated to a non-finite element centroid in this rank's slice";
            ensureOnAllRanks(failure.empty() ? nullptr : failure.c_str());
        }

        // The box above is the control net's; a rational patch can still
        // evaluate to a non-finite point from a finite net, and such a
        // centroid would break both the curve quantisation and the RCB
        // median split. One O(n_loc) pass, next to nothing beside the evaluation.
        GISMO_ENSURE(m_centroids.allFinite(),
                     "gsGeometricPartitioner: the geometry evaluated to a "
                     "non-finite element centroid; refusing to partition.");
    }

    // ------------------------------------------------------------------
    // Label strategy: recursive coordinate bisection
    // ------------------------------------------------------------------

    /**
       @brief RCB labels, replicated: serial constructor and communicators of
       size <= 1 only (a distributed partitioner uses computeOwnElements()).
       The recursion rcbSplit() over all leaves, O(N log nparts).
    */
    void labelsByRcb(std::vector<index_t>& labels) const
    {
        GISMO_ASSERT(m_comm.size() <= 1, "gsGeometricPartitioner::labelsByRcb: replicated RCB only.");
        const index_t N = static_cast<index_t>(labels.size());
        std::vector<index_t> order(N);
        for (index_t i = 0; i != N; ++i) order[i] = i;

        rcbSplit(m_centroids, m_weights, m_opts.weightByDofs, labels, order,
                 0, N, this->nparts(), 0);
    }

    // ------------------------------------------------------------------
    // Label strategy: Hilbert / Morton space-filling curve
    // ------------------------------------------------------------------

    void labelsByCurve(std::vector<index_t>& labels,
                       gsSpaceFillingCurve::Curve curve) const
    {
        std::vector<uint64_t> keys;
        curveKeys(curve, keys);
        curveLabels(keys, m_weights, this->nparts(), labels);
    }

private:

    const gsMultiPatch<T>&    m_mp;
    const Options             m_opts;
    typename gsDomain<T>::Ptr m_dom;
    index_t                   m_numElements; ///< element count of m_dom, cached (the composite domain recounts in O(nPieces))
    gsMpiComm                 m_comm;        ///< copy of the pass-A communicator; size() <= 1 means serial

    gsMatrix<T>          m_bbox;      ///< geoDim x 2 (lower, upper corner)
    gsMatrix<T>          m_centroids; ///< geoDim x numElements on the serial path, geoDim x (hi - lo) when distributed (column i = element lo + i, lo = floor(N rank/P)); released after partition()
    std::vector<index_t> m_weights;   ///< one weight per element, indexed like the columns of m_centroids

}; // class gsGeometricPartitioner

} // namespace gismo
