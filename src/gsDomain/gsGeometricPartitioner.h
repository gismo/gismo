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
    subdomainForRank() are shared, and a caller can switch between them.

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
#include <gsDomain/gsIndexSubDomain.h>
#include <gsDomain/gsPartitionedDofMapper.h>
#include <gsDomain/gsPartitionerBase.h>
#include <gsUtils/gsSpaceFillingCurve.h>

#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
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

   \note Every rank is expected to run this partitioner independently and to
   obtain bit-identical labels (there is no broadcast of the labelling), so
   every ordering used internally is a strict \em total order: the primary
   key is compared first and every tie is broken on the element id.

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
        bool     keepCentroids = false; ///< retain centroids() after partition() (geoDim*8 B/element)
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
    : Base(mb, mapper, nparts), m_mp(mp), m_opts(opts), m_dom(mb.domain())
    { }

    /// @brief "rcb" | "hilbert" | "morton"; fails on anything else.
    static Strategy strategyFromString(const std::string& s)
    {
        if ("rcb"     == s) return rcb;
        if ("hilbert" == s) return hilbert;
        if ("morton"  == s) return morton;
        GISMO_ERROR("gsGeometricPartitioner: unknown strategy \""<<s<<"\" "
                    "(expected \"rcb\", \"hilbert\" or \"morton\").");
    }

    /// @brief Per-element weight from pass A: all 1 by default, and the
    /// element's free-DOF count when Options::weightByDofs is set.
    /// Call partition() first.
    const std::vector<index_t>& elementWeights() const
    {
        this->requireLabels();
        return m_weights;
    }

    /// @brief geoDim x numElements() physical centroids from pass A.
    /// Only available if Options::keepCentroids was set: the centroids are
    /// released at the end of partition() otherwise (they are geoDim*8
    /// bytes per element, which is the memory this class is built to save).
    const gsMatrix<T>& centroids() const
    {
        GISMO_ENSURE(m_opts.keepCentroids,
                     "gsGeometricPartitioner::centroids(): the centroids are "
                     "released at the end of partition(); construct with "
                     "Options::keepCentroids = true to retain them.");
        return m_centroids;
    }

    /// @brief Number of elements of the partitioned domain.
    index_t numElements() const { return static_cast<index_t>(m_dom->numElements()); }

    /**
       @brief Pass B: graph-free DOF ownership -> gsPartitionedDofMapper.

       Streams the element mesh a second time and folds, for every free DOF,
       the minimum and maximum partition label of the elements active on it.
       That reduction is order-independent and idempotent under duplicates,
       so it reproduces gsPartitionedDofMapper's own (element -> DOF
       incidence) reduction exactly -- without materializing the incidence.
    */
    gsPartitionedDofMapper makeDofMapper(index_t nranks) const override
    {
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
       and the label count/range afterwards. Idempotent: every buffer used
       here is (re)sized from scratch.
    */
    void computeLabels() override
    {
        checkTensorPatches();
        passA();

        const index_t N = static_cast<index_t>(m_dom->numElements());
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
    */
    void passA()
    {
        const index_t N      = static_cast<index_t>(m_dom->numElements());
        const short_t parDim = m_mp.parDim();
        const short_t geoDim = m_mp.geoDim();
        const index_t chunk  = math::min( math::max((index_t)1, m_opts.chunkSize),
                                          math::max((index_t)1, N) );

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

        m_centroids.resize(geoDim, N);   // column e = centroid of element e
        m_weights.assign(N, 1);

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
            // bufElem is contiguous by construction: `expected` increments by
            // one per element and the GISMO_ENSURE(e == expected) below rejects
            // any gap, while a chunk is flushed before it can span two patches.
            // So bufElem[k] == bufElem[0] + k and the whole copy is one block
            // assignment. Measured 2026-08-13: the per-element form put ~38% of
            // pass A into Eigen Block construction and dense-assignment loops
            // (perf, self time, pass A isolated).
            m_centroids.middleCols(bufElem[0], nBuf) = phys;
            nBuf = 0;
        };

        auto       it  = m_dom->beginAll();
        const auto end = m_dom->endAll();
        index_t expected = 0;
        for (; it != end; ++it, ++expected)
        {
            const index_t e = static_cast<index_t>(it.id());
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
                m_weights[e] = cnt;
            }
            ++nBuf;
        }
        flush();

        // The box above is the control net's; a rational patch can still
        // evaluate to a non-finite point from a finite net, and such a
        // centroid would break both the curve quantisation and the RCB
        // median split. One O(N) pass, next to nothing beside the evaluation.
        GISMO_ENSURE(m_centroids.allFinite(),
                     "gsGeometricPartitioner: the geometry evaluated to a "
                     "non-finite element centroid; refusing to partition.");
    }

    // ------------------------------------------------------------------
    // Label strategy: recursive coordinate bisection
    // ------------------------------------------------------------------

    void labelsByRcb(std::vector<index_t>& labels) const
    {
        const index_t N = static_cast<index_t>(labels.size());
        std::vector<index_t> order(N);
        for (index_t i = 0; i != N; ++i) order[i] = i;
        rcbSplit(labels, order, 0, N, this->nparts(), 0);
    }

    /**
       @brief Recursive coordinate bisection of order[first,last) into \a k
       parts labelled label0 .. label0+k-1.

       Complexity: O(n) expected per level (std::nth_element, or ~16 O(n)
       histogram passes in the weighted case) and O(log k) levels, i.e.
       O(N log k) expected overall. No sort per level.
    */
    void rcbSplit(std::vector<index_t>& labels, std::vector<index_t>& order,
                  index_t first, index_t last, index_t k, index_t label0) const
    {
        if (1 == k)
        {
            for (index_t i = first; i != last; ++i) labels[order[i]] = label0;
            return;
        }

        // ceil/floor, so nparts need not be a power of two
        const index_t kL = (k + 1) / 2;
        const index_t kR = k - kL;

        // Empty range: parts label0..label0+k-1 simply stay empty, which is
        // legal (gsIndexSubDomain and PETSc both accept an empty part).
        if (first == last) return;

        const index_t n = last - first;

        // Longest axis of the bounding box of this range; ties go to the
        // lowest axis index (strict >).
        const short_t geoDim = static_cast<short_t>(m_centroids.rows());
        short_t axis  = 0;
        T       bestE = -1, axLo = 0, axHi = 0;
        for (short_t d = 0; d != geoDim; ++d)
        {
            T lo = m_centroids(d, order[first]), hi = lo;
            for (index_t i = first + 1; i != last; ++i)
            {
                const T v = m_centroids(d, order[i]);
                if (v < lo) lo = v;
                if (v > hi) hi = v;
            }
            if (hi - lo > bestE) { bestE = hi - lo; axis = d; axLo = lo; axHi = hi; }
        }

        index_t mid = first;
        bool    exactSplit = !m_opts.weightByDofs;

        if (!exactSplit)
        {
            // Coordinate-histogram bisection on the weights: a bounded number
            // of O(n) passes, instead of an O(n log n) sort per level.
            int64_t totalW = 0;
            for (index_t i = first; i != last; ++i)
                totalW += static_cast<int64_t>(m_weights[order[i]]);
            const int64_t target = (totalW * static_cast<int64_t>(kL))
                                 / static_cast<int64_t>(k);

            T lo = axLo, hi = axHi;
            for (int iter = 0; iter != 16; ++iter)
            {
                const T midc = T(0.5) * (lo + hi);
                int64_t wL = 0;
                for (index_t i = first; i != last; ++i)
                    if (m_centroids(axis, order[i]) < midc)
                        wL += static_cast<int64_t>(m_weights[order[i]]);
                if (wL < target) lo = midc; else hi = midc;
            }
            const T cut = T(0.5) * (lo + hi);

            // Order-preserving partition: the outcome then depends on the data
            // only, never on the incoming order of `order`.
            typename std::vector<index_t>::iterator pivot =
                std::stable_partition(order.begin() + first, order.begin() + last,
                    [&](index_t a) { return m_centroids(axis, a) < cut; });
            mid = static_cast<index_t>(pivot - order.begin());

            // Degenerate split (e.g. every centroid shares this coordinate):
            // fall back to the exact count-based split, otherwise one part
            // would swallow the whole range.
            if ((mid == first || mid == last) && n >= 2) exactSplit = true;
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
            // element id. std::nth_element is not stable, and there is no
            // broadcast of the labels -- every rank recomputes this partition
            // independently and must agree bit for bit, so a merely weak
            // ordering would be a silent correctness bug, not a quality nit.
            const gsMatrix<T>& C = m_centroids;
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

        rcbSplit(labels, order, first, mid,  kL, label0);
        rcbSplit(labels, order, mid,   last, kR, label0 + kL);
    }

    // ------------------------------------------------------------------
    // Label strategy: Hilbert / Morton space-filling curve
    // ------------------------------------------------------------------

    void labelsByCurve(std::vector<index_t>& labels,
                       gsSpaceFillingCurve::Curve curve) const
    {
        const index_t N      = static_cast<index_t>(labels.size());
        const short_t geoDim = static_cast<short_t>(m_centroids.rows());

        // gsSpaceFillingCurve is not templated: it works in real_t, so the
        // box and the points cross a T -> real_t boundary here.
        gsMatrix<real_t> box = m_bbox.template cast<real_t>();
        // The curve owns the bit budget: bits() == bitsPerAxis(curveDim()),
        // keyed on the number of NON-degenerate axes (a flat plate in 3D gets
        // 31 bits per axis, not 21). Do not recompute it from geoDim.
        const gsSpaceFillingCurve sfc(box, curve);

        std::vector<uint64_t> keys(N, 0);
        gsVector<real_t>      pt(geoDim);
        for (index_t e = 0; e != N; ++e)
        {
            for (short_t d = 0; d != geoDim; ++d)
                pt[d] = static_cast<real_t>(m_centroids(d, e));
            keys[e] = sfc.encode(pt);
        }

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
        // (acc * nparts reaches ~2.2e12 at plan scale, hence int64_t).
        int64_t totalW = 0;
        for (index_t e = 0; e != N; ++e) totalW += static_cast<int64_t>(m_weights[e]);

        const int64_t nparts64 = static_cast<int64_t>(this->nparts());
        int64_t acc = 0;
        index_t j   = 0;
        for (index_t pos = 0; pos != N; ++pos)
        {
            const index_t e = order[pos];
            while (j + 1 < this->nparts() &&
                   acc * nparts64 >= totalW * static_cast<int64_t>(j + 1))
                ++j;
            labels[e] = j;                       // non-decreasing along the curve
            acc += static_cast<int64_t>(m_weights[e]);
        }
    }

private:

    const gsMultiPatch<T>&    m_mp;
    const Options             m_opts;
    typename gsDomain<T>::Ptr m_dom;

    gsMatrix<T>          m_bbox;      ///< geoDim x 2 (lower, upper corner)
    gsMatrix<T>          m_centroids; ///< geoDim x numElements (released after partition())
    std::vector<index_t> m_weights;   ///< one weight per element

}; // class gsGeometricPartitioner

} // namespace gismo
