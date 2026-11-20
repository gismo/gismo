/** @file gsHPatchRule.h

    @brief Provides hierarchical patch-wise quadrature rule (per-atomic-cell).

    For each leaf (atomic cell) of the hierarchical domain, a tensorized
    Johannessen Newton quadrature rule is built in the parameter domain of
    that cell.  The \ref mapTo method looks up the atomic cell containing a
    queried element and returns the subset of quadrature points inside it.

    The class is dimension-agnostic (templated only on the scalar type \a T,
    like \ref gsPatchRule) and is constructed from a \ref gsDomain that must be
    a hierarchical domain (\ref gsHDomain).  This lets it be created by the
    \ref gsQuadrature factory and used directly through \ref gsExprAssembler.

    Construction parameters:
      - \a degree    : target degree of the integration space. By convention this
        should equal 2*p, where p is the degree of the underlying hierarchical
        basis in the corresponding direction (use 2*p for mass); the class does
        not enforce this, so a caller that violates it silently changes which
        bilinear form the rule is exact for.
      - \a regularity: target regularity **at a simple interior knot** (use p-1
        for max-continuity mass; p-2 for stiffness). At an interior breakpoint
        of multiplicity mu_j > 1 in the level knot vector, the effective target
        regularity used for the leaf's target knot vector is offset to
        \a regularity + 1 - mu_j (floored at -1, i.e. target multiplicity capped
        at \a degree + 1), so that the requested space is exact regardless of
        how smooth the underlying basis actually is there.
      - \a overintegrate : not yet implemented (must be false)

    The choice between the Johannessen rule and element-wise Gauss is made per
    leaf and per direction (other than the fixed direction). A direction uses
    element-wise Gauss with floor(\a degree / 2) + 1 points per element, the
    smallest count that is exact to degree \a degree, whenever:
      - its candidate Johannessen node count is >= the element-wise Gauss node
        count of the leaf in that direction (in effect a tie, which costs the
        same and skips the Newton solve), or
      - its Newton solve fails to converge.
    The other directions of the same leaf keep their Johannessen rule. A tensor
    product of such per-direction rules is exact on the leaf's target space.
    numFallbackCells() counts the leaves with a direction in Gauss because its
    Newton solve failed; a tie is not counted.

    Each 1D rule is solved on the leaf's target knot vector mapped to [0,1] and
    scaled to the leaf (node a + L xi, weight L omega). The solve on [0,1] stops at
    1e-14/min(L,1): this is the physical stop 1e-14 for a leaf direction of length
    L <= 1 and the stop of a length-1 leaf for L > 1 (see the tol parameter of
    gsHPatchRule::_compute1D). A @ref gsHPatchRuleCache shares the
    solved shapes between leaves, directions and builds.

    Leaf merging (gsHPatchRule::MergeMode). The rule can replace the leaves of each level
    by fewer boxes that tile the same region. Every such box is a union of leaves of its
    level, so the rule of the box is exact on the tensor space of that level on the box by
    the argument of the class documentation, and it has fewer nodes than the rules of the
    leaves it replaces. Everything above holds per cell, where a cell is a leaf (no merging)
    or such a box.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsAssembler/gsQuadRule.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsDomain/gsHDomain.h>

#include <array>
#include <functional>
#include <map>
#include <utility>
#include <vector>

namespace gismo
{

template<short_t d, class T, bool Trunc> class gsTHBSplineBasis;

template<short_t d, class T> class gsHTensorBasis;

template<class T> class gsHPatchRule;

/**
    @brief Caller-owned store of solved 1D quadrature shapes for @ref gsHPatchRule.

    A shape is the Johannessen rule of one even target knot vector mapped to [0,1]. A
    @ref gsHPatchRule solves every leaf direction on [0,1] and scales the result to the leaf
    (node \f$a + L\xi\f$, weight \f$L\omega\f$, \f$L = b-a\f$), so leaf directions with
    the same shape need one solve only. Passing the same cache to several builds shares
    the shapes between them.

    **Key.** The clamped length \f$\min(L,1)\f$ together with the degree and the full
    repeated-knot sequence of the normalised even target knot vector, all compared
    exactly (no tolerance, no quantisation). This is exactly the input of the solve: the
    knot vector on [0,1] and the stopping tolerance \f$10^{-14}/\min(L,1)\f$ on the
    2-norm of the moment residual. A hit therefore returns bitwise what a fresh solve
    of the same key returns. Leaf directions of length at least 1 with the same normalised
    shape share an entry; a direction shorter than 1 shares with directions of the same
    length only. Two leaves whose normalised knots or lengths differ by one ulp get two
    entries, which costs a solve and never a wrong rule.

    **Value.** The nodes and weights on [0,1], or a recorded failure (a solve that did not
    converge). A failed entry holds no nodes and weights; the partial iterates of a
    failed solve are never stored.

    **Counters.** For every lookup of a patch direction (a non-fixed direction whose
    Johannessen node count is below the Gauss node count) a build increments either
    numHits() (the key is present, converged or failed) or numSolves() (the key is
    absent, so a solve is made and stored), and a failed solve also increments
    numFailures() and appends its key to failedKeys(). After any number of builds:
    numHits() + numSolves() is the number of lookups, numFailures() == failedKeys().size()
    <= numSolves(), and size() == numSolves(). Directions that use Gauss because of a tie
    never touch the cache.

    **Thread safety and lifetime.** The cache is not thread-safe: two builds that run
    concurrently must not share one cache. A rule does not keep a pointer to the cache,
    so the cache only has to outlive the constructor call it is passed to. A factory from
    gsHPatchRule::makeFactory is the exception to the lifetime: it keeps the pointer, uses
    the cache in makeFactory and, under its own lock, in every later request that builds a
    boundary rule, so the cache must outlive the factory and all its copies. No two builds
    that use one cache may run at the same time; in particular two factories that share a
    cache must not be called concurrently, and a constructor call with that cache must
    not run while one of the factories is being called.

    **Complexity.** A build makes S solves instead of one per patch direction, where S is
    the number of distinct keys. A lookup costs O(log S) key comparisons of cost
    O(knot vector size) each. Memory is O(sum of the sizes of the stored shapes) and
    grows across builds until clear().

    @ingroup Assembler
*/
template<class T>
class GISMO_EXPORT gsHPatchRuleCache
{
public:

    /// Empty cache with all counters zero.
    gsHPatchRuleCache();

    /// Number of lookups answered from the cache (no solve).
    index_t numHits() const;

    /// Number of 1D solves made through the cache.
    index_t numSolves() const;

    /// Number of solves that did not converge; equals failedKeys().size().
    index_t numFailures() const;

    /// Number of stored shapes (converged and failed).
    index_t size() const;

    /// (normalised even target knot vector, min(L,1)) of the failed solves, in order of failure.
    const std::vector< std::pair<gsKnotVector<T>, T> > & failedKeys() const;

    /// Drops all shapes and the failed-key list and resets the counters.
    void clear();

private:

    friend class gsHPatchRule<T>;

    /// Rule of one shape on [0,1]; \a nodes and \a weights are empty if !converged.
    struct Entry
    {
        bool        converged;
        gsVector<T> nodes;
        gsVector<T> weights;
    };

    /// (min(L,1), (degree, repeated knots of the normalised even target knot vector))
    typedef std::pair< T, std::pair< index_t, std::vector<T> > > Key;

    std::map<Key, Entry> m_shapes;
    std::vector< std::pair<gsKnotVector<T>, T> > m_failedKeys;
    index_t m_hits, m_solves, m_failures;
};

/**
    \brief Per-atomic-cell hierarchical patch quadrature rule.

    For each leaf (atomic cell = connected macro region at a fixed level) of
    the hierarchical k-d-tree, a 1D Johannessen Newton rule is computed per
    parametric direction and the results are tensorized.  The resulting
    quadrature points are stored in per-direction sorted node→weight arrays
    (one per atomic cell) and are filtered by element bounds in \ref mapTo —
    exactly as \ref gsPatchRule stores and queries a single patch.

    Exactness is guaranteed for every atomic cell, truncated or not: any
    active (T)HB basis function restricted to a leaf D_r of level m coincides
    there with a plain level-m tensor B-spline (truncation only removes
    finer-level functions whose support lies in \f$\Omega^{m+1}\f$, which
    vanishes identically on \f$D_r\f$). Hence a rule targeting the cell's
    level-m breakpoints integrates all active-function products on that cell
    exactly, and summing over the leaf partition gives the exact global
    matrix -- \ref makeIntegrationBasis is not required for exactness. It
    remains available as an optional utility for callers that want every
    element to lie under at most one finer truncated function for other
    reasons (e.g. simplifying downstream bookkeeping).

    **Leaf merging.** With a MergeMode other than \c none, the leaves of each level are replaced
    by boxes of the same level that tile the same region (the cells of the rule), and one
    rule is built per box. A box is a union of leaves of its level and is itself a box, so
    on it every active function coincides with a level-m tensor B-spline combination and
    the exactness argument above carries over. The kd-tree of the domain is not modified.
    The kd leaves stay the unit of the lookup: see @ref _findCell.

    **Handles.** The leaf data (the sorted node-weight arrays of all cells and the index
    from kd leaves to cells) is filled once by
    the constructor and is immutable afterwards. Copying a rule (copy construction or
    assignment) gives a handle: it shares the leaf data and owns its lookup position
    (see @ref m_lastLeaf) by value. Copy construction costs O(1): one reference-count
    increment and a few scalars; it allocates nothing. Copy assignment may release the
    previous leaf data, which costs O(L) (L = number of cells) when the target was its last
    owner. The leaf data lives as long as the last handle. A handle
    must be used by one thread at a time, but handles of one rule can be used concurrently
    by different threads. Copy a rule only while no other thread calls mapTo on it, because
    the copy reads the lookup position.

    \ingroup Assembler
*/
template<class T>
class gsHPatchRule GISMO_FINAL : public gsQuadRule<T>
{
public:

    typedef memory::unique_ptr<gsHPatchRule> uPtr;

    /// Same type as gsExprAssembler<T>::QuadratureFactory.
    typedef std::function<typename gsQuadRule<T>::uPtr(const gsBasis<T> & basis,
                                                        const gsOptionList & options,
                                                        index_t patch,
                                                        short_t fixedDirection)> QuadratureFactory;

    /**
     * @brief How the constructor partitions the region of each level into the cells of the rule.
     *
     * The cells of a level are boxes (in the level-local unique-knot indices) whose union is
     * the union of the kd-tree leaves of that level. Merging acts per level, because
     * different levels have different target spaces, and per rule, hence per patch.
     * A candidate is scored with the node count of the rule itself: per direction other
     * than the fixed one, the number of nodes of the even target knot vector of the box,
     * multiplied over the directions. The score assumes that the 1D solves converge; a box
     * whose solve fails is integrated by element-wise Gauss in that direction, see
     * numFallbackCells().
     *
     * - \c none: the cells are the kd leaves, in leaf-iterator order.
     * - \c pairwise: repeated merging of boxes whose union is a box, starting from the kd
     *   leaves. Never more nodes than \c none. Complexity O(L^2 log L) in the worst case
     *   (L = leaves of the level), a few O(L log L) passes in practice. It depends on the
     *   kd leaves, hence on the refinement history of the domain.
     * - \c sweep: the canonical slab decomposition of the region, with the axis order of
     *   fewest nodes per level (ties to the lexicographically smallest order). The result
     *   depends on the region only, not on the kd-tree: it is history-independent. It can have
     *   more nodes than \c none. Complexity O(d! L^d log L) in the worst case, far less
     *   in practice.
     * - \c best: both candidates, and per level the one with fewer model nodes (a tie
     *   goes to \c sweep). Hence best <= pairwise <= none in node count, and best <= sweep.
     *   As \c pairwise, \c best is not history-independent in general; it is wherever
     *   \c sweep wins or ties on every level.
     *
     * With every mode other than \c none the cells are sorted by (level, lower corner, upper
     * corner), so equal partitions give bitwise equal rules.
     */
    enum class MergeMode { none, pairwise, sweep, best };

    /// Partition of one build; shared by all handles of the rule. Node counts are model
    /// counts: the sum over the boxes of the product of the 1D node counts of the
    /// directions other than the fixed one. When no 1D solve fails, numPointsTotal() equals
    /// the model count of the partition that was kept: nodesLeaves for MergeMode::none,
    /// nodesPairwise for \c pairwise, nodesSweep for \c sweep, and for \c best the sum over
    /// the levels of the smaller of the two candidate counts of each level. That sum is at
    /// most min(nodesPairwise, nodesSweep), and equal to it when levelsPairwise == 0.
    struct PartitionInfo
    {
        PartitionInfo()
        : numLeaves(0), numCells(0), nodesLeaves(0), nodesPairwise(-1), nodesSweep(-1),
          levelsPairwise(0), seconds(0) {}

        index_t numLeaves;      ///< kd-tree leaves of the domain
        index_t numCells;       ///< cells of the rule
        index_t nodesLeaves;    ///< model node count of the kd-tree leaves
        index_t nodesPairwise;  ///< model node count of the pairwise candidate, -1 if not computed
        index_t nodesSweep;     ///< model node count of the sweep candidate, -1 if not computed
        index_t levelsPairwise; ///< levels on which the pairwise candidate was kept
        double  seconds;        ///< wall time of candidates + choice (gsStopwatch), 0 for none
    };

    /// Default empty constructor
    gsHPatchRule()
    : m_cells(memory::make_shared(new std::vector<CellData>)),
      m_index(memory::make_shared(new LeafIndex)),
      m_deg(0), m_reg(0), m_over(false), m_fixDir(-1), m_dim(0), m_maxLevel(0),
      m_merge(MergeMode::none), m_lastLeaf(0) {}

    /**
     * @brief Construct a per-cell hierarchical patch rule (a cell is a kd leaf, or with merging a box of leaves of one level).
     *
     * @param domain       A hierarchical domain (must be a \ref gsHDomain, as
     *                     returned by gsHTensorBasis::domain()).
     * @param degree       Target integration degree, by convention 2*p where p
     *                     is the underlying basis degree (e.g. 2*p for mass).
     * @param regularity   Target regularity at a *simple* interior knot (e.g.
     *                     p-1 for max-continuity mass, p-2 for stiffness). At
     *                     an interior breakpoint of multiplicity mu_j > 1 the
     *                     effective target regularity used is offset to
     *                     regularity + 1 - mu_j (floored at -1).
     * @param overintegrate Reserved (must be false).
     * @param fixDir       Direction to fix (-1 = none).
     * @param cache        Optional, caller-owned store of solved 1D shapes; not stored by
     *                     the rule, so it only has to outlive this call. \c nullptr uses
     *                     a cache local to this build. Must not be shared by builds that
     *                     run concurrently (see @ref gsHPatchRuleCache).
     * @param merge        How the leaves of each level are merged into cells (see
     *                     @ref MergeMode for the guarantees, the history dependence and
     *                     the complexity of each mode). The default \c none keeps the
     *                     kd leaves as cells. The partition is built before the cells,
     *                     and its wall time is reported by partitionInfo().
     */
    gsHPatchRule(const gsDomain<T> & domain,
                 index_t  degree,
                 index_t  regularity,
                 bool     overintegrate,
                 short_t  fixDir = -1,
                 gsHPatchRuleCache<T> * cache = nullptr,
                 MergeMode merge = MergeMode::none);

    /// Same arguments as the constructor; returns the rule in a unique pointer.
    static uPtr make(const gsDomain<T> & domain,
                     index_t  degree,
                     index_t  regularity,
                     bool     overintegrate,
                     short_t  fixDir = -1,
                     gsHPatchRuleCache<T> * cache = nullptr,
                     MergeMode merge = MergeMode::none)
    { return uPtr(new gsHPatchRule(domain,degree,regularity,overintegrate,fixDir,cache,merge)); }

    /**
     * @brief Quadrature factory that hands out handles of rules built once.
     *
     * Returns a callable of the type of gsExprAssembler<T>::QuadratureFactory, for
     * gsExprAssembler::setQuadratureFactory. It owns one volume rule per patch of
     * \a integration and one boundary rule per (patch, direction); a call returns a handle
     * (see the class documentation) of the requested rule, so the Newton solves of a
     * build are paid once and not once per thread and per assembly.
     *
     * A call with (\a basis, \a options, \a patch, \a fixedDirection) returns
     * - for \a fixedDirection == -1 a handle to the volume rule of \a patch,
     * - for \a fixedDirection >= 0 a handle to the boundary rule of
     *   (\a patch, \a fixedDirection), which serves both sides of that direction.
     *
     * It ignores \a basis and \a options (apart from a dimension assertion): the rule
     * parameters are fixed at creation, and the rule is selected by \a patch and
     * \a fixedDirection only. A \a patch outside the patches of \a integration, or a
     * \a fixedDirection outside [-1, dim), raises an error.
     *
     * The volume rules are built eagerly and serially in makeFactory (a patch that is not
     * hierarchical makes the constructor throw there). A boundary rule is built on the
     * first request of its (patch, direction), under a lock, and kept. Volume requests
     * take no lock.
     *
     * **Thread safety.** The callable is safe to call concurrently, and copies of it share
     * one state. Each returned handle belongs to one thread at a time. The lock is
     * an OpenMP lock (the stub of gsOpenMP when OpenMP is off) taken through a scope guard,
     * so an error raised inside a lazy build does not leave it set.
     *
     * **Lifetime.** The factory clones \a integration, so the argument only has to live
     * during this call; all rules, volume and boundary, are built on the clone. Refining the
     * caller's multi-basis afterwards does not update the factory: build a new one. A
     * non-null \a cache is kept as a pointer and used for every build of the factory, so it
     * must outlive the factory and all its copies; no two builds that use it may run
     * at the same time (see @ref gsHPatchRuleCache). With \c nullptr the factory owns a
     * cache of its own. Handles do not depend on the factory or on the cache: they stay
     * valid after every copy of the factory is destroyed.
     *
     * **Complexity.** Creation makes one rule build per patch. A volume call is O(1) and
     * allocates one object. A boundary call is O(1) plus the lock, and one rule build on
     * the first request of a slot. A handle shares O(number of nodes) data and owns O(1).
     *
     * @param integration   Multi-basis of hierarchical bases; patch \a i is built from
     *                      its basis \a i.
     * @param degree        Target integration degree (see the constructor).
     * @param regularity    Target regularity (see the constructor).
     * @param overintegrate Reserved (must be false), as in the constructor.
     * @param cache         Optional, caller-owned store of solved 1D shapes (see above).
     * @param merge         Leaf merging of every volume and boundary rule of the factory
     *                      (see @ref MergeMode). The default \c none keeps the kd leaves.
     */
    static QuadratureFactory makeFactory(const gsMultiBasis<T> & integration,
                                         index_t  degree,
                                         index_t  regularity,
                                         bool     overintegrate,
                                         gsHPatchRuleCache<T> * cache = nullptr,
                                         MergeMode merge = MergeMode::none);

    // User-declared, so that no implicit move constructor or move assignment exists: a move
    // would leave the shared pointers m_cells and m_index of the source null, and every
    // handle must keep both non-null.
    ~gsHPatchRule() {}

public:

    /**
     * @brief Map quadrature rule to the element \a [lower, upper].
     *
     * Looks up the cell (a kd leaf, or a merged box) whose parametric box contains
     * \a [lower, upper], filters the cell's per-direction node arrays to the element
     * bounds, and tensorizes the result (direction 0 fastest). For a merged cell the
     * elements of the cell share its nodes: the nodes returned for an element are those
     * of the cell inside the element, and the elements of the cell together get every
     * node once.
     *
     * An element outside every cell is a precondition violation and raises an
     * error.
     *
     * Complexity: the cell lookup of @ref _findCell, plus O(d log M) for the
     * range search (M = nodes per direction of the cell) plus O(d * nPts) to
     * write the tensor product. No heap allocation except the resize of
     * \a nodes and \a weights when the point count changes.
     */
    void mapTo(const gsVector<T>& lower, const gsVector<T>& upper,
               gsMatrix<T> & nodes, gsVector<T> & weights) const;

    /// 1D overload (requires dim() == 1); same behaviour and complexity as
    /// the d-dimensional overload.
    void mapTo(T startVal, T endVal,
               gsMatrix<T> & nodes, gsVector<T> & weights) const;

    void mapToAll(const std::vector<T> &,
                  gsMatrix<T> &, gsVector<T> &) const
    { GISMO_NO_IMPLEMENTATION }

    index_t dim() const { return m_dim; }

    /// Total number of quadrature points across all cells, merged or not.
    index_t numPointsTotal() const;

    /**
     * @brief The partition of the build: leaves, cells, model node counts of the candidates.
     *
     * Shared by all handles of the rule, valid as long as one handle lives. The node
     * counts are model counts that assume convergent 1D solves. With numFallbackCells() == 0,
     * numPointsTotal() equals the count of the partition that was kept (see PartitionInfo);
     * for \c best that is the per-level minimum summed over the levels, at most
     * min(nodesPairwise, nodesSweep). Guarantees (see @ref MergeMode):
     * best <= pairwise <= leaves in model count, and best <= sweep; \c sweep is
     * history-independent while \c pairwise and \c best are not in general. The
     * member \c seconds is the wall time of the candidates and the choice and is
     * 0 for MergeMode::none. Complexity O(1).
     */
    const PartitionInfo & partitionInfo() const;

    /**
     * @brief Number of cells (merged or not) with at least one direction (other than the
     *        fixed one) in element-wise Gauss because its 1D solve failed.
     *
     * A direction with equally many Johannessen and Gauss nodes uses Gauss but is
     * not counted: it costs no node and loses nothing. The Johannessen rule never has
     * more nodes than Gauss: a block of e elements between two knots of multiplicity
     * deg+1 has at most ceil((e*deg+1)/2) nodes, which is at most (deg/2+1)*e. Complexity O(number of cells).
     */
    index_t numFallbackCells() const;

    /**
     * @brief Refine \a basis so every element lies under at most one finer
     *        truncated function.
     *
     * Optional utility: exactness of \ref gsHPatchRule does not depend on
     * this refinement (see the class documentation), but some callers may
     * still want the one-truncated-function-per-element property for other
     * purposes.
     */
    template<short_t d>
    static gsTHBSplineBasis<d,T,true>
    makeIntegrationBasis(const gsTHBSplineBasis<d,T,true> & basis);

    /// True if \a domain is a hierarchical domain (\ref gsHDomain) of dim 1..4.
    static bool isHierarchical(const gsDomain<T> & domain);

    /**
     * @brief Returns the finest-level tensor domain of a hierarchical \a domain.
     *
     * Used by the \ref gsQuadrature factory to let \ref gsPatchRule operate on
     * a hierarchical basis by defaulting to the finest level.  Returns a null
     * pointer if \a domain is not hierarchical.
     */
    static typename gsDomain<T>::Ptr finestTensorDomain(const gsDomain<T> & domain);

private:

    /// Per-cell quadrature data (sorted node→weight arrays per direction); a cell is a
    /// kd leaf or a merged box of leaves of one level.
    struct CellData
    {
        gsVector<T>                lowerParam;  ///< Parametric lower corner
        gsVector<T>                upperParam;  ///< Parametric upper corner
        gsVector<index_t>          lowerIndex;  ///< Unique-knot index lower corner (at cell level)
        gsVector<index_t>          upperIndex;  ///< Unique-knot index upper corner (at cell level)
        index_t                    level;       ///< Hierarchical level
        /// dirMaps[dd]: ascending-by-node list of (node,weight) in direction dd
        std::vector<std::vector<std::pair<T,T>>> dirMaps;
        /// True if a non-fixed direction is in element-wise Gauss because its solve failed
        /// (a tie is not counted)
        bool                       usedFallback;
    };

    /// Shared state of a factory returned by makeFactory.
    struct FactoryState;

    /// Box of elements of one level, in level-local unique-knot indices: the elements
    /// with lo[dd] <= i_dd < up[dd] in every direction dd.
    template<short_t d>
    struct IndexBox
    {
        index_t lo[d];
        index_t up[d];
    };

    /// Lookup structure shared by all handles; immutable after construction.
    /// The kd leaves in leaf-iterator order (the order of volume iteration) with their
    /// parametric boxes, and for each leaf the cells that meet it in a set of positive volume.
    struct LeafIndex
    {
        std::vector<T>       lo, up;  ///< Parametric kd-leaf boxes, d entries per leaf
        std::vector<index_t> first;   ///< CSR offsets into \a cell, numLeaves + 1 entries
        std::vector<index_t> cell;    ///< Cells that meet leaf i: cell[first[i]] .. cell[first[i+1]-1]
        PartitionInfo        info;
    };

    // ----------------------------------------------------------------
    // Private helpers
    // ----------------------------------------------------------------

    /// Dimension-specific construction (downcasts \a domain to gsHDomain<d,...>);
    /// the 1D rules of the patch directions come from \a cache.
    template<short_t d>
    void _build(const gsDomain<T> & domain, gsHPatchRuleCache<T> & cache);

    /**
     * @brief Even target knot vector of a cell in one direction, and its node count.
     *
     * The target knot vector has multiplicity m_deg+1 at both ends of the unique-knot index
     * range [\a lo, \a up] of \a kvLev and, at an interior unique knot of multiplicity
     * mu_j in \a kvLev, multiplicity min(m_deg - m_reg - 1 + mu_j, m_deg + 1) (target
     * regularity m_reg + 1 - mu_j, floored at -1). It is made even block by block by
     * @ref _ensureEven1D. The 1D rule has \a n1d = (dimension)/2 nodes. The build of the
     * cells and the scoring of the candidate partitions both use this function, so the
     * model count and the rule cannot differ.
     *
     * Complexity: O((up - lo) * m_deg) for the knots, plus _ensureEven1D.
     */
    gsKnotVector<T> _evenTargetKnots1D(const gsKnotVector<T> & kvLev,
                                       index_t lo, index_t up, index_t & n1d) const;

    /// Model node count of \a boxes of level \a level: the sum over the boxes of the
    /// product of the 1D node counts of the non-fixed directions. Never touches a
    /// gsHPatchRuleCache. Complexity O(sum over boxes of the extents of the box).
    template<short_t d>
    index_t _modelNodes(const gsHTensorBasis<d,T> & basis, index_t level,
                        const std::vector<IndexBox<d> > & boxes) const;

    /// Strict weak order of boxes: lower corner lexicographically, then upper corner.
    template<short_t d>
    static bool _boxLess(const IndexBox<d> & x, const IndexBox<d> & y);

    /**
     * @brief Repeatedly merges boxes whose union is a box, in place.
     *
     * Two boxes merge when they agree in every axis but one and touch along it. Passes
     * over the axes repeat until nothing merges. The result depends on the set of
     * boxes only, and a merge never raises the node count, since the merged target space
     * has at most the dimension of the sum of the two (the interior multiplicity at the
     * shared knot is added and one end block of m_deg+1 is dropped).
     *
     * Complexity: a pass over one axis sorts O(L log L) and merges runs in one scan; every
     * productive pass removes a box, so O(L^2 log L) in the worst case and a few passes in practice.
     */
    template<short_t d>
    static void _mergePairwise(std::vector<IndexBox<d> > & boxes);

    /**
     * @brief Canonical slab decomposition of the union of the disjoint boxes \a B for the
     *        axis order \a sigma, from depth \a k on.
     *
     * At depth k the cut set holds every face coordinate of \a B along a = sigma[k]. Each
     * slab between consecutive cuts is decomposed by the recursive call (axes sigma[k+1..])
     * on the boxes spanning the slab, clipped to it; a piece is glued to the piece of the
     * previous slab that has the same extents in the remaining axes and touches it. The
     * last axis merges touching intervals directly.
     *
     * Why this is canonical: the cut set contains every face coordinate of the kd leaves,
     * including interior faces that depend on the refinement history. At a cut where the
     * cross-section of the region is the same on both sides, the recursive results on both
     * sides are identical, so every piece continues across the cut and the cut leaves no
     * trace. Every other cut is a face of the region itself. The output therefore depends on
     * the region and \a sigma, not on the kd-tree. Merging happens inner axes first; a second
     * pass would find nothing, because pieces touching along sigma[k] with equal other
     * extents are already glued inside the slabs.
     *
     * Complexity: a depth has at most 2L-1 slabs for L boxes and the pieces before gluing
     * number at most (2L)^(d-1) L = O(L^d), so the result can contain more boxes than L (a
     * T-shaped region made of 2 leaves gives 3 pieces for one axis order). The slab filter costs O(L) per slab
     * and the ordered maps add a log factor: O(L^d log L) per axis order. Ordered containers
     * only, so the result does not depend on hash order.
     */
    template<short_t d>
    static std::vector<IndexBox<d> > _sweepRec(const std::vector<IndexBox<d> > & B,
                                               const std::array<short_t,d> & sigma, short_t k);

    /// Sweep candidate of one level: runs _sweepRec for every axis order (lexicographic
    /// order from 0,..,d-1) and keeps the order with fewest model nodes, ties to the
    /// lexicographically smallest order. The pieces are sorted by _boxLess; \a nodes is
    /// their model count. Complexity: d! times _sweepRec and _modelNodes.
    template<short_t d>
    std::vector<IndexBox<d> > _sweepLevel(const gsHTensorBasis<d,T> & basis, index_t level,
                                          const std::vector<IndexBox<d> > & leaves,
                                          index_t & nodes) const;

    /**
     * @brief Cells of every kd leaf: the CSR lists of LeafIndex.
     *
     * Leaf i is listed with every cell of its level whose box meets the leaf in a set of
     * positive volume, in ascending cell order. The cells of a level are indexed by a
     * bounding-volume tree, so a leaf is queried in O(log C + f) for C cells of its level
     * and f cells that meet it. Total O(C log C + L log C + sum of f). In debug builds
     * it checks that the intersections with the listed cells cover every leaf exactly
     * and that the cells of every level cover the leaves of that level exactly.
     */
    template<short_t d>
    static void _leafCells(const std::vector<IndexBox<d> > & leafBox,
                           const std::vector<index_t> & leafLev,
                           const std::vector<IndexBox<d> > & cellBox,
                           const std::vector<index_t> & cellLev,
                           std::vector<index_t> & first, std::vector<index_t> & cell);

    /**
     * @brief Chooses the cells of the rule from the kd leaves according to m_merge.
     *
     * For MergeMode::none the cells are the leaves in leaf-iterator order. Otherwise the
     * leaves of each level are replaced by the pairwise and/or the sweep candidate (per
     * m_merge; for \c best the one with fewer model nodes, a tie to the sweep), and the cells
     * are sorted by (level, _boxLess). Fills the candidate fields of \a info.
     *
     * Complexity: see MergeMode.
     */
    template<short_t d>
    void _partition(const gsHTensorBasis<d,T> & basis,
                    const std::vector<IndexBox<d> > & leafBox,
                    const std::vector<index_t> & leafLev,
                    std::vector<IndexBox<d> > & cellBox,
                    std::vector<index_t> & cellLev,
                    PartitionInfo & info) const;

    /**
     * @brief Rule of a normalised even target space for a leaf direction of clamped length.
     *
     * Looks \a refKv (the even target knot vector mapped to [0,1]) and \a Lc =
     * min(L,1) up in \a cache; on a miss solves it with the tolerance
     * \f$10^{-14}/L_c\f$ and stores the result, converged or failed. The tolerance is
     * the physical stop \f$10^{-14}\f$ for a leaf direction of length \f$L \le 1\f$ (a
     * rule scaled by \f$L\f$ has residual \f$L\f$ times that on [0,1]) and the stop of a
     * length-1 leaf for \f$L > 1\f$. Updates the counters of \a cache.
     *
     * Complexity: O(log S) key comparisons for the lookup, plus one _compute1D on a miss.
     *
     * @return the stored entry; valid until \a cache is cleared or destroyed
     */
    static const typename gsHPatchRuleCache<T>::Entry &
    _referenceRule(const gsKnotVector<T> & refKv, T Lc, gsHPatchRuleCache<T> & cache);

    /**
     * @brief Ensure that the basis dimension of \a kv is even, block by block.
     *
     * A knot of multiplicity \f$q+1\f$ makes the space discontinuous there and splits it
     * into independent blocks, one between each pair of consecutive such knots (the end
     * knots included). A block with \f$u\f$ interior knots of multiplicities \f$m_j\f$ has
     * dimension \f$q+1+\sum_j m_j\f$. A node serves one block only, so a rule with
     * \f$\text{dim}/2\f$ nodes exists only if every block has an even dimension. One simple
     * knot is inserted in the middle of the longest span of every block of odd dimension. If
     * \a kv has no interior knot of multiplicity \f$q+1\f$ it is a single block, and one knot
     * is inserted if its dimension is odd.
     *
     * @param kv knot vector of the target space
     * @return \a kv with the knots inserted
     */
    static gsKnotVector<T> _ensureEven1D(const gsKnotVector<T> & kv);

    /**
     * @brief Compute Greville abscissae and exact B-spline integrals for \a kv.
     *
     * Uses the exact formula: integral N_{i,q} = (kv[i+q+1] - kv[i]) / (q+1).
     */
    static std::pair<gsMatrix<T>, gsVector<T>> _integrate1D(const gsKnotVector<T> & kv);

    /**
     * @brief 1D Johannessen Newton solver. Requires even basis dimension.
     *
     * Let \f$N_i\f$, \f$i=0,\dots,\text{size}-1\f$, be the B-splines of degree
     * \f$q\f$ on \a evenKv, \f$\text{size}=2n\f$, and \f$I_i=\int N_i\f$. The
     * \f$n\f$ weights \f$w_k\f$ and nodes \f$\xi_k\f$ solve the moment equations
     * \f[ F_i(w,\xi) = \sum_{k=0}^{n-1} w_k N_i(\xi_k) - I_i = 0,
     *     \qquad i = 0,\dots,2n-1, \f]
     * by Newton's method started from the midpoints of consecutive Greville
     * abscissae. If the Newton iteration fails (a node leaves the open interval of
     * \a evenKv, a pivot is zero or not finite, or \a itMax iterations are used up), the
     * solve restarts from the same starting guess with a Levenberg-Marquardt iteration
     * on \f$\tfrac12\|F\|^2\f$. Each step solves
     * \f$(J^T J + \mu\,\mathrm{diag}(J^T J))\,d = -J^T F\f$ with the damping \f$\mu\f$
     * lowered by 3 after an accepted step and raised by 4 after a rejected one. A step is
     * accepted only if all nodes stay in the open interval and \f$\|F\|_2\f$ decreases. The
     * full Newton step is large when \f$J\f$ is ill-conditioned (at high degree and low
     * regularity), and this damping bounds it. A converged Newton iteration is not
     * changed by the restart. The Jacobian has the entries
     * \f$\partial F_i/\partial w_k = N_i(\xi_k)\f$ and
     * \f$\partial F_i/\partial \xi_k = w_k N_i'(\xi_k)\f$.
     *
     * The unknowns are ordered interleaved, \f$(w_0,\xi_0,w_1,\xi_1,\dots)\f$.
     * The pair of columns \f$(2k,2k+1)\f$ is nonzero only in the \f$q+1\f$ rows
     * \f$a_k,\dots,a_k+q\f$ of the B-splines active at \f$\xi_k\f$. At the
     * starting guess \f$\xi_k\f$ lies between the Greville abscissae
     * \f$g_{2k}\f$ and \f$g_{2k+1}\f$ of gsKnotVector::greville_into, which
     * shifts the first of two equal abscissae (at a knot of multiplicity
     * \f$q+1\f$) down by \f$10^{-9}\f$; hence \f$2k-q \le a_k \le 2k+1\f$, and
     * the Jacobian is banded with lower bandwidth \f$kl \le q+1\f$ and upper
     * bandwidth \f$ku \le q+1\f$ (\f$ku = q+1\f$ when such a knot starts at
     * index \f$2k+1\f$; the Jacobian is then singular). Nodes may move
     * during the iteration, so \f$kl\f$ and \f$ku\f$ are recomputed from the
     * active indices of the current iterate in every iteration. The Jacobian
     * is stored in band form and solved by _bandLUSolve().
     *
     * \f$J^T J\f$ is banded, since the columns of two nodes couple only if their active
     * B-spline ranges overlap, and is stored in band form. Its half-bandwidth is
     * \f$hb = O(q)\f$ while the nodes stay in ascending order; if an iterate reorders
     * nodes, \f$hb\f$ can grow up to O(n).
     *
     * Complexity: O(size q^2) time and O(size q) memory per Newton iteration. Per
     * Levenberg-Marquardt iteration (the path of every target whose Newton iteration
     * fails) O(n q^2) to find the coupled node pairs and fill \f$J^T J\f$ when the nodes
     * are in ascending order, O(n^2) otherwise, plus O(size hb^2) per damping trial.
     *
     * @param[in]     evenKv  knot vector of the target space (even dimension)
     * @param[out]    nodes   the n nodes
     * @param[out]    weights the n weights
     * @param[in]     tol     absolute convergence tolerance on the 2-norm of the
     *                        residual \f$F\f$ of the system it is given (not scaled by the
     *                        norm of the integrals). The constructor passes the target space
     *                        mapped to [0,1] with \f$\text{tol} = 10^{-14}/\min(L,1)\f$:
     *                        the physical stop \f$10^{-14}\f$ for a leaf direction of length
     *                        \f$L \le 1\f$, and the stop of a length-1 leaf for \f$L > 1\f$.
     * @param[in]     itMax   maximum number of iterations of each of the two phases
     * @return true on convergence of either phase, false if both fail (the
     *         Levenberg-Marquardt phase finds no step that lowers the residual, produces
     *         a non-finite residual, or \a itMax is reached). A non-finite residual or
     *         Newton update ends the Newton phase at once; the nodes and weights of a
     *         call that returns false are partial iterates and must not be used.
     */
    static bool _compute1D(const gsKnotVector<T> & evenKv,
                            gsVector<T> & nodes,
                            gsVector<T> & weights,
                            T       tol   = T(1e-14),
                            index_t itMax = 200);

    /**
     * @brief Solves A x = b in place by banded Gaussian elimination with row partial
     *        pivoting (the algorithm of LAPACK DGBTF2/DGBTRS).
     *
     * \a AB holds A in LAPACK band storage: (2*kl+ku+1) x N, with A(i,j) at
     * AB(kl+ku+i-j, j) for max(0,j-ku) <= i <= min(N-1,j+kl); the top kl rows are
     * workspace for the fill-in of the pivoting and must be zero on entry. On return
     * \a AB holds the factors and \a b the solution.
     *
     * Complexity: O(N kl (kl+ku)) time, (2kl+ku+1) N memory.
     *
     * @param[in,out] AB  band storage of A, overwritten by the LU factors
     * @param[in]     kl  number of sub-diagonals of A
     * @param[in]     ku  number of super-diagonals of A
     * @param[in,out] b   right-hand side on entry, solution on return
     * @return false on a zero or non-finite pivot (A singular or not finite); \a b is then undefined.
     */
    static bool _bandLUSolve(gsMatrix<T> & AB, index_t kl, index_t ku, gsVector<T> & b);

    /**
     * @brief Fill \a dirMap with element-wise Gauss nodes/weights for one
     *        direction over the unique-knot index range [lowerIdx, upperIdx).
     */
    static void _fillGaussFallback1D(std::vector<std::pair<T,T>> & dirMap,
                                     const gsKnotVector<T> & kv,
                                     index_t lowerIdx,
                                     index_t upperIdx,
                                     index_t nq);

    /// Sort \a dirMap ascending by node (and assert uniqueness in debug).
    static void _finalizeDirMap(std::vector<std::pair<T,T>> & dirMap);

    /**
     * @brief Map the quadrature rule to the element \a [lower, upper]
     *        (arrays of length dim()); the implementation of both @ref mapTo
     *        overloads.
     *
     * Complexity: see @ref mapTo. Allocation-free apart from the resize of
     * \a nodes and \a weights; per-direction ranges live in fixed-size
     * scratch (dim() <= 4).
     */
    void _mapTo(const T * lower, const T * upper,
                gsMatrix<T> & nodes, gsVector<T> & weights) const;

    /**
     * @brief Index of the cell whose parametric box contains the
     *        element \a [lower, upper] (arrays of length dim()), or -1.
     *
     * Containment: cell contains [lower,upper] iff
     * cell.lower[dd] <= lower[dd] and upper[dd] <= cell.upper[dd] for all dd
     * (exact comparison). Two levels: probes the kd leaves, starting at the last hit
     * \a m_lastLeaf, then its successor, then cyclically; the leaf found becomes the new
     * \a m_lastLeaf, and the first of the cells that meet it that contains the element is
     * returned. If no leaf contains the element (an element of a coarser or foreign
     * domain), or none of its cells does, all cells are scanned in order.
     * Handles degenerate (boundary) elements where lower[dd] == upper[dd].
     * Never throws and never allocates; a containing cell is found whenever one exists,
     * in any iteration order.
     *
     * Complexity, with L_kd = number of kd leaves, L = number of cells and f_i = number of
     * cells that meet leaf i (1 when the leaf lies in one cell):
     * - O(1) amortised plus O(f_i d) per element in volume iteration order
     *   (gsHDomainIterator visits the leaves in the order of LeafIndex);
     * - O(nBdrElem + L_kd) per sweep for boundary or interface iteration with a
     *   fixed direction (the boundary iterator skips leaves not on the side);
     * - one O(L_kd) scan at the start of each OpenMP chunk;
     * - O(L_kd + L) per miss when elements arrive in an order unrelated to the leaf
     *   order (shuffled, or a domain other than the rule's own);
     * - returns -1 if no cell contains the element.
     */
    index_t _findCell(const T * lower, const T * upper) const;

    /// True if the element [lower, upper] (arrays of length \a dim) is
    /// contained in \a cell. O(dim).
    static bool _contains(const CellData & cell,
                          const T * lower, const T * upper, index_t dim);

    /// Raise the error for an element that no cell contains. Cold path of
    /// @ref mapTo; kept out of line so the hit path stays small.
#if defined(__GNUC__)
    [[noreturn]] __attribute__((noinline))
#elif defined(_MSC_VER)
    [[noreturn]] __declspec(noinline)
#else
    [[noreturn]]
#endif
    static void _missError(const T * lower, const T * upper, index_t dim);

    // ----------------------------------------------------------------
    // Data
    // ----------------------------------------------------------------
    /// One entry per cell (a kd leaf, or a merged box). Never null; filled once by the
    /// constructor and immutable afterwards, shared by all handles (copies) of the rule.
    memory::shared_ptr<const std::vector<CellData> > m_cells;

    /// The kd leaves, the cells that meet each, and the partition info. Never null;
    /// filled once by the constructor and immutable afterwards, shared by all handles.
    memory::shared_ptr<const LeafIndex> m_index;

    index_t m_deg;
    index_t m_reg;
    bool    m_over;
    short_t m_fixDir;
    index_t m_dim;
    index_t m_maxLevel;
    MergeMode m_merge;

    /// Position of the last-hit kd leaf of *this handle* in the leaf list of m_index, in
    /// [0, L_kd) when L_kd > 0. Mutable lookup cache. Handles share the leaf data and never
    /// share this field, so a handle needs no synchronisation as long as one thread uses it.
    mutable index_t m_lastLeaf;

}; // class gsHPatchRule

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHPatchRule.hpp)
#endif
