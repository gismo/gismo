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

#include <functional>
#include <map>
#include <utility>
#include <vector>

namespace gismo
{

template<short_t d, class T, bool Trunc> class gsTHBSplineBasis;

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

    **Handles.** The leaf data (the sorted node-weight arrays of all leaves) is filled once by
    the constructor and is immutable afterwards. Copying a rule (copy construction or
    assignment) gives a handle: it shares the leaf data and owns its lookup position
    (see @ref m_lastCell) by value. Copy construction costs O(1): one reference-count
    increment and a few scalars; it allocates nothing. Copy assignment may release the
    previous leaf data, which costs O(L) (L = number of leaves) when the target was its last
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

    /// Default empty constructor
    gsHPatchRule()
    : m_cells(memory::make_shared(new std::vector<CellData>)),
      m_deg(0), m_reg(0), m_over(false), m_fixDir(-1), m_dim(0), m_maxLevel(0),
      m_lastCell(0) {}

    /**
     * @brief Construct a per-atomic-cell hierarchical patch rule.
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
     */
    gsHPatchRule(const gsDomain<T> & domain,
                 index_t  degree,
                 index_t  regularity,
                 bool     overintegrate,
                 short_t  fixDir = -1,
                 gsHPatchRuleCache<T> * cache = nullptr);

    /// Same arguments as the constructor; returns the rule in a unique pointer.
    static uPtr make(const gsDomain<T> & domain,
                     index_t  degree,
                     index_t  regularity,
                     bool     overintegrate,
                     short_t  fixDir = -1,
                     gsHPatchRuleCache<T> * cache = nullptr)
    { return uPtr(new gsHPatchRule(domain,degree,regularity,overintegrate,fixDir,cache)); }

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
     */
    static QuadratureFactory makeFactory(const gsMultiBasis<T> & integration,
                                         index_t  degree,
                                         index_t  regularity,
                                         bool     overintegrate,
                                         gsHPatchRuleCache<T> * cache = nullptr);

    ~gsHPatchRule() {}

public:

    /**
     * @brief Map quadrature rule to the element \a [lower, upper].
     *
     * Looks up the atomic cell whose parametric box contains \a [lower, upper],
     * filters the cell's per-direction node arrays to the element bounds, and
     * tensorizes the result (direction 0 fastest).
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

    /// Total number of quadrature points across all atomic cells.
    index_t numPointsTotal() const;

    /**
     * @brief Number of atomic cells with at least one direction (other than the
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

    /// Per-atomic-cell quadrature data (sorted node→weight arrays per direction).
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

    // ----------------------------------------------------------------
    // Private helpers
    // ----------------------------------------------------------------

    /// Dimension-specific construction (downcasts \a domain to gsHDomain<d,...>);
    /// the 1D rules of the patch directions come from \a cache.
    template<short_t d>
    void _build(const gsDomain<T> & domain, gsHPatchRuleCache<T> & cache);

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
     * @brief Index of the atomic cell whose parametric box contains the
     *        element \a [lower, upper] (arrays of length dim()), or -1.
     *
     * Containment: cell contains [lower,upper] iff
     * cell.lower[dd] <= lower[dd] and upper[dd] <= cell.upper[dd] for all dd
     * (exact comparison). Probes the last hit \a m_lastCell, then its
     * successor, then scans the remaining cells cyclically starting at
     * \a m_lastCell + 2; the cell found becomes the new \a m_lastCell.
     * Handles degenerate (boundary) elements where lower[dd] == upper[dd].
     * Never throws; a containing cell is found whenever one exists, in any
     * iteration order.
     *
     * Complexity, with L = number of cells:
     * - O(1) amortised per element in volume iteration order
     *   (gsHDomainIterator visits the leaves in the order of m_cells);
     * - O(nBdrElem + L) per sweep for boundary or interface iteration with a
     *   fixed direction (the boundary iterator skips leaves not on the side);
     * - one O(L) scan at the start of each OpenMP chunk;
     * - O(L) per miss when elements arrive in an order unrelated to the leaf
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
    /// One entry per leaf of the HTree. Never null; filled once by the constructor and
    /// immutable afterwards, shared by all handles (copies) of the rule.
    memory::shared_ptr<const std::vector<CellData> > m_cells;

    index_t m_deg;
    index_t m_reg;
    bool    m_over;
    short_t m_fixDir;
    index_t m_dim;
    index_t m_maxLevel;

    /// Index of the last-hit cell of *this handle*, in [0, L) when L = m_cells->size() > 0.
    /// Mutable lookup cache. Handles share the leaf data and never share this field, so a
    /// handle needs no synchronisation as long as one thread uses it.
    mutable index_t m_lastCell;

}; // class gsHPatchRule

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHPatchRule.hpp)
#endif
