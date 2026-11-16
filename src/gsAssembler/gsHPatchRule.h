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

    A Gauss fallback is used per cell whenever:
      - the candidate Johannessen point count >= element-wise Gauss count, or
      - the Newton solver fails to converge.

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

namespace gismo
{

template<short_t d, class T, bool Trunc> class gsTHBSplineBasis;

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

    \ingroup Assembler
*/
template<class T>
class gsHPatchRule GISMO_FINAL : public gsQuadRule<T>
{
public:

    typedef memory::unique_ptr<gsHPatchRule> uPtr;

    /// Default empty constructor
    gsHPatchRule()
    : m_deg(0), m_reg(0), m_over(false), m_fixDir(-1), m_dim(0), m_maxLevel(0),
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
     */
    gsHPatchRule(const gsDomain<T> & domain,
                 index_t  degree,
                 index_t  regularity,
                 bool     overintegrate,
                 short_t  fixDir = -1);

    static uPtr make(const gsDomain<T> & domain,
                     index_t  degree,
                     index_t  regularity,
                     bool     overintegrate,
                     short_t  fixDir = -1)
    { return uPtr(new gsHPatchRule(domain,degree,regularity,overintegrate,fixDir)); }

    ~gsHPatchRule() {}

public:

    /**
     * @brief Map quadrature rule to the element \a [lower, upper].
     *
     * Looks up the atomic cell containing the element centre, filters the
     * cell's per-direction node arrays to the element bounds, and tensorizes
     * the result.
     */
    void mapTo(const gsVector<T>& lower, const gsVector<T>& upper,
               gsMatrix<T> & nodes, gsVector<T> & weights) const;

    /// 1D overload (requires dim() == 1).
    void mapTo(T startVal, T endVal,
               gsMatrix<T> & nodes, gsVector<T> & weights) const;

    void mapToAll(const std::vector<T> &,
                  gsMatrix<T> &, gsVector<T> &) const
    { GISMO_NO_IMPLEMENTATION }

    index_t dim() const { return m_dim; }

    /// Total number of quadrature points across all atomic cells.
    index_t numPointsTotal() const;

    /// Number of atomic cells that fell back to element-wise Gauss.
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
        bool                       usedFallback;///< True if Gauss fallback was applied
    };

    // ----------------------------------------------------------------
    // Private helpers
    // ----------------------------------------------------------------

    /// Dimension-specific construction (downcasts \a domain to gsHDomain<d,...>).
    template<short_t d>
    void _build(const gsDomain<T> & domain);

    /**
     * @brief Ensure the basis dimension of \a kv is even by inserting one knot.
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
     * @return true on convergence, false on failure.
     */
    static bool _compute1D(const gsKnotVector<T> & evenKv,
                            gsVector<T> & nodes,
                            gsVector<T> & weights,
                            T       tol   = T(1e-14),
                            index_t itMax = 200);

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
     * @brief Find the atomic cell whose parametric box contains the element
     *        \a [lower, upper].
     *
     * Uses element containment: cell contains [lower,upper] iff
     * cell.lower[dd] <= lower[dd] and upper[dd] <= cell.upper[dd] for all dd.
     * Uses a last-hit cache (the standard element iteration order is
     * cell-coherent, so this is amortized O(1)); falls back to a linear scan.
     * Handles degenerate (boundary) elements where lower[dd] == upper[dd].
     */
    const CellData & _findCell(const gsVector<T> & lower,
                               const gsVector<T> & upper) const;

    /// True if the element [lower, upper] is contained in \a cell.
    static bool _contains(const CellData & cell,
                           const gsVector<T> & lower,
                           const gsVector<T> & upper);

    // ----------------------------------------------------------------
    // Data
    // ----------------------------------------------------------------
    std::vector<CellData> m_cells;    ///< One entry per leaf of the HTree

    index_t m_deg;
    index_t m_reg;
    bool    m_over;
    short_t m_fixDir;
    index_t m_dim;
    index_t m_maxLevel;

    mutable index_t m_lastCell;       ///< Last-hit cell cache (per-thread rule => safe)

}; // class gsHPatchRule

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHPatchRule.hpp)
#endif
