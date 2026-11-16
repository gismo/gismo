/** @file gsHPatchRule.hpp

    @brief Implementation of the per-atomic-cell hierarchical patch quadrature rule.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsAssembler/gsHPatchRule.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsNurbs/gsBSplineBasis.h>
#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsHSplines/gsHElement.h>
#include <gsHSplines/gsHElementHelper.h>
#include <gsDomain/gsHDomainIterator.h>
#include <gsMatrix/gsSparseSolver.h>

#include <algorithm>
#include <map>

namespace gismo
{

// ============================================================
// Constructor + dimension dispatch
// ============================================================

template<class T>
gsHPatchRule<T>::gsHPatchRule(const gsDomain<T> & domain,
                              const index_t  degree,
                              const index_t  regularity,
                              const bool     overintegrate,
                              const short_t  fixDir)
:
m_deg(degree),
m_reg(regularity),
m_over(overintegrate),
m_fixDir(fixDir),
m_dim(domain.dim()),
m_maxLevel(0),
m_lastCell(0)
{
    GISMO_ENSURE(m_reg < m_deg, "regularity must be strictly less than degree");
    GISMO_ENSURE(!m_over, "gsHPatchRule: over-integration is not implemented.");

    switch (m_dim)
    {
        case 1: _build<1>(domain); break;
        case 2: _build<2>(domain); break;
        case 3: _build<3>(domain); break;
        case 4: _build<4>(domain); break;
        default:
            GISMO_ERROR("gsHPatchRule: unsupported domain dimension " << m_dim);
    }
}

// ============================================================
// Dimension-specific construction
// ============================================================

template<class T>
template<short_t d>
void gsHPatchRule<T>::_build(const gsDomain<T> & domain)
{
    const gsHDomain<d,T,index_t> * hd =
        dynamic_cast<const gsHDomain<d,T,index_t>*>(&domain);
    GISMO_ENSURE(hd != nullptr,
                 "gsHPatchRule requires a hierarchical domain (gsHDomain).");

    const gsHTensorBasis<d,T> & basis = hd->basis();
    GISMO_ASSERT(!basis.manualLevels(),
                 "gsHPatchRule: manualLevels bases are not supported.");

    m_maxLevel = basis.maxLevel();

    auto leaf = hd->tree().beginLeafIterator();
    for (; leaf.good(); leaf.next())
    {
        const index_t l = leaf.level();

        CellData cell;
        cell.level = l;
        cell.lowerParam.resize(d);
        cell.upperParam.resize(d);
        cell.lowerIndex.resize(d);
        cell.upperIndex.resize(d);
        cell.dirMaps.resize(d);
        cell.usedFallback = false;

        // Convert index corners -> parametric corners (per direction)
        std::vector<gsKnotVector<T>> evenKvs(d);
        std::vector<index_t>         n1d(d, 0);
        index_t nCandidateTensor = 1;
        index_t nGaussTotal      = 1;

        for (short_t dd = 0; dd < d; ++dd)
        {
            const gsKnotVector<T> & kvLev =
                basis.tensorLevel(l).component(dd).knots();

            cell.lowerIndex[dd] = leaf.lowerCorner()[dd];
            cell.upperIndex[dd] = leaf.upperCorner()[dd];
            cell.lowerParam[dd] = kvLev(cell.lowerIndex[dd]);
            cell.upperParam[dd] = kvLev(cell.upperIndex[dd]);

            const index_t nElemDD = cell.upperIndex[dd] - cell.lowerIndex[dd];
            nGaussTotal *= nElemDD * (basis.degree(dd) + 1);

            if (short_t(dd) == m_fixDir)
            {
                // Fixed direction: single node at the cell midpoint, weight = span
                n1d[dd] = 1;
                continue;
            }

            // Target knot vector on the leaf: degree m_deg, and at every interior
            // breakpoint the multiplicity implied by the level-l knot vector's own
            // multiplicity mu_j.  Target regularity there is r_j - s = m_reg + 1 - mu_j
            // (Cor. 2.1), floored at -1, i.e. multiplicity capped at m_deg + 1.
            std::vector<T> repKnots;
            for (index_t i = cell.lowerIndex[dd]; i <= cell.upperIndex[dd]; ++i)
            {
                const T    z    = kvLev(i);                 // unique index -> knot value
                const bool endK = (i == cell.lowerIndex[dd] || i == cell.upperIndex[dd]);
                index_t    mult;
                if (endK)
                    mult = m_deg + 1;
                else
                {
                    const index_t muj = kvLev.multiplicity(z);
                    mult = math::min(m_deg - m_reg - 1 + muj, m_deg + 1);
                }
                repKnots.insert(repKnots.end(), mult, z);
            }
            gsKnotVector<T> targetKv(repKnots, static_cast<short_t>(m_deg));
            evenKvs[dd] = _ensureEven1D(targetKv);

            const index_t N = gsBSplineBasis<T>(evenKvs[dd]).size();
            GISMO_ENSURE(N % 2 == 0,
                "gsHPatchRule: even target basis construction failed in direction " << dd);
            n1d[dd] = N / 2;
            nCandidateTensor *= n1d[dd];
        }

        // --- Decide: Johannessen or Gauss fallback ---
        bool useFallback = (nCandidateTensor >= nGaussTotal);

        std::vector<gsVector<T>> dirNodes(d), dirWeights(d);
        if (!useFallback)
        {
            for (short_t dd = 0; dd < d && !useFallback; ++dd)
            {
                if (short_t(dd) == m_fixDir) continue;
                if (!_compute1D(evenKvs[dd], dirNodes[dd], dirWeights[dd]))
                    useFallback = true;
            }
        }

        if (useFallback)
        {
            cell.usedFallback = true;
            for (short_t dd = 0; dd < d; ++dd)
            {
                cell.dirMaps[dd].clear();
                if (short_t(dd) == m_fixDir)
                {
                    // Store both cell boundaries so degenerate boundary elements
                    // (lower[dd]==upper[dd]) can pick the correct boundary node.
                    cell.dirMaps[dd].push_back(std::make_pair(cell.lowerParam[dd], T(1)));
                    cell.dirMaps[dd].push_back(std::make_pair(cell.upperParam[dd], T(1)));
                    continue;
                }
                const gsKnotVector<T> & kvLev =
                    basis.tensorLevel(l).component(dd).knots();
                _fillGaussFallback1D(cell.dirMaps[dd], kvLev,
                                     cell.lowerIndex[dd], cell.upperIndex[dd],
                                     basis.degree(dd) + 1);
                _finalizeDirMap(cell.dirMaps[dd]);
            }
            m_cells.push_back(give(cell));
            continue;
        }

        // --- Store Johannessen results (sorted by node) ---
        for (short_t dd = 0; dd < d; ++dd)
        {
            if (short_t(dd) == m_fixDir)
            {
                // Store both cell boundaries so degenerate boundary elements
                // (lower[dd]==upper[dd]) can pick the correct boundary node.
                cell.dirMaps[dd].push_back(std::make_pair(cell.lowerParam[dd], T(1)));
                cell.dirMaps[dd].push_back(std::make_pair(cell.upperParam[dd], T(1)));
                continue;
            }
            cell.dirMaps[dd].reserve(n1d[dd]);
            for (index_t k = 0; k < n1d[dd]; ++k)
                cell.dirMaps[dd].push_back(
                    std::make_pair(dirNodes[dd][k], dirWeights[dd][k]));
            _finalizeDirMap(cell.dirMaps[dd]);
        }
        m_cells.push_back(give(cell));
    }

    m_lastCell = 0;
}

// ============================================================
// mapTo (d-dimensional)
// ============================================================

template<class T>
void gsHPatchRule<T>::mapTo(const gsVector<T> & lower,
                            const gsVector<T> & upper,
                            gsMatrix<T> & nodes,
                            gsVector<T> & weights) const
{
    const CellData & cell = _findCell(lower, upper);

#ifndef NDEBUG
    for (short_t dd = 0; dd < m_dim; ++dd)
        GISMO_ASSERT(lower[dd] >= cell.lowerParam[dd] &&
                     upper[dd] <= cell.upperParam[dd],
                     "gsHPatchRule::mapTo: integration element ["
                     << lower.transpose() << "],[" << upper.transpose()
                     << "] is not contained in selected leaf cell ["
                     << cell.lowerParam.transpose() << "],["
                     << cell.upperParam.transpose() << "]");
#endif

    // node-only comparator for lower_bound on sorted (node,weight) arrays
    struct NodeLess
    { bool operator()(const std::pair<T,T> & p, T v) const { return p.first < v; } };

    index_t size = 1;
    std::vector<gsVector<T>> elNodes(m_dim), elWeights(m_dim);

    for (short_t dd = 0; dd < m_dim; ++dd)
    {
        const std::vector<std::pair<T,T>> & dm = cell.dirMaps[dd];
        const bool isEnd  = (upper[dd] >= cell.upperParam[dd]);
        const bool isDegen = (lower[dd] >= upper[dd]);  // degenerate: boundary element

        // The "pick the single nearest node" branch below is only correct in
        // the fixed direction, whose dirMap stores the two cell-boundary
        // sentinels exactly at cell.lowerParam/upperParam. In any other
        // direction dirMap holds interior Johannessen (or Gauss fallback)
        // nodes, so a degenerate query there would silently pick an
        // arbitrary interior node instead of the boundary value.
        GISMO_ASSERT(!isDegen || dd == m_fixDir,
                     "gsHPatchRule::mapTo: degenerate element in direction " << dd
                     << " (m_fixDir=" << m_fixDir << "); only the fixed direction "
                     << "supports degenerate (boundary) queries.");

        auto itBegin = std::lower_bound(dm.begin(), dm.end(), lower[dd], NodeLess());
        // For degenerate elements (lower==upper, e.g. from boundary iteration):
        // include exactly one node — the one at the boundary value.
        auto itEnd = isDegen
                   ? (itBegin != dm.end() ? std::next(itBegin) : dm.end())
                   : (isEnd ? dm.end()
                            : std::lower_bound(dm.begin(), dm.end(), upper[dd], NodeLess()));

        const index_t cnt = static_cast<index_t>(std::distance(itBegin, itEnd));
        elNodes[dd].resize(cnt);
        elWeights[dd].resize(cnt);
        index_t k = 0;
        for (auto it = itBegin; it != itEnd; ++it, ++k)
        {
            elNodes[dd][k]   = it->first;
            elWeights[dd][k] = it->second;
        }
        size *= cnt;
    }

    nodes.resize(m_dim, size);
    weights.resize(size);
    if (size == 0) return;

    // Tensorize: direction 0 varies fastest (same convention as gsPatchRule)
    index_t nPrev = elNodes[0].size();
    nodes.row(0).head(nPrev)  = elNodes[0].transpose();
    weights.head(nPrev)        = elWeights[0];

    for (short_t dd = 1; dd < m_dim; ++dd)
    {
        const index_t nd   = elNodes[dd].size();
        const index_t nNew = nPrev * nd;

        const gsMatrix<T> prevNodes   = nodes.topRows(dd).leftCols(nPrev); // snapshot
        const gsVector<T> prevWeights = weights.head(nPrev);               // snapshot

        nodes.topRows(dd).leftCols(nNew) = prevNodes.replicate(1, nd);
        for (index_t k = 0; k < nd; ++k)
        {
            nodes.row(dd).segment(k * nPrev, nPrev).setConstant(elNodes[dd][k]);
            weights.segment(k * nPrev, nPrev) = elWeights[dd][k] * prevWeights;
        }
        nPrev = nNew;
    }
}

template<class T>
void gsHPatchRule<T>::mapTo(T startVal, T endVal,
                            gsMatrix<T> & nodes, gsVector<T> & weights) const
{
    GISMO_ASSERT(m_dim == 1, "gsHPatchRule::mapTo(T,T,...) requires dim()==1");
    gsVector<T> lo(1), up(1);
    lo[0] = startVal; up[0] = endVal;
    mapTo(lo, up, nodes, weights);
}

// ============================================================
// Diagnostics
// ============================================================

template<class T>
index_t gsHPatchRule<T>::numPointsTotal() const
{
    index_t total = 0;
    for (const CellData & cell : m_cells)
    {
        index_t nPts = 1;
        for (short_t dd = 0; dd < m_dim; ++dd)
        {
            // dirMaps[m_fixDir] stores the two cell-boundary sentinels used by
            // mapTo to pick a single boundary node for degenerate elements,
            // not two distinct quadrature points; the fixed direction always
            // contributes exactly one point.
            if (dd == m_fixDir) continue;
            nPts *= static_cast<index_t>(cell.dirMaps[dd].size());
        }
        total += nPts;
    }
    return total;
}

template<class T>
index_t gsHPatchRule<T>::numFallbackCells() const
{
    index_t cnt = 0;
    for (const CellData & cell : m_cells)
        if (cell.usedFallback) ++cnt;
    return cnt;
}

// ============================================================
// Hierarchical-domain helpers (dimension dispatch)
// ============================================================

template<class T>
bool gsHPatchRule<T>::isHierarchical(const gsDomain<T> & domain)
{
    switch (domain.dim())
    {
        case 1: return dynamic_cast<const gsHDomain<1,T,index_t>*>(&domain) != nullptr;
        case 2: return dynamic_cast<const gsHDomain<2,T,index_t>*>(&domain) != nullptr;
        case 3: return dynamic_cast<const gsHDomain<3,T,index_t>*>(&domain) != nullptr;
        case 4: return dynamic_cast<const gsHDomain<4,T,index_t>*>(&domain) != nullptr;
        default: return false;
    }
}

template<class T>
typename gsDomain<T>::Ptr
gsHPatchRule<T>::finestTensorDomain(const gsDomain<T> & domain)
{
#define GISMO_HPATCH_FINEST(D)                                                 \
    {                                                                          \
        const gsHDomain<D,T,index_t> * hd =                                    \
            dynamic_cast<const gsHDomain<D,T,index_t>*>(&domain);              \
        if (!hd) return typename gsDomain<T>::Ptr();                           \
        const gsHTensorBasis<D,T> & b = hd->basis();                          \
        return b.tensorLevel(b.maxLevel()).domain();                          \
    }
    switch (domain.dim())
    {
        case 1: GISMO_HPATCH_FINEST(1)
        case 2: GISMO_HPATCH_FINEST(2)
        case 3: GISMO_HPATCH_FINEST(3)
        case 4: GISMO_HPATCH_FINEST(4)
        default: return typename gsDomain<T>::Ptr();
    }
#undef GISMO_HPATCH_FINEST
}

// ============================================================
// makeIntegrationBasis
// ============================================================

template<class T>
template<short_t d>
gsTHBSplineBasis<d,T,true>
gsHPatchRule<T>::makeIntegrationBasis(const gsTHBSplineBasis<d,T,true> & basis)
{
    gsTHBSplineBasis<d,T,true> result(basis);
    gsHElementHelper<d,T> helper(result);

    auto overlap = [](const gsMatrix<index_t,d,2> & a,
                      const gsMatrix<index_t,d,2> & b) -> bool
    {
        for (short_t k = 0; k != d; ++k)
            if (a(k,0) >= b(k,1) || a(k,1) <= b(k,0))
                return false;
        return true;
    };

    // Offline construction cost: O(nElem * nTruncated * nnz(representation)).
    // This keeps the quadrature hot path simple and guarantees that partial
    // overlap with representation-level pieces is not missed by center tests.
    std::vector<index_t> refboxes;
    for (auto & elem : basis.domain()->allElements())
    {
        const index_t elemLevel =
            static_cast<const gsHDomainIterator<T,d,index_t> *>(&elem)->getLevel();
        index_t targetLevel = elemLevel;

        const gsVector<T,d> low = elem.lowerCorner();
        const gsVector<T,d> upp = elem.upperCorner();

        for (typename std::map<index_t, gsSparseVector<T> >::const_iterator it =
                 basis.truncatedBegin(); it != basis.truncatedEnd(); ++it)
        {
            const index_t funIndex = it->first;
            const index_t repLevel = basis.repLevel(funIndex);
            if (repLevel <= targetLevel)
                continue;

            const typename gsHElementHelper<d,T>::element_t repElement =
                helper.toElement(low, upp, repLevel);
            gsMatrix<index_t,d,2> repElementBox;
            repElementBox.col(0) = repElement.lowerCorner();
            repElementBox.col(1) = repElement.upperCorner();

            const gsTensorBSplineBasis<d,T> & repBasis =
                basis.tensorLevel(repLevel);
            const gsSparseVector<T> & coefs = it->second;

            gsMatrix<index_t,d,2> supp;
            for (typename gsSparseVector<T>::InnerIterator coefIt(coefs);
                 coefIt; ++coefIt)
            {
                if (coefIt.value() == T(0))
                    continue;

                repBasis.elementSupport_into(coefIt.index(), supp);
                if (overlap(supp, repElementBox))
                {
                    targetLevel = repLevel;
                    break;
                }
            }
        }

        if (targetLevel <= elemLevel)
            continue;

        auto element = helper.toElement(low, upp, targetLevel);
        const auto & lowCorner = element.lowerCorner();
        const auto & uppCorner = element.upperCorner();

        refboxes.push_back(targetLevel);
        for (short_t k = 0; k != d; ++k) refboxes.push_back(lowCorner(k));
        for (short_t k = 0; k != d; ++k) refboxes.push_back(uppCorner(k));
    }

    if (!refboxes.empty())
        result.refineElements(refboxes);
    return result;
}

// ============================================================
// Private: _findCell + _inside
// ============================================================

template<class T>
bool gsHPatchRule<T>::_contains(const CellData & cell,
                                const gsVector<T> & lower,
                                const gsVector<T> & upper)
{
    for (index_t dd = 0; dd < lower.size(); ++dd)
        if (lower[dd] < cell.lowerParam[dd] || upper[dd] > cell.upperParam[dd])
            return false;
    return true;
}

template<class T>
const typename gsHPatchRule<T>::CellData &
gsHPatchRule<T>::_findCell(const gsVector<T> & lower,
                           const gsVector<T> & upper) const
{
    // Cache: the standard element iteration order is cell-coherent, so a hit
    // on the previously found cell is the common case (amortized O(1)).
    if (m_lastCell < static_cast<index_t>(m_cells.size()) &&
        _contains(m_cells[m_lastCell], lower, upper))
        return m_cells[m_lastCell];

    for (index_t c = 0; c < static_cast<index_t>(m_cells.size()); ++c)
        if (_contains(m_cells[c], lower, upper))
        {
            m_lastCell = c;
            return m_cells[c];
        }

    GISMO_ERROR("gsHPatchRule::_findCell: no cell contains element ["
                << lower.transpose() << "," << upper.transpose() << "]");
}

// ============================================================
// Private static: _ensureEven1D
// ============================================================

template<class T>
gsKnotVector<T> gsHPatchRule<T>::_ensureEven1D(const gsKnotVector<T> & kv)
{
    gsBSplineBasis<T> basis(kv);
    const index_t N = basis.size();
    if (N % 2 == 0) return kv;

    gsKnotVector<T> result = kv;

    typedef typename gsKnotVector<T>::const_iterator kIt;
    std::vector<T> diff;
    {
        kIt prev = result.begin();
        for (kIt it = result.begin() + 1; it != result.end(); ++it)
        {
            diff.push_back(*it - *prev);
            prev = it;
        }
    }

    const T maxSpan = *std::max_element(diff.begin(), diff.end());
    std::vector<index_t> maxIdx;
    for (index_t k = 0; k < (index_t)diff.size(); ++k)
        if (math::abs(diff[k] - maxSpan) / maxSpan < T(1e-15))
            maxIdx.push_back(k);

    const index_t i = maxIdx.at(
        static_cast<index_t>(std::ceil(maxIdx.size() / 2.0)) - 1);
    result.insert((result.at(i) + result.at(i + 1)) / T(2));
    return result;
}

// ============================================================
// Private static: _integrate1D
// ============================================================

template<class T>
std::pair<gsMatrix<T>, gsVector<T>>
gsHPatchRule<T>::_integrate1D(const gsKnotVector<T> & kv)
{
    const index_t q = kv.degree();
    gsBSplineBasis<T> basis(kv);
    const index_t N = basis.size();

    gsMatrix<T> greville;
    kv.greville_into(greville); // 1 x N matrix

    // Exact B-spline integrals: integral N_{i,q} = (kv[i+q+1] - kv[i])/(q+1)
    gsVector<T> integrals(N);
    for (index_t i = 0; i < N; ++i)
        integrals[i] = (kv[i + q + 1] - kv[i]) / T(q + 1);

    return std::make_pair(greville, integrals);
}

// ============================================================
// Private static: _compute1D  (dense Newton)
// ============================================================

template<class T>
bool gsHPatchRule<T>::_compute1D(const gsKnotVector<T> & evenKv,
                                 gsVector<T> & nodes,
                                 gsVector<T> & weights,
                                 T       tol,
                                 index_t itMax)
{
    gsBSplineBasis<T> basis(evenKv);
    const index_t size = basis.size();
    GISMO_ENSURE(size % 2 == 0, "gsHPatchRule::_compute1D: basis dimension must be even");

    const index_t n = size / 2;

    gsMatrix<T> greville;
    gsVector<T> integrals;
    std::tie(greville, integrals) = _integrate1D(evenKv);

    nodes.resize(n);
    weights.resize(n);
    for (index_t k = 0; k < n; ++k)
    {
        nodes[k]   = T(0.5) * (greville(0, 2*k) + greville(0, 2*k+1));
        weights[k] = integrals[2*k] + integrals[2*k+1];
    }

    gsMatrix<index_t> actives;
    gsMatrix<T>       vals(size, n), dvals(size, n);
    gsMatrix<T>       vals_tmp, dvals_tmp;
    gsVector<T>       res, update(size);
    gsMatrix<T>       dres(size, size);            // allocate once (F3)

    for (index_t it = 0; it < itMax; ++it)
    {
        vals.setZero();
        dvals.setZero();

        basis.active_into(nodes.transpose(), actives);
        basis.eval_into  (nodes.transpose(), vals_tmp);
        basis.deriv_into (nodes.transpose(), dvals_tmp);

        for (index_t act = 0; act < actives.rows(); ++act)
            for (index_t pt = 0; pt < actives.cols(); ++pt)
            {
                vals (actives(act, pt), pt) = vals_tmp (act, pt);
                dvals(actives(act, pt), pt) = dvals_tmp(act, pt);
            }

        res = vals * weights - integrals;
        if (res.norm() < tol) return true;

        // Jacobian [dF/dw | dF/dxi]  (dense; small, square)
        dres.leftCols(n)  = vals;
        dres.rightCols(n) = dvals * weights.asDiagonal();

        update = dres.partialPivLu().solve(-res);

        weights += update.head(n);
        nodes   += update.tail(n);

        if (nodes.minCoeff() <= evenKv.first()) return false;
        if (nodes.maxCoeff() >= evenKv.last())  return false;
    }
    return false; // max iterations reached
}

// ============================================================
// Private static: _fillGaussFallback1D + _finalizeDirMap
// ============================================================

template<class T>
void gsHPatchRule<T>::_fillGaussFallback1D(std::vector<std::pair<T,T>> & dirMap,
                                           const gsKnotVector<T> & kv,
                                           index_t lowerIdx,
                                           index_t upperIdx,
                                           index_t nq)
{
    gsGaussRule<T> gaussRule1D(nq);
    for (index_t e = lowerIdx; e < upperIdx; ++e)
    {
        const T a = kv(e);     // unique knot at position e
        const T b = kv(e + 1); // unique knot at position e+1

        gsMatrix<T> elNodes;
        gsVector<T> elWeights;
        gaussRule1D.mapTo(a, b, elNodes, elWeights);

        for (index_t k = 0; k < nq; ++k)
            dirMap.push_back(std::make_pair(elNodes(0,k), elWeights[k]));
    }
}

template<class T>
void gsHPatchRule<T>::_finalizeDirMap(std::vector<std::pair<T,T>> & dirMap)
{
    std::sort(dirMap.begin(), dirMap.end(),
              [](const std::pair<T,T> & a, const std::pair<T,T> & b)
              { return a.first < b.first; });
}

} // namespace gismo
