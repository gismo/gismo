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
                              const short_t  fixDir,
                              gsHPatchRuleCache<T> * cache)
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

    gsHPatchRuleCache<T> local;
    gsHPatchRuleCache<T> & c = cache ? *cache : local;

    switch (m_dim)
    {
        case 1: _build<1>(domain, c); break;
        case 2: _build<2>(domain, c); break;
        case 3: _build<3>(domain, c); break;
        case 4: _build<4>(domain, c); break;
        default:
            GISMO_ERROR("gsHPatchRule: unsupported domain dimension " << m_dim);
    }
}

// ============================================================
// Dimension-specific construction
// ============================================================

template<class T>
template<short_t d>
void gsHPatchRule<T>::_build(const gsDomain<T> & domain, gsHPatchRuleCache<T> & cache)
{
    const gsHDomain<d,T,index_t> * hd =
        dynamic_cast<const gsHDomain<d,T,index_t>*>(&domain);
    GISMO_ENSURE(hd != nullptr,
                 "gsHPatchRule requires a hierarchical domain (gsHDomain).");

    const gsHTensorBasis<d,T> & basis = hd->basis();
    GISMO_ASSERT(!basis.manualLevels(),
                 "gsHPatchRule: manualLevels bases are not supported.");

    m_maxLevel = basis.maxLevel();

    // Element-wise Gauss with nq points per element is exact to degree 2nq-1 >= m_deg.
    const index_t nq = m_deg / 2 + 1;

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
        std::vector<bool>            gaussDir(d, false); // direction uses element-wise Gauss

        for (short_t dd = 0; dd < d; ++dd)
        {
            const gsKnotVector<T> & kvLev =
                basis.tensorLevel(l).component(dd).knots();

            cell.lowerIndex[dd] = leaf.lowerCorner()[dd];
            cell.upperIndex[dd] = leaf.upperCorner()[dd];
            cell.lowerParam[dd] = kvLev(cell.lowerIndex[dd]);
            cell.upperParam[dd] = kvLev(cell.upperIndex[dd]);

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

            // The patch rule is used in this direction only if it has strictly fewer
            // nodes than element-wise Gauss; a tie costs the same and skips the solve.
            // It never has more: a block of e elements between two knots of multiplicity
            // m_deg+1 has at most e*m_deg+1 functions, hence after even-ing at most
            // ceil((e*m_deg+1)/2) <= e*nq nodes.
            const index_t nElemDD = cell.upperIndex[dd] - cell.lowerIndex[dd];
            GISMO_ASSERT(n1d[dd] <= nElemDD * nq,
                         "gsHPatchRule: the Johannessen rule has more nodes than element-wise Gauss in direction " << dd);
            gaussDir[dd] = (n1d[dd] >= nElemDD * nq);
        }

        // Per direction either the Johannessen rule on evenKvs[dd] or nq-point Gauss on
        // every level-l element of the leaf.  On the leaf the target space is the tensor
        // product S_1 x ... x S_d of univariate spaces S_k; Q_k is exact on S_k in both
        // cases (the even-ized space contains S_k; piecewise polynomials of degree
        // <= 2nq-1 >= m_deg on the leaf's elements contain S_k).  For separable
        // f = f_1...f_d with f_k in S_k the tensor rule gives prod Q_k(f_k) = prod int f_k
        // = int f, and by linearity it is exact on S, whichever direction uses which rule.
        // The Johannessen rule of a direction is solved on the even target knot vector mapped
        // affinely to [0,1] and scaled back to the leaf [a,a+L] (node a + L*xi, weight
        // L*omega): B-splines are invariant under the map and the moments I_i scale with L,
        // so a rule that is exact on [0,1] is exact on the leaf.
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

            bool newtonFailed = false;
            if (!gaussDir[dd])
            {
                // Even target knot vector mapped to [0,1]: division per knot, monotone and
                // exact at both ends (z = a gives 0, z = b gives 1).
                const T a = cell.lowerParam[dd];
                const T L = cell.upperParam[dd] - a;
                const std::vector<T> & physKnots = evenKvs[dd];
                std::vector<T> refKnots(physKnots.size());
                for (size_t i = 0; i != physKnots.size(); ++i)
                    refKnots[i] = (physKnots[i] - a) / L;
                const gsKnotVector<T> refKv(refKnots, evenKvs[dd].degree());

                const typename gsHPatchRuleCache<T>::Entry & rule =
                    _referenceRule(refKv, (std::min)(L, T(1)), cache);
                newtonFailed = !rule.converged;
                if (!newtonFailed)
                {
                    // Johannessen results, sorted by node
                    cell.dirMaps[dd].reserve(n1d[dd]);
                    for (index_t k = 0; k < n1d[dd]; ++k)
                        cell.dirMaps[dd].push_back(
                            std::make_pair(a + L * rule.nodes[k], L * rule.weights[k]));
                }
            }

            if (gaussDir[dd] || newtonFailed)
            {
                // A tie (no solve) is not a fallback: same node count, no loss.
                if (newtonFailed)
                    cell.usedFallback = true;
                const gsKnotVector<T> & kvLev =
                    basis.tensorLevel(l).component(dd).knots();
                _fillGaussFallback1D(cell.dirMaps[dd], kvLev,
                                     cell.lowerIndex[dd], cell.upperIndex[dd], nq);
            }
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
    GISMO_ASSERT(lower.size() == m_dim && upper.size() == m_dim,
                 "gsHPatchRule::mapTo: element corners must have dim() entries");
    _mapTo(lower.data(), upper.data(), nodes, weights);
}

template<class T>
void gsHPatchRule<T>::mapTo(T startVal, T endVal,
                            gsMatrix<T> & nodes, gsVector<T> & weights) const
{
    GISMO_ASSERT(m_dim == 1, "gsHPatchRule::mapTo(T,T,...) requires dim()==1");
    _mapTo(&startVal, &endVal, nodes, weights);
}

template<class T>
void gsHPatchRule<T>::_mapTo(const T * lower, const T * upper,
                             gsMatrix<T> & nodes, gsVector<T> & weights) const
{
    const index_t c = _findCell(lower, upper);
    if (c < 0)
        _missError(lower, upper, m_dim);
    const CellData & cell = m_cells[c];

#ifndef NDEBUG
    for (index_t dd = 0; dd < m_dim; ++dd)
        GISMO_ASSERT(lower[dd] >= cell.lowerParam[dd] &&
                     upper[dd] <= cell.upperParam[dd],
                     "gsHPatchRule::mapTo: integration element is not contained "
                     "in selected leaf cell in direction " << dd << ": ["
                     << lower[dd] << "," << upper[dd] << "] not in ["
                     << cell.lowerParam[dd] << "," << cell.upperParam[dd] << "]");
#endif

    // node-only comparator for lower_bound on sorted (node,weight) arrays
    struct NodeLess
    { bool operator()(const std::pair<T,T> & p, T v) const { return p.first < v; } };

    // The constructor accepts dim() = 1..4 only.
    const std::pair<T,T> * first[4];  // start of the selected range in direction dd
    index_t                cnt[4];    // its length
    index_t size = 1;

    for (index_t dd = 0; dd < m_dim; ++dd)
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

        cnt[dd]   = static_cast<index_t>(std::distance(itBegin, itEnd));
        first[dd] = (cnt[dd] > 0 ? &*itBegin : nullptr);
        size *= cnt[dd];
    }

    nodes.resize(m_dim, size);
    weights.resize(size);
    if (size == 0) return;

    // Tensorize: direction 0 varies fastest (same convention as gsPatchRule),
    // column i <-> (idx[0],...,idx[m_dim-1]) with i = idx[0] + cnt[0]*(idx[1] + ...).
    // The weight is accumulated as w_{d-1} * (... * (w_1 * w_0)).
    index_t idx[4] = {0, 0, 0, 0};
    for (index_t i = 0; i < size; ++i)
    {
        T w = first[0][idx[0]].second;
        nodes(0, i) = first[0][idx[0]].first;
        for (index_t dd = 1; dd < m_dim; ++dd)
        {
            nodes(dd, i) = first[dd][idx[dd]].first;
            w = first[dd][idx[dd]].second * w;
        }
        weights[i] = w;

        for (index_t dd = 0; dd < m_dim; ++dd)
        {
            if (++idx[dd] < cnt[dd]) break;
            idx[dd] = 0;
        }
    }
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
// Shape cache and reference-interval solve
// ============================================================

template<class T>
gsHPatchRuleCache<T>::gsHPatchRuleCache()
: m_hits(0), m_solves(0), m_failures(0)
{ }

template<class T>
index_t gsHPatchRuleCache<T>::numHits() const { return m_hits; }

template<class T>
index_t gsHPatchRuleCache<T>::numSolves() const { return m_solves; }

template<class T>
index_t gsHPatchRuleCache<T>::numFailures() const { return m_failures; }

template<class T>
index_t gsHPatchRuleCache<T>::size() const
{ return static_cast<index_t>(m_shapes.size()); }

template<class T>
const std::vector< std::pair<gsKnotVector<T>, T> > &
gsHPatchRuleCache<T>::failedKeys() const
{ return m_failedKeys; }

template<class T>
void gsHPatchRuleCache<T>::clear()
{
    m_shapes.clear();
    m_failedKeys.clear();
    m_hits = m_solves = m_failures = 0;
}

template<class T>
const typename gsHPatchRuleCache<T>::Entry &
gsHPatchRule<T>::_referenceRule(const gsKnotVector<T> & refKv, T Lc,
                                gsHPatchRuleCache<T> & cache)
{
    typedef typename gsHPatchRuleCache<T>::Key   Key;
    typedef typename gsHPatchRuleCache<T>::Entry Entry;

    const Key key(Lc, std::make_pair(static_cast<index_t>(refKv.degree()),
                                     static_cast<const std::vector<T> &>(refKv)));
    typename std::map<Key, Entry>::const_iterator it = cache.m_shapes.find(key);
    if (it != cache.m_shapes.end())
    {
        ++cache.m_hits;
        return it->second;
    }

    ++cache.m_solves;
    Entry entry;
    entry.converged = _compute1D(refKv, entry.nodes, entry.weights, T(1e-14) / Lc);
    if (!entry.converged)
    {
        entry.nodes.resize(0);
        entry.weights.resize(0);
        ++cache.m_failures;
        cache.m_failedKeys.push_back(std::make_pair(refKv, Lc));
    }
    return cache.m_shapes.insert(std::make_pair(key, give(entry))).first->second;
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
                                const T * lower, const T * upper, index_t dim)
{
    for (index_t dd = 0; dd < dim; ++dd)
        if (lower[dd] < cell.lowerParam[dd] || upper[dd] > cell.upperParam[dd])
            return false;
    return true;
}

template<class T>
index_t gsHPatchRule<T>::_findCell(const T * lower, const T * upper) const
{
    const index_t L = static_cast<index_t>(m_cells.size());
    if (L == 0)
        return -1;

    // Probe the last hit (k = 0), then its successor (k = 1), then continue
    // cyclically: the leaf order of volume iteration makes k <= 1 the common
    // case, and the wrap-around finds leaf 0 right after the last leaf.
    const index_t c = m_lastCell;
    for (index_t k = 0; k < L; ++k)
    {
        index_t j = c + k;
        if (j >= L) j -= L;
        if (_contains(m_cells[j], lower, upper, m_dim))
        {
            m_lastCell = j;
            return j;
        }
    }
    return -1;
}

template<class T>
void gsHPatchRule<T>::_missError(const T * lower, const T * upper, index_t dim)
{
    gsVector<T> lo(dim), up(dim);
    for (index_t dd = 0; dd < dim; ++dd)
    {
        lo[dd] = lower[dd];
        up[dd] = upper[dd];
    }
    GISMO_ERROR("gsHPatchRule::mapTo: no cell contains element ["
                << lo.transpose() << "," << up.transpose() << "]");
}

// ============================================================
// Private static: _ensureEven1D
// ============================================================

template<class T>
gsKnotVector<T> gsHPatchRule<T>::_ensureEven1D(const gsKnotVector<T> & kv)
{
    const index_t q = kv.degree();
    const std::vector<T>       uk   = kv.unique();
    const std::vector<index_t> mult = kv.multiplicities();
    const index_t nU = static_cast<index_t>(uk.size());

    // A knot of multiplicity q+1 splits the space into independent blocks, and a node
    // cannot serve two blocks. A block between two such knots (the end knots included)
    // with u unique knots and interior multiplicities m_j has dimension
    // q+1 + sum_j m_j, and needs an even dimension of its own.
    std::vector<T> newKnots;
    index_t first = 0;
    for (index_t u = 1; u < nU; ++u)
    {
        if (u != nU - 1 && mult[u] < q + 1)
            continue;

        index_t dim = q + 1;
        for (index_t j = first + 1; j < u; ++j)
            dim += mult[j];

        if (dim % 2 != 0)
        {
            // Split the middle one of the longest spans of the block
            T maxSpan = T(0);
            for (index_t j = first; j < u; ++j)
                maxSpan = (std::max)(maxSpan, uk[j + 1] - uk[j]);
            std::vector<index_t> maxIdx;
            for (index_t j = first; j < u; ++j)
                if (math::abs(uk[j + 1] - uk[j] - maxSpan) / maxSpan < T(1e-15))
                    maxIdx.push_back(j);
            const index_t i = maxIdx.at(
                static_cast<index_t>(std::ceil(maxIdx.size() / 2.0)) - 1);
            newKnots.push_back((uk[i] + uk[i + 1]) / T(2));
        }
        first = u;
    }

    gsKnotVector<T> result = kv;
    for (size_t i = 0; i != newKnots.size(); ++i)
        result.insert(newKnots[i]);
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
// Private static: _compute1D  (Newton, banded Jacobian) + _bandLUSolve
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

    const index_t q = evenKv.degree();
    const gsVector<T> nodes0 = nodes, weights0 = weights;

    // Local evaluations: (q+1) x n blocks holding, per node, the q+1 active
    // B-splines (consecutive, ascending indices) and their values/derivatives.
    gsMatrix<index_t> actives;
    gsMatrix<T>       vals, dvals;
    gsVector<T>       res(size), update(size);
    gsMatrix<T>       AB; // LAPACK band storage of the Jacobian

    for (index_t it = 0; it < itMax; ++it)
    {
        basis.active_into(nodes.transpose(), actives);
        basis.eval_into  (nodes.transpose(), vals);
        basis.deriv_into (nodes.transpose(), dvals);

        // Residual F_i = sum_k w_k N_i(xi_k) - I_i
        res = -integrals;
        for (index_t k = 0; k < n; ++k)
            for (index_t j = 0; j <= q; ++j)
                res[actives(j, k)] += vals(j, k) * weights[k];

        const T resNormN = res.norm();
        if (!math::isfinite(resNormN)) break;
        if (resNormN < tol) return true;

        // Bandwidths of the Jacobian with the unknowns ordered (w_0, xi_0, w_1, xi_1, ...)
        index_t kl = 0, ku = 0;
        for (index_t k = 0; k < n; ++k)
        {
            kl = (std::max)(kl, actives(0, k) + q - 2*k);
            ku = (std::max)(ku, 2*k + 1 - actives(0, k));
        }

        // Jacobian columns 2k: dF/dw_k = N_i(xi_k), 2k+1: dF/dxi_k = w_k N_i'(xi_k)
        const index_t kv = kl + ku;
        AB.setZero(2*kl + ku + 1, size);
        for (index_t k = 0; k < n; ++k)
            for (index_t j = 0; j <= q; ++j)
            {
                const index_t i = actives(j, k);
                AB(kv + i - 2*k    , 2*k    ) = vals (j, k);
                AB(kv + i - 2*k - 1, 2*k + 1) = dvals(j, k) * weights[k];
            }

        update = -res;
        if (!_bandLUSolve(AB, kl, ku, update)) break;
        if (!math::isfinite(update.sum())) break;

        for (index_t k = 0; k < n; ++k)
        {
            weights[k] += update[2*k];
            nodes[k]   += update[2*k + 1];
        }

        if (nodes.minCoeff() <= evenKv.first()) break;
        if (nodes.maxCoeff() >= evenKv.last())  break;
    }

    // Levenberg-Marquardt restart from the starting guess: minimises ||F||^2 with the
    // damped normal equations (J^T J + mu diag(J^T J)) d = -J^T F, which stay solvable
    // when J is nearly singular and bound the step where the full Newton step leaves the
    // basin of attraction. J^T J inherits the band structure of J (half-bandwidth hb).
    nodes = nodes0;
    weights = weights0;
    const T lo = evenKv.first(), hi = evenKv.last();
    gsMatrix<T> Hb, ABm;
    gsVector<T> grad(size), dH(size), delta(size), resNew(size), nodesNew(n), weightsNew(n);
    T mu = T(0);
    for (index_t it = 0; it < itMax; ++it)
    {
        basis.active_into(nodes.transpose(), actives);
        basis.eval_into  (nodes.transpose(), vals);
        basis.deriv_into (nodes.transpose(), dvals);

        res = -integrals;
        for (index_t k = 0; k < n; ++k)
            for (index_t j = 0; j <= q; ++j)
                res[actives(j, k)] += vals(j, k) * weights[k];

        const T resNorm = res.norm();
        if (!math::isfinite(resNorm)) return false;
        if (resNorm < tol) return true;

        // Gradient J^T F; the column of xi_k carries the factor w_k
        for (index_t k = 0; k < n; ++k)
        {
            T gw = T(0), gx = T(0);
            for (index_t j = 0; j <= q; ++j)
            {
                gw += vals (j, k) * res[actives(j, k)];
                gx += dvals(j, k) * res[actives(j, k)];
            }
            grad[2*k]     = gw;
            grad[2*k + 1] = weights[k] * gx;
        }

        // Half-bandwidth of J^T J: columns of nodes j and k couple iff their active ranges overlap
        // For nodes in ascending order of their first active index the coupled nodes of k form
        // a contiguous run, so the scans stop at its end (O(n q)); otherwise all pairs are scanned.
        bool ordered = true;
        for (index_t k = 1; k < n; ++k)
            if (actives(0, k) < actives(0, k - 1)) { ordered = false; break; }

        index_t hb = 1;
        for (index_t k = 0; k < n; ++k)
            for (index_t j = k + 1; j < n; ++j)
            {
                if (math::abs(actives(0, j) - actives(0, k)) <= q)
                    hb = (std::max)(hb, 2 * (j - k) + 1);
                else if (ordered)
                    break;
            }

        // J^T J in band storage (LAPACK layout of _bandLUSolve with kl = ku = hb)
        const index_t kvb = 2 * hb;
        Hb.setZero(3 * hb + 1, size);
        for (index_t k = 0; k < n; ++k)
            for (index_t j = k; j < n; ++j)
            {
                const index_t ak = actives(0, k), aj = actives(0, j);
                if (math::abs(aj - ak) > q)
                {
                    if (ordered) break;
                    continue;
                }
                T b00 = T(0), b01 = T(0), b10 = T(0), b11 = T(0);
                const index_t rLo = (std::max)(ak, aj), rHi = (std::min)(ak, aj) + q;
                for (index_t r = rLo; r <= rHi; ++r)
                {
                    const T vk = vals(r - ak, k), dk = dvals(r - ak, k) * weights[k];
                    const T vj = vals(r - aj, j), dj = dvals(r - aj, j) * weights[j];
                    b00 += vk * vj;  b01 += vk * dj;
                    b10 += dk * vj;  b11 += dk * dj;
                }
                // entry (row, col) of J^T J lives at Hb(kvb + row - col, col)
                Hb(kvb + 2*k     - 2*j,     2*j    ) = b00;
                Hb(kvb + 2*k     - (2*j+1), 2*j + 1) = b01;
                Hb(kvb + 2*k + 1 - 2*j,     2*j    ) = b10;
                Hb(kvb + 2*k + 1 - (2*j+1), 2*j + 1) = b11;
                Hb(kvb + 2*j     - 2*k,     2*k    ) = b00;
                Hb(kvb + 2*j + 1 - 2*k,     2*k    ) = b01;
                Hb(kvb + 2*j     - (2*k+1), 2*k + 1) = b10;
                Hb(kvb + 2*j + 1 - (2*k+1), 2*k + 1) = b11;
            }
        for (index_t i = 0; i < size; ++i)
            dH[i] = (std::max)(Hb(kvb, i), std::numeric_limits<T>::min());
        if (it == 0)
            mu = T(1e-3) * dH.maxCoeff();

        bool accepted = false;
        for (index_t trial = 0; trial < 40 && !accepted; ++trial)
        {
            ABm = Hb;
            for (index_t i = 0; i < size; ++i)
                ABm(kvb, i) += mu * dH[i];
            delta = -grad;
            if (_bandLUSolve(ABm, hb, hb, delta))
            {
                for (index_t k = 0; k < n; ++k)
                {
                    weightsNew[k] = weights[k] + delta[2*k];
                    nodesNew[k]   = nodes[k]   + delta[2*k + 1];
                }
                if (nodesNew.allFinite() && weightsNew.allFinite() &&
                    nodesNew.minCoeff() > lo && nodesNew.maxCoeff() < hi)
                {
                    basis.active_into(nodesNew.transpose(), actives);
                    basis.eval_into  (nodesNew.transpose(), vals);
                    resNew = -integrals;
                    for (index_t k = 0; k < n; ++k)
                        for (index_t j = 0; j <= q; ++j)
                            resNew[actives(j, k)] += vals(j, k) * weightsNew[k];
                    if (resNew.allFinite() && resNew.norm() < resNorm)
                    {
                        nodes = nodesNew;
                        weights = weightsNew;
                        mu = (std::max)(mu / T(3), T(1e-14));
                        accepted = true;
                        break;
                    }
                }
            }
            mu *= T(4);
        }
        if (!accepted) return false; // no damped step lowers the residual
    }
    return false; // max iterations reached
}

template<class T>
bool gsHPatchRule<T>::_bandLUSolve(gsMatrix<T> & AB, index_t kl, index_t ku, gsVector<T> & b)
{
    const index_t N  = AB.cols();
    const index_t kv = kl + ku; // row of the diagonal in AB; U has upper bandwidth kv
    std::vector<index_t> piv(N);

    // Factorization P A = L U, row partial pivoting
    for (index_t j = 0; j < N; ++j)
    {
        const index_t iMax = (std::min)(N - 1, j + kl);
        const index_t cMax = (std::min)(N - 1, j + kv);

        index_t p = j;
        T pmax = math::abs(AB(kv, j));
        for (index_t i = j + 1; i <= iMax; ++i)
            if (math::abs(AB(kv + i - j, j)) > pmax)
            {
                pmax = math::abs(AB(kv + i - j, j));
                p = i;
            }
        if (!(pmax > T(0)) || !(pmax <= (std::numeric_limits<T>::max)())) // zero, NaN or inf pivot
            return false;
        piv[j] = p;

        if (p != j)
            for (index_t c = j; c <= cMax; ++c)
                std::swap(AB(kv + j - c, c), AB(kv + p - c, c));

        const T d = AB(kv, j);
        for (index_t i = j + 1; i <= iMax; ++i)
            AB(kv + i - j, j) /= d;

        for (index_t c = j + 1; c <= cMax; ++c)
        {
            const T u = AB(kv + j - c, c);
            for (index_t i = j + 1; i <= iMax; ++i)
                AB(kv + i - c, c) -= AB(kv + i - j, j) * u;
        }
    }

    // Forward substitution L y = P b
    for (index_t j = 0; j < N; ++j)
    {
        if (piv[j] != j)
            std::swap(b[j], b[piv[j]]);
        const index_t iMax = (std::min)(N - 1, j + kl);
        for (index_t i = j + 1; i <= iMax; ++i)
            b[i] -= AB(kv + i - j, j) * b[j];
    }

    // Back substitution U x = y
    for (index_t j = N - 1; j >= 0; --j)
    {
        b[j] /= AB(kv, j);
        for (index_t i = (std::max)(index_t(0), j - kv); i < j; ++i)
            b[i] -= AB(kv + i - j, j) * b[j];
    }
    return true;
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
