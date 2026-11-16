/** @file hierarchical_macro_mass_rule_tensor.cpp

    @brief Tensorized Johannessen quadrature rule for each atomic cell of a
    hierarchical THB-spline domain.

    For each gsHTree leaf (= one atomic cell = one macro region):
      1. Builds a local open target knot vector per direction:
           degree q = 2p, regularity p-1 (internal multiplicity p+1),
           using the unique breakpoints of the leaf at its level.
      2. Applies the 1D Johannessen Newton rule per direction independently,
         obtaining (nodes, weights) in each direction.
      3. Tensorizes the two 1D rules into a 2D quadrature rule with
           n_tensor = ceil(Nx/2) * ceil(Ny/2)
         points, where Nx, Ny are the actual 1D target basis dimensions
         (not derived from formula, but from basis.size() after even-fixing).
      4. Falls back to element-wise Gauss only when n_tensor >= n_Gauss
         (i.e., on very small leaves with no efficiency gain).
      5. Validates the rule against level tensor-basis products on the leaf.

    This replaces the 2D nonlinear LM solver in hierarchical_macro_mass_rule_single.cpp.
    Asymptotically (large leaves), the tensor rule uses ~1/4 the points of standard Gauss.

    NOTE on future extension:
      The local target knot vector uses regularity p-1 (internal multiplicity p+1).
      This is correct for the current simple-knot (max-continuity) THB test.
      For solution knots with internal multiplicity m_sol > 1, the target regularity
      should be (p - m_sol), giving internal target multiplicity (p + m_sol).
      A hook for this is marked in buildTensorRule().

    Usage:
        ./hierarchical_macro_mass_rule_tensor [-c N] [-o overQ] [--plot] [--all]

    Options:
        -c N        Atomic cell to validate in detail (-1 = auto, default)
        -o N        Extra Gauss points per direction for reference integration (default: 4)
        --plot      Write ParaView mesh file
        --all       Run detailed validation on every cell

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <unordered_map>

using namespace gismo;

// ============================================================
// Shared infrastructure
// ============================================================

struct ElementBox
{
    gsVector<real_t> lower;
    gsVector<real_t> upper;
};

struct AtomicCell
{
    index_t            level;
    gsVector<index_t>  lowerIndex;
    gsVector<index_t>  upperIndex;
    gsVector<real_t>   lowerParam;
    gsVector<real_t>   upperParam;
    std::vector<ElementBox> elements;
};

template<short_t d, class T>
std::vector<AtomicCell> collectAtomicCells(const gsHTensorBasis<d,T>& thb)
{
    GISMO_ASSERT(!thb.manualLevels(),
        "collectAtomicCells: manualLevels bases not supported here.");

    std::vector<AtomicCell> cells;
    auto leaf = thb.tree().beginLeafIterator();

    for (; leaf.good(); leaf.next())
    {
        AtomicCell cell;
        cell.level      = leaf.level();
        cell.lowerIndex = leaf.lowerCorner();
        cell.upperIndex = leaf.upperCorner();
        cell.lowerParam.resize(d);
        cell.upperParam.resize(d);

        for (short_t dd = 0; dd < d; ++dd)
        {
            const gsKnotVector<T>& kv =
                thb.tensorLevel(cell.level).component(dd).knots();
            cell.lowerParam[dd] = kv(cell.lowerIndex[dd]);
            cell.upperParam[dd] = kv(cell.upperIndex[dd]);
        }

        gsVector<index_t> alpha(d);
        alpha = cell.lowerIndex;

        bool more = true;
        while (more)
        {
            ElementBox elem;
            elem.lower.resize(d);
            elem.upper.resize(d);
            for (short_t dd = 0; dd < d; ++dd)
            {
                const gsKnotVector<T>& kv =
                    thb.tensorLevel(cell.level).component(dd).knots();
                elem.lower[dd] = kv(alpha[dd]);
                elem.upper[dd] = kv(alpha[dd] + 1);
            }
            cell.elements.push_back(elem);

            more = false;
            for (short_t dd = d - 1; dd >= 0; --dd)
            {
                if (alpha[dd] < cell.upperIndex[dd] - 1)
                { ++alpha[dd]; more = true; break; }
                alpha[dd] = cell.lowerIndex[dd];
            }
        }

        cells.push_back(cell);
    }
    return cells;
}

// ============================================================
// 1D Johannessen Newton rule
// (logic extracted from gsPatchRule::_init/_integrate/_compute)
// ============================================================

/**
 * @brief Ensure the 1D basis has even dimension by inserting one knot.
 *
 * Checks N = gsBSplineBasis<T>(kv).size() (NOT kv.size() or kv.uSize()).
 * If N is odd, inserts a knot at the midpoint of the centered-largest span,
 * matching the tie-break in gsPatchRule::_init.
 */
template<class T>
gsKnotVector<T> ensureEven1D(const gsKnotVector<T>& kv)
{
    gsBSplineBasis<T> basis(kv);
    const index_t N = basis.size();
    if (N % 2 == 0) return kv;

    gsKnotVector<T> result = kv;

    // Compute gaps between all consecutive knot values (including repeats)
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

    // Find centered index among all maximum-span positions
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

/**
 * @brief Compute Greville abscissae and exact B-spline integrals.
 *
 * Uses the exact formula:  integral N_{i,q}(x) dx = (tau[i+q+1] - tau[i]) / (q+1)
 * where tau[i] = kv[i] is the repeated knot vector.
 *
 * @return {greville (1 x N matrix), integrals (N-vector)}
 */
template<class T>
std::pair<gsMatrix<T>, gsVector<T> > integrate1D(const gsKnotVector<T>& kv)
{
    const index_t q = kv.degree();
    gsBSplineBasis<T> basis(kv);
    const index_t N = basis.size();

    gsMatrix<T> greville;
    kv.greville_into(greville);  // returns 1 x N matrix

    // Exact B-spline integral via knot differences (repeated knot access: kv[i])
    gsVector<T> integrals(N);
    for (index_t i = 0; i < N; ++i)
        integrals[i] = (kv[i + q + 1] - kv[i]) / T(q + 1);

    return std::make_pair(greville, integrals);
}

/**
 * @brief 1D Johannessen Newton solver.
 *
 * Solves:  sum_{q=1}^{n}  w_q * N_i(x_q) = integral N_i,   i = 1..N
 * for n = N/2 unknown (node, weight) pairs.
 *
 * Initial guess: Greville midpoints as nodes, paired integral sums as weights.
 * Convergence check on F = vals*weights - integrals BEFORE applying each update.
 * vals and dvals are reset to zero every iteration to avoid stale entries when
 * a node crosses a support boundary.
 *
 * Returns false (recoverable) if nodes leave the domain or Newton does not converge.
 * Returns true on success.
 *
 * @param kv      Target knot vector (basis dimension N must be even — programmer error)
 * @param nodes   Output: n-vector of quadrature nodes
 * @param weights Output: n-vector of quadrature weights
 * @param tol     Newton residual tolerance on ||F||
 * @param itMax   Maximum Newton iterations
 */
template<class T>
bool compute1D(const gsKnotVector<T>& kv,
               gsVector<T>&           nodes,
               gsVector<T>&           weights,
               T                      tol   = T(1e-10),
               index_t                itMax = 100)
{
    gsBSplineBasis<T> basis(kv);
    const index_t size = basis.size();
    GISMO_ENSURE(size % 2 == 0, "Basis dimension must be even in compute1D");

    const index_t n = size / 2;

    // Initial guess: Greville midpoints and paired integrals.
    // greville_into returns a 1 x N matrix, so access as greville(0, k).
    gsMatrix<T> greville;
    gsVector<T> integrals;
    std::tie(greville, integrals) = integrate1D(kv);

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
    gsMatrix<T>       dres;

    typename gsSparseSolver<T>::QR solver;
    index_t it = 0;

    for (; it < itMax; ++it)
    {
        vals.setZero();   // reset to avoid stale entries from previous iteration
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

        // Compute residual for current iterate
        res = vals * weights - integrals;

        // Check convergence BEFORE computing the update
        if (res.norm() < tol) return true;

        // Jacobian: [dF/dw | dF/dxi] in R^{size x size}
        dres = gsMatrix<T>::Zero(size, size);
        dres.block(0,    0, size,    n) = vals;
        dres.block(0,    n, size,    n) = dvals * weights.asDiagonal();

        solver.compute(dres.sparseView());
        update = solver.solve(-res);

        weights += update.head(n);
        nodes   += update.tail(n);

        if (nodes.minCoeff() <= kv.first()) return false;
        if (nodes.maxCoeff() >= kv.last())  return false;
    }
    return false; // max iterations reached
}

/**
 * @brief Run 1D Johannessen Newton on an already-even target knot vector.
 *
 * @param evenKv  Even-dimensioned open target KV (from ensureEven1D)
 * @param nodes   Output: quadrature nodes
 * @param weights Output: quadrature weights
 * @return true on success, false on Newton failure
 */
template<class T>
bool johannessen1DFromEvenKv(const gsKnotVector<T>& evenKv,
                              gsVector<T>&           nodes,
                              gsVector<T>&           weights)
{
    return compute1D(evenKv, nodes, weights);
}

// ============================================================
// Tensor rule builder
// ============================================================

template<class T>
struct TensorRuleResult
{
    gsMatrix<T>              pts;              // d x nUsed
    gsVector<T>              wts;              // nUsed
    index_t                  nCandidateTensor; // Johannessen tensor point count (before decision)
    index_t                  nUsed;            // actual rule points = pts.cols()
    index_t                  nGauss;           // standard element-wise Gauss point count
    std::vector<index_t>     n1d;              // per-direction 1D Johannessen point counts
    bool                     usedFallback;     // true if Gauss fallback was applied
    bool                     newtonFailed;     // true if fallback was due to Newton failure
    std::string              fallbackReason;   // "none", "count", or "newton"
};

/**
 * @brief Build element-wise Gauss quadrature rule (fallback for buildTensorRule).
 *
 * Used when nCandidateTensor > nGauss (count fallback) or Newton fails.
 */
template<short_t d, class T>
TensorRuleResult<T>
buildGaussFallback(const AtomicCell&            cell,
                   const std::vector<index_t>&  pPerDir,
                   index_t                      nCandidateTensor,
                   index_t                      nGauss,
                   const std::vector<index_t>&  n1d_in,
                   bool                         newtonFailed,
                   const std::string&           reason)
{
    TensorRuleResult<T> result;
    result.nCandidateTensor = nCandidateTensor;
    result.nUsed            = nGauss;
    result.nGauss           = nGauss;
    result.n1d              = n1d_in;
    result.usedFallback     = true;
    result.newtonFailed     = newtonFailed;
    result.fallbackReason   = reason;

    result.pts.resize(d, nGauss);
    result.wts.resize(nGauss);

    gsVector<index_t> numPts(d);
    for (short_t dd = 0; dd < d; ++dd) numPts[dd] = pPerDir[dd] + 1;
    gsGaussRule<T> gaussRule(numPts);

    gsMatrix<T> elNodes;
    gsVector<T> elWeights;
    index_t k = 0;
    for (const ElementBox& e : cell.elements)
    {
        gaussRule.mapTo(e.lower, e.upper, elNodes, elWeights);
        const index_t nq = elWeights.size();
        result.pts.block(0, k, d, nq) = elNodes;
        result.wts.segment(k, nq)     = elWeights;
        k += nq;
    }
    return result;
}

/**
 * @brief Build a tensorized Johannessen quadrature rule for one atomic cell.
 *
 * For each direction dd, builds a local open target knot vector:
 *   - degree q_dd = 2 * p_dd   (p_dd = solution degree in direction dd)
 *   - regularity p_dd - 1      (internal multiplicity p_dd + 1)
 *   - breakpoints = unique leaf breakpoints from the level-l knot vector
 *
 * NOTE for future extension (non-simple knots):
 *   For solution knots with internal multiplicity m_sol > 1, the correct
 *   target regularity is (p_dd - m_sol), giving internal target multiplicity
 *   (p_dd + m_sol). Replace "p_dd - 1" below with "p_dd - m_sol" in that case.
 *
 * Decision policy:
 *   1. Build even target KVs and compute candidate counts without running Newton.
 *   2. Fall back to Gauss (reason="count")   if nCandidateTensor >= nGauss.
 *   3. Run 1D Newton per direction.
 *   4. Fall back to Gauss (reason="newton")  if Newton fails in any direction.
 *   5. Otherwise tensorize with a multi-index loop (direction 0 varies fastest).
 *
 * Works for any d >= 1.
 */
template<short_t d, class T>
TensorRuleResult<T> buildTensorRule(const AtomicCell&          cell,
                                     const gsHTensorBasis<d,T>& thb)
{
    // Per-direction solution degrees (anisotropic support)
    std::vector<index_t> pPerDir(d);
    for (short_t dd = 0; dd < d; ++dd) pPerDir[dd] = thb.degree(dd);

    // --- Step 1: build even target KVs and count candidates (no Newton yet) ---
    std::vector<gsKnotVector<T> > evenKv(d);
    std::vector<index_t>          n1d(d);

    for (short_t dd = 0; dd < d; ++dd)
    {
        const index_t p_dd = pPerDir[dd];
        const index_t q_dd = 2 * p_dd;

        const gsKnotVector<T>& kv_lev =
            thb.tensorLevel(cell.level).component(dd).knots();

        // Unique breakpoints inside the leaf (unique-knot access: kv_lev(i))
        std::vector<T> uKnots;
        for (index_t i = cell.lowerIndex[dd]; i <= cell.upperIndex[dd]; ++i)
            uKnots.push_back(kv_lev(i));

        // Open target KV: degree q_dd, C^{p_dd - 1} interior (mult p_dd + 1 at internals)
        // NOTE: for non-simple-knot case, replace (p_dd - 1) with (p_dd - m_sol).
        gsKnotVector<T> targetKv(uKnots, q_dd, p_dd - 1);
        evenKv[dd] = ensureEven1D(targetKv);

        const gsBSplineBasis<T> basisDD(evenKv[dd]);
        GISMO_ENSURE(basisDD.size() % 2 == 0,
            "Even target basis construction failed in direction " << dd);
        n1d[dd] = basisDD.size() / 2;
    }

    index_t nCandidateTensor = 1;
    for (short_t dd = 0; dd < d; ++dd) nCandidateTensor *= n1d[dd];

    // Standard Gauss: prod(p_dd + 1) points per element
    index_t ppow = 1;
    for (short_t dd = 0; dd < d; ++dd) ppow *= pPerDir[dd] + 1;
    const index_t nGauss = ppow * (index_t)cell.elements.size();

    // --- Step 2: fallback on count ---
    if (nCandidateTensor >= nGauss)
        return buildGaussFallback<d,T>(cell, pPerDir, nCandidateTensor, nGauss,
                                       n1d, false, "count");

    // --- Step 3: run Newton per direction ---
    std::vector<gsVector<T> > dir_nodes(d), dir_weights(d);
    for (short_t dd = 0; dd < d; ++dd)
    {
        const bool ok =
            johannessen1DFromEvenKv<T>(evenKv[dd], dir_nodes[dd], dir_weights[dd]);
        if (!ok)
            return buildGaussFallback<d,T>(cell, pPerDir, nCandidateTensor, nGauss,
                                           n1d, true, "newton");
    }

    // --- Step 4: tensorize with multi-index loop (direction 0 varies fastest) ---
    TensorRuleResult<T> result;
    result.nCandidateTensor = nCandidateTensor;
    result.nUsed            = nCandidateTensor;
    result.nGauss           = nGauss;
    result.n1d              = n1d;
    result.usedFallback     = false;
    result.newtonFailed     = false;
    result.fallbackReason   = "none";

    result.pts.resize(d, nCandidateTensor);
    result.wts.resize(nCandidateTensor);

    gsVector<index_t> idx(d);
    idx.setZero();
    for (index_t k = 0; k < nCandidateTensor; ++k)
    {
        T w = T(1);
        for (short_t dd = 0; dd < d; ++dd)
        {
            result.pts(dd, k) = dir_nodes[dd][idx[dd]];
            w *= dir_weights[dd][idx[dd]];
        }
        result.wts[k] = w;

        // Increment multi-index (direction 0 varies fastest)
        for (short_t dd = 0; dd < d; ++dd)
        {
            ++idx[dd];
            if (idx[dd] < n1d[dd]) break;
            idx[dd] = 0;
        }
    }

    return result;
}

// ============================================================
// Validation
// ============================================================

template <short_t d, class T>
struct ValidationResult
{
    index_t numPairs;
    T       maxAbsErr;
    T       sumSqErr;
    T       cellArea;
    TensorRuleResult<T> ruleResult;
};

template<short_t d, class T>
void printValidation(const ValidationResult<d,T>& valResult)
{
    gsInfo << "    Pairs (B_i B_j):  " << valResult.numPairs << "\n";
    gsInfo << "    Max |err|:        " << std::scientific << valResult.maxAbsErr << "\n";
    gsInfo << "    RMS |err|:        " << std::scientific
           << std::sqrt(valResult.sumSqErr / T(valResult.numPairs)) << "\n";
    gsInfo << "    Sum of weights:   " << std::fixed << std::setprecision(8)
           << valResult.ruleResult.wts.sum()
           << "  (cell area = " << valResult.cellArea << ")\n";
    gsInfo << "    Min/max weight:   "
           << valResult.ruleResult.wts.minCoeff() << " / " << valResult.ruleResult.wts.maxCoeff() << "\n";
    gsInfo << "    n1dX / n1dY:      " << valResult.ruleResult.n1d[0]
           << " / " << valResult.ruleResult.n1d[1] << "\n";
    gsInfo << "    Candidate tensor: " << valResult.ruleResult.nCandidateTensor << "\n";
    gsInfo << "    n_used / nGauss:  " << valResult.ruleResult.nUsed
           << " / " << valResult.ruleResult.nGauss << "  (ratio "
           << std::fixed << std::setprecision(3)
           << (valResult.ruleResult.nGauss > 0 ?
               (real_t)valResult.ruleResult.nUsed / valResult.ruleResult.nGauss : 0.0)
           << ")\n";
    gsInfo << "    Fallback reason:  " << valResult.ruleResult.fallbackReason << "\n";
}

/**
 * @brief Validate a quadrature rule against level tensor-basis products on the leaf.
 *
 * Computes reference moments  integral_{D_r} B_i * B_j  via high-order Gauss,
 * then compares with  sum_q w_q * B_i(x_q) * B_j(x_q)  from the rule.
 *
 * Uses the tensor basis at cell.level (not the global THB basis): active THB
 * functions restricted to a level-l leaf are representable as linear combinations
 * of level-l tensor B-splines. Therefore exactness for all level-l tensor-basis
 * product pairs implies exactness for THB mass-product restrictions on this leaf.
 *
 * Prints clearly: "Validation: level tensor-basis products on the leaf"
 */
template<short_t d, class T>
ValidationResult<d,T>
validateRule(const TensorRuleResult<T>& ruleResult,
                  const AtomicCell&          cell,
                  const gsHTensorBasis<d,T>& thb,
                  index_t overQuad = 4)
{
    // gsInfo << "\n  Validation: level tensor-basis products on the leaf\n";

    struct PairHash {
        size_t operator()(std::pair<index_t,index_t> p) const noexcept {
            return std::hash<size_t>()(
                (size_t(uint32_t(p.first)) << 32) | size_t(uint32_t(p.second)));
        }
    };
    std::unordered_map<std::pair<index_t,index_t>, index_t, PairHash> rowOfPair;
    rowOfPair.reserve(1024);
    std::vector<T> refMoments;

    // --- 1. Reference moments via high-order Gauss ---
    gsOptionList opts;
    opts.addInt   ("quRule", "Quadrature rule", gsQuadrature::GaussLegendre);
    opts.addReal  ("quA",    "quA",             T(1));
    opts.addInt   ("quB",    "quB",             (int)overQuad);
    opts.addSwitch("overInt","Over-integration", false);
    typename gsQuadRule<T>::uPtr gauss =
        gsQuadrature::getPtr(thb.tensorLevel(cell.level), opts);

    const gsTensorBSplineBasis<d,T>& tensor =
        dynamic_cast<const gsTensorBSplineBasis<d,T>&>(thb.tensorLevel(cell.level));

    gsFuncData<T> fd;
    fd.flags = NEED_VALUE | NEED_ACTIVE | SAME_ELEMENT;
    gsMatrix<T> nodes;
    gsVector<T> weights;

    for (const ElementBox& e : cell.elements)
    {
        gauss->mapTo(e.lower, e.upper, nodes, weights);
        tensor.compute(nodes, fd);

        const index_t na = fd.allActives().rows();
        const index_t nq = (index_t)weights.size();

        std::vector<index_t> rowCache;
        rowCache.reserve((na * (na + 1)) / 2);

        for (index_t a = 0; a < na; ++a)
        {
            const index_t bi = fd.allActives()(a, 0);
            for (index_t c = a; c < na; ++c)
            {
                const index_t bj = fd.allActives()(c, 0);
                const index_t ii = std::min(bi, bj), jj = std::max(bi, bj);
                auto key = std::make_pair(ii, jj);
                auto it  = rowOfPair.find(key);
                index_t row;
                if (it == rowOfPair.end())
                {
                    row = (index_t)rowOfPair.size();
                    rowOfPair[key] = row;
                    refMoments.push_back(T(0));
                }
                else
                    row = it->second;
                rowCache.push_back(row);
            }
        }

        for (index_t q = 0; q < nq; ++q)
        {
            const T w = weights[q];
            index_t k = 0;
            for (index_t a = 0; a < na; ++a)
                for (index_t c = a; c < na; ++c, ++k)
                    refMoments[rowCache[k]] +=
                        w * fd.allValues()(a, q) * fd.allValues()(c, q);
        }
    }

    // --- 2. Rule moments ---
    const index_t numPairs = (index_t)rowOfPair.size();
    std::vector<T> ruleMoments(numPairs, T(0));

    // Rule pts span multiple elements, so evaluate without SAME_ELEMENT
    gsFuncData<T> fd2;
    fd2.flags = NEED_VALUE | NEED_ACTIVE;
    tensor.compute(ruleResult.pts, fd2);

    const index_t na2 = fd2.allActives().rows();
    for (index_t q = 0; q < (index_t)ruleResult.wts.size(); ++q)
    {
        const T wq = ruleResult.wts[q];
        for (index_t a = 0; a < na2; ++a)
        {
            const index_t bi = fd2.allActives()(a, q);
            for (index_t c = a; c < na2; ++c)
            {
                const index_t bj = fd2.allActives()(c, q);
                const index_t ii = std::min(bi, bj), jj = std::max(bi, bj);
                auto it = rowOfPair.find(std::make_pair(ii, jj));
                if (it != rowOfPair.end())
                    ruleMoments[it->second] +=
                        wq * fd2.allValues()(a, q) * fd2.allValues()(c, q);
            }
        }
    }

    // --- 3. Error report ---
    T maxAbsErr = T(0), sumSqErr = T(0);
    for (index_t r = 0; r < numPairs; ++r)
    {
        const T err = std::abs(ruleMoments[r] - refMoments[r]);
        maxAbsErr = std::max(maxAbsErr, err);
        sumSqErr += err * err;
    }

    const T cellArea =
        (cell.upperParam[0] - cell.lowerParam[0]) *
        (cell.upperParam[1] - cell.lowerParam[1]);

    ValidationResult<d,T> valResult;
    valResult.numPairs = numPairs;
    valResult.maxAbsErr = maxAbsErr;
    valResult.sumSqErr = sumSqErr;
    valResult.cellArea = cellArea;
    valResult.ruleResult = ruleResult;
    return valResult;
}

// ============================================================
// main
// ============================================================

template<short_t d, class T>
void run_main(bool plot, index_t numLevels, index_t overQuad)
{
    // ------------------------------------------------------------------
    // Build THB-spline basis with L-shaped local refinement
    // (same geometry as hierarchical_macro_mass_rule_single.cpp)
    // ------------------------------------------------------------------
    gsKnotVector<T> kv(0, 1.0, 0, 3, 1);
    kv.uniformRefine();
    kv.uniformRefine();

    std::vector<gsKnotVector<T> > kvs;
    for (short_t dd = 0; dd < d; ++dd) kvs.push_back(kv);

    gsTensorBSplineBasis<d,T> tbsb(kvs);
    gsTHBSplineBasis<d,T>     thb(tbsb);

    {
        thb.uniformRefine();
        for (index_t lvl = 1; lvl <= numLevels; ++lvl)
        {
            const index_t N = 4 * (1 << lvl);
            std::vector<index_t> boxes;
            for (short_t DD = 0; DD < d; ++DD)
            {
                boxes.push_back(lvl);      
                for (short_t dd = 0; dd < d; ++dd)
                    boxes.push_back(0);

                for (short_t dd = 0; dd < d; ++dd)
                    if (dd == DD)   
                        boxes.push_back(N);
                    else
                         boxes.push_back(4);
            }
            thb.refineElements(boxes);
        }
    }

    const index_t deg = thb.degree(0);
    gsInfo << "Basis: degree=" << deg
           << ", maxLevel=" << thb.maxLevel()
           << ", DOF=" << thb.size() << "\n\n";

    if (plot)
    {
        gsMesh<T> mesh(thb);
        gsWriteParaview(mesh, "hierarchical_tensor_mesh");
    }

    // ------------------------------------------------------------------
    // Collect atomic cells
    // ------------------------------------------------------------------
    const std::vector<AtomicCell> cells = collectAtomicCells<d,T>(thb);
    gsInfo << "Atomic cells: " << cells.size() << "\n\n";

    // ------------------------------------------------------------------
    // Build tensor rules for all cells
    // ------------------------------------------------------------------
    std::vector<TensorRuleResult<T>> results;
    results.reserve(cells.size());
    for (index_t r = 0; r < (index_t)cells.size(); ++r)
        results.push_back(buildTensorRule<d,T>(cells[r], thb));

    // ------------------------------------------------------------------
    // Summary table
    // ------------------------------------------------------------------
    gsInfo << std::setw(6)  << "cell"
           << std::setw(6)  << "lev"
           << std::setw(7)  << "elems"
           << std::setw(7)  << "n1dX"
           << std::setw(7)  << "n1dY"
           << std::setw(11) << "candTensor"
           << std::setw(9)  << "used"
           << std::setw(9)  << "nGauss"
           << std::setw(9)  << "ratio"
           << std::setw(10) << "fallback"
           << "\n";
    gsInfo << std::string(81, '-') << "\n";

    index_t sumUsed = 0, sumGauss = 0;
    index_t nFallbackCount = 0, nFallbackNewton = 0;

    for (index_t r = 0; r < (index_t)cells.size(); ++r)
    {
        const AtomicCell&               cell = cells[r];
        const TensorRuleResult<T>& res  = results[r];

        gsInfo << std::setw(6) << r
               << std::setw(6) << cell.level
               << std::setw(7) << (index_t)cell.elements.size()
               << std::setw(7) << res.n1d[0]
               << std::setw(7) << res.n1d[1]
               << std::setw(11) << res.nCandidateTensor
               << std::setw(9) << res.nUsed
               << std::setw(9) << res.nGauss
               << std::setw(9) << std::fixed << std::setprecision(3)
                               << (res.nGauss > 0 ?
                                   (real_t)res.nUsed / res.nGauss : 0.0)
               << std::setw(10) << res.fallbackReason
               << "\n";

        sumUsed  += res.nUsed;
        sumGauss += res.nGauss;
        if (res.fallbackReason == "count")  ++nFallbackCount;
        if (res.fallbackReason == "newton") ++nFallbackNewton;
    }

    // gsInfo << "\nGlobal summary\n==============\n";
    // gsInfo << "  Total cells:        " << cells.size()      << "\n";
    // gsInfo << "  Fallback (count):   " << nFallbackCount    << "\n";
    // gsInfo << "  Fallback (newton):  " << nFallbackNewton   << "\n";
    // gsInfo << "  Sum used:           " << sumUsed           << "\n";
    // gsInfo << "  Sum nGauss:         " << sumGauss          << "\n";
    // gsInfo << "  Global ratio:       " << std::fixed << std::setprecision(4)
    //        << (sumGauss > 0 ? (real_t)sumUsed / sumGauss : 0.0) << "\n\n";

    // if ((real_t)sumUsed / sumGauss < 1.0)
    //     gsInfo << ">>> Viability: PROMISING — used pts < Gauss pts globally.\n\n";
    // else
    //     gsInfo << ">>> Viability: NOT PROMISING — used pts >= Gauss pts globally.\n\n";


    // ------------------------------------------------------------------
    // Detailed validation
    // ------------------------------------------------------------------
    const index_t first = 0;
    const index_t last  = (index_t)cells.size();
    gsInfo<<"|------------|-------|-------|-----------|----------|------------|----------|----------|----------------|\n";
    gsInfo<<"| Cell index | Level | Elems | maxAbsErr | sumSqErr | sumWeights | cellArea | N/NGauss | fallbackReason |\n";
    gsInfo<<"|------------|-------|-------|-----------|----------|------------|----------|----------|----------------|\n";
    for (index_t r = first; r < last; ++r)
    {
        const AtomicCell& cell = cells[r];
        index_t overQ = 4;
        auto result = validateRule<d,T>(results[r], cell, thb, overQ);
        gsInfo<<"| " << std::setw(10) << r
              << " | " << std::setw(5)  << cell.level
              << " | " << std::setw(5)  << (index_t)cell.elements.size()
              << " | " << std::setw(10) << std::scientific << result.maxAbsErr
              << " | " << std::setw(8) << std::scientific << result.sumSqErr
              << " | " << std::setw(10) << std::fixed << std::setprecision(8) << result.ruleResult.wts.sum()
              << " | " << std::setw(9) << std::fixed << std::setprecision(8) << result.cellArea
              << " | " << std::setw(8) << std::fixed << std::setprecision(3)
                              << (result.ruleResult.nGauss > 0 ?
                                  (real_t)result.ruleResult.nUsed / result.ruleResult.nGauss : 0.0)
              << " | " << std::setw(14) << result.ruleResult.fallbackReason
              << "\n";
        // gsInfo << "\nCell " << r
        //        << "  level=" << cell.level
        //        << "  elems=" << cell.elements.size()
        //        << "  [" << cell.lowerParam.transpose()
        //        << " , " << cell.upperParam.transpose() << "]\n";
        // auto result = validateRule<d,T>(results[r], cell, thb, overQ);
        // printValidation(result);
    }
    gsInfo<<"|------------|-------|-------|-----------|----------|------------|----------|----------|----------------|\n";
    gsInfo<<"| " << std::setw(10) << "TOTAL "
        << " | " << std::setw(5)  << "-"
        << " | " << std::setw(5)  << (index_t)cells.size()
        << " | " << std::setw(10) << std::scientific << "-"
        << " | " << std::setw(8) << std::scientific << "-"
        << " | " << std::setw(10) << std::fixed << std::setprecision(8) << "-"
        << " | " << std::setw(9) << std::fixed << std::setprecision(8) << "-"
        << " | " << std::setw(8) << std::fixed << std::setprecision(3)
                        << (sumGauss > 0 ?
                            (real_t)sumUsed / sumGauss : 0.0)
        << " | " << std::setw(14) << "-"
        << "\n";

}

int main(int argc, char* argv[])
{
    // index_t overQ  = 4;
    bool    plot   = false;
    short_t d      = 2;
    index_t numLevels = 5;

    gsCmdLine cmd("Tensorized Johannessen quadrature rule for hierarchical domains.");
    // cmd.addInt   ("o", "overQ",  "Extra Gauss pts per direction (reference)",     overQ);
    cmd.addInt   ("d", "dim",    "Problem dimension (default: 2)",                d);
    cmd.addInt   ("l", "levels", "Number of levels (default: 5)",                 numLevels);
    cmd.addSwitch("plot",        "Write ParaView mesh file",                       plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (d==2)
        run_main<2, real_t>(plot, numLevels, /*overQ=*/4);
    else if (d==3)
        run_main<3, real_t>(plot, numLevels, /*overQ=*/4);
    else if (d==4)
        run_main<4, real_t>(plot, numLevels, /*overQ=*/4);
    else
    {
        gsInfo << "Dimension " << d << " not supported.\n";
        return 1;
    }
    return 0;
}
