/** @file hierarchical_patch_rule.cpp

    @brief Hierarchical patch rule for THB-splines.

    Core idea (derived from discussion on nested quadrature, May 2026):

      1. Compute ONE shared set of quadrature points from the finest-level
         Johannessen patch rule (target space: degree 2p, regularity p-1).

      2. For each level l = 0..L, solve the minimum-norm linear system
             N^l * w^l = I^l
         to obtain level-l directional weights at those shared points.
         This system always has a solution because S^p_{tau^l} ⊂ S^p_{tau^L}
         (knot refinement), so the finest-level points are sufficient for
         every coarser level's target space.

      3. Assembly uses *rep-level dispatch*: for a basis-function pair (i,j)
         with max representationLevel = Lp, use the shared points with the
         Lp-specific weight vector.

    Benefit: geometry map, Jacobian, and material-law evaluations are performed
    ONCE per quadrature point per element, regardless of how many levels have
    active basis functions there.  Only the scalar weight changes per level.

    Relationship to existing Cases:
      - Case 1  (existing file): finest-level patch rule for all elements.
        Equivalent to our rule with Lp = maxLevel for every pair — correct
        by the inclusion argument, but slightly over-integrates coarser pairs.
      - Case 10 (existing file): rep-level dispatch with SEPARATE (non-nested)
        gsPatchRule objects per level — requires multiple basis evaluations per
        element.  This file replaces those separate rules with a single shared
        point set and level-specific weights.

    Author: H.M. Verhelst
*/

#include <gismo.h>
#include <gsAssembler/gsPatchRule.h>

using namespace gismo;

// =============================================================================
// Free-function helpers
// =============================================================================

/**
 * @brief Replicate gsPatchRule::_init to obtain the target integration knot
 *        vector for a given physical knot vector.
 *
 * Matches gsPatchRule::_init behaviour exactly (including the known second-
 * branch typo — only rdiff>0 is applied).
 *
 * @param kv      Physical knot vector of one parametric direction
 * @param intDeg  Target integration degree  (e.g. 2*p for mass matrix)
 * @param reg     Target regularity          (e.g. p-1 for max-cont splines)
 * @param overInt Over-integrate boundary elements
 */
template<class T>
gsKnotVector<T> makeTargetKv(const gsKnotVector<T>& kv,
                              index_t intDeg, index_t reg,
                              bool overInt = false)
{
    gsKnotVector<T> result = kv;

    // Compute diffs from the ORIGINAL knot vector (before any modification),
    // exactly as in gsPatchRule::_init.
    const index_t pdiff = intDeg - result.degree();
    std::vector<index_t> mults = result.multiplicities();
    const index_t rmin  = *std::min_element(mults.begin(), mults.end());
    const index_t rdiff = (intDeg - reg) - rmin;

    // Apply degree change
    if      (pdiff > 0) result.degreeIncrease(pdiff);
    else if (pdiff < 0) result.degreeDecrease(-pdiff);

    // Apply regularity change (increase only — matching gsPatchRule behaviour)
    if (rdiff > 0)
        result.increaseMultiplicity(rdiff);

    index_t size = result.size() - result.degree() - 1;

    // Over-integration: insert (deg-1) knots in each boundary element
    if (overInt && result.degree() > 1)
    {
        const index_t numOver = result.degree() - 1;
        const T lLen = (result[1]              - result.first()) / T(numOver + 1);
        const T uLen = (result.last() - result[result.uSize()-2]) / T(numOver + 1);
        for (index_t k = 0; k < numOver; ++k)
        {
            result.insert(result.first() + T(k+1) * lLen);
            result.insert(result.last()  - T(k+1) * uLen);
        }
        size += 2 * numOver;
    }

    // Make size even: insert one knot at the centre of the largest knot span
    if (size % 2 == 1)
    {
        T maxSpan = T(0), mid = T(0);
        for (index_t i = 0; i + 1 < result.size(); ++i)
        {
            const T span = result[i+1] - result[i];
            if (span > maxSpan)
            {
                maxSpan = span;
                mid     = (result[i] + result[i+1]) / T(2);
            }
        }
        result.insert(mid);
    }

    return result;
}

/**
 * @brief Compute 1-D quadrature weights at the given fixed points @p pts for
 *        the spline target space defined by @p kvTarget.
 *
 * Solves the minimum-norm (underdetermined) or least-squares (overdetermined)
 * linear system
 *         N^l * w^l = I^l
 * where
 *   N^l[j,k] = N_j(pts[k])   (target-space basis evaluated at shared points)
 *   I^l[j]   = (kv[j+p+1] - kv[j]) / (p+1)   (exact integrals, de Boor formula)
 *
 * For l < L (finest): system is underdetermined  — returns min-norm solution.
 * For l = L (finest): system is overdetermined   — returns least-squares sol.
 *   (At the Johannessen points the overdetermined system has an exact solution;
 *    the SVD recovers it when the residual is at machine precision.)
 *
 * @param kvTarget  Target integration knot vector for level l
 * @param pts       Shared finest-level 1-D quadrature points
 */
template<class T>
gsVector<T> computeLevelWeights(const gsKnotVector<T>& kvTarget,
                                 const gsVector<T>& pts)
{
    gsBSplineBasis<T> basis(kvTarget);
    const index_t n = basis.size();   // number of equations  (basis fns at level l)
    const index_t K = pts.size();     // number of unknowns   (shared finest pts)
    const index_t p = kvTarget.degree();

    // Build N^l ∈ ℝ^{n × K}
    gsMatrix<T> N = gsMatrix<T>::Zero(n, K);
    gsMatrix<T>      ptsRow = pts.transpose();   // 1 × K
    gsMatrix<index_t> acts;
    gsMatrix<T>       vals;
    basis.active_into(ptsRow, acts);
    basis.eval_into  (ptsRow, vals);
    for (index_t k = 0; k < K; ++k)
        for (index_t r = 0; r < acts.rows(); ++r)
            N(acts(r,k), k) = vals(r,k);

    // Build exact integral vector I^l (de Boor formula, eq. (6) of Johannessen 2017)
    gsVector<T> I(n);
    for (index_t j = 0; j < n; ++j)
        I(j) = (kvTarget[j + p + 1] - kvTarget[j]) / T(p + 1);

    // Min-norm / least-squares solve via thin SVD
    return N.bdcSvd(gsEigen::ComputeThinU | gsEigen::ComputeThinV).solve(I);
}

// =============================================================================
// gsHierarchicalPatchRule
// =============================================================================

/**
 * @brief Hierarchical / embedded patch rule for 2-D THB-splines.
 *
 * This version is deliberately different from the previous "finest points
 * everywhere" experiment.
 *
 * For every evaluation level m = 0..L it stores the ordinary level-m
 * Johannessen point set X_m.  For every target level l <= m it computes a
 * weight vector w^{l->m} on X_m by solving
 *
 *      N^l(X_m) w^{l->m} = I^l.
 *
 * Hence an element can be evaluated at a local maximum level m only once, and
 * all lower-level contributions on that same element can reuse those function
 * evaluations with different weights.  The point sets X_m are not forced to be
 * geometrically nested; the nestedness here is in the embedding of lower-level
 * target spaces onto higher-level evaluation point sets.
 */
template<class T>
class gsHierarchicalPatchRule
{
public:
    /**
     * @param basis    THB-spline basis (provides per-level knot vectors)
     * @param intDeg   Integration degree  (e.g. 2*p for mass)
     * @param reg      Integration regularity
     * @param overInt  Over-integrate boundary elements
     */
    gsHierarchicalPatchRule(const gsHTensorBasis<2,T>& basis,
                             index_t intDeg, index_t reg,
                             bool overInt = false)
        : m_L(basis.maxLevel())
    {
        m_end = basis.support().col(1);

        m_nodes.resize(2);
        m_nodeMap.resize(2);
        m_levelW.resize(2);

        gsOptionList opts;
        opts.addInt   ("quRule", "", gsQuadrature::PatchRule);
        opts.addReal  ("quA",    "", T(intDeg));
        opts.addInt   ("quB",    "", reg);
        opts.addSwitch("overInt","", overInt);

        for (short_t d = 0; d < 2; ++d)
        {
            m_nodes[d].resize(m_L + 1);
            m_nodeMap[d].resize(m_L + 1);
            m_levelW[d].resize(m_L + 1); // target level
            for (index_t l = 0; l <= m_L; ++l)
                m_levelW[d][l].resize(m_L + 1); // evaluation level
        }

        // Build one ordinary patch point set X_m for every level m.
        for (index_t m = 0; m <= m_L; ++m)
        {
            auto ruleM = gsQuadrature::getPtr(basis.tensorLevel(m), opts);
            const gsPatchRule<T>* prM =
                dynamic_cast<const gsPatchRule<T>*>(ruleM.get());
            GISMO_ENSURE(prM,
                "gsHierarchicalPatchRule: gsQuadrature did not produce a gsPatchRule. "
                "Check that PatchRule is registered for the given basis type.");

            for (short_t d = 0; d < 2; ++d)
            {
                m_nodes[d][m] = prM->nodes(d);
                for (index_t k = 0; k < m_nodes[d][m].size(); ++k)
                    m_nodeMap[d][m][m_nodes[d][m][k]] = k;
            }
        }

        // For every target level l and every evaluation level m >= l, compute
        // weights w^{l->m} on the level-m point set.
        for (short_t d = 0; d < 2; ++d)
        {
            for (index_t l = 0; l <= m_L; ++l)
            {
                const gsKnotVector<T>* kvPhys =
                    static_cast<const gsKnotVector<T>*>(&basis.tensorLevel(l).knots(d));
                gsKnotVector<T> kvT = makeTargetKv(*kvPhys, intDeg, reg, overInt);

                for (index_t m = l; m <= m_L; ++m)
                    m_levelW[d][l][m] = computeLevelWeights(kvT, m_nodes[d][m]);
            }
        }
    }

    /**
     * @brief Quadrature nodes and weights on an element.
     *
     * @param targetLevel  Level of the target space to integrate.
     * @param evalLevel    Level of the point set used for evaluation. Must be
     *                     >= targetLevel.
     */
    void mapTo(const gsVector<T>& lower, const gsVector<T>& upper,
               index_t targetLevel, index_t evalLevel,
               gsMatrix<T>& nodes, gsVector<T>& weights) const
    {
        GISMO_ASSERT(targetLevel >= 0 && targetLevel <= m_L,
            "targetLevel " << targetLevel << " out of range [0," << m_L << "]");
        GISMO_ASSERT(evalLevel >= 0 && evalLevel <= m_L,
            "evalLevel " << evalLevel << " out of range [0," << m_L << "]");
        GISMO_ASSERT(evalLevel >= targetLevel,
            "Need evalLevel >= targetLevel for embedded weights.");

        std::vector<std::vector<index_t>> idxs(2);
        for (short_t d = 0; d < 2; ++d)
        {
            const bool isLast = (upper[d] >= m_end[d] - T(1e-14));
            const std::map<T,index_t>& mp = m_nodeMap[d][evalLevel];
            auto itLo = mp.lower_bound(lower[d]);
            auto itHi = isLast ? mp.end() : mp.lower_bound(upper[d]);
            for (auto it = itLo; it != itHi; ++it)
                idxs[d].push_back(it->second);
        }
        _tensorProduct(idxs, targetLevel, evalLevel, nodes, weights);
    }

    /// Convenience overload: target and evaluation levels coincide.
    void mapTo(const gsVector<T>& lower, const gsVector<T>& upper,
               index_t level,
               gsMatrix<T>& nodes, gsVector<T>& weights) const
    {
        mapTo(lower, upper, level, level, nodes, weights);
    }

    /**
     * @brief Only return weights for element [lower, upper] with target level
     *        @p targetLevel on evaluation point set @p evalLevel.
     */
    void getWeights(const gsVector<T>& lower, const gsVector<T>& upper,
                    index_t targetLevel, index_t evalLevel,
                    gsVector<T>& weights) const
    {
        gsMatrix<T> nodes;
        mapTo(lower, upper, targetLevel, evalLevel, nodes, weights);
    }

    /// Convenience overload: target and evaluation levels coincide.
    void getWeights(const gsVector<T>& lower, const gsVector<T>& upper,
                    index_t level, gsVector<T>& weights) const
    {
        getWeights(lower, upper, level, level, weights);
    }

    index_t maxLevel() const { return m_L; }

    const gsVector<T>& nodes(short_t d, index_t evalLevel) const
    {
        GISMO_ASSERT(d >= 0 && d < 2, "Direction out of range.");
        GISMO_ASSERT(evalLevel >= 0 && evalLevel <= m_L, "Evaluation level out of range.");
        return m_nodes[d][evalLevel];
    }

private:
    // Build tensor-product nodes and weights from per-direction index sets.
    void _tensorProduct(const std::vector<std::vector<index_t>>& idxs,
                        index_t targetLevel, index_t evalLevel,
                        gsMatrix<T>& nodes, gsVector<T>& weights) const
    {
        const index_t nx = static_cast<index_t>(idxs[0].size());
        const index_t ny = static_cast<index_t>(idxs[1].size());
        const index_t sz = nx * ny;

        nodes.resize(2, sz);
        weights.resize(sz);

        if (sz == 0) return;

        const gsVector<T>& wx = m_levelW[0][targetLevel][evalLevel];
        const gsVector<T>& wy = m_levelW[1][targetLevel][evalLevel];
        const gsVector<T>& x  = m_nodes[0][evalLevel];
        const gsVector<T>& y  = m_nodes[1][evalLevel];

        index_t cnt = 0;
        for (index_t jy = 0; jy < ny; ++jy)
        {
            const index_t ky = idxs[1][jy];
            for (index_t jx = 0; jx < nx; ++jx, ++cnt)
            {
                const index_t kx = idxs[0][jx];
                nodes(0, cnt) = x[kx];
                nodes(1, cnt) = y[ky];
                weights[cnt]  = wx[kx] * wy[ky];
            }
        }
    }

    index_t     m_L;
    gsVector<T> m_end;   // parametric domain upper corner

    // m_nodes[d][m]: ordinary level-m 1-D patch nodes in direction d
    std::vector<std::vector<gsVector<T>>>          m_nodes;
    std::vector<std::vector<std::map<T,index_t>>>  m_nodeMap;

    // m_levelW[d][l][m]: weights for target level l on evaluation points X_m
    std::vector<std::vector<std::vector<gsVector<T>>>> m_levelW;
};

// =============================================================================
// main
// =============================================================================

int main(int argc, char* argv[])
{
    index_t order      = 2;
    index_t regularity = 1;
    bool    plot       = false;

    gsCmdLine cmd("Hierarchical patch rule: shared points, per-level weights.");
    cmd.addInt   ("P","deg","Polynomial order of basis",  order);
    cmd.addInt   ("R","reg","Regularity of basis",        regularity);
    cmd.addSwitch("plot","Export mesh to Paraview",       plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(regularity < order, "Require regularity < order.");

    // ── Build THB-spline basis ─────────────────────────────────────────────
    gsKnotVector<> kv(0.0, 1.0, 4, order+1, 1);  // 4 elems, C^{order-1}
    gsTensorBSplineBasis<2,real_t> tbsb(kv, kv);
    gsTHBSplineBasis<2,real_t,true> thb(tbsb);

    thb.uniformRefine();
    thb.refineElements({1, 0, 0, 4, 4});   // lower-left quarter at level 1
    thb.refineElements({2, 0, 0, 4, 4});   // finer corner at level 2

    const index_t L = thb.maxLevel();
    gsInfo << "THB basis: " << thb.size() << " functions, maxLevel=" << L << "\n";
    {
        std::map<index_t,index_t> cnt;
        for (index_t i = 0; i != thb.size(); ++i) ++cnt[thb.levelOf(i)];
        gsInfo << "  per-level counts:";
        for (auto& p : cnt) gsInfo << "  L" << p.first << ":" << p.second;
        gsInfo << "\n";
    }

    if (plot) { gsMesh<> m(thb); gsWriteParaview(m,"hrule_mesh"); }

    // ── Integration target degrees ─────────────────────────────────────────
    // Mass matrix: B_i*B_j uses the univariate target S^{2p}_{r}.
    // Plane stiffness: grad(B_i)·grad(B_j) contains mixed tensor targets:
    //   d/dxi contribution  : S^{2p-2}_{r-1} in xi  × S^{2p}_{r} in eta
    //   d/deta contribution : S^{2p}_{r} in xi      × S^{2p-2}_{r-1} in eta
    // This example uses the safe tensor-product superspace S^{2p}_{r-1} in
    // both directions, so that the stiffness test is not under-integrated.
    const index_t intDegM = 2 * order;
    const index_t intRegM = regularity;
    const index_t intDegK = 2 * order;
    const index_t intRegK = regularity - 1;

    // ── Reference: Gauss-Legendre on finest-level mesh ─────────────────────
    gsOptionList gOpts;
    gOpts.addInt   ("quRule","", gsQuadrature::GaussLegendre);
    gOpts.addReal  ("quA",   "", 1.0);
    gOpts.addInt   ("quB",   "", 1);
    gOpts.addSwitch("overInt","",false);
    auto gaussRule = gsQuadrature::getPtr(thb.tensorLevel(L), gOpts);

    gsVector<real_t> Iref(thb.size()); Iref.setZero();
    gsMatrix<real_t> Mref = gsMatrix<real_t>::Zero(thb.size(), thb.size());
    gsMatrix<real_t> Kref = gsMatrix<real_t>::Zero(thb.size(), thb.size());
    index_t gaussTotalPts = 0;
    {
        gsMatrix<real_t> nodes; gsVector<real_t> weights;
        gsMatrix<index_t> acts; gsMatrix<real_t> vals, derivs;
        for (auto& elem : thb.domain()->allElements())
        {
            gaussRule->mapTo(elem.lowerCorner(), elem.upperCorner(), nodes, weights);
            gaussTotalPts += weights.size();
            thb.active_into(nodes, acts);
            thb.eval_into  (nodes, vals);
            thb.deriv_into (nodes, derivs);
            for (index_t k = 0; k < weights.size(); ++k)
                for (index_t r = 0; r < acts.rows(); ++r)
                {
                    Iref[acts(r,0)] += weights[k] * vals(r,k);
                    for (index_t c = r; c < acts.rows(); ++c)
                    {
                        const index_t i = acts(r,0), j = acts(c,0);
                        const real_t sm = weights[k] * vals(r,k) * vals(c,k);
                        Mref(i,j) += sm; if (i!=j) Mref(j,i) += sm;
                        const real_t sk = weights[k] *
                            (derivs(2*r+0,k)*derivs(2*c+0,k) +
                             derivs(2*r+1,k)*derivs(2*c+1,k));
                        Kref(i,j) += sk; if (i!=j) Kref(j,i) += sk;
                    }
                }
        }
    }

    // ── Hierarchical patch rules (mass and stiffness) ──────────────────────
    gsHierarchicalPatchRule<real_t> hprM(thb, intDegM, intRegM);
    gsHierarchicalPatchRule<real_t> hprK(thb, intDegK, intRegK);

    // ── Assembly with local evaluation-level dispatch ───────────────────────
    //
    // For each leaf element we choose evalLevel = max repLevel among the local
    // active pair set.  The basis is evaluated once on the level-evalLevel
    // point set.  Contributions with target level Lp <= evalLevel reuse these
    // evaluations with weights w^{Lp->evalLevel}.
    //
    // This is an embedded-weight experiment.  It is cheaper than evaluating
    // unrelated level-wise point clouds on the same element, but it is not the
    // same as a globally exact hierarchical moment-fitting rule.
    //
    gsVector<real_t> IA(thb.size()); IA.setZero();
    gsMatrix<real_t> MA = gsMatrix<real_t>::Zero(thb.size(), thb.size());
    gsMatrix<real_t> KA = gsMatrix<real_t>::Zero(thb.size(), thb.size());
    index_t hprMassPts  = 0;
    index_t hprStiffPts = 0;
    {
        gsMatrix<real_t> nodesM, nodesK;
        gsVector<real_t> wEval, wLp, wLpK;
        gsMatrix<index_t> acts, actsK;
        gsMatrix<real_t>  vals, derivs;
        gsMatrix<real_t>  ctr(2,1);
        gsMatrix<index_t> actCtr;

        for (auto& elem : thb.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,2>*>(&elem);
            GISMO_ASSERT(helem, "Element iterator is not gsHDomainIterator.");
            const auto& low = helem->lowerCorner();
            const auto& upp = helem->upperCorner();

            // Determine active basis functions via element centre.
            ctr.col(0) = 0.5 * (low + upp);
            thb.active_into(ctr, actCtr);
            if (actCtr.rows() == 0) continue;

            // Distinct target levels for local matrix entries, and the local
            // evaluation level used for the shared point set.
            std::set<index_t> Lp_set;
            index_t evalLevel = 0;
            for (index_t r = 0; r < actCtr.rows(); ++r)
            for (index_t c = r; c < actCtr.rows(); ++c)
            {
                const index_t Lp = std::max(thb.repLevel(actCtr(r,0)),
                                            thb.repLevel(actCtr(c,0)));
                Lp_set.insert(Lp);
                evalLevel = std::max(evalLevel, Lp);
            }

            // Shared mass-rule nodes at the local evaluation level.
            hprM.mapTo(low, upp, evalLevel, evalLevel, nodesM, wEval);
            if (wEval.size() == 0) continue;
            hprMassPts += wEval.size();

            thb.active_into(nodesM, acts);
            thb.eval_into  (nodesM, vals);

            // Shared stiffness-rule nodes at the same local evaluation level.
            // hprK may have a different target space and therefore a different
            // point set from hprM.
            hprK.mapTo(low, upp, evalLevel, evalLevel, nodesK, wEval);
            if (nodesK.cols() > 0)
            {
                hprStiffPts += wEval.size();
                thb.active_into(nodesK, actsK);
                thb.deriv_into (nodesK, derivs);
            }

            // Apply target-level-specific weights on the local evalLevel nodes.
            for (index_t Lp : Lp_set)
            {
                // Mass weights w^{Lp->evalLevel}.
                hprM.getWeights(low, upp, Lp, evalLevel, wLp);
                if (wLp.size() == 0) continue;

                // Per-basis integrals: dispatch by the function's own repLevel.
                for (index_t r = 0; r < acts.rows(); ++r)
                {
                    const index_t i = acts(r,0);
                    if (thb.repLevel(i) != Lp) continue;
                    for (index_t k = 0; k < wLp.size(); ++k)
                        IA[i] += wLp[k] * vals(r,k);
                }

                // Mass matrix pairs with max(repLevel) == Lp.
                for (index_t r = 0; r < acts.rows(); ++r)
                for (index_t c = r; c < acts.rows(); ++c)
                {
                    const index_t i = acts(r,0), j = acts(c,0);
                    if (std::max(thb.repLevel(i), thb.repLevel(j)) != Lp) continue;
                    real_t sm = 0;
                    for (index_t k = 0; k < wLp.size(); ++k)
                        sm += wLp[k] * vals(r,k) * vals(c,k);
                    MA(i,j) += sm; if (i != j) MA(j,i) += sm;
                }

                // Stiffness matrix pairs with max(repLevel) == Lp.
                if (nodesK.cols() > 0)
                {
                    hprK.getWeights(low, upp, Lp, evalLevel, wLpK);
                    for (index_t r = 0; r < actsK.rows(); ++r)
                    for (index_t c = r; c < actsK.rows(); ++c)
                    {
                        const index_t i = actsK(r,0), j = actsK(c,0);
                        if (std::max(thb.repLevel(i), thb.repLevel(j)) != Lp) continue;
                        real_t sk = 0;
                        for (index_t k = 0; k < wLpK.size(); ++k)
                            sk += wLpK[k] * (derivs(2*r+0,k)*derivs(2*c+0,k) +
                                             derivs(2*r+1,k)*derivs(2*c+1,k));
                        KA(i,j) += sk; if (i != j) KA(j,i) += sk;
                    }
                }
            }
        }
    }

    // ── Verification: finest-level weights already integrate all levels ─────
    //
    // Since S^p_{tau^l} ⊂ S^p_{tau^L}, the finest-level weights w^L give
    // exact integration for ALL levels' basis functions — no level dispatch
    // needed for scalar integrals.  This verifies the theoretical claim.
    //
    gsVector<real_t> I_finest(thb.size()); I_finest.setZero();
    {
        gsMatrix<real_t> nodes; gsVector<real_t> weights;
        gsMatrix<index_t> acts; gsMatrix<real_t> vals;
        for (auto& elem : thb.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,2>*>(&elem);
            hprM.mapTo(helem->lowerCorner(), helem->upperCorner(), L,
                       nodes, weights);
            if (weights.size() == 0) continue;
            thb.active_into(nodes, acts);
            thb.eval_into  (nodes, vals);
            for (index_t r = 0; r < acts.rows(); ++r)
                for (index_t k = 0; k < weights.size(); ++k)
                    I_finest[acts(r,0)] += weights[k] * vals(r,k);
        }
    }

    // ── Results ────────────────────────────────────────────────────────────
    gsInfo << "\n========== Hierarchical Patch Rule Results ==========\n";

    gsInfo << "\n--- Scalar integrals ---\n";
    gsInfo << "||I_hierarchical - I_ref||_inf = "
           << (IA - Iref).cwiseAbs().maxCoeff() << "\n";
    gsInfo << "||I_finest_weights - I_ref||_inf = "
           << (I_finest - Iref).cwiseAbs().maxCoeff()
           << "  (should be ≈ machine eps: finest weights suffice for all levels)\n";
    {
        real_t eU = 0, eT = 0;
        for (index_t i = 0; i < thb.size(); ++i)
        {
            const real_t e = std::abs(IA[i] - Iref[i]);
            if (thb.isTruncated(i)) eT = std::max(eT, e);
            else                    eU = std::max(eU, e);
        }
        gsInfo << "  untruncated=" << eU << "  truncated=" << eT << "\n";
    }

    gsInfo << "\n--- Mass matrix ---\n";
    gsInfo << "||M_hierarchical - M_ref||_inf = "
           << (MA - Mref).cwiseAbs().maxCoeff() << "\n";
    {
        real_t eUU=0, eUT=0, eTT=0;
        for (index_t i = 0; i < thb.size(); ++i)
        for (index_t j = i; j < thb.size(); ++j)
        {
            const real_t e  = std::abs(MA(i,j)-Mref(i,j));
            const bool ti = thb.isTruncated(i), tj = thb.isTruncated(j);
            if      (!ti&&!tj) eUU=std::max(eUU,e);
            else if ( ti&& tj) eTT=std::max(eTT,e);
            else               eUT=std::max(eUT,e);
        }
        gsInfo << "  U-U=" << eUU << "  U-T=" << eUT << "  T-T=" << eTT << "\n";
    }

    gsInfo << "\n--- Stiffness matrix ---\n";
    gsInfo << "||K_hierarchical - K_ref||_inf = "
           << (KA - Kref).cwiseAbs().maxCoeff() << "\n";

    gsInfo << "\n--- Quadrature point counts ---\n";
    gsInfo << "Gauss-Legendre (all elements):       " << gaussTotalPts << " pts\n";
    gsInfo << "Hierarchical patch, mass eval pts:  " << hprMassPts  << " pts"
           << "  (ratio " << real_t(hprMassPts)/gaussTotalPts << "x vs Gauss)\n";
    gsInfo << "Hierarchical patch, stiff eval pts: " << hprStiffPts << " pts"
           << "  (ratio " << real_t(hprStiffPts)/gaussTotalPts << "x vs Gauss)\n";
    gsInfo << "Note: points are now selected with a LOCAL evalLevel per element,\n";
    gsInfo << "      not with the global finest level everywhere.\n";

    gsInfo << "\n--- Embedded weight residuals ||N^target(X_eval) w^{target->eval} - I^target||_2 ---\n";
    gsInfo << "(Mass target only; verifies linear solve quality per direction.)\n";
    {
        for (short_t d = 0; d < 2; ++d)
        {
            for (index_t eval = 0; eval <= L; ++eval)
            {
                const gsVector<real_t>& pts = hprM.nodes(d, eval);
                for (index_t target = 0; target <= eval; ++target)
                {
                    const gsKnotVector<real_t>* kvPhys =
                        static_cast<const gsKnotVector<real_t>*>(&thb.tensorLevel(target).knots(d));
                    gsKnotVector<real_t> kvT = makeTargetKv(*kvPhys, intDegM, intRegM);
                    const gsVector<real_t> wl = computeLevelWeights(kvT, pts);

                    // Reconstruct N^target(X_eval) and I^target.
                    gsBSplineBasis<real_t> basis(kvT);
                    const index_t n = basis.size(), K = pts.size();
                    gsMatrix<real_t> N = gsMatrix<real_t>::Zero(n,K);
                    gsMatrix<real_t>      pr = pts.transpose();
                    gsMatrix<index_t>     a; gsMatrix<real_t> v;
                    basis.active_into(pr,a); basis.eval_into(pr,v);
                    for (index_t k=0;k<K;++k)
                        for (index_t r=0;r<a.rows();++r)
                            N(a(r,k),k)=v(r,k);
                    gsVector<real_t> I(n);
                    const index_t p = kvT.degree();
                    for (index_t j=0;j<n;++j)
                        I(j)=(kvT[j+p+1]-kvT[j])/real_t(p+1);

                    const real_t res = (N*wl - I).norm();
                    gsInfo << "  d=" << d
                           << "  target=" << target
                           << "  eval=" << eval
                           << "  n_target=" << n
                           << "  K_eval=" << K
                           << "  residual=" << res << "\n";
                }
            }
        }
    }

    // Plot mass-rule evaluation points using the same local evalLevel logic as
    // in the assembly loop.  The exported point cloud should no longer be the
    // full finest tensor grid over the whole patch.
    gsMatrix<> allNodes(2, hprMassPts);
    gsVector<> allWeights(hprMassPts);
    {
        gsMatrix<> tmpNodes;
        gsVector<> tmpWeights;
        gsMatrix<> ctr(2,1);
        gsMatrix<index_t> actCtr;
        index_t cnt = 0;
        for (auto& elem : thb.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,2>*>(&elem);
            const auto& low = helem->lowerCorner();
            const auto& upp = helem->upperCorner();

            ctr.col(0) = 0.5 * (low + upp);
            thb.active_into(ctr, actCtr);
            if (actCtr.rows() == 0) continue;

            index_t evalLevel = 0;
            for (index_t r = 0; r < actCtr.rows(); ++r)
            for (index_t c = r; c < actCtr.rows(); ++c)
                evalLevel = std::max(evalLevel,
                    std::max(thb.repLevel(actCtr(r,0)), thb.repLevel(actCtr(c,0))));

            hprM.mapTo(low, upp, evalLevel, evalLevel, tmpNodes, tmpWeights);
            if (tmpWeights.size() == 0) continue;
            cnt += tmpWeights.size();
            allNodes.block(0, cnt - tmpWeights.size(), 2, tmpWeights.size()) = tmpNodes;
            allWeights.segment(cnt - tmpWeights.size(), tmpWeights.size()) = tmpWeights;
        }
        GISMO_ENSURE(cnt == hprMassPts, "Mismatch in total point count.");
    }
    gsInfo<< "Exporting " << hprMassPts << " local-eval hierarchical mass points to Paraview...\n";
    gsMatrix<> X = allNodes.row(0);
    gsMatrix<> Y = allNodes.row(1);
    gsMatrix<> Z(1, hprMassPts); Z.setZero();
    gsMatrix<> W = allWeights.transpose();

    gsInfo<<"X = "<<X.rows()<<" x "<<X.cols()<<"\n";
    gsInfo<<"Y = "<<Y.rows()<<" x "<<Y.cols()<<"\n";
    gsInfo<<"Z = "<<Z.rows()<<" x "<<Z.cols()<<"\n";
    gsInfo<<"W = "<<W.rows()<<" x "<<W.cols()<<"\n";

    gsWriteParaviewPoints(X, Y, Z, W, "hrule_points");

    // Plot gauss 
    gsMatrix<> allNodesG(2, gaussTotalPts);
    gsVector<> allWeightsG(gaussTotalPts);
     {
        gsMatrix<> tmpNodes;
        gsVector<> tmpWeights;
        gsMatrix<> ctr(2,1);
        gsMatrix<index_t> actCtr;
        index_t cnt = 0;
        for (auto& elem : thb.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,2>*>(&elem);
            const auto& low = helem->lowerCorner();
            const auto& upp = helem->upperCorner();

            gaussRule->mapTo(low, upp, tmpNodes, tmpWeights);
            if (tmpWeights.size() == 0) continue;
            cnt += tmpWeights.size();
            allNodesG.block(0, cnt - tmpWeights.size(), 2, tmpWeights.size()) = tmpNodes;
            allWeightsG.segment(cnt - tmpWeights.size(), tmpWeights.size()) = tmpWeights;
        }
        GISMO_ENSURE(cnt == gaussTotalPts, "Mismatch in total point count.");
    }
    gsInfo<< "Exporting " << gaussTotalPts << " Gauss-Legendre points to Paraview...\n";
    gsMatrix<> XG = allNodesG.row(0);
    gsMatrix<> YG = allNodesG.row(1); 
    gsMatrix<> ZG(1, gaussTotalPts); ZG.setZero();
    gsMatrix<> WG = allWeightsG.transpose();
    gsWriteParaviewPoints(XG, YG, ZG, WG, "gauss_points");

    // Plot mesh
    gsMesh<> mesh(thb);
    gsWriteParaview(mesh, "hrule_mesh");
        
    return 0;
}
