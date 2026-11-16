/** @file hierarchical_patch_rule_topdown.cpp

    @brief Top-down nested hierarchical patch rule for THB-splines.

    Core idea:

      1. Construct the finest-level 1D Johannessen patch-rule point set X_L.

      2. Construct coarser point sets by nested subset selection

             X_0 subset X_1 subset ... subset X_L,

         going from L down to 0.  Coordinates are not moved.  At level l we
         select a subset of X_{l+1}; afterwards only the weights are solved.

         This is intentionally not an optimal generalized Gaussian rule at the
         coarse levels.  It is a fixed-node/moment-fitting rule.  The benefit is
         true nestedness, which lets one evaluate once at X_m and reuse the same
         values with different weight vectors w^{l->m} for all l <= m.

      3. For every target level l and every evaluation level m >= l, solve

             N^l(X_m) w^{l->m} = I^l

         in the least-squares/SVD sense.  During assembly, an element chooses an
         evaluation level m = max repLevel of the local entry group/element, and
         contributions with target level l use w^{l->m} on the same nodes X_m.

    Notes:

      * This file is a research prototype.  The nested subset selection is a
        greedy backward pruning algorithm.  It tries to remove points while the
        1D moment residual remains below a requested tolerance.  If exactness
        cannot be preserved, it keeps more points.

      * The stiffness rule below uses a safe tensor-product target
            S^{2p}_{r-1}
        rather than S^{2p-2}_{r-1} in both directions.  The latter
        under-integrates the undifferentiated transverse direction in 2D.

    Author: H.M. Verhelst, adapted prototype
*/

#include <gismo.h>
#include <gsAssembler/gsPatchRule.h>

#include <algorithm>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <vector>

using namespace gismo;

// =============================================================================
// Free-function helpers
// =============================================================================

template<class T>
gsKnotVector<T> makeTargetKv(const gsKnotVector<T>& kv,
                              index_t intDeg, index_t reg,
                              bool overInt = false)
{
    gsKnotVector<T> result = kv;

    const index_t pdiff = intDeg - result.degree();
    std::vector<index_t> mults = result.multiplicities();
    const index_t rmin  = *std::min_element(mults.begin(), mults.end());
    const index_t rdiff = (intDeg - reg) - rmin;

    if      (pdiff > 0) result.degreeIncrease(pdiff);
    else if (pdiff < 0) result.degreeDecrease(-pdiff);

    if (rdiff > 0)
        result.increaseMultiplicity(rdiff);

    index_t size = result.size() - result.degree() - 1;

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

template<class T>
gsMatrix<T> buildMomentMatrix(const gsKnotVector<T>& kvTarget,
                              const gsVector<T>& pts)
{
    gsBSplineBasis<T> basis(kvTarget);
    const index_t n = basis.size();
    const index_t K = pts.size();

    gsMatrix<T> N = gsMatrix<T>::Zero(n, K);
    if (K == 0)
        return N;

    gsMatrix<T>      ptsRow = pts.transpose();
    gsMatrix<index_t> acts;
    gsMatrix<T>       vals;
    basis.active_into(ptsRow, acts);
    basis.eval_into  (ptsRow, vals);

    for (index_t k = 0; k < K; ++k)
        for (index_t r = 0; r < acts.rows(); ++r)
            N(acts(r,k), k) = vals(r,k);

    return N;
}

template<class T>
gsVector<T> buildMomentRhs(const gsKnotVector<T>& kvTarget)
{
    gsBSplineBasis<T> basis(kvTarget);
    const index_t n = basis.size();
    const index_t p = kvTarget.degree();

    gsVector<T> I(n);
    for (index_t j = 0; j < n; ++j)
        I(j) = (kvTarget[j + p + 1] - kvTarget[j]) / T(p + 1);
    return I;
}

template<class T>
gsVector<T> computeLevelWeights(const gsKnotVector<T>& kvTarget,
                                const gsVector<T>& pts)
{
    const gsMatrix<T> N = buildMomentMatrix(kvTarget, pts);
    const gsVector<T> I = buildMomentRhs(kvTarget);

    if (pts.size() == 0)
        return gsVector<T>();

    return N.bdcSvd(gsEigen::ComputeThinU | gsEigen::ComputeThinV).solve(I);
}

template<class T>
T momentResidual(const gsKnotVector<T>& kvTarget,
                 const gsVector<T>& pts,
                 gsVector<T>* weightsOut = nullptr)
{
    const gsMatrix<T> N = buildMomentMatrix(kvTarget, pts);
    const gsVector<T> I = buildMomentRhs(kvTarget);

    if (pts.size() == 0)
        return I.norm();

    const gsVector<T> w =
        N.bdcSvd(gsEigen::ComputeThinU | gsEigen::ComputeThinV).solve(I);

    if (weightsOut)
        *weightsOut = w;

    return (N*w - I).norm();
}

template<class T>
gsVector<T> subsetByIndices(const gsVector<T>& candidates,
                            const std::vector<index_t>& idx)
{
    gsVector<T> out(idx.size());
    for (index_t k = 0; k < static_cast<index_t>(idx.size()); ++k)
        out[k] = candidates[idx[k]];
    return out;
}

/**
 * @brief Greedy top-down nested subset selection.
 *
 * Given candidates X_{l+1}, select X_l subset X_{l+1}.  Coordinates are fixed;
 * the weights are recomputed after every tentative removal.  A point is removed
 * only if the relative residual remains below @p relTol.  The algorithm stops
 * once @p minCount is reached or no admissible removal exists.
 */
template<class T>
gsVector<T> selectNestedSubsetGreedy(const gsKnotVector<T>& kvTarget,
                                     const gsVector<T>& candidates,
                                     index_t minCount,
                                     T relTol,
                                     T absTol,
                                     bool verbose,
                                     short_t dir,
                                     index_t level)
{
    const gsVector<T> I = buildMomentRhs(kvTarget);
    const T rhsNorm = std::max<T>(I.norm(), T(1));

    std::vector<index_t> idx(candidates.size());
    std::iota(idx.begin(), idx.end(), index_t(0));

    minCount = std::min<index_t>(minCount, candidates.size());

    T currentRes = momentResidual(kvTarget, candidates);
    if (verbose)
        gsInfo << "    d=" << dir << " l=" << level
               << " start K=" << candidates.size()
               << " minK=" << minCount
               << " relRes=" << currentRes/rhsNorm << "\n";

    // If even the full candidate set cannot hit the moments, keep all points.
    if (currentRes > absTol && currentRes/rhsNorm > relTol)
    {
        if (verbose)
            gsInfo << "      full candidate residual too large; keeping all points\n";
        return candidates;
    }

    while (static_cast<index_t>(idx.size()) > minCount)
    {
        T bestRes = std::numeric_limits<T>::max();
        index_t bestPos = -1;

        for (index_t pos = 0; pos < static_cast<index_t>(idx.size()); ++pos)
        {
            std::vector<index_t> trial = idx;
            trial.erase(trial.begin() + pos);

            const gsVector<T> trialPts = subsetByIndices(candidates, trial);
            const T res = momentResidual(kvTarget, trialPts);

            if (res < bestRes)
            {
                bestRes = res;
                bestPos = pos;
            }
        }

        if (bestPos < 0)
            break;

        if (bestRes <= absTol || bestRes/rhsNorm <= relTol)
        {
            idx.erase(idx.begin() + bestPos);
            currentRes = bestRes;
        }
        else
        {
            break;
        }
    }

    const gsVector<T> selected = subsetByIndices(candidates, idx);

    if (verbose)
        gsInfo << "      selected K=" << selected.size()
               << " relRes=" << momentResidual(kvTarget, selected)/rhsNorm
               << "\n";

    return selected;
}

// =============================================================================
// gsTopDownNestedPatchRule
// =============================================================================

template<class T>
class gsTopDownNestedPatchRule
{
public:
    /**
     * @param basis        THB/HB tensor hierarchy
     * @param intDeg       Target integration degree
     * @param reg          Target integration regularity
     * @param overInt      Boundary over-integration flag
     * @param relTol       Greedy subset residual tolerance
     * @param safety       Extra points kept beyond the nominal fixed-node count
     * @param verbose      Print subset construction diagnostics
     */
    gsTopDownNestedPatchRule(const gsHTensorBasis<2,T>& basis,
                             index_t intDeg, index_t reg,
                             bool overInt = false,
                             T relTol = T(1e-12),
                             index_t safety = 0,
                             bool verbose = true)
        : m_L(basis.maxLevel()),
          m_intDeg(intDeg),
          m_reg(reg),
          m_overInt(overInt)
    {
        m_end = basis.support().col(1);
        m_nodes.resize(2);
        m_nodeMap.resize(2);
        m_weights.resize(2);

        gsOptionList opts;
        opts.addInt   ("quRule", "", gsQuadrature::PatchRule);
        opts.addReal  ("quA",    "", T(intDeg));
        opts.addInt   ("quB",    "", reg);
        opts.addSwitch("overInt", "", overInt);

        auto finestRule = gsQuadrature::getPtr(basis.tensorLevel(m_L), opts);
        const gsPatchRule<T>* fpr =
            dynamic_cast<const gsPatchRule<T>*>(finestRule.get());
        GISMO_ENSURE(fpr,
            "gsTopDownNestedPatchRule: gsQuadrature did not produce a gsPatchRule.");

        for (short_t d = 0; d < 2; ++d)
        {
            m_nodes[d].resize(m_L + 1);
            m_nodeMap[d].resize(m_L + 1);
            m_weights[d].resize(m_L + 1);

            // Finest level: ordinary Johannessen nodes.
            m_nodes[d][m_L] = fpr->nodes(d);

            // Top-down nested subsets: X_l subset X_{l+1}.
            for (index_t lev = m_L; lev > 0; --lev)
            {
                const index_t l = lev - 1;
                const gsKnotVector<T>* kvPhys =
                    static_cast<const gsKnotVector<T>*>(&basis.tensorLevel(l).knots(d));
                const gsKnotVector<T> kvT = makeTargetKv(*kvPhys, intDeg, reg, overInt);
                gsBSplineBasis<T> targetBasis(kvT);

                // Fixed-node rules normally need around dim(T_l) points.
                // If X_{l+1} has fewer points than dim(T_l), we keep all of it.
                const index_t nominal = targetBasis.size() + safety;
                const index_t minCount = std::min<index_t>(m_nodes[d][lev].size(), nominal);

                m_nodes[d][l] = selectNestedSubsetGreedy(kvT,
                                                         m_nodes[d][lev],
                                                         minCount,
                                                         relTol,
                                                         T(1e-14),
                                                         verbose,
                                                         d,
                                                         l);
            }

            // Build maps for each eval level.
            for (index_t m = 0; m <= m_L; ++m)
            {
                for (index_t k = 0; k < m_nodes[d][m].size(); ++k)
                    m_nodeMap[d][m][m_nodes[d][m][k]] = k;
            }

            // Build embedded weights w^{target -> eval} for all target <= eval.
            for (index_t target = 0; target <= m_L; ++target)
            {
                m_weights[d][target].resize(m_L + 1);

                const gsKnotVector<T>* kvPhys =
                    static_cast<const gsKnotVector<T>*>(&basis.tensorLevel(target).knots(d));
                const gsKnotVector<T> kvT = makeTargetKv(*kvPhys, intDeg, reg, overInt);

                for (index_t eval = target; eval <= m_L; ++eval)
                    m_weights[d][target][eval] = computeLevelWeights(kvT, m_nodes[d][eval]);
            }
        }
    }

    /**
     * @brief Return nodes X_eval and weights w^{targetLevel -> evalLevel} on
     *        the element [lower,upper].
     */
    void mapTo(const gsVector<T>& lower, const gsVector<T>& upper,
               index_t targetLevel, index_t evalLevel,
               gsMatrix<T>& nodes, gsVector<T>& weights) const
    {
        GISMO_ASSERT(targetLevel >= 0 && targetLevel <= m_L,
            "targetLevel " << targetLevel << " out of range [0," << m_L << "]");
        GISMO_ASSERT(evalLevel >= 0 && evalLevel <= m_L,
            "evalLevel " << evalLevel << " out of range [0," << m_L << "]");
        GISMO_ASSERT(targetLevel <= evalLevel,
            "Need targetLevel <= evalLevel for embedded weights.");

        std::vector<std::vector<index_t> > idxs(2);
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

    void getWeights(const gsVector<T>& lower, const gsVector<T>& upper,
                    index_t targetLevel, index_t evalLevel,
                    gsVector<T>& weights) const
    {
        gsMatrix<T> nodes;
        mapTo(lower, upper, targetLevel, evalLevel, nodes, weights);
    }

    index_t maxLevel() const { return m_L; }

    const gsVector<T>& nodes(short_t d, index_t evalLevel) const
    {
        return m_nodes[d][evalLevel];
    }

    index_t numNodes(short_t d, index_t evalLevel) const
    {
        return m_nodes[d][evalLevel].size();
    }

    bool checkNested(short_t d, T tol = T(1e-12)) const
    {
        for (index_t l = 0; l < m_L; ++l)
        {
            const gsVector<T>& coarse = m_nodes[d][l];
            const gsVector<T>& fine   = m_nodes[d][l+1];
            for (index_t i = 0; i < coarse.size(); ++i)
            {
                bool found = false;
                for (index_t j = 0; j < fine.size(); ++j)
                    if (math::abs(coarse[i] - fine[j]) <= tol)
                    {
                        found = true;
                        break;
                    }
                if (!found)
                    return false;
            }
        }
        return true;
    }

private:
    void _tensorProduct(const std::vector<std::vector<index_t> >& idxs,
                        index_t targetLevel, index_t evalLevel,
                        gsMatrix<T>& nodes, gsVector<T>& weights) const
    {
        const index_t nx = static_cast<index_t>(idxs[0].size());
        const index_t ny = static_cast<index_t>(idxs[1].size());
        const index_t sz = nx * ny;

        nodes.resize(2, sz);
        weights.resize(sz);

        if (sz == 0) return;

        const gsVector<T>& wx = m_weights[0][targetLevel][evalLevel];
        const gsVector<T>& wy = m_weights[1][targetLevel][evalLevel];

        index_t cnt = 0;
        for (index_t jy = 0; jy < ny; ++jy)
        {
            const index_t ky = idxs[1][jy];
            for (index_t jx = 0; jx < nx; ++jx, ++cnt)
            {
                const index_t kx = idxs[0][jx];
                nodes(0, cnt) = m_nodes[0][evalLevel][kx];
                nodes(1, cnt) = m_nodes[1][evalLevel][ky];
                weights[cnt]  = wx[kx] * wy[ky];
            }
        }
    }

private:
    index_t     m_L;
    index_t     m_intDeg;
    index_t     m_reg;
    bool        m_overInt;
    gsVector<T> m_end;

    // [direction][evalLevel]
    std::vector<std::vector<gsVector<T> > > m_nodes;

    // [direction][evalLevel] : coordinate -> local index in m_nodes[direction][evalLevel]
    std::vector<std::vector<std::map<T,index_t> > > m_nodeMap;

    // [direction][targetLevel][evalLevel]
    std::vector<std::vector<std::vector<gsVector<T> > > > m_weights;
};

// =============================================================================
// main
// =============================================================================

int main(int argc, char* argv[])
{
    index_t order      = 2;
    index_t regularity = 1;
    bool    plot       = false;
    bool    verbose    = false;
    index_t safety     = 0;
    real_t  relTol     = 1e-12;

    gsCmdLine cmd("Top-down nested hierarchical patch rule: fixed coordinates, embedded weights.");
    cmd.addInt   ("P", "deg",     "Polynomial order of basis", order);
    cmd.addInt   ("R", "reg",     "Regularity of basis", regularity);
    cmd.addInt   ("s", "safety",  "Extra points kept in greedy nested subset", safety);
    cmd.addReal  ("t", "tol",     "Relative residual tolerance for subset pruning", relTol);
    cmd.addSwitch("plot",          "Export mesh and quadrature points to Paraview", plot);
    cmd.addSwitch("verbose",       "Verbose nested subset construction", verbose);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(regularity < order, "Require regularity < order.");

    // Build THB-spline basis.
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

    // Integration targets.
    const index_t intDegM = 2 * order;
    const index_t intRegM = regularity;

    // Safe stiffness target for 2D gradient products:
    // derivative direction needs regularity-1, transverse direction needs degree 2p.
    // This safe tensor target uses degree 2p and regularity r-1 in every direction.
    const index_t intDegK = 2 * order;
    const index_t intRegK = regularity > 0 ? regularity - 1 : 0;

    // Reference: Gauss-Legendre on the leaf hierarchy.
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

    // Top-down nested rules.
    gsTopDownNestedPatchRule<real_t> hprM(thb, intDegM, intRegM,
                                          false, relTol, safety, verbose);
    gsTopDownNestedPatchRule<real_t> hprK(thb, intDegK, intRegK,
                                          false, relTol, safety, verbose);

    gsInfo << "\n--- 1D nested node counts, mass rule ---\n";
    for (short_t d = 0; d < 2; ++d)
    {
        gsInfo << "  d=" << d << ":";
        for (index_t l = 0; l <= L; ++l)
            gsInfo << "  K" << l << "=" << hprM.numNodes(d,l);
        gsInfo << "  nested=" << (hprM.checkNested(d) ? "yes" : "NO") << "\n";
    }
    gsInfo << "--- 1D nested node counts, stiffness rule ---\n";
    for (short_t d = 0; d < 2; ++d)
    {
        gsInfo << "  d=" << d << ":";
        for (index_t l = 0; l <= L; ++l)
            gsInfo << "  K" << l << "=" << hprK.numNodes(d,l);
        gsInfo << "  nested=" << (hprK.checkNested(d) ? "yes" : "NO") << "\n";
    }

    // Assembly with local evaluation-level dispatch.
    gsVector<real_t> IA(thb.size()); IA.setZero();
    gsMatrix<real_t> MA = gsMatrix<real_t>::Zero(thb.size(), thb.size());
    gsMatrix<real_t> KA = gsMatrix<real_t>::Zero(thb.size(), thb.size());
    index_t hprMassPts = 0;
    index_t hprStiffPts = 0;
    {
        gsMatrix<real_t> nodesM, nodesK;
        gsVector<real_t> wEval, wLp, wEvalK, wLpK;
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

            ctr.col(0) = 0.5 * (low + upp);
            thb.active_into(ctr, actCtr);

            std::set<index_t> Lp_set;
            for (index_t r = 0; r < actCtr.rows(); ++r)
                for (index_t c = r; c < actCtr.rows(); ++c)
                    Lp_set.insert(std::max(thb.repLevel(actCtr(r,0)),
                                           thb.repLevel(actCtr(c,0))));

            if (Lp_set.empty()) continue;
            const index_t evalLevel = *Lp_set.rbegin();

            // Shared mass nodes on local evaluation level.
            hprM.mapTo(low, upp, evalLevel, evalLevel, nodesM, wEval);
            if (wEval.size() == 0) continue;
            hprMassPts += wEval.size();

            thb.active_into(nodesM, acts);
            thb.eval_into  (nodesM, vals);

            // Shared stiffness nodes on the same local evaluation level, but for
            // the stiffness target rule.
            hprK.mapTo(low, upp, evalLevel, evalLevel, nodesK, wEvalK);
            if (wEvalK.size() > 0)
            {
                hprStiffPts += wEvalK.size();
                thb.active_into(nodesK, actsK);
                thb.deriv_into (nodesK, derivs);
            }

            for (index_t Lp : Lp_set)
            {
                // Mass weights w^{Lp -> evalLevel}.
                hprM.getWeights(low, upp, Lp, evalLevel, wLp);
                if (wLp.size() == 0) continue;

                for (index_t r = 0; r < acts.rows(); ++r)
                {
                    const index_t i = acts(r,0);
                    if (thb.repLevel(i) != Lp) continue;
                    for (index_t k = 0; k < wLp.size(); ++k)
                        IA[i] += wLp[k] * vals(r,k);
                }

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

                // Stiffness weights w^{Lp -> evalLevel}.
                if (wEvalK.size() > 0)
                {
                    hprK.getWeights(low, upp, Lp, evalLevel, wLpK);
                    for (index_t r = 0; r < actsK.rows(); ++r)
                        for (index_t c = r; c < actsK.rows(); ++c)
                        {
                            const index_t i = actsK(r,0), j = actsK(c,0);
                            if (std::max(thb.repLevel(i), thb.repLevel(j)) != Lp) continue;
                            real_t sk = 0;
                            for (index_t k = 0; k < wLpK.size(); ++k)
                                sk += wLpK[k] *
                                      (derivs(2*r+0,k)*derivs(2*c+0,k) +
                                       derivs(2*r+1,k)*derivs(2*c+1,k));
                            KA(i,j) += sk; if (i != j) KA(j,i) += sk;
                        }
                }
            }
        }
    }

    // Results.
    gsInfo << "\n========== Top-Down Nested Hierarchical Patch Rule Results ==========" << "\n";

    gsInfo << "\n--- Scalar integrals ---\n";
    gsInfo << "||I_nested - I_ref||_inf = "
           << (IA - Iref).cwiseAbs().maxCoeff() << "\n";
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
    gsInfo << "||M_nested - M_ref||_inf = "
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
    gsInfo << "||K_nested - K_ref||_inf = "
           << (KA - Kref).cwiseAbs().maxCoeff() << "\n";

    gsInfo << "\n--- Quadrature point counts ---\n";
    gsInfo << "Gauss-Legendre (all elements):         " << gaussTotalPts << " pts\n";
    gsInfo << "Nested patch mass eval points:        " << hprMassPts << " pts"
           << "  (ratio " << real_t(hprMassPts)/gaussTotalPts << "x vs Gauss)\n";
    gsInfo << "Nested patch stiffness eval points:   " << hprStiffPts << " pts"
           << "  (ratio " << real_t(hprStiffPts)/gaussTotalPts << "x vs Gauss)\n";
    gsInfo << "Note: these are evaluation points with local evalLevel dispatch, not finest-everywhere points.\n";

    gsInfo << "\n--- Embedded 1D weight residuals for mass rule ---\n";
    for (short_t d = 0; d < 2; ++d)
    {
        for (index_t target = 0; target <= L; ++target)
        {
            const gsKnotVector<real_t>* kvPhys =
                static_cast<const gsKnotVector<real_t>*>(&thb.tensorLevel(target).knots(d));
            const gsKnotVector<real_t> kvT = makeTargetKv(*kvPhys, intDegM, intRegM);
            const gsVector<real_t> rhs = buildMomentRhs(kvT);
            const real_t rhsNorm = std::max<real_t>(rhs.norm(), 1.0);

            for (index_t eval = target; eval <= L; ++eval)
            {
                gsVector<real_t> w;
                const real_t res = momentResidual(kvT, hprM.nodes(d, eval), &w);
                gsInfo << "  d=" << d
                       << "  target=" << target
                       << "  eval=" << eval
                       << "  K_eval=" << hprM.numNodes(d, eval)
                       << "  relRes=" << res/rhsNorm << "\n";
            }
        }
    }

    // Plot mass evaluation points with local evalLevel.
    if (plot)
    {
        gsMatrix<> allNodes(2, hprMassPts);
        gsVector<> allWeights(hprMassPts);
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

            index_t evalLevel = 0;
            for (index_t r = 0; r < actCtr.rows(); ++r)
                evalLevel = std::max(evalLevel, thb.repLevel(actCtr(r,0)));

            hprM.mapTo(low, upp, evalLevel, evalLevel, tmpNodes, tmpWeights);
            if (tmpWeights.size() == 0) continue;

            allNodes.block(0, cnt, 2, tmpWeights.size()) = tmpNodes;
            allWeights.segment(cnt, tmpWeights.size()) = tmpWeights;
            cnt += tmpWeights.size();
        }
        GISMO_ENSURE(cnt == hprMassPts, "Mismatch in total point count.");

        gsInfo<< "Exporting " << hprMassPts
              << " nested hierarchical patch rule points to Paraview...\n";
        gsMatrix<> X = allNodes.row(0);
        gsMatrix<> Y = allNodes.row(1);
        gsMatrix<> Z(1, hprMassPts); Z.setZero();
        gsMatrix<> W = allWeights.transpose();
        gsWriteParaviewPoints(X, Y, Z, W, "hrule_points_nested");
    }

    return 0;
}
