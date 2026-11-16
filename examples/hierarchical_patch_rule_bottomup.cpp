/** @file hierarchical_patch_rule_bottomup.cpp

    @brief Bottom-up nested hierarchical patch rule experiment for THB-splines.

    Core idea:

      Construct one nested family of 1-D quadrature point sets per parametric
      direction,

          X_0 subset X_1 subset ... subset X_L,

      by going from coarse to fine.  X_0 is the ordinary Johannessen patch rule.
      To construct X_{l+1}, keep all coordinates in X_l fixed, insert additional
      coordinates in the largest gaps, and solve a constrained nonlinear moment
      system for

          - all weights on X_{l+1}, and
          - only the newly inserted coordinates.

      The old coordinates are constrained; the old weights are not.

      After all coordinates have been constructed, embedded weights

          w^{targetLevel -> evalLevel}

      are computed by linear moment fitting on X_evalLevel.  During assembly one
      evaluates once on X_evalLevel, where evalLevel is the largest level needed
      by the local entries, and then reuses the same values with different
      weight vectors.

    Notes:

      - This file is an experimental prototype.  I kept the structure close to
        the original standalone example.
      - The nonlinear solve is deliberately damped and conservative.  If it does
        not converge, the file falls back to the current coordinates and reports
        the residual instead of hanging in a greedy subset-selection loop.
*/

#include <gismo.h>
#include <gsAssembler/gsPatchRule.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
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
        const T lLen = (result[1] - result.first()) / T(numOver + 1);
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
gsVector<T> sortUnique(const gsVector<T>& x, T tol = T(1e-13))
{
    std::vector<T> v(x.size());
    for (index_t i = 0; i < x.size(); ++i) v[i] = x[i];
    std::sort(v.begin(), v.end());

    std::vector<T> u;
    for (typename std::vector<T>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
        if (u.empty() || std::abs(*it - u.back()) > tol)
            u.push_back(*it);
    }

    gsVector<T> y(u.size());
    for (index_t i = 0; i < static_cast<index_t>(u.size()); ++i)
        y[i] = u[i];
    return y;
}

template<class T>
gsVector<T> concatenate(const gsVector<T>& a, const gsVector<T>& b)
{
    gsVector<T> c(a.size() + b.size());
    for (index_t i = 0; i < a.size(); ++i) c[i] = a[i];
    for (index_t i = 0; i < b.size(); ++i) c[a.size()+i] = b[i];
    return c;
}

template<class T>
bool validPointSet(const gsVector<T>& x, T a, T b, T tol = T(1e-12))
{
    if (x.size() == 0) return true;
    for (index_t i = 0; i < x.size(); ++i)
        if (!(x[i] > a + tol && x[i] < b - tol)) return false;

    gsVector<T> xs = sortUnique(x, tol);
    return xs.size() == x.size();
}

/**
 * @brief Insert @p numNew points in the largest gaps of @p oldPts.
 *
 * The endpoints a,b are treated as fixed artificial gap boundaries, but no new
 * point is inserted exactly at an endpoint.
 */
template<class T>
gsVector<T> insertInLargestGaps(const gsVector<T>& oldPts,
                                index_t numNew,
                                T a,
                                T b)
{
    gsVector<T> old = sortUnique(oldPts);

    struct Gap { T left; T right; };
    std::vector<Gap> gaps;

    if (old.size() == 0)
        gaps.push_back(Gap{a,b});
    else
    {
        gaps.push_back(Gap{a, old[0]});
        for (index_t i = 0; i + 1 < old.size(); ++i)
            gaps.push_back(Gap{old[i], old[i+1]});
        gaps.push_back(Gap{old[old.size()-1], b});
    }

    std::vector<T> newPts;
    for (index_t k = 0; k < numNew; ++k)
    {
        index_t best = 0;
        T bestLen = T(-1);
        for (index_t g = 0; g < static_cast<index_t>(gaps.size()); ++g)
        {
            const T len = gaps[g].right - gaps[g].left;
            if (len > bestLen)
            {
                bestLen = len;
                best = g;
            }
        }

        const T mid = T(0.5) * (gaps[best].left + gaps[best].right);
        newPts.push_back(mid);

        const Gap oldGap = gaps[best];
        gaps.erase(gaps.begin() + best);
        gaps.push_back(Gap{oldGap.left, mid});
        gaps.push_back(Gap{mid, oldGap.right});
    }

    gsVector<T> y(newPts.size());
    for (index_t i = 0; i < static_cast<index_t>(newPts.size()); ++i)
        y[i] = newPts[i];
    return sortUnique(y);
}

/**
 * @brief Build moment matrix A, exact moment vector I, and optionally dA/dx.
 */
template<class T>
void buildMomentSystem(const gsKnotVector<T>& kvTarget,
                       const gsVector<T>& pts,
                       gsMatrix<T>& A,
                       gsVector<T>& I,
                       gsMatrix<T>* dA = nullptr)
{
    gsBSplineBasis<T> basis(kvTarget);
    const index_t n = basis.size();
    const index_t K = pts.size();
    const index_t p = kvTarget.degree();

    A = gsMatrix<T>::Zero(n, K);
    if (dA) *dA = gsMatrix<T>::Zero(n, K);

    gsMatrix<T> ptsRow = pts.transpose();
    gsMatrix<index_t> acts;
    gsMatrix<T> vals, dvals;

    basis.active_into(ptsRow, acts);
    basis.eval_into  (ptsRow, vals);
    if (dA) basis.deriv_into(ptsRow, dvals);

    for (index_t k = 0; k < K; ++k)
    {
        for (index_t r = 0; r < acts.rows(); ++r)
        {
            A(acts(r,k), k) = vals(r,k);
            if (dA) (*dA)(acts(r,k), k) = dvals(r,k);
        }
    }

    I.resize(n);
    for (index_t j = 0; j < n; ++j)
        I(j) = (kvTarget[j + p + 1] - kvTarget[j]) / T(p + 1);
}

template<class T>
gsVector<T> computeLevelWeights(const gsKnotVector<T>& kvTarget,
                                 const gsVector<T>& pts)
{
    gsMatrix<T> A;
    gsVector<T> I;
    buildMomentSystem(kvTarget, pts, A, I, static_cast<gsMatrix<T>*>(nullptr));
    return A.bdcSvd(gsEigen::ComputeThinU | gsEigen::ComputeThinV).solve(I);
}

template<class T>
T momentResidual(const gsKnotVector<T>& kvTarget,
                 const gsVector<T>& pts,
                 const gsVector<T>& weights)
{
    gsMatrix<T> A;
    gsVector<T> I;
    buildMomentSystem(kvTarget, pts, A, I, static_cast<gsMatrix<T>*>(nullptr));
    return (A * weights - I).norm();
}

/**
 * @brief Construct X_{l+1} from X_l by constrained nonlinear moment fitting.
 *
 * Old coordinates are fixed.  New coordinates and all weights are free.  The
 * returned point vector is sorted; weights are intentionally not returned here,
 * because all embedded weights are recomputed afterwards by linear moment fits.
 */
template<class T>
gsVector<T> constructEmbeddedNextLevel(const gsVector<T>& oldPtsIn,
                                       const gsKnotVector<T>& kvTargetNext,
                                       index_t safety,
                                       T tol,
                                       index_t maxIter)
{
    const T a = kvTargetNext.first();
    const T b = kvTargetNext.last();
    gsVector<T> oldPts = sortUnique(oldPtsIn);

    gsBSplineBasis<T> basis(kvTargetNext);
    const index_t n = basis.size();
    const index_t Kold = oldPts.size();

    index_t numNew = 0;
    if (n > Kold)
        numNew = static_cast<index_t>(std::ceil(T(n - Kold) / T(2))) + safety;
    else
        numNew = safety;

    // Keep the nonlinear systems moderate.  A small safety usually helps;
    // excessive safety makes the SVD solve more expensive and is not needed.
    numNew = std::max<index_t>(numNew, 1);

    gsVector<T> newPts = insertInLargestGaps(oldPts, numNew, a, b);
    if (newPts.size() != numNew)
        gsInfo << "WARNING: duplicate new points removed during initialization. "
               << "Requested " << numNew << ", got " << newPts.size() << "\n";
    numNew = newPts.size();

    gsVector<T> allPts = concatenate(oldPts, newPts);
    gsVector<T> weights = computeLevelWeights(kvTargetNext, allPts);

    T res0 = momentResidual(kvTargetNext, allPts, weights);
    T res  = res0;

    gsInfo << "    embedded solve: n=" << n
           << ", Kold=" << Kold
           << ", Knew=" << numNew
           << ", Ktotal=" << allPts.size()
           << ", initial residual=" << res0 << "\n";

    for (index_t it = 0; it < maxIter; ++it)
    {
        gsMatrix<T> A, dA;
        gsVector<T> I;
        buildMomentSystem(kvTargetNext, allPts, A, I, &dA);
        gsVector<T> F = A * weights - I;
        res = F.norm();

        if (res < tol)
        {
            gsInfo << "      converged at it=" << it << ", residual=" << res << "\n";
            return sortUnique(allPts);
        }

        const index_t K = allPts.size();
        gsMatrix<T> J = gsMatrix<T>::Zero(n, K + numNew);
        J.block(0, 0, n, K) = A;

        for (index_t s = 0; s < numNew; ++s)
        {
            const index_t colInA = Kold + s;
            J.col(K + s) = dA.col(colInA) * weights[colInA];
        }

        gsVector<T> dz = J.bdcSvd(gsEigen::ComputeThinU | gsEigen::ComputeThinV).solve(-F);

        bool accepted = false;
        T alpha = T(1);
        for (index_t ls = 0; ls < 25; ++ls)
        {
            gsVector<T> wTrial = weights + alpha * dz.head(K);
            gsVector<T> yTrial = newPts  + alpha * dz.tail(numNew);
            gsVector<T> xTrial = concatenate(oldPts, yTrial);

            if (!validPointSet(xTrial, a, b))
            {
                alpha *= T(0.5);
                continue;
            }

            const T trialRes = momentResidual(kvTargetNext, xTrial, wTrial);
            if (trialRes < res)
            {
                weights = wTrial;
                newPts  = yTrial;
                allPts  = xTrial;
                accepted = true;
                break;
            }

            alpha *= T(0.5);
        }

        if (!accepted)
        {
            gsInfo << "      WARNING: line search failed at it=" << it
                   << ", residual=" << res << ". Returning current points.\n";
            break;
        }
    }

    weights = computeLevelWeights(kvTargetNext, allPts);
    res = momentResidual(kvTargetNext, allPts, weights);
    gsInfo << "      WARNING: maxIter reached/fallback, final residual=" << res << "\n";
    return sortUnique(allPts);
}

// =============================================================================
// gsBottomUpNestedPatchRule
// =============================================================================

template<class T>
class gsBottomUpNestedPatchRule
{
public:
    /**
     * @param basis     THB/HB tensor basis.
     * @param intDeg    Target integration degree.
     * @param reg       Target integration regularity.
     * @param overInt   Boundary over-integration flag.
     * @param safety    Extra new points added at each embedded refinement step.
     * @param tol       Nonlinear moment residual tolerance.
     * @param maxIter   Maximum damped Newton iterations per embedded step.
     */
    gsBottomUpNestedPatchRule(const gsHTensorBasis<2,T>& basis,
                              index_t intDeg, index_t reg,
                              bool overInt = false,
                              index_t safety = 1,
                              T tol = T(1e-12),
                              index_t maxIter = 80)
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

            std::vector<gsKnotVector<T> > targetKv(m_L + 1);
            for (index_t l = 0; l <= m_L; ++l)
            {
                const gsKnotVector<T>* kvPhys =
                    static_cast<const gsKnotVector<T>*>(&basis.tensorLevel(l).knots(d));
                targetKv[l] = makeTargetKv(*kvPhys, intDeg, reg, overInt);
            }

            // Level 0: ordinary Johannessen rule.
            auto rule0 = gsQuadrature::getPtr(basis.tensorLevel(0), opts);
            const gsPatchRule<T>* pr0 = dynamic_cast<const gsPatchRule<T>*>(rule0.get());
            GISMO_ENSURE(pr0, "Bottom-up nested patch rule: PatchRule construction failed at level 0.");
            m_nodes[d][0] = sortUnique(pr0->nodes(d));

            gsInfo << "  d=" << d << " L0 nodes=" << m_nodes[d][0].size() << "\n";

            // Levels 1..L: embedded construction, preserving old coordinates.
            for (index_t l = 0; l < m_L; ++l)
            {
                gsInfo << "  d=" << d << " constructing embedded level " << (l+1) << "\n";
                m_nodes[d][l+1] = constructEmbeddedNextLevel(
                    m_nodes[d][l], targetKv[l+1], safety, tol, maxIter);
                gsInfo << "  d=" << d << " L" << (l+1)
                       << " nodes=" << m_nodes[d][l+1].size() << "\n";
            }

            // Build point maps for element extraction.
            for (index_t l = 0; l <= m_L; ++l)
            {
                for (index_t k = 0; k < m_nodes[d][l].size(); ++k)
                    m_nodeMap[d][l][m_nodes[d][l][k]] = k;
            }

            // Embedded weights: target level ell integrated on eval level m.
            m_levelW[d].resize(m_L + 1);
            for (index_t ell = 0; ell <= m_L; ++ell)
            {
                m_levelW[d][ell].resize(m_L + 1);
                for (index_t m = ell; m <= m_L; ++m)
                {
                    m_levelW[d][ell][m] = computeLevelWeights(targetKv[ell], m_nodes[d][m]);
                }
            }
        }
    }

    void mapTo(const gsVector<T>& lower, const gsVector<T>& upper,
               index_t targetLevel, index_t evalLevel,
               gsMatrix<T>& nodes, gsVector<T>& weights) const
    {
        GISMO_ASSERT(targetLevel >= 0 && targetLevel <= m_L,
            "targetLevel " << targetLevel << " out of range [0," << m_L << "]");
        GISMO_ASSERT(evalLevel >= 0 && evalLevel <= m_L,
            "evalLevel " << evalLevel << " out of range [0," << m_L << "]");
        GISMO_ASSERT(targetLevel <= evalLevel,
            "targetLevel must be <= evalLevel for embedded weights.");

        std::vector<std::vector<index_t> > idxs(2);
        for (short_t d = 0; d < 2; ++d)
        {
            const bool isLast = (upper[d] >= m_end[d] - T(1e-14));
            const std::map<T,index_t>& mp = m_nodeMap[d][evalLevel];
            typename std::map<T,index_t>::const_iterator itLo = mp.lower_bound(lower[d]);
            typename std::map<T,index_t>::const_iterator itHi = isLast ? mp.end() : mp.lower_bound(upper[d]);
            for (typename std::map<T,index_t>::const_iterator it = itLo; it != itHi; ++it)
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

    index_t numNodes1D(short_t d, index_t l) const { return m_nodes[d][l].size(); }

    bool checkNested(short_t d, T tol = T(1e-12)) const
    {
        for (index_t l = 0; l < m_L; ++l)
        {
            for (index_t i = 0; i < m_nodes[d][l].size(); ++i)
            {
                bool found = false;
                for (index_t j = 0; j < m_nodes[d][l+1].size(); ++j)
                    if (std::abs(m_nodes[d][l][i] - m_nodes[d][l+1][j]) < tol)
                    {
                        found = true;
                        break;
                    }
                if (!found) return false;
            }
        }
        return true;
    }

private:
    void _tensorProduct(const std::vector<std::vector<index_t> >& idxs,
                        index_t targetLevel,
                        index_t evalLevel,
                        gsMatrix<T>& nodes,
                        gsVector<T>& weights) const
    {
        const index_t nx = static_cast<index_t>(idxs[0].size());
        const index_t ny = static_cast<index_t>(idxs[1].size());
        const index_t sz = nx * ny;

        nodes.resize(2, sz);
        weights.resize(sz);
        if (sz == 0) return;

        const gsVector<T>& wx = m_levelW[0][targetLevel][evalLevel];
        const gsVector<T>& wy = m_levelW[1][targetLevel][evalLevel];

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
    gsVector<T> m_end;

    // m_nodes[d][evalLevel]
    std::vector<std::vector<gsVector<T> > > m_nodes;

    // m_nodeMap[d][evalLevel]: coordinate -> index in m_nodes[d][evalLevel]
    std::vector<std::vector<std::map<T,index_t> > > m_nodeMap;

    // m_levelW[d][targetLevel][evalLevel]
    std::vector<std::vector<std::vector<gsVector<T> > > > m_levelW;
};

// =============================================================================
// main
// =============================================================================

int main(int argc, char* argv[])
{
    index_t order      = 2;
    index_t regularity = 1;
    index_t safety     = 1;
    bool    plot       = false;

    gsCmdLine cmd("Bottom-up nested hierarchical patch rule: embedded coordinates.");
    cmd.addInt   ("P","deg",   "Polynomial order of basis",              order);
    cmd.addInt   ("R","reg",   "Regularity of basis",                    regularity);
    cmd.addInt   ("s","safety","Extra new points per embedded level",    safety);
    cmd.addSwitch("plot",       "Export mesh and points to Paraview",      plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(regularity < order, "Require regularity < order.");

    // ── Build THB-spline basis ─────────────────────────────────────────────
    gsKnotVector<> kv(0.0, 1.0, 4, order+1, 1);  // 4 elems, max continuity
    gsTensorBSplineBasis<2,real_t> tbsb(kv, kv);
    gsTHBSplineBasis<2,real_t,true> thb(tbsb);

    thb.uniformRefine();
    thb.refineElements({1, 0, 0, 4, 4});
    thb.refineElements({2, 0, 0, 4, 4});

    const index_t L = thb.maxLevel();
    gsInfo << "THB basis: " << thb.size() << " functions, maxLevel=" << L << "\n";
    {
        std::map<index_t,index_t> cnt;
        for (index_t i = 0; i != thb.size(); ++i) ++cnt[thb.levelOf(i)];
        gsInfo << "  per-level counts:";
        for (std::map<index_t,index_t>::const_iterator it = cnt.begin(); it != cnt.end(); ++it)
            gsInfo << "  L" << it->first << ":" << it->second;
        gsInfo << "\n";
    }

    if (plot) { gsMesh<> m(thb); gsWriteParaview(m,"hrule_mesh"); }

    // Mass target: S^{2p}_r.  Plane stiffness safe target: S^{2p}_{r-1}.
    const index_t intDegM = 2 * order;
    const index_t intRegM = regularity;
    const index_t intDegK = 2 * order;
    const index_t intRegK = std::max<index_t>(regularity - 1, 0);

    // ── Reference: Gauss-Legendre on hierarchical leaf mesh ────────────────
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

    // ── Bottom-up nested hierarchical patch rules ─────────────────────────
    gsInfo << "\nConstructing bottom-up nested mass rule\n";
    gsBottomUpNestedPatchRule<real_t> hprM(thb, intDegM, intRegM, false, safety);
    gsInfo << "\nConstructing bottom-up nested stiffness-safe rule\n";
    gsBottomUpNestedPatchRule<real_t> hprK(thb, intDegK, intRegK, false, safety);

    gsInfo << "\n--- 1D node counts / nestedness, mass rule ---\n";
    for (short_t d = 0; d < 2; ++d)
    {
        gsInfo << "  d=" << d;
        for (index_t l = 0; l <= L; ++l)
            gsInfo << "  L" << l << ":" << hprM.numNodes1D(d,l);
        gsInfo << "  nested=" << (hprM.checkNested(d) ? "yes" : "NO") << "\n";
    }

    // ── Assembly with local evaluation-level dispatch ─────────────────────
    gsVector<real_t> IA(thb.size()); IA.setZero();
    gsMatrix<real_t> MA = gsMatrix<real_t>::Zero(thb.size(), thb.size());
    gsMatrix<real_t> KA = gsMatrix<real_t>::Zero(thb.size(), thb.size());
    index_t hprTotalPtsM = 0;
    index_t hprTotalPtsK = 0;
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

            ctr.col(0) = 0.5 * (low + upp);
            thb.active_into(ctr, actCtr);

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

            // Evaluate once on the local maximum level.
            hprM.mapTo(low, upp, evalLevel, evalLevel, nodesM, wEval);
            if (wEval.size() == 0) continue;
            hprTotalPtsM += wEval.size();
            thb.active_into(nodesM, acts);
            thb.eval_into  (nodesM, vals);

            hprK.mapTo(low, upp, evalLevel, evalLevel, nodesK, wEval);
            if (wEval.size() > 0)
            {
                hprTotalPtsK += wEval.size();
                thb.active_into(nodesK, actsK);
                thb.deriv_into (nodesK, derivs);
            }

            for (std::set<index_t>::const_iterator itLp = Lp_set.begin(); itLp != Lp_set.end(); ++itLp)
            {
                const index_t Lp = *itLp;

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

    // ── Results ────────────────────────────────────────────────────────────
    gsInfo << "\n========== Bottom-up Nested Hierarchical Patch Rule Results ==========\n";

    gsInfo << "\n--- Scalar integrals ---\n";
    gsInfo << "||I_hierarchical - I_ref||_inf = "
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
    gsInfo << "||M_hierarchical - M_ref||_inf = "
           << (MA - Mref).cwiseAbs().maxCoeff() << "\n";

    gsInfo << "\n--- Stiffness matrix ---\n";
    gsInfo << "||K_hierarchical - K_ref||_inf = "
           << (KA - Kref).cwiseAbs().maxCoeff() << "\n";

    gsInfo << "\n--- Quadrature point counts ---\n";
    gsInfo << "Gauss-Legendre (all elements):       " << gaussTotalPts << " pts\n";
    gsInfo << "Nested patch mass eval points:      " << hprTotalPtsM << " pts"
           << "  (ratio " << real_t(hprTotalPtsM)/gaussTotalPts << "x vs Gauss)\n";
    gsInfo << "Nested patch stiffness eval points: " << hprTotalPtsK << " pts"
           << "  (ratio " << real_t(hprTotalPtsK)/gaussTotalPts << "x vs Gauss)\n";

    // ── Plot locally evaluated mass-rule points ────────────────────────────
    if (plot)
    {
        gsMatrix<> allNodes(2, hprTotalPtsM);
        gsVector<> allWeights(hprTotalPtsM);
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
            for (index_t c = r; c < actCtr.rows(); ++c)
                evalLevel = std::max(evalLevel,
                    std::max(thb.repLevel(actCtr(r,0)), thb.repLevel(actCtr(c,0))));

            hprM.mapTo(low, upp, evalLevel, evalLevel, tmpNodes, tmpWeights);
            if (tmpWeights.size() == 0) continue;
            allNodes.block(0, cnt, 2, tmpWeights.size()) = tmpNodes;
            allWeights.segment(cnt, tmpWeights.size()) = tmpWeights;
            cnt += tmpWeights.size();
        }
        GISMO_ENSURE(cnt == hprTotalPtsM, "Mismatch in total point count.");

        gsInfo << "Exporting " << hprTotalPtsM
               << " locally evaluated nested points to Paraview...\n";
        gsMatrix<> X = allNodes.row(0);
        gsMatrix<> Y = allNodes.row(1);
        gsMatrix<> Z(1, hprTotalPtsM); Z.setZero();
        gsMatrix<> W = allWeights.transpose();
        gsWriteParaviewPoints(X, Y, Z, W, "hrule_bottomup_points");
    }

    return 0;
}
