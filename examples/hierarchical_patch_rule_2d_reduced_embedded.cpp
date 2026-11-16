/** @file hierarchical_patch_rule_2d_reduced_embedded.cpp

    @brief 2-D THB experiment with a reduced embedded nested patch rule and
           zero-extended weights.

    Construction in each parametric direction:

        X_0 = ordinary Johannessen/gsPatchRule nodes on level 0.

        Given X_{l-1}, construct X_l = X_{l-1} union Y_l, where the old
        coordinates X_{l-1} are fixed.  The initial candidates Y_l are chosen
        as a reduced subset of the ordinary level-l patch-rule nodes by a
        projected greedy rank-revealing selection.  Then a variable-projection
        Gauss--Newton solve moves ONLY the new coordinates Y_l; all weights are
        recomputed freely by the linear moment fit A(X_l) w = b.

    This is the reduced version of the reduced-embedded rule.  It keeps exact
    nestedness X_0 subset X_1 subset ... subset X_L, but aims for point counts
    closer to the embedded lower bound

        |X_l| >= ceil((dim(T_l)+|X_{l-1}|)/2).

    In 2-D the 1-D rule is tensorized.  On every leaf element, we choose

        evalLevel = max repLevel(i) over active THB functions on the element,

    evaluate the local basis once at X_evalLevel x X_evalLevel, and assemble
    the main ZERO-EXTENDED policy:

        targetLevel = repLevel(i) for scalar integrals,
        targetLevel = max(repLevel(i),repLevel(j)) for matrix entries.

    If targetLevel < evalLevel, self weights of X_targetLevel are copied onto
    the matching nested nodes in X_evalLevel and all added nodes receive zero.

    This file is an experimental prototype; it prioritizes diagnostics over
    performance.
*/

#include <gismo.h>
#include <gsAssembler/gsPatchRule.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace gismo;

// =============================================================================
// 1-D helper functions
// =============================================================================

template<class T>
gsKnotVector<T> makeTargetKv(const gsKnotVector<T>& kv,
                             index_t intDeg, index_t reg,
                             bool overInt = false)
{
    gsKnotVector<T> result = kv;

    // Mirrors gsPatchRule::_init as closely as possible.
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

    // gsPatchRule wants an even target-space dimension.
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
        if (u.empty() || std::abs(*it - u.back()) > tol)
            u.push_back(*it);

    gsVector<T> y(u.size());
    for (index_t i = 0; i < static_cast<index_t>(u.size()); ++i) y[i] = u[i];
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
gsVector<T> solveMomentWeights(const gsKnotVector<T>& kvTarget,
                               const gsVector<T>& pts)
{
    gsMatrix<T> A;
    gsVector<T> I;
    buildMomentSystem(kvTarget, pts, A, I);
    return A.bdcSvd(gsEigen::ComputeThinU | gsEigen::ComputeThinV).solve(I);
}

template<class T>
T momentResidual(const gsKnotVector<T>& kvTarget,
                 const gsVector<T>& pts,
                 const gsVector<T>& weights)
{
    gsMatrix<T> A;
    gsVector<T> I;
    buildMomentSystem(kvTarget, pts, A, I);
    return (A * weights - I).norm();
}

template<class T>
gsVector<T> ordinaryPatchNodes(const gsKnotVector<T>& physKv,
                               index_t intDeg, index_t reg,
                               bool overInt)
{
    gsBSplineBasis<T> physBasis(physKv);
    gsPatchRule<T> rule(*physBasis.domain(), intDeg, reg, overInt);
    return sortUnique(rule.nodes(0));
}

template<class T>
bool checkNested(const std::vector<gsVector<T> >& X, T tol = T(1e-12))
{
    for (index_t l = 0; l + 1 < static_cast<index_t>(X.size()); ++l)
    {
        for (index_t i = 0; i < X[l].size(); ++i)
        {
            bool found = false;
            for (index_t j = 0; j < X[l+1].size(); ++j)
            {
                if (std::abs(X[l][i] - X[l+1][j]) < tol)
                {
                    found = true;
                    break;
                }
            }
            if (!found) return false;
        }
    }
    return true;
}


template<class T>
T minDistanceToSet(T x, const gsVector<T>& pts)
{
    T d = std::numeric_limits<T>::max();
    for (index_t i = 0; i < pts.size(); ++i)
        d = std::min<T>(d, std::abs(x - pts[i]));
    return d;
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

template<class T>
gsVector<T> candidateFinePoints(const gsVector<T>& oldPtsIn,
                                const gsVector<T>& ordinaryFineIn,
                                T a, T b,
                                T tol = T(1e-11))
{
    gsVector<T> oldPts = sortUnique(oldPtsIn, tol);
    gsVector<T> ordinaryFine = sortUnique(ordinaryFineIn, tol);
    std::vector<T> cand;
    for (index_t i = 0; i < ordinaryFine.size(); ++i)
    {
        const T x = ordinaryFine[i];
        if (x <= a + tol || x >= b - tol) continue;
        if (minDistanceToSet(x, oldPts) > tol)
            cand.push_back(x);
    }
    gsVector<T> y(cand.size());
    for (index_t i = 0; i < static_cast<index_t>(cand.size()); ++i) y[i] = cand[i];
    return sortUnique(y, tol);
}

template<class T>
gsVector<T> insertLargestGapExtra(const gsVector<T>& fixedPtsIn,
                                  index_t numExtra,
                                  T a, T b,
                                  T tol = T(1e-13))
{
    gsVector<T> fixedPts = sortUnique(fixedPtsIn, tol);
    std::vector<T> pts;
    for (index_t i = 0; i < fixedPts.size(); ++i) pts.push_back(fixedPts[i]);

    std::vector<T> extras;
    for (index_t k = 0; k < numExtra; ++k)
    {
        std::sort(pts.begin(), pts.end());
        T bestLeft = a, bestRight = b, bestLen = T(-1);

        if (!pts.empty())
        {
            if (pts.front() - a > bestLen)
            {
                bestLen = pts.front() - a;
                bestLeft = a;
                bestRight = pts.front();
            }
            for (index_t i = 0; i + 1 < static_cast<index_t>(pts.size()); ++i)
            {
                const T len = pts[i+1] - pts[i];
                if (len > bestLen)
                {
                    bestLen = len;
                    bestLeft = pts[i];
                    bestRight = pts[i+1];
                }
            }
            if (b - pts.back() > bestLen)
            {
                bestLen = b - pts.back();
                bestLeft = pts.back();
                bestRight = b;
            }
        }

        const T x = T(0.5) * (bestLeft + bestRight);
        pts.push_back(x);
        extras.push_back(x);
    }

    gsVector<T> y(extras.size());
    for (index_t i = 0; i < static_cast<index_t>(extras.size()); ++i) y[i] = extras[i];
    return sortUnique(y, tol);
}

// Modified Gram--Schmidt helper for greedy candidate selection.
template<class T>
bool addNormalizedColumn(const gsVector<T>& a,
                         std::vector<gsVector<T> >& Q,
                         T tol)
{
    gsVector<T> r = a;
    for (typename std::vector<gsVector<T> >::const_iterator it = Q.begin(); it != Q.end(); ++it)
        r -= (*it) * ((*it).dot(r));
    const T nrm = r.norm();
    if (nrm <= tol) return false;
    Q.push_back(r / nrm);
    return true;
}

/**
 * @brief Select a reduced subset of ordinary fine nodes.
 *
 * The old nodes are fixed.  We project candidate columns against the span of
 * A(old) and greedily add candidates with the largest remaining column norm.
 * This is a cheap rank-revealing selection; it avoids the previous expensive
 * remove-one-point-at-a-time strategy.
 */
template<class T>
gsVector<T> selectReducedNewPoints(const gsKnotVector<T>& kvTarget,
                                   const gsVector<T>& oldPts,
                                   const gsVector<T>& candidates,
                                   index_t numNew,
                                   T tol = T(1e-12))
{
    if (numNew <= 0 || candidates.size() == 0)
        return gsVector<T>(0);

    const index_t take = std::min<index_t>(numNew, candidates.size());

    gsMatrix<T> Aold, Acand;
    gsVector<T> Iold, Icand;
    buildMomentSystem(kvTarget, oldPts, Aold, Iold);
    buildMomentSystem(kvTarget, candidates, Acand, Icand);

    std::vector<gsVector<T> > Q;
    for (index_t j = 0; j < Aold.cols(); ++j)
        addNormalizedColumn<T>(Aold.col(j), Q, tol);

    std::vector<bool> used(candidates.size(), false);
    std::vector<index_t> selected;
    selected.reserve(take);

    for (index_t s = 0; s < take; ++s)
    {
        index_t best = -1;
        T bestNorm = T(-1);
        gsVector<T> bestResidual;

        for (index_t j = 0; j < candidates.size(); ++j)
        {
            if (used[j]) continue;
            gsVector<T> r = Acand.col(j);
            for (typename std::vector<gsVector<T> >::const_iterator it = Q.begin(); it != Q.end(); ++it)
                r -= (*it) * ((*it).dot(r));
            const T nrm = r.norm();
            if (nrm > bestNorm)
            {
                bestNorm = nrm;
                best = j;
                bestResidual = r;
            }
        }

        if (best < 0) break;
        used[best] = true;
        selected.push_back(best);

        const T nrm = bestResidual.norm();
        if (nrm > tol)
            Q.push_back(bestResidual / nrm);
        else
            addNormalizedColumn<T>(Acand.col(best), Q, tol);
    }

    std::vector<T> x;
    for (typename std::vector<index_t>::const_iterator it = selected.begin(); it != selected.end(); ++it)
        x.push_back(candidates[*it]);
    std::sort(x.begin(), x.end());

    gsVector<T> y(x.size());
    for (index_t i = 0; i < static_cast<index_t>(x.size()); ++i) y[i] = x[i];
    return sortUnique(y, tol);
}

template<class T>
struct EmbeddedSolveResult
{
    gsVector<T> points;
    gsVector<T> weights;
    T residual;
    index_t iterations;
    bool converged;
};

template<class T>
T bestMomentResidual(const gsKnotVector<T>& kvTarget,
                     const gsVector<T>& pts,
                     gsVector<T>* weightsOut = nullptr)
{
    gsVector<T> w = solveMomentWeights(kvTarget, pts);
    if (weightsOut) *weightsOut = w;
    return momentResidual(kvTarget, pts, w);
}

/**
 * @brief Variable-projection solve for X_next = oldPts union newPts.
 *
 * Old coordinates are fixed. New coordinates are moved. At every trial point
 * set, all weights are recomputed by the linear least-squares/moment fit.
 */
template<class T>
EmbeddedSolveResult<T> constructEmbeddedNextLevelVPFromInitial(
        const gsVector<T>& oldPtsIn,
        const gsVector<T>& newPtsInitialIn,
        const gsKnotVector<T>& kvTargetNext,
        T tol,
        index_t maxIter)
{
    EmbeddedSolveResult<T> out;
    out.iterations = 0;
    out.converged = false;

    const T a = kvTargetNext.first();
    const T b = kvTargetNext.last();
    gsVector<T> oldPts = sortUnique(oldPtsIn);
    gsVector<T> newPts = sortUnique(newPtsInitialIn);

    gsVector<T> allPts = sortUnique(concatenate(oldPts, newPts));
    gsVector<T> weights;
    T res = bestMomentResidual(kvTargetNext, allPts, &weights);

    if (res < tol)
    {
        out.points = allPts;
        out.weights = weights;
        out.residual = res;
        out.converged = true;
        return out;
    }

    const index_t Kold = oldPts.size();
    index_t numNew = newPts.size();
    if (numNew == 0)
    {
        out.points = allPts;
        out.weights = weights;
        out.residual = res;
        out.converged = false;
        return out;
    }

    for (index_t it = 0; it < maxIter; ++it)
    {
        allPts = concatenate(oldPts, newPts);

        gsMatrix<T> A, dA;
        gsVector<T> I;
        buildMomentSystem(kvTargetNext, allPts, A, I, &dA);

        weights = A.bdcSvd(gsEigen::ComputeThinU | gsEigen::ComputeThinV).solve(I);
        gsVector<T> F = A * weights - I;
        res = F.norm();

        if (res < tol)
        {
            out.points = sortUnique(allPts);
            out.weights = solveMomentWeights(kvTargetNext, out.points);
            out.residual = momentResidual(kvTargetNext, out.points, out.weights);
            out.iterations = it;
            out.converged = true;
            return out;
        }

        // Variable-projection Jacobian: projected derivative of A(y)w.
        gsMatrix<T> J(A.rows(), numNew);
        J.setZero();
        for (index_t s = 0; s < numNew; ++s)
        {
            const index_t col = Kold + s;
            gsVector<T> g = dA.col(col);
            g *= weights[col];
            const gsVector<T> z = A.bdcSvd(gsEigen::ComputeThinU | gsEigen::ComputeThinV).solve(g);
            J.col(s) = g - A * z;
        }

        gsVector<T> dy = J.bdcSvd(gsEigen::ComputeThinU | gsEigen::ComputeThinV).solve(-F);

        const T maxAbs = dy.size() ? dy.cwiseAbs().maxCoeff() : T(0);
        if (maxAbs > T(0.25) * (b - a))
            dy *= (T(0.25) * (b - a) / maxAbs);

        bool accepted = false;
        T alpha = T(1);
        for (index_t ls = 0; ls < 35; ++ls)
        {
            gsVector<T> yTrial = newPts + alpha * dy;
            gsVector<T> xTrial = concatenate(oldPts, yTrial);

            if (!validPointSet(xTrial, a, b))
            {
                alpha *= T(0.5);
                continue;
            }

            gsVector<T> wTrial;
            const T trialRes = bestMomentResidual(kvTargetNext, xTrial, &wTrial);
            if (trialRes < res * (T(1) - T(1e-4) * alpha) || trialRes < tol)
            {
                newPts = yTrial;
                weights = wTrial;
                res = trialRes;
                accepted = true;
                break;
            }
            alpha *= T(0.5);
        }

        out.iterations = it + 1;
        if (!accepted) break;
    }

    allPts = sortUnique(concatenate(oldPts, newPts));
    out.points = allPts;
    out.weights = solveMomentWeights(kvTargetNext, out.points);
    out.residual = momentResidual(kvTargetNext, out.points, out.weights);
    out.converged = (out.residual < tol);
    return out;
}

/**
 * @brief Construct one reduced embedded level.
 *
 * Tries increasing point counts.  If the reduced rule does not converge, the
 * caller can optionally fall back to the reduced embedded by setting maxExtra large
 * enough to include all ordinary fine candidates.
 */
template<class T>
EmbeddedSolveResult<T> constructReducedEmbeddedLevel(
        const gsVector<T>& oldPts,
        const gsVector<T>& ordinaryFinePts,
        const gsKnotVector<T>& kvTarget,
        index_t safety,
        index_t retryStep,
        index_t maxExtra,
        T tol,
        index_t maxIter)
{
    EmbeddedSolveResult<T> best;
    best.points = oldPts;
    best.weights = solveMomentWeights(kvTarget, oldPts);
    best.residual = momentResidual(kvTarget, oldPts, best.weights);
    best.iterations = 0;
    best.converged = false;

    const T a = kvTarget.first();
    const T b = kvTarget.last();
    gsBSplineBasis<T> tb(kvTarget);
    const index_t n = tb.size();
    const index_t Kold = oldPts.size();

    index_t minNew = 1;
    if (n > Kold)
        minNew = static_cast<index_t>(std::ceil(T(n - Kold) / T(2))) + safety;
    else
        minNew = std::max<index_t>(safety, 1);

    gsVector<T> cand = candidateFinePoints(oldPts, ordinaryFinePts, a, b);
    for (index_t extra = 0; extra <= maxExtra; extra += std::max<index_t>(retryStep,1))
    {
        index_t Knew = minNew + extra;
        gsVector<T> selected;
        if (cand.size() > 0)
            selected = selectReducedNewPoints(kvTarget, oldPts, cand, std::min<index_t>(Knew, cand.size()));
        else
            selected.resize(0);

        if (selected.size() < Knew)
        {
            gsVector<T> fixed = sortUnique(concatenate(oldPts, selected));
            gsVector<T> extraPts = insertLargestGapExtra(fixed, Knew - selected.size(), a, b);
            selected = sortUnique(concatenate(selected, extraPts));
        }

        gsVector<T> initialPts = sortUnique(concatenate(oldPts, selected));
        gsVector<T> wInit;
        const T initRes = bestMomentResidual(kvTarget, initialPts, &wInit);
        gsInfo << "      try Knew=" << selected.size()
               << " Ktotal=" << initialPts.size()
               << " initial residual=" << initRes << "\n";

        EmbeddedSolveResult<T> r = constructEmbeddedNextLevelVPFromInitial(
            oldPts, selected, kvTarget, tol, maxIter);

        gsInfo << "        final residual=" << r.residual
               << " iterations=" << r.iterations
               << " converged=" << (r.converged ? "yes" : "no") << "\n";

        if (r.residual < best.residual || best.points.size() == oldPts.size())
            best = r;

        if (r.converged)
            return r;
    }

    return best;
}

// =============================================================================
// Reduced embedded nested tensor-product rule for 2-D
// =============================================================================

template<class T>
class gsReducedEmbeddedPatchRule2D
{
public:
    gsReducedEmbeddedPatchRule2D(const gsHTensorBasis<2,T>& basis,
                             index_t intDeg,
                             index_t reg,
                             bool overInt = false,
                             const std::string& name = "rule",
                             index_t safety = 0,
                             index_t retryStep = 2,
                             index_t maxExtra = 12,
                             T tol = T(1e-12),
                             index_t maxIter = 100)
        : m_L(basis.maxLevel()),
          m_name(name)
    {
        m_end = basis.support().col(1);
        m_nodes.resize(2);
        m_targetKv.resize(2);
        m_weights.resize(2);
        m_zeroWeights.resize(2);

        for (short_t d = 0; d < 2; ++d)
        {
            m_nodes[d].resize(m_L + 1);
            m_targetKv[d].resize(m_L + 1);

            gsInfo << "\nConstructing reduced-embedded " << m_name
                   << " rule, direction d=" << d << "\n";

            std::vector<gsVector<T> > ordinaryNodes(m_L + 1);
            for (index_t l = 0; l <= m_L; ++l)
            {
                const gsKnotVector<T>* kvPhys =
                    static_cast<const gsKnotVector<T>*>(&basis.tensorLevel(l).knots(d));

                m_targetKv[d][l] = makeTargetKv(*kvPhys, intDeg, reg, overInt);
                ordinaryNodes[l] = ordinaryPatchNodes(*kvPhys, intDeg, reg, overInt);

                gsVector<T> wOrd = solveMomentWeights(m_targetKv[d][l], ordinaryNodes[l]);
                const T rOrd = momentResidual(m_targetKv[d][l], ordinaryNodes[l], wOrd);
                gsBSplineBasis<T> tb(m_targetKv[d][l]);
                gsInfo << "  ordinary L" << l
                       << ": target dim=" << tb.size()
                       << ", nodes=" << ordinaryNodes[l].size()
                       << ", residual=" << rOrd << "\n";
            }

            m_nodes[d][0] = ordinaryNodes[0];
            {
                gsVector<T> wSelf = solveMomentWeights(m_targetKv[d][0], m_nodes[d][0]);
                const T resSelf = momentResidual(m_targetKv[d][0], m_nodes[d][0], wSelf);
                gsBSplineBasis<T> tb(m_targetKv[d][0]);
                gsInfo << "  embedded L0"
                       << ": target dim=" << tb.size()
                       << ", nodes=" << m_nodes[d][0].size()
                       << ", self residual=" << resSelf << "\n";
            }

            for (index_t l = 1; l <= m_L; ++l)
            {
                gsBSplineBasis<T> tb(m_targetKv[d][l]);
                const index_t n = tb.size();
                const index_t Kold = m_nodes[d][l-1].size();
                const index_t minNew = (n > Kold)
                    ? static_cast<index_t>(std::ceil(T(n - Kold) / T(2))) + safety
                    : std::max<index_t>(safety, 1);

                gsInfo << "  constructing embedded L" << l
                       << " from L" << (l-1)
                       << ": n=" << n
                       << " Kold=" << Kold
                       << " minNew=" << minNew << "\n";

                EmbeddedSolveResult<T> r = constructReducedEmbeddedLevel(
                    m_nodes[d][l-1], ordinaryNodes[l], m_targetKv[d][l],
                    safety, retryStep, maxExtra, tol, maxIter);

                m_nodes[d][l] = r.points;
                gsInfo << "  embedded L" << l
                       << ": target dim=" << n
                       << ", ordinary nodes=" << ordinaryNodes[l].size()
                       << ", nested nodes=" << m_nodes[d][l].size()
                       << ", self residual=" << r.residual
                       << ", converged=" << (r.converged ? "yes" : "no")
                       << "\n";
            }

            gsInfo << "  nested=" << (checkNested(m_nodes[d]) ? "yes" : "NO") << "\n";

            // weights[targetLevel][evalLevel], valid for targetLevel <= evalLevel.
            // These are embedded fixed-node moment-fit weights w^{ell->m}.
            m_weights[d].resize(m_L + 1);
            for (index_t ell = 0; ell <= m_L; ++ell)
            {
                m_weights[d][ell].resize(m_L + 1);
                for (index_t m = ell; m <= m_L; ++m)
                {
                    m_weights[d][ell][m] = solveMomentWeights(m_targetKv[d][ell], m_nodes[d][m]);
                    const T res = momentResidual(m_targetKv[d][ell], m_nodes[d][m], m_weights[d][ell][m]);
                    gsInfo << "    embedded weight residual ell=" << ell
                           << " eval=" << m
                           << " K=" << m_nodes[d][m].size()
                           << " res=" << res
                           << " minW=" << m_weights[d][ell][m].minCoeff()
                           << " maxW=" << m_weights[d][ell][m].maxCoeff()
                           << "\n";
                }
            }

            // zeroWeights[targetLevel][evalLevel], valid for targetLevel <= evalLevel.
            // These keep the self weights on X_targetLevel and set all added
            // nodes in X_evalLevel \ X_targetLevel to zero.
            m_zeroWeights[d].resize(m_L + 1);
            for (index_t ell = 0; ell <= m_L; ++ell)
            {
                m_zeroWeights[d][ell].resize(m_L + 1);
                for (index_t m = ell; m <= m_L; ++m)
                {
                    m_zeroWeights[d][ell][m] = zeroExtend1DWeights(d, ell, m);
                    const T zres = momentResidual(m_targetKv[d][ell], m_nodes[d][m], m_zeroWeights[d][ell][m]);
                    gsInfo << "    zero-extended residual ell=" << ell
                           << " eval=" << m
                           << " K=" << m_nodes[d][m].size()
                           << " res=" << zres
                           << " minW=" << m_zeroWeights[d][ell][m].minCoeff()
                           << " maxW=" << m_zeroWeights[d][ell][m].maxCoeff()
                           << "\n";
                }
            }
        }
    }

    index_t maxLevel() const { return m_L; }

    index_t numNodes1D(short_t d, index_t evalLevel) const
    {
        return m_nodes[d][evalLevel].size();
    }

    /**
     * @brief Return tensor-product nodes and weights on one leaf element.
     *
     * @param targetLevel  target/moment level ell
     * @param evalLevel    coordinate/evaluation level m >= ell
     */
    void mapTo(const gsVector<T>& lower, const gsVector<T>& upper,
               index_t targetLevel, index_t evalLevel,
               gsMatrix<T>& nodes, gsVector<T>& weights) const
    {
        GISMO_ASSERT(targetLevel >= 0 && targetLevel <= m_L,
            "targetLevel " << targetLevel << " out of range");
        GISMO_ASSERT(evalLevel >= targetLevel && evalLevel <= m_L,
            "Need targetLevel <= evalLevel <= maxLevel");

        std::vector<std::vector<index_t> > idxs(2);
        for (short_t d = 0; d < 2; ++d)
            idxs[d] = _indicesInElement(d, evalLevel, lower[d], upper[d]);

        _tensorProduct(idxs, targetLevel, evalLevel, nodes, weights);
    }

    void getWeights(const gsVector<T>& lower, const gsVector<T>& upper,
                    index_t targetLevel, index_t evalLevel,
                    gsVector<T>& weights) const
    {
        gsMatrix<T> nodes;
        mapTo(lower, upper, targetLevel, evalLevel, nodes, weights);
    }

    /**
     * @brief Return tensor-product nodes at evalLevel and zero-extended
     *        self weights of targetLevel on those nodes.
     *
     * For targetLevel < evalLevel, weights are nonzero only at nodes that
     * already belong to X_targetLevel.  The added nodes in
     * X_evalLevel \ X_targetLevel get zero weight.
     */
    void mapToZeroExtended(const gsVector<T>& lower, const gsVector<T>& upper,
                           index_t targetLevel, index_t evalLevel,
                           gsMatrix<T>& nodes, gsVector<T>& weights) const
    {
        GISMO_ASSERT(targetLevel >= 0 && targetLevel <= m_L,
            "targetLevel " << targetLevel << " out of range");
        GISMO_ASSERT(evalLevel >= targetLevel && evalLevel <= m_L,
            "Need targetLevel <= evalLevel <= maxLevel");

        std::vector<std::vector<index_t> > idxs(2);
        for (short_t d = 0; d < 2; ++d)
            idxs[d] = _indicesInElement(d, evalLevel, lower[d], upper[d]);

        _tensorProductZeroExtended(idxs, targetLevel, evalLevel, nodes, weights);
    }

    void getZeroExtendedWeights(const gsVector<T>& lower, const gsVector<T>& upper,
                                index_t targetLevel, index_t evalLevel,
                                gsVector<T>& weights) const
    {
        gsMatrix<T> nodes;
        mapToZeroExtended(lower, upper, targetLevel, evalLevel, nodes, weights);
    }

private:
    gsVector<T> zeroExtend1DWeights(short_t d, index_t targetLevel, index_t evalLevel) const
    {
        GISMO_ASSERT(targetLevel <= evalLevel, "targetLevel must be <= evalLevel");

        const gsVector<T>& xTarget = m_nodes[d][targetLevel];
        const gsVector<T>& xEval   = m_nodes[d][evalLevel];
        const gsVector<T>& wTarget = m_weights[d][targetLevel][targetLevel];

        gsVector<T> w = gsVector<T>::Zero(xEval.size());
        const T tol = T(1e-11);

        for (index_t ke = 0; ke < xEval.size(); ++ke)
        {
            for (index_t kt = 0; kt < xTarget.size(); ++kt)
            {
                if (std::abs(xEval[ke] - xTarget[kt]) < tol)
                {
                    w[ke] = wTarget[kt];
                    break;
                }
            }
        }
        return w;
    }

    std::vector<index_t> _indicesInElement(short_t d, index_t evalLevel,
                                           T lower, T upper) const
    {
        std::vector<index_t> idx;
        const bool isLast = (upper >= m_end[d] - T(1e-14));
        const gsVector<T>& x = m_nodes[d][evalLevel];

        for (index_t k = 0; k < x.size(); ++k)
        {
            const bool inside = (x[k] >= lower - T(1e-14)) &&
                                (x[k] <  upper - T(1e-14) ||
                                 (isLast && x[k] <= upper + T(1e-14)));
            if (inside) idx.push_back(k);
        }
        return idx;
    }

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
                nodes(0,cnt) = m_nodes[0][evalLevel][kx];
                nodes(1,cnt) = m_nodes[1][evalLevel][ky];
                weights[cnt] = wx[kx] * wy[ky];
            }
        }
    }

    void _tensorProductZeroExtended(const std::vector<std::vector<index_t> >& idxs,
                                    index_t targetLevel, index_t evalLevel,
                                    gsMatrix<T>& nodes, gsVector<T>& weights) const
    {
        const index_t nx = static_cast<index_t>(idxs[0].size());
        const index_t ny = static_cast<index_t>(idxs[1].size());
        const index_t sz = nx * ny;

        nodes.resize(2, sz);
        weights.resize(sz);
        if (sz == 0) return;

        const gsVector<T>& wx = m_zeroWeights[0][targetLevel][evalLevel];
        const gsVector<T>& wy = m_zeroWeights[1][targetLevel][evalLevel];

        index_t cnt = 0;
        for (index_t jy = 0; jy < ny; ++jy)
        {
            const index_t ky = idxs[1][jy];
            for (index_t jx = 0; jx < nx; ++jx, ++cnt)
            {
                const index_t kx = idxs[0][jx];
                nodes(0,cnt) = m_nodes[0][evalLevel][kx];
                nodes(1,cnt) = m_nodes[1][evalLevel][ky];
                weights[cnt] = wx[kx] * wy[ky];
            }
        }
    }

private:
    index_t m_L;
    std::string m_name;
    gsVector<T> m_end;

    // [d][evalLevel]
    std::vector<std::vector<gsVector<T> > > m_nodes;

    // [d][targetLevel]
    std::vector<std::vector<gsKnotVector<T> > > m_targetKv;

    // [d][targetLevel][evalLevel], embedded fixed-node moment-fit weights
    std::vector<std::vector<std::vector<gsVector<T> > > > m_weights;

    // [d][targetLevel][evalLevel], self weights of targetLevel zero-extended to evalLevel
    std::vector<std::vector<std::vector<gsVector<T> > > > m_zeroWeights;
};

// =============================================================================
// Main experiment
// =============================================================================

int main(int argc, char* argv[])
{
    index_t order      = 2;
    index_t regularity = 1;
    bool    plot       = false;
    bool    overInt    = false;
    std::string csv    = "reduced_embedded_2d_points.csv";
    index_t safety     = 0;
    index_t retryStep  = 2;
    index_t maxExtra   = 12;
    index_t maxIter    = 100;
    real_t  tol        = 1e-12;

    gsCmdLine cmd("2-D THB reduced embedded nested patch-rule experiment with zero-extended weights.");
    cmd.addInt   ("P", "deg", "Polynomial degree of physical basis", order);
    cmd.addInt   ("R", "reg", "Physical regularity", regularity);
    cmd.addSwitch("",  "plot",     "Export point cloud to CSV/Paraview", plot);
    cmd.addSwitch("",  "over",     "Use boundary over-integration", overInt);
    cmd.addString("", "csv",       "CSV output file for evaluation points", csv);
    cmd.addInt   ("s", "safety", "Extra new points beyond embedded dimension-count minimum", safety);
    cmd.addInt   ("q", "retryStep", "Retry increment for extra new points", retryStep);
    cmd.addInt   ("x", "maxExtra", "Maximum extra new points to try", maxExtra);
    cmd.addInt   ("i", "maxIter", "Maximum variable-projection iterations per level", maxIter);
    cmd.addReal  ("t", "tol", "Moment residual tolerance", tol);
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
        for (std::map<index_t,index_t>::const_iterator it = cnt.begin(); it != cnt.end(); ++it)
            gsInfo << "  L" << it->first << ":" << it->second;
        gsInfo << "\n";
    }

    if (plot)
    {
        gsMesh<> m(thb);
        gsWriteParaview(m,"reduced_embedded_mesh");
    }

    // ── Target spaces ─────────────────────────────────────────────────────
    // Mass/scalar: S^{2p}_{r}.  For max-cont p=2, r=1 -> S^4_C1.
    const index_t intDegM = 2 * order;
    const index_t intRegM = regularity;

    // Safe 2-D stiffness target: S^{2p}_{r-1} in both directions.
    // This is more expensive than mixed derivative tensor rules, but avoids
    // the previous under-integration of the transverse direction.
    const index_t intDegK = 2 * order;
    const index_t intRegK = std::max<index_t>(regularity - 1, 0);

    gsInfo << "\nMass target:      degree=" << intDegM << " regularity=" << intRegM << "\n";
    gsInfo << "Stiffness target: degree=" << intDegK << " regularity=" << intRegK
           << "  [safe tensor target]\n";

    gsReducedEmbeddedPatchRule2D<real_t> ruleM(thb, intDegM, intRegM, overInt, "mass",
                                                safety, retryStep, maxExtra, tol, maxIter);
    gsReducedEmbeddedPatchRule2D<real_t> ruleK(thb, intDegK, intRegK, overInt, "stiffness-safe",
                                                safety, retryStep, maxExtra, tol, maxIter);

    // ── Reference: ordinary Gauss on all leaf elements ─────────────────────
    gsOptionList gOpts;
    gOpts.addInt   ("quRule", "", gsQuadrature::GaussLegendre);
    gOpts.addReal  ("quA",    "", 1.0);
    gOpts.addInt   ("quB",    "", 1);
    gOpts.addSwitch("overInt", "", false);
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
            {
                for (index_t r = 0; r < acts.rows(); ++r)
                {
                    const index_t i = acts(r,0);
                    Iref[i] += weights[k] * vals(r,k);
                    for (index_t c = r; c < acts.rows(); ++c)
                    {
                        const index_t j = acts(c,0);
                        const real_t sm = weights[k] * vals(r,k) * vals(c,k);
                        Mref(i,j) += sm; if (i != j) Mref(j,i) += sm;

                        const real_t sk = weights[k] *
                            (derivs(2*r+0,k) * derivs(2*c+0,k) +
                             derivs(2*r+1,k) * derivs(2*c+1,k));
                        Kref(i,j) += sk; if (i != j) Kref(j,i) += sk;
                    }
                }
            }
        }
    }

    // ── Reduced embedded assembly: zero-extended, embedded, and eval weights ─
    gsVector<real_t> Izero(thb.size()); Izero.setZero();
    gsVector<real_t> Idisp(thb.size()); Idisp.setZero();
    gsVector<real_t> Ieval(thb.size()); Ieval.setZero();
    gsMatrix<real_t> Mzero = gsMatrix<real_t>::Zero(thb.size(), thb.size());
    gsMatrix<real_t> Mdisp = gsMatrix<real_t>::Zero(thb.size(), thb.size());
    gsMatrix<real_t> Meval = gsMatrix<real_t>::Zero(thb.size(), thb.size());
    gsMatrix<real_t> Kzero = gsMatrix<real_t>::Zero(thb.size(), thb.size());
    gsMatrix<real_t> Kdisp = gsMatrix<real_t>::Zero(thb.size(), thb.size());
    gsMatrix<real_t> Keval = gsMatrix<real_t>::Zero(thb.size(), thb.size());

    index_t massEvalPts = 0;
    index_t stiffEvalPts = 0;
    std::vector<index_t> elemEvalCount(L+1, 0);

    {
        gsMatrix<real_t> ctr(2,1);
        gsMatrix<index_t> actCtr;

        gsMatrix<real_t> nodesM, nodesK;
        gsVector<real_t> wEvalM, wLpM, wZeroM, wEvalK, wLpK, wZeroK;
        gsMatrix<index_t> actsM, actsK;
        gsMatrix<real_t> vals, derivs;

        for (auto& elem : thb.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,2>*>(&elem);
            GISMO_ASSERT(helem, "Element iterator is not gsHDomainIterator.");
            const gsVector<real_t>& low = helem->lowerCorner();
            const gsVector<real_t>& upp = helem->upperCorner();

            // Local evaluation level: max representation level of active functions.
            ctr.col(0) = real_t(0.5) * (low + upp);
            thb.active_into(ctr, actCtr);

            index_t evalLevel = 0;
            std::set<index_t> Lp_set;
            for (index_t r = 0; r < actCtr.rows(); ++r)
            {
                const index_t i = actCtr(r,0);
                evalLevel = std::max<index_t>(evalLevel, thb.repLevel(i));
            }
            for (index_t r = 0; r < actCtr.rows(); ++r)
                for (index_t c = r; c < actCtr.rows(); ++c)
                    Lp_set.insert(std::max(thb.repLevel(actCtr(r,0)),
                                           thb.repLevel(actCtr(c,0))));

            GISMO_ASSERT(evalLevel <= L, "evalLevel out of range");
            ++elemEvalCount[evalLevel];

            // ── Mass/scalar nodes, evaluated once at evalLevel ─────────────
            ruleM.mapTo(low, upp, evalLevel, evalLevel, nodesM, wEvalM);
            if (wEvalM.size() > 0)
            {
                massEvalPts += wEvalM.size();
                thb.active_into(nodesM, actsM);
                thb.eval_into  (nodesM, vals);

                // eval-weights policy: use targetLevel = evalLevel for everything.
                for (index_t r = 0; r < actsM.rows(); ++r)
                {
                    const index_t i = actsM(r,0);
                    for (index_t k = 0; k < wEvalM.size(); ++k)
                        Ieval[i] += wEvalM[k] * vals(r,k);
                }
                for (index_t r = 0; r < actsM.rows(); ++r)
                for (index_t c = r; c < actsM.rows(); ++c)
                {
                    const index_t i = actsM(r,0), j = actsM(c,0);
                    real_t sm = 0;
                    for (index_t k = 0; k < wEvalM.size(); ++k)
                        sm += wEvalM[k] * vals(r,k) * vals(c,k);
                    Meval(i,j) += sm; if (i != j) Meval(j,i) += sm;
                }

                // zero-extended policy: use the self weights of targetLevel
                // and zero out all added nodes in X_evalLevel \ X_targetLevel.
                for (index_t Lp = 0; Lp <= evalLevel; ++Lp)
                {
                    ruleM.getZeroExtendedWeights(low, upp, Lp, evalLevel, wZeroM);
                    if (wZeroM.size() == 0) continue;

                    for (index_t r = 0; r < actsM.rows(); ++r)
                    {
                        const index_t i = actsM(r,0);
                        if (thb.repLevel(i) != Lp) continue;
                        for (index_t k = 0; k < wZeroM.size(); ++k)
                            Izero[i] += wZeroM[k] * vals(r,k);
                    }

                    for (index_t r = 0; r < actsM.rows(); ++r)
                    for (index_t c = r; c < actsM.rows(); ++c)
                    {
                        const index_t i = actsM(r,0), j = actsM(c,0);
                        if (std::max(thb.repLevel(i), thb.repLevel(j)) != Lp) continue;
                        real_t sm = 0;
                        for (index_t k = 0; k < wZeroM.size(); ++k)
                            sm += wZeroM[k] * vals(r,k) * vals(c,k);
                        Mzero(i,j) += sm; if (i != j) Mzero(j,i) += sm;
                    }
                }

                // embedded-fit policy: solve weights w^{Lp->evalLevel} on the
                // whole eval point set.  Kept only as a diagnostic.
                for (index_t Lp = 0; Lp <= evalLevel; ++Lp)
                {
                    ruleM.getWeights(low, upp, Lp, evalLevel, wLpM);
                    if (wLpM.size() == 0) continue;

                    for (index_t r = 0; r < actsM.rows(); ++r)
                    {
                        const index_t i = actsM(r,0);
                        if (thb.repLevel(i) != Lp) continue;
                        for (index_t k = 0; k < wLpM.size(); ++k)
                            Idisp[i] += wLpM[k] * vals(r,k);
                    }

                    for (index_t r = 0; r < actsM.rows(); ++r)
                    for (index_t c = r; c < actsM.rows(); ++c)
                    {
                        const index_t i = actsM(r,0), j = actsM(c,0);
                        if (std::max(thb.repLevel(i), thb.repLevel(j)) != Lp) continue;
                        real_t sm = 0;
                        for (index_t k = 0; k < wLpM.size(); ++k)
                            sm += wLpM[k] * vals(r,k) * vals(c,k);
                        Mdisp(i,j) += sm; if (i != j) Mdisp(j,i) += sm;
                    }
                }
            }

            // ── Stiffness-safe nodes, evaluated once at evalLevel ──────────
            ruleK.mapTo(low, upp, evalLevel, evalLevel, nodesK, wEvalK);
            if (wEvalK.size() > 0)
            {
                stiffEvalPts += wEvalK.size();
                thb.active_into(nodesK, actsK);
                thb.deriv_into (nodesK, derivs);

                // eval-weights policy.
                for (index_t r = 0; r < actsK.rows(); ++r)
                for (index_t c = r; c < actsK.rows(); ++c)
                {
                    const index_t i = actsK(r,0), j = actsK(c,0);
                    real_t sk = 0;
                    for (index_t k = 0; k < wEvalK.size(); ++k)
                        sk += wEvalK[k] *
                              (derivs(2*r+0,k) * derivs(2*c+0,k) +
                               derivs(2*r+1,k) * derivs(2*c+1,k));
                    Keval(i,j) += sk; if (i != j) Keval(j,i) += sk;
                }

                // zero-extended policy.
                for (index_t Lp = 0; Lp <= evalLevel; ++Lp)
                {
                    ruleK.getZeroExtendedWeights(low, upp, Lp, evalLevel, wZeroK);
                    if (wZeroK.size() == 0) continue;

                    for (index_t r = 0; r < actsK.rows(); ++r)
                    for (index_t c = r; c < actsK.rows(); ++c)
                    {
                        const index_t i = actsK(r,0), j = actsK(c,0);
                        if (std::max(thb.repLevel(i), thb.repLevel(j)) != Lp) continue;
                        real_t sk = 0;
                        for (index_t k = 0; k < wZeroK.size(); ++k)
                            sk += wZeroK[k] *
                                  (derivs(2*r+0,k) * derivs(2*c+0,k) +
                                   derivs(2*r+1,k) * derivs(2*c+1,k));
                        Kzero(i,j) += sk; if (i != j) Kzero(j,i) += sk;
                    }
                }

                // embedded-fit policy, kept only as a diagnostic.
                for (index_t Lp = 0; Lp <= evalLevel; ++Lp)
                {
                    ruleK.getWeights(low, upp, Lp, evalLevel, wLpK);
                    if (wLpK.size() == 0) continue;

                    for (index_t r = 0; r < actsK.rows(); ++r)
                    for (index_t c = r; c < actsK.rows(); ++c)
                    {
                        const index_t i = actsK(r,0), j = actsK(c,0);
                        if (std::max(thb.repLevel(i), thb.repLevel(j)) != Lp) continue;
                        real_t sk = 0;
                        for (index_t k = 0; k < wLpK.size(); ++k)
                            sk += wLpK[k] *
                                  (derivs(2*r+0,k) * derivs(2*c+0,k) +
                                   derivs(2*r+1,k) * derivs(2*c+1,k));
                        Kdisp(i,j) += sk; if (i != j) Kdisp(j,i) += sk;
                    }
                }
            }
        }
    }

    // ── Results ────────────────────────────────────────────────────────────
    gsInfo << "\n========== 2-D Reduced Embedded Hierarchical Patch Rule Results ==========" << "\n";

    gsInfo << "\n--- Element evaluation-level counts ---\n";
    for (index_t l = 0; l <= L; ++l)
        gsInfo << "  evalLevel L" << l << ": " << elemEvalCount[l] << " elements\n";

    gsInfo << "\n--- Scalar integrals ---\n";
    gsInfo << "zero-extended:    ||I - I_ref||_inf = "
           << (Izero - Iref).cwiseAbs().maxCoeff() << "\n";
    gsInfo << "dispatch weights: ||I - I_ref||_inf = "
           << (Idisp - Iref).cwiseAbs().maxCoeff() << "\n";
    gsInfo << "eval weights:     ||I - I_ref||_inf = "
           << (Ieval - Iref).cwiseAbs().maxCoeff() << "\n";

    gsInfo << "\n--- Mass matrix ---\n";
    gsInfo << "zero-extended:    ||M - M_ref||_inf = "
           << (Mzero - Mref).cwiseAbs().maxCoeff() << "\n";
    gsInfo << "dispatch weights: ||M - M_ref||_inf = "
           << (Mdisp - Mref).cwiseAbs().maxCoeff() << "\n";
    gsInfo << "eval weights:     ||M - M_ref||_inf = "
           << (Meval - Mref).cwiseAbs().maxCoeff() << "\n";

    gsInfo << "\n--- Stiffness matrix ---\n";
    gsInfo << "zero-extended:    ||K - K_ref||_inf = "
           << (Kzero - Kref).cwiseAbs().maxCoeff() << "\n";
    gsInfo << "dispatch weights: ||K - K_ref||_inf = "
           << (Kdisp - Kref).cwiseAbs().maxCoeff() << "\n";
    gsInfo << "eval weights:     ||K - K_ref||_inf = "
           << (Keval - Kref).cwiseAbs().maxCoeff() << "\n";

    gsInfo << "\n--- Quadrature point counts ---\n";
    gsInfo << "Gauss-Legendre all leaf elements: " << gaussTotalPts << " pts\n";
    gsInfo << "Reduced embedded mass eval points:   " << massEvalPts
           << " pts  (ratio " << real_t(massEvalPts)/gaussTotalPts << "x vs Gauss)\n";
    gsInfo << "Reduced embedded stiff eval points:  " << stiffEvalPts
           << " pts  (ratio " << real_t(stiffEvalPts)/gaussTotalPts << "x vs Gauss)\n";

    // CSV export of actual mass evaluation point cloud, with evalLevel per point.
    if (plot)
    {
        std::ofstream os(csv.c_str());
        os << std::setprecision(17);
        os << "x,y,evalLevel\n";

        gsMatrix<real_t> ctr(2,1);
        gsMatrix<index_t> actCtr;
        gsMatrix<real_t> nodes;
        gsVector<real_t> weights;

        for (auto& elem : thb.domain()->allElements())
        {
            auto helem = dynamic_cast<gsHDomainIterator<real_t,2>*>(&elem);
            const gsVector<real_t>& low = helem->lowerCorner();
            const gsVector<real_t>& upp = helem->upperCorner();
            ctr.col(0) = real_t(0.5) * (low + upp);
            thb.active_into(ctr, actCtr);

            index_t evalLevel = 0;
            for (index_t r = 0; r < actCtr.rows(); ++r)
                evalLevel = std::max<index_t>(evalLevel, thb.repLevel(actCtr(r,0)));

            ruleM.mapTo(low, upp, evalLevel, evalLevel, nodes, weights);
            for (index_t k = 0; k < weights.size(); ++k)
                os << nodes(0,k) << "," << nodes(1,k) << "," << evalLevel << "\n";
        }
        os.close();
        gsInfo << "\nWrote mass evaluation point cloud to " << csv << "\n";
    }

    return 0;
}
