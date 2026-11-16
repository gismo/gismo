/** @file hierarchical_patch_rule_1d_bottomup_vp.cpp

    @brief 1-D bottom-up nested embedded patch-rule experiment.

    This file is intentionally 1-D.  It is meant to debug the quadrature-rule
    construction before going back to the 2-D/THB assembler.

    Construction:

        X_0 subset X_1 subset ... subset X_L

    Level 0 is the ordinary Johannessen/gsPatchRule rule.  To construct
    X_{l+1}, all coordinates in X_l are kept fixed.  New coordinates are
    initialized from the ordinary level-(l+1) patch-rule nodes, with optional
    extra points added if needed.  Then a variable-projection nonlinear solve
    moves ONLY the new coordinates.  At every trial coordinate set, the weights
    are recomputed by the linear least-squares/minimum-norm moment solve

        A(X) w = b.

    Thus old coordinates are constrained, but all weights are free.

    Output focuses on 1-D diagnostics:
      - target-space dimension n_l,
      - ordinary vs nested point counts,
      - moment residuals,
      - span coverage,
      - nestedness.

    Author: H.M. Verhelst / experimental prototype
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
#include <vector>

using namespace gismo;

// =============================================================================
// Basic helpers
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

    // gsPatchRule requires an even target-space dimension.
    if (size % 2 == 1)
    {
        T maxSpan = T(0), mid = T(0);
        for (index_t i = 0; i + 1 < static_cast<index_t>(result.size()); ++i)
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

// =============================================================================
// Moment system utilities
// =============================================================================

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

template<class T>
T bestMomentResidual(const gsKnotVector<T>& kvTarget,
                     const gsVector<T>& pts,
                     gsVector<T>* weightsOut = nullptr)
{
    gsVector<T> w = solveMomentWeights(kvTarget, pts);
    if (weightsOut) *weightsOut = w;
    return momentResidual(kvTarget, pts, w);
}

// Ordinary Johannessen rule, used for level 0 and as initialization guidance.
template<class T>
gsVector<T> ordinaryPatchNodes(const gsKnotVector<T>& physKv,
                               index_t intDeg, index_t reg,
                               bool overInt)
{
    gsBSplineBasis<T> physBasis(physKv);
    gsPatchRule<T> rule(*physBasis.domain(), intDeg, reg, overInt);
    return sortUnique(rule.nodes(0));
}

// =============================================================================
// Diagnostics
// =============================================================================

template<class T>
void coverageDiagnostics(const gsKnotVector<T>& kvTarget,
                         const gsVector<T>& pts,
                         const std::string& prefix)
{
    gsMatrix<T> A;
    gsVector<T> I;
    buildMomentSystem(kvTarget, pts, A, I, static_cast<gsMatrix<T>*>(nullptr));

    index_t zeroRows = 0;
    T minRowNorm = std::numeric_limits<T>::max();
    for (index_t r = 0; r < A.rows(); ++r)
    {
        const T rn = A.row(r).norm();
        minRowNorm = std::min<T>(minRowNorm, rn);
        if (rn < T(1e-14)) ++zeroRows;
    }

    index_t zeroSpans = 0;
    index_t minPerSpan = std::numeric_limits<index_t>::max();
    index_t maxPerSpan = 0;
    for (index_t i = 0; i + 1 < static_cast<index_t>(kvTarget.size()); ++i)
    {
        const T a = kvTarget[i], b = kvTarget[i+1];
        if (!(b > a)) continue;
        index_t cnt = 0;
        for (index_t k = 0; k < pts.size(); ++k)
        {
            const bool lastSpan = (b >= kvTarget.last() - T(1e-14));
            if (pts[k] >= a - T(1e-14) && (pts[k] < b - T(1e-14) || (lastSpan && pts[k] <= b + T(1e-14))))
                ++cnt;
        }
        if (cnt == 0) ++zeroSpans;
        minPerSpan = std::min<index_t>(minPerSpan, cnt);
        maxPerSpan = std::max<index_t>(maxPerSpan, cnt);
    }
    if (minPerSpan == std::numeric_limits<index_t>::max()) minPerSpan = 0;

    gsInfo << prefix
           << " coverage: zero basis rows=" << zeroRows
           << ", min row norm=" << minRowNorm
           << ", zero knot spans=" << zeroSpans
           << ", min/max pts per nonzero span=" << minPerSpan << "/" << maxPerSpan
           << "\n";
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

// =============================================================================
// Initialization of new points
// =============================================================================

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

/**
 * @brief Initialize new points for X_{l+1} from the ordinary fine rule.
 *
 * By default, this keeps all ordinary level-(l+1) patch nodes that are not
 * already close to X_l.  If this still gives fewer than minNew points, extra
 * midpoint points are inserted in the largest gaps of X_l union candidates.
 */
template<class T>
gsVector<T> initializeNewPointsFromFineRule(const gsVector<T>& oldPtsIn,
                                            const gsVector<T>& ordinaryFineIn,
                                            index_t minNew,
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
    y = sortUnique(y, tol);

    if (y.size() < minNew)
    {
        gsVector<T> fixed = sortUnique(concatenate(oldPts, y), tol);
        gsVector<T> extra = insertLargestGapExtra(fixed, minNew - y.size(), a, b, tol);
        y = sortUnique(concatenate(y, extra), tol);
    }

    return y;
}

// =============================================================================
// Variable-projection embedded solve
// =============================================================================

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
EmbeddedSolveResult<T> constructEmbeddedNextLevelVP(const gsVector<T>& oldPtsIn,
                                                     const gsVector<T>& ordinaryFinePts,
                                                     const gsKnotVector<T>& kvTargetNext,
                                                     index_t safety,
                                                     T tol,
                                                     index_t maxIter)
{
    EmbeddedSolveResult<T> out;
    out.iterations = 0;
    out.converged = false;

    const T a = kvTargetNext.first();
    const T b = kvTargetNext.last();
    gsVector<T> oldPts = sortUnique(oldPtsIn);

    gsBSplineBasis<T> basis(kvTargetNext);
    const index_t n = basis.size();
    const index_t Kold = oldPts.size();

    index_t minNew = 1;
    if (n > Kold)
        minNew = static_cast<index_t>(std::ceil(T(n - Kold) / T(2))) + safety;
    else
        minNew = safety;
    minNew = std::max<index_t>(minNew, 1);

    gsVector<T> newPts = initializeNewPointsFromFineRule(oldPts, ordinaryFinePts, minNew, a, b);
    index_t numNew = newPts.size();
    gsVector<T> allPts = concatenate(oldPts, newPts);

    gsVector<T> weights;
    T res = bestMomentResidual(kvTargetNext, allPts, &weights);

    gsInfo << "    VP embedded solve: n=" << n
           << ", Kold=" << Kold
           << ", minNew=" << minNew
           << ", Knew=" << numNew
           << ", Ktotal=" << allPts.size()
           << ", initial best residual=" << res << "\n";
    coverageDiagnostics(kvTargetNext, allPts, "      initial");

    // If the fixed-node system already solves the moments, do not move points.
    if (res < tol)
    {
        out.points = sortUnique(allPts);
        out.weights = weights;
        out.residual = res;
        out.converged = true;
        return out;
    }

    for (index_t it = 0; it < maxIter; ++it)
    {
        allPts = concatenate(oldPts, newPts);

        gsMatrix<T> A, dA;
        gsVector<T> I;
        buildMomentSystem(kvTargetNext, allPts, A, I, &dA);

        // Variable projection: recompute optimal weights at current coordinates.
        weights = A.bdcSvd(gsEigen::ComputeThinU | gsEigen::ComputeThinV).solve(I);
        gsVector<T> F = A * weights - I;
        res = F.norm();

        if (res < tol)
        {
            gsInfo << "      converged at it=" << it << ", residual=" << res << "\n";
            out.points = sortUnique(allPts);
            out.weights = solveMomentWeights(kvTargetNext, out.points);
            out.residual = momentResidual(kvTargetNext, out.points, out.weights);
            out.iterations = it;
            out.converged = true;
            return out;
        }

        // Jacobian of variable-projection residual.  We project dA*w onto the
        // orthogonal complement of range(A) by g - A*A^+*g.
        //
        // Important G+Smo/Eigen compatibility note:
        // Do not instantiate gsEigen::BDCSVD<gsMatrix<T> > explicitly here.
        // In some G+Smo versions this triggers evaluator/const-cast template
        // errors.  Use the matrix member interface A.bdcSvd(...).solve(...),
        // which is the interface used elsewhere in this file and in G+Smo.
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

        // Gauss-Newton step for new coordinates only.
        gsVector<T> dy = J.bdcSvd(gsEigen::ComputeThinU | gsEigen::ComputeThinV).solve(-F);

        // Avoid enormous trial steps.
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

            // Mild Armijo-style decrease; pure decrease can stagnate on noise.
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

        if (!accepted)
        {
            gsInfo << "      WARNING: line search failed at it=" << it
                   << ", residual=" << res << "\n";
            break;
        }

        if (it % 10 == 0 || res < T(1e-8))
            gsInfo << "      it=" << it << " residual=" << res << "\n";
    }

    allPts = concatenate(oldPts, newPts);
    out.points = sortUnique(allPts);
    out.weights = solveMomentWeights(kvTargetNext, out.points);
    out.residual = momentResidual(kvTargetNext, out.points, out.weights);
    out.converged = (out.residual < tol);

    if (!out.converged)
        gsInfo << "      WARNING: final residual=" << out.residual
               << " after " << out.iterations << " iterations\n";
    else
        gsInfo << "      converged/final residual=" << out.residual << "\n";

    coverageDiagnostics(kvTargetNext, out.points, "      final");
    return out;
}

// =============================================================================
// Main 1-D experiment
// =============================================================================

int main(int argc, char* argv[])
{
    index_t order      = 2;
    index_t regularity = 1;
    index_t levels     = 2;
    index_t elements0  = 4;
    index_t safety     = 0;
    index_t maxIter    = 100;
    real_t  tol        = 1e-12;
    bool    overInt    = false;
    bool    stiffness  = false;
    std::string csv    = "nested_1d_points.csv";

    gsCmdLine cmd("1-D bottom-up nested embedded patch-rule experiment.");
    cmd.addInt   ("P", "deg",      "Polynomial degree of the physical spline basis", order);
    cmd.addInt   ("R", "reg",      "Physical regularity", regularity);
    cmd.addInt   ("L", "levels",   "Finest level", levels);
    cmd.addInt   ("E", "elements", "Number of level-0 elements", elements0);
    cmd.addInt   ("s", "safety",   "Extra new points beyond dimension-count minimum", safety);
    cmd.addInt   ("i", "maxIter",  "Maximum variable-projection iterations per level", maxIter);
    cmd.addReal  ("t", "tol",      "Moment residual tolerance", tol);
    cmd.addSwitch("",  "over",     "Use boundary over-integration", overInt);
    cmd.addSwitch("",  "stiffness","Use 1-D stiffness target S^{2p-2}_{r-1} instead of mass target S^{2p}_r", stiffness);
    cmd.addString("", "csv",       "CSV output file for nested points", csv);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    GISMO_ENSURE(regularity < order, "Require regularity < order.");
    GISMO_ENSURE(levels >= 0, "levels must be nonnegative.");

    index_t intDeg, intReg;
    if (stiffness)
    {
        intDeg = std::max<index_t>(2*(order-1), 1);
        intReg = std::max<index_t>(regularity-1, 0);
    }
    else
    {
        intDeg = 2 * order;
        intReg = regularity;
    }

    gsInfo << "\n========== 1-D nested embedded patch-rule experiment ==========" << "\n";
    gsInfo << "basis degree=" << order << ", regularity=" << regularity
           << ", levels=" << levels << ", level-0 elements=" << elements0 << "\n";
    gsInfo << "target degree=" << intDeg << ", target regularity=" << intReg
           << (stiffness ? "  [stiffness]\n" : "  [mass]\n");

    // Build physical knot vectors for levels 0..L by uniform dyadic refinement.
    std::vector<gsKnotVector<real_t> > physKv(levels + 1);
    physKv[0] = gsKnotVector<real_t>(0.0, 1.0, elements0, order + 1, 1);
    for (index_t l = 1; l <= levels; ++l)
    {
        physKv[l] = physKv[l-1];
        physKv[l].uniformRefine();
    }

    std::vector<gsKnotVector<real_t> > targetKv(levels + 1);
    std::vector<gsVector<real_t> > ordinaryNodes(levels + 1);
    for (index_t l = 0; l <= levels; ++l)
    {
        targetKv[l] = makeTargetKv(physKv[l], intDeg, intReg, overInt);
        ordinaryNodes[l] = ordinaryPatchNodes(physKv[l], intDeg, intReg, overInt);
        gsBSplineBasis<real_t> tb(targetKv[l]);
        gsVector<real_t> wOrd = solveMomentWeights(targetKv[l], ordinaryNodes[l]);
        const real_t rOrd = momentResidual(targetKv[l], ordinaryNodes[l], wOrd);
        gsInfo << "ordinary L" << l
               << ": target dim=" << tb.size()
               << ", nodes=" << ordinaryNodes[l].size()
               << ", residual=" << rOrd << "\n";
    }

    // Bottom-up embedded construction.
    std::vector<gsVector<real_t> > X(levels + 1);
    std::vector<gsVector<real_t> > Wself(levels + 1);
    std::vector<real_t> residual(levels + 1, 0.0);

    X[0] = ordinaryNodes[0];
    Wself[0] = solveMomentWeights(targetKv[0], X[0]);
    residual[0] = momentResidual(targetKv[0], X[0], Wself[0]);

    gsInfo << "\nembedded L0: nodes=" << X[0].size()
           << ", residual=" << residual[0] << "\n";
    coverageDiagnostics(targetKv[0], X[0], "  L0");

    for (index_t l = 0; l < levels; ++l)
    {
        gsInfo << "\nconstructing embedded L" << (l+1) << " from L" << l << "\n";
        EmbeddedSolveResult<real_t> r = constructEmbeddedNextLevelVP(
            X[l], ordinaryNodes[l+1], targetKv[l+1], safety, tol, maxIter);

        X[l+1] = r.points;
        Wself[l+1] = r.weights;
        residual[l+1] = r.residual;

        gsInfo << "embedded L" << (l+1)
               << ": nodes=" << X[l+1].size()
               << ", residual=" << residual[l+1]
               << ", nested=" << (checkNested(std::vector<gsVector<real_t> >(X.begin(), X.begin()+l+2)) ? "yes" : "NO")
               << "\n";
    }

    // Embedded weights diagnostics: integrate every target level ell using every eval level m>=ell.
    gsInfo << "\n--- Embedded weight residuals ||A_ell(X_m) w^{ell->m}-b_ell||_2 ---\n";
    for (index_t ell = 0; ell <= levels; ++ell)
    {
        for (index_t m = ell; m <= levels; ++m)
        {
            gsVector<real_t> w = solveMomentWeights(targetKv[ell], X[m]);
            const real_t r = momentResidual(targetKv[ell], X[m], w);
            gsInfo << "  ell=" << ell << " eval=" << m
                   << "  K=" << X[m].size()
                   << "  residual=" << r
                   << "  minW=" << w.minCoeff()
                   << "  maxW=" << w.maxCoeff()
                   << "\n";
        }
    }

    gsInfo << "\n--- Summary ---\n";
    for (index_t l = 0; l <= levels; ++l)
    {
        gsBSplineBasis<real_t> tb(targetKv[l]);
        gsInfo << "  L" << l
               << ": target dim=" << tb.size()
               << ", ordinary nodes=" << ordinaryNodes[l].size()
               << ", nested nodes=" << X[l].size()
               << ", self residual=" << residual[l]
               << "\n";
    }
    gsInfo << "  nested family = " << (checkNested(X) ? "yes" : "NO") << "\n";

    // Write CSV with points and self weights for plotting/debugging.
    std::ofstream os(csv.c_str());
    os << std::setprecision(17);
    os << "level,x,w_self\n";
    for (index_t l = 0; l <= levels; ++l)
    {
        gsVector<real_t> w = solveMomentWeights(targetKv[l], X[l]);
        for (index_t k = 0; k < X[l].size(); ++k)
            os << l << "," << X[l][k] << "," << w[k] << "\n";
    }
    os.close();
    gsInfo << "\nWrote nested points to " << csv << "\n";

    return 0;
}
