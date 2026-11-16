/** @file hierarchical_patch_rule_2d_zero_extended.cpp

    @brief 2-D THB experiment with nested union patch rules and zero-extended weights.

    Construction in each parametric direction:

        X_0 = ordinary Johannessen/gsPatchRule nodes on level 0,
        X_l = sort_unique( X_{l-1} union ordinary Johannessen nodes on level l ).

    Hence

        X_0 subset X_1 subset ... subset X_L.

    For every target level ell and evaluation level m >= ell, we compute fixed-node
    moment-fitting weights

        A_ell(X_m) w^{ell->m} = b_ell.

    In 2-D the rule is tensorized.  On every leaf element, we choose

        evalLevel = max repLevel(i) over active THB functions on the element.

    We then evaluate basis values / derivatives ONCE at the tensor product nodes
    X_evalLevel x X_evalLevel.

    The main policy tested here is ZERO-EXTENDED WEIGHTS:

        targetLevel = repLevel(i)                 for scalar integrals,
        targetLevel = max(repLevel(i),repLevel(j)) for matrix entries.

    If targetLevel < evalLevel, we keep the self weights of X_targetLevel at
    the nodes that belong to X_targetLevel, and set the weights of the added
    nodes X_evalLevel \ X_targetLevel to zero.  Thus a level-l contribution is
    integrated by the same level-l global patch rule on every element, while we
    still evaluate once at the finest representation level required locally.

    For diagnostics, the file also keeps the previous two policies:
      - embedded-fit weights w^{ell->m}, and
      - eval-level weights w^{m->m} for everything.

    This is an experimental prototype.  It prioritizes transparent diagnostics
    over performance.
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
                       gsVector<T>& I)
{
    gsBSplineBasis<T> basis(kvTarget);
    const index_t n = basis.size();
    const index_t K = pts.size();
    const index_t p = kvTarget.degree();

    A = gsMatrix<T>::Zero(n, K);

    gsMatrix<T> ptsRow = pts.transpose();
    gsMatrix<index_t> acts;
    gsMatrix<T> vals;
    basis.active_into(ptsRow, acts);
    basis.eval_into  (ptsRow, vals);

    for (index_t k = 0; k < K; ++k)
        for (index_t r = 0; r < acts.rows(); ++r)
            A(acts(r,k), k) = vals(r,k);

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

// =============================================================================
// Nested union tensor-product rule for 2-D
// =============================================================================

template<class T>
class gsNestedUnionPatchRule2D
{
public:
    gsNestedUnionPatchRule2D(const gsHTensorBasis<2,T>& basis,
                             index_t intDeg,
                             index_t reg,
                             bool overInt = false,
                             const std::string& name = "rule")
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

            gsInfo << "\nConstructing nested-union " << m_name
                   << " rule, direction d=" << d << "\n";

            for (index_t l = 0; l <= m_L; ++l)
            {
                const gsKnotVector<T>* kvPhys =
                    static_cast<const gsKnotVector<T>*>(&basis.tensorLevel(l).knots(d));

                m_targetKv[d][l] = makeTargetKv(*kvPhys, intDeg, reg, overInt);
                gsVector<T> ordinary = ordinaryPatchNodes(*kvPhys, intDeg, reg, overInt);

                if (l == 0)
                    m_nodes[d][l] = ordinary;
                else
                    m_nodes[d][l] = sortUnique(concatenate(m_nodes[d][l-1], ordinary));

                gsVector<T> wSelf = solveMomentWeights(m_targetKv[d][l], m_nodes[d][l]);
                const T resSelf = momentResidual(m_targetKv[d][l], m_nodes[d][l], wSelf);

                gsBSplineBasis<T> tb(m_targetKv[d][l]);
                gsInfo << "  L" << l
                       << ": target dim=" << tb.size()
                       << ", ordinary nodes=" << ordinary.size()
                       << ", nested nodes=" << m_nodes[d][l].size()
                       << ", self residual=" << resSelf << "\n";
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
    std::string csv    = "nested_2d_points.csv";

    gsCmdLine cmd("2-D THB nested-union patch-rule experiment with zero-extended weights.");
    cmd.addInt   ("P", "deg", "Polynomial degree of physical basis", order);
    cmd.addInt   ("R", "reg", "Physical regularity", regularity);
    cmd.addSwitch("",  "plot",     "Export point cloud to CSV/Paraview", plot);
    cmd.addSwitch("",  "over",     "Use boundary over-integration", overInt);
    cmd.addString("", "csv",       "CSV output file for evaluation points", csv);
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
        gsWriteParaview(m,"nested_union_mesh");
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

    gsNestedUnionPatchRule2D<real_t> ruleM(thb, intDegM, intRegM, overInt, "mass");
    gsNestedUnionPatchRule2D<real_t> ruleK(thb, intDegK, intRegK, overInt, "stiffness-safe");

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

    // ── Nested-union assembly: zero-extended, embedded, and eval weights ─
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
    gsInfo << "\n========== 2-D Nested-Union Hierarchical Patch Rule Results ==========" << "\n";

    gsInfo << "\n--- Element evaluation-level counts ---\n";
    for (index_t l = 0; l <= L; ++l)
        gsInfo << "  evalLevel L" << l << ": " << elemEvalCount[l] << " elements\n";

    gsInfo << "\n--- Scalar integrals ---\n";
    gsInfo << "dispatch weights: ||I - I_ref||_inf = "
           << (Idisp - Iref).cwiseAbs().maxCoeff() << "\n";
    gsInfo << "eval weights:     ||I - I_ref||_inf = "
           << (Ieval - Iref).cwiseAbs().maxCoeff() << "\n";

    gsInfo << "\n--- Mass matrix ---\n";
    gsInfo << "dispatch weights: ||M - M_ref||_inf = "
           << (Mdisp - Mref).cwiseAbs().maxCoeff() << "\n";
    gsInfo << "eval weights:     ||M - M_ref||_inf = "
           << (Meval - Mref).cwiseAbs().maxCoeff() << "\n";

    gsInfo << "\n--- Stiffness matrix ---\n";
    gsInfo << "dispatch weights: ||K - K_ref||_inf = "
           << (Kdisp - Kref).cwiseAbs().maxCoeff() << "\n";
    gsInfo << "eval weights:     ||K - K_ref||_inf = "
           << (Keval - Kref).cwiseAbs().maxCoeff() << "\n";

    gsInfo << "\n--- Quadrature point counts ---\n";
    gsInfo << "Gauss-Legendre all leaf elements: " << gaussTotalPts << " pts\n";
    gsInfo << "Nested-union mass eval points:   " << massEvalPts
           << " pts  (ratio " << real_t(massEvalPts)/gaussTotalPts << "x vs Gauss)\n";
    gsInfo << "Nested-union stiff eval points:  " << stiffEvalPts
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
