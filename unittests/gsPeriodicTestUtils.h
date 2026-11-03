/** @file gsPeriodicTestUtils.h

    @brief Reference helpers for the periodic B-spline unit tests.

    Every expected value is built from a period description (degree,
    domain, breakpoints with multiplicity) and never from the object under
    test.  Include this header after gismo_unittest.h.  No helper calls a
    CHECK macro: each comparison returns a Verdict.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
**/

#pragma once

#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <utility>

namespace gismo {
namespace periodic_test {

/// Result of a comparison: ok, a magnitude (max error, residual, ...)
/// and a diagnostic that is non-empty whenever ok is false.
struct Verdict
{
    bool        ok;
    real_t      value;
    std::string why;
};

inline Verdict mkOk(real_t v = 0)
{
    Verdict r; r.ok = true; r.value = v; return r;
}

inline Verdict mkBad(real_t v, const std::string & why)
{
    Verdict r; r.ok = false; r.value = v; r.why = why; return r;
}

inline std::string num(real_t v)
{
    std::ostringstream os;
    os.precision(12);
    os << v;
    return os.str();
}

inline std::string vecStr(const std::vector<real_t> & v)
{
    std::ostringstream os;
    os.precision(12);
    os << "[";
    for (size_t i = 0; i != v.size(); ++i)
        os << (i ? ", " : "") << v[i];
    os << "]";
    return os.str();
}

// ---------------------------------------------------------------------
// Falsification switch
// ---------------------------------------------------------------------

/// True iff the environment variable GISMO_PERIODIC_FALSIFY names the
/// test that is currently running.
inline bool falsifying()
{
    const char * e = std::getenv("GISMO_PERIODIC_FALSIFY");
    if (NULL == e)
        return false;
    const UnitTest::TestDetails * td = UnitTest::CurrentTest::Details();
    return (NULL != td) && (0 == std::strcmp(e, td->testName));
}

/// Expected integer, shifted by one under falsification.
inline index_t fzN(index_t v) { return falsifying() ? v + 1 : v; }

/// Expected boolean, negated under falsification.
inline bool fzB(bool v) { return falsifying() ? !v : v; }

/// Entry (0,0) shifted by rel*max(1,max|e|) under falsification.
inline gsMatrix<real_t> perturbed(const gsMatrix<real_t> & e, real_t rel = 1e-6)
{
    if (!falsifying() || 0 == e.size())
        return e;
    gsMatrix<real_t> r = e;
    r(0, 0) += rel * (std::max)(real_t(1), e.cwiseAbs().maxCoeff());
    return r;
}

inline std::vector<real_t> perturbedVec(const std::vector<real_t> & e)
{
    if (!falsifying() || e.empty())
        return e;
    std::vector<real_t> r = e;
    real_t mx = 1;
    for (size_t i = 0; i != e.size(); ++i)
        mx = (std::max)(mx, math::abs(e[i]));
    r[0] += 1e-6 * mx;
    return r;
}

// ---------------------------------------------------------------------
// Period description and the stretched-knot generator
// ---------------------------------------------------------------------

struct PeriodDesc
{
    std::string         tag;
    short_t             p;
    real_t              a, b;
    std::vector<real_t> s; ///< breakpoints of [a,b) with multiplicity, s[0]==a

    index_t n()  const { return (index_t)s.size(); }
    index_t mu() const
    {
        index_t c = 0;
        for (size_t i = 0; i != s.size(); ++i)
            if (s[i] == a) ++c;
        return c;
    }
    index_t m()  const { return p + 1 - mu(); }
    real_t  L()  const { return b - a; }
    real_t  scale() const
    {
        return (std::max)((std::max)(math::abs(a), math::abs(b)), L());
    }
};

/// A periodic basis exists iff 1 <= mu <= p and n >= p+1.
inline bool constructible(const PeriodDesc & d)
{
    return d.mu() >= 1 && d.mu() <= d.p && d.n() >= d.p + 1;
}

inline index_t floorDiv(index_t q, index_t n)
{
    index_t r = q / n;
    if ((q % n != 0) && ((q < 0) != (n < 0)))
        --r;
    return r;
}

inline index_t floorMod(index_t q, index_t n)
{
    return q - n * floorDiv(q, n);
}

/// Stretched knots t_k = s[(k-m) mod n] + L*floor((k-m)/n), k=0..n+m+p.
/// The copies of b are stored as b itself.
inline std::vector<real_t> stretchedKnots(const PeriodDesc & d)
{
    const index_t n = d.n(), m = d.m(), p = d.p;
    std::vector<real_t> t(n + m + p + 1);
    for (index_t k = 0; k <= n + m + p; ++k)
    {
        const index_t q = k - m, f = floorDiv(q, n), r = floorMod(q, n);
        t[k] = (f == 1 && d.s[r] == d.a) ? d.b : d.s[r] + d.L() * f;
    }
    return t;
}

inline gsKnotVector<real_t> stretchedKV(const PeriodDesc & d)
{
    return gsKnotVector<real_t>(stretchedKnots(d), d.p);
}

/// Structural check of a raw stretched vector against a description.
inline Verdict selfCheck(const PeriodDesc & d, const std::vector<real_t> & t)
{
    const index_t n = d.n(), m = d.m(), p = d.p;
    const real_t tol = 1e-14 * d.scale();
    if ((index_t)t.size() != n + m + p + 1)
        return mkBad(0, "size " + num(t.size()) + " != " + num(n + m + p + 1));
    for (size_t k = 1; k < t.size(); ++k)
        if (t[k] < t[k - 1])
            return mkBad(0, "not non-decreasing at " + num(k));
    if (math::abs(t[p] - d.a) > tol)
        return mkBad(0, "t[p] != a");
    if (math::abs(t[n + m] - d.b) > tol)
        return mkBad(0, "t[n+m] != b");
    for (index_t k = 0; k <= m + p; ++k)
        if (math::abs(t[k + n] - (t[k] + d.L())) > tol)
            return mkBad(0, "t[" + num(k + n) + "] != t[" + num(k) + "]+L");
    if (!(t[m - 1] < d.a))
        return mkBad(0, "t[m-1] not below a");
    return mkOk();
}

// ---------------------------------------------------------------------
// Reference bases
// ---------------------------------------------------------------------

/// R1: non-periodic basis on the stretched knots (n+m functions).
inline gsBSplineBasis<real_t> makeR1(const PeriodDesc & d)
{
    return gsBSplineBasis<real_t>(stretchedKV(d));
}

/// R2: non-periodic basis on enough translated copies of the breakpoints
/// (built from s only); folded modulo n it spans the periodic space.
inline gsBSplineBasis<real_t> makeR2(const PeriodDesc & d)
{
    const index_t n = d.n();
    const index_t P = (d.p + n) / n + 1;
    std::vector<real_t> kn;
    for (index_t k = -P; k <= P; ++k)
        for (index_t i = 0; i < n; ++i)
            kn.push_back(d.s[i] + k * d.L());
    kn.push_back(d.a + (P + 1) * d.L());
    return gsBSplineBasis<real_t>(gsKnotVector<real_t>(kn, d.p));
}

/// Dense evaluation (nRows x u.cols()) of derivative order der, summing
/// the active functions of B into row (fold>0 ? index mod fold : index).
inline gsMatrix<real_t> denseEval(const gsBasis<real_t> & B,
                                  const gsMatrix<real_t> & u,
                                  index_t nRows, index_t fold, int der,
                                  bool & inRange)
{
    inRange = true;
    gsMatrix<real_t> vals;
    gsMatrix<index_t> act;
    if (0 == der)      B.eval_into(u, vals);
    else if (1 == der) B.deriv_into(u, vals);
    else if (2 == der) B.deriv2_into(u, vals);
    else
    {
        std::vector<gsMatrix<real_t> > all;
        B.evalAllDers_into(u, der, all);
        vals = all[der];
    }
    B.active_into(u, act);

    gsMatrix<real_t> out = gsMatrix<real_t>::Zero(nRows, u.cols());
    for (index_t c = 0; c != act.cols(); ++c)
        for (index_t r = 0; r != act.rows(); ++r)
        {
            const index_t idx = fold > 0 ? floorMod(act(r, c), fold) : act(r, c);
            if (idx < 0 || idx >= nRows)
            {
                inRange = false;
                continue;
            }
            out(idx, c) += vals(r, c);
        }
    return out;
}

inline std::vector<real_t> uniqueSorted(std::vector<real_t> v, real_t tol = 0)
{
    std::sort(v.begin(), v.end());
    std::vector<real_t> r;
    for (size_t i = 0; i != v.size(); ++i)
        if (r.empty() || v[i] - r.back() > tol)
            r.push_back(v[i]);
    return r;
}

/// Sample points of [a,b] (1 x npts): p+1 Gauss points per span of the
/// breakpoints of d (plus extraBreaks inside [a,b]), the breakpoints
/// themselves, and points approaching both domain ends.
inline gsMatrix<real_t> samplePoints(const PeriodDesc & d,
                                     const std::vector<real_t> & extraBreaks = std::vector<real_t>())
{
    std::vector<real_t> br(d.s.begin(), d.s.end());
    br.push_back(d.b);
    for (size_t i = 0; i != extraBreaks.size(); ++i)
        if (extraBreaks[i] >= d.a && extraBreaks[i] <= d.b)
            br.push_back(extraBreaks[i]);
    br = uniqueSorted(br, 1e-13 * d.scale());

    gsGaussRule<real_t> gr(d.p + 1);
    gsMatrix<real_t> nodes;
    gsVector<real_t> w;
    gr.mapToAll(br, nodes, w);

    std::vector<real_t> pts;
    for (index_t i = 0; i != nodes.cols(); ++i)
        pts.push_back(nodes(0, i));
    for (size_t i = 0; i != br.size(); ++i)
        pts.push_back(br[i]);
    const real_t h = (std::min)(br[1] - br[0], br[br.size() - 1] - br[br.size() - 2]);
    for (int k = 4; k <= 8; ++k)
    {
        const real_t e = h / real_t(1 << k);
        pts.push_back(d.a + e);
        pts.push_back(d.b - e);
    }
    gsMatrix<real_t> U(1, pts.size());
    for (size_t i = 0; i != pts.size(); ++i)
        U(0, i) = pts[i];
    return U;
}

inline gsMatrix<real_t> pointMat(real_t x)
{
    gsMatrix<real_t> u(1, 1);
    u(0, 0) = x;
    return u;
}

/// Folded R1 values (der-th derivative) of the description at U.
inline gsMatrix<real_t> r1Dense(const PeriodDesc & d, const gsMatrix<real_t> & U,
                                int der, bool & ok)
{
    const gsBSplineBasis<real_t> r1 = makeR1(d);
    return denseEval(r1, U, d.n(), d.n(), der, ok);
}

/// Left limit at b of the folded der-th derivatives: the polynomial of
/// the last span is interpolated at p+1 interior points and evaluated
/// exactly at b (for bases whose evaluation at b returns the right limit).
inline gsMatrix<real_t> denseEvalLeftAtB(const gsBasis<real_t> & B, const PeriodDesc & d,
                                         index_t nRows, index_t fold, int der)
{
    std::vector<real_t> br(d.s.begin(), d.s.end());
    br.push_back(d.b);
    br = uniqueSorted(br);
    const real_t h = br[br.size() - 1] - br[br.size() - 2];
    const index_t q = d.p + 1;
    gsMatrix<real_t> X(1, q);
    for (index_t j = 0; j != q; ++j)
        X(0, j) = d.b - h * real_t(j + 1) / real_t(q + 1);
    bool ok;
    const gsMatrix<real_t> V = denseEval(B, X, nRows, fold, der, ok);
    gsVector<real_t> wt(q);
    for (index_t j = 0; j != q; ++j)
    {
        real_t w = 1;
        for (index_t i = 0; i != q; ++i)
            if (i != j)
                w *= (d.b - X(0, i)) / (X(0, j) - X(0, i));
        wt[j] = w;
    }
    return V * wt;
}

// ---------------------------------------------------------------------
// Comparison helpers (expected input is perturbed when falsifying)
// ---------------------------------------------------------------------

inline Verdict valuesClose(const gsMatrix<real_t> & got, const gsMatrix<real_t> & expected)
{
    if (got.rows() != expected.rows() || got.cols() != expected.cols())
        return mkBad(0, "shape " + num(got.rows()) + "x" + num(got.cols()) +
                     " vs " + num(expected.rows()) + "x" + num(expected.cols()));
    const gsMatrix<real_t> e = perturbed(expected);
    const real_t err = (got - e).cwiseAbs().maxCoeff();
    if (!(err <= 1e-12))
        return mkBad(err, "max |got-expected| = " + num(err));
    return mkOk(err);
}

inline Verdict derivsClose(const gsMatrix<real_t> & got, const gsMatrix<real_t> & expected)
{
    if (got.rows() != expected.rows() || got.cols() != expected.cols())
        return mkBad(0, "shape " + num(got.rows()) + "x" + num(got.cols()) +
                     " vs " + num(expected.rows()) + "x" + num(expected.cols()));
    const gsMatrix<real_t> e = perturbed(expected);
    const real_t err = (got - e).cwiseAbs().maxCoeff();
    const real_t tol = 1e-11 * (std::max)(real_t(1), e.cwiseAbs().maxCoeff());
    if (!(err <= tol))
        return mkBad(err, "max |got-expected| = " + num(err) + " > " + num(tol));
    return mkOk(err);
}

/// Finite-difference sanity check, entrywise
/// |fd-deriv| <= max(absTol, relTol*|deriv|).
inline Verdict fdClose(const gsMatrix<real_t> & fd, const gsMatrix<real_t> & deriv,
                       real_t absTol, real_t relTol)
{
    if (fd.rows() != deriv.rows() || fd.cols() != deriv.cols())
        return mkBad(0, "shape mismatch");
    const gsMatrix<real_t> e = perturbed(deriv, 1e-3);
    real_t worst = 0;
    for (index_t c = 0; c != e.cols(); ++c)
        for (index_t r = 0; r != e.rows(); ++r)
        {
            const real_t err = math::abs(fd(r, c) - e(r, c));
            const real_t tol = (std::max)(absTol, relTol * math::abs(e(r, c)));
            if (!(err <= tol))
                return mkBad(err, "fd mismatch " + num(err) + " > " + num(tol) +
                             " at (" + num(r) + "," + num(c) + ")");
            worst = (std::max)(worst, err);
        }
    return mkOk(worst);
}

/// Central finite differences of order 1 (h=1e-5) or 2 (h=1e-4) of the
/// folded dense evaluation.
inline gsMatrix<real_t> fdMatrix(const gsBasis<real_t> & B, const gsMatrix<real_t> & U,
                                 index_t nRows, index_t fold, int order)
{
    const real_t h = (1 == order) ? 1e-5 : 1e-4;
    bool ok;
    gsMatrix<real_t> Up = U, Um = U;
    Up.array() += h;
    Um.array() -= h;
    const gsMatrix<real_t> vp = denseEval(B, Up, nRows, fold, 0, ok);
    const gsMatrix<real_t> vm = denseEval(B, Um, nRows, fold, 0, ok);
    if (1 == order)
        return (vp - vm) / (2 * h);
    const gsMatrix<real_t> v0 = denseEval(B, U, nRows, fold, 0, ok);
    return (vp - 2 * v0 + vm) / (h * h);
}

inline Verdict partitionOfUnity(const gsMatrix<real_t> & vals)
{
    const real_t target = falsifying() ? real_t(1) + 1e-6 : real_t(1);
    real_t worst = 0;
    for (index_t c = 0; c != vals.cols(); ++c)
        worst = (std::max)(worst, math::abs(vals.col(c).sum() - target));
    if (!(worst <= 1e-12))
        return mkBad(worst, "column sum deviates from 1 by " + num(worst));
    return mkOk(worst);
}

typedef std::vector<std::pair<gsMatrix<real_t>, gsMatrix<real_t> > > SeamData;

/// D[k] = (order-k derivatives at a, the same at b), k = 0..p-mu+1.
/// Continuous up to order p-mu, discontinuous at order p-mu+1.
inline Verdict seamContinuity(const SeamData & D0, short_t p, index_t mu)
{
    SeamData D = D0;
    if (falsifying() && !D.empty() && D[0].second.size() > 0)
        D[0].second(0, 0) += 1e-6 * (std::max)(real_t(1), D[0].second.cwiseAbs().maxCoeff());
    const index_t kmax = p - mu + 1;
    if ((index_t)D.size() < kmax + 1)
        return mkBad(0, "too few derivative orders");
    for (index_t k = 0; k <= kmax; ++k)
    {
        const real_t scale = (std::max)(real_t(1),
            (std::max)(D[k].first.cwiseAbs().maxCoeff(), D[k].second.cwiseAbs().maxCoeff()));
        const real_t jump = (D[k].first - D[k].second).cwiseAbs().maxCoeff();
        if (k < kmax)
        {
            if (!(jump <= 1e-10 * scale))
                return mkBad(jump, "order " + num(k) + " jump " + num(jump));
        }
        else if (!(jump > 1e-6 * scale))
            return mkBad(jump, "order " + num(k) + " has no jump (" + num(jump) + ")");
    }
    return mkOk();
}

/// Seam data of a folded basis, evaluating b directly (left limit for a
/// basis whose domain end is b) or by exact extrapolation from the last span.
inline SeamData seamData(const gsBasis<real_t> & B, const PeriodDesc & d,
                         index_t nRows, index_t fold, bool extrapolateB)
{
    SeamData D;
    bool ok;
    for (index_t k = 0; k <= d.p - d.mu() + 1; ++k)
    {
        const gsMatrix<real_t> A = denseEval(B, pointMat(d.a), nRows, fold, (int)k, ok);
        const gsMatrix<real_t> Bv = extrapolateB
            ? denseEvalLeftAtB(B, d, nRows, fold, (int)k)
            : denseEval(B, pointMat(d.b), nRows, fold, (int)k, ok);
        D.push_back(std::make_pair(A, Bv));
    }
    return D;
}

/// Equality of the spaces spanned by the rows of VA (nA x npts) and VB.
inline Verdict spanEqual(const gsMatrix<real_t> & VA, const gsMatrix<real_t> & VB0)
{
    gsMatrix<real_t> VB = VB0;
    if (falsifying() && VB.size() > 0)
        VB(0, 0) += 1e-6 * (std::max)(real_t(1), VB.cwiseAbs().maxCoeff());
    const gsMatrix<real_t> A = VA.transpose(), B = VB.transpose();
    if (A.rows() < A.cols() || B.rows() < B.cols())
        return mkBad(0, "fewer sample points than functions");
    gsMatrix<real_t>::JacobiSVD svdA = A.jacobiSvd(gsEigen::ComputeThinU | gsEigen::ComputeThinV);
    gsMatrix<real_t>::JacobiSVD svdB = B.jacobiSvd(gsEigen::ComputeThinU | gsEigen::ComputeThinV);
    const real_t sA = svdA.singularValues()[svdA.singularValues().size() - 1] / svdA.singularValues()[0];
    const real_t sB = svdB.singularValues()[svdB.singularValues().size() - 1] / svdB.singularValues()[0];
    if (!(sA > 1e-10))
        return mkBad(sA, "rank: first system, sigma_min/sigma_max = " + num(sA));
    if (!(sB > 1e-10))
        return mkBad(sB, "rank: second system, sigma_min/sigma_max = " + num(sB));
    const real_t r1 = (B - A * svdA.solve(B)).norm() / B.norm();
    const real_t r2 = (A - B * svdB.solve(A)).norm() / A.norm();
    if (!(r1 <= 1e-10))
        return mkBad(r1, "second space not inside first, relative residual " + num(r1));
    if (!(r2 <= 1e-10))
        return mkBad(r2, "first space not inside second, relative residual " + num(r2));
    return mkOk((std::max)(r1, r2));
}

/// Every active index lies in [0,n) and none repeats within a column.
inline Verdict activeUniqueInRange(const gsMatrix<index_t> & act, index_t n)
{
    for (index_t c = 0; c != act.cols(); ++c)
        for (index_t r = 0; r != act.rows(); ++r)
        {
            if (act(r, c) < 0 || act(r, c) >= n)
                return mkBad(0, "index " + num(act(r, c)) + " out of [0," + num(n) + ")");
            for (index_t r2 = r + 1; r2 != act.rows(); ++r2)
                if (act(r2, c) == act(r, c))
                    return mkBad(0, "duplicate index " + num(act(r, c)) + " in column " + num(c));
        }
    return mkOk();
}

/// Vector comparison of knots, tolerance 1e-14*scale; size is compared exactly.
inline Verdict knotsClose(const std::vector<real_t> & got, const std::vector<real_t> & expected,
                          real_t scale)
{
    const std::vector<real_t> e = perturbedVec(expected);
    if (got.size() != e.size())
        return mkBad(0, "knot count " + num(got.size()) + " vs " + num(e.size()) +
                     ": got " + vecStr(got) + " expected " + vecStr(expected));
    real_t err = 0;
    for (size_t i = 0; i != got.size(); ++i)
        err = (std::max)(err, math::abs(got[i] - e[i]));
    if (!(err <= 1e-14 * scale))
        return mkBad(err, "knots differ by " + num(err) + ": got " + vecStr(got) +
                     " expected " + vecStr(expected));
    return mkOk(err);
}

/// Scalar comparison, tolerance 1e-14*scale on the (perturbed) expectation.
inline Verdict scalarClose(real_t got, real_t expected, real_t scale)
{
    const real_t e = falsifying() ? expected + 1e-6 * (std::max)(real_t(1), math::abs(expected))
                                  : expected;
    if (!(math::abs(got - e) <= 1e-14 * scale))
        return mkBad(math::abs(got - e), "got " + num(got) + " expected " + num(expected));
    return mkOk();
}

/// Basis under test is periodic with n functions (and, when ghost, m
/// crossing functions).  Not falsified.
inline bool guardPeriodic(const gsBSplineBasis<real_t> & b, const PeriodDesc & d, bool ghost)
{
    return b.isPeriodic() && b.size() == d.n() &&
           (!ghost || b.numCrossingFunctions() == d.m());
}

inline std::string guardMsg(const char * what, const gsBSplineBasis<real_t> & b, const PeriodDesc & d)
{
    std::ostringstream os;
    os << what << ": periodic=" << (b.isPeriodic() ? 1 : 0) << " size=" << b.size()
       << " expected n=" << d.n();
    return os.str();
}

inline std::vector<real_t> knotsOf(const gsBSplineBasis<real_t> & b)
{
    const gsKnotVector<real_t> & kv = b.knots();
    return std::vector<real_t>(kv.begin(), kv.end());
}

/// Tested basis (fold 0, d.n() rows) against folded R1 at samplePoints(d),
/// derivative order der (values: valuesClose, der>0: derivsClose).
inline Verdict basisVsR1(const gsBasis<real_t> & b, const PeriodDesc & d, int der = 0)
{
    const gsMatrix<real_t> U = samplePoints(d);
    bool ok1, ok2;
    const gsMatrix<real_t> got = denseEval(b, U, d.n(), 0, der, ok1);
    const gsMatrix<real_t> ref = r1Dense(d, U, der, ok2);
    if (!ok1)
        return mkBad(0, "active index outside [0,n)");
    return 0 == der ? valuesClose(got, ref) : derivsClose(got, ref);
}

// ---------------------------------------------------------------------
// Snapshot of a basis for refusal tests
// ---------------------------------------------------------------------

struct BasisSnapshot
{
    std::vector<real_t> knots;
    index_t             size;
    short_t             degree;
    bool                periodic;
};

inline BasisSnapshot snapshot(const gsBSplineBasis<real_t> & b)
{
    BasisSnapshot s;
    s.knots    = knotsOf(b);
    s.size     = b.size();
    s.degree   = b.degree();
    s.periodic = b.isPeriodic();
    return s;
}

/// Exact comparison of knots, size, degree and periodicity.
inline Verdict unchanged(const BasisSnapshot & before0, const gsBSplineBasis<real_t> & after)
{
    BasisSnapshot before = before0;
    if (falsifying() && !before.knots.empty())
        before.knots[0] += 1e-6 * (std::max)(real_t(1), math::abs(before.knots[0]));
    const std::vector<real_t> k = knotsOf(after);
    if (k != before.knots)
        return mkBad(0, "knots changed: " + vecStr(before.knots) + " -> " + vecStr(k));
    if (after.size() != before.size)
        return mkBad(0, "size changed: " + num(before.size) + " -> " + num(after.size()));
    if (after.degree() != before.degree)
        return mkBad(0, "degree changed");
    if (after.isPeriodic() != before.periodic)
        return mkBad(0, "periodic flag changed");
    return mkOk();
}

// ---------------------------------------------------------------------
// Operations on period descriptions (expected results of mutations)
// ---------------------------------------------------------------------

inline index_t countEq(const std::vector<real_t> & s, real_t x)
{
    index_t c = 0;
    for (size_t i = 0; i != s.size(); ++i)
        if (s[i] == x) ++c;
    return c;
}

inline PeriodDesc withTag(PeriodDesc d, const std::string & suffix)
{
    d.tag += suffix;
    return d;
}

inline void pushMult(std::vector<real_t> & s, real_t x, index_t mult)
{
    for (index_t i = 0; i < mult; ++i)
        s.push_back(x);
}

inline PeriodDesc opInsert(PeriodDesc d, real_t x)
{
    d.s.push_back(x);
    std::sort(d.s.begin(), d.s.end());
    return withTag(d, "+ins");
}

inline PeriodDesc opRemove(PeriodDesc d, real_t x)
{
    std::vector<real_t>::iterator it = std::find(d.s.begin(), d.s.end(), x);
    if (it != d.s.end())
        d.s.erase(it);
    return withTag(d, "-rem");
}

inline PeriodDesc opUniformRefine(PeriodDesc d, index_t numKnots, index_t mul)
{
    std::vector<real_t> ub(d.s.begin(), d.s.end());
    ub.push_back(d.b);
    ub = uniqueSorted(ub);
    for (size_t i = 0; i + 1 < ub.size(); ++i)
        for (index_t j = 1; j <= numKnots; ++j)
            pushMult(d.s, ub[i] + (ub[i + 1] - ub[i]) * real_t(j) / real_t(numKnots + 1), mul);
    std::sort(d.s.begin(), d.s.end());
    return withTag(d, "+ref");
}

inline PeriodDesc opDegreeElevate(PeriodDesc d, index_t i)
{
    std::vector<real_t> u = uniqueSorted(d.s), s;
    for (size_t j = 0; j != u.size(); ++j)
        pushMult(s, u[j], countEq(d.s, u[j]) + i);
    d.s = s;
    d.p = d.p + (short_t)i;
    return withTag(d, "+de");
}

inline PeriodDesc opDegreeIncrease(PeriodDesc d, index_t i)
{
    d.p = d.p + (short_t)i;
    return withTag(d, "+di");
}

inline PeriodDesc opDegreeReduce(PeriodDesc d, index_t i)
{
    std::vector<real_t> u = uniqueSorted(d.s), s;
    for (size_t j = 0; j != u.size(); ++j)
    {
        const index_t mult = countEq(d.s, u[j]);
        pushMult(s, u[j], u[j] == d.a ? (std::max)(index_t(1), mult - i)
                                      : (std::max)(index_t(0), mult - i));
    }
    d.s = s;
    d.p = d.p - (short_t)i;
    return withTag(d, "-dr");
}

inline PeriodDesc opDegreeDecrease(PeriodDesc d, index_t i)
{
    d.p = d.p - (short_t)i;
    return withTag(d, "-dd");
}

inline PeriodDesc opElevateContinuity(PeriodDesc d, index_t i)
{
    std::vector<real_t> u = uniqueSorted(d.s), s;
    for (size_t j = 0; j != u.size(); ++j)
    {
        const index_t mult = countEq(d.s, u[j]);
        pushMult(s, u[j], u[j] == d.a ? (std::max)(index_t(1), mult - i)
                                      : (std::max)(index_t(0), mult - i));
    }
    d.s = s;
    return withTag(d, "+ec");
}

inline PeriodDesc opReduceContinuity(PeriodDesc d, index_t i)
{
    std::vector<real_t> u = uniqueSorted(d.s), s;
    for (size_t j = 0; j != u.size(); ++j)
    {
        const index_t mult = countEq(d.s, u[j]);
        pushMult(s, u[j], u[j] == d.a ? (std::min)(index_t(d.p), mult + i)
                                      : (std::min)(index_t(d.p) + 1, mult + i));
    }
    d.s = s;
    return withTag(d, "-rc");
}

/// Walks the distinct breakpoints in (a,b): removes knotRemove (all
/// copies), keeps knotSkip, and repeats.
inline PeriodDesc opCoarsen(PeriodDesc d, index_t knotRemove, index_t knotSkip = 1)
{
    std::vector<real_t> u = uniqueSorted(d.s), inner, gone;
    for (size_t j = 0; j != u.size(); ++j)
        if (u[j] != d.a)
            inner.push_back(u[j]);
    size_t i = 0;
    while (i < inner.size())
    {
        for (index_t c = 0; c != knotRemove && i < inner.size(); ++c, ++i)
            gone.push_back(inner[i]);
        i += knotSkip;
    }
    std::vector<real_t> s;
    for (size_t j = 0; j != d.s.size(); ++j)
        if (std::find(gone.begin(), gone.end(), d.s[j]) == gone.end())
            s.push_back(d.s[j]);
    d.s = s;
    return withTag(d, "-co");
}

inline PeriodDesc opReverse(PeriodDesc d)
{
    std::vector<real_t> s, rest;
    pushMult(s, d.a, d.mu());
    for (size_t j = 0; j != d.s.size(); ++j)
        if (d.s[j] != d.a)
            rest.push_back(d.a + d.b - d.s[j]);
    std::sort(rest.begin(), rest.end());
    s.insert(s.end(), rest.begin(), rest.end());
    d.s = s;
    return withTag(d, "rev");
}

// ---------------------------------------------------------------------
// Cases
// ---------------------------------------------------------------------

inline PeriodDesc mkDesc(const std::string & tag, short_t p, real_t a, real_t b,
                         const std::vector<real_t> & s)
{
    PeriodDesc d;
    d.tag = tag; d.p = p; d.a = a; d.b = b; d.s = s;
    return d;
}

inline std::vector<real_t> vec(const real_t * v, size_t k)
{
    return std::vector<real_t>(v, v + k);
}

inline std::string tagPN(const char * prefix, int p, int N)
{
    std::ostringstream os;
    os << prefix << "P" << p << "N" << N;
    return os.str();
}

inline PeriodDesc caseC2u(short_t p, index_t N)
{
    std::vector<real_t> s;
    for (index_t i = 0; i != N; ++i)
        s.push_back(real_t(i) / real_t(N));
    return mkDesc(tagPN("C2u", p, (int)N), p, 0, 1, s);
}

inline PeriodDesc caseC2a()
{
    std::vector<real_t> s;
    for (int i = 0; i != 4; ++i)
        s.push_back(-1 + real_t(i) * 7 / 8);
    return mkDesc("C2a", 2, -1, 2.5, s);
}

inline PeriodDesc caseC3n()
{
    const real_t v[] = {0, .5, 1.5, 2, 3, 3.5};
    return mkDesc("C3n", 2, 0, 4, vec(v, 6));
}

inline PeriodDesc caseC3y()
{
    const real_t v[] = {0, .5, 1.5, 2.5, 3, 3.5};
    return mkDesc("C3y", 2, 0, 4, vec(v, 6));
}

inline PeriodDesc caseC4i()
{
    const real_t v[] = {0, 1, 1, 2, 3};
    return mkDesc("C4i", 2, 0, 4, vec(v, 5));
}

inline PeriodDesc caseC4o()
{
    const real_t v[] = {0, 1, 2, 3, 3, 4, 5};
    return mkDesc("C4o", 2, 0, 6, vec(v, 7));
}

inline PeriodDesc caseE2()
{
    const real_t v[] = {0, 0, .25, .25, .5, .5, .75, .75};
    return mkDesc("E2", 3, 0, 1, vec(v, 8));
}

inline PeriodDesc caseC5a()
{
    const real_t v[] = {0, .5};
    return mkDesc("C5a", 3, 0, 1, vec(v, 2));
}

inline PeriodDesc caseC5b()
{
    const real_t v[] = {0, 1, 2};
    return mkDesc("C5b", 3, 0, 3, vec(v, 3));
}

/// Clamped conversion: p copies of 0 followed by i/N, i=1..N-1.
inline PeriodDesc caseC1k(short_t p, index_t N)
{
    std::vector<real_t> s(p, real_t(0));
    for (index_t i = 1; i != N; ++i)
        s.push_back(real_t(i) / real_t(N));
    return mkDesc(tagPN("C1k", p, (int)N), p, 0, 1, s);
}

inline PeriodDesc caseC1kL4(short_t p)
{
    std::vector<real_t> s(p, real_t(0));
    s.push_back(1); s.push_back(2); s.push_back(3);
    std::ostringstream os;
    os << "C1kL4P" << p;
    return mkDesc(os.str(), p, 0, 4, s);
}

/// From-scratch constructor: i/N with multiplicity mult.
inline PeriodDesc caseC1c(short_t p, index_t N, index_t mult = 1)
{
    std::vector<real_t> s;
    for (index_t i = 0; i != N; ++i)
        pushMult(s, real_t(i) / real_t(N), mult);
    std::string tag = tagPN("C1c", p, (int)N);
    if (mult != 1)
        tag += "M" + num(mult);
    return mkDesc(tag, p, 0, 1, s);
}

inline std::vector<PeriodDesc> constructibleCases()
{
    std::vector<PeriodDesc> c;
    const int pn[7][2] = {{1,4},{2,4},{3,4},{1,8},{2,8},{3,8},{4,8}};
    for (int i = 0; i != 7; ++i)
        c.push_back(caseC2u((short_t)pn[i][0], pn[i][1]));
    c.push_back(caseC2a());
    c.push_back(caseC3n());
    c.push_back(caseC3y());
    c.push_back(caseC4i());
    c.push_back(caseC4o());
    c.push_back(caseE2());
    return c;
}

inline std::vector<PeriodDesc> c2uCases()
{
    std::vector<PeriodDesc> all = constructibleCases(), c;
    for (size_t i = 0; i != 7; ++i)
        c.push_back(all[i]);
    return c;
}

inline std::vector<PeriodDesc> c1cCases()
{
    std::vector<PeriodDesc> c;
    const int pn[7][2] = {{1,4},{2,4},{3,4},{1,8},{2,8},{3,8},{4,8}};
    for (int i = 0; i != 7; ++i)
        c.push_back(caseC1c((short_t)pn[i][0], pn[i][1]));
    return c;
}

inline std::vector<PeriodDesc> c1kCases()
{
    std::vector<PeriodDesc> c;
    for (short_t p = 1; p <= 4; ++p)
    {
        c.push_back(caseC1k(p, 4));
        c.push_back(caseC1k(p, 8));
    }
    for (short_t p = 1; p <= 4; ++p)
        c.push_back(caseC1kL4(p));
    return c;
}

/// Non-periodic clamped quadratic basis on [0,.25,.5,.75,1].
inline gsKnotVector<real_t> c6Knots()
{
    const real_t v[] = {0, 0, 0, .25, .5, .75, 1, 1, 1};
    return gsKnotVector<real_t>(vec(v, 9), 2);
}

inline gsBSplineBasis<real_t> makeC6()
{
    return gsBSplineBasis<real_t>(c6Knots());
}

// ---------------------------------------------------------------------
// Test drivers
// ---------------------------------------------------------------------

/// Runs f on every case; exceptions become failure text.  Returns the
/// concatenated failures ("" if none).
inline std::string runCases(const std::vector<PeriodDesc> & cases,
                            const std::function<std::string(const PeriodDesc &)> & f)
{
    std::string fails;
    for (size_t i = 0; i != cases.size(); ++i)
    {
        std::string w;
        try { w = f(cases[i]); }
        catch (std::exception & e) { w = std::string("exception: ") + e.what(); }
        if (!w.empty())
            fails += cases[i].tag + ": " + w + "; ";
    }
    return fails;
}

/// Diagnostic of a verdict: empty iff ok.
inline std::string whyOf(const Verdict & v)
{
    return v.ok ? std::string() : (v.why.empty() ? std::string("failed") : v.why);
}

/// Prints what the non-periodic twin of a periodic refinement test does:
/// R1 (same stretched knots, flag off) with the same operation applied.
inline void printTwin(const char * name, const PeriodDesc & d,
                      std::function<void(gsBSplineBasis<real_t> &)> op)
{
    try
    {
        gsBSplineBasis<real_t> b = makeR1(d);
        bool threw = false;
        std::string what;
        try { op(b); }
        catch (std::exception & e) { threw = true; what = e.what(); }
        gsInfo << "info_twin of " << name << " on " << d.tag << ": threw=" << (threw ? 1 : 0);
        if (threw)
            gsInfo << " (" << what << ")";
        gsInfo << " size=" << b.size() << " knots=" << vecStr(knotsOf(b)) << "\n";
    }
    catch (std::exception & e)
    {
        gsInfo << "info_twin of " << name << ": exception: " << e.what() << "\n";
    }
}

} // namespace periodic_test
} // namespace gismo
