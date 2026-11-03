/** @file gsPeriodicBSpline_test.cpp

    @brief Periodic B-spline curves, tensor bases and geometries, multipatch
    topology of periodic patches and XML round trips.

    Expected values come from period descriptions and from non-periodic
    reference bases (see gsPeriodicTestUtils.h), never from the object under
    test.  Curve and tensor-geometry references sum the index-aligned
    reference functions with the coefficients of the n reduced rows.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
**/

#include "gismo_unittest.h"
#include "gsPeriodicTestUtils.h"

#include <cstring>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <new>

using namespace gismo::periodic_test;

namespace {

typedef gsBSpline<real_t>               Crv;
typedef gsBSplineBasis<real_t>          Basis;
typedef gsTensorBSplineBasis<2, real_t> TB;
typedef gsTensorBSpline<2, real_t>      TG;
typedef std::function<void(Crv &)>      COp;
typedef std::function<void(TG &)>       GOp;

const real_t PI = 3.14159265358979323846;

/// Number of functions of the clamped, non-periodic direction of the tensor cases.
const index_t NNP = 5;

// ---------------------------------------------------------------------
// Coefficients and references (a) curve, (b) tensor geometry
// ---------------------------------------------------------------------

/// n reduced rows of a generic closed planar curve.
gsMatrix<real_t> crvCoefs(const PeriodDesc & d)
{
    const index_t n = d.n();
    gsMatrix<real_t> c(n, 2);
    for (index_t j = 0; j != n; ++j)
    {
        const real_t th = 2 * PI * real_t(j) / real_t(n);
        c(j, 0) = math::cos(th) + real_t(0.1) * j * j;
        c(j, 1) = math::sin(th) - real_t(0.05) * j;
    }
    return c;
}

/// Rows 0..n-1 of cn followed by m ghost rows, row n+j equal to row j.
gsMatrix<real_t> expandRows(const gsMatrix<real_t> & cn, index_t n, index_t m)
{
    gsMatrix<real_t> e(n + m, cn.cols());
    for (index_t i = 0; i != n + m; ++i)
        e.row(i) = cn.row(i % n);
    return e;
}

real_t maxAbs(const gsMatrix<real_t> & c)
{
    return (0 == c.size()) ? real_t(1) : c.cwiseAbs().maxCoeff();
}

/// Entrywise comparison, tolerance rel*scale; the expectation is perturbed
/// under falsification.
std::string geoClose(const gsMatrix<real_t> & got, const gsMatrix<real_t> & expected,
                     real_t scale, real_t rel = 1e-12)
{
    if (got.rows() != expected.rows() || got.cols() != expected.cols())
        return "shape " + num(got.rows()) + "x" + num(got.cols()) + " vs " +
               num(expected.rows()) + "x" + num(expected.cols());
    if (0 == got.size())
        return "empty";
    const gsMatrix<real_t> e = perturbed(expected);
    const real_t err = (got - e).cwiseAbs().maxCoeff();
    if (!(err <= rel * scale))
        return "max |got-expected| = " + num(err) + " > " + num(rel * scale);
    return "";
}

/// Reference (a): x(u) = sum_i c_{i mod n} B^{R1}_i(u), geoDim x npts, for the
/// der-th derivative.
gsMatrix<real_t> curveRef(const PeriodDesc & d, const gsMatrix<real_t> & cn,
                          const gsMatrix<real_t> & U, int der = 0)
{
    bool ok;
    const gsMatrix<real_t> B = r1Dense(d, U, der, ok);
    return cn.transpose() * B;
}

gsMatrix<real_t> crvDer(const Crv & c, const gsMatrix<real_t> & U, int k)
{
    gsMatrix<real_t> r;
    if (0 == k)      c.eval_into(U, r);
    else if (1 == k) c.deriv_into(U, r);
    else             c.deriv2_into(U, r);
    return r;
}

std::string curveVsRef(const Crv & c, const PeriodDesc & d, const gsMatrix<real_t> & cn,
                       const gsMatrix<real_t> & U, int der = 0, real_t rel = 1e-12)
{
    const gsMatrix<real_t> got = crvDer(c, U, der);
    const gsMatrix<real_t> ref = curveRef(d, cn, U, der);
    return 0 == der ? geoClose(got, ref, maxAbs(cn), rel) : whyOf(derivsClose(got, ref));
}

std::vector<real_t> knotsVec(const Crv & c)
{
    return std::vector<real_t>(c.knots().begin(), c.knots().end());
}

/// Periodic, n basis functions and n coefficient rows.
std::string crvGuard(const Crv & c, const PeriodDesc & d)
{
    if (c.isPeriodic() && c.basis().size() == d.n() && c.coefs().rows() == d.n())
        return "";
    std::ostringstream os;
    os << "guard: periodic=" << (c.isPeriodic() ? 1 : 0) << " size=" << c.basis().size()
       << " rows=" << c.coefs().rows() << " expected n=" << d.n();
    return os.str();
}

#define GUARD_CRV(c, d) \
    { const std::string guard_ = crvGuard(c, d); \
      CHECK_EQUAL(std::string(""), guard_); \
      if (!guard_.empty()) return; }

/// Main route: periodic basis on the stretched knots, then the n reduced rows.
Crv makeCurve(const PeriodDesc & d, const gsMatrix<real_t> & cn)
{
    Basis b(stretchedKV(d));
    b.setPeriodic();
    return Crv(b, cn);
}

bool throwsRT(const std::function<void()> & f)
{
    try { f(); }
    catch (std::runtime_error &) { return true; }
    catch (...) { return false; }
    return false;
}

struct CrvSnap
{
    bool                per;
    index_t             size;
    short_t             deg;
    std::vector<real_t> knots;
    gsMatrix<real_t>    coefs;
};

CrvSnap crvSnap(const Crv & c)
{
    CrvSnap s;
    s.per = c.isPeriodic();
    s.size = c.basis().size();
    s.deg = c.degree();
    s.knots = knotsVec(c);
    s.coefs = c.coefs();
    return s;
}

/// Knots, size, degree, flag and coefficients equal those of the snapshot.
std::string crvUnchanged(const CrvSnap & s0, const Crv & c)
{
    CrvSnap s = s0;
    s.coefs = perturbed(s0.coefs);
    if (knotsVec(c) != s.knots)
        return "knots changed";
    if (c.basis().size() != s.size)
        return "size changed: " + num(s.size) + " -> " + num(c.basis().size());
    if (c.degree() != s.deg)
        return "degree changed";
    if (c.isPeriodic() != s.per)
        return "periodic flag changed";
    if (c.coefs().rows() != s.coefs.rows() || c.coefs().cols() != s.coefs.cols())
        return "coefficient shape changed";
    if (c.coefs() != s.coefs)
        return "coefficients changed";
    return "";
}

/// op must throw std::runtime_error and leave c as it was.
void refuseCrv(Crv & c, const COp & op)
{
    const CrvSnap s0 = crvSnap(c);
    CHECK_THROW(op(c), std::runtime_error);
    CHECK_EQUAL(std::string(""), crvUnchanged(s0, c));
}

/// Periodic curve of d, op, then the periodic curve of dp with the values of
/// the curve before the operation.
void crvRow(const PeriodDesc & d, const COp & op, const PeriodDesc & dp)
{
    const gsMatrix<real_t> cn = crvCoefs(d);
    CHECK(constructible(dp) && selfCheck(dp, stretchedKnots(dp)).ok);
    Crv c = makeCurve(d, cn);
    GUARD_CRV(c, d)
    try { op(c); }
    catch (std::exception & e)
    {
        CHECK_EQUAL(std::string(""), std::string("operation threw: ") + e.what());
        return;
    }
    GUARD_CRV(c, dp)
    CHECK_EQUAL(fzN(dp.p), (index_t)c.degree());
    CHECK_EQUAL(fzN(dp.m()), (index_t)c.basis().numCrossingFunctions());
    CHECK_EQUAL(std::string(""), whyOf(knotsClose(knotsVec(c), stretchedKnots(dp), dp.scale())));
    const gsMatrix<real_t> U = samplePoints(dp);
    gsMatrix<real_t> got;
    c.eval_into(U, got);
    CHECK_EQUAL(std::string(""), geoClose(got, curveRef(d, cn, U), maxAbs(cn)));
}

/// Like crvRow for operations that approximate: the new curve interpolates the
/// old one at the anchors of its new basis.
void crvApprox(const PeriodDesc & d, const COp & op, const PeriodDesc & dp)
{
    const gsMatrix<real_t> cn = crvCoefs(d);
    CHECK(constructible(dp) && selfCheck(dp, stretchedKnots(dp)).ok);
    Crv c = makeCurve(d, cn);
    GUARD_CRV(c, d)
    try { op(c); }
    catch (std::exception & e)
    {
        CHECK_EQUAL(std::string(""), std::string("operation threw: ") + e.what());
        return;
    }
    GUARD_CRV(c, dp)
    CHECK_EQUAL(fzN(dp.p), (index_t)c.degree());
    CHECK_EQUAL(fzN(dp.m()), (index_t)c.basis().numCrossingFunctions());
    const gsMatrix<real_t> A = c.basis().anchors();
    CHECK_EQUAL(fzN(dp.n()), (index_t)A.cols());
    gsMatrix<real_t> got;
    c.eval_into(A, got);
    CHECK_EQUAL(std::string(""), geoClose(got, curveRef(d, cn, A), maxAbs(cn), 1e-10));
}

COp cInsert(real_t x)    { return [=](Crv & c) { c.insertKnot(x); }; }
COp cInsertKnots(real_t x)
{
    return [=](Crv & c) { std::vector<real_t> v(1, x); c.insertKnots(v.begin(), v.end()); };
}
COp cUniformRefine()     { return [](Crv & c) { c.uniformRefine(); }; }
COp cElevate()           { return [](Crv & c) { c.degreeElevate((short_t)1); }; }
COp cIncrease()          { return [](Crv & c) { c.degreeIncrease((short_t)1); }; }
COp cDecrease()          { return [](Crv & c) { c.degreeDecrease((short_t)1); }; }
COp cReduce()            { return [](Crv & c) { c.degreeReduce((short_t)1); }; }
COp cReverse()           { return [](Crv & c) { c.reverse(); }; }

/// Prints what the non-periodic twin of a curve test does: the same stretched
/// knots with n+m rows, flag off, the same operation, error against the
/// unchanged (or, with reflect, the reflected) reference curve.  The line
/// starts with B32 for refinement twins and with "twin" otherwise.
void twinCrv(const char * name, const PeriodDesc & d, const COp & op,
             bool reflect = false, bool refinement = true)
{
    try
    {
        const gsMatrix<real_t> cn = crvCoefs(d);
        Crv c(stretchedKV(d), expandRows(cn, d.n(), d.m()));
        bool threw = false;
        std::string what;
        try { op(c); }
        catch (std::exception & e) { threw = true; what = e.what(); }
        std::ostringstream os;
        os << (refinement ? "B32: " : "twin: ") << name << " on " << d.tag
           << ": threw=" << (threw ? 1 : 0);
        if (threw)
            os << " (" << what << ")";
        os << " size=" << c.basis().size() << " rows=" << c.coefs().rows();
        if (!threw)
        {
            const gsMatrix<real_t> U = samplePoints(d);
            gsMatrix<real_t> Ur = U;
            if (reflect)
                for (index_t i = 0; i != U.cols(); ++i)
                    Ur(0, i) = d.a + d.b - U(0, i);
            gsMatrix<real_t> got;
            c.eval_into(U, got);
            os << " maxValueError=" << (got - curveRef(d, cn, Ur)).cwiseAbs().maxCoeff();
        }
        gsInfo << os.str() << "\n";
    }
    catch (std::exception & e)
    {
        gsInfo << (refinement ? "B32: " : "twin: ") << name << ": exception: " << e.what() << "\n";
    }
}

/// Left limit at b of the k-th derivative (k <= 2) of a curve: the last span
/// is interpolated at p+1 interior points and evaluated at b.
gsMatrix<real_t> crvLeftAtB(const Crv & c, const PeriodDesc & d, int k)
{
    std::vector<real_t> br(d.s.begin(), d.s.end());
    br.push_back(d.b);
    br = uniqueSorted(br);
    const real_t h = br[br.size() - 1] - br[br.size() - 2];
    const index_t q = d.p + 1;
    gsMatrix<real_t> X(1, q);
    for (index_t j = 0; j != q; ++j)
        X(0, j) = d.b - h * real_t(j + 1) / real_t(q + 1);
    const gsMatrix<real_t> V = crvDer(c, X, k);
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

/// C^{p-mu} at the seam, jump in the derivative of order p-mu+1 (orders <= 2).
std::string crvSeam(const Crv & c, const PeriodDesc & d)
{
    if (d.p - d.mu() + 1 > 2)
        return "orders above 2 are not supported by this helper";
    SeamData D;
    for (int k = 0; k <= d.p - d.mu() + 1; ++k)
        D.push_back(std::make_pair(crvDer(c, pointMat(d.a), k), crvLeftAtB(c, d, k)));
    return whyOf(seamContinuity(D, d.p, d.mu()));
}

// ---------------------------------------------------------------------
// Tensor references: (b) geometry, basis; clamped direction fixed
// ---------------------------------------------------------------------

/// Clamped quadratic knot vector of the non-periodic direction, NNP functions.
gsKnotVector<real_t> npKnots()
{
    const real_t v[] = {0, 0, 0, 0.25, 0.5, 1, 1, 1};
    return gsKnotVector<real_t>(vec(v, 8), 2);
}

/// Direction 0 runs fastest: periodic direction 0 gives iper + inp*nper,
/// periodic direction 1 gives inp + iper*nnp.
index_t gidx(int dir, index_t iper, index_t inp, index_t nper, index_t nnp)
{
    return 0 == dir ? iper + inp * nper : inp + iper * nnp;
}

/// n*NNP reduced rows of a generic planar surface.
gsMatrix<real_t> tgCoefs(const PeriodDesc & d, int dir)
{
    const index_t n = d.n();
    gsMatrix<real_t> c(n * NNP, 2);
    for (index_t ip = 0; ip != n; ++ip)
        for (index_t iq = 0; iq != NNP; ++iq)
        {
            const real_t th = 2 * PI * real_t(ip) / real_t(n);
            c(gidx(dir, ip, iq, n, NNP), 0) = math::cos(th) + real_t(0.1) * ip * ip + real_t(0.3) * iq;
            c(gidx(dir, ip, iq, n, NNP), 1) = math::sin(th) - real_t(0.05) * ip + real_t(0.2) * iq * iq;
        }
    return c;
}

/// (n+m) slices along the periodic direction, ghost slice n+j equal to slice j.
gsMatrix<real_t> expandTG(const gsMatrix<real_t> & cn, const PeriodDesc & d, int dir)
{
    const index_t n = d.n(), m = d.m();
    gsMatrix<real_t> e((n + m) * NNP, cn.cols());
    for (index_t ip = 0; ip != n + m; ++ip)
        for (index_t iq = 0; iq != NNP; ++iq)
            e.row(gidx(dir, ip, iq, n + m, NNP)) = cn.row(gidx(dir, ip % n, iq, n, NNP));
    return e;
}

/// 2 x npts tensor grid: the sample points of d along the periodic direction,
/// seven values along the clamped direction on [0,1].
gsMatrix<real_t> tgPoints(const PeriodDesc & d, int dir)
{
    const gsMatrix<real_t> P = samplePoints(d);
    const real_t v[] = {0, 0.1, 0.25, 0.375, 0.5, 0.8, 1};
    gsMatrix<real_t> U(2, P.cols() * 7);
    index_t k = 0;
    for (int j = 0; j != 7; ++j)
        for (index_t i = 0; i != P.cols(); ++i, ++k)
        {
            U(dir, k) = P(0, i);
            U(1 - dir, k) = v[j];
        }
    return U;
}

/// Reference (b): the product of the folded R1 values of the periodic
/// direction and the clamped basis of the other direction, summed with the
/// n*nnp reduced coefficient rows.
gsMatrix<real_t> tgRef(const PeriodDesc & d, int dir, const gsMatrix<real_t> & cn,
                       const gsMatrix<real_t> & U)
{
    bool ok;
    const gsMatrix<real_t> up = U.row(dir), un = U.row(1 - dir);
    const gsMatrix<real_t> Bp = r1Dense(d, up, 0, ok);
    const Basis np(npKnots());
    const gsMatrix<real_t> Bn = denseEval(np, un, NNP, 0, 0, ok);
    gsMatrix<real_t> out = gsMatrix<real_t>::Zero(cn.cols(), U.cols());
    for (index_t pt = 0; pt != U.cols(); ++pt)
        for (index_t ip = 0; ip != d.n(); ++ip)
            for (index_t iq = 0; iq != NNP; ++iq)
                out.col(pt) += cn.row(gidx(dir, ip, iq, d.n(), NNP)).transpose() * (Bp(ip, pt) * Bn(iq, pt));
    return out;
}

/// All n*nnp basis function values: product of the folded R1 values and the
/// values of the clamped basis on knots np.
gsMatrix<real_t> tbRef(const PeriodDesc & d, int dir, const gsMatrix<real_t> & U,
                       const gsKnotVector<real_t> & np)
{
    bool ok;
    const gsMatrix<real_t> up = U.row(dir), un = U.row(1 - dir);
    const gsMatrix<real_t> Bp = r1Dense(d, up, 0, ok);
    const Basis npb(np);
    const index_t nnp = npb.size();
    const gsMatrix<real_t> Bn = denseEval(npb, un, nnp, 0, 0, ok);
    gsMatrix<real_t> out = gsMatrix<real_t>::Zero(d.n() * nnp, U.cols());
    for (index_t pt = 0; pt != U.cols(); ++pt)
        for (index_t ip = 0; ip != d.n(); ++ip)
            for (index_t iq = 0; iq != nnp; ++iq)
                out(gidx(dir, ip, iq, d.n(), nnp), pt) = Bp(ip, pt) * Bn(iq, pt);
    return out;
}

/// Periodic tensor basis, periodic in direction dir.
TB makeTB(const PeriodDesc & d, int dir)
{
    TB tb = (0 == dir) ? TB(stretchedKV(d), npKnots()) : TB(npKnots(), stretchedKV(d));
    tb.setPeriodic(dir);
    return tb;
}

/// Default route: periodic tensor basis, then the n*NNP reduced rows.
TG makeTG(const PeriodDesc & d, int dir, const gsMatrix<real_t> & cn)
{
    return TG(makeTB(d, dir), cn);
}

std::string tbGuardN(const TB & b, int dir, index_t nper, index_t nnp)
{
    if (b.isPeriodic() && b.periodicDirection() == dir && b.size() == nper * nnp)
        return "";
    std::ostringstream os;
    os << "guard: periodic=" << (b.isPeriodic() ? 1 : 0) << " periodicDirection=" << b.periodicDirection()
       << " size=" << b.size() << " expected " << nper * nnp << " (direction " << dir << ")";
    return os.str();
}

std::string tgGuardN(const TG & g, int dir, index_t nper, index_t nnp)
{
    std::string w = tbGuardN(g.basis(), dir, nper, nnp);
    if (w.empty() && g.coefs().rows() != nper * nnp)
        w = "guard: coefficient rows " + num(g.coefs().rows()) + " != size " + num(nper * nnp);
    return w;
}

#define GUARD_TB(b, dir, nper, nnp) \
    { const std::string guard_ = tbGuardN(b, dir, nper, nnp); \
      CHECK_EQUAL(std::string(""), guard_); \
      if (!guard_.empty()) return; }

#define GUARD_TG(g, dir, nper, nnp) \
    { const std::string guard_ = tgGuardN(g, dir, nper, nnp); \
      CHECK_EQUAL(std::string(""), guard_); \
      if (!guard_.empty()) return; }

std::vector<real_t> tgKnots(const TG & g, int dir)
{
    const gsKnotVector<real_t> & kv = g.basis().knots(dir);
    return std::vector<real_t>(kv.begin(), kv.end());
}

/// Periodic tensor geometry of d, op, then the periodic geometry of dp
/// (nnp functions along the clamped direction) with the values of the
/// geometry before the operation.
void tgRow(const PeriodDesc & d, int dir, const GOp & op, const PeriodDesc & dp, index_t nnp)
{
    const gsMatrix<real_t> cn = tgCoefs(d, dir);
    CHECK(constructible(dp) && selfCheck(dp, stretchedKnots(dp)).ok);
    TG g = makeTG(d, dir, cn);
    GUARD_TG(g, dir, d.n(), NNP)
    try { op(g); }
    catch (std::exception & e)
    {
        CHECK_EQUAL(std::string(""), std::string("operation threw: ") + e.what());
        return;
    }
    GUARD_TG(g, dir, dp.n(), nnp)
    CHECK_EQUAL(fzN(dp.p), (index_t)g.basis().degree(dir));
    CHECK_EQUAL(std::string(""), whyOf(knotsClose(tgKnots(g, dir), stretchedKnots(dp), dp.scale())));
    const gsMatrix<real_t> U = tgPoints(dp, dir);
    gsMatrix<real_t> got;
    g.eval_into(U, got);
    CHECK_EQUAL(std::string(""), geoClose(got, tgRef(d, dir, cn, U), maxAbs(cn)));
}

/// Non-periodic twin of a tensor-geometry test: the same stretched knots
/// with (n+m)*NNP rows, flag off, the same operation.  With reflect the
/// reference is evaluated at the periodic coordinate reflected about a+b.
void twinTG(const char * name, const PeriodDesc & d, int dir, const GOp & op, bool reflect = false)
{
    try
    {
        const gsMatrix<real_t> cn = tgCoefs(d, dir);
        const gsMatrix<real_t> ex = expandTG(cn, d, dir);
        TG g = (0 == dir) ? TG(stretchedKV(d), npKnots(), ex) : TG(npKnots(), stretchedKV(d), ex);
        bool threw = false;
        std::string what;
        try { op(g); }
        catch (std::exception & e) { threw = true; what = e.what(); }
        std::ostringstream os;
        os << "B32: " << name << " on " << d.tag << " direction " << dir << ": threw=" << (threw ? 1 : 0);
        if (threw)
            os << " (" << what << ")";
        os << " size=" << g.basis().size() << " rows=" << g.coefs().rows();
        if (!threw)
        {
            const gsMatrix<real_t> U = tgPoints(d, dir);
            gsMatrix<real_t> Ur = U;
            if (reflect)
                for (index_t i = 0; i != U.cols(); ++i)
                    Ur(dir, i) = d.a + d.b - U(dir, i);
            gsMatrix<real_t> got;
            g.eval_into(U, got);
            os << " maxValueError=" << (got - tgRef(d, dir, cn, Ur)).cwiseAbs().maxCoeff();
        }
        gsInfo << os.str() << "\n";
    }
    catch (std::exception & e)
    {
        gsInfo << "B32: " << name << ": exception: " << e.what() << "\n";
    }
}

std::string xmlPath(const char * name)
{
    return gsFileManager::getTempPath() + "/gsPeriodicBSpline_test_" + name + ".xml";
}

} // namespace

SUITE(gsPeriodicBSpline_test)
{

// =====================================================================
// G1: curve storage and constructors
// =====================================================================

void kvCtorExpanded(const PeriodDesc & d)
{
    const gsMatrix<real_t> cn = crvCoefs(d);
    Crv c(stretchedKV(d), expandRows(cn, d.n(), d.m()), true);
    CHECK_EQUAL(fzB(true), c.isPeriodic());
    if (!c.isPeriodic() || c.basis().size() != d.n())
    {
        CHECK_EQUAL(std::string(""), crvGuard(c, d));
        return;
    }
    CHECK_EQUAL(fzN(d.n()), (index_t)c.basis().size());
    CHECK_EQUAL(fzN(d.n()), (index_t)c.coefs().rows());
    CHECK_EQUAL(std::string(""), curveVsRef(c, d, cn, samplePoints(d)));
}

TEST(space_B16_crvKVCtorExpandedRowsCuni)  { kvCtorExpanded(caseC2u(2, 4)); }
TEST(space_B16_crvKVCtorExpandedRowsCnonu) { kvCtorExpanded(caseC3n()); }

void kvCtorReduced(const PeriodDesc & d)
{
    const gsMatrix<real_t> cn = crvCoefs(d);
    Crv c(stretchedKV(d), cn, true);
    CHECK_EQUAL(fzB(true), c.isPeriodic());
    if (!c.isPeriodic() || c.basis().size() != d.n())
    {
        CHECK_EQUAL(std::string(""), crvGuard(c, d));
        return;
    }
    CHECK_EQUAL(fzN(d.n()), (index_t)c.basis().size());
    CHECK_EQUAL(fzN(d.n()), (index_t)c.coefs().rows());
    CHECK_EQUAL(std::string(""), curveVsRef(c, d, cn, samplePoints(d)));
}

TEST(space_B16_crvKVCtorReducedRowsCuni)  { kvCtorReduced(caseC2u(2, 4)); }
TEST(space_B16_crvKVCtorReducedRowsCnonu) { kvCtorReduced(caseC3n()); }

void basisCtorRows(const PeriodDesc & d)
{
    const gsMatrix<real_t> cn = crvCoefs(d);
    Basis b(stretchedKV(d));
    b.setPeriodic();
    if (!guardPeriodic(b, d, false))
    {
        CHECK_EQUAL(std::string(""), guardMsg("pre-op guard", b, d));
        return;
    }
    Crv c(b, cn);
    CHECK_EQUAL(fzB(true), c.isPeriodic());
    CHECK_EQUAL(fzN(d.n()), (index_t)c.coefs().rows());
    if (c.isPeriodic() && c.basis().size() == d.n())
        CHECK_EQUAL(std::string(""), curveVsRef(c, d, cn, samplePoints(d)));
}

TEST(space_B16_crvBasisCtorRowsCuni)  { basisCtorRows(caseC2u(2, 4)); }
TEST(space_B16_crvBasisCtorRowsCnonu) { basisCtorRows(caseC3n()); }

TEST(space_NEW_crvSetPeriodicStorageCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = crvCoefs(d);
    Crv c(stretchedKV(d), expandRows(cn, d.n(), d.m()));
    CHECK_EQUAL(fzB(false), c.isPeriodic());
    c.setPeriodic();
    GUARD_CRV(c, d)
    CHECK_EQUAL(std::string(""), curveVsRef(c, d, cn, samplePoints(d)));
}

TEST(space_NEW_crvExpandedRowsMismatchRefusedCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    gsMatrix<real_t> bad = expandRows(crvCoefs(d), d.n(), d.m());
    bad(d.n(), 0) += 1e-3;
    CHECK_EQUAL(fzB(true), throwsRT([&] { (void)Crv(stretchedKV(d), bad, true); }));
}

TEST(space_NEW_crvCtorPeriodicKVFlagFalseCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = crvCoefs(d);
    Basis bPer(stretchedKV(d));
    bPer.setPeriodic();
    if (!guardPeriodic(bPer, d, false))
    {
        CHECK_EQUAL(std::string(""), guardMsg("pre-op guard", bPer, d));
        return;
    }
    gsKnotVector<real_t> kv = bPer.knots();
    Crv c(kv, cn, false);
    GUARD_CRV(c, d)
    CHECK_EQUAL(std::string(""), curveVsRef(c, d, cn, samplePoints(d)));
}

TEST(space_NEW_crvUniformCtorCctor)
{
    const std::string fails = runCases(c1cCases(), [](const PeriodDesc & d) -> std::string
    {
        const gsMatrix<real_t> cn = crvCoefs(d);
        Crv c(real_t(0), real_t(1), (unsigned)(d.n() - 1), d.p, cn, 1, true);
        if (!c.isPeriodic() || c.basis().size() != d.n())
            return "guard: periodic=" + num(c.isPeriodic()) + " size=" + num(c.basis().size());
        std::string w;
        if (fzN(d.n()) != (index_t)c.coefs().rows())
            w += "rows " + num(c.coefs().rows()) + " != " + num(d.n()) + "; ";
        return w + curveVsRef(c, d, cn, samplePoints(d));
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_crvUniformCtorMultRefusedCctor)
{
    gsMatrix<real_t> cf = gsMatrix<real_t>::Random(12, 2);
    CHECK_EQUAL(fzB(true), throwsRT([&] { (void)Crv(real_t(0), real_t(1), 3u, 2, cf, 3u, true); }));
}

TEST(space_B1_crvClampedConversionCclamp)
{
    const PeriodDesc dk = caseC1k(2, 4);
    const index_t n = dk.n();
    gsMatrix<real_t> closed(n + 1, 2);
    closed.topRows(n) = crvCoefs(dk);
    closed.row(n) = closed.row(0);
    const gsKnotVector<real_t> kv(dk.a, dk.b, (unsigned)3, dk.p + 1, 1, dk.p);
    const Crv c0(kv, closed);
    Crv c = c0;
    c.setPeriodic();
    GUARD_CRV(c, dk)
    const gsMatrix<real_t> U = samplePoints(dk);
    gsMatrix<real_t> got, ref;
    c.eval_into(U, got);
    c0.eval_into(U, ref);
    CHECK_EQUAL(std::string(""), geoClose(got, ref, maxAbs(closed)));
}

void seamCurve(const PeriodDesc & d)
{
    Crv c = makeCurve(d, crvCoefs(d));
    GUARD_CRV(c, d)
    CHECK_EQUAL(std::string(""), crvSeam(c, d));
}

TEST(space_NEW_crvSeamContinuityCuni)  { seamCurve(caseC2u(2, 4)); }
TEST(space_NEW_crvSeamContinuityCnonu) { seamCurve(caseC3n()); }

TEST(space_NEW_crvUnwrapCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = crvCoefs(d);
    Crv c = makeCurve(d, cn);
    GUARD_CRV(c, d)
    c.setPeriodic(false);
    CHECK_EQUAL(fzB(false), c.isPeriodic());
    CHECK_EQUAL(fzN(d.n() + d.m()), (index_t)c.basis().size());
    CHECK_EQUAL(fzN(d.n() + d.m()), (index_t)c.coefs().rows());
    CHECK_EQUAL(std::string(""), curveVsRef(c, d, cn, samplePoints(d)));
}

TEST(space_B9_crvDerivReducedRowsCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = crvCoefs(d);
    Crv c(stretchedKV(d), expandRows(cn, d.n(), d.m()), true);
    GUARD_CRV(c, d)
    const gsMatrix<real_t> U = samplePoints(d);
    CHECK_EQUAL(std::string(""), curveVsRef(c, d, cn, U, 1));
    CHECK_EQUAL(std::string(""), curveVsRef(c, d, cn, U, 2));
}

TEST(space_B19_crvCoordCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = crvCoefs(d);
    Crv c(stretchedKV(d), expandRows(cn, d.n(), d.m()), true);
    GUARD_CRV(c, d)
    try
    {
        const gsGeometry<real_t>::uPtr c0 = c.coord(0);
        CHECK_EQUAL(fzN(d.n()), (index_t)c0->coefs().rows());
        const gsMatrix<real_t> U = samplePoints(d);
        gsMatrix<real_t> got;
        c0->eval_into(U, got);
        const gsMatrix<real_t> ref = curveRef(d, cn, U);
        CHECK_EQUAL(std::string(""), geoClose(got, ref.topRows(1), maxAbs(cn)));
    }
    catch (std::exception & e)
    {
        CHECK_EQUAL(std::string(""), std::string("coord threw: ") + e.what());
    }
}

TEST(space_B35_crvControlNetCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    Crv c = makeCurve(d, crvCoefs(d));
    GUARD_CRV(c, d)
    try
    {
        gsMesh<real_t> m;
        c.controlNet(m);
        CHECK_EQUAL(fzN(d.n()), (index_t)m.numVertices());
    }
    catch (std::exception & e)
    {
        CHECK_EQUAL(std::string(""), std::string("controlNet threw: ") + e.what());
    }
}

// =====================================================================
// G2: curve refinement and degree
// =====================================================================

TEST(space_B6_crvUniformRefineCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    crvRow(d, cUniformRefine(), opUniformRefine(d, 1, 1));
}
TEST(info_twin_crvUniformRefineCuni) { twinCrv("crvUniformRefine", caseC2u(2, 4), cUniformRefine()); }

TEST(space_B6_crvUniformRefineCnonu)
{
    const PeriodDesc d = caseC3n();
    crvRow(d, cUniformRefine(), opUniformRefine(d, 1, 1));
}
TEST(info_twin_crvUniformRefineCnonu) { twinCrv("crvUniformRefine", caseC3n(), cUniformRefine()); }

TEST(space_B6_crvUniformRefineCdblin)
{
    const PeriodDesc d = caseC4i();
    crvRow(d, cUniformRefine(), opUniformRefine(d, 1, 1));
}
TEST(info_twin_crvUniformRefineCdblin) { twinCrv("crvUniformRefine", caseC4i(), cUniformRefine()); }

TEST(space_B6_crvUniformRefineCdblout)
{
    const PeriodDesc d = caseC4o();
    crvRow(d, cUniformRefine(), opUniformRefine(d, 1, 1));
}
TEST(info_twin_crvUniformRefineCdblout) { twinCrv("crvUniformRefine", caseC4o(), cUniformRefine()); }

TEST(space_B6_crvTransferShapeCuni)
{
    const PeriodDesc d = caseC2u(2, 4), dp = opUniformRefine(d, 1, 1);
    Basis b(stretchedKV(d));
    b.setPeriodic();
    if (!guardPeriodic(b, d, false))
    {
        CHECK_EQUAL(std::string(""), guardMsg("pre-op guard", b, d));
        return;
    }
    gsSparseMatrix<real_t, RowMajor> T;
    b.uniformRefine_withTransfer(T, 1, 1);
    CHECK_EQUAL(fzN(dp.n()), (index_t)T.rows());
    CHECK_EQUAL(fzN(d.n()), (index_t)T.cols());
    if (T.rows() != dp.n() || T.cols() != d.n())
        return;
    const gsMatrix<real_t> cn = crvCoefs(d);
    const gsMatrix<real_t> cp = T * cn;
    const gsMatrix<real_t> U = samplePoints(dp);
    CHECK_EQUAL(std::string(""), geoClose(curveRef(dp, cp, U), curveRef(d, cn, U), maxAbs(cn)));
}

TEST(info_twin_crvTransferShapeCuni)
{
    try
    {
        const PeriodDesc d = caseC2u(2, 4);
        Basis b(stretchedKV(d));
        std::ostringstream os;
        os << "B32: crvTransferShape on " << d.tag << ": ";
        try
        {
            gsSparseMatrix<real_t, RowMajor> T;
            b.uniformRefine_withTransfer(T, 1, 1);
            os << "threw=0 transfer " << T.rows() << "x" << T.cols() << " size=" << b.size();
        }
        catch (std::exception & e)
        {
            os << "threw=1 (" << e.what() << ")";
        }
        gsInfo << os.str() << "\n";
    }
    catch (std::exception & e)
    {
        gsInfo << "B32: crvTransferShape: exception: " << e.what() << "\n";
    }
}

TEST(space_B17_crvInsertKnotGhostWindowCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    crvRow(d, cInsert(0.125), opInsert(d, 0.125));
}
TEST(info_twin_crvInsertKnotGhostWindowCuni) { twinCrv("crvInsertKnotGhostWindow", caseC2u(2, 4), cInsert(0.125)); }

TEST(space_B17_crvInsertKnotGhostWindowCnonu)
{
    const PeriodDesc d = caseC3n();
    crvRow(d, cInsert(0.25), opInsert(d, 0.25));
}
TEST(info_twin_crvInsertKnotGhostWindowCnonu) { twinCrv("crvInsertKnotGhostWindow", caseC3n(), cInsert(0.25)); }

TEST(space_B17_crvInsertKnotWrapCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    crvRow(d, cInsert(1.125), opInsert(d, 0.125));
}
TEST(info_twin_crvInsertKnotWrapCuni) { twinCrv("crvInsertKnotWrap", caseC2u(2, 4), cInsert(1.125)); }

TEST(space_NEW_crvInsertKnotInteriorCdblout)
{
    const PeriodDesc d = caseC4o();
    crvRow(d, cInsert(3.5), opInsert(d, 3.5));
}
TEST(info_twin_crvInsertKnotInteriorCdblout) { twinCrv("crvInsertKnotInterior", caseC4o(), cInsert(3.5)); }

TEST(space_NEW_crvInsertKnotSeamRefusedCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    Crv c = makeCurve(d, crvCoefs(d));
    GUARD_CRV(c, d)
    refuseCrv(c, cInsert(d.a));
    refuseCrv(c, cInsert(d.b));
}
TEST(info_twin_crvInsertKnotSeamRefusedCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    twinCrv("crvInsertKnotSeamA", d, cInsert(d.a), false, false);
    twinCrv("crvInsertKnotSeamB", d, cInsert(d.b), false, false);
}

TEST(space_NEW_crvInsertKnotFarRefusedCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    Crv c = makeCurve(d, crvCoefs(d));
    GUARD_CRV(c, d)
    refuseCrv(c, cInsert(d.a - real_t(1.5) * d.L()));
}
TEST(info_twin_crvInsertKnotFarRefusedCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    twinCrv("crvInsertKnotFar", d, cInsert(d.a - real_t(1.5) * d.L()), false, false);
}

TEST(space_NEW_crvInsertKnotsSeamARefusedCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    Crv c = makeCurve(d, crvCoefs(d));
    GUARD_CRV(c, d)
    refuseCrv(c, cInsertKnots(d.a));
}
TEST(info_twin_crvInsertKnotsSeamARefusedCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    twinCrv("crvInsertKnotsSeamA", d, cInsertKnots(d.a), false, false);
}

TEST(space_B18_crvInsertKnotsSeamBRefusedCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    Crv c = makeCurve(d, crvCoefs(d));
    GUARD_CRV(c, d)
    refuseCrv(c, cInsertKnots(d.b));
}
TEST(info_twin_crvInsertKnotsSeamBRefusedCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    twinCrv("crvInsertKnotsSeamB", d, cInsertKnots(d.b), false, false);
}

TEST(space_B26_crvInsertKnotsGhostWindowCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    crvRow(d, cInsertKnots(0.125), opInsert(d, 0.125));
}
TEST(info_twin_crvInsertKnotsGhostWindowCuni) { twinCrv("crvInsertKnotsGhostWindow", caseC2u(2, 4), cInsertKnots(0.125)); }

TEST(space_B26_crvInsertKnotsGhostWindowCnonu)
{
    const PeriodDesc d = caseC3n();
    crvRow(d, cInsertKnots(0.25), opInsert(d, 0.25));
}
TEST(info_twin_crvInsertKnotsGhostWindowCnonu) { twinCrv("crvInsertKnotsGhostWindow", caseC3n(), cInsertKnots(0.25)); }

TEST(space_B27_crvDegreeElevateCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    crvRow(d, cElevate(), opDegreeElevate(d, 1));
}
TEST(info_twin_crvDegreeElevateCuni) { twinCrv("crvDegreeElevate", caseC2u(2, 4), cElevate()); }

TEST(space_B27_crvDegreeElevateCnonu)
{
    const PeriodDesc d = caseC3n();
    crvRow(d, cElevate(), opDegreeElevate(d, 1));
}
TEST(info_twin_crvDegreeElevateCnonu) { twinCrv("crvDegreeElevate", caseC3n(), cElevate()); }

void reverseCurve(const PeriodDesc & d)
{
    const gsMatrix<real_t> cn = crvCoefs(d);
    Crv c = makeCurve(d, cn);
    GUARD_CRV(c, d)
    c.reverse();
    GUARD_CRV(c, d)
    const gsMatrix<real_t> pr = c.parameterRange();
    CHECK_EQUAL(std::string(""), whyOf(scalarClose(pr(0, 0), d.a, d.scale())));
    CHECK_EQUAL(std::string(""), whyOf(scalarClose(pr(0, 1), d.b, d.scale())));
    const gsMatrix<real_t> U = samplePoints(d);
    gsMatrix<real_t> Ur = U;
    for (index_t i = 0; i != U.cols(); ++i)
        Ur(0, i) = d.a + d.b - U(0, i);
    gsMatrix<real_t> got;
    c.eval_into(U, got);
    CHECK_EQUAL(std::string(""), geoClose(got, curveRef(d, cn, Ur), maxAbs(cn)));
}

TEST(space_B27_crvReverseCuni)  { reverseCurve(caseC2u(2, 4)); }
TEST(info_twin_crvReverseCuni)  { twinCrv("crvReverse", caseC2u(2, 4), cReverse(), true); }
TEST(space_B27_crvReverseCnonu) { reverseCurve(caseC3n()); }
TEST(info_twin_crvReverseCnonu) { twinCrv("crvReverse", caseC3n(), cReverse(), true); }

TEST(space_NEW_crvDegreeIncreaseCuni)
{
    const PeriodDesc d = caseC2u(2, 8);
    crvApprox(d, cIncrease(), opDegreeIncrease(d, 1));
}
TEST(info_twin_crvDegreeIncreaseCuni) { twinCrv("crvDegreeIncrease", caseC2u(2, 8), cIncrease(), false, false); }

TEST(space_NEW_crvDegreeIncreaseRefusedCuni)
{
    const PeriodDesc d = caseC2u(3, 4);
    Crv c = makeCurve(d, crvCoefs(d));
    GUARD_CRV(c, d)
    refuseCrv(c, cIncrease());
}
TEST(info_twin_crvDegreeIncreaseRefusedCuni) { twinCrv("crvDegreeIncreaseRefused", caseC2u(3, 4), cIncrease(), false, false); }

TEST(space_NEW_crvDegreeDecreaseCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    crvApprox(d, cDecrease(), opDegreeDecrease(d, 1));
}
TEST(info_twin_crvDegreeDecreaseCuni) { twinCrv("crvDegreeDecrease", caseC2u(2, 4), cDecrease(), false, false); }

TEST(space_NEW_crvDegreeReduceRefusedCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    Crv c = makeCurve(d, crvCoefs(d));
    GUARD_CRV(c, d)
    refuseCrv(c, cReduce());
}
TEST(info_twin_crvDegreeReduceRefusedCuni) { twinCrv("crvDegreeReduceRefused", caseC2u(2, 4), cReduce(), false, false); }

TEST(space_NEW_crvDegreeReduceCuni)
{
    const PeriodDesc d = caseE2();
    crvApprox(d, cReduce(), opDegreeReduce(d, 1));
}
TEST(info_twin_crvDegreeReduceCuni) { twinCrv("crvDegreeReduce", caseE2(), cReduce(), false, false); }

TEST(space_B10_crvInterpolateAtAnchorsCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = crvCoefs(d);
    Basis b(stretchedKV(d));
    b.setPeriodic();
    if (!guardPeriodic(b, d, false))
    {
        CHECK_EQUAL(std::string(""), guardMsg("pre-op guard", b, d));
        return;
    }
    const gsMatrix<real_t> A = b.anchors();
    CHECK_EQUAL(fzN(d.n()), (index_t)A.cols());
    if (A.cols() != d.n())
        return;
    const gsMatrix<real_t> vals = curveRef(d, cn, A);
    try
    {
        const gsGeometry<real_t>::uPtr g = b.interpolateAtAnchors(vals);
        const Crv * pc = dynamic_cast<const Crv *>(g.get());
        CHECK_EQUAL(fzB(true), NULL != pc);
        if (NULL == pc)
            return;
        GUARD_CRV(*pc, d)
        gsMatrix<real_t> got;
        pc->eval_into(A, got);
        CHECK_EQUAL(std::string(""), geoClose(got, vals, maxAbs(cn), 1e-10));
    }
    catch (std::exception & e)
    {
        CHECK_EQUAL(std::string(""), std::string("interpolateAtAnchors threw: ") + e.what());
    }
}

// =====================================================================
// G3: curve policies
// =====================================================================

TEST(space_B36_crvCornerFunctionsRefusedCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = crvCoefs(d);
    Crv c = makeCurve(d, cn);
    GUARD_CRV(c, d)
    const CrvSnap s0 = crvSnap(c);
    const gsMatrix<real_t> v = cn.row(0).transpose();
    gsVector<index_t, 1> curr;
    CHECK_THROW(c.isPatchCorner(v), std::runtime_error);
    CHECK_THROW(c.findCorner(v, curr), std::runtime_error);
    CHECK_THROW(c.setOriginCorner(v), std::runtime_error);
    CHECK_THROW(c.setFurthestCorner(v), std::runtime_error);
    CHECK_EQUAL(std::string(""), crvUnchanged(s0, c));
}

TEST(space_NEW_crvNoPeriodicMeaningRefusedCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = crvCoefs(d);
    Crv c = makeCurve(d, cn);
    GUARD_CRV(c, d)
    Crv other = c;
    const CrvSnap s0 = crvSnap(c);
    Crv l, r;
    CHECK_THROW(c.segmentFromTo(0.25, 0.75), std::runtime_error);
    CHECK_THROW(c.splitAt(0.5, l, r), std::runtime_error);
    CHECK_THROW(c.toBezier(), std::runtime_error);
    CHECK_THROW(c.merge(&other), std::runtime_error);
    CHECK_THROW(c.intersect(other), std::runtime_error);
    CHECK_EQUAL(std::string(""), crvUnchanged(s0, c));
}

// =====================================================================
// G4: tensor basis, periodic in direction 0 (D0) or 1 (D1)
// =====================================================================

void tbActive(int dir)
{
    const PeriodDesc d = caseC2u(2, 4);
    TB tb = makeTB(d, dir);
    GUARD_TB(tb, dir, d.n(), NNP)
    gsMatrix<index_t> act;
    tb.active_into(tgPoints(d, dir), act);
    CHECK_EQUAL(fzN(9), (index_t)act.rows());
    CHECK_EQUAL(std::string(""), whyOf(activeUniqueInRange(act, tb.size())));
}
TEST(space_B11_tbActiveD0) { tbActive(0); }
TEST(space_B11_tbActiveD1) { tbActive(1); }

void tbEval(int dir)
{
    const PeriodDesc d = caseC2u(2, 4);
    TB tb = makeTB(d, dir);
    GUARD_TB(tb, dir, d.n(), NNP)
    const gsMatrix<real_t> U = tgPoints(d, dir);
    bool ok;
    const gsMatrix<real_t> got = denseEval(tb, U, d.n() * NNP, 0, 0, ok);
    CHECK_EQUAL(fzB(true), ok);
    CHECK_EQUAL(std::string(""), whyOf(valuesClose(got, tbRef(d, dir, U, npKnots()))));
}
TEST(space_B11_tbEvalD0) { tbEval(0); }
TEST(space_B11_tbEvalD1) { tbEval(1); }

void tbAllBoundary(int dir)
{
    const PeriodDesc d = caseC2u(2, 4);
    TB tb = makeTB(d, dir);
    GUARD_TB(tb, dir, d.n(), NNP)
    std::vector<index_t> expected;
    for (index_t ip = 0; ip != d.n(); ++ip)
    {
        expected.push_back(gidx(dir, ip, 0, d.n(), NNP));
        expected.push_back(gidx(dir, ip, NNP - 1, d.n(), NNP));
    }
    std::sort(expected.begin(), expected.end());
    if (falsifying())
        expected[0] += 1;
    const gsMatrix<index_t> ab = tb.allBoundary();
    std::vector<index_t> got(ab.data(), ab.data() + ab.size());
    std::sort(got.begin(), got.end());
    CHECK(got == expected);
}
TEST(space_B12_tbAllBoundaryD0) { tbAllBoundary(0); }
TEST(space_B12_tbAllBoundaryD1) { tbAllBoundary(1); }

void tbBoundaryOffsetRefused(int dir)
{
    const PeriodDesc d = caseC2u(2, 4);
    TB tb = makeTB(d, dir);
    GUARD_TB(tb, dir, d.n(), NNP)
    const boxSide s0 = (0 == dir) ? boxSide(boundary::west) : boxSide(boundary::south);
    const boxSide s1 = (0 == dir) ? boxSide(boundary::east) : boxSide(boundary::north);
    CHECK_THROW(tb.boundaryOffset(s0, 0), std::runtime_error);
    CHECK_EQUAL(fzB(true), throwsRT([&] { (void)tb.boundaryOffset(s1, 0); }));
}
TEST(space_B12_tbBoundaryOffsetRefusedD0) { tbBoundaryOffsetRefused(0); }
TEST(space_B12_tbBoundaryOffsetRefusedD1) { tbBoundaryOffsetRefused(1); }

TEST(space_NEW_tbFunctionAtCornerRefusedD0)
{
    const PeriodDesc d = caseC2u(2, 4);
    TB tb = makeTB(d, 0);
    GUARD_TB(tb, 0, d.n(), NNP)
    CHECK_THROW(tb.functionAtCorner(boxCorner::getFirst(2)), std::runtime_error);
    CHECK_EQUAL(fzB(true), throwsRT([&] { (void)tb.functionAtCorner(boxCorner::getLast(2)); }));
}

// ---------------------------------------------------------------------
// ghost_ block: internals of the tensor storage (diagnostic)
// ---------------------------------------------------------------------

void tbPerCoefs(int dir)
{
    const PeriodDesc d = caseC2u(2, 4);
    TB tb = makeTB(d, dir);
    GUARD_TB(tb, dir, d.n(), NNP)
    const index_t n = d.n(), m = d.m();
    gsMatrix<real_t> in((n + m) * NNP, 2);
    for (index_t ip = 0; ip != n + m; ++ip)
        for (index_t iq = 0; iq != NNP; ++iq)
        {
            in(gidx(dir, ip, iq, n + m, NNP), 0) = 100 * ip + iq + real_t(0.5);
            in(gidx(dir, ip, iq, n + m, NNP), 1) = 7 * ip - iq;
        }
    gsMatrix<real_t> expected = in;
    for (index_t j = 0; j != m; ++j)
        for (index_t iq = 0; iq != NNP; ++iq)
            expected.row(gidx(dir, n + j, iq, n + m, NNP)) = in.row(gidx(dir, j, iq, n + m, NNP));
    const gsMatrix<real_t> got = tb.perCoefs(in, (short_t)dir);
    CHECK_EQUAL(std::string(""), geoClose(got, expected, maxAbs(in)));
}
TEST(ghost_B13_tbPerCoefsD0) { tbPerCoefs(0); }
TEST(ghost_B13_tbPerCoefsD1) { tbPerCoefs(1); }

// ---------------------------------------------------------------------

TEST(space_B14_tbDefaultCtorNotPeriodic)
{
    alignas(TB) unsigned char buf[sizeof(TB)];
    std::memset(buf, 0, sizeof(buf));
    TB * b = new (buf) TB();
    CHECK_EQUAL(fzB(false), b->isPeriodic());
    CHECK_EQUAL(fzN(-1), (index_t)b->periodicDirection());
    b->~TB();
}

void tbSwapDirections(int dir)
{
    const PeriodDesc d = caseC2u(2, 4);
    TB tb = makeTB(d, dir);
    GUARD_TB(tb, dir, d.n(), NNP)
    tb.swapDirections(0, 1);
    CHECK_EQUAL(fzB(true), tb.isPeriodic());
    CHECK_EQUAL(fzN(1 - dir), (index_t)tb.periodicDirection());
    const gsMatrix<real_t> U = tgPoints(d, 1 - dir);
    gsMatrix<index_t> act;
    tb.active_into(U, act);
    CHECK_EQUAL(std::string(""), whyOf(activeUniqueInRange(act, tb.size())));
    bool ok;
    const gsMatrix<real_t> got = denseEval(tb, U, d.n() * NNP, 0, 0, ok);
    CHECK_EQUAL(std::string(""), whyOf(valuesClose(got, tbRef(d, 1 - dir, U, npKnots()))));
}
TEST(space_B15_tbSwapDirectionsD0) { tbSwapDirections(0); }
TEST(space_B15_tbSwapDirectionsD1) { tbSwapDirections(1); }

TEST(space_NEW_tbCopyKeepsPeriodicD0)
{
    const PeriodDesc d = caseC2u(2, 4);
    TB tb = makeTB(d, 0);
    GUARD_TB(tb, 0, d.n(), NNP)
    const TB copy(tb);
    CHECK_EQUAL(fzB(true), copy.isPeriodic());
    CHECK_EQUAL(fzN(0), (index_t)copy.periodicDirection());
    CHECK_EQUAL(fzN(d.n() * NNP), (index_t)copy.size());
}

TEST(space_NEW_tbSwapKeepsPeriodicD0)
{
    const PeriodDesc d = caseC2u(2, 4);
    TB tb = makeTB(d, 0);
    GUARD_TB(tb, 0, d.n(), NNP)
    TB other(npKnots(), npKnots());
    tb.swap(other);
    CHECK_EQUAL(fzB(true), other.isPeriodic());
    CHECK_EQUAL(fzN(0), (index_t)other.periodicDirection());
    CHECK_EQUAL(fzN(d.n() * NNP), (index_t)other.size());
    CHECK_EQUAL(fzB(false), tb.isPeriodic());
    CHECK_EQUAL(fzN(NNP * NNP), (index_t)tb.size());
}

TEST(space_NEW_tbCreatePeriodicKVRefused)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis bPer(stretchedKV(d));
    bPer.setPeriodic();
    if (!guardPeriodic(bPer, d, false))
    {
        CHECK_EQUAL(std::string(""), guardMsg("pre-op guard", bPer, d));
        return;
    }
    CHECK_THROW(Basis::create(bPer.knots(), 2), std::runtime_error);
    CHECK_EQUAL(fzB(true), throwsRT([&] { (void)Basis::create(bPer.knots(), 2); }));
}

gsKnotVector<real_t> npKnotsRefined()
{
    const real_t v[] = {0, 0, 0, 0.125, 0.25, 0.375, 0.5, 0.75, 1, 1, 1};
    return gsKnotVector<real_t>(vec(v, 11), 2);
}

void tbUniformRefine(int dir)
{
    const PeriodDesc d = caseC2u(2, 4), dp = opUniformRefine(d, 1, 1);
    TB tb = makeTB(d, dir);
    GUARD_TB(tb, dir, d.n(), NNP)
    tb.uniformRefine();
    GUARD_TB(tb, dir, dp.n(), NNP + 3)
    const gsMatrix<real_t> U = tgPoints(dp, dir);
    gsMatrix<index_t> act;
    tb.active_into(U, act);
    CHECK_EQUAL(std::string(""), whyOf(activeUniqueInRange(act, tb.size())));
    bool ok;
    const gsMatrix<real_t> got = denseEval(tb, U, dp.n() * (NNP + 3), 0, 0, ok);
    CHECK_EQUAL(std::string(""), whyOf(valuesClose(got, tbRef(dp, dir, U, npKnotsRefined()))));
}
TEST(space_NEW_tbUniformRefineD0) { tbUniformRefine(0); }
TEST(space_NEW_tbUniformRefineD1) { tbUniformRefine(1); }

TEST(space_NEW_tbInsertKnotGhostWindowD0)
{
    const PeriodDesc d = caseC2u(2, 4), dp = opInsert(d, 0.125);
    TB tb = makeTB(d, 0);
    GUARD_TB(tb, 0, d.n(), NNP)
    tb.insertKnot(0.125, 0);
    GUARD_TB(tb, 0, dp.n(), NNP)
    CHECK_EQUAL(fzN(dp.n() * NNP), (index_t)tb.size());
    CHECK_EQUAL(std::string(""), whyOf(knotsClose(
        std::vector<real_t>(tb.knots(0).begin(), tb.knots(0).end()), stretchedKnots(dp), dp.scale())));
}

// =====================================================================
// G5: tensor geometry
// =====================================================================

void tgSetPeriodic(int dir)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = tgCoefs(d, dir), ex = expandTG(cn, d, dir);
    TG g = (0 == dir) ? TG(stretchedKV(d), npKnots(), ex) : TG(npKnots(), stretchedKV(d), ex);
    g.setPeriodic(dir);
    GUARD_TG(g, dir, d.n(), NNP)
    gsMatrix<real_t> got;
    const gsMatrix<real_t> U = tgPoints(d, dir);
    g.eval_into(U, got);
    CHECK_EQUAL(std::string(""), geoClose(got, tgRef(d, dir, cn, U), maxAbs(cn)));
}
TEST(space_B20_tgSetPeriodicD0) { tgSetPeriodic(0); }
TEST(space_B20_tgSetPeriodicD1) { tgSetPeriodic(1); }

void tgEvalFromBasis(int dir)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = tgCoefs(d, dir);
    TG g = makeTG(d, dir, cn);
    GUARD_TG(g, dir, d.n(), NNP)
    gsMatrix<real_t> got;
    const gsMatrix<real_t> U = tgPoints(d, dir);
    g.eval_into(U, got);
    CHECK_EQUAL(std::string(""), geoClose(got, tgRef(d, dir, cn, U), maxAbs(cn)));
}
TEST(space_B11_tgEvalFromBasisD0) { tgEvalFromBasis(0); }
TEST(space_B11_tgEvalFromBasisD1) { tgEvalFromBasis(1); }

GOp gUniformRefine() { return [](TG & g) { g.uniformRefine(); }; }
GOp gInsertKnot(real_t x, int dir) { return [=](TG & g) { g.insertKnot(x, dir); }; }
GOp gElevate(int dir) { return [=](TG & g) { g.degreeElevate((short_t)1, (short_t)dir); }; }
GOp gReverse(int dir) { return [=](TG & g) { g.reverse((unsigned)dir); }; }

TEST(space_NEW_tgUniformRefineD0)
{
    const PeriodDesc d = caseC2u(2, 4);
    tgRow(d, 0, gUniformRefine(), opUniformRefine(d, 1, 1), NNP + 3);
}
TEST(info_twin_tgUniformRefineD0) { twinTG("tgUniformRefine", caseC2u(2, 4), 0, gUniformRefine()); }

TEST(space_NEW_tgUniformRefineD1)
{
    const PeriodDesc d = caseC2u(2, 4);
    tgRow(d, 1, gUniformRefine(), opUniformRefine(d, 1, 1), NNP + 3);
}
TEST(info_twin_tgUniformRefineD1) { twinTG("tgUniformRefine", caseC2u(2, 4), 1, gUniformRefine()); }

TEST(space_NEW_tgInsertKnotGhostWindowD0)
{
    const PeriodDesc d = caseC2u(2, 4);
    tgRow(d, 0, gInsertKnot(0.125, 0), opInsert(d, 0.125), NNP);
}
TEST(info_twin_tgInsertKnotGhostWindowD0) { twinTG("tgInsertKnotGhostWindow", caseC2u(2, 4), 0, gInsertKnot(0.125, 0)); }

TEST(space_NEW_tgDegreeElevateD0)
{
    const PeriodDesc d = caseC2u(2, 4);
    tgRow(d, 0, gElevate(0), opDegreeElevate(d, 1), NNP);
}
TEST(info_twin_tgDegreeElevateD0) { twinTG("tgDegreeElevate", caseC2u(2, 4), 0, gElevate(0)); }

void tgReverse(int dir)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = tgCoefs(d, dir);
    const TG g = makeTG(d, dir, cn);
    GUARD_TG(g, dir, d.n(), NNP)
    const gsMatrix<real_t> U = tgPoints(d, dir);

    TG gp = g;
    gp.reverse((unsigned)dir);
    GUARD_TG(gp, dir, d.n(), NNP)
    const gsMatrix<real_t> pr = gp.parameterRange();
    CHECK_EQUAL(std::string(""), whyOf(scalarClose(pr(dir, 0), d.a, d.scale())));
    CHECK_EQUAL(std::string(""), whyOf(scalarClose(pr(dir, 1), d.b, d.scale())));
    gsMatrix<real_t> Ur = U, got;
    for (index_t i = 0; i != U.cols(); ++i)
        Ur(dir, i) = d.a + d.b - U(dir, i);
    gp.eval_into(U, got);
    CHECK_EQUAL(std::string(""), geoClose(got, tgRef(d, dir, cn, Ur), maxAbs(cn)));

    TG gq = g;
    gq.reverse((unsigned)(1 - dir));
    GUARD_TG(gq, dir, d.n(), NNP)
    Ur = U;
    for (index_t i = 0; i != U.cols(); ++i)
        Ur(1 - dir, i) = 1 - U(1 - dir, i);
    gq.eval_into(U, got);
    CHECK_EQUAL(std::string(""), geoClose(got, tgRef(d, dir, cn, Ur), maxAbs(cn)));
}
TEST(space_NEW_tgReverseD0) { tgReverse(0); }
TEST(info_twin_tgReverseD0) { twinTG("tgReverse", caseC2u(2, 4), 0, gReverse(0), true); }
TEST(space_NEW_tgReverseD1) { tgReverse(1); }
TEST(info_twin_tgReverseD1) { twinTG("tgReverse", caseC2u(2, 4), 1, gReverse(1), true); }

TEST(space_NEW_tgSwapDirectionsD0)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = tgCoefs(d, 0);
    TG g = makeTG(d, 0, cn);
    GUARD_TG(g, 0, d.n(), NNP)
    g.swapDirections(0, 1);
    GUARD_TG(g, 1, d.n(), NNP)
    const gsMatrix<real_t> U = tgPoints(d, 0);
    gsMatrix<real_t> Us(2, U.cols());
    Us.row(0) = U.row(1);
    Us.row(1) = U.row(0);
    gsMatrix<real_t> got;
    g.eval_into(Us, got);
    CHECK_EQUAL(std::string(""), geoClose(got, tgRef(d, 0, cn, U), maxAbs(cn)));
}

TEST(space_NEW_tgSliceNonPeriodicParamD0)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = tgCoefs(d, 0);
    const TG g = makeTG(d, 0, cn);
    GUARD_TG(g, 0, d.n(), NNP)
    Crv crv;
    g.slice(1, real_t(0.375), crv);
    GUARD_CRV(crv, d)
    const gsMatrix<real_t> P = samplePoints(d);
    gsMatrix<real_t> U(2, P.cols());
    U.row(0) = P.row(0);
    U.row(1).setConstant(0.375);
    gsMatrix<real_t> got;
    crv.eval_into(P, got);
    CHECK_EQUAL(std::string(""), geoClose(got, tgRef(d, 0, cn, U), maxAbs(cn)));
}

TEST(space_NEW_tgSlicePeriodicParamD0)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = tgCoefs(d, 0);
    const TG g = makeTG(d, 0, cn);
    GUARD_TG(g, 0, d.n(), NNP)
    Crv crv;
    g.slice(0, real_t(0.375), crv);
    CHECK_EQUAL(fzB(false), crv.isPeriodic());
    CHECK_EQUAL(fzN(NNP), (index_t)crv.coefs().rows());
    const real_t v[] = {0, 0.1, 0.25, 0.4, 0.5, 0.8, 1};
    gsMatrix<real_t> V(1, 7), U(2, 7);
    for (int i = 0; i != 7; ++i)
    {
        V(0, i) = v[i];
        U(0, i) = 0.375;
        U(1, i) = v[i];
    }
    gsMatrix<real_t> got;
    crv.eval_into(V, got);
    CHECK_EQUAL(std::string(""), geoClose(got, tgRef(d, 0, cn, U), maxAbs(cn)));
}

TEST(space_NEW_tgSplitAtPeriodicRefusedD0)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = tgCoefs(d, 0);
    TG g = makeTG(d, 0, cn);
    GUARD_TG(g, 0, d.n(), NNP)
    const TG g0 = g;
    TG l, r;
    CHECK_THROW(g.splitAt(0, real_t(0.5), l, r), std::runtime_error);
    CHECK_EQUAL(std::string(""), geoClose(g.coefs(), g0.coefs(), maxAbs(cn), 0));
    CHECK_EQUAL(std::string(""), tgGuardN(g, 0, d.n(), NNP));
    CHECK_EQUAL(std::string(""), whyOf(knotsClose(tgKnots(g, 0), tgKnots(g0, 0), d.scale())));
    g.splitAt(1, real_t(0.5), l, r);
    CHECK_EQUAL(fzB(true), l.basis().isPeriodic() && r.basis().isPeriodic());
    CHECK_EQUAL(fzN(0), (index_t)l.basis().periodicDirection());
    CHECK_EQUAL(fzN(0), (index_t)r.basis().periodicDirection());
    const gsMatrix<real_t> P = samplePoints(d);
    gsMatrix<real_t> U(2, P.cols()), got;
    U.row(0) = P.row(0);
    U.row(1).setConstant(0.25);
    l.eval_into(U, got);
    CHECK_EQUAL(std::string(""), geoClose(got, tgRef(d, 0, cn, U), maxAbs(cn)));
    U.row(1).setConstant(0.75);
    r.eval_into(U, got);
    CHECK_EQUAL(std::string(""), geoClose(got, tgRef(d, 0, cn, U), maxAbs(cn)));
}

// =====================================================================
// G6: multipatch topology of a periodic cylinder
// =====================================================================

/// Cylinder rows (cos theta, sin theta, z) in the index order of direction 0.
gsMatrix<real_t> cylCoefs(const PeriodDesc & d)
{
    const index_t n = d.n();
    gsMatrix<real_t> c(n * NNP, 3);
    for (index_t ip = 0; ip != n; ++ip)
        for (index_t iq = 0; iq != NNP; ++iq)
        {
            const real_t th = 2 * PI * real_t(ip) / real_t(n);
            c.row(gidx(0, ip, iq, n, NNP)) << math::cos(th), math::sin(th), real_t(iq) / real_t(NNP - 1);
        }
    return c;
}

TEST(space_NEW_mpCylinderCtor)
{
    const PeriodDesc d = caseC2u(2, 4);
    const TG g = makeTG(d, 0, cylCoefs(d));
    GUARD_TG(g, 0, d.n(), NNP)
    gsMultiPatch<real_t> mp(g);
    CHECK_EQUAL((size_t)fzN(2), mp.nBoundary());
    CHECK_EQUAL((size_t)0, mp.nInterfaces());
    CHECK_EQUAL(fzB(true), mp.isBoundary(0, boxSide(boundary::south)));
    CHECK_EQUAL(fzB(true), mp.isBoundary(0, boxSide(boundary::north)));
    CHECK_EQUAL(fzB(false), mp.isBoundary(0, boxSide(boundary::west)));
    CHECK_EQUAL(fzB(false), mp.isBoundary(0, boxSide(boundary::east)));
}

TEST(space_NEW_mpCylinderComputeTopology)
{
    const PeriodDesc d = caseC2u(2, 4);
    const TG g = makeTG(d, 0, cylCoefs(d));
    GUARD_TG(g, 0, d.n(), NNP)
    gsMultiPatch<real_t> mp(g);
    mp.computeTopology();
    CHECK_EQUAL((size_t)fzN(2), mp.nBoundary());
    CHECK_EQUAL((size_t)fzN(0), mp.nInterfaces());
}

TEST(info_clamped_mpCylinderTopology)
{
    try
    {
        const index_t n0 = 6;
        gsMatrix<real_t> c(n0 * NNP, 3);
        for (index_t i0 = 0; i0 != n0; ++i0)
            for (index_t i1 = 0; i1 != NNP; ++i1)
            {
                const real_t th = 2 * PI * real_t(i0 % (n0 - 1)) / real_t(n0 - 1);
                c.row(i0 + i1 * n0) << math::cos(th), math::sin(th), real_t(i1) / real_t(NNP - 1);
            }
        const gsKnotVector<real_t> kv0(real_t(0), real_t(1), 3u, 3, 1, 2);
        const TG g(kv0, npKnots(), c);
        gsMultiPatch<real_t> mp(g);
        mp.computeTopology();
        gsInfo << "clamped closed cylinder: nBoundary=" << mp.nBoundary()
               << " nInterfaces=" << mp.nInterfaces() << "\n";
    }
    catch (std::exception & e)
    {
        gsInfo << "clamped closed cylinder: exception: " << e.what() << "\n";
    }
}

// =====================================================================
// G7: XML round trips
// =====================================================================

void xmlBasisRoundTrip(const PeriodDesc & d, const char * name)
{
    Basis b(stretchedKV(d));
    b.setPeriodic();
    if (!guardPeriodic(b, d, false))
    {
        CHECK_EQUAL(std::string(""), guardMsg("pre-op guard", b, d));
        return;
    }
    const std::string path = xmlPath(name);
    {
        gsFileData<> fd;
        fd << b;
        fd.save(path);
    }
    gsFileData<> in(path);
    memory::unique_ptr<Basis> r = in.getFirst<Basis>();
    std::remove(path.c_str());
    CHECK_EQUAL(fzB(true), NULL != r.get());
    if (!r)
        return;
    CHECK_EQUAL(fzB(true), r->isPeriodic());
    CHECK_EQUAL(fzN(d.n()), (index_t)r->size());
    CHECK_EQUAL(std::string(""), whyOf(knotsClose(knotsOf(*r), stretchedKnots(d), d.scale())));
    if (r->isPeriodic() && r->size() == d.n())
        CHECK_EQUAL(std::string(""), whyOf(basisVsR1(*r, d)));
}
TEST(space_B21_xmlBasisRoundTripCuni)  { xmlBasisRoundTrip(caseC2u(2, 4), "basisCuni"); }
TEST(space_B21_xmlBasisRoundTripCnonu) { xmlBasisRoundTrip(caseC3n(), "basisCnonu"); }

TEST(space_B21_xmlCurveRoundTripCuni)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = crvCoefs(d);
    Crv c(stretchedKV(d), expandRows(cn, d.n(), d.m()), true);
    GUARD_CRV(c, d)
    const std::string path = xmlPath("curveCuni");
    {
        gsFileData<> fd;
        fd << c;
        fd.save(path);
    }
    gsFileData<> in(path);
    memory::unique_ptr<Crv> r = in.getFirst<Crv>();
    std::remove(path.c_str());
    CHECK_EQUAL(fzB(true), NULL != r.get());
    if (!r)
        return;
    CHECK_EQUAL(fzB(true), r->isPeriodic());
    CHECK_EQUAL(fzN(d.n()), (index_t)r->coefs().rows());
    CHECK_EQUAL(std::string(""), whyOf(knotsClose(knotsVec(*r), stretchedKnots(d), d.scale())));
    CHECK_EQUAL(std::string(""), curveVsRef(*r, d, cn, samplePoints(d)));
}

void xmlTensorRoundTrip(int dir, const char * name)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = tgCoefs(d, dir);
    const TG g = makeTG(d, dir, cn);
    GUARD_TG(g, dir, d.n(), NNP)
    const std::string path = xmlPath(name);
    {
        gsFileData<> fd;
        fd << g;
        fd.save(path);
    }
    gsFileData<> in(path);
    memory::unique_ptr<TG> r = in.getFirst<TG>();
    std::remove(path.c_str());
    CHECK_EQUAL(fzB(true), NULL != r.get());
    if (!r)
        return;
    CHECK_EQUAL(fzB(true), r->basis().isPeriodic());
    CHECK_EQUAL(fzN(dir), (index_t)r->basis().periodicDirection());
    CHECK_EQUAL(fzN(d.n() * NNP), (index_t)r->coefs().rows());
    CHECK_EQUAL(std::string(""), whyOf(knotsClose(tgKnots(*r, dir), stretchedKnots(d), d.scale())));
    const std::vector<real_t> kn = tgKnots(*r, 1 - dir);
    const gsKnotVector<real_t> npk = npKnots();
    const std::vector<real_t> np(npk.begin(), npk.end());
    CHECK_EQUAL(std::string(""), whyOf(knotsClose(kn, np, 1)));
    gsMatrix<real_t> got;
    const gsMatrix<real_t> U = tgPoints(d, dir);
    r->eval_into(U, got);
    CHECK_EQUAL(std::string(""), geoClose(got, tgRef(d, dir, cn, U), maxAbs(cn)));
}
TEST(space_B21_xmlTensorGeoRoundTripD0) { xmlTensorRoundTrip(0, "tensorD0"); }
TEST(space_B21_xmlTensorGeoRoundTripD1) { xmlTensorRoundTrip(1, "tensorD1"); }

TEST(space_NEW_xmlNoAttributeBasis)
{
    const std::string path = xmlPath("noAttrBasis");
    {
        std::ofstream f(path.c_str());
        f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<xml>\n"
          << "<Basis type=\"BSplineBasis\"><KnotVector degree=\"2\">"
          << "-0.5 -0.25 0 0.25 0.5 0.75 1 1.25 1.5</KnotVector></Basis>\n</xml>\n";
    }
    gsFileData<> in(path);
    memory::unique_ptr<Basis> r = in.getFirst<Basis>();
    std::remove(path.c_str());
    CHECK_EQUAL(fzB(true), NULL != r.get());
    if (!r)
        return;
    CHECK_EQUAL(fzB(false), r->isPeriodic());
    CHECK_EQUAL(fzN(6), (index_t)r->size());
    const PeriodDesc d = caseC2u(2, 4);
    CHECK_EQUAL(std::string(""), whyOf(knotsClose(knotsOf(*r), stretchedKnots(d), d.scale())));
}

TEST(space_NEW_xmlNoAttributeLegacyCurve)
{
    const PeriodDesc d = caseC2u(2, 4);
    const gsMatrix<real_t> cn = crvCoefs(d), ex = expandRows(cn, d.n(), d.m());
    const std::string path = xmlPath("noAttrCurve");
    {
        std::ofstream f(path.c_str());
        f << std::setprecision(17);
        f << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<xml>\n<Geometry type=\"BSpline\">\n"
          << "<Basis type=\"BSplineBasis\"><KnotVector degree=\"2\">"
          << "-0.5 -0.25 0 0.25 0.5 0.75 1 1.25 1.5</KnotVector></Basis>\n<coefs geoDim=\"2\">\n";
        for (index_t i = 0; i != ex.rows(); ++i)
            f << ex(i, 0) << " " << ex(i, 1) << "\n";
        f << "</coefs>\n</Geometry>\n</xml>\n";
    }
    gsFileData<> in(path);
    memory::unique_ptr<Crv> r = in.getFirst<Crv>();
    std::remove(path.c_str());
    CHECK_EQUAL(fzB(true), NULL != r.get());
    if (!r)
        return;
    CHECK_EQUAL(fzB(false), r->isPeriodic());
    CHECK_EQUAL(fzN(d.n() + d.m()), (index_t)r->coefs().rows());
    CHECK_EQUAL(std::string(""), curveVsRef(*r, d, cn, samplePoints(d)));
}

// =====================================================================
// Controls: the references (a) and (b) and the guards
// =====================================================================

TEST(control_helperA_crvVsNonPeriodic)
{
    std::vector<PeriodDesc> cs;
    cs.push_back(caseC2u(2, 4));
    cs.push_back(caseC2a());
    cs.push_back(caseC3n());
    cs.push_back(caseC4i());
    cs.push_back(caseC4o());
    const std::string fails = runCases(cs, [](const PeriodDesc & d) -> std::string
    {
        const gsMatrix<real_t> cn = crvCoefs(d);
        const Crv c(stretchedKV(d), expandRows(cn, d.n(), d.m()));
        const gsMatrix<real_t> U = samplePoints(d);
        return curveVsRef(c, d, cn, U) + curveVsRef(c, d, cn, U, 1) + curveVsRef(c, d, cn, U, 2);
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(control_helperA_crvNegPerturbedRow)
{
    const PeriodDesc d = caseC2u(2, 4);
    gsMatrix<real_t> cn = crvCoefs(d);
    const Crv c(stretchedKV(d), expandRows(cn, d.n(), d.m()));
    cn(1, 0) += 1e-3;
    CHECK(!curveVsRef(c, d, cn, samplePoints(d)).empty());
}

TEST(control_helperA_crvNegGuardAndUnchanged)
{
    const PeriodDesc d = caseC2u(2, 4);
    Crv c(stretchedKV(d), expandRows(crvCoefs(d), d.n(), d.m()));
    CHECK(!crvGuard(c, d).empty());
    const CrvSnap s0 = crvSnap(c);
    CHECK(crvUnchanged(s0, c).empty());
    c.insertKnot(0.125);
    CHECK(!crvUnchanged(s0, c).empty());
}

TEST(control_helperB_tgVsNonPeriodic)
{
    const PeriodDesc d = caseC2u(2, 4);
    for (int dir = 0; dir != 2; ++dir)
    {
        const gsMatrix<real_t> cn = tgCoefs(d, dir), ex = expandTG(cn, d, dir);
        const TG g = (0 == dir) ? TG(stretchedKV(d), npKnots(), ex) : TG(npKnots(), stretchedKV(d), ex);
        gsMatrix<real_t> got;
        const gsMatrix<real_t> U = tgPoints(d, dir);
        g.eval_into(U, got);
        CHECK_EQUAL(std::string(""), geoClose(got, tgRef(d, dir, cn, U), maxAbs(cn)));
    }
}

TEST(control_helperB_tgNegPerturbedRow)
{
    const PeriodDesc d = caseC2u(2, 4);
    for (int dir = 0; dir != 2; ++dir)
    {
        gsMatrix<real_t> cn = tgCoefs(d, dir);
        const gsMatrix<real_t> ex = expandTG(cn, d, dir);
        const TG g = (0 == dir) ? TG(stretchedKV(d), npKnots(), ex) : TG(npKnots(), stretchedKV(d), ex);
        cn(3, 1) += 1e-3;
        gsMatrix<real_t> got;
        const gsMatrix<real_t> U = tgPoints(d, dir);
        g.eval_into(U, got);
        CHECK(!geoClose(got, tgRef(d, dir, cn, U), maxAbs(cn)).empty());
    }
}

TEST(control_helperB_tbPartitionOfUnity)
{
    const PeriodDesc d = caseC2u(2, 4);
    for (int dir = 0; dir != 2; ++dir)
    {
        const Verdict v = partitionOfUnity(tbRef(d, dir, tgPoints(d, dir), npKnots()));
        CHECK_EQUAL(std::string(""), whyOf(v));
    }
    gsMatrix<real_t> bad = tbRef(d, 0, tgPoints(d, 0), npKnots());
    bad.row(0) *= 1.01;
    CHECK(!partitionOfUnity(bad).ok);
}

TEST(control_helperB_tbGuardsNegative)
{
    const PeriodDesc d = caseC2u(2, 4);
    const TB tb(stretchedKV(d), npKnots());
    CHECK(!tbGuardN(tb, 0, d.n(), NNP).empty());
}

} // SUITE
