/** @file gsPeriodicBSplineBasis_test.cpp

    @brief Public behaviour of one-dimensional periodic B-spline bases.

    Expected values come from period descriptions and from non-periodic
    reference bases (see gsPeriodicTestUtils.h), never from the basis under
    test.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
**/

#include "gismo_unittest.h"
#include "gsPeriodicTestUtils.h"

using namespace gismo::periodic_test;

namespace {

typedef gsBSplineBasis<real_t> Basis;
typedef std::function<void(Basis &)> Op;

/// Periodic basis from the stretched knots of d.
Basis convert(const PeriodDesc & d)
{
    Basis b(stretchedKV(d));
    b.setPeriodic();
    return b;
}

/// Clamped knot vector of d (one span per unique breakpoint), converted.
Basis clampedConv(const PeriodDesc & d)
{
    std::vector<real_t> ub(d.s.begin(), d.s.end());
    ub.push_back(d.b);
    const index_t N = uniqueSorted(ub).size() - 1;
    gsKnotVector<real_t> kv(d.a, d.b, (unsigned)(N - 1), d.p + 1, 1, d.p);
    Basis b(kv);
    b.setPeriodic();
    return b;
}

/// Gauss points only (no knots, no domain ends), p+1 per span.
gsMatrix<real_t> gaussOnly(const PeriodDesc & d)
{
    std::vector<real_t> br(d.s.begin(), d.s.end());
    br.push_back(d.b);
    br = uniqueSorted(br);
    gsGaussRule<real_t> gr(d.p + 1);
    gsMatrix<real_t> nodes;
    gsVector<real_t> w;
    gr.mapToAll(br, nodes, w);
    return nodes;
}

std::string literalKnots(const PeriodDesc & d, const real_t * lit, size_t k)
{
    return whyOf(knotsClose(stretchedKnots(d), vec(lit, k), d.scale()));
}

std::string sEquals(const PeriodDesc & got, const real_t * lit, size_t k, short_t p)
{
    std::string w = whyOf(knotsClose(got.s, vec(lit, k), got.scale()));
    if (got.p != p)
        w += " degree " + num(got.p) + " != " + num(p);
    return w.empty() ? w : got.tag + ": " + w + "; ";
}

// --- operations -----------------------------------------------------

Op uRef(int k, int mul)   { return [=](Basis & b) { b.uniformRefine(k, mul); }; }
Op insK(real_t x, int m = 1) { return [=](Basis & b) { b.insertKnot(x, m); }; }
Op remK(real_t x)         { return [=](Basis & b) { b.removeKnot(x); }; }
Op dElev()                { return [](Basis & b) { b.degreeElevate((short_t)1); }; }
Op dInc()                 { return [](Basis & b) { b.degreeIncrease((short_t)1); }; }
Op dRed()                 { return [](Basis & b) { b.degreeReduce((short_t)1); }; }
Op dDec()                 { return [](Basis & b) { b.degreeDecrease((short_t)1); }; }
Op elevCont()             { return [](Basis & b) { b.elevateContinuity(1); }; }
Op redCont()              { return [](Basis & b) { b.reduceContinuity(1); }; }
Op coarsen1()             { return [](Basis & b) { b.uniformCoarsen(1); }; }
Op rev()                  { return [](Basis & b) { b.reverse(); }; }

Op insKnotsVec(real_t x, real_t y)
{
    return [=](Basis & b)
    {
        std::vector<std::vector<real_t> > kk(1);
        kk[0].push_back(x);
        kk[0].push_back(y);
        b.insertKnots(kk);
    };
}

Op refineElem0()
{
    return [](Basis & b) { b.refineElements(std::vector<index_t>(1, 0)); };
}

Op directInsert(real_t x)
{
    return [=](Basis & b) { b.knots().insert(x); };
}

Op refineK(const Basis & other)
{
    return [&other](Basis & b) { b.refine_k(other, 1); };
}

// --- drivers ----------------------------------------------------------

/// Converts d; returns "" or the pre-operation guard message.
std::string preGuard(const PeriodDesc & d, Basis & b)
{
    try { b.setPeriodic(); }
    catch (std::exception & e)
    { return std::string("pre-op guard: conversion threw: ") + e.what(); }
    return guardPeriodic(b, d, false) ? std::string() : guardMsg("pre-op guard", b, d);
}

/// Refinement row: convert d, apply op, expect the periodic space of dp.
void rowSuccess(const PeriodDesc & d, const Op & op, const PeriodDesc & dp)
{
    Basis b(stretchedKV(d));
    const std::string pre = preGuard(d, b);
    CHECK_EQUAL(std::string(""), pre);
    if (!pre.empty())
        return;
    try { op(b); }
    catch (std::exception & e)
    {
        CHECK_EQUAL(std::string(""), std::string("operation threw: ") + e.what());
        return;
    }
    const std::string post = guardPeriodic(b, dp, false) ? std::string()
                                                          : guardMsg("post-op guard", b, dp);
    CHECK_EQUAL(std::string(""), post);
    if (!post.empty())
        return;
    CHECK_EQUAL(fzN(dp.p), (index_t)b.degree());
    try
    {
        const gsMatrix<real_t> U = samplePoints(dp);
        bool ok1, ok2;
        const gsMatrix<real_t> got = denseEval(b, U, dp.n(), 0, 0, ok1);
        const gsMatrix<real_t> ref = denseEval(makeR1(dp), U, dp.n(), dp.n(), 0, ok2);
        CHECK(ok1);
        CHECK_EQUAL(std::string(""), whyOf(valuesClose(got, ref)));
    }
    catch (std::exception & e)
    {
        CHECK_EQUAL(std::string(""), std::string("evaluation threw: ") + e.what());
    }
}

/// Refusal of op on b: throws std::runtime_error and leaves b unchanged.
void refuse(Basis & b, const Op & op)
{
    const BasisSnapshot s0 = snapshot(b);
    CHECK_THROW(op(b), std::runtime_error);
    CHECK_EQUAL(std::string(""), whyOf(unchanged(s0, b)));
}

/// Refusal row: each op is tried on a fresh periodic basis of d.
void rowRefusal(const PeriodDesc & d, const Op & op1, const Op & op2 = Op())
{
    for (int i = 0; i != 2; ++i)
    {
        const Op & op = (0 == i) ? op1 : op2;
        if (!op)
            continue;
        Basis b(stretchedKV(d));
        const std::string pre = preGuard(d, b);
        CHECK_EQUAL(std::string(""), pre);
        if (!pre.empty())
            return;
        refuse(b, op);
    }
}

/// Reverse(): function j at u is the R1 function (m-1-j) mod n at a+b-u,
/// and the same as the R1 functions of the reversed description.
void reverseValues(const PeriodDesc & d)
{
    Basis b(stretchedKV(d));
    const std::string pre = preGuard(d, b);
    CHECK_EQUAL(std::string(""), pre);
    if (!pre.empty())
        return;
    try { b.reverse(); }
    catch (std::exception & e)
    {
        CHECK_EQUAL(std::string(""), std::string("operation threw: ") + e.what());
        return;
    }
    const std::string post = guardPeriodic(b, d, false) ? std::string()
                                                         : guardMsg("post-op guard", b, d);
    CHECK_EQUAL(std::string(""), post);
    if (!post.empty())
        return;
    const gsMatrix<real_t> U = samplePoints(d);
    gsMatrix<real_t> Ur = U;
    for (index_t i = 0; i != U.cols(); ++i)
        Ur(0, i) = d.a + d.b - U(0, i);
    bool ok;
    const gsMatrix<real_t> got = denseEval(b, U, d.n(), 0, 0, ok);
    const gsMatrix<real_t> ref = r1Dense(d, Ur, 0, ok);
    gsMatrix<real_t> expected(d.n(), U.cols());
    for (index_t j = 0; j != d.n(); ++j)
        expected.row(j) = ref.row(floorMod(d.m() - 1 - j, d.n()));
    CHECK_EQUAL(std::string(""), whyOf(valuesClose(got, expected)));

    const PeriodDesc dr = opReverse(d);
    CHECK_EQUAL(fzB(true), constructible(dr));
    const gsMatrix<real_t> refRev = r1Dense(dr, U, 0, ok);
    CHECK_EQUAL(std::string(""), whyOf(valuesClose(got, refRev)));
}

/// Maximum R1 value error and structural data of a probe case.
void probe(const char * name, const PeriodDesc & d)
{
    try
    {
        Basis b(stretchedKV(d));
        b.setPeriodic();
        gsInfo << name << ": periodic=" << b.isPeriodic() << " size=" << b.size()
               << " numCrossingFunctions=" << b.numCrossingFunctions();
        const Verdict v = basisVsR1(b, d);
        gsInfo << " maxR1ValueError=" << v.value
               << " selfCheck=" << (selfCheck(d, stretchedKnots(d)).ok ? "ok" : "bad")
               << (v.ok ? "" : " (" + v.why + ")") << "\n";
    }
    catch (std::exception & e)
    {
        gsInfo << name << ": exception: " << e.what() << "\n";
    }
}

} // namespace

SUITE(gsPeriodicBSplineBasis_test)
{

// =====================================================================
// Controls
// =====================================================================

TEST(control_genExamples_all)
{
    std::string fails;

    // Generator output, literal vectors.
    const real_t c2uP2N4[] = {-.5, -.25, 0, .25, .5, .75, 1, 1.25, 1.5};
    fails += literalKnots(caseC2u(2, 4), c2uP2N4, 9);
    const real_t c2a[] = {-2.75, -1.875, -1, -0.125, 0.75, 1.625, 2.5, 3.375, 4.25};
    fails += literalKnots(caseC2a(), c2a, 9);
    const real_t c3n[] = {-1, -.5, 0, .5, 1.5, 2, 3, 3.5, 4, 4.5, 5.5};
    fails += literalKnots(caseC3n(), c3n, 11);
    const real_t c3y[] = {-1, -.5, 0, .5, 1.5, 2.5, 3, 3.5, 4, 4.5, 5.5};
    fails += literalKnots(caseC3y(), c3y, 11);
    const real_t c4i[] = {-2, -1, 0, 1, 1, 2, 3, 4, 5, 5};
    fails += literalKnots(caseC4i(), c4i, 10);
    const real_t c4o[] = {-2, -1, 0, 1, 2, 3, 3, 4, 5, 6, 7, 8};
    fails += literalKnots(caseC4o(), c4o, 12);
    const real_t e2[] = {-.25, -.25, 0, 0, .25, .25, .5, .5, .75, .75, 1, 1, 1.25, 1.25};
    fails += literalKnots(caseE2(), e2, 14);
    const real_t c1kL4P2[] = {-1, 0, 0, 1, 2, 3, 4, 4, 5};
    fails += literalKnots(caseC1kL4(2), c1kL4P2, 9);
    const real_t c1kP2N4[] = {-.25, 0, 0, .25, .5, .75, 1, 1, 1.25};
    fails += literalKnots(caseC1k(2, 4), c1kP2N4, 9);
    const real_t c5a[] = {-1.5, -1, -.5, 0, .5, 1, 1.5, 2, 2.5};
    fails += literalKnots(caseC5a(), c5a, 9);
    const real_t c5b[] = {-3, -2, -1, 0, 1, 2, 3, 4, 5, 6};
    fails += literalKnots(caseC5b(), c5b, 10);

    // Worked example with seam multiplicity 2: p=2, s=[0,0,1,2,3] on [0,4].
    {
        const real_t s[] = {0, 0, 1, 2, 3};
        const PeriodDesc d = mkDesc("C1kL4P2alt", 2, 0, 4, vec(s, 5));
        fails += literalKnots(d, c1kL4P2, 9);
    }

    // Uniform cases: t_k = (k-p)/N, for C2u and C1c (mult 1).
    {
        std::vector<PeriodDesc> u = c2uCases(), c = c1cCases();
        u.insert(u.end(), c.begin(), c.end());
        for (size_t i = 0; i != u.size(); ++i)
        {
            const index_t N = u[i].n(), p = u[i].p;
            std::vector<real_t> t;
            for (index_t k = 0; k <= N + p + p; ++k)
                t.push_back(real_t(k - p) / real_t(N));
            fails += whyOf(knotsClose(stretchedKnots(u[i]), t, u[i].scale()));
        }
    }

    // selfCheck is ok on every generated vector, also the refused ones.
    {
        std::vector<PeriodDesc> all = constructibleCases(), k = c1kCases(), c = c1cCases();
        all.insert(all.end(), k.begin(), k.end());
        all.insert(all.end(), c.begin(), c.end());
        all.push_back(caseC5a());
        all.push_back(caseC5b());
        for (size_t i = 0; i != all.size(); ++i)
        {
            const Verdict v = selfCheck(all[i], stretchedKnots(all[i]));
            if (!v.ok)
                fails += all[i].tag + " selfCheck: " + v.why + "; ";
        }
    }
    CHECK_EQUAL(std::string(""), fails);

    CHECK_EQUAL(fzN(-1), floorDiv(-1, 5));
    CHECK_EQUAL(fzN(4), floorMod(-1, 5));
    CHECK_EQUAL(fzN(-2), floorDiv(-6, 5));
    CHECK_EQUAL(fzN(4), floorMod(-6, 5));
    CHECK_EQUAL(fzN(1), floorDiv(5, 5));

    CHECK_EQUAL(fzB(false), constructible(caseC5a()));
    CHECK_EQUAL(fzB(false), constructible(caseC5b()));
    CHECK_EQUAL(fzB(false), constructible(caseC2u(4, 4)));
    {
        std::vector<PeriodDesc> cs = constructibleCases();
        bool all = true;
        for (size_t i = 0; i != cs.size(); ++i)
            all = all && constructible(cs[i]);
        CHECK_EQUAL(fzB(true), all);
    }

    // Expected breakpoint sets of the op... helpers.
    std::string ops;
    {
        const PeriodDesc u24 = caseC2u(2, 4);
        const real_t r11[] = {0, .125, .25, .375, .5, .625, .75, .875};
        ops += sEquals(opUniformRefine(u24, 1, 1), r11, 8, 2);
        const real_t r12[] = {0, .125, .125, .25, .375, .375, .5, .625, .625, .75, .875, .875};
        ops += sEquals(opUniformRefine(u24, 1, 2), r12, 12, 2);
        const real_t c3nRef[] = {0, .25, .5, 1, 1.5, 1.75, 2, 2.5, 3, 3.25, 3.5, 3.75};
        ops += sEquals(opUniformRefine(caseC3n(), 1, 1), c3nRef, 12, 2);
        const real_t c4iRef[] = {0, .5, 1, 1, 1.5, 2, 2.5, 3, 3.5};
        ops += sEquals(opUniformRefine(caseC4i(), 1, 1), c4iRef, 9, 2);

        const real_t i25[] = {0, .5, 1.5, 2, 2.5, 3, 3.5};
        ops += sEquals(opInsert(caseC3n(), 2.5), i25, 7, 2);
        const real_t i025[] = {0, .25, .5, 1.5, 2, 3, 3.5};
        ops += sEquals(opInsert(caseC3n(), 0.25), i025, 7, 2);
        const real_t i375[] = {0, .5, 1.5, 2, 3, 3.5, 3.75};
        ops += sEquals(opInsert(caseC3n(), 3.75), i375, 7, 2);
        const real_t c4oIns[] = {0, 1, 2, 3, 3, 3.5, 4, 5};
        ops += sEquals(opInsert(caseC4o(), 3.5), c4oIns, 8, 2);
        const real_t i2 [] = {0, .25, .5, 1.5, 2, 2.5, 3, 3.5};
        ops += sEquals(opInsert(opInsert(caseC3n(), 0.25), 2.5), i2, 8, 2);

        const real_t rm20[] = {0, .5, 1.5, 3, 3.5};
        ops += sEquals(opRemove(caseC3n(), 2.0), rm20, 5, 2);
        const real_t rm05[] = {0, 1.5, 2, 3, 3.5};
        ops += sEquals(opRemove(caseC3n(), 0.5), rm05, 5, 2);
        const real_t rm1[] = {0, 1, 2, 3};
        ops += sEquals(opRemove(caseC4i(), 1.0), rm1, 4, 2);

        const real_t de[] = {0, 0, .25, .25, .5, .5, .75, .75};
        ops += sEquals(opDegreeElevate(u24, 1), de, 8, 3);
        const real_t deC3n[] = {0, 0, .5, .5, 1.5, 1.5, 2, 2, 3, 3, 3.5, 3.5};
        ops += sEquals(opDegreeElevate(caseC3n(), 1), deC3n, 12, 3);
        const real_t quarters[] = {0, .25, .5, .75};
        ops += sEquals(opDegreeIncrease(u24, 1), quarters, 4, 3);
        ops += sEquals(opDegreeReduce(caseE2(), 1), quarters, 4, 2);
        const real_t only0[] = {0};
        ops += sEquals(opDegreeReduce(u24, 1), only0, 1, 1);
        ops += sEquals(opDegreeDecrease(u24, 1), quarters, 4, 1);
        ops += sEquals(opElevateContinuity(caseE2(), 1), quarters, 4, 3);
        ops += sEquals(opElevateContinuity(u24, 1), only0, 1, 2);
        ops += sEquals(opReduceContinuity(u24, 1), de, 8, 2);
        const real_t half[] = {0, .5};
        ops += sEquals(opCoarsen(caseC2u(2, 8), 1), quarters, 4, 2);
        ops += sEquals(opCoarsen(u24, 1), half, 2, 2);
        const real_t revC3n[] = {0, .5, 1, 2, 2.5, 3.5};
        ops += sEquals(opReverse(caseC3n()), revC3n, 6, 2);
    }
    CHECK_EQUAL(std::string(""), ops);
}

TEST(control_genSelfCheckNeg_C3x)
{
    const real_t c3x[] = {-1.5, -.5, 0, .5, 1.5, 2, 3, 3.5, 4, 4.5, 5.5};
    CHECK_EQUAL(fzB(false), selfCheck(caseC3n(), vec(c3x, 11)).ok);

    const real_t s31[] = {-1, 0, 0, 1, 2, 3, 4, 5, 6};
    const real_t s[] = {0, 0, 1, 2, 3};
    const PeriodDesc d = mkDesc("S31desc", 2, 0, 4, vec(s, 5));
    CHECK_EQUAL(fzB(false), selfCheck(d, vec(s31, 9)).ok);
}

TEST(control_R1vsR2_allCases)
{
    std::vector<PeriodDesc> cs = constructibleCases(), k = c1kCases();
    cs.insert(cs.end(), k.begin(), k.end());
    const std::string fails = runCases(cs, [](const PeriodDesc & d) -> std::string
    {
        const index_t n = d.n();
        const gsMatrix<real_t> U = samplePoints(d);
        bool ok;
        const Basis r1 = makeR1(d), r2 = makeR2(d);
        const gsMatrix<real_t> V1 = denseEval(r1, U, n, n, 0, ok);
        const gsMatrix<real_t> V2 = denseEval(r2, U, n, n, 0, ok);
        Verdict v = spanEqual(V1, V2);
        if (!v.ok) return "span: " + v.why;
        v = partitionOfUnity(V1);
        if (!v.ok) return "PoU: " + v.why;
        v = seamContinuity(seamData(r1, d, n, n, false), d.p, d.mu());
        if (!v.ok) return "seam: " + v.why;
        return "";
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(control_R2seam_allCases)
{
    const std::string fails = runCases(constructibleCases(), [](const PeriodDesc & d) -> std::string
    {
        const Basis r2 = makeR2(d);
        const Verdict v = seamContinuity(seamData(r2, d, d.n(), d.n(), true), d.p, d.mu());
        return whyOf(v);
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(control_R1anchors_p1hats_C2uP1N4)
{
    const PeriodDesc d = caseC2u(1, 4);
    const gsMatrix<real_t> U = samplePoints(d);
    const real_t N = 4;
    gsMatrix<real_t> hats(4, U.cols());
    for (index_t j = 0; j != 4; ++j)
        for (index_t c = 0; c != U.cols(); ++c)
        {
            real_t v = 0;
            for (int k = -1; k <= 1; ++k)
                v = (std::max)(v, (std::max)(real_t(0),
                        1 - N * math::abs(U(0, c) - real_t(j) / N - real_t(k))));
            hats(j, c) = v;
        }
    bool ok;
    CHECK_EQUAL(std::string(""), whyOf(valuesClose(r1Dense(d, U, 0, ok), hats)));
}

TEST(control_R1anchors_p2uniform_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    gsMatrix<real_t> U(1, 3);
    U << 0.25, 0.125, 0.875;
    gsMatrix<real_t> E(4, 3);
    E << 0,    0.125, 0.75,
         0.5,  0.75,  0.125,
         0.5,  0.125, 0,
         0,    0,     0.125;
    bool ok;
    CHECK_EQUAL(std::string(""), whyOf(valuesClose(r1Dense(d, U, 0, ok), E)));
}

TEST(control_helpersC6_nonPeriodic)
{
    const Basis c6 = makeC6();
    std::vector<real_t> br;
    br.push_back(0); br.push_back(.25); br.push_back(.5); br.push_back(.75); br.push_back(1);
    gsGaussRule<real_t> gr(3);
    gsMatrix<real_t> U;
    gsVector<real_t> w;
    gr.mapToAll(br, U, w);
    bool ok;
    CHECK_EQUAL(std::string(""), whyOf(partitionOfUnity(denseEval(c6, U, 6, 0, 0, ok))));
    CHECK_EQUAL(std::string(""), whyOf(fdClose(fdMatrix(c6, U, 6, 0, 1),
                                               denseEval(c6, U, 6, 0, 1, ok), 1e-6, 1e-6)));
    CHECK_EQUAL(std::string(""), whyOf(fdClose(fdMatrix(c6, U, 6, 0, 2),
                                               denseEval(c6, U, 6, 0, 2, ok), 1e-4, 1e-4)));

    const gsMatrix<real_t> u = pointMat(0.1);
    gsMatrix<real_t> d0(1, 1), d1(1, 1);
    d0 << 0.36;
    d1 << -4.8;
    CHECK_EQUAL(std::string(""), whyOf(valuesClose(denseEval(c6, u, 6, 0, 0, ok).row(0), d0)));
    CHECK_EQUAL(std::string(""), whyOf(derivsClose(denseEval(c6, u, 6, 0, 1, ok).row(0), d1)));
}

TEST(control_fdR1_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    const Basis r1 = makeR1(d);
    const gsMatrix<real_t> U = gaussOnly(d);
    bool ok;
    CHECK_EQUAL(std::string(""), whyOf(fdClose(fdMatrix(r1, U, 4, 4, 1),
                                               denseEval(r1, U, 4, 4, 1, ok), 1e-6, 1e-6)));
    CHECK_EQUAL(std::string(""), whyOf(fdClose(fdMatrix(r1, U, 4, 4, 2),
                                               denseEval(r1, U, 4, 4, 2, ok), 1e-4, 1e-4)));
}

TEST(control_negSpan_dropKnot_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    const real_t s[] = {0, .25, .75};
    const PeriodDesc dd = mkDesc("dropped", 2, 0, 1, vec(s, 3));
    const gsMatrix<real_t> U = samplePoints(d);
    bool ok;
    const gsMatrix<real_t> VA = r1Dense(d, U, 0, ok);
    const gsMatrix<real_t> VB = denseEval(makeR2(dd), U, 3, 3, 0, ok);
    CHECK_EQUAL(fzB(false), spanEqual(VA, VB).ok);
}

TEST(control_negSeam_C6)
{
    const Basis c6 = makeC6();
    SeamData D;
    bool ok;
    for (int k = 0; k <= 2; ++k)
        D.push_back(std::make_pair(denseEval(c6, pointMat(0), 6, 0, k, ok),
                                   denseEval(c6, pointMat(1), 6, 0, k, ok)));
    CHECK_EQUAL(fzB(false), seamContinuity(D, 2, 1).ok);
}

TEST(control_negPoU_scaled_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    bool ok;
    gsMatrix<real_t> V = r1Dense(d, samplePoints(d), 0, ok);
    V.row(0) *= 1.01;
    CHECK_EQUAL(fzB(false), partitionOfUnity(V).ok);
}

TEST(control_negR1_shiftedKnot_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    std::vector<real_t> t = stretchedKnots(d);
    t[4] += 1.0 / 1024;
    const Basis shifted((gsKnotVector<real_t>(t, 2)));
    const gsMatrix<real_t> U = samplePoints(d);
    bool ok;
    CHECK_EQUAL(fzB(false), valuesClose(denseEval(shifted, U, 4, 4, 0, ok),
                                        r1Dense(d, U, 0, ok)).ok);
    CHECK_EQUAL(fzB(false), derivsClose(denseEval(shifted, U, 4, 4, 1, ok),
                                        r1Dense(d, U, 1, ok)).ok);
}

TEST(control_negFd_C6)
{
    const Basis c6 = makeC6();
    const gsMatrix<real_t> U = gaussOnly(mkDesc("c6", 2, 0, 1,
        std::vector<real_t>(caseC2u(2, 4).s)));
    bool ok;
    gsMatrix<real_t> deriv = denseEval(c6, U, 6, 0, 1, ok);
    deriv *= 1.01;
    CHECK_EQUAL(fzB(false), fdClose(fdMatrix(c6, U, 6, 0, 1), deriv, 1e-6, 1e-6).ok);
}

TEST(control_negUnchanged_C6)
{
    const Basis c6 = makeC6();
    const BasisSnapshot snap = snapshot(c6);
    CHECK(unchanged(snap, c6).ok);
    Basis copy = c6;
    copy.uniformRefine();
    CHECK_EQUAL(fzB(false), unchanged(snap, copy).ok);
}

TEST(control_negActive_dup)
{
    gsMatrix<index_t> dup(3, 1), out(3, 1);
    dup << 0, 0, 1;
    out << 0, 1, 4;
    CHECK_EQUAL(fzB(false), activeUniqueInRange(dup, 4).ok);
    CHECK_EQUAL(fzB(false), activeUniqueInRange(out, 4).ok);
}

TEST(control_negGuard_R1)
{
    const PeriodDesc d = caseC2u(2, 4);
    CHECK_EQUAL(fzB(false), guardPeriodic(makeR1(d), d, false));
}

// =====================================================================
// Constructor and conversion contract
// =====================================================================

TEST(space_NEW_ctorValues_C1cAll)
{
    const std::string fails = runCases(c1cCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b(0.0, 1.0, (unsigned)(d.n() - 1), (int)d.p, 1u, true);
        if (!guardPeriodic(b, d, false)) return "guard";
        return whyOf(basisVsR1(b, d));
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_ctorSeam_C1cAll)
{
    const std::string fails = runCases(c1cCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b(0.0, 1.0, (unsigned)(d.n() - 1), (int)d.p, 1u, true);
        if (!guardPeriodic(b, d, false)) return "guard";
        return whyOf(seamContinuity(seamData(b, d, d.n(), 0, false), d.p, d.mu()));
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_ctorSpanR2_C1cAll)
{
    const std::string fails = runCases(c1cCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b(0.0, 1.0, (unsigned)(d.n() - 1), (int)d.p, 1u, true);
        if (!guardPeriodic(b, d, false)) return "guard";
        const gsMatrix<real_t> U = samplePoints(d, knotsOf(b));
        bool ok;
        return whyOf(spanEqual(denseEval(b, U, d.n(), 0, 0, ok),
                               denseEval(makeR2(d), U, d.n(), d.n(), 0, ok)));
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_ctorMult2_C1cP3N4)
{
    const PeriodDesc d = caseC1c(3, 4, 2);
    Basis b(0.0, 1.0, 3u, 3, 2u, true);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    CHECK_EQUAL(std::string(""), whyOf(basisVsR1(b, d)));
}

TEST(space_NEW_ctorRefusedSmallN_C1cP4N4)
{
    CHECK_THROW((gsBSplineBasis<>(0.0, 1.0, 3, 4, 1, true)), std::runtime_error);
    CHECK(!fzB(false));
}

TEST(space_NEW_ctorRefusedMultP1_C1cP2N4)
{
    CHECK_THROW((gsBSplineBasis<>(0.0, 1.0, 3, 2, 3, true)), std::runtime_error);
    CHECK(!fzB(false));
}

TEST(space_B1_clampedConvValues_C1kAll)
{
    const std::string fails = runCases(c1kCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b = clampedConv(d);
        if (!guardPeriodic(b, d, false)) return "guard";
        return whyOf(basisVsR1(b, d));
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_B1_clampedConvDomain_C1kAll)
{
    const std::string fails = runCases(c1kCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b = clampedConv(d);
        if (!guardPeriodic(b, d, false)) return "guard";
        std::string w = whyOf(scalarClose(b.domainStart(), d.a, d.scale()));
        w += whyOf(scalarClose(b.domainEnd(), d.b, d.scale()));
        return w;
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_B1_clampedConvSpanR2_C1kAll)
{
    const std::string fails = runCases(c1kCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b = clampedConv(d);
        if (!guardPeriodic(b, d, false)) return "guard";
        const gsMatrix<real_t> U = samplePoints(d, knotsOf(b));
        bool ok;
        return whyOf(spanEqual(denseEval(b, U, d.n(), 0, 0, ok),
                               denseEval(makeR2(d), U, d.n(), d.n(), 0, ok)));
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_convValues_C2uAll)
{
    const std::string fails = runCases(c2uCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b = convert(d);
        if (!guardPeriodic(b, d, false)) return "guard";
        return whyOf(basisVsR1(b, d));
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_convValues_C2a)
{
    const PeriodDesc d = caseC2a();
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    CHECK_EQUAL(std::string(""), whyOf(basisVsR1(b, d)));
}

TEST(space_NEW_convValues_C3n)
{
    const PeriodDesc d = caseC3n();
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    CHECK_EQUAL(std::string(""), whyOf(basisVsR1(b, d)));
}

TEST(space_NEW_convValues_C4i)
{
    const PeriodDesc d = caseC4i();
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    CHECK_EQUAL(std::string(""), whyOf(basisVsR1(b, d)));
}

TEST(space_NEW_convValues_C4o)
{
    const PeriodDesc d = caseC4o();
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    CHECK_EQUAL(std::string(""), whyOf(basisVsR1(b, d)));
}

TEST(space_NEW_convValues_E2)
{
    const PeriodDesc d = caseE2();
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    CHECK_EQUAL(std::string(""), whyOf(basisVsR1(b, d)));
}

TEST(space_NEW_convCtorFlag_C3n)
{
    const PeriodDesc d = caseC3n();
    Basis b(stretchedKV(d), true);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    CHECK_EQUAL(std::string(""), whyOf(basisVsR1(b, d)));
}

TEST(space_NEW_convIntOverload_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b(stretchedKV(d));
    b.setPeriodic(0);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    CHECK_EQUAL(std::string(""), whyOf(basisVsR1(b, d)));
}

TEST(space_NEW_spanR2_allCases)
{
    const std::string fails = runCases(constructibleCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b = convert(d);
        if (!guardPeriodic(b, d, false)) return "guard";
        const gsMatrix<real_t> U = samplePoints(d, knotsOf(b));
        bool ok;
        return whyOf(spanEqual(denseEval(b, U, d.n(), 0, 0, ok),
                               denseEval(makeR2(d), U, d.n(), d.n(), 0, ok)));
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_partitionOfUnity_allCases)
{
    const std::string fails = runCases(constructibleCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b = convert(d);
        if (!guardPeriodic(b, d, false)) return "guard";
        bool ok;
        return whyOf(partitionOfUnity(denseEval(b, samplePoints(d), d.n(), 0, 0, ok)));
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_seamContinuity_allCases)
{
    const std::string fails = runCases(constructibleCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b = convert(d);
        if (!guardPeriodic(b, d, false)) return "guard";
        return whyOf(seamContinuity(seamData(b, d, d.n(), 0, false), d.p, d.mu()));
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_activeUnique_allCases)
{
    const std::string fails = runCases(constructibleCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b = convert(d);
        if (!guardPeriodic(b, d, false)) return "guard";
        gsMatrix<index_t> act;
        b.active_into(samplePoints(d), act);
        std::string w = whyOf(activeUniqueInRange(act, d.n()));
        if (b.numActive() != fzN(d.p + 1))
            w += " numActive " + num(b.numActive()) + " != " + num(d.p + 1);
        return w;
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_derivsR1_allCases)
{
    const std::string fails = runCases(constructibleCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b = convert(d);
        if (!guardPeriodic(b, d, false)) return "guard";
        std::string w;
        for (int k = 1; k <= (std::max)(2, (int)d.p); ++k)
        {
            const std::string v = whyOf(basisVsR1(b, d, k));
            if (!v.empty())
                w += "order " + num(k) + ": " + v + " ";
        }
        return w;
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_zeroOutside_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    gsMatrix<real_t> u(1, 2), vals;
    u << -0.1, 1.1;
    b.eval_into(u, vals);
    CHECK_EQUAL(fzB(true), vals.size() > 0);
    CHECK_EQUAL(std::string(""), whyOf(valuesClose(vals,
        gsMatrix<real_t>::Zero(vals.rows(), vals.cols()))));
}

TEST(space_B22_freeKnotAccepted_C3y)
{
    const PeriodDesc d = caseC3y();
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    CHECK_EQUAL(std::string(""), whyOf(basisVsR1(b, d)));
}

TEST(space_B22_firstSpanRefused_C3x)
{
    const real_t c3x[] = {-1.5, -.5, 0, .5, 1.5, 2, 3, 3.5, 4, 4.5, 5.5};
    Basis b(gsKnotVector<real_t>(vec(c3x, 11), 2));
    refuse(b, [](Basis & x) { x.setPeriodic(); });
}

TEST(space_B23_smallNRefused_C5a)
{
    Basis b(stretchedKV(caseC5a()));
    refuse(b, [](Basis & x) { x.setPeriodic(); });
}

TEST(space_B23_smallNRefused_C5b)
{
    Basis b(stretchedKV(caseC5b()));
    refuse(b, [](Basis & x) { x.setPeriodic(); });
}

TEST(space_B23_smallNRefused_C2uP4N4)
{
    Basis b(stretchedKV(caseC2u(4, 4)));
    refuse(b, [](Basis & x) { x.setPeriodic(); });
}

TEST(space_B31_asymSeamRefused_S31)
{
    const real_t s31[] = {-1, 0, 0, 1, 2, 3, 4, 5, 6};
    Basis b(gsKnotVector<real_t>(vec(s31, 9), 2));
    refuse(b, [](Basis & x) { x.setPeriodic(); });
}

TEST(space_B29_setPeriodicTwice_C2uAll)
{
    const std::string fails = runCases(c2uCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b = convert(d);
        b.setPeriodic();
        if (!guardPeriodic(b, d, false)) return "guard";
        return whyOf(basisVsR1(b, d));
    });
    CHECK_EQUAL(std::string(""), fails);
}

// =====================================================================
// Single-function routines
// =====================================================================

TEST(space_B4_evalSingle_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    const gsMatrix<real_t> U = samplePoints(d);
    bool ok;
    const gsMatrix<real_t> ref = r1Dense(d, U, 0, ok);
    std::string fails;
    for (index_t i = 0; i != 4; ++i)
    {
        gsMatrix<real_t> r;
        b.evalSingle_into(i, U, r);
        const Verdict v = valuesClose(r, ref.row(i));
        if (!v.ok) fails += "i=" + num(i) + ": " + v.why + "; ";
    }
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_B4_evalSingleTwoPoints_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    gsMatrix<real_t> u(1, 2), r;
    u << 0.875, 0.125;
    b.evalSingle_into(0, u, r);
    bool ok;
    const gsMatrix<real_t> ref = r1Dense(d, u, 0, ok);
    CHECK_EQUAL(std::string(""), whyOf(valuesClose(r, ref.row(0))));
}

TEST(space_B4_evalDerSingle_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    const gsMatrix<real_t> U = samplePoints(d);
    bool ok;
    const gsMatrix<real_t> ref = r1Dense(d, U, 1, ok);
    std::string fails;
    for (index_t i = 0; i != 4; ++i)
    {
        gsMatrix<real_t> r;
        b.evalDerSingle_into(i, U, 1, r);
        const Verdict v = derivsClose(r, ref.row(i));
        if (!v.ok) fails += "i=" + num(i) + ": " + v.why + "; ";
    }
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_B4_evalAllDersSingle_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    const gsMatrix<real_t> U = samplePoints(d);
    bool ok;
    std::string fails;
    for (index_t i = 0; i != 4; ++i)
    {
        gsMatrix<real_t> r;
        b.evalAllDersSingle_into(i, U, 2, r);
        for (int k = 0; k <= 2; ++k)
        {
            const gsMatrix<real_t> ref = r1Dense(d, U, k, ok);
            const Verdict v = 0 == k ? valuesClose(r.row(k), ref.row(i))
                                     : derivsClose(r.row(k), ref.row(i));
            if (!v.ok) fails += "i=" + num(i) + " k=" + num(k) + ": " + v.why + "; ";
        }
    }
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_B5_derivSingle_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    const gsMatrix<real_t> U = samplePoints(d);
    bool ok;
    const gsMatrix<real_t> ref = r1Dense(d, U, 1, ok);
    std::string fails;
    for (index_t i = 0; i != 4; ++i)
    {
        gsMatrix<real_t> r;
        b.derivSingle_into(i, U, r);
        const Verdict v = derivsClose(r, ref.row(i));
        if (!v.ok) fails += "i=" + num(i) + ": " + v.why + "; ";
    }
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_B5_deriv2Single_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    const gsMatrix<real_t> U = samplePoints(d);
    bool ok;
    const gsMatrix<real_t> ref = r1Dense(d, U, 2, ok);
    std::string fails;
    for (index_t i = 0; i != 4; ++i)
    {
        gsMatrix<real_t> r;
        b.deriv2Single_into(i, U, r);
        const Verdict v = derivsClose(r, ref.row(i));
        if (!v.ok) fails += "i=" + num(i) + ": " + v.why + "; ";
    }
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_B3_evalDerSingleNonPer_C6)
{
    const Basis c6 = makeC6();
    gsMatrix<real_t> u(1, 2), r, e(1, 2);
    u << 0.9, 0.1;
    e << 0, -4.8;
    c6.evalDerSingle_into(0, u, 1, r);
    CHECK_EQUAL(std::string(""), whyOf(derivsClose(r, e)));

    gsMatrix<real_t> u1(1, 1), r1, e1(1, 1);
    u1 << 1.0;
    e1 << 0;
    c6.evalDerSingle_into(0, u1, 1, r1);
    CHECK_EQUAL(std::string(""), whyOf(derivsClose(r1, e1)));
}

TEST(space_B3_evalAllDersSingleNonPer_C6)
{
    const Basis c6 = makeC6();
    gsMatrix<real_t> u(1, 3), r, e(2, 3);
    u << 0.9, 0.1, 1.0;
    e << 0, 0.36, 0,
         0, -4.8, 0;
    c6.evalAllDersSingle_into(0, u, 1, r);
    CHECK_EQUAL(std::string(""), whyOf(derivsClose(r, e)));
}

// =====================================================================
// Supports and boundary
// =====================================================================

TEST(space_B24_isActive_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    std::vector<real_t> us;
    for (int i = 0; i != 4; ++i)
        us.push_back((i + 0.5) / 4);
    for (int k = 4; k <= 8; ++k)
        us.push_back(1 - 0.25 / real_t(1 << k));
    std::string fails;
    for (size_t c = 0; c != us.size(); ++c)
    {
        bool ok;
        const gsMatrix<real_t> ref = r1Dense(d, pointMat(us[c]), 0, ok);
        gsVector<real_t> uv(1);
        uv << us[c];
        for (index_t j = 0; j != 4; ++j)
            if (b.isActive(j, uv) != fzB(ref(j, 0) > 1e-14))
                fails += "j=" + num(j) + " u=" + num(us[c]) + "; ";
    }
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_B25_supportHull_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    const std::vector<real_t> t = stretchedKnots(d);
    std::string fails;
    for (index_t j = 0; j != 4; ++j)
    {
        gsMatrix<real_t> e(1, 2);
        if (j < d.m())
            e << d.a, d.b;
        else
            e << (std::max)(d.a, t[j]), (std::min)(d.b, t[j + d.p + 1]);
        const Verdict v = valuesClose(b.support(j), e);
        if (!v.ok) fails += "j=" + num(j) + ": " + v.why + "; ";
    }
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_B25_elementSupport_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    const Basis r1 = makeR1(d);
    std::string fails;
    for (index_t j = 0; j != 4; ++j)
    {
        gsMatrix<index_t, 1, 2> es, ref;
        b.elementSupport_into(j, es);
        if (j < d.m())
            ref << 0, 4;
        else
            r1.elementSupport_into(j, ref);
        for (int c = 0; c != 2; ++c)
            if (es(0, c) != fzN(ref(0, c)))
                fails += "j=" + num(j) + " col " + num(c) + ": got " + num(es(0, c)) +
                         " expected " + num(ref(0, c)) + "; ";
    }
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_B25_elementInSupportOf_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    std::string fails;
    for (index_t j = 0; j != 4; ++j)
    {
        const gsMatrix<real_t> e = b.elementInSupportOf(j);
        const real_t mid = 0.5 * (e(0, 0) + e(0, 1));
        bool ok;
        const gsMatrix<real_t> ref = r1Dense(d, pointMat(mid), 0, ok);
        if (!fzB(ref(j, 0) > 1e-14))
            fails += "j=" + num(j) + " midpoint " + num(mid) + " outside the support; ";
    }
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_allBoundaryEmpty_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    const gsMatrix<index_t> ab = b.allBoundary();
    CHECK_EQUAL(fzN(0), (index_t)ab.size());
}

TEST(space_B8_boundaryOffsetRefused_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    refuse(b, [](Basis & x) { x.boundaryOffset(boundary::west, 0); });
    refuse(b, [](Basis & x) { x.boundaryOffset(boundary::east, 0); });
}

TEST(space_NEW_functionAtCornerRefused_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    refuse(b, [](Basis & x) { x.functionAtCorner(boxCorner(1)); });
    refuse(b, [](Basis & x) { x.functionAtCorner(boxCorner(2)); });
}

TEST(space_NEW_boundaryBasisRefused_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    refuse(b, [](Basis & x) { x.boundaryBasis(boundary::west); });
}

TEST(space_B28_matchWithRefused_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, false);
    CHECK(g);
    if (!g) return;
    const Basis c6 = makeC6();
    boundaryInterface bi(patchSide(0, boundary::east), patchSide(1, boundary::west),
                         static_cast<short_t>(1));
    refuse(b, [&](Basis & x)
    {
        gsMatrix<index_t> bt, bo;
        x.matchWith(bi, c6, bt, bo);
    });
}

// =====================================================================
// Refinement table
// =====================================================================

TEST(space_B2_uniformRefine_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    rowSuccess(d, uRef(1, 1), opUniformRefine(d, 1, 1));
}
TEST(info_twin_B2_uniformRefine_C2uP2N4)
{
    printTwin("space_B2_uniformRefine_C2uP2N4", caseC2u(2, 4), uRef(1, 1));
}

TEST(space_B2_uniformRefineMul2_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    rowSuccess(d, uRef(1, 2), opUniformRefine(d, 1, 2));
}
TEST(info_twin_B2_uniformRefineMul2_C2uP2N4)
{
    printTwin("space_B2_uniformRefineMul2_C2uP2N4", caseC2u(2, 4), uRef(1, 2));
}

TEST(space_B2_uniformRefine_C3n)
{
    const PeriodDesc d = caseC3n();
    rowSuccess(d, uRef(1, 1), opUniformRefine(d, 1, 1));
}
TEST(info_twin_B2_uniformRefine_C3n)
{
    printTwin("space_B2_uniformRefine_C3n", caseC3n(), uRef(1, 1));
}

TEST(space_B2_uniformRefine_C4i)
{
    const PeriodDesc d = caseC4i();
    rowSuccess(d, uRef(1, 1), opUniformRefine(d, 1, 1));
}
TEST(info_twin_B2_uniformRefine_C4i)
{
    printTwin("space_B2_uniformRefine_C4i", caseC4i(), uRef(1, 1));
}

TEST(space_NEW_insertKnotInside_C3n)
{
    const PeriodDesc d = caseC3n();
    rowSuccess(d, insK(2.5), opInsert(d, 2.5));
}
TEST(info_twin_NEW_insertKnotInside_C3n)
{
    printTwin("space_NEW_insertKnotInside_C3n", caseC3n(), insK(2.5));
}

TEST(space_NEW_insertKnotLeftWin_C3n)
{
    const PeriodDesc d = caseC3n();
    rowSuccess(d, insK(0.25), opInsert(d, 0.25));
}
TEST(info_twin_NEW_insertKnotLeftWin_C3n)
{
    printTwin("space_NEW_insertKnotLeftWin_C3n", caseC3n(), insK(0.25));
}

TEST(space_NEW_insertKnotRightWin_C3n)
{
    const PeriodDesc d = caseC3n();
    rowSuccess(d, insK(3.75), opInsert(d, 3.75));
}
TEST(info_twin_NEW_insertKnotRightWin_C3n)
{
    printTwin("space_NEW_insertKnotRightWin_C3n", caseC3n(), insK(3.75));
}

TEST(space_NEW_insertKnotOutside_C4o)
{
    const PeriodDesc d = caseC4o();
    rowSuccess(d, insK(3.5), opInsert(d, 3.5));
}
TEST(info_twin_NEW_insertKnotOutside_C4o)
{
    printTwin("space_NEW_insertKnotOutside_C4o", caseC4o(), insK(3.5));
}

TEST(space_NEW_insertKnotBelowA_C3n)
{
    const PeriodDesc d = caseC3n();
    rowSuccess(d, insK(-3.75), opInsert(d, 0.25));
}
TEST(info_twin_NEW_insertKnotBelowA_C3n)
{
    printTwin("space_NEW_insertKnotBelowA_C3n", caseC3n(), insK(-3.75));
}

TEST(space_NEW_insertKnotAboveB_C3n)
{
    const PeriodDesc d = caseC3n();
    rowSuccess(d, insK(4.25), opInsert(d, 0.25));
}
TEST(info_twin_NEW_insertKnotAboveB_C3n)
{
    printTwin("space_NEW_insertKnotAboveB_C3n", caseC3n(), insK(4.25));
}

TEST(space_NEW_insertKnotsVec_C3n)
{
    const PeriodDesc d = caseC3n();
    rowSuccess(d, insKnotsVec(0.25, 2.5), opInsert(opInsert(d, 0.25), 2.5));
}
TEST(info_twin_NEW_insertKnotsVec_C3n)
{
    printTwin("space_NEW_insertKnotsVec_C3n", caseC3n(), insKnotsVec(0.25, 2.5));
}

TEST(space_NEW_directKnotsEdit_C3n)
{
    const PeriodDesc d = caseC3n();
    rowSuccess(d, directInsert(0.25), opInsert(d, 0.25));
}
TEST(info_twin_NEW_directKnotsEdit_C3n)
{
    printTwin("space_NEW_directKnotsEdit_C3n", caseC3n(), directInsert(0.25));
}

TEST(space_NEW_insertKnotSeamRefused_C3n)
{
    const PeriodDesc d = caseC3n();
    rowRefusal(d, insK(0.0));
    rowRefusal(d, insK(4.0));
    rowRefusal(d, insK(-4.0));
    rowRefusal(d, insK(8.0));
}
TEST(info_twin_NEW_insertKnotSeamRefused_C3n)
{
    const real_t xs[] = {0.0, 4.0, -4.0, 8.0};
    for (int i = 0; i != 4; ++i)
        printTwin("space_NEW_insertKnotSeamRefused_C3n", caseC3n(), insK(xs[i]));
}

TEST(space_NEW_insertKnotFarRefused_C3n)
{
    const PeriodDesc d = caseC3n();
    rowRefusal(d, insK(-4.5));
    rowRefusal(d, insK(8.5));
}
TEST(info_twin_NEW_insertKnotFarRefused_C3n)
{
    printTwin("space_NEW_insertKnotFarRefused_C3n", caseC3n(), insK(-4.5));
    printTwin("space_NEW_insertKnotFarRefused_C3n", caseC3n(), insK(8.5));
}

TEST(space_NEW_insertKnotMultRefused_C3n)
{
    rowRefusal(caseC3n(), insK(2.0, 3));
}
TEST(info_twin_NEW_insertKnotMultRefused_C3n)
{
    printTwin("space_NEW_insertKnotMultRefused_C3n", caseC3n(), insK(2.0, 3));
}

TEST(space_NEW_refineElementsDomainIndex_C3n)
{
    const PeriodDesc d = caseC3n();
    rowSuccess(d, refineElem0(), opInsert(d, 0.25));
}
TEST(info_twin_NEW_refineElementsDomainIndex_C3n)
{
    printTwin("space_NEW_refineElementsDomainIndex_C3n", caseC3n(), refineElem0());
}

TEST(space_B2_removeKnot_C3n)
{
    const PeriodDesc d = caseC3n();
    rowSuccess(d, remK(2.0), opRemove(d, 2.0));
}
TEST(info_twin_B2_removeKnot_C3n)
{
    printTwin("space_B2_removeKnot_C3n", caseC3n(), remK(2.0));
}

TEST(space_B2_removeKnotWin_C3n)
{
    const PeriodDesc d = caseC3n();
    rowSuccess(d, remK(0.5), opRemove(d, 0.5));
}
TEST(info_twin_B2_removeKnotWin_C3n)
{
    printTwin("space_B2_removeKnotWin_C3n", caseC3n(), remK(0.5));
}

TEST(space_B2_removeKnotWin_C4i)
{
    const PeriodDesc d = caseC4i();
    rowSuccess(d, remK(1.0), opRemove(d, 1.0));
}
TEST(info_twin_B2_removeKnotWin_C4i)
{
    printTwin("space_B2_removeKnotWin_C4i", caseC4i(), remK(1.0));
}

TEST(space_B2_removeSeamRefused_C3n)
{
    rowRefusal(caseC3n(), remK(0.0));
}
TEST(info_twin_B2_removeSeamRefused_C3n)
{
    printTwin("space_B2_removeSeamRefused_C3n", caseC3n(), remK(0.0));
}

TEST(space_B2_removeSeamRefused_E2)
{
    rowRefusal(caseE2(), remK(0.0));
}
TEST(info_twin_B2_removeSeamRefused_E2)
{
    printTwin("space_B2_removeSeamRefused_E2", caseE2(), remK(0.0));
}

TEST(space_B2_removeSmallNRefused_C2uP3N4)
{
    rowRefusal(caseC2u(3, 4), remK(0.25));
}
TEST(info_twin_B2_removeSmallNRefused_C2uP3N4)
{
    printTwin("space_B2_removeSmallNRefused_C2uP3N4", caseC2u(3, 4), remK(0.25));
}

TEST(space_B2_degreeElevate_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    rowSuccess(d, dElev(), opDegreeElevate(d, 1));
}
TEST(info_twin_B2_degreeElevate_C2uP2N4)
{
    printTwin("space_B2_degreeElevate_C2uP2N4", caseC2u(2, 4), dElev());
}

TEST(space_B2_degreeElevate_C3n)
{
    const PeriodDesc d = caseC3n();
    rowSuccess(d, dElev(), opDegreeElevate(d, 1));
}
TEST(info_twin_B2_degreeElevate_C3n)
{
    printTwin("space_B2_degreeElevate_C3n", caseC3n(), dElev());
}

TEST(space_B2_degreeIncrease_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    rowSuccess(d, dInc(), opDegreeIncrease(d, 1));
}
TEST(info_twin_B2_degreeIncrease_C2uP2N4)
{
    printTwin("space_B2_degreeIncrease_C2uP2N4", caseC2u(2, 4), dInc());
}

TEST(space_B2_degreeIncRefused_C2uP3N4)
{
    rowRefusal(caseC2u(3, 4), dInc());
}
TEST(info_twin_B2_degreeIncRefused_C2uP3N4)
{
    printTwin("space_B2_degreeIncRefused_C2uP3N4", caseC2u(3, 4), dInc());
}

TEST(space_B2_degreeReduce_E2)
{
    const PeriodDesc d = caseE2();
    rowSuccess(d, dRed(), opDegreeReduce(d, 1));
}
TEST(info_twin_B2_degreeReduce_E2)
{
    printTwin("space_B2_degreeReduce_E2", caseE2(), dRed());
}

TEST(space_B2_degreeRedRefused_C2uP2N4)
{
    rowRefusal(caseC2u(2, 4), dRed());
}
TEST(info_twin_B2_degreeRedRefused_C2uP2N4)
{
    printTwin("space_B2_degreeRedRefused_C2uP2N4", caseC2u(2, 4), dRed());
}

TEST(space_B2_degreeDecrease_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    rowSuccess(d, dDec(), opDegreeDecrease(d, 1));
}
TEST(info_twin_B2_degreeDecrease_C2uP2N4)
{
    printTwin("space_B2_degreeDecrease_C2uP2N4", caseC2u(2, 4), dDec());
}

TEST(space_B2_degreeDecRefused_C2uP1N4)
{
    rowRefusal(caseC2u(1, 4), dDec());
}
TEST(info_twin_B2_degreeDecRefused_C2uP1N4)
{
    printTwin("space_B2_degreeDecRefused_C2uP1N4", caseC2u(1, 4), dDec());
}

TEST(space_B2_elevateContinuity_E2)
{
    const PeriodDesc d = caseE2();
    rowSuccess(d, elevCont(), opElevateContinuity(d, 1));
}
TEST(info_twin_B2_elevateContinuity_E2)
{
    printTwin("space_B2_elevateContinuity_E2", caseE2(), elevCont());
}

TEST(space_B2_elevContRefused_C2uP2N4)
{
    rowRefusal(caseC2u(2, 4), elevCont());
}
TEST(info_twin_B2_elevContRefused_C2uP2N4)
{
    printTwin("space_B2_elevContRefused_C2uP2N4", caseC2u(2, 4), elevCont());
}

TEST(space_B2_reduceContinuity_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    rowSuccess(d, redCont(), opReduceContinuity(d, 1));
}
TEST(info_twin_B2_reduceContinuity_C2uP2N4)
{
    printTwin("space_B2_reduceContinuity_C2uP2N4", caseC2u(2, 4), redCont());
}

TEST(space_B2_uniformCoarsen_C2uP2N8)
{
    const PeriodDesc d = caseC2u(2, 8);
    rowSuccess(d, coarsen1(), opCoarsen(d, 1));
}
TEST(info_twin_B2_uniformCoarsen_C2uP2N8)
{
    printTwin("space_B2_uniformCoarsen_C2uP2N8", caseC2u(2, 8), coarsen1());
}

TEST(space_B2_uniformCoarsenRefused_C2uP2N4)
{
    rowRefusal(caseC2u(2, 4), coarsen1());
}
TEST(info_twin_B2_uniformCoarsenRefused_C2uP2N4)
{
    printTwin("space_B2_uniformCoarsenRefused_C2uP2N4", caseC2u(2, 4), coarsen1());
}

TEST(space_B2_refineKRefused_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4), dOther = caseC2u(1, 4);
    const Basis other = convert(dOther);
    const bool go = guardPeriodic(other, dOther, false);
    CHECK(go);
    if (!go)
        return;
    rowRefusal(d, refineK(other));
}
TEST(info_twin_B2_refineKRefused_C2uP2N4)
{
    printTwin("space_B2_refineKRefused_C2uP2N4", caseC2u(2, 4), [](Basis & b)
    {
        const Basis other = makeR1(caseC2u(1, 4));
        b.refine_k(other, 1);
    });
}

TEST(space_NEW_reverseDomain_C3n)
{
    const PeriodDesc d = caseC3n();
    Basis b(stretchedKV(d));
    const std::string pre = preGuard(d, b);
    CHECK_EQUAL(std::string(""), pre);
    if (!pre.empty())
        return;
    try { b.reverse(); }
    catch (std::exception & e)
    {
        CHECK_EQUAL(std::string(""), std::string("operation threw: ") + e.what());
        return;
    }
    const std::string post = guardPeriodic(b, d, false) ? std::string()
                                                         : guardMsg("post-op guard", b, d);
    CHECK_EQUAL(std::string(""), post);
    if (!post.empty())
        return;
    CHECK_EQUAL(std::string(""), whyOf(scalarClose(b.domainStart(), 0, d.scale())));
    CHECK_EQUAL(std::string(""), whyOf(scalarClose(b.domainEnd(), 4, d.scale())));
}
TEST(info_twin_NEW_reverseDomain_C3n)
{
    printTwin("space_NEW_reverseDomain_C3n", caseC3n(), rev());
}

TEST(space_NEW_reverseValues_C3n)
{
    reverseValues(caseC3n());
}
TEST(info_twin_NEW_reverseValues_C3n)
{
    printTwin("space_NEW_reverseValues_C3n", caseC3n(), rev());
}

TEST(space_NEW_reverseValues_C2uP2N4)
{
    reverseValues(caseC2u(2, 4));
}
TEST(info_twin_NEW_reverseValues_C2uP2N4)
{
    printTwin("space_NEW_reverseValues_C2uP2N4", caseC2u(2, 4), rev());
}

// =====================================================================
// Probes (never fail)
// =====================================================================

TEST(info_probe_nonDyadic_P03)
{
    std::vector<real_t> s;
    for (int i = 0; i != 4; ++i)
        s.push_back(i * 0.075);
    probe("info_probe_nonDyadic_P03", mkDesc("P03", 2, 0, 0.3, s));
}

TEST(info_probe_tinyDomain_E9)
{
    std::vector<real_t> s;
    for (int i = 0; i != 4; ++i)
        s.push_back(i * 2.5e-10);
    probe("info_probe_tinyDomain_E9", mkDesc("E9", 2, 0, 1e-9, s));
}

TEST(info_probe_farPeriod_E6)
{
    std::vector<real_t> s;
    for (int i = 0; i != 4; ++i)
        s.push_back(1e6 + i / 4.0);
    probe("info_probe_farPeriod_E6", mkDesc("E6", 2, 1e6, 1e6 + 1, s));
}

}
