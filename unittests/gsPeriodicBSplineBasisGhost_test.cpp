/** @file gsPeriodicBSplineBasisGhost_test.cpp

    @brief Internals contract of one-dimensional periodic B-spline bases:
    stretched knots, number of crossing functions, trueSize, twin and the
    coefficient helpers.

    Expected knot vectors come from the stretched-knot generator of
    gsPeriodicTestUtils.h, which uses a period description only.

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

Basis convert(const PeriodDesc & d)
{
    Basis b(stretchedKV(d));
    b.setPeriodic();
    return b;
}

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

/// "" iff b is periodic with n functions and m crossing functions and its
/// knot vector is the stretched vector of d.
std::string ghostCheck(const Basis & b, const PeriodDesc & d)
{
    if (!guardPeriodic(b, d, true))
        return guardMsg("guard", b, d) + " crossing=" + num(b.numCrossingFunctions()) +
               " expected m=" + num(d.m());
    const std::vector<real_t> got = knotsOf(b), exp = stretchedKnots(d);
    if ((index_t)got.size() != fzN((index_t)exp.size()))
        return "knot count " + num(got.size()) + " vs " + num(exp.size()) +
               ": got " + vecStr(got);
    return whyOf(knotsClose(got, exp, d.scale()));
}

std::string ghostCheckLiteral(const Basis & b, const PeriodDesc & d, const real_t * lit, size_t k)
{
    if (!guardPeriodic(b, d, true))
        return guardMsg("guard", b, d) + " crossing=" + num(b.numCrossingFunctions()) +
               " expected m=" + num(d.m());
    return whyOf(knotsClose(knotsOf(b), vec(lit, k), d.scale()));
}

/// Convert d, apply op, and compare the stretched knots with those of dp.
void ghostRow(const PeriodDesc & d, const Op & op, const PeriodDesc & dp)
{
    Basis b(stretchedKV(d));
    std::string pre;
    try { b.setPeriodic(); }
    catch (std::exception & e)
    { pre = std::string("pre-op guard: conversion threw: ") + e.what(); }
    if (pre.empty() && !guardPeriodic(b, d, true))
        pre = guardMsg("pre-op guard", b, d);
    CHECK_EQUAL(std::string(""), pre);
    if (!pre.empty())
        return;
    try { op(b); }
    catch (std::exception & e)
    {
        CHECK_EQUAL(std::string(""), std::string("operation threw: ") + e.what());
        return;
    }
    CHECK_EQUAL(std::string(""), ghostCheck(b, dp));
}

Op uRef()      { return [](Basis & b) { b.uniformRefine(1, 1); }; }
Op dElev()     { return [](Basis & b) { b.degreeElevate((short_t)1); }; }
Op dInc()      { return [](Basis & b) { b.degreeIncrease((short_t)1); }; }
Op dRed()      { return [](Basis & b) { b.degreeReduce((short_t)1); }; }
Op dDec()      { return [](Basis & b) { b.degreeDecrease((short_t)1); }; }
Op elevCont()  { return [](Basis & b) { b.elevateContinuity(1); }; }
Op redCont()   { return [](Basis & b) { b.reduceContinuity(1); }; }
Op coarsen1()  { return [](Basis & b) { b.uniformCoarsen(1); }; }
Op remK(real_t x)   { return [=](Basis & b) { b.removeKnot(x); }; }
Op insK(real_t x)   { return [=](Basis & b) { b.insertKnot(x); }; }
Op directInsert(real_t x) { return [=](Basis & b) { b.knots().insert(x); }; }
Op rev()       { return [](Basis & b) { b.reverse(); }; }

} // namespace

SUITE(gsPeriodicBSplineBasisGhost_test)
{

TEST(ghost_NEW_stretchedKnots_C2uAll)
{
    const std::string fails = runCases(c2uCases(), [](const PeriodDesc & d) -> std::string
    {
        return ghostCheck(convert(d), d);
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(ghost_NEW_stretchedKnots_C2a)
{
    CHECK_EQUAL(std::string(""), ghostCheck(convert(caseC2a()), caseC2a()));
}

TEST(ghost_NEW_stretchedKnots_C3n)
{
    CHECK_EQUAL(std::string(""), ghostCheck(convert(caseC3n()), caseC3n()));
}

TEST(ghost_NEW_stretchedKnots_C3y)
{
    CHECK_EQUAL(std::string(""), ghostCheck(convert(caseC3y()), caseC3y()));
}

TEST(ghost_NEW_stretchedKnots_C4i)
{
    CHECK_EQUAL(std::string(""), ghostCheck(convert(caseC4i()), caseC4i()));
}

TEST(ghost_NEW_stretchedKnots_C4o)
{
    CHECK_EQUAL(std::string(""), ghostCheck(convert(caseC4o()), caseC4o()));
}

TEST(ghost_NEW_stretchedKnots_E2)
{
    CHECK_EQUAL(std::string(""), ghostCheck(convert(caseE2()), caseE2()));
}

TEST(ghost_B1_clampedConvKnots_C1kAll)
{
    const std::string fails = runCases(c1kCases(), [](const PeriodDesc & d) -> std::string
    {
        return ghostCheck(clampedConv(d), d);
    });
    CHECK_EQUAL(std::string(""), fails);

    const real_t l4p2[] = {-1, 0, 0, 1, 2, 3, 4, 4, 5};
    CHECK_EQUAL(std::string(""), ghostCheckLiteral(clampedConv(caseC1kL4(2)), caseC1kL4(2), l4p2, 9));
}

TEST(ghost_NEW_ctorKnots_C1cAll)
{
    const std::string fails = runCases(c1cCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b(0.0, 1.0, (unsigned)(d.n() - 1), (int)d.p, 1u, true);
        return ghostCheck(b, d);
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(ghost_NEW_trueSize_allCases)
{
    const std::string fails = runCases(constructibleCases(), [](const PeriodDesc & d) -> std::string
    {
        Basis b = convert(d);
        if (!guardPeriodic(b, d, true)) return "guard";
        if (b.trueSize() != fzN(d.n() + d.m()))
            return "trueSize " + num(b.trueSize()) + " != " + num(d.n() + d.m());
        return "";
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(ghost_B7_twinAll_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, true);
    CHECK(g);
    if (!g) return;
    const index_t n = d.n(), m = d.m();
    for (index_t i = 0; i != n + m; ++i)
    {
        const index_t exp = i < m ? i + n : (i >= n ? i - n : i);
        CHECK_EQUAL(fzN(exp), b.twin(i));
    }
    CHECK_EQUAL(fzN(0), b.twin(b.size()));
}

TEST(ghost_NEW_perCoefs_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, true);
    CHECK(g);
    if (!g) return;
    const index_t n = d.n(), m = d.m();
    gsMatrix<real_t> C(n + m, 2), E(n + m, 2);
    for (index_t r = 0; r != n + m; ++r)
    {
        C(r, 0) = r + 1;
        C(r, 1) = 10 * (r + 1);
        const index_t src = r < n ? r : r - n;
        E(r, 0) = src + 1;
        E(r, 1) = 10 * (src + 1);
    }
    CHECK_EQUAL(std::string(""), whyOf(valuesClose(b.perCoefs(C), E)));
}

TEST(ghost_NEW_expandCoefs_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, true);
    CHECK(g);
    if (!g) return;
    const index_t n = d.n(), m = d.m();
    gsMatrix<real_t> C(n, 2), E(n + m, 2);
    for (index_t r = 0; r != n; ++r)
    {
        C(r, 0) = r + 1;
        C(r, 1) = 10 * (r + 1);
    }
    for (index_t r = 0; r != n + m; ++r)
    {
        const index_t src = r < n ? r : r - n;
        E(r, 0) = src + 1;
        E(r, 1) = 10 * (src + 1);
    }
    b.expandCoefs(C);
    CHECK_EQUAL(fzN(n + m), (index_t)C.rows());
    CHECK_EQUAL(std::string(""), whyOf(valuesClose(C, E)));
}

TEST(ghost_NEW_trimCoefs_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, true);
    CHECK(g);
    if (!g) return;
    const index_t n = d.n(), m = d.m();
    gsMatrix<real_t> C(n + m, 2), E(n, 2);
    for (index_t r = 0; r != n + m; ++r)
    {
        C(r, 0) = r + 1;
        C(r, 1) = 10 * (r + 1);
        if (r < n)
        {
            E(r, 0) = r + 1;
            E(r, 1) = 10 * (r + 1);
        }
    }
    b.trimCoefs(C);
    CHECK_EQUAL(fzN(n), (index_t)C.rows());
    CHECK_EQUAL(std::string(""), whyOf(valuesClose(C, E)));
}

TEST(ghost_B2_uniformRefineKnots_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    ghostRow(d, uRef(), opUniformRefine(d, 1, 1));
}

TEST(ghost_B2_uniformRefineKnots_C4i)
{
    const PeriodDesc d = caseC4i();
    ghostRow(d, uRef(), opUniformRefine(d, 1, 1));
}

TEST(ghost_B2_degreeElevateKnots_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    ghostRow(d, dElev(), opDegreeElevate(d, 1));
}

TEST(ghost_B2_degreeIncreaseKnots_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    ghostRow(d, dInc(), opDegreeIncrease(d, 1));
}

TEST(ghost_B2_degreeReduceKnots_E2)
{
    const PeriodDesc d = caseE2();
    ghostRow(d, dRed(), opDegreeReduce(d, 1));
}

TEST(ghost_B2_degreeDecreaseKnots_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    ghostRow(d, dDec(), opDegreeDecrease(d, 1));
}

TEST(ghost_B2_elevContKnots_E2)
{
    const PeriodDesc d = caseE2();
    ghostRow(d, elevCont(), opElevateContinuity(d, 1));
}

TEST(ghost_B2_reduceContKnots_C2uP2N4)
{
    const PeriodDesc d = caseC2u(2, 4);
    ghostRow(d, redCont(), opReduceContinuity(d, 1));
}

TEST(ghost_B2_uniformCoarsenKnots_C2uP2N8)
{
    const PeriodDesc d = caseC2u(2, 8);
    ghostRow(d, coarsen1(), opCoarsen(d, 1));
}

TEST(ghost_B2_removeKnotKnots_C3n)
{
    const PeriodDesc d = caseC3n();
    ghostRow(d, remK(2.0), opRemove(d, 2.0));
}

TEST(ghost_NEW_insertKnotKnots_C3n)
{
    const PeriodDesc d = caseC3n();
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, true);
    CHECK(g);
    if (!g) return;
    insK(0.25)(b);
    const real_t exp[] = {-1, -.5, 0, .25, .5, 1.5, 2, 3, 3.5, 4, 4.25, 4.5};
    CHECK_EQUAL(std::string(""), ghostCheckLiteral(b, opInsert(d, 0.25), exp, 12));
}

TEST(ghost_NEW_directKnotsEditI1_C3n)
{
    const PeriodDesc d = caseC3n();
    Basis b = convert(d);
    const bool g = guardPeriodic(b, d, true);
    CHECK(g);
    if (!g) return;
    directInsert(0.25)(b);
    const real_t exp[] = {-1, -.5, 0, .25, .5, 1.5, 2, 3, 3.5, 4, 4.25, 4.5};
    CHECK_EQUAL(std::string(""), ghostCheckLiteral(b, opInsert(d, 0.25), exp, 12));
}

TEST(ghost_NEW_reverseKnots_C3n)
{
    const PeriodDesc d = caseC3n();
    ghostRow(d, rev(), opReverse(d));
}

}
