/** @file gsPeriodicKnotVector_test.cpp

    @brief Public behaviour of periodic one-dimensional knot vectors.

    A periodic gsKnotVector stores the stretched knots of a period
    description (degree, domain [a,b], breakpoints of [a,b) with
    multiplicity).  Expected values come from such descriptions (see
    gsPeriodicTestUtils.h), never from the knot vector under test; the only
    exceptions are the unchanged/carry checks, which compare an object with
    its own earlier snapshot.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
**/

#include "gismo_unittest.h"
#include "gsPeriodicTestUtils.h"

#include <cmath>
#include <iterator>
#include <memory>

using namespace gismo::periodic_test;

namespace {

typedef gsKnotVector<real_t>               KV;
typedef std::function<void(KV &)>          KOp;

// --- state helpers -------------------------------------------------------

std::vector<real_t> kvKnots(const KV & kv)
{
    return std::vector<real_t>(kv.begin(), kv.end());
}

/// Observable state of a knot vector, for "unchanged" checks.
struct KvSnap
{
    std::vector<real_t> knots;
    size_t              size;
    short_t             deg;
    bool                per;
    index_t             m, n;
};

KvSnap kvSnap(const KV & kv)
{
    KvSnap s;
    s.knots = kvKnots(kv);
    s.size  = kv.size();
    s.deg   = kv.degree();
    s.per   = kv.isPeriodic();
    s.m     = kv.numCrossingFunctions();
    s.n     = s.per ? kv.periodicSize() : 0;
    return s;
}

/// Exact comparison of a vector with its earlier snapshot.
Verdict kvUnchanged(const KvSnap & before0, const KV & after)
{
    KvSnap before = before0;
    if (falsifying() && !before.knots.empty())
        before.knots[0] += 1e-6 * (std::max)(real_t(1), math::abs(before.knots[0]));
    const std::vector<real_t> k = kvKnots(after);
    if (k != before.knots)
        return mkBad(0, "knots changed: " + vecStr(before.knots) + " -> " + vecStr(k));
    if (after.size() != before.size)
        return mkBad(0, "size changed: " + num(before.size) + " -> " + num(after.size()));
    if (after.degree() != before.deg)
        return mkBad(0, "degree changed");
    if (after.isPeriodic() != before.per)
        return mkBad(0, "periodic flag changed");
    if (after.numCrossingFunctions() != before.m)
        return mkBad(0, "numCrossingFunctions changed");
    if (after.isPeriodic() && after.periodicSize() != before.n)
        return mkBad(0, "periodicSize changed");
    return mkOk();
}

/// Prefixes a non-empty failure text with the case that produced it.
std::string tagged(const std::string & who, const std::string & w)
{
    return w.empty() ? w : who + ": " + w;
}

/// Refusal of op on kv: throws std::runtime_error and leaves kv unchanged.
void kvRefuse(KV & kv, const KOp & op, const std::string & who = std::string())
{
    const KvSnap s0 = kvSnap(kv);
    std::string failure;
    try { op(kv); failure = "did not throw"; }
    catch (std::runtime_error &) {}
    catch (std::exception & e) { failure = std::string("threw a different type: ") + e.what(); }
    CHECK_EQUAL(std::string(""), tagged(who, failure));
    CHECK_EQUAL(std::string(""), tagged(who, whyOf(kvUnchanged(s0, kv))));
}

/// "" if f runs without throwing, else the exception text.
std::string thrown(const std::function<void()> & f)
{
    try { f(); }
    catch (std::exception & e) { return std::string("threw: ") + e.what(); }
    return std::string();
}

KV kvConvert(const PeriodDesc & d)
{
    KV kv = stretchedKV(d);
    kv.setPeriodic(true);
    return kv;
}

/// With fz false the expectations are not perturbed under falsification.
bool fB(bool fz, bool v)        { return fz ? fzB(v) : v; }
index_t fN(bool fz, index_t v)  { return fz ? fzN(v) : v; }

std::string plainKnots(const std::vector<real_t> & got, const std::vector<real_t> & e, real_t scale)
{
    if (got.size() != e.size())
        return "knot count " + num(got.size()) + " vs " + num(e.size()) +
               ": got " + vecStr(got) + " expected " + vecStr(e);
    real_t err = 0;
    for (size_t i = 0; i != got.size(); ++i)
        err = (std::max)(err, math::abs(got[i] - e[i]));
    return err <= 1e-14 * scale ? std::string()
        : "knots differ by " + num(err) + ": got " + vecStr(got) + " expected " + vecStr(e);
}

/// "" iff kv is the periodic knot vector of dp: flag, m, n, L, degree,
/// size, stretched knots, a and (bit-exactly) b.
std::string kvState(const KV & kv, const PeriodDesc & dp, bool fz = true)
{
    std::string w;
    if (!fB(fz, kv.isPeriodic()))
        return "not periodic";
    if (kv.numCrossingFunctions() != fN(fz, dp.m()))
        w += "m " + num(kv.numCrossingFunctions()) + " != " + num(dp.m()) + "; ";
    if (kv.periodicSize() != fN(fz, dp.n()))
        w += "n " + num(kv.periodicSize()) + " != " + num(dp.n()) + "; ";
    if (fz)
        w += whyOf(scalarClose(kv.period(), dp.L(), dp.scale()));
    else if (!(math::abs(kv.period() - dp.L()) <= 1e-14 * dp.scale()))
        w += "period " + num(kv.period()) + " != " + num(dp.L());
    if (kv.degree() != dp.p)
        w += "degree " + num(kv.degree()) + " != " + num(dp.p) + "; ";
    if ((index_t)kv.size() != dp.n() + dp.m() + dp.p + 1)
        w += "size " + num(kv.size()) + " != " + num(dp.n() + dp.m() + dp.p + 1) + "; ";
    w += fz ? whyOf(knotsClose(kvKnots(kv), stretchedKnots(dp), dp.scale()))
            : plainKnots(kvKnots(kv), stretchedKnots(dp), dp.scale());
    if (!fB(fz, *kv.domainBegin() == dp.a))
        w += "a " + num(*kv.domainBegin()) + " != " + num(dp.a) + "; ";
    if (!fB(fz, *kv.domainEnd() == dp.b))
        w += "b " + num(*kv.domainEnd()) + " != " + num(dp.b) + " (exact); ";
    return w;
}

// --- cases ---------------------------------------------------------------

PeriodDesc caseByTag(const std::string & t)
{
    if ("C2uP1N4" == t) return caseC2u(1, 4);
    if ("C2uP2N4" == t) return caseC2u(2, 4);
    if ("C2uP2N8" == t) return caseC2u(2, 8);
    if ("C2uP3N4" == t) return caseC2u(3, 4);
    if ("C2uP4N4" == t) return caseC2u(4, 4);
    if ("C3n" == t)     return caseC3n();
    if ("C3y" == t)     return caseC3y();
    if ("C4i" == t)     return caseC4i();
    if ("C4o" == t)     return caseC4o();
    if ("E2" == t)      return caseE2();
    GISMO_ERROR("unknown case tag " << t);
}

std::vector<PeriodDesc> fourCases()
{
    std::vector<PeriodDesc> c;
    c.push_back(caseC2u(2, 4));
    c.push_back(caseC3n());
    c.push_back(caseC4i());
    c.push_back(caseC4o());
    return c;
}

std::vector<PeriodDesc> fiveCases()
{
    std::vector<PeriodDesc> c = fourCases();
    c.push_back(caseE2());
    return c;
}

/// Knots of the literal vectors that are not period descriptions.
std::vector<real_t> knotsC3x()
{
    const real_t v[] = {-1.5, -.5, 0, .5, 1.5, 2, 3, 3.5, 4, 4.5, 5.5};
    return vec(v, 11);
}

std::vector<real_t> knotsS31()
{
    const real_t v[] = {-1, 0, 0, 1, 2, 3, 4, 5, 6};
    return vec(v, 9);
}

PeriodDesc caseE6()
{
    std::vector<real_t> s;
    for (int i = 0; i != 4; ++i)
        s.push_back(1e6 + real_t(i) / 4);
    return mkDesc("E6", 2, 1e6, 1e6 + 1, s);
}

/// d with all breakpoints and the domain shifted by amount.
PeriodDesc shifted(PeriodDesc d, real_t amount)
{
    d.a += amount;
    d.b += amount;
    for (size_t i = 0; i != d.s.size(); ++i)
        d.s[i] += amount;
    return d;
}

/// d with every breakpoint above a limited to multiplicity cap.
PeriodDesc capMultiplicity(PeriodDesc d, index_t cap)
{
    std::vector<real_t> u = uniqueSorted(d.s), s;
    for (size_t j = 0; j != u.size(); ++j)
    {
        const index_t mult = countEq(d.s, u[j]);
        pushMult(s, u[j], u[j] == d.a ? mult : (std::min)(mult, cap));
    }
    d.s = s;
    return d;
}

PeriodDesc reducedInsert(PeriodDesc d, real_t x)
{
    return opInsert(d, x < d.a ? x + d.L() : x - d.L());
}

// --- scenarios: operation, case and expected description -----------------

struct Scen
{
    PeriodDesc d;
    std::string label;
    KOp op;
    PeriodDesc dp;
    std::shared_ptr<std::vector<real_t> > got; ///< coarsen: returned knots
    std::vector<real_t> removed;               ///< coarsen: expected returned knots
    bool hasRet;
};

Scen mkScen(const PeriodDesc & d, const std::string & label, const KOp & op,
            const PeriodDesc & dp)
{
    Scen s;
    s.d = d; s.label = label; s.op = op; s.dp = dp; s.hasRet = false;
    return s;
}

std::string labelOf(real_t x) { return "x=" + num(x); }

/// Row driver: converts d, applies op; the result is the periodic vector
/// of dp when constructible(dp) and a refusal (unchanged vector) otherwise.
/// Returns true when the operation was expected to be accepted.
bool runScen(const Scen & sc)
{
    const std::string who = sc.d.tag + (sc.label.empty() ? "" : " " + sc.label);
    KV kv;
    try { kv = kvConvert(sc.d); }
    catch (std::exception & e)
    {
        CHECK_EQUAL(std::string(""), who + ": pre-op guard: conversion threw: " + e.what());
        return false;
    }
    const std::string pre = kvState(kv, sc.d, false);
    CHECK_EQUAL(std::string(""), tagged(who, pre.empty() ? pre : "pre-op guard: " + pre));
    if (!pre.empty())
        return false;
    if (!constructible(sc.dp))
    {
        kvRefuse(kv, sc.op, who);
        return false;
    }
    try { sc.op(kv); }
    catch (std::exception & e)
    {
        CHECK_EQUAL(std::string(""), who + ": operation threw: " + e.what());
        return true;
    }
    CHECK_EQUAL(std::string(""), tagged(who, kvState(kv, sc.dp)));
    if (sc.hasRet)
        CHECK_EQUAL(std::string(""),
                    tagged(who, whyOf(knotsClose(*sc.got, sc.removed, sc.d.scale()))));
    return true;
}

/// Number of scenarios expected to be accepted.
int runAll(const std::vector<Scen> & v)
{
    int accepted = 0;
    for (size_t i = 0; i != v.size(); ++i)
        if (runScen(v[i]))
            ++accepted;
    return accepted;
}

/// Prints what the flag-off twin (same stretched knots) does with op.
void printKvTwin(const std::string & name, const PeriodDesc & d, const KOp & op)
{
    try
    {
        KV kv = stretchedKV(d);
        bool threw = false;
        std::string what;
        try { op(kv); }
        catch (std::exception & e) { threw = true; what = e.what(); }
        gsInfo << "info_twin of " << name << " on " << d.tag << ": threw=" << (threw ? 1 : 0);
        if (threw)
            gsInfo << " (" << what << ")";
        gsInfo << " size=" << kv.size() << " knots=" << vecStr(kvKnots(kv)) << "\n";
    }
    catch (std::exception & e)
    {
        gsInfo << "info_twin of " << name << ": exception: " << e.what() << "\n";
    }
}

void twinAll(const char * name, const std::vector<Scen> & v)
{
    try
    {
        for (size_t i = 0; i != v.size(); ++i)
            printKvTwin(std::string(name) + " " + v[i].label, v[i].d, v[i].op);
    }
    catch (std::exception & e)
    {
        gsInfo << "info_twin of " << name << ": exception: " << e.what() << "\n";
    }
}

// --- scenario lists (shared by the tests and their twins) ----------------

struct TagX { const char * tag; real_t x; };

KOp insertOp(real_t x, index_t mult = 1)
{
    return [=](KV & k) { k.insert(x, mult); };
}

std::vector<Scen> scInsertInside()
{
    const TagX t[] = {{"C2uP2N4", .125}, {"C2uP2N4", .625}, {"C3n", .25}, {"C4i", .5},
                      {"C4o", 1.5}, {"C3n", 2.5}, {"C4i", 1.5}, {"C4o", 2.5}};
    std::vector<Scen> v;
    for (size_t i = 0; i != 8; ++i)
    {
        const PeriodDesc d = caseByTag(t[i].tag);
        v.push_back(mkScen(d, labelOf(t[i].x), insertOp(t[i].x), opInsert(d, t[i].x)));
    }
    return v;
}

std::vector<Scen> scInsertMult()
{
    const TagX t[] = {{"C2uP2N4", .375}, {"C3n", 2.5}, {"C4i", 2.5}, {"C4o", 2.5}};
    std::vector<Scen> v;
    for (size_t i = 0; i != 4; ++i)
    {
        const PeriodDesc d = caseByTag(t[i].tag);
        v.push_back(mkScen(d, labelOf(t[i].x), insertOp(t[i].x, 2),
                           opInsert(opInsert(d, t[i].x), t[i].x)));
    }
    return v;
}

std::vector<Scen> scInsertWrap()
{
    const TagX t[] = {{"C2uP2N4", -.875}, {"C2uP2N4", 1.125}, {"C3n", -3.75}, {"C3n", 4.25},
                      {"C4i", -3.5}, {"C4i", 4.5}, {"C4o", -4.5}, {"C4o", 7.5}};
    std::vector<Scen> v;
    for (size_t i = 0; i != 8; ++i)
    {
        const PeriodDesc d = caseByTag(t[i].tag);
        v.push_back(mkScen(d, labelOf(t[i].x), insertOp(t[i].x), reducedInsert(d, t[i].x)));
    }
    return v;
}

/// Sorted batch {window, outside, wrap} per case.
std::vector<std::vector<real_t> > rangeBatches()
{
    const real_t b0[] = {.125, .625, 1.125}, b1[] = {-3.75, .25, 2.5},
                 b2[] = {.5, 1.5, 4.5},      b3[] = {-4.5, 1.5, 2.5};
    std::vector<std::vector<real_t> > b;
    b.push_back(vec(b0, 3)); b.push_back(vec(b1, 3));
    b.push_back(vec(b2, 3)); b.push_back(vec(b3, 3));
    return b;
}

PeriodDesc insertAll(PeriodDesc d, const std::vector<real_t> & v)
{
    const PeriodDesc d0 = d;
    for (size_t i = 0; i != v.size(); ++i)
        d = (v[i] < d0.a || v[i] > d0.b) ? reducedInsert(d, v[i]) : opInsert(d, v[i]);
    return d;
}

std::vector<Scen> scInsertRange()
{
    const std::vector<PeriodDesc> c = fourCases();
    const std::vector<std::vector<real_t> > b = rangeBatches();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
    {
        const std::vector<real_t> x = b[i];
        v.push_back(mkScen(c[i], vecStr(x),
                           [=](KV & k) { k.insert(x.begin(), x.end()); },
                           insertAll(c[i], x)));
    }
    return v;
}

std::vector<Scen> scInsertVec()
{
    const std::vector<PeriodDesc> c = fourCases();
    const std::vector<std::vector<real_t> > b = rangeBatches();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
    {
        const std::vector<real_t> x = b[i];
        v.push_back(mkScen(c[i], vecStr(x),
                           [=](KV & k) { k.insert(x, 1); },
                           insertAll(c[i], x)));
    }
    return v;
}

const TagX removeValues[] = {{"C2uP2N4", .5}, {"C3n", .5}, {"C3n", 2}, {"C4i", 1},
                             {"C4i", 2}, {"C4o", 3}, {"C4o", 1}};

std::vector<Scen> scRemoveValue()
{
    std::vector<Scen> v;
    for (size_t i = 0; i != 7; ++i)
    {
        const PeriodDesc d = caseByTag(removeValues[i].tag);
        const real_t x = removeValues[i].x;
        v.push_back(mkScen(d, labelOf(x), [=](KV & k) { k.remove(x, 1); }, opRemove(d, x)));
    }
    return v;
}

std::vector<Scen> scRemoveWrap()
{
    const PeriodDesc d = caseC3n();
    std::vector<Scen> v;
    v.push_back(mkScen(d, "x=-0.5", [](KV & k) { k.remove(real_t(-0.5), 1); },
                       opRemove(d, 3.5)));
    return v;
}

std::vector<Scen> scRemoveIter()
{
    std::vector<Scen> v;
    for (size_t i = 0; i != 7; ++i)
    {
        const PeriodDesc d = caseByTag(removeValues[i].tag);
        const real_t x = removeValues[i].x;
        v.push_back(mkScen(d, labelOf(x),
                           [=](KV & k) { k.remove(std::lower_bound(k.ubegin(), k.uend(), x), 1); },
                           opRemove(d, x)));
    }
    return v;
}

std::vector<Scen> scUniformRefine(index_t mul)
{
    const std::vector<PeriodDesc> c = fourCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
        v.push_back(mkScen(c[i], "", [=](KV & k) { k.uniformRefine(1, mul); },
                           opUniformRefine(c[i], 1, mul)));
    return v;
}

struct TagE { const char * tag; index_t e0, e1; real_t x0, x1; };

const TagE spanData[] = {{"C2uP2N4", 0, 3, .125, .875}, {"C3n", 0, 5, .25, 3.75},
                         {"C4i", 0, 3, .5, 3.5},        {"C4o", 0, 5, .5, 5.5}};

std::vector<Scen> scRefineSpansMult()
{
    std::vector<Scen> v;
    for (size_t i = 0; i != 4; ++i)
    {
        const TagE e = spanData[i];
        const PeriodDesc d = caseByTag(e.tag);
        KV::multContainer idx;
        idx.push_back(e.e0); idx.push_back(e.e1);
        v.push_back(mkScen(d, vecStr(std::vector<real_t>(1, real_t(e.e0))) + vecStr(std::vector<real_t>(1, real_t(e.e1))),
                           [=](KV & k) { k.refineSpans(idx, 1); },
                           opInsert(opInsert(d, e.x0), e.x1)));
    }
    return v;
}

std::vector<Scen> scRefineSpansUnsigned()
{
    std::vector<Scen> v;
    for (size_t i = 0; i != 4; ++i)
    {
        const TagE e = spanData[i];
        const PeriodDesc d = caseByTag(e.tag);
        std::vector<unsigned> idx;
        idx.push_back((unsigned)e.e0); idx.push_back((unsigned)e.e1);
        v.push_back(mkScen(d, vecStr(std::vector<real_t>(1, real_t(e.e0))) + vecStr(std::vector<real_t>(1, real_t(e.e1))),
                           [=](KV & k) { k.refineSpans(idx, 1); },
                           opInsert(opInsert(d, e.x0), e.x1)));
    }
    return v;
}

std::vector<Scen> scRefineSpansIter()
{
    const PeriodDesc d = caseC3n();
    std::vector<Scen> v;
    v.push_back(mkScen(d, "e={0},k=3",
                       [](KV & k)
                       {
                           std::vector<index_t> e(1, 0);
                           k.refineSpans(e.begin(), e.end(), 3);
                       },
                       opInsert(opInsert(opInsert(d, .125), .25), .375)));
    return v;
}

std::vector<Scen> scDegreeElevate()
{
    const std::vector<PeriodDesc> c = fourCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
        v.push_back(mkScen(c[i], "", [](KV & k) { k.degreeElevate((short_t)1); },
                           opDegreeElevate(c[i], 1)));
    return v;
}

std::vector<Scen> scDegreeIncrease()
{
    const std::vector<PeriodDesc> c = fourCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
        v.push_back(mkScen(c[i], "", [](KV & k) { k.degreeIncrease(1); },
                           opDegreeIncrease(c[i], 1)));
    return v;
}

std::vector<Scen> scDegreeReduce()
{
    const std::vector<PeriodDesc> c = fiveCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
        v.push_back(mkScen(c[i], "", [](KV & k) { k.degreeReduce((short_t)1); },
                           opDegreeReduce(c[i], 1)));
    return v;
}

std::vector<Scen> scDegreeDecrease()
{
    const std::vector<PeriodDesc> c = fourCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
        v.push_back(mkScen(c[i], "", [](KV & k) { k.degreeDecrease(1); },
                           opDegreeDecrease(c[i], 1)));
    return v;
}

std::vector<Scen> scDegreeDecInterior()
{
    const std::vector<PeriodDesc> c = fourCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
        v.push_back(mkScen(c[i], "", [](KV & k) { k.degreeDecrease(1, true); },
                           capMultiplicity(opDegreeDecrease(c[i], 1), c[i].p - 1)));
    return v;
}

std::vector<Scen> scIncMult()
{
    const std::vector<PeriodDesc> c = fourCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
        for (int bnd = 0; bnd != 2; ++bnd)
            v.push_back(mkScen(c[i], bnd ? "boundary=1" : "boundary=0",
                               [=](KV & k) { k.increaseMultiplicity(1, 0 != bnd); },
                               opReduceContinuity(c[i], 1)));
    return v;
}

std::vector<Scen> scRedMult()
{
    const std::vector<PeriodDesc> c = fiveCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
        for (int bnd = 0; bnd != 2; ++bnd)
            v.push_back(mkScen(c[i], bnd ? "boundary=1" : "boundary=0",
                               [=](KV & k) { k.reduceMultiplicity(1, 0 != bnd); },
                               opElevateContinuity(c[i], 1)));
    return v;
}

std::vector<Scen> scCoarsen(const std::vector<PeriodDesc> & c)
{
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
    {
        const PeriodDesc dp = opCoarsen(c[i], 1, 1);
        std::shared_ptr<std::vector<real_t> > got(new std::vector<real_t>());
        Scen s = mkScen(c[i], "", [got](KV & k) { *got = k.coarsen(1, 1); }, dp);
        s.got = got;
        s.hasRet = true;
        std::set_difference(c[i].s.begin(), c[i].s.end(), dp.s.begin(), dp.s.end(),
                            std::back_inserter(s.removed));
        v.push_back(s);
    }
    return v;
}

std::vector<Scen> scSetDegree()
{
    const std::vector<PeriodDesc> c = fourCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
    {
        const short_t up = (short_t)(c[i].p + 1), down = (short_t)(c[i].p - 1);
        v.push_back(mkScen(c[i], "p+1", [=](KV & k) { k.set_degree(up); },
                           opDegreeIncrease(c[i], 1)));
        v.push_back(mkScen(c[i], "p-1", [=](KV & k) { k.set_degree(down); },
                           opDegreeDecrease(c[i], 1)));
    }
    return v;
}

std::vector<Scen> scReverse()
{
    const std::vector<PeriodDesc> c = fourCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
        v.push_back(mkScen(c[i], "", [](KV & k) { k.reverse(); }, opReverse(c[i])));
    return v;
}

std::vector<Scen> scAffine()
{
    const std::vector<PeriodDesc> c = fourCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
    {
        PeriodDesc dp = c[i];
        const real_t nb = -2 + 2 * c[i].L();
        for (size_t j = 0; j != dp.s.size(); ++j)
            dp.s[j] = -2 + 2 * (c[i].s[j] - c[i].a);
        dp.a = -2;
        dp.b = nb;
        v.push_back(mkScen(c[i], "", [=](KV & k) { k.affineTransformTo(-2, nb); }, dp));
    }
    return v;
}

std::vector<Scen> scTransform()
{
    const std::vector<PeriodDesc> c = fourCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
    {
        const real_t na = c[i].a + .5, nb = c[i].b + .5;
        v.push_back(mkScen(c[i], "", [=](KV & k) { k.transform(na, nb); },
                           shifted(c[i], .5)));
    }
    return v;
}

std::vector<Scen> scAddConstant()
{
    const std::vector<PeriodDesc> c = fourCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
        v.push_back(mkScen(c[i], "", [](KV & k) { k.addConstant(real_t(.5)); },
                           shifted(c[i], .5)));
    return v;
}

std::vector<Scen> scAddConstantStart()
{
    const std::vector<PeriodDesc> c = fourCases();
    std::vector<Scen> v;
    for (size_t i = 0; i != c.size(); ++i)
    {
        const real_t t0 = stretchedKnots(c[i])[0];
        v.push_back(mkScen(c[i], "start=t0", [=](KV & k) { k.addConstant(t0, real_t(.5)); },
                           shifted(c[i], .5)));
        v.push_back(mkScen(c[i], "start=t0-1", [=](KV & k) { k.addConstant(t0 - 1, real_t(.5)); },
                           shifted(c[i], .5)));
    }
    return v;
}

/// Flattened values of a matrix, column by column.
std::vector<real_t> flat(const gsMatrix<real_t> & m)
{
    return std::vector<real_t>(m.data(), m.data() + m.size());
}

/// Every const query of greville on kv, flattened, including the
/// pointwise greville(int); records whether anything threw.
std::vector<real_t> grevilleAll(const KV & kv, index_t count, std::string & threw,
                                index_t & cols)
{
    std::vector<real_t> all;
    try
    {
        const gsMatrix<real_t> g0 = kv.greville();
        gsMatrix<real_t> g1, g2;
        kv.greville_into(g1, true);
        kv.greville_into(g2, false);
        cols = g0.cols();
        const std::vector<real_t> f0 = flat(g0), f1 = flat(g1), f2 = flat(g2);
        all.insert(all.end(), f0.begin(), f0.end());
        all.insert(all.end(), f1.begin(), f1.end());
        all.insert(all.end(), f2.begin(), f2.end());
        for (index_t i = 0; i != count; ++i)
            all.push_back(kv.greville((int)i));
    }
    catch (std::exception & e) { threw = e.what(); }
    return all;
}

void grevilleOn(const PeriodDesc & d)
{
    const KV per = kvConvert(d), twin = stretchedKV(d);
    const KvSnap s0 = kvSnap(per);
    std::string threwPer, threwTwin;
    index_t cols = -1, colsTwin = -1;
    const index_t count = d.n() + d.m();
    const std::vector<real_t> got = grevilleAll(per, count, threwPer, cols);
    const std::vector<real_t> ref = grevilleAll(twin, count, threwTwin, colsTwin);
    CHECK_EQUAL(std::string(""), tagged(d.tag, threwPer));
    CHECK_EQUAL(std::string(""), tagged(d.tag, threwTwin));
    CHECK_EQUAL(std::string(""), tagged(d.tag, whyOf(knotsClose(got, ref, d.scale()))));
    CHECK_EQUAL(fzN(count), cols);
    CHECK_EQUAL(std::string(""), tagged(d.tag, whyOf(kvUnchanged(s0, per))));
}

/// Expected ghost-free reference result of an init* call on a vector that
/// was periodic: the same call on a default-constructed vector.
void initResets(const KOp & init)
{
    KV kv = kvConvert(caseC3n()), ref;
    init(kv);
    init(ref);
    CHECK(fzB(!kv.isPeriodic()));
    CHECK_EQUAL(fzN(0), kv.numCrossingFunctions());
    CHECK(fzB(kv == ref));
}

} // namespace

SUITE(gsPeriodicKnotVector_test)
{

// ===========================================================================
// Conversion with setPeriodic
// ===========================================================================

TEST(space_NEW_kvConvertQueries_ConvAll)
{
    const std::string fails = runCases(constructibleCases(), [](const PeriodDesc & d) -> std::string
    {
        std::string w;
        for (int k = 0; k != 2; ++k)
        {
            KV kv = stretchedKV(d);
            if (kv.isPeriodic() || 0 != kv.numCrossingFunctions())
                w += "flag set before conversion; ";
            if (0 == k) kv.setPeriodic(); else kv.setPeriodic(true);
            if (!fzB(kv.isPeriodic()))
                w += "not periodic after conversion; ";
            if (kv.periodicSize() != fzN(d.n()))
                w += "n " + num(kv.periodicSize()) + " != " + num(d.n()) + "; ";
            w += whyOf(scalarClose(kv.period(), d.L(), d.scale()));
        }
        return w;
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_kvUnwrap_ConvAll)
{
    const std::string fails = runCases(constructibleCases(), [](const PeriodDesc & d) -> std::string
    {
        std::string w;
        KV kv = kvConvert(d);
        const std::vector<real_t> knots = kvKnots(kv);
        const short_t deg = kv.degree();
        kv.setPeriodic(false);
        if (!fzB(!kv.isPeriodic()) || fzN(0) != kv.numCrossingFunctions())
            w += "flag or m after unwrap; ";
        if (kvKnots(kv) != knots)
            w += "knots changed by unwrap; ";
        if (kv.degree() != deg)
            w += "degree changed by unwrap; ";
        if (!fzB(kv == stretchedKV(d)))
            w += "unwrapped vector differs from the raw stretched vector; ";
        kv.setPeriodic(true);
        w += kvState(kv, d);
        return w;
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(space_NEW_kvSeamTolFar_E6)
{
    const PeriodDesc d = caseE6();
    KV kv = kvConvert(d);
    CHECK_EQUAL(std::string(""), kvState(kv, d));

    // The seam tolerance is 100*eps*max(|a|,|b|,L) ~ 2.2e-8 here, far above
    // 100*eps*L (below one ulp of 1e6): a deviation of 1e-9 in a constrained
    // ghost knot is accepted and regenerated, 1e-7 and 1e-4 are refused.
    std::vector<real_t> t = stretchedKnots(d);
    t[0] += 1e-9;
    KV near(t, 2);
    near.setPeriodic(true);
    CHECK_EQUAL(std::string(""), kvState(near, d));

    t = stretchedKnots(d);
    t[0] += 1e-7;
    KV mid(t, 2);
    kvRefuse(mid, [](KV & k) { k.setPeriodic(true); }, "delta=1e-7");

    t = stretchedKnots(d);
    t[0] += 1e-4;
    KV far(t, 2);
    kvRefuse(far, [](KV & k) { k.setPeriodic(true); }, "delta=1e-4");
}

TEST(space_B23_kvClampedSmallNRefused_C1kP3N1)
{
    KV kv(real_t(0), real_t(1), 0u, 4, 1, 3);
    kvRefuse(kv, [](KV & k) { k.setPeriodic(true); }, "C1kP3N1");
}

TEST(space_B22_kvFirstSpanRefused_C3x)
{
    KV kv(knotsC3x(), 2);
    kvRefuse(kv, [](KV & k) { k.setPeriodic(true); }, "C3x");
}

TEST(space_B22_kvFreeKnotAccepted_C3y)
{
    const PeriodDesc d = caseC3y();
    KV kv = stretchedKV(d);
    CHECK_EQUAL(std::string(""), thrown([&]() { kv.setPeriodic(true); }));
    CHECK_EQUAL(std::string(""), kvState(kv, d));
}

TEST(space_B23_kvSmallNRefused_C5a)
{
    KV kv = stretchedKV(caseC5a());
    kvRefuse(kv, [](KV & k) { k.setPeriodic(true); }, "C5a");
}

TEST(space_B23_kvSmallNRefused_C5b)
{
    KV kv = stretchedKV(caseC5b());
    kvRefuse(kv, [](KV & k) { k.setPeriodic(true); }, "C5b");
}

TEST(space_B23_kvSmallNRefused_C2uP4N4)
{
    KV kv = stretchedKV(caseC2u(4, 4));
    kvRefuse(kv, [](KV & k) { k.setPeriodic(true); }, "C2uP4N4");
}

TEST(space_B31_kvAsymSeamRefused_S31)
{
    KV kv(knotsS31(), 2);
    kvRefuse(kv, [](KV & k) { k.setPeriodic(true); }, "S31");
}

TEST(space_NEW_kvNonPeriodicQueries_C6)
{
    const KV kv = c6Knots();
    CHECK(fzB(!kv.isPeriodic()));
    CHECK_EQUAL(fzN(0), kv.numCrossingFunctions());
    // 6 functions of degree 2 on [0,.25,.5,.75,1]
    CHECK_EQUAL(fzN(6), kv.periodicSize());
    CHECK_EQUAL(std::string(""), whyOf(scalarClose(kv.period(), 1, 1)));
}

TEST(space_NEW_kvSeamRoundoffRefused_ClampedRandom)
{
    const real_t a = -7.312715117751976, b = 1.1631479142509775, mid = (a + b) / 2;
    const real_t k[] = {a, a, a, mid, b, b, b};
    KV kv(vec(k, 7), 2);
    kv.setPeriodic();
    // s = {a, a, mid}: mu = p = 2 copies of a, one crossing function.
    CHECK_EQUAL(fzN(3), kv.periodicSize());
    CHECK_EQUAL(fzN(1), kv.numCrossingFunctions());

    const real_t nearB = std::nextafter(b, real_t(-1e300));
    kvRefuse(kv, [=](KV & q) { q.insert(nearB); }, "insert(nextafter(b))");

    const KvSnap s0 = kvSnap(kv);
    CHECK_EQUAL(std::string(""), thrown([&]() { kv.insert(real_t(0)); }));
    CHECK_EQUAL(fzN(s0.n + 1), kv.periodicSize());
    CHECK_EQUAL(fzN(1), kv.numCrossingFunctions());
    CHECK_EQUAL(s0.size + 1, kv.size());
}

// ===========================================================================
// ghost_ tests: the stored stretched knots and the crossing count m.
// Diagnostic only; they describe the current stretched-knot design.
// ===========================================================================

TEST(ghost_NEW_kvConvertKnots_ConvAll)
{
    const std::string fails = runCases(constructibleCases(), [](const PeriodDesc & d)
    {
        return kvState(kvConvert(d), d);
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(ghost_NEW_kvConvertKnots_C1cAll)
{
    std::vector<PeriodDesc> cases = c1cCases();
    cases.push_back(caseC1c(3, 4, 2));
    const std::string fails = runCases(cases, [](const PeriodDesc & d)
    {
        return kvState(kvConvert(d), d);
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(ghost_B1_kvClampedConv_C1kAll)
{
    const std::string fails = runCases(c1kCases(), [](const PeriodDesc & d) -> std::string
    {
        std::vector<real_t> ub(d.s.begin(), d.s.end());
        ub.push_back(d.b);
        const index_t N = (index_t)uniqueSorted(ub).size() - 1;
        KV kv(d.a, d.b, (unsigned)(N - 1), d.p + 1, 1, d.p);
        kv.setPeriodic(true);
        std::string w = kvState(kv, d);
        if (kv.numCrossingFunctions() != fzN(1))
            w += "m != 1; ";
        if (kv.periodicSize() != fzN(N + d.p - 1))
            w += "n " + num(kv.periodicSize()) + " != N+p-1; ";
        return w;
    });
    CHECK_EQUAL(std::string(""), fails);
}

TEST(ghost_B29_kvIdempotent_ConvAll)
{
    const std::string fails = runCases(constructibleCases(), [](const PeriodDesc & d) -> std::string
    {
        std::string w;
        KV kv = kvConvert(d);
        const std::vector<real_t> before = kvKnots(kv);
        kv.setPeriodic(true);
        if (!fzB(kvKnots(kv) == before))
            w += "knots changed by a second setPeriodic; ";
        w += kvState(kv, d);
        return w;
    });
    CHECK_EQUAL(std::string(""), fails);
}

// --- mutator rows: act on the period part, then regenerate -----------------

TEST(ghost_NEW_kvInsertInside_Four)    { CHECK(runAll(scInsertInside()) > 0); }
TEST(ghost_NEW_kvInsertMult_Four)      { CHECK(runAll(scInsertMult()) > 0); }
TEST(ghost_NEW_kvInsertWrap_Four)      { CHECK(runAll(scInsertWrap()) > 0); }
TEST(ghost_NEW_kvInsertRange_Four)     { CHECK(runAll(scInsertRange()) > 0); }
TEST(ghost_NEW_kvInsertVec_Four)       { CHECK(runAll(scInsertVec()) > 0); }
TEST(ghost_NEW_kvRemoveValue_Four)     { CHECK(runAll(scRemoveValue()) > 0); }
TEST(ghost_NEW_kvRemoveWrap_C3n)       { CHECK(runAll(scRemoveWrap()) > 0); }
TEST(ghost_NEW_kvRemoveIter_Four)      { CHECK(runAll(scRemoveIter()) > 0); }
TEST(ghost_NEW_kvUniformRefine_Four)   { CHECK(runAll(scUniformRefine(1)) > 0); }
TEST(ghost_NEW_kvUniformRefineMul2_Four) { CHECK(runAll(scUniformRefine(2)) > 0); }
TEST(ghost_NEW_kvRefineSpansMult_Four) { CHECK(runAll(scRefineSpansMult()) > 0); }
TEST(ghost_NEW_kvRefineSpansUnsigned_Four) { CHECK(runAll(scRefineSpansUnsigned()) > 0); }
TEST(ghost_NEW_kvRefineSpansIter_C3n)  { CHECK(runAll(scRefineSpansIter()) > 0); }
TEST(ghost_NEW_kvDegreeElevate_Four)   { CHECK(runAll(scDegreeElevate()) > 0); }
TEST(ghost_NEW_kvDegreeIncrease_Four)  { CHECK(runAll(scDegreeIncrease()) > 0); }
TEST(ghost_NEW_kvDegreeReduce_Five)    { CHECK(runAll(scDegreeReduce()) > 0); }
TEST(ghost_NEW_kvDegreeDecrease_Four)  { CHECK(runAll(scDegreeDecrease()) > 0); }
TEST(ghost_NEW_kvDegreeDecInterior_Four) { CHECK(runAll(scDegreeDecInterior()) > 0); }
TEST(ghost_NEW_kvIncMult_Four)         { CHECK(runAll(scIncMult()) > 0); }
TEST(ghost_NEW_kvRedMult_Five)         { CHECK(runAll(scRedMult()) > 0); }
TEST(ghost_NEW_kvCoarsen_Four)         { CHECK(runAll(scCoarsen(fourCases())) > 0); }
TEST(ghost_NEW_kvCoarsen_C2uP2N8)      { CHECK(runAll(scCoarsen(std::vector<PeriodDesc>(1, caseC2u(2, 8)))) > 0); }
TEST(ghost_NEW_kvSetDegree_Four)       { CHECK(runAll(scSetDegree()) > 0); }
TEST(ghost_NEW_kvReverse_Four)         { CHECK(runAll(scReverse()) > 0); }
TEST(ghost_NEW_kvAffine_Four)          { CHECK(runAll(scAffine()) > 0); }
TEST(ghost_NEW_kvTransform_Four)       { CHECK(runAll(scTransform()) > 0); }
TEST(ghost_NEW_kvAddConstant_Four)     { CHECK(runAll(scAddConstant()) > 0); }
TEST(ghost_NEW_kvAddConstantStart_Four) { CHECK(runAll(scAddConstantStart()) > 0); }

TEST(ghost_NEW_kvRemoveLeftGhostIter_Clamped)
{
    // clamped quadratic on [0,1]; periodic: s = {0,0,.1,.5,.9}, stretched
    // knots -.1 0 0 .1 .5 .9 1 1 1.1.
    const real_t k[] = {0, 0, 0, .1, .5, .9, 1, 1, 1};
    const real_t s0[] = {0, 0, .1, .5, .9}, s1[] = {0, 0, .1, .5};
    const PeriodDesc d = mkDesc("Clamped", 2, 0, 1, vec(s0, 5));
    const PeriodDesc dp = mkDesc("Clamped-ghost", 2, 0, 1, vec(s1, 4));
    KV kv(vec(k, 9), 2);
    kv.setPeriodic();
    CHECK_EQUAL(std::string(""), kvState(kv, d, false));
    CHECK_EQUAL(std::string(""), thrown([&]() { kv.remove(kv.ubegin(), 1); }));
    CHECK_EQUAL(std::string(""), kvState(kv, dp));
}

// ===========================================================================
// Const members on a periodic vector
// ===========================================================================

TEST(space_NEW_kvGreville_C2uP2N4)  { grevilleOn(caseC2u(2, 4)); }
TEST(space_NEW_kvGreville_C3n)      { grevilleOn(caseC3n()); }

// ===========================================================================
// Explicit refusals
// ===========================================================================

TEST(space_NEW_kvInsertSeamRefused_Four)
{
    const std::vector<PeriodDesc> c = fourCases();
    for (size_t i = 0; i != c.size(); ++i)
    {
        const real_t a = c[i].a, b = c[i].b, L = c[i].L();
        const real_t xs[] = {a, b, a - L, b + L};
        for (size_t j = 0; j != 4; ++j)
        {
            KV kv = kvConvert(c[i]);
            kvRefuse(kv, insertOp(xs[j]), c[i].tag + " " + labelOf(xs[j]));
        }
    }
}

TEST(space_NEW_kvInsertFarRefused_Four)
{
    const std::vector<PeriodDesc> c = fourCases();
    for (size_t i = 0; i != c.size(); ++i)
    {
        const real_t xs[] = {c[i].a - 1.5 * c[i].L(), c[i].b + 1.5 * c[i].L()};
        for (size_t j = 0; j != 2; ++j)
        {
            KV kv = kvConvert(c[i]);
            kvRefuse(kv, insertOp(xs[j]), c[i].tag + " " + labelOf(xs[j]));
        }
    }
}

TEST(space_NEW_kvInsertRangeAtomic_C3n)
{
    KV kv = kvConvert(caseC3n());
    const real_t x[] = {.25, 4};
    const std::vector<real_t> v = vec(x, 2);
    kvRefuse(kv, [=](KV & k) { k.insert(v.begin(), v.end()); }, "{.25,4}");
}

TEST(space_NEW_kvRemoveSeamRefused_C3n)
{
    const real_t xs[] = {0, 4};
    for (size_t j = 0; j != 2; ++j)
    {
        KV kv = kvConvert(caseC3n());
        const real_t x = xs[j];
        kvRefuse(kv, [=](KV & k) { k.remove(x, 1); }, labelOf(x));
    }
}

TEST(space_NEW_kvRemoveSeamRefused_E2)
{
    KV kv = kvConvert(caseE2());
    kvRefuse(kv, [](KV & k) { k.remove(real_t(0), 1); }, "x=0");
}

TEST(space_NEW_kvRemoveSmallNRefused_C2uP3N4)
{
    KV kv = kvConvert(caseC2u(3, 4));
    kvRefuse(kv, [](KV & k) { k.remove(real_t(.25), 1); }, "x=.25");
}

TEST(space_NEW_kvRemoveNonKnotIgnored_C3n)
{
    // .25 lies inside the period and is not a knot.
    KV kv = kvConvert(caseC3n());
    const KvSnap s0 = kvSnap(kv);
    CHECK_EQUAL(std::string(""), thrown([&]() { kv.remove(real_t(.25), 1); }));
    CHECK_EQUAL(std::string(""), whyOf(kvUnchanged(s0, kv)));
}

TEST(space_NEW_kvDegreeIncRefused_C2uP3N4)
{
    KV kv = kvConvert(caseC2u(3, 4));
    kvRefuse(kv, [](KV & k) { k.degreeIncrease(1); }, "C2uP3N4");
}

TEST(space_NEW_kvDegreeDecRefused_C2uP1N4)
{
    KV kv = kvConvert(caseC2u(1, 4));
    kvRefuse(kv, [](KV & k) { k.degreeDecrease(1); }, "C2uP1N4");
}

TEST(space_NEW_kvSetDegreeRefused_Four)
{
    const std::vector<PeriodDesc> c = fourCases();
    for (size_t i = 0; i != c.size(); ++i)
    {
        KV kv = kvConvert(c[i]);
        kvRefuse(kv, [](KV & k) { k.set_degree(0); }, c[i].tag);
    }
}

TEST(space_NEW_kvAddConstantStartRefused_Four)
{
    const std::vector<PeriodDesc> c = fourCases();
    for (size_t i = 0; i != c.size(); ++i)
    {
        KV kv = kvConvert(c[i]);
        const real_t a = c[i].a;
        kvRefuse(kv, [=](KV & k) { k.addConstant(a, real_t(.5)); }, c[i].tag);
    }
}

TEST(space_NEW_kvInsertMultRefused_C4i)
{
    // multiplicity of 1 would become 2+2 = 4 > p+1 = 3.
    KV kv = kvConvert(caseC4i());
    kvRefuse(kv, [](KV & k) { k.insert(real_t(1), (index_t)2); }, "x=1,mult=2");
}

// --- refused raw surgery ---------------------------------------------------

TEST(space_NEW_kvTrimLeftRefused_C3n)
{
    KV kv = kvConvert(caseC3n());
    kvRefuse(kv, [](KV & k) { k.trimLeft(1); }, "C3n");
}

TEST(space_NEW_kvTrimRightRefused_C3n)
{
    KV kv = kvConvert(caseC3n());
    kvRefuse(kv, [](KV & k) { k.trimRight(1); }, "C3n");
}

TEST(space_NEW_kvAppendRefused_C3n)
{
    KV kv = kvConvert(caseC3n());
    const std::vector<real_t> v(1, real_t(2.5));
    kvRefuse(kv, [=](KV & k) { k.append(v.begin(), v.end()); }, "C3n");
}

TEST(space_NEW_kvEraseRefused_C3n)
{
    KV kv = kvConvert(caseC3n());
    kvRefuse(kv, [](KV & k) { k.erase(0, 1); }, "C3n");
}

// ===========================================================================
// Flag reset and carry
// ===========================================================================

TEST(space_NEW_kvInitResets_C3n)
{
    std::vector<KOp> inits;
    inits.push_back([](KV & k) { k.initUniform(real_t(0), real_t(1), 3u, 3u, 1u, (short_t)2); });
    inits.push_back([](KV & k) { k.initUniform(5u, 3u, 1u, (short_t)2); });
    inits.push_back([](KV & k) { k.initClamped(2, 4u, 1u); });
    inits.push_back([](KV & k) { k.initClamped(real_t(0), real_t(1), 2, 3u, 1u); });
    inits.push_back([](KV & k) { k.initGraded(5u, 2, real_t(0.5), 1u); });
    inits.push_back([](KV & k) { k.initGraded(real_t(0), real_t(1), 3u, 2, real_t(0.5), 1u); });
    for (size_t i = 0; i != inits.size(); ++i)
        initResets(inits[i]);
}

TEST(space_NEW_kvRawCtorNotPeriodic_C3n)
{
    const KV per = kvConvert(caseC3n());
    const KV fromKnots(kvKnots(per), 2), empty, uniform(real_t(0), real_t(1), 3u, 3u, 1u, 2);
    CHECK(fzB(!fromKnots.isPeriodic()));
    CHECK(!empty.isPeriodic());
    CHECK(!uniform.isPeriodic());
    CHECK_EQUAL(fzN(0), fromKnots.numCrossingFunctions());
    CHECK_EQUAL(index_t(0), empty.numCrossingFunctions());
    CHECK_EQUAL(index_t(0), uniform.numCrossingFunctions());
}

TEST(space_NEW_kvCopyCarries_C3n)
{
    const PeriodDesc d = caseC3n();
    const KV src = kvConvert(d);

    const KV copy(src);
    KV assigned = c6Knots();
    assigned = src;
    const memory::unique_ptr<KV> cloned(src.clone());

    CHECK_EQUAL(std::string(""), tagged("copy ctor", kvState(copy, d)));
    CHECK_EQUAL(std::string(""), tagged("copy assignment", kvState(assigned, d)));
    CHECK_EQUAL(std::string(""), tagged("clone", kvState(*cloned, d)));
    CHECK(fzB(copy == src));
    CHECK(assigned == src);
    CHECK(*cloned == src);
}

TEST(space_NEW_kvMoveCarries_C3n)
{
    const PeriodDesc d = caseC3n();
    KV s1 = kvConvert(d), s2 = kvConvert(d);

    const KV moved(std::move(s1));
    KV assigned = c6Knots();
    assigned = std::move(s2);

    CHECK_EQUAL(std::string(""), tagged("move ctor", kvState(moved, d)));
    CHECK_EQUAL(std::string(""), tagged("move assignment", kvState(assigned, d)));
    CHECK(fzB(moved == assigned));
}

TEST(space_NEW_kvSwapCarries_C3n)
{
    const PeriodDesc d = caseC3n();
    KV p = kvConvert(d), q = c6Knots();
    p.swap(q);
    CHECK_EQUAL(std::string(""), kvState(q, d));
    CHECK(fzB(p == c6Knots()));
    CHECK(!p.isPeriodic());
}

TEST(space_NEW_kvEqualityFlag_C3n)
{
    const PeriodDesc d = caseC3n();
    const KV a = kvConvert(d), b = kvConvert(d), raw = stretchedKV(d);
    CHECK(fzB(a == b));
    CHECK(!(a == raw));
    CHECK(fzB(a != raw));
    KV unwrapped = kvConvert(d);
    unwrapped.setPeriodic(false);
    CHECK(fzB(unwrapped == raw));
}

// ===========================================================================
// Non-periodic twins: what the same operation does on the flag-off vector
// ===========================================================================

TEST(info_twin_kvInsertInside_Four)    { twinAll("kvInsertInside", scInsertInside()); }
TEST(info_twin_kvInsertMult_Four)      { twinAll("kvInsertMult", scInsertMult()); }
TEST(info_twin_kvInsertWrap_Four)      { twinAll("kvInsertWrap", scInsertWrap()); }
TEST(info_twin_kvInsertRange_Four)     { twinAll("kvInsertRange", scInsertRange()); }
TEST(info_twin_kvInsertVec_Four)       { twinAll("kvInsertVec", scInsertVec()); }
TEST(info_twin_kvRemoveValue_Four)     { twinAll("kvRemoveValue", scRemoveValue()); }
TEST(info_twin_kvRemoveWrap_C3n)       { twinAll("kvRemoveWrap", scRemoveWrap()); }
TEST(info_twin_kvRemoveIter_Four)      { twinAll("kvRemoveIter", scRemoveIter()); }
TEST(info_twin_kvUniformRefine_Four)   { twinAll("kvUniformRefine", scUniformRefine(1)); }
TEST(info_twin_kvUniformRefineMul2_Four) { twinAll("kvUniformRefineMul2", scUniformRefine(2)); }
TEST(info_twin_kvRefineSpansMult_Four) { twinAll("kvRefineSpansMult", scRefineSpansMult()); }
TEST(info_twin_kvRefineSpansUnsigned_Four) { twinAll("kvRefineSpansUnsigned", scRefineSpansUnsigned()); }
TEST(info_twin_kvRefineSpansIter_C3n)  { twinAll("kvRefineSpansIter", scRefineSpansIter()); }
TEST(info_twin_kvDegreeElevate_Four)   { twinAll("kvDegreeElevate", scDegreeElevate()); }
TEST(info_twin_kvDegreeIncrease_Four)  { twinAll("kvDegreeIncrease", scDegreeIncrease()); }
TEST(info_twin_kvDegreeReduce_Five)    { twinAll("kvDegreeReduce", scDegreeReduce()); }
TEST(info_twin_kvDegreeDecrease_Four)  { twinAll("kvDegreeDecrease", scDegreeDecrease()); }
TEST(info_twin_kvDegreeDecInterior_Four) { twinAll("kvDegreeDecInterior", scDegreeDecInterior()); }
TEST(info_twin_kvIncMult_Four)         { twinAll("kvIncMult", scIncMult()); }
TEST(info_twin_kvRedMult_Five)         { twinAll("kvRedMult", scRedMult()); }
TEST(info_twin_kvCoarsen_Four)         { twinAll("kvCoarsen", scCoarsen(fourCases())); }
TEST(info_twin_kvCoarsen_C2uP2N8)      { twinAll("kvCoarsen", scCoarsen(std::vector<PeriodDesc>(1, caseC2u(2, 8)))); }
TEST(info_twin_kvSetDegree_Four)       { twinAll("kvSetDegree", scSetDegree()); }

} // SUITE
