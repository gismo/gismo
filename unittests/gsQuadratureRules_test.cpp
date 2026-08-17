/** @file gsQuadratureRules_test.cpp

    @brief Tests available gsQuadratureRules (containing gsGaussRules 1-21)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Vogl, A. Mantzaflaris
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework

#include <gsDomain/gsImplicitTrimmedDomain.h>

namespace {

/// RAII capture of std::cout (gsWarn/gsInfo both write there, gsDebug.h:43,51).
struct CoutCapture
{
    std::ostringstream oss;
    std::streambuf *   old;
    CoutCapture() : old(std::cout.rdbuf(oss.rdbuf())) {}
    ~CoutCapture() { std::cout.rdbuf(old); }
    std::string str() const { return oss.str(); }
};

/// A fully-interior gsImplicitTrimmedDomain over the unit square, built with
/// m_deg == 1 (constant level set -1 everywhere -> no trimming), used by the
/// explicit-degrees tests (1 and 2) to isolate the "degrees" parameter from
/// domain.degree() -- both would otherwise be 1 and the tests would not
/// distinguish the explicit-degrees path from the fallback.
inline gsImplicitTrimmedDomain<2,real_t> unitSquareDomainDeg1()
{
    gsFunctionExpr<real_t> phi("-1", 2);
    gsMatrix<real_t> bbox(2,2);
    bbox << 0.0, 0.0,     // lower corners
            1.0, 1.0;     // upper corners
    gsVector<index_t,2> nc; nc << 4, 4;
    return gsImplicitTrimmedDomain<2,real_t>(phi, bbox, nc, 5, /*deg=*/1);
}

} // anonymous namespace

SUITE(gsQuadratureRules_test)                 // The suite should have the same name as the file
{

void testWork(const index_t nodes[], size_t dim);

real_t calcAntiDerivative(gsVector<index_t> const &deg, const index_t dim);

real_t calcPoly(gsVector<index_t> const &deg,
                gsQuadRule<real_t> const &gr,
                const index_t dim);

const char *addPlus(const index_t d, index_t i);

TEST(tensor_quad_3)
{
    index_t array[] = {3};
    testWork(array, 1);
}

TEST(tensor_quad_5)
{
    index_t array[] = {5};
    testWork(array, 1);
}

TEST(tensor_quad_1_5_9)
{
    index_t array[] = {1, 5, 9};
    testWork(array, 3);
}

TEST(tensor_quad_2_6_7)
{
    index_t array[] = {2, 6, 7};
    testWork(array, 3);
}

TEST(tensor_quad_3_4_8)
{
    index_t array[] = {3, 4, 8};
    testWork(array, 3);
}

TEST(tensor_quad_10_17)
{
    index_t array[] = {10, 17};
    testWork(array, 2);
}

TEST(tensor_quad_11_16)
{
    index_t array[] = {11, 16};
    testWork(array, 2);
}

TEST(tensor_quad_12_15)
{
    index_t array[] = {12, 15};
    testWork(array, 2);
}

TEST(tensor_quad_13_14)
{
    index_t array[] = {13, 14};
    testWork(array, 2);
}

TEST(tensor_quad_18)
{
    index_t array[] = {18};
    testWork(array, 1);
}

TEST(tensor_quad_19)
{
    index_t array[] = {19};
    testWork(array, 1);
}

TEST(tensor_quad_20)
{
    index_t array[] = {20};
    testWork(array, 1);
}

TEST(tensor_quad_21)
{
    index_t array[] = {21};
    testWork(array, 1);
}

TEST(tensor_quad_22)
{
    index_t array[] = {22};
    testWork(array, 1);
}

/// Test 1 -- gsQuadrature::numNodes(domain, quA, quB, fixDir, degrees) uses
/// the explicitly supplied per-direction degrees and falls back to
/// domain.degree(i) only when degrees is empty.
///
/// The domain is built with m_deg == 1 (domain.degree(i) == 1 for both
/// directions), so any node count that used the explicit degrees {3,3} or
/// {3,2} instead of the fallback is unambiguous.
TEST(explicitDegrees_numNodes)
{
    gsImplicitTrimmedDomain<2,real_t> domain = unitSquareDomainDeg1();

    // Sanity: the domain's own guessed degree is 1 in both directions.
    CHECK_EQUAL((short_t)1, domain.degree(0));
    CHECK_EQUAL((short_t)1, domain.degree(1));

    // Explicit isotropic degrees {3,3}: quA*deg+quB+0.5 = 1*3+1 = 4.
    gsVector<short_t> degs33(2); degs33 << 3, 3;
    gsVector<index_t> nn33 = gsQuadrature::numNodes(domain, 1.0, 1, -1, degs33);
    CHECK_EQUAL(4, nn33[0]);
    CHECK_EQUAL(4, nn33[1]);

    // Empty degrees: falls back to domain.degree(i) == 1 -> 1*1+1 = 2.
    gsVector<index_t> nnEmpty = gsQuadrature::numNodes(domain, 1.0, 1, -1, gsVector<short_t>());
    CHECK_EQUAL(2, nnEmpty[0]);
    CHECK_EQUAL(2, nnEmpty[1]);

    // Anisotropic guard: {3,2} -> (4,3). Catches a direction collapse or a
    // maxCoeff() shortcut leaking into numNodes.
    gsVector<short_t> degs32(2); degs32 << 3, 2;
    gsVector<index_t> nn32 = gsQuadrature::numNodes(domain, 1.0, 1, -1, degs32);
    CHECK_EQUAL(4, nn32[0]);
    CHECK_EQUAL(3, nn32[1]);
}

/// Test 2 -- gsQuadrature::getPtr(domain, options, fixDir, degrees) forwards
/// the same degrees into the constructed rule (node count).
TEST(explicitDegrees_getPtr)
{
    gsImplicitTrimmedDomain<2,real_t> domain = unitSquareDomainDeg1();

    gsOptionList opt;
    opt.addInt ("quRule", "", gsQuadrature::GaussLegendre);
    opt.addReal("quA",    "", 1.0);
    opt.addInt ("quB",    "", 1);

    gsVector<short_t> degs33(2); degs33 << 3, 3;
    typename gsQuadRule<real_t>::uPtr ruleExplicit =
        gsQuadrature::getPtr<real_t>(domain, opt, -1, degs33);
    CHECK_EQUAL(16, ruleExplicit->numNodes()); // 4 x 4

    // Contrast: with empty degrees, the fallback domain.degree(i) == 1 gives
    // 2 nodes/direction, i.e. 4 total.
    typename gsQuadRule<real_t>::uPtr ruleFallback =
        gsQuadrature::getPtr<real_t>(domain, opt, -1, gsVector<short_t>());
    CHECK_EQUAL(4, ruleFallback->numNodes()); // 2 x 2
}

/// Test 3 -- numNodes(basis, ...) and numNodes(*basis.domain(), ...) agree
/// for both gsTensorBSplineBasis and gsTHBSplineBasis: guards the reroute of
/// the basis overload away from a naive *basis.domain() forward that would
/// drop the (accurate) per-basis degrees. Anisotropic degrees (2 in x, 3 in
/// y) so a direction collapse or a maxDegree() shortcut cannot hide.
TEST(basisDomainDegreeAgreement)
{
    gsKnotVector<real_t> kvx(0, 1, 3, 3); // degree 2
    gsKnotVector<real_t> kvy(0, 1, 3, 4); // degree 3
    gsTensorBSplineBasis<2,real_t> tb(kvx, kvy);
    gsTHBSplineBasis<2,real_t> thb(tb);

    CHECK_EQUAL((short_t)2, tb.degree(0));
    CHECK_EQUAL((short_t)3, tb.degree(1));
    CHECK_EQUAL((short_t)2, thb.degree(0));
    CHECK_EQUAL((short_t)3, thb.degree(1));

    {
        gsVector<index_t> nnBasis  = gsQuadrature::numNodes(tb, 1.0, 1, -1);
        gsVector<index_t> nnDomain = gsQuadrature::numNodes(*tb.domain(), 1.0, 1, -1);
        CHECK_EQUAL(nnDomain[0], nnBasis[0]);
        CHECK_EQUAL(nnDomain[1], nnBasis[1]);
        CHECK_EQUAL(tb.degree(0) + 1, nnBasis[0]);
        CHECK_EQUAL(tb.degree(1) + 1, nnBasis[1]);
    }
    {
        gsVector<index_t> nnBasis  = gsQuadrature::numNodes(thb, 1.0, 1, -1);
        gsVector<index_t> nnDomain = gsQuadrature::numNodes(*thb.domain(), 1.0, 1, -1);
        CHECK_EQUAL(nnDomain[0], nnBasis[0]);
        CHECK_EQUAL(nnDomain[1], nnBasis[1]);
        CHECK_EQUAL(thb.degree(0) + 1, nnBasis[0]);
        CHECK_EQUAL(thb.degree(1) + 1, nnBasis[1]);
    }
}

/// Test 5 -- the moment-fitting exactness guard's verdict (warn / no-warn)
/// is decided by _maxDegree(domain, degrees) alone: the two arms below use
/// IDENTICAL quA/quB and differ ONLY in the `degrees` argument, and they are
/// constructed to land on opposite sides of the `n < 2*deg+1` inequality
/// because of the degree source alone. If `_maxDegree` ignored `degrees`
/// (i.e. always returned domain.degree() == 1), both arms would compute the
/// SAME deg/n and give the SAME verdict, and arm 2's assertion below would
/// fail -- see the falsification requirement in the task report.
///
/// Arithmetic oracle (deg = _maxDegree(domain, degrees), quA=1.0, quB=2,
/// n = lround(quA*deg) + quB, warn iff n < 2*deg+1):
///
/// | arm | degrees | deg | n = lround(1.0*deg)+2 | 2*deg+1 | verdict       |
/// |-----|---------|-----|------------------------|---------|---------------|
/// | 1   | empty   |   1 | 1+2 = 3                | 3       | no warning    |
/// | 2   | {3,3}   |   3 | 3+2 = 5                | 7       | warns         |
///
/// Neither arm throws: gsMomentRule's output grid is n per direction
/// (3 in arm 1, 5 in arm 2), both >= 2, so gsMomentRule.h:228's
/// `at least 2 output points per direction` guard is satisfied in both
/// arms and the verdict below is decided by the warning assertion alone,
/// not by an exception.
TEST(momentFittingExactnessGuard)
{
    gsFunctionExpr<real_t> phi("sqrt(x^2+y^2)-1", 2);
    gsMatrix<real_t> bbox(2,2);
    bbox << -1.2, -1.0,
             1.2,  1.0;
    gsVector<index_t,2> nc; nc << 5, 4;
    gsImplicitTrimmedDomain<2,real_t> domain(phi, bbox, nc, 5, /*deg=*/1);

    gsVector<short_t> degs(2); degs << 3, 3;
    const std::string needle = "moment fitting is exact only up to degree n-1";

    // NOTE on ordering: gsQuadrature::makeMomentFittingPtr guards its warning
    // behind a one-shot `static bool warned` (gsQuadrature.h:806, accessed via
    // the public gsQuadrature::momentExactnessWarned() test hook). We reset
    // that flag immediately before EACH arm below so neither arm depends on
    // what any other test in the binary did first, but we ALSO keep the
    // no-warning case FIRST and the warning case SECOND in this same TEST
    // body: asserting "no warning" only proves something if a warning was
    // still possible at that point. If the positive (warning) case ran
    // first, the "no warning" assertion afterwards would be checking a
    // permanently spent flag and could never fail -- do not reorder these
    // two blocks.

    // 1) FIRST, must not warn: empty degrees -> deg = domain.degree() = 1,
    //    n = 3, 2*deg+1 = 3, 3 < 3 is false.
    {
        gsQuadrature::momentExactnessWarned() = false;

        gsOptionList opt;
        opt.addInt ("quRule",             "", gsQuadrature::MomentFittingRule);
        opt.addInt ("quMomentUnderlying", "", gsQuadrature::CutCellRule);
        opt.addReal("quA",                "", 1.0);
        opt.addInt ("quB",                "", 2);

        CoutCapture cap;
        typename gsQuadRule<real_t>::uPtr rule =
            gsQuadrature::getPtr<real_t>(domain, opt, -1, gsVector<short_t>());
        CHECK(cap.str().find(needle) == std::string::npos);
    }

    // 2) SECOND, must warn: explicit degrees {3,3} -> deg = maxCoeff = 3,
    //    n = 5, 2*deg+1 = 7, 5 < 7 is true. Same quA/quB as arm 1 -- the
    //    ONLY difference is the `degrees` argument.
    {
        gsQuadrature::momentExactnessWarned() = false;

        gsOptionList opt;
        opt.addInt ("quRule",             "", gsQuadrature::MomentFittingRule);
        opt.addInt ("quMomentUnderlying", "", gsQuadrature::CutCellRule);
        opt.addReal("quA",                "", 1.0);
        opt.addInt ("quB",                "", 2);

        CoutCapture cap;
        typename gsQuadRule<real_t>::uPtr rule =
            gsQuadrature::getPtr<real_t>(domain, opt, -1, degs);
        CHECK(cap.str().find(needle) != std::string::npos);
    }
}

/// Test 6 -- numNodes(domain, quA, quB, fixDir, degrees) and
/// getPtr(domain, options, fixDir, degrees) with a NON-NEGATIVE fixDir and
/// an anisotropic explicit `degrees`, pinning that:
///  (a) the fixed direction is dropped to 1 node, and
///  (b) each SURVIVING direction uses ITS OWN degree from `degrees`, not
///      degrees.maxCoeff() and not domain.degree(i).
///
/// Domain: unitSquareDomainDeg1(), domain.degree(i) == 1 in both directions.
/// degs32 = {3,2} is anisotropic on purpose, and degs32[0] == 3 is the
/// maximum, so fixDir = 0 (which drops direction 0 and keeps direction 1,
/// whose own degree is 2) separates "per-direction degrees" from
/// "the maximum": a maxCoeff() shortcut leaking into numNodes would report
/// 4 nodes in direction 1 instead of 3.
TEST(explicitDegrees_fixDir)
{
    gsImplicitTrimmedDomain<2,real_t> domain = unitSquareDomainDeg1();

    gsVector<short_t> degs32(2); degs32 << 3, 2;

    // numNodes: quA=1.0, quB=1, nnodes[i] = lround(quA*deg_i + quB + 0.5)...
    // (cast<Real,index_t>(quA*deg_i + quB + 0.5)); fixed direction -> 1 node.

    // fixDir = 0: direction 0 fixed -> 1 node; direction 1 survives with its
    // OWN degree 2 -> 1*2+1+0.5 -> 3. NOT 4 (which would be maxCoeff()==3
    // used instead), NOT 2 (which would be the domain.degree() fallback).
    {
        gsVector<index_t> nn =
            gsQuadrature::numNodes(domain, 1.0, 1, /*fixDir=*/0, degs32);
        CHECK_EQUAL(1, nn[0]);
        CHECK_EQUAL(3, nn[1]);
    }

    // fixDir = 1: direction 1 fixed -> 1 node; direction 0 survives with its
    // OWN degree 3 -> 1*3+1+0.5 -> 4.
    {
        gsVector<index_t> nn =
            gsQuadrature::numNodes(domain, 1.0, 1, /*fixDir=*/1, degs32);
        CHECK_EQUAL(4, nn[0]);
        CHECK_EQUAL(1, nn[1]);
    }

    // Contrast: fixDir = 0 with empty degrees -> direction 1 falls back to
    // domain.degree(1) == 1 -> 1*1+1+0.5 -> 2 (not 3).
    {
        gsVector<index_t> nn =
            gsQuadrature::numNodes(domain, 1.0, 1, /*fixDir=*/0, gsVector<short_t>());
        CHECK_EQUAL(1, nn[0]);
        CHECK_EQUAL(2, nn[1]);
    }

    // Same, through getPtr -- both quA and quB keys must be present since
    // getPtr reads them with getReal("quA")/getInt("quB").
    gsOptionList opt;
    opt.addInt ("quRule", "", gsQuadrature::GaussLegendre);
    opt.addReal("quA",    "", 1.0);
    opt.addInt ("quB",    "", 1);

    CHECK_EQUAL(3, gsQuadrature::getPtr<real_t>(domain, opt, /*fixDir=*/0, degs32)->numNodes()); // 1 x 3
    CHECK_EQUAL(4, gsQuadrature::getPtr<real_t>(domain, opt, /*fixDir=*/1, degs32)->numNodes()); // 4 x 1
    CHECK_EQUAL(2, gsQuadrature::getPtr<real_t>(domain, opt, /*fixDir=*/0, gsVector<short_t>())->numNodes()); // 1 x 2 fallback
}

void testWork(const index_t nodes[], const size_t dim)
{
    gsVector<index_t> numNodes = gsAsConstVector<index_t>(nodes, dim);

    // ---------------- Lobatto

    gsGaussRule<real_t> legendreRule(numNodes);
    gsGaussRule<real_t> legendreRuleComp(numNodes, REAL_DIG);
    gsQuadRule<real_t> quadLeg_lookup = gsQuadrature::get<real_t>(gsQuadrature::GaussLegendre, numNodes, 0);
    gsQuadRule<real_t> quadLeg_compute = gsQuadrature::get<real_t>(gsQuadrature::GaussLegendre, numNodes, REAL_DIG);
    CHECK(legendreRule.referenceNodes() == quadLeg_lookup.referenceNodes());
    CHECK_MATRIX_CLOSE(quadLeg_compute.referenceNodes(), quadLeg_lookup.referenceNodes(), EPSILON);
    CHECK_MATRIX_CLOSE(quadLeg_compute.referenceWeights(), quadLeg_lookup.referenceWeights(), EPSILON);

    gsVector<index_t> legVec = (2 * numNodes - 1 * gsVector<index_t>::Ones(dim)).cwiseMax(gsVector<index_t>::Zero(dim));
    real_t expectedLeg = calcAntiDerivative(legVec, dim);
    real_t lookupLeg = calcPoly(legVec, quadLeg_lookup, dim);
    real_t computeLeg = calcPoly(legVec, quadLeg_compute, dim);
    CHECK_CLOSE(expectedLeg, lookupLeg, EPSILON);
    CHECK_CLOSE(expectedLeg, computeLeg, EPSILON);

    // ---------------- Lobatto

    gsLobattoRule<real_t> lobattoRule(numNodes);
    gsLobattoRule<real_t> lobattoRuleComp(numNodes, REAL_DIG);
    gsQuadRule<real_t> quadLob_lookup = gsQuadrature::get<real_t>(gsQuadrature::GaussLobatto, numNodes, 0);
    gsQuadRule<real_t> quadLob_compute = gsQuadrature::get<real_t>(gsQuadrature::GaussLobatto, numNodes, REAL_DIG);
    CHECK(lobattoRule.referenceNodes() == quadLob_lookup.referenceNodes());
    CHECK_MATRIX_CLOSE(quadLob_lookup.referenceNodes(), quadLob_compute.referenceNodes(), EPSILON);
    CHECK_MATRIX_CLOSE(quadLob_compute.referenceWeights(), quadLob_lookup.referenceWeights(), EPSILON);

    gsVector<index_t> lobVec = (2 * numNodes - 3 * gsVector<index_t>::Ones(dim)).cwiseMax(gsVector<index_t>::Zero(dim));
    real_t expectedLob = calcAntiDerivative(lobVec, dim);
    real_t lookupLob = calcPoly(lobVec, quadLob_lookup, dim);
    real_t computeLob = calcPoly(lobVec, quadLob_compute, dim);
    CHECK_CLOSE(expectedLob, lookupLob, EPSILON);
    CHECK_CLOSE(expectedLob, computeLob, EPSILON);
}

real_t calcAntiDerivative(gsVector<index_t> const &deg, const index_t dim)
{
    // Test index_tegration
    gsVector<real_t> u;
    u.setConstant(dim, 1.0123);

    // Construct polynomial
    std::string var = std::string("xyzwuv").substr(0, dim);    // defining variable names
    std::stringstream poly;
    for (index_t i = 0; i < dim; ++i)
    {
        // Construct the anti-derivative
        std::string tmp(var);
        tmp.erase(i, 1); // cut out current variable, that will be to the power of deg[i]
        for (index_t j = dim - 1; j; --j) tmp.insert(j, "*");   // add * between variables
        poly << "(1.0/" << deg[i] + 1 << ")*" << tmp << var[i] << "^" << deg[i] + 1 << addPlus(dim, i);
    }

    gsFunctionExpr<real_t> antideriv(poly.str(), dim);
    return antideriv.eval(u).at(0);
}

real_t calcPoly(gsVector<index_t> const &deg,
                gsQuadRule<real_t> const &gr,
                const index_t dim)
{
    // Test index_tegration
    gsVector<real_t> u;
    u.setConstant(dim, 1.0123);

    const index_t d = gr.dim();
    CHECK_EQUAL(d, dim);

    const gsVector<real_t> l = gsVector<real_t>::Zero(d);

    gsMatrix<real_t> ngrid;
    gsVector<real_t> wgrid;

    // Map rule to index_tegration domain
    gr.mapTo(l, u, ngrid, wgrid);

    // Construct polynomial
    std::string var = std::string("xyzwuv").substr(0, d);    // defining variable names
    std::stringstream poly;
    for (index_t i = 0; i < d; ++i)
    {
        // Make a polynomial of requested degree
        poly << var[i] << "^" << deg[i] << addPlus(d, i);
    }

    gsFunctionExpr<real_t> index_tegrand(poly.str(), d);

    ngrid = index_tegrand.eval(ngrid);
    return wgrid.dot(ngrid.row(0));
}

const char *addPlus(const index_t d, index_t i)
{ return (i == d - 1 ? "" : "+"); }

}
