/** @file gsQuadratureRules_test.cpp

    @brief Tests available gsQuadratureRules (containing gsGaussRules 1-21)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Vogl, A. Mantzaflaris
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework

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

// Integrate f element-by-element over domain using gsPatchRule
real_t patchIntegrate(const gsDomain<real_t>& domain,
                      gsPatchRule<real_t>& rule,
                      const gsFunctionExpr<real_t>& f)
{
    gsMatrix<real_t> nodes;
    gsVector<real_t> weights;
    real_t result = 0;
    auto it    = domain.beginAll();
    auto itEnd = domain.endAll();
    for (; it < itEnd; ++it)
    {
        rule.mapTo(it.lowerCorner(), it.upperCorner(), nodes, weights);
        if (weights.size() == 0) continue;
        gsMatrix<real_t> vals = f.eval(nodes);
        result += weights.dot(vals.row(0).transpose());
    }
    return result;
}

// Integrate f element-by-element over domain using a Gauss-Legendre rule
real_t gaussIntegrate(const gsDomain<real_t>& domain,
                      gsGaussRule<real_t>& rule,
                      const gsFunctionExpr<real_t>& f)
{
    gsMatrix<real_t> nodes;
    gsVector<real_t> weights;
    real_t result = 0;
    auto it    = domain.beginAll();
    auto itEnd = domain.endAll();
    for (; it < itEnd; ++it)
    {
        rule.mapTo(it.lowerCorner(), it.upperCorner(), nodes, weights);
        gsMatrix<real_t> vals = f.eval(nodes);
        result += weights.dot(vals.row(0).transpose());
    }
    return result;
}

TEST(patchrule_1d_weight_sum)
{
    gsKnotVector<real_t> kv(0, 1, 3, 3, 1); // degree 3, 4 elements
    gsTensorBSplineBasis<1,real_t> basis(kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);
    gsFunctionExpr<real_t> f("1", 1);
    CHECK_CLOSE(1.0, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_1d_polynomial_exactness)
{
    gsKnotVector<real_t> kv(0, 1, 3, 3, 1);
    gsTensorBSplineBasis<1,real_t> basis(kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);
    gsFunctionExpr<real_t> f("x^3", 1);
    CHECK_CLOSE(0.25, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_2d_weight_sum)
{
    gsKnotVector<real_t> kv(0, 1, 3, 3, 1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);
    gsFunctionExpr<real_t> f("1", 2);
    CHECK_CLOSE(1.0, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_2d_polynomial_exactness)
{
    gsKnotVector<real_t> kv(0, 1, 3, 3, 1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);
    gsFunctionExpr<real_t> f("x^2*y^2", 2);
    // integral of x^2*y^2 over [0,1]^2 = (1/3)*(1/3) = 1/9
    CHECK_CLOSE(1.0/9.0, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_3d_weight_sum)
{
    gsKnotVector<real_t> kv(0, 1, 2, 3, 1); // degree 3, 3 elements
    gsTensorBSplineBasis<3,real_t> basis(kv, kv, kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);
    gsFunctionExpr<real_t> f("1", 3);
    CHECK_CLOSE(1.0, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_3d_polynomial_exactness)
{
    gsKnotVector<real_t> kv(0, 1, 2, 3, 1);
    gsTensorBSplineBasis<3,real_t> basis(kv, kv, kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);
    gsFunctionExpr<real_t> f("x^2*y^2*z^2", 3);
    // integral of x^2*y^2*z^2 over [0,1]^3 = (1/3)^3 = 1/27
    CHECK_CLOSE(1.0/27.0, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_overintegration_2d)
{
    gsKnotVector<real_t> kv(0, 1, 3, 3, 1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, true);
    gsFunctionExpr<real_t> f("x^3*y^3", 2);
    // integral of x^3*y^3 over [0,1]^2 = (1/4)*(1/4) = 1/16
    CHECK_CLOSE(1.0/16.0, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_matches_gauss_on_polynomial)
{
    gsKnotVector<real_t> kv(0, 1, 3, 3, 1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);

    // Gauss rule with 4 pts/direction is exact for degree 7, well above degree 3
    gsVector<index_t> numPts(2); numPts << 4, 4;
    gsGaussRule<real_t> gauss(numPts);

    gsFunctionExpr<real_t> f("x^3*y", 2);
    // analytic: (1/4)*(1/2) = 1/8
    real_t patchResult = patchIntegrate(*basis.domain(), rule, f);
    real_t gaussResult = gaussIntegrate(*basis.domain(), gauss, f);
    CHECK_CLOSE(1.0/8.0, patchResult, EPSILON);
    CHECK_CLOSE(gaussResult, patchResult, EPSILON);
}

} // SUITE gsQuadratureRules_test
