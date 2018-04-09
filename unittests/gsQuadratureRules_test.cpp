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

void testWork(const int dim, const int nodes[]);

real_t calcAntiDerivative(gsVector<index_t> const &deg, const int dim);

real_t calcPoly(gsVector<index_t> const &deg,
                gsQuadRule<real_t> const &gr,
                const int dim);

gsVector<index_t> noneMinus(gsVector<index_t> inVec);

const char *addPlus(const int d, int i);

TEST(R3_1)
{
    int array[] = {1, 5, 9};
    testWork(3, array);
}

TEST(R3_2)
{
    int array[] = {2, 6, 7};
    testWork(3, array);
}

TEST(R3_3)
{
    int array[] = {3, 4, 8};
    testWork(3, array);
}

TEST(R2_1)
{
    int array[] = {10, 17};
    testWork(2, array);
}

TEST(R2_2)
{
    int array[] = {11, 16};
    testWork(2, array);
}

TEST(R2_3)
{
    int array[] = {12, 15};
    testWork(2, array);
}

TEST(R2_4)
{
    int array[] = {13, 14};
    testWork(2, array);
}

TEST(R1_1)
{
    int array[] = {18};
    testWork(1, array);
}

TEST(R1_2)
{
    int array[] = {19};
    testWork(1, array);
}

TEST(R1_3)
{
    int array[] = {20};
    testWork(1, array);
}

TEST(R1_4)
{
    int array[] = {21};
    testWork(1, array);
}

void testWork(const int dim, const int nodes[])
{
    std::vector<index_t> qNodes;    // row vector
    qNodes.reserve(dim);

    for (int i = 0; i < dim; ++i)
    {
        qNodes.push_back(nodes[i]);
    }

    // Dimension of the rule
    const int d = qNodes.size();
    CHECK_EQUAL(dim, d);

    // Number of quadrature points
    gsVector<index_t> numNodes = gsAsMatrix<index_t>(qNodes).transpose();   // column vector
    CHECK_EQUAL(dim, numNodes.size());

    // Setup the reference rule
    gsGaussRule<real_t> legendreRule = gsGaussRule<real_t>(numNodes);
    gsGaussRule<real_t> legendreRuleComp = gsGaussRule<real_t>(numNodes, REAL_DIG);

    gsQuadRule<real_t> quadLeg_lookup = gsQuadrature::get<real_t>(gsQuadrature::GaussLegendre, numNodes, 0);
    gsQuadRule<real_t> quadLeg_compute = gsQuadrature::get<real_t>(gsQuadrature::GaussLegendre, numNodes, REAL_DIG);

    CHECK(legendreRule.referenceNodes() == quadLeg_lookup.referenceNodes());
    CHECK(legendreRuleComp.referenceNodes() == quadLeg_compute.referenceNodes());

    gsLobattoRule<real_t> lobattoRule = gsLobattoRule<real_t>(numNodes);
    gsLobattoRule<real_t> lobattoRuleComp = gsLobattoRule<real_t>(numNodes, REAL_DIG);

    gsQuadRule<real_t> quadLob_lookup = gsQuadrature::get<real_t>(gsQuadrature::GaussLobatto, numNodes, 0);
    gsQuadRule<real_t> quadLob_compute = gsQuadrature::get<real_t>(gsQuadrature::GaussLobatto, numNodes, REAL_DIG);

    CHECK(lobattoRule.referenceNodes() == quadLob_lookup.referenceNodes());
    CHECK(lobattoRuleComp.referenceNodes() == quadLob_compute.referenceNodes());

    gsVector<index_t> legVec = noneMinus(2 * numNodes - 1 * gsVector<index_t>::Ones(d));
    gsVector<index_t> lobVec = noneMinus(2 * numNodes - 3 * gsVector<index_t>::Ones(d));

    real_t expectedLeg = calcAntiDerivative(legVec, dim);
    real_t lookupLeg = calcPoly(legVec, quadLeg_lookup, dim);
    real_t computeLeg = calcPoly(legVec, quadLeg_compute, dim);

    real_t expectedLob = calcAntiDerivative(lobVec, dim);
    real_t lookupLob = calcPoly(lobVec, quadLob_lookup, dim);
    real_t computeLob = calcPoly(lobVec, quadLob_compute, dim);

    CHECK_CLOSE(expectedLeg, lookupLeg, EPSILON);
    CHECK_CLOSE(expectedLeg, computeLeg, EPSILON);
    //CHECK_CLOSE(computeLeg,  lookupLeg,  EPSILON);

    CHECK_CLOSE(expectedLob, lookupLob, EPSILON);
    CHECK_CLOSE(expectedLob, computeLob, EPSILON);
    //CHECK_CLOSE(computeLob,  lookupLob, EPSILON);
}

real_t calcAntiDerivative(gsVector<index_t> const &deg, const int dim)
{
    // Test integration
    gsVector<real_t> u;
    u.setConstant(dim, 1.0123);

    // Construct both polynoms (integrand and anti-derivative)
    std::string var = std::string("xyzwuv").substr(0, dim);    // defining variable names
    std::stringstream poly;
    for (int i = 0; i < dim; ++i)
    {
        // Compute the anti-derivative
        std::string tmp(var);
        tmp.erase(i, 1); // cut out current variable, that will be to the power of deg[i]
        for (int j = dim - 1; j; --j) tmp.insert(j, "*");   // add * between variables
        poly << "(1.0/" << deg[i] + 1 << ")*" << tmp << var[i] << "^" << deg[i] + 1 << addPlus(dim, i);
    }

    gsFunctionExpr<real_t> antideriv(poly.str(), dim);
    return antideriv.eval(u).at(0);
}

real_t calcPoly(gsVector<index_t> const &deg,
                gsQuadRule<real_t> const &gr,
                const int dim)
{
    // Test integration
    gsVector<real_t> u;
    u.setConstant(dim, 1.0123);

    const int d = gr.dim();
    CHECK_EQUAL(d, dim);

    const gsVector<real_t> l = gsVector<real_t>::Zero(d);

    gsMatrix<real_t> ngrid;
    gsVector<real_t> wgrid;

    // Map rule to integration domain
    gr.mapTo(l, u, ngrid, wgrid);

    // Construct polynomial
    std::string var = std::string("xyzwuv").substr(0, d);    // defining variable names
    std::stringstream poly;
    for (int i = 0; i < d; ++i)
    {
        // Make a polynomial of requested degree
        poly << var[i] << "^" << deg[i] << addPlus(d, i);
    }

    gsFunctionExpr<real_t> integrand(poly.str(), d);

    ngrid = integrand.eval(ngrid);
    return wgrid.dot(ngrid.row(0));
}

gsVector<index_t> noneMinus(gsVector<index_t> inVec)
{
    const int size = inVec.size();
    for (int i = 0; i < size; ++i)
    {
        if (inVec[i] < 0)
            inVec[i] = 0;
    }
    return inVec;
}

const char *addPlus(const int d, int i)
{ return (i == d - 1 ? "" : "+"); }

}
