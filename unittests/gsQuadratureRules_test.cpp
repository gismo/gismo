/** @file gsQuadratureRules_test.cpp

    @brief Tutorial for setting up unittests

    This is an example unit test that doesn't really do anything useful.
    It is here as a reference for you when creating additional unit tests.
    For additional reference information, see the "test.h" header.

    == BASIC REFERENCE ==
         - TEST(NAME_OF_TEST) { body_of_test }
         - TEST_FIXTURE(NAME_OF_FIXTURE,NAME_OF_TEST){ body_of_test }

    == CHECK MACRO REFERENCE ==
         - CHECK(EXPR);
         - CHECK_EQUAL(EXPECTED,ACTUAL);
         - CHECK_CLOSE(EXPECTED,ACTUAL,EPSILON);
         - CHECK_ARRAY_EQUAL(EXPECTED,ACTUAL,LENGTH);
         - CHECK_ARRAY_CLOSE(EXPECTED,ACTUAL,LENGTH,EPSILON);
         - CHECK_ARRAY2D_EQUAL(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT);
         - CHECK_ARRAY2D_CLOSE(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT,EPSILON);
         - CHECK_THROW(EXPR,EXCEPTION_TYPE_EXPECTED);

    == TIME CONSTRAINTS ==
         - UNITTEST_TIME_CONSTRAINT(TIME_IN_MILLISECONDS);
         - UNITTEST_TIME_CONSTRAINT_EXEMPT();

    == MORE INFO ==
         See: https://unittest-cpp.github.io/

    Author(s): J. Vogl, A. Mantzaflaris
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework

SUITE(gsQuadratureRules_test)                 // The suite should have the same name as the file
{

void testWork(const int dim, const int nodes[], const double epsilon);

void testPoly(gsVector<index_t> const &deg,
              gsQuadRule<real_t> const &gr,
              gsVector<real_t> const &u,
              const int dim,
              const double epsilon);

const char *addPlus(const int d, int i);

TEST(R3_1)
{
    int array[] = {1, 5, 9};
    testWork(3, array, 1e-16);
}

TEST(R3_2)
{
    int array[] = {2, 6, 7};
    testWork(3, array, 1e-16);
}

TEST(R3_3)
{
    int array[] = {3, 4, 8};
    testWork(3, array, 2e-16);
}

TEST(R2_1)
{
    int array[] = {10, 17};
    testWork(2, array, 1e-16);
}

TEST(R2_2)
{
    int array[] = {11, 16};
    testWork(2, array, 3e-16);
}

TEST(R2_3)
{
    int array[] = {12, 15};
    testWork(2, array, 1e-16);
}

TEST(R2_4)
{
    int array[] = {13, 14};
    testWork(2, array, 1e-16);
}

TEST(R1_1)
{
    int array[] = {18};
    testWork(1, array, 2e-16);
}

TEST(R1_2)
{
    int array[] = {19};
    testWork(1, array, 1e-16);
}

TEST(R1_3)
{
    int array[] = {20};
    testWork(1, array, 1e-16);
}

TEST(R1_4)
{
    int array[] = {21};
    testWork(1, array, 1e-16);
}

void testWork(const int dim, const int nodes[], const double epsilon)
{
    std::vector<index_t> qNodes;    // row vector
    qNodes.reserve(dim);

    for (int i = 0; i < dim; ++i)
    {
        qNodes.push_back(nodes[i]);
    }

    // Dimension of the rule
    int d = qNodes.size();
    CHECK_EQUAL(dim, d);

    // Number of quadrature points
    gsVector<index_t> numNodes = gsAsMatrix<index_t>(qNodes).transpose();   // column vector
    CHECK_EQUAL(dim, numNodes.size());

    // Setup the reference rule
    gsGaussRule<real_t> GaussRule;
    GaussRule.setNodes(numNodes);

    // Test integration
    gsVector<real_t> u;
    u.setConstant(d, 1.0123);


    testPoly(2 * numNodes - gsVector<index_t>::Ones(d), GaussRule, u, dim, epsilon);
}

void testPoly(gsVector<index_t> const &deg,
              gsQuadRule<real_t> const &gr,
              gsVector<real_t> const &u, const int dim, const double epsilon)
{
    const int d = gr.dim();
    CHECK_EQUAL(d, dim);

    const gsVector<real_t> l = gsVector<real_t>::Zero(d);

    gsMatrix<real_t> ngrid;
    gsVector<real_t> wgrid;

    // Map rule to integration domain
    gr.mapTo(l, u, ngrid, wgrid);

    // Construct both polynoms (integrand and anti-derivative)
    std::string var = std::string("xyzwuv").substr(0, d);    // defining variable names
    std::stringstream poly1, poly2;
    for (int i = 0; i < d; ++i)
    {
        // Make a polynomial of requested degree
        poly1 << var[i] << "^" << deg[i] << addPlus(d, i);

        // Compute the anti-derivative
        std::string tmp(var);
        tmp.erase(i, 1); // cut out current variable, that will be to the power of deg[i]
        for (int j = d - 1; j; --j) tmp.insert(j, "*");   // add * between variables
        poly2 << "(1.0/" << deg[i] + 1 << ")*" << tmp << var[i] << "^" << deg[i] + 1 << addPlus(d, i);
    }

    gsFunctionExpr<real_t> integrand(poly1.str(), d);
    gsFunctionExpr<real_t> antideriv(poly2.str(), d);

    ngrid = integrand.eval(ngrid);
    real_t quadRes = wgrid.dot(ngrid.row(0));
    real_t aderRes = antideriv.eval(u).at(0);

    CHECK_CLOSE(aderRes, quadRes, epsilon);
}

const char *addPlus(const int d, int i)
{ return (i == d - 1 ? "" : "+"); }

}
