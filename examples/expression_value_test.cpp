/** @file expression_value_test.cpp

    @brief Testing the ExpressionResult container

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gismo.h>
#include <gsNewExpressions/ExpressionUtils.h>
#include <gsNewExpressions/ExpressionResult.h>

using namespace gismo;
using namespace Expr;

int main(int argc, char *argv[])
{
    gsCmdLine cmd("Testing ExpressionResult container.");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "=== Testing ExpressionResult ===\n\n";

    // Test 1: Scalar expression (1,1)
    gsInfo << "Test 1: Scalar expression (cardinality 1,1)\n";
    ExpressionResult<real_t> scalar_val = makeExpressionResult<real_t>(SpaceType::None);
    scalar_val().resize(2, 3);
    scalar_val().setConstant(5.0);
    gsInfo << "  Cardinality: (" << scalar_val.rowCardinality() << ", " << scalar_val.colCardinality() << ")\n";
    gsInfo << "  Matrix(0,0):\n" << scalar_val(0, 0) << "\n\n";
    gsInfo << "  Matrix() [shorthand]:\n" << scalar_val() << "\n\n";

    // Test 2: Test space expression (N,1)
    gsInfo << "Test 2: Test space expression (cardinality N,1) with N=3\n";
    index_t N = 3;
    ExpressionResult<real_t> test_val = makeExpressionResult<real_t>(SpaceType::Test, N, 1, 2, 2);
    for (index_t i = 0; i < N; ++i)
    {
        test_val(i, 0).setConstant(i + 1.0);
    }
    gsInfo << "  Cardinality: (" << test_val.rowCardinality() << ", " << test_val.colCardinality() << ")\n";
    for (index_t i = 0; i < N; ++i)
    {
        gsInfo << "  Matrix(" << i << ",0):\n" << test_val(i, 0) << "\n";
    }
    gsInfo << "\n";

    // Test 3: Trial space expression (1,M)
    gsInfo << "Test 3: Trial space expression (cardinality 1,M) with M=4\n";
    index_t M = 4;
    ExpressionResult<real_t> trial_val = makeExpressionResult<real_t>(SpaceType::Trial, 1, M, 3, 1);
    for (index_t j = 0; j < M; ++j)
    {
        trial_val(0, j).setConstant(10.0 * (j + 1.0));
    }
    gsInfo << "  Cardinality: (" << trial_val.rowCardinality() << ", " << trial_val.colCardinality() << ")\n";
    for (index_t j = 0; j < M; ++j)
    {
        gsInfo << "  Matrix(0," << j << "):\n" << trial_val(0, j).transpose() << "\n";
    }
    gsInfo << "\n";

    // Test 4: Bilinear expression (N,M)
    gsInfo << "Test 4: Bilinear expression (cardinality N,M) with N=2, M=3\n";
    N = 2;
    M = 3;
    ExpressionResult<real_t> bilinear_val = makeExpressionResult<real_t>(SpaceType::Both, N, M, 1, 1);
    for (index_t i = 0; i < N; ++i)
    {
        for (index_t j = 0; j < M; ++j)
        {
            bilinear_val(i, j)(0, 0) = 100.0 * (i + 1) + (j + 1);
        }
    }
    gsInfo << "  Cardinality: (" << bilinear_val.rowCardinality() << ", " << bilinear_val.colCardinality() << ")\n";
    gsInfo << "  Matrix values:\n";
    for (index_t i = 0; i < N; ++i)
    {
        for (index_t j = 0; j < M; ++j)
        {
            gsInfo << "    (" << i << "," << j << "): " << bilinear_val(i, j)(0, 0) << "\n";
        }
    }
    gsInfo << "\n";

    // Test 5: Iteration
    gsInfo << "Test 5: Iterating over all matrices\n";
    ExpressionResult<real_t> test_iter(2, 2, 1, 1);
    index_t counter = 0;
    for (auto& mat : test_iter)
    {
        mat(0, 0) = counter++;
    }
    gsInfo << "  Using range-based for loop:\n";
    counter = 0;
    for (const auto& mat : test_iter)
    {
        gsInfo << "    Matrix[" << counter++ << "]: " << mat(0, 0) << "\n";
    }
    gsInfo << "\n";

    // Test 6: Resize
    gsInfo << "Test 6: Resizing cardinality\n";
    ExpressionResult<real_t> resize_val(2, 2, 1, 1);
    resize_val(0, 0)(0, 0) = 1.0;
    resize_val(1, 1)(0, 0) = 2.0;
    gsInfo << "  Original cardinality: (" << resize_val.rowCardinality() << ", " << resize_val.colCardinality() << ")\n";
    gsInfo << "  Matrix(0,0): " << resize_val(0, 0)(0, 0) << "\n";
    gsInfo << "  Matrix(1,1): " << resize_val(1, 1)(0, 0) << "\n";
    
    resize_val.resize(3, 3, true);
    gsInfo << "  After resize to (3,3) with preserve=true:\n";
    gsInfo << "    New cardinality: (" << resize_val.rowCardinality() << ", " << resize_val.colCardinality() << ")\n";
    gsInfo << "    Matrix(0,0): " << resize_val(0, 0)(0, 0) << "\n";
    gsInfo << "    Matrix(1,1): " << resize_val(1, 1)(0, 0) << "\n";
    gsInfo << "\n";

    gsInfo << "=== All tests completed successfully! ===\n";

    return EXIT_SUCCESS;
}
