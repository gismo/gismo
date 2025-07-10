/** @file gsExpressions_test.cpp

    @brief Testing integral computation using the expression evaluator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

#include <gsNewExpressions/NewExpressions.h>
using namespace gismo;
using namespace Expr;







int main(int argc, char *argv[])
{

    bool verbose = false;
    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    ScalarExpression<real_t> scalar_expr(5.0);
    gsVector<real_t> vector_data = gsVector<real_t>::LinSpaced(3, 1.0, 3.0);
    VectorExpression<real_t> vector_expr(vector_data);
    gsMatrix<real_t> matrix_data(2, 2);
    matrix_data << 1.0, 2.0,
                   3.0, 4.0;
    MatrixExpression<real_t> matrix_expr(matrix_data);

    // Test addition of scalar + scalar
    auto scalar_plus_scalar = scalar_expr + scalar_expr;
    gsMatrix<real_t> result_scalar = scalar_plus_scalar.eval(0);
    std::cout << "Result of scalar + scalar:\n" << result_scalar << std::endl;

    // // Test addition of scalar + vector
    // auto scalar_plus_vector = scalar_expr + vector_expr;
    // gsMatrix<real_t> result_vector = scalar_plus_vector.eval(0);
    // std::cout << "Result of scalar + vector:\n" << result_vector << std::endl;

    // // Test addition of vector + scalar
    // auto vector_plus_scalar = vector_expr + scalar_expr;
    // gsMatrix<real_t> result_vector2 = vector_plus_scalar.eval(0);
    // std::cout << "Result of vector + scalar:\n" << result_vector2 << std::endl;

    // Test addition of two vectors
    auto vector_plus_vector = vector_expr + vector_expr;
    gsMatrix<real_t> result_vector3 = vector_plus_vector.eval(0);
    std::cout << "Result of vector + vector:\n" << result_vector3 << std::endl;

    // // Test addition of scalar + matrix
    // auto scalar_plus_matrix = scalar_expr + matrix_expr;
    // gsMatrix<real_t> result_matrix = scalar_plus_matrix.eval(0);
    // std::cout << "Result of scalar + matrix:\n" << result_matrix << std::endl;

    // // Test addition of matrix + scalar
    // auto matrix_plus_scalar = matrix_expr + scalar_expr;
    // gsMatrix<real_t> result_matrix2 = matrix_plus_scalar.eval(0);
    // std::cout << "Result of matrix + scalar:\n" << result_matrix2 << std::endl;

    return EXIT_SUCCESS;
}
