/** @file linearSolveAdjointExample.cpp

    @brief Example showing adjoint computation for a linear solve

    Demonstrates how to compute the adjoint (sensitivity) of a scalar loss
    that depends on the solution of a linear system x = A^{-1} b by solving
    the transposed linear system A^T * lambda = dL/dx (i.e. using solve_transpose).

    This mirrors the approach used in `codipack_example.cpp`'s `solveSystem_rev`.

*/

#include <gismo.h>

using namespace gismo;

int main()
{
    std::cout << "=== linearSolveAdjointExample ===" << std::endl;

    // Small 2x2 example
    gsMatrix<double> A(2,2);
    A << 2.0, 1.0,
         1.0, 3.0;
    gsMatrix<double> b(2,1);
    b << 5.0,
         7.0;

    // Solve Ax = b
    gsMatrix<double> x = A.fullPivLu().solve(b);
    std::cout << "Solution x:\n" << x << std::endl;

    // Define a simple scalar loss L = 0.5 * ||x - x_target||^2
    gsMatrix<double> x_target(2,1);
    x_target << 1.0, 2.0;
    gsMatrix<double> diff = x - x_target;
    double L = 0.5 * diff.squaredNorm();
    std::cout << "Loss L = " << L << std::endl;

    // dL/dx = x - x_target
    gsMatrix<double> dL_dx = diff;

    // To compute adjoint on b: solve A^T * lambda = dL/dx
    // This is equivalent to calling solver.solve_transpose(dL_dx) in a solver wrapper.
    gsMatrix<double> lambda;
    {
        // use LU on A^T
        auto solver = A.transpose().fullPivLu();
        lambda = solver.solve(dL_dx);
    }

    std::cout << "Adjoint (lambda) solving A^T * lambda = dL/dx:\n" << lambda << std::endl;

    // Adjoint on b is lambda (since x = A^{-1} b -> dx/db = A^{-1})
    gsMatrix<double> adj_b = lambda;
    std::cout << "Adjoint w.r.t. b (dL/db):\n" << adj_b << std::endl;

    // Adjoint on A: for small perturbation dA, dx = -A^{-1} dA x
    // so dL/dA = - lambda * x^T (outer product)
    gsMatrix<double> adj_A = - (lambda * x.transpose());
    std::cout << "Adjoint w.r.t. A (dL/dA) as matrix:\n" << adj_A << std::endl;

    // Flattened view (entrywise) if needed
    std::cout << "Adjoint entries (dL/dA) flattened:\n";
    for (index_t i=0;i<adj_A.rows();++i)
        for (index_t j=0;j<adj_A.cols();++j)
            std::cout << "dL/dA("<<i<<","<<j<<") = "<< adj_A(i,j) << std::endl;

    return 0;
}
