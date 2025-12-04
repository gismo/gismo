/** @file solveExprUsage.cpp

    @brief Example showing how to use `gismo::solve` from `gsAutoDiffUtils.h`

    Demonstrates creating expression pointers for matrix `A` and vector `b`,
    calling the convenience `solve(a_expr,b_expr,solver)` to produce a
    `SolveExpr` instance and explains how the adjoint (backward) step is
    performed via `solver.solve_transpose(...)` (as used in `SolveExpr::propagate`).

*/

#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsAutoDiff/gsAutoDiffUtils.h>

using namespace gismo;
using namespace autodiff::reverse::detail; // for ExprPtr and constant()

// Simple dense solver wrapper that holds A and implements solve and solve_transpose
struct DenseSolver
{
    gsMatrix<double> A;
    DenseSolver() {}
    DenseSolver(const gsMatrix<double> & A_) : A(A_) {}

    // Forward solve: x = A^{-1} b
    gsMatrix<double> solve(const gsMatrix<double> & b) const
    {
        return A.fullPivLu().solve(b);
    }

    // Solve transposed system: solves A^T y = wprime
    gsMatrix<double> solve_transpose(const gsMatrix<double> & wprime) const
    {
        return A.transpose().fullPivLu().solve(wprime);
    }
};

int main()
{
    std::cout << "=== solveExprUsage example ===" << std::endl;

    // Build numeric A and b
    gsMatrix<double> A(2,2);
    A << 2.0, 1.0,
         1.0, 3.0;
    gsMatrix<double> b(2,1);
    b << 5.0,
         7.0;

    // Create ExprPtr constants (these are leaf expressions in the autodiff graph)
    auto a_expr = constant<gsMatrix<double>>(A);
    auto b_expr = constant<gsMatrix<double>>(b);

    // Create solver initialized with A's numeric value
    DenseSolver solver(A);

    // Use the convenience solve function from gsAutoDiffUtils.h
    // This returns an ExprPtr<gsMatrix<double>> whose concrete type is SolveExpr<gsMatrix<double>,DenseSolver>
    auto x_expr = gismo::solve<gsMatrix<double>, DenseSolver>(a_expr, b_expr, solver);

    // Forward value: x_expr->val
    std::cout << "Forward: x = A^{-1} b =\n" << x_expr->val << std::endl;

    // Suppose we have a scalar loss L that depends on x; its upstream gradient w.r.t. x is dL/dx
    gsMatrix<double> dL_dx(2,1);
    dL_dx << 1.0, 2.0; // example upstream gradient

    // In SolveExpr::propagate the following happens:
    //   adjoint_b = solver.solve_transpose(wprime);
    //   r->propagate(adjoint_b);
    // So the solver must be able to compute A^T \ wprime.
    gsMatrix<double> adjoint_b = solver.solve_transpose(dL_dx);
    std::cout << "Simulated backward: adjoint on b (A^T \\ dL/dx) =\n" << adjoint_b << std::endl;

    // If A were variable, its adjoint contribution would be:
    //   dL/dA = - adjoint_b * x^T  (outer product)
    gsMatrix<double> adj_A = - (adjoint_b * x_expr->val.transpose());
    std::cout << "Adjoint w.r.t. A (entrywise) =\n" << adj_A << std::endl;

    std::cout << "\nNotes:\n"
              << " - `gismo::solve` constructs a SolveExpr node that uses the provided solver.\n"
              << " - During reverse-mode, SolveExpr::propagate calls `solver.solve_transpose(wprime)`.\n"
              << " - The example simulates the backward step by calling `solver.solve_transpose` directly.\n";

    return 0;
}
