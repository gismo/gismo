/** @file linearSolveAutoDiff.cpp

    @brief Example demonstrating automatic differentiation through linear solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst

*/

//! [Include namespace]
#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>

using namespace gismo;
//! [Include namespace]

using namespace autodiff;

int main()
{
    std::cout << "=== Linear Solve with Automatic Differentiation ===" << std::endl;
    std::cout << "Solving Ax = b and computing derivatives w.r.t. parameters\n" << std::endl;
    
    // ========================================================================
    // Example 1: LHS depends on parameter, RHS does not
    // ========================================================================
    std::cout << "=== Example 1: LHS depends on parameter alpha, RHS is constant ===" << std::endl;
    std::cout << "Solving: A(alpha) * x = b, where A(alpha) = [[2+alpha, 1], [1, 3]]" << std::endl;
    
    dual alpha = 1.0;
    autodiff::detail::seed<1>(alpha, 1.0);  // Seed alpha for differentiation
    
    // Construct LHS matrix that depends on alpha
    gsMatrix<dual> A1(2, 2);
    A1 << 2.0 + alpha, 1.0,
          1.0,         3.0;
    
    // Constant RHS vector
    gsMatrix<dual> b1(2, 1);
    b1 << 5.0,
          7.0;
    
    // Solve the linear system
    gsMatrix<dual> x1 = A1.fullPivLu().solve(b1);
    
    // Extract values and derivatives
    gsMatrix<double> x1_val(2, 1);
    gsMatrix<double> dx1_dalpha(2, 1);
    for (index_t i = 0; i < 2; ++i)
    {
        x1_val(i, 0) = autodiff::detail::derivative<0>(x1(i, 0));
        dx1_dalpha(i, 0) = autodiff::detail::derivative<1>(x1(i, 0));
    }
    
    std::cout << "Solution x = \n" << x1_val << std::endl;
    std::cout << "Derivative dx/d(alpha) = \n" << dx1_dalpha << std::endl;
    
    // ========================================================================
    // Example 2: RHS depends on parameter, LHS does not
    // ========================================================================
    std::cout << "\n=== Example 2: RHS depends on parameter beta, LHS is constant ===" << std::endl;
    std::cout << "Solving: A * x = b(beta), where b(beta) = [5+beta, 7]" << std::endl;
    
    dual beta = 2.0;
    autodiff::detail::seed<1>(beta, 1.0);  // Seed beta for differentiation
    
    // Constant LHS matrix
    gsMatrix<dual> A2(2, 2);
    A2 << 3.0, 1.0,
          1.0, 2.0;
    
    // RHS vector that depends on beta
    gsMatrix<dual> b2(2, 1);
    b2 << 5.0 + beta,
          7.0;
    
    // Solve the linear system
    gsMatrix<dual> x2 = A2.fullPivLu().solve(b2);
    
    // Extract values and derivatives
    gsMatrix<double> x2_val(2, 1);
    gsMatrix<double> dx2_dbeta(2, 1);
    for (index_t i = 0; i < 2; ++i)
    {
        x2_val(i, 0) = autodiff::detail::derivative<0>(x2(i, 0));
        dx2_dbeta(i, 0) = autodiff::detail::derivative<1>(x2(i, 0));
    }
    
    std::cout << "Solution x = \n" << x2_val << std::endl;
    std::cout << "Derivative dx/d(beta) = \n" << dx2_dbeta << std::endl;
    
    // ========================================================================
    // Example 3: Both LHS and RHS depend on parameter
    // ========================================================================
    std::cout << "\n=== Example 3: Both LHS and RHS depend on parameter gamma ===" << std::endl;
    std::cout << "Solving: A(gamma) * x = b(gamma)" << std::endl;
    std::cout << "where A(gamma) = [[2+gamma, 1], [1, 3+gamma]]" << std::endl;
    std::cout << "and   b(gamma) = [5+gamma, 7-gamma]" << std::endl;
    
    dual gamma = 0.5;
    autodiff::detail::seed<1>(gamma, 1.0);  // Seed gamma for differentiation
    
    // LHS matrix that depends on gamma
    gsMatrix<dual> A3(2, 2);
    A3 << 2.0 + gamma, 1.0,
          1.0,         3.0 + gamma;
    
    // RHS vector that depends on gamma
    gsMatrix<dual> b3(2, 1);
    b3 << 5.0 + gamma,
          7.0 - gamma;
    
    // Solve the linear system
    gsMatrix<dual> x3 = A3.fullPivLu().solve(b3);
    
    // Extract values and derivatives
    gsMatrix<double> x3_val(2, 1);
    gsMatrix<double> dx3_dgamma(2, 1);
    for (index_t i = 0; i < 2; ++i)
    {
        x3_val(i, 0) = autodiff::detail::derivative<0>(x3(i, 0));
        dx3_dgamma(i, 0) = autodiff::detail::derivative<1>(x3(i, 0));
    }
    
    std::cout << "Solution x = \n" << x3_val << std::endl;
    std::cout << "Derivative dx/d(gamma) = \n" << dx3_dgamma << std::endl;
    
    // ========================================================================
    // Verification using implicit function theorem
    // ========================================================================
    std::cout << "\n=== Verification (Example 2) ===" << std::endl;
    std::cout << "Using implicit function theorem: dx/dbeta = A^(-1) * db/dbeta" << std::endl;
    
    // For Example 2: A * x = b(beta), where b(beta) = [5+beta, 7]
    gsMatrix<double> A2_check(2, 2);
    A2_check << 3.0, 1.0,
                1.0, 2.0;
    
    gsMatrix<double> db_dbeta(2, 1);
    db_dbeta << 1.0,  // d(5+beta)/dbeta = 1
                0.0;  // d(7)/dbeta = 0
    
    gsMatrix<double> dx_dbeta_check = A2_check.fullPivLu().solve(db_dbeta);
    
    std::cout << "Analytical dx/dbeta = \n" << dx_dbeta_check << std::endl;
    std::cout << "AutoDiff dx/dbeta   = \n" << dx2_dbeta << std::endl;
    std::cout << "Difference          = \n" << (dx_dbeta_check - dx2_dbeta).norm() << std::endl;
    
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "The automatic differentiation correctly propagates derivatives" << std::endl;
    std::cout << "through linear solvers using the implicit function theorem:" << std::endl;
    std::cout << "  If A(p)*x(p) = b(p), then dx/dp = -A^(-1) * (dA/dp * x - db/dp)" << std::endl;
    
    return 0;
}
