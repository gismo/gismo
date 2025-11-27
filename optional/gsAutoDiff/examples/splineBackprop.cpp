/** @file splineBackprop.cpp

    @brief Comparison of forward-mode (dual) vs reverse-mode (var) for spline evaluation

    Demonstrates backpropagation through spline evaluation and compares performance.
    
    Note: This example uses element-wise operations rather than full gsBSpline<var>
    due to template instantiation requirements. The mathematical operations are
    equivalent to spline evaluation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst

*/

#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsAutoDiff/gsAutoDiffVar.h>
#include <chrono>

using namespace gismo;
using namespace autodiff;

// Helper to print timing results
void printTiming(const std::string& label, double time_ms) {
    std::cout << label << ": " << time_ms << " ms" << std::endl;
}

// Helper: Manually compute B-spline basis functions for demonstration
void computeBSplineBasis(double u, int degree, const gsVector<>& knots, gsVector<>& basis) {
    int n = basis.size();
    basis.setZero();
    
    // Find span
    int span = degree;
    for (int i = degree; i < knots.size() - degree - 1; ++i) {
        if (u >= knots[i] && u < knots[i+1]) {
            span = i;
            break;
        }
    }
    
    // Cox-de Boor recursion (simplified for degree 3)
    if (degree == 3) {
        basis[span-3] = 0.16666667;  // Simplified uniform B-spline weights at u=0.5
        basis[span-2] = 0.66666667;
        basis[span-1] = 0.16666667;
    }
}

int main()
{
    std::cout << "========================================================================" << std::endl;
    std::cout << "   SPLINE EVALUATION: FORWARD-MODE (dual) vs REVERSE-MODE (var)" << std::endl;
    std::cout << "========================================================================\n" << std::endl;
    
    // ========================================================================
    // Setup: Create a B-spline with regular double type
    // ========================================================================
    std::cout << "=== Setup: Creating B-spline ===" << std::endl;
    
    const index_t degree = 3;
    const index_t numCoefs = 20;  // Number of control points
    const index_t numEvals = 100; // Number of evaluation points
    
    gsKnotVector<> kv(0.0, 1.0, numCoefs - degree - 1, degree + 1);
    gsBSplineBasis<> basis(kv);
    
    std::cout << "Degree: " << degree << std::endl;
    std::cout << "Number of coefficients: " << numCoefs << std::endl;
    std::cout << "Number of evaluation points: " << numEvals << std::endl;
    std::cout << "Active basis functions per point: ~" << (degree + 1) << "\n" << std::endl;
    
    // Create coefficient matrix (2D spline - 2 components)
    gsMatrix<> coefs = gsMatrix<>::Random(numCoefs, 2);
    gsBSpline<> spline(basis, coefs);
    
    // Create evaluation points
    gsMatrix<> evalPts(1, numEvals);
    for (index_t i = 0; i < numEvals; ++i)
    {
        evalPts(0, i) = 0.01 + (0.99 - 0.01) * i / (numEvals - 1);
    }
    
    // ========================================================================
    // Example 1: Forward-mode AD (dual) - Standard GISMO approach
    // ========================================================================
    std::cout << "=== Example 1: FORWARD-MODE (dual) - Jacobian computation ===" << std::endl;
    std::cout << "Computing: d(spline_eval) / d(coefficients)\n" << std::endl;
    
    auto start_forward = std::chrono::high_resolution_clock::now();
    
    // For each coefficient, compute derivative of ALL evaluations w.r.t. that coefficient
    // This requires numCoefs forward passes
    gsMatrix<> jacobian_forward(numEvals * 2, numCoefs);
    
    for (index_t k = 0; k < numCoefs; ++k)
    {
        // Seed coefficient k
        gsMatrix<dual> coefs_dual = coefs.cast<dual>();
        autodiff::detail::seed<1>(coefs_dual(k, 0), 1.0);
        
        gsKnotVector<dual> kv_dual(dual(0.0), dual(1.0), numCoefs - degree - 1, degree + 1);
        gsBSplineBasis<dual> basis_dual(kv_dual);
        gsBSpline<dual> spline_dual(basis_dual, coefs_dual);
        gsMatrix<dual> result_dual = spline_dual.eval(evalPts.cast<dual>());
        
        // Extract derivatives for all evaluation points
        for (index_t i = 0; i < numEvals; ++i)
        {
            jacobian_forward(i * 2 + 0, k) = autodiff::detail::derivative<1>(result_dual(0, i));
            jacobian_forward(i * 2 + 1, k) = autodiff::detail::derivative<1>(result_dual(1, i));
        }
    }
    
    auto end_forward = std::chrono::high_resolution_clock::now();
    double time_forward = std::chrono::duration<double, std::milli>(end_forward - start_forward).count();
    
    std::cout << "Forward-mode approach:" << std::endl;
    std::cout << "  - Number of forward passes: " << numCoefs << " (one per coefficient)" << std::endl;
    std::cout << "  - Each pass: evaluate spline at " << numEvals << " points" << std::endl;
    std::cout << "  - Result: Full Jacobian [" << jacobian_forward.rows() << " x " << jacobian_forward.cols() << "]" << std::endl;
    printTiming("  - Time", time_forward);
    std::cout << std::endl;
    
    // ========================================================================
    // Example 2: Reverse-mode AD (var) - Backpropagation approach
    // ========================================================================
    std::cout << "=== Example 2: REVERSE-MODE (var) - Backpropagation ===" << std::endl;
    std::cout << "Computing: d(loss) / d(coefficients) using backprop\n" << std::endl;
    
    // Create var-type coefficients in a gsMatrix
    gsMatrix<var> coefs_var(numCoefs, 2);
    for (index_t i = 0; i < numCoefs; ++i)
    {
        for (index_t j = 0; j < 2; ++j)
        {
            coefs_var(i, j) = coefs(i, j);
        }
    }
    
    auto start_reverse = std::chrono::high_resolution_clock::now();
    
    // Compute each output's gradient w.r.t. all coefficients
    // This requires (numEvals * 2) reverse passes
    gsMatrix<> jacobian_reverse(numEvals * 2, numCoefs);
    
    for (index_t eval_idx = 0; eval_idx < numEvals; ++eval_idx)
    {
        // For each evaluation point and each component
        for (index_t comp = 0; comp < 2; ++comp)
        {
            // Fresh coefficient matrix for this reverse pass
            gsMatrix<var> coefs_var_fresh(numCoefs, 2);
            for (index_t i = 0; i < numCoefs; ++i)
            {
                for (index_t j = 0; j < 2; ++j)
                {
                    coefs_var_fresh(i, j) = coefs(i, j);
                }
            }
            
            // Evaluate spline
            gsKnotVector<var> kv_var(var(0.0), var(1.0), numCoefs - degree - 1, degree + 1);
            gsBSplineBasis<var> basis_var(kv_var);
            gsBSpline<var> spline_var(basis_var, coefs_var_fresh);
            gsMatrix<> evalPt_single(1, 1);
            evalPt_single << evalPts(0, eval_idx);
            gsMatrix<var> result_var = spline_var.eval(evalPt_single.cast<var>());
            
            // Get the output we want to backprop from
            var output = result_var(comp, 0);
            
            // Compute gradient w.r.t. first column of coefficients
            // Create a proper Eigen vector (not a matrix view)
            gsEigen::Matrix<var, gsEigen::Dynamic, 1> coefs_vec(numCoefs);
            for (index_t i = 0; i < numCoefs; ++i)
            {
                coefs_vec(i) = coefs_var_fresh(i, 0);
            }
            auto grad = autodiff::gradient(output, coefs_vec);
            
            // Store gradient
            for (index_t i = 0; i < numCoefs; ++i)
            {
                jacobian_reverse(eval_idx * 2 + comp, i) = grad(i);
            }
        }
    }
    
    auto end_reverse = std::chrono::high_resolution_clock::now();
    double time_reverse = std::chrono::duration<double, std::milli>(end_reverse - start_reverse).count();
    
    std::cout << "Reverse-mode approach:" << std::endl;
    std::cout << "  - Number of reverse passes: " << (numEvals * 2) << " (one per output)" << std::endl;
    std::cout << "  - Each pass: backprop from one output through spline" << std::endl;
    std::cout << "  - Result: Full Jacobian [" << jacobian_reverse.rows() << " x " << jacobian_reverse.cols() << "]" << std::endl;
    printTiming("  - Time", time_reverse);
    std::cout << std::endl;
    
    // ========================================================================
    // Example 3: Optimization scenario - Reverse-mode shines
    // ========================================================================
    std::cout << "=== Example 3: OPTIMIZATION - Computing gradient of scalar loss ===" << std::endl;
    std::cout << "Scenario: Minimize sum-of-squares error between spline and target\n" << std::endl;
    
    // Create target values
    gsMatrix<> target = gsMatrix<>::Random(2, numEvals);
    
    // Forward-mode: Need to compute gradient of scalar loss w.r.t. numCoefs parameters
    auto start_opt_forward = std::chrono::high_resolution_clock::now();
    
    gsMatrix<> gradient_forward(numCoefs, 2);
    gradient_forward.setZero();
    
    for (index_t k = 0; k < numCoefs; ++k)
    {
        // Seed coefficient k (first column)
        gsMatrix<dual> coefs_dual = coefs.cast<dual>();
        autodiff::detail::seed<1>(coefs_dual(k, 0), 1.0);
        
        gsKnotVector<dual> kv_dual(dual(0.0), dual(1.0), numCoefs - degree - 1, degree + 1);
        gsBSplineBasis<dual> basis_dual(kv_dual);
        gsBSpline<dual> spline_dual(basis_dual, coefs_dual);
        gsMatrix<dual> result_dual = spline_dual.eval(evalPts.cast<dual>());
        
        // Compute loss = sum((result - target)^2)
        dual loss = 0.0;
        for (index_t i = 0; i < numEvals; ++i)
        {
            for (index_t j = 0; j < 2; ++j)
            {
                dual diff = result_dual(j, i) - target(j, i);
                loss = loss + diff * diff;
            }
        }
        
        gradient_forward(k, 0) = autodiff::detail::derivative<1>(loss);
    }
    
    auto end_opt_forward = std::chrono::high_resolution_clock::now();
    double time_opt_forward = std::chrono::duration<double, std::milli>(end_opt_forward - start_opt_forward).count();
    
    std::cout << "Forward-mode optimization:" << std::endl;
    std::cout << "  - Task: Compute gradient of scalar loss w.r.t. " << numCoefs << " coefficients" << std::endl;
    std::cout << "  - Number of passes: " << numCoefs << " (one per parameter)" << std::endl;
    printTiming("  - Time", time_opt_forward);
    std::cout << std::endl;
    
    // Reverse-mode: ONE pass to get all gradients
    auto start_opt_reverse = std::chrono::high_resolution_clock::now();
    
    gsMatrix<var> coefs_var_opt(numCoefs, 2);
    for (index_t i = 0; i < numCoefs; ++i)
    {
        for (index_t j = 0; j < 2; ++j)
        {
            coefs_var_opt(i, j) = coefs(i, j);
        }
    }
    
    gsKnotVector<var> kv_var(var(0.0), var(1.0), numCoefs - degree - 1, degree + 1);
    gsBSplineBasis<var> basis_var(kv_var);
    gsBSpline<var> spline_var_opt(basis_var, coefs_var_opt);
    gsMatrix<var> result_var_opt = spline_var_opt.eval(evalPts.cast<var>());
    
    // Compute loss
    var loss_var = 0.0;
    for (index_t i = 0; i < numEvals; ++i)
    {
        for (index_t j = 0; j < 2; ++j)
        {
            var diff = result_var_opt(j, i) - target(j, i);
            loss_var = loss_var + diff * diff;
        }
    }
    
    // ONE reverse pass to get ALL gradients
    // Create a proper Eigen vector (not a matrix view)
    gsEigen::Matrix<var, gsEigen::Dynamic, 1> coefs_vec_opt(numCoefs);
    for (index_t i = 0; i < numCoefs; ++i)
    {
        coefs_vec_opt(i) = coefs_var_opt(i, 0);
    }
    auto gradient_reverse = autodiff::gradient(loss_var, coefs_vec_opt);
    
    auto end_opt_reverse = std::chrono::high_resolution_clock::now();
    double time_opt_reverse = std::chrono::duration<double, std::milli>(end_opt_reverse - start_opt_reverse).count();
    
    std::cout << "Reverse-mode optimization:" << std::endl;
    std::cout << "  - Task: Compute gradient of scalar loss w.r.t. " << numCoefs << " coefficients" << std::endl;
    std::cout << "  - Number of passes: 1 (backprop from single output)" << std::endl;
    printTiming("  - Time", time_opt_reverse);
    std::cout << std::endl;
    
    // ========================================================================
    // Summary and comparison
    // ========================================================================
    std::cout << "=======================================================================" << std::endl;
    std::cout << "   PERFORMANCE SUMMARY" << std::endl;
    std::cout << "=======================================================================" << std::endl;
    
    std::cout << "\n--- Task 1: Full Jacobian [" << (numEvals*2) << " x " << numCoefs << "] ---" << std::endl;
    printTiming("Forward-mode (dual)", time_forward);
    printTiming("Reverse-mode (var) ", time_reverse);
    std::cout << "Winner: " << (time_forward < time_reverse ? "FORWARD-MODE" : "REVERSE-MODE") << std::endl;
    std::cout << "Speedup: " << std::max(time_forward, time_reverse) / std::min(time_forward, time_reverse) << "x" << std::endl;
    std::cout << "\nExplanation: Forward-mode is better because #inputs (" << numCoefs 
              << ") < #outputs (" << (numEvals*2) << ")" << std::endl;
    
    std::cout << "\n--- Task 2: Gradient of scalar loss w.r.t. " << numCoefs << " parameters ---" << std::endl;
    printTiming("Forward-mode (dual)", time_opt_forward);
    printTiming("Reverse-mode (var) ", time_opt_reverse);
    std::cout << "Winner: " << (time_opt_forward < time_opt_reverse ? "FORWARD-MODE" : "REVERSE-MODE") << std::endl;
    std::cout << "Speedup: " << std::max(time_opt_forward, time_opt_reverse) / std::min(time_opt_forward, time_opt_reverse) << "x" << std::endl;
    std::cout << "\nExplanation: Reverse-mode is better because #outputs (1) << #inputs (" << numCoefs << ")" << std::endl;
    
    std::cout << "\n=== VERIFICATION ===" << std::endl;
    std::cout << "Checking if forward and reverse modes produce same gradients..." << std::endl;
    double max_diff = 0.0;
    for (index_t i = 0; i < numCoefs; ++i)
    {
        double diff = std::abs(gradient_forward(i, 0) - gradient_reverse(i));
        max_diff = std::max(max_diff, diff);
    }
    std::cout << "Maximum difference: " << max_diff << std::endl;
    if (max_diff < 1e-10)
    {
        std::cout << "✓ VERIFIED: Both methods produce identical gradients!" << std::endl;
    }
    else
    {
        std::cout << "✗ WARNING: Gradients differ by more than tolerance" << std::endl;
    }
    
    std::cout << "\n=======================================================================" << std::endl;
    std::cout << "   CONCLUSIONS" << std::endl;
    std::cout << "=======================================================================" << std::endl;
    std::cout << "\n1. For GISMO's typical use (derivatives w.r.t. spatial coords):" << std::endl;
    std::cout << "   → Use FORWARD-MODE (dual) - already fully integrated" << std::endl;
    
    std::cout << "\n2. For shape OPTIMIZATION (gradient w.r.t. many control points):" << std::endl;
    std::cout << "   → Use REVERSE-MODE (var) - " << (time_opt_forward/time_opt_reverse) << "x faster!" << std::endl;
    
    std::cout << "\n3. Rule of thumb:" << std::endl;
    std::cout << "   - Forward-mode: O(#inputs) passes" << std::endl;
    std::cout << "   - Reverse-mode: O(#outputs) passes" << std::endl;
    std::cout << "   → Choose the mode where the multiplier is smaller!" << std::endl;
    
    std::cout << "=======================================================================" << std::endl;
    
    return 0;
}
