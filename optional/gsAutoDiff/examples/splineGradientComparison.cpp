/** @file splineGradientComparison.cpp

    @brief Performance comparison: forward-mode (dual) vs reverse-mode (var)

    Demonstrates the key difference between forward and reverse mode AD by computing
    the gradient of a scalar loss function with respect to many parameters.
    
    Scenario: Fitting B-spline curve to target points by minimizing least-squares error.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst

*/

#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsAutoDiff/gsAutoDiffVar.h>

using namespace gismo;
using namespace autodiff;

int main(int argc, char *argv[])
{
    index_t numCoefs = 50;   // Number of control points (parameters)
    index_t numEvals = 100;  // Number of evaluation points
    bool validate = false;   // Validate gradients against analytical computation
    gsCmdLine cmd("Compares forward-mode vs reverse-mode AD performance for spline gradient computation.");
    cmd.addInt("C","numCoefs", "Number of control points (parameters)", numCoefs);
    cmd.addInt("E","numEvals", "Number of evaluation points", numEvals);
    cmd.addSwitch("validate", "Validate gradients against analytical basis function computation", validate);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    std::cout << "======================================================================" << std::endl;
    std::cout << "   SPLINE GRADIENT: FORWARD-MODE (dual) vs REVERSE-MODE (var)" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    std::cout << "Problem: Minimize least-squares error for spline curve fitting" << std::endl;
    std::cout << "  - Number of control points (parameters): " << numCoefs << std::endl;
    std::cout << "  - Number of evaluation points: " << numEvals << std::endl;
    std::cout << "  - Task: Compute gradient of scalar loss w.r.t. all parameters" << std::endl;
    std::cout << "  - Validation: " << (validate ? "enabled" : "disabled") << "\n" << std::endl;
    
    // Create B-spline with double type
    const index_t degree = 3;
    gsKnotVector<> kv(0.0, 1.0, numCoefs - degree - 1, degree + 1);
    gsBSplineBasis<> basis(kv);
    gsMatrix<> coefs = gsMatrix<>::Random(numCoefs, 2);
    gsBSpline<> spline(basis, coefs);
    
    // Create evaluation points
    gsMatrix<> evalPts(1, numEvals);
    for (index_t i = 0; i < numEvals; ++i)
    {
        evalPts(0, i) = 0.01 + (0.99 - 0.01) * i / (numEvals - 1);
    }
    
    // Create target values (what we're fitting to)
    gsMatrix<> target = gsMatrix<>::Random(2, numEvals);
    
    // Compute analytical gradient using basis functions (for validation)
    gsMatrix<> gradient_analytical(numCoefs, 2);
    gradient_analytical.setZero();
    if (validate)
    {
        // Evaluate ALL basis functions densely
        gsMatrix<> basisValues_analytical(numCoefs, numEvals);
        for (index_t k = 0; k < numCoefs; ++k)
        {
            gsMatrix<> singleBasisResult;
            basis.evalSingle_into(k, evalPts, singleBasisResult);
            basisValues_analytical.row(k) = singleBasisResult.row(0);
        }
        gsMatrix<> splineValues = spline.eval(evalPts);
        
        for (index_t k = 0; k < numCoefs; ++k)
        {
            for (index_t dim = 0; dim < 2; ++dim)
            {
                double grad = 0.0;
                for (index_t i = 0; i < numEvals; ++i)
                {
                    // dLoss/dCoef_k = 2 * sum_i( (spline(u_i) - target_i) * dSpline/dCoef_k )
                    // where dSpline/dCoef_k = basis_k(u_i)
                    double residual = splineValues(dim, i) - target(dim, i);
                    grad += 2.0 * residual * basisValues_analytical(k, i);
                }
                gradient_analytical(k, dim) = grad;
            }
        }
    }
    
    std::cout << "======================================================================" << std::endl;
    std::cout << "   METHOD 1: FORWARD-MODE (dual)" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    std::cout << "Strategy: Seed one parameter at a time, evaluate full forward pass" << std::endl;
    std::cout << "  - Number of forward passes needed: " << numCoefs << " (one per parameter)" << std::endl;
    std::cout << "  - Each pass: evaluate spline + compute loss at all " << numEvals << " points\n" << std::endl;
    
    gsStopwatch timer;
    
    gsMatrix<> gradient_forward(numCoefs, 2);
    gradient_forward.setZero();
    
    for (index_t k = 0; k < numCoefs; ++k)
    {
        for (index_t dim = 0; dim < 2; ++dim)
        {
            // Seed coefficient k for dimension dim
            gsMatrix<dual> coefs_dual = coefs.cast<dual>();
            autodiff::detail::seed<1>(coefs_dual(k, dim), 1.0);
        
            // Create spline with seeded coefficients
            gsKnotVector<dual> kv_dual(dual(0.0), dual(1.0), numCoefs - degree - 1, degree + 1);
            gsBSplineBasis<dual> basis_dual(kv_dual);
            gsBSpline<dual> spline_dual(basis_dual, coefs_dual);
            
            // Evaluate spline at all points
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
            
            // Extract gradient component for parameter k, dimension dim
            gradient_forward(k, dim) = autodiff::detail::derivative<1>(loss);
        }
    }
    
    double time_forward = timer.stop() * 1000.0;  // Convert to milliseconds
    
    // Estimate memory consumption for forward-mode
    size_t mem_forward = numCoefs * sizeof(dual) * 2 +  // dual coefficients
                         numEvals * sizeof(dual) * 2 +   // dual evaluation results
                         numCoefs * 2 * sizeof(double);  // gradient storage
    
    std::cout << "Result: Successfully computed gradient with " << numCoefs << " components" << std::endl;
    std::cout << "Time: " << time_forward << " ms" << std::endl;
    std::cout << "Memory (estimated): " << mem_forward / 1024.0 << " KB" << std::endl;
    std::cout << "Computational complexity: O(numCoefs × numEvals) = O(" << numCoefs << " × " << numEvals << ")" << std::endl;
    
    if (validate)
    {
        double maxDiff = (gradient_forward - gradient_analytical).cwiseAbs().maxCoeff();
        std::cout << "Validation: max difference = " << maxDiff << " (analytical vs forward-mode)" << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << "======================================================================" << std::endl;
    std::cout << "   METHOD 2: REVERSE-MODE (var) with element-wise operations" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    std::cout << "Strategy: Build computation graph, then backpropagate from single loss" << std::endl;
    std::cout << "  - Number of backward passes needed: 1 (from single output)" << std::endl;
    std::cout << "  - Single pass: backprop through all " << numEvals << " evaluation points\n" << std::endl;
    
    timer.restart();
    
    // Evaluate ALL basis functions at evaluation points (dense matrix: numCoefs x numEvals)
    gsMatrix<> basisValues(numCoefs, numEvals);
    basisValues.setZero();
    for (index_t k = 0; k < numCoefs; ++k)
    {
        gsMatrix<> singleBasisResult;
        basis.evalSingle_into(k, evalPts, singleBasisResult);
        // singleBasisResult is (1 x numEvals) - copy to row k
        basisValues.row(k) = singleBasisResult.row(0);
    }
    
    gsMatrix<> gradient_reverse(numCoefs, 2);
    gradient_reverse.setZero();
    
    // Compute gradients for each dimension separately (autodiff::gradient requires vectors)
    for (index_t dim = 0; dim < 2; ++dim)
    {
        // Create var-type coefficient vector for this dimension
        gsEigen::Matrix<var, gsEigen::Dynamic, 1> coefs_var(numCoefs);
        for (index_t i = 0; i < numCoefs; ++i)
        {
            coefs_var(i) = coefs(i, dim);
        }
        // Build computation graph using actual basis function values
        var loss_var = 0.0;
        for (index_t i = 0; i < numEvals; ++i)
        {
            // Compute spline value at evaluation point i using basis functions
            var spline_value = 0.0;
            for (index_t k = 0; k < numCoefs; ++k)
            {
                spline_value = spline_value + coefs_var(k) * basisValues(k, i);
            }
            var diff = spline_value - target(dim, i);
            loss_var = loss_var + diff * diff;
        }
        
        // ONE backward pass to get ALL gradients for this dimension
        gsEigen::VectorXd grad_dim = autodiff::gradient(loss_var, coefs_var);
        for (index_t k = 0; k < numCoefs; ++k)
        {
            gradient_reverse(k, dim) = grad_dim[k];
        }
    }
    
    double time_reverse = timer.stop() * 1000.0;  // Convert to milliseconds
    
    // Estimate memory consumption for reverse-mode
    // Each var node stores: value (8 bytes) + expression pointer (8 bytes) + gradient (8 bytes)
    size_t nodes_in_graph = numCoefs * 2 +                    // coefficient nodes
                            numEvals * 2 +                     // spline value nodes
                            numEvals * 2 +                     // difference nodes
                            numEvals * 2 +                     // squared difference nodes
                            1;                                 // loss node
    size_t mem_reverse = nodes_in_graph * 24 +                // var nodes (approximate)
                         numCoefs * 2 * sizeof(double) +       // gradient storage
                         numEvals * numCoefs * sizeof(double); // basis values
    
    std::cout << "Result: Successfully computed gradient with " << numCoefs << " components" << std::endl;
    std::cout << "Time: " << time_reverse << " ms" << std::endl;
    std::cout << "Memory (estimated): " << mem_reverse / 1024.0 << " KB" << std::endl;
    std::cout << "  - Computation graph nodes: " << nodes_in_graph << std::endl;
    std::cout << "Computational complexity: O(numEvals) + O(graph_size) ≈ O(" << numEvals << ")" << std::endl;
    
    if (validate)
    {
        double maxDiff = (gradient_reverse - gradient_analytical).cwiseAbs().maxCoeff();
        std::cout << "Validation: max difference = " << maxDiff << " (analytical vs reverse-mode)" << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << "======================================================================" << std::endl;
    std::cout << "   PERFORMANCE SUMMARY" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    std::cout << "Task: Compute gradient of scalar loss w.r.t. " << numCoefs << " parameters\n" << std::endl;
    
    std::cout << "Forward-mode (dual): " << time_forward << " ms, "
              << mem_forward / 1024.0 << " KB memory" << std::endl;
    std::cout << "  - Requires " << numCoefs << " separate evaluations" << std::endl;
    std::cout << "  - Each evaluation processes " << numEvals << " points" << std::endl;
    std::cout << "  - Best for: few parameters, many outputs\n" << std::endl;
    
    std::cout << "Reverse-mode (var):  " << time_reverse << " ms, "
              << mem_reverse / 1024.0 << " KB memory" << std::endl;
    std::cout << "  - Requires only 1 backpropagation" << std::endl;
    std::cout << "  - Processes all " << numEvals << " points in single pass" << std::endl;
    std::cout << "  - Best for: many parameters, few outputs\n" << std::endl;
    
    if (time_forward > time_reverse)
    {
        double speedup = time_forward / time_reverse;
        double mem_ratio = (double)mem_reverse / mem_forward;
        std::cout << "🚀 REVERSE-MODE WINS: " << speedup << "x faster!" << std::endl;
        std::cout << "   Memory overhead: " << mem_ratio << "x (" 
                  << (mem_reverse - mem_forward) / 1024.0 << " KB extra)" << std::endl;
        std::cout << "\nFor optimization with " << numCoefs << " parameters:" << std::endl;
        std::cout << "  → Use reverse-mode (var) for gradient computation" << std::endl;
        std::cout << "  → Speedup increases with more parameters!" << std::endl;
    }
    else
    {
        std::cout << "⚠️  Forward-mode faster (due to overhead in small example)" << std::endl;
        std::cout << "  → Reverse-mode advantage grows with parameter count" << std::endl;
    }
    
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "   CONCLUSIONS" << std::endl;
    std::cout << "======================================================================\n" << std::endl;
    
    std::cout << "1. GISMO typical use (geometry evaluation with spatial derivatives):" << std::endl;
    std::cout << "   → Use FORWARD-MODE (dual)" << std::endl;
    std::cout << "   → Few inputs (x,y,z) → many outputs (function values)" << std::endl;
    std::cout << "   → Already fully integrated in GISMO\n" << std::endl;
    
    std::cout << "2. Shape OPTIMIZATION (fitting splines to data):" << std::endl;
    std::cout << "   → Use REVERSE-MODE (var)" << std::endl;
    std::cout << "   → Many inputs (control points) → one output (loss)" << std::endl;
    std::cout << "   → Up to " << (time_forward/time_reverse) << "x faster for this problem size\n" << std::endl;
    
    std::cout << "3. Rule of thumb:" << std::endl;
    std::cout << "   → Forward-mode: Efficient when #inputs ≤ #outputs" << std::endl;
    std::cout << "   → Reverse-mode: Efficient when #outputs ≤ #inputs" << std::endl;
    std::cout << "   → For optimization (1 loss → N params): use reverse-mode!" << std::endl;
    
    std::cout << "======================================================================" << std::endl;
    
    return 0;
}
