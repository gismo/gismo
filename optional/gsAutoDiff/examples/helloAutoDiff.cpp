/** @file helloMUMPS.cpp

    @brief First example of submodule

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

std::ostream &operator<<(std::ostream &out, const gsMatrix<dual> &s)
{
    return operator<<(out,s.template cast<double>());
}

dual testFunction(dual x)
{
    return 1 + x + x*x + 1/x + log(x);
}

int main()
{
    std::cout << "========================================" << std::endl;
    std::cout << "   FORWARD-MODE AD (dual numbers)" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Example 1: Computing derivatives directly from dual numbers
    std::cout << "\n=== Example 1: Direct derivative extraction ===" << std::endl;
    dual x = 1.0;
    autodiff::detail::seed<1>(x, 1.0);  // Seed x for differentiation
    dual u = testFunction(x);
    
    double value = autodiff::detail::derivative<0>(u);  // Get function value
    double dudx = autodiff::detail::derivative<1>(u);   // Get first derivative
    
    std::cout << "f(x) = " << value << std::endl;
    std::cout << "f'(x) = " << dudx << std::endl;

    // Example 2: B-spline evaluation with derivatives
    std::cout << "\n=== Example 2: B-spline eval() derivatives ===" << std::endl;
    gsKnotVector<dual> kv(0.0, 1.0, 2, 2);
    gsBSplineBasis<dual> basis(kv);
    gsMatrix<dual> coefs = gsMatrix<dual>::Random(basis.size(), 2);
    std::cout << "Coefficients:\n" << coefs << std::endl;

    gsBSpline<dual> spline(basis, coefs);
    gsMatrix<dual> pts(1,1);
    pts << 0.5;
    
    std::cout << "\nEvaluating spline at point " << pts(0,0) << std::endl;
    
    // Evaluate spline without seeding - just get values
    gsMatrix<dual> result = spline.eval(pts);
    std::cout << "Spline values (component 0): " << autodiff::detail::derivative<0>(result(0,0)) << std::endl;
    std::cout << "Spline values (component 1): " << autodiff::detail::derivative<0>(result(1,0)) << std::endl;

    // Example 3: Computing derivatives of eval() w.r.t. coefficients
    std::cout << "\n=== Example 3: Derivatives of eval() w.r.t. coefficients ===" << std::endl;
    
    std::cout << "Computing Jacobian (derivatives of each output w.r.t. each input):" << std::endl;
    
    // Compute full Jacobian: d eval() / d coefs
    // For each coefficient, compute how it affects each output
    gsMatrix<double> jacobian(result.rows(), coefs.rows());
    
    for (index_t k=0; k<coefs.rows(); ++k)
    {
        // Seed coefficient k (first column only for simplicity)
        gsMatrix<dual> coefs_seeded = coefs;
        autodiff::detail::seed<1>(coefs_seeded(k,0), 1.0);
        gsBSpline<dual> spline_k(basis, coefs_seeded);
        gsMatrix<dual> res_k = spline_k.eval(pts);
        
        // Extract derivatives for all outputs at once using derivative<>
        jacobian(0, k) = autodiff::detail::derivative<1>(res_k(0,0));
        jacobian(1, k) = autodiff::detail::derivative<1>(res_k(1,0));
    }
    
    std::cout << "Jacobian d eval() / d coefs (first column):\n" << jacobian << std::endl;
    
    std::cout << "\nNote: Only coefficients 1 and 2 have non-zero derivatives" << std::endl;
    std::cout << "because only basis functions 1 and 2 are active at x=0.5" << std::endl;

    std::cout << "\nBasis functions at evaluation point:\n" << basis.eval(pts) << std::endl;
    
    // Example 4: Verify derivatives match basis functions
    std::cout << "\n=== Example 4: Verification ===" << std::endl;
    std::cout << "The derivatives should match the basis function values" << std::endl;
    std::cout << "since spline.eval() = sum(coef[i] * basis[i])" << std::endl;

    std::cout << "\n\n========================================" << std::endl;
    std::cout << "   DUAL vs VAR: Conceptual Comparison" << std::endl;
    std::cout << "========================================" << std::endl;
    
    std::cout << "\n=== Example 5: Understanding Forward vs Reverse Mode ===" << std::endl;
    std::cout << "\nFORWARD MODE (autodiff::dual):" << std::endl;
    std::cout << "  - Implemented and fully supported in GISMO" << std::endl;
    std::cout << "  - Works with gsMatrix<dual>, gsBSpline<dual>, etc." << std::endl;
    std::cout << "  - Computes derivatives by propagating dual numbers forward" << std::endl;
    std::cout << "  - One pass per input variable to get all outputs' derivatives" << std::endl;
    std::cout << "  - Efficient when #inputs <= #outputs" << std::endl;
    std::cout << "\n  Example: For a spline with " << coefs.rows() << " coefficients:" << std::endl;
    std::cout << "    - To get Jacobian (all outputs w.r.t. all inputs): " << coefs.rows() << " forward passes" << std::endl;
    std::cout << "    - Each pass: seed one coefficient, evaluate spline" << std::endl;
    std::cout << "    - Result: one column of Jacobian per pass" << std::endl;
    
    std::cout << "\nREVERSE MODE (autodiff::var):" << std::endl;
    std::cout << "  - Available in autodiff library" << std::endl;
    std::cout << "  - NOT YET fully integrated with GISMO types" << std::endl;
    std::cout << "  - Requires NumTraits and other Eigen support" << std::endl;
    std::cout << "  - Computes derivatives by backpropagation" << std::endl;
    std::cout << "  - One pass per output variable to get all inputs' derivatives" << std::endl;
    std::cout << "  - Efficient when #outputs <= #inputs" << std::endl;
    std::cout << "\n  Example: For the same spline with " << coefs.rows() << " coefficients:" << std::endl;
    std::cout << "    - To get gradient of one output w.r.t. all inputs: 1 reverse pass" << std::endl;
    std::cout << "    - Each pass: evaluate forward, then backpropagate from one output" << std::endl;
    std::cout << "    - Result: one row of Jacobian per pass" << std::endl;
    
    std::cout << "\n=== Example 6: When to use which mode ===" << std::endl;
    std::cout << "\nUse FORWARD MODE (dual) when:" << std::endl;
    std::cout << "  1. Computing derivatives w.r.t. spatial coordinates (x, y, z)" << std::endl;
    std::cout << "     Example: df/dx, df/dy for a function f(x,y)" << std::endl;
    std::cout << "  2. Computing Jacobians with few inputs, many outputs" << std::endl;
    std::cout << "  3. Working with GISMO geometry objects (splines, patches, etc.)" << std::endl;
    std::cout << "  4. Need is ALREADY SUPPORTED in GISMO!" << std::endl;
    
    std::cout << "\nUse REVERSE MODE (var) when:" << std::endl;
    std::cout << "  1. Computing gradients for optimization (many parameters, one objective)" << std::endl;
    std::cout << "     Example: gradient of error function w.r.t. 1000s of coefficients" << std::endl;
    std::cout << "  2. Computing Jacobians with many inputs, few outputs" << std::endl;
    std::cout << "  3. Machine learning / neural network applications" << std::endl;
    std::cout << "  4. Would require additional GISMO integration work" << std::endl;
    
    std::cout << "\n=== Example 7: Efficiency comparison ===" << std::endl;
    std::cout << "\nScenario: Spline with " << coefs.rows() << " coefficients, " << result.rows() << " output components" << std::endl;
    std::cout << "\nCompute full Jacobian (all outputs w.r.t. all inputs):" << std::endl;
    std::cout << "  - Forward mode: " << coefs.rows() << " passes (one per input)" << std::endl;
    std::cout << "  - Reverse mode: " << result.rows() << " passes (one per output)" << std::endl;
    std::cout << "  - Winner: " << (coefs.rows() < result.rows() ? "Forward" : "Reverse") << " mode" << std::endl;
    
    std::cout << "\nCompute gradient of scalar objective w.r.t. all " << coefs.rows() << " coefficients:" << std::endl;
    std::cout << "  - Forward mode: " << coefs.rows() << " passes" << std::endl;
    std::cout << "  - Reverse mode: 1 pass" << std::endl;
    std::cout << "  - Winner: Reverse mode (much faster!)" << std::endl;
    
    std::cout << "\nCurrent example (derivatives w.r.t. spatial point):" << std::endl;
    std::cout << "  - Inputs: 1-3 (x, y, z coordinates)" << std::endl;
    std::cout << "  - Outputs: " << result.rows() << " (spline components)" << std::endl;
    std::cout << "  - Winner: Forward mode (as used in GISMO)" << std::endl;
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "For GISMO development, forward mode (dual) is the primary choice" << std::endl;
    std::cout << "because most use cases involve derivatives w.r.t. few spatial variables." << std::endl;
    std::cout << "Reverse mode (var) could be beneficial for future shape optimization" << std::endl;
    std::cout << "or parameter fitting applications with many unknowns." << std::endl;
    std::cout << "========================================\n" << std::endl;
}


// int main(int argc, char *argv[])
// {
//     // Size of global sparse matrix
//     index_t mat_size = 10;
//     std::string spm(""); // sparse matrix from a file

//     gsCmdLine cmd("Testing the use of sparse linear solvers.");
//     cmd.addInt("n", "size", "Size of the matrices", mat_size);
//     cmd.addString("m", "matrix", "Filename to read sparse matrix and right hand side", spm);

//     try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


//     return EXIT_SUCCESS;
// }
