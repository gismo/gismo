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
    // Example 1: Computing derivatives directly from dual numbers
    std::cout << "=== Example 1: Direct derivative extraction ===" << std::endl;
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
