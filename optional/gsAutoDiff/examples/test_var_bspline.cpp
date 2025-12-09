/** @file test_var_bspline.cpp

    @brief Test if gsTensorBSpline<2,var> works with autodiff.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Test
*/

#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>

using namespace gismo;
using namespace autodiff;

int main(int argc, char* argv[])
{
    // Create a simple B-spline surface with double coefficients
    gsKnotVector<> kv1(0, 1, 2, 2);  // 3 basis functions, degree 1
    gsKnotVector<> kv2(0, 1, 2, 2);
    gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);
    
    gsMatrix<> coefs(basis.size(), 3);
    coefs.setRandom();
    coefs.col(2).setZero();  // z = 0 for 2D surface
    
    // Create surface with double coefficients
    gsTensorBSpline<2, real_t> surfaceDouble(basis, coefs);
    
    gsInfo << "Created double surface with " << surfaceDouble.coefs().rows() << " control points\n";
    
    // Now try to create a surface with var coefficients
    // Key: coefficients CAN be var, but basis must be double and we evaluate with var coefs
    try {
        gsMatrix<var> coefsVar = coefs.cast<var>();
        gsInfo << "Successfully cast coefficients to var\n";
        
        // Keep basis as double, just use var coefficients in the surface
        gsTensorBSpline<2, var> surfaceVar(basis.clone(), coefsVar);
        gsInfo << "Successfully created surface with var coefficients\n";
        
        // Evaluate at a point
        gsMatrix<var> pt(2, 1);
        pt(0, 0) = 0.5;
        pt(1, 0) = 0.5;
        
        gsMatrix<var> result = surfaceVar.eval(pt);
        gsInfo << "Successfully evaluated surface with var\n";
        gsInfo << "Result has " << result.rows() << " rows and " << result.cols() << " cols\n";
    }
    catch (const std::exception& e) {
        gsInfo << "Error: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
