/** @file test_var_basis.cpp

    @brief Minimal test to verify gsTensorBSpline works with autodiff::var

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Test
*/

#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff.h>

using namespace gismo;
using namespace autodiff;

int main(int argc, char* argv[])
{
    // Create a simple tensor B-spline basis
    gsKnotVector<> kv1(0, 1, 2, 3);  // degree 2
    gsKnotVector<> kv2(0, 1, 2, 3);
    gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);

    // Create control point coefficients
    gsMatrix<> coefs = basis.anchors().transpose();
    coefs.col(2) = gsMatrix<>::Zero(coefs.rows(), 1);

    // Try to create a surface with var type
    gsInfo << "Testing gsTensorBSpline<2, var>...\n";

    try
    {
        // Create basis with var
        gsTensorBSplineBasis<2, var> basis_var(kv1, kv2);
        gsInfo << "✓ Basis created successfully with var\n";

        // Try to evaluate basis functions
        gsMatrix<var> pts(2, 1);
        pts(0, 0) = var(0.5);
        pts(1, 0) = var(0.5);

        gsMatrix<var> values;
        basis_var.eval_into(pts, values);
        gsInfo << "✓ Basis evaluation successful\n";
        gsInfo << "Number of basis functions: " << basis_var.nBases() << "\n";

        // Try with coefficients
        gsMatrix<var> coefs_var(coefs.rows(), coefs.cols());
        for (index_t i = 0; i < coefs.rows(); ++i)
            for (index_t j = 0; j < coefs.cols(); ++j)
                coefs_var(i, j) = var(coefs(i, j));

        gsTensorBSpline<2, var> surface_var(basis_var, coefs_var);
        gsInfo << "✓ Tensor BSpline surface created with var\n";

        // Evaluate the surface
        gsMatrix<var> eval_result;
        surface_var.eval_into(pts, eval_result);
        gsInfo << "✓ Surface evaluation successful\n";
        gsInfo << "Evaluation result: " << eval_result << "\n";
    }
    catch (const std::exception& e)
    {
        gsInfo << "✗ Error: " << e.what() << "\n";
        return 1;
    }

    gsInfo << "\nAll tests passed!\n";
    return 0;
}
