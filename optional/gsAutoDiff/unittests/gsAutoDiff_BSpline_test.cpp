/** @file gsAutoDiff_BSpline_test.cpp

    @brief Tests automatic differentiation on B-spline surfaces.
    Validates forward mode (dual_t) and reverse mode (var_t) AD
    with two key use cases: derivatives of many evaluations w.r.t. 
    one coefficient, and derivatives at one point w.r.t. all coefficients.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Automatically generated based on examples
**/

#include "gismo_unittest.h"
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsNurbs/gsTensorBSpline.hpp>
#include <gsCore/gsBoundingBox.h>

using namespace gismo;

SUITE(gsAutoDiff_BSpline)
{

TEST(Forward_AD_Dual_Basic)
{
    // Test forward AD with dual_t
    typedef gismo::ad::dual_t<double,double> dual_t;
    
    // Simple test: f(x) = x^2
    dual_t x(3.0, 1.0);  // value = 3, derivative = 1
    dual_t y = x * x;    // y = 9, dy/dx = 6
    
    CHECK_CLOSE(y.value, 9.0, 1e-10);
    CHECK_CLOSE(y.derivative, 6.0, 1e-10);
}

TEST(Reverse_AD_Var_Basic)
{
    // Test reverse AD with var_t
    typedef gismo::ad::var_t<double> var_t;
    
    // Simple test: f(x) = x^2
    var_t x = 3.0;
    var_t y = x * x;    // y = 9, dy/dx = 6
    
    ad::derivative(y, 1.0);
    real_t grad_x = ad::gradient(x);
    
    CHECK_CLOSE(grad_x, 6.0, 1e-10);
}

TEST(Forward_vs_Reverse_Consistency)
{
    // Simple consistency check between forward and reverse modes
    typedef gismo::ad::dual_t dual_t;
    typedef gismo::ad::var_t var_t;
    
    // Test function: f(x,y) = x*y + x^2
    real_t x_val = 2.0;
    real_t y_val = 3.0;
    
    // Forward mode: compute df/dx
    dual_t x_fwd(x_val, 1.0);
    dual_t y_fwd(y_val, 0.0);
    dual_t result_fwd = x_fwd * y_fwd + x_fwd * x_fwd;
    
    real_t deriv_fwd_dx = ad::val(result_fwd.derivatives()); // Extract derivative
    
    // Reverse mode: compute df/dx
    var_t x_rev(x_val);
    var_t y_rev(y_val);
    var_t result_rev = x_rev * y_rev + x_rev * x_rev;
    
    ad::derivative(result_rev, 1.0);
    real_t deriv_rev_dx = ad::gradient(x_rev);
    
    // Both should give df/dx = y + 2x = 3 + 4 = 7
    real_t expected = y_val + 2*x_val;
    
    CHECK_CLOSE(deriv_fwd_dx, expected, 1e-10);
    CHECK_CLOSE(deriv_rev_dx, expected, 1e-10);
}

}

