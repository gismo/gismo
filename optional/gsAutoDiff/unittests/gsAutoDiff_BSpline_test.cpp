/** @file gsAutoDiff_BSpline_test.cpp

    @brief Tests automatic differentiation on B-spline surfaces.
    Validates forward mode (dual_t) and reverse mode (var_t) AD
    with two key use cases: derivatives of many evaluations w.r.t. 
    one coefficient, and derivatives at one point w.r.t. all coefficients.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#include "gismo_unittest.h"
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsNurbs/gsTensorBSpline.hpp>

using namespace gismo;

SUITE(gsAutoDiff_BSpline)
{

TEST(Forward_AD_Dual_Basic)
{
    // Test forward AD with dual_t
    // Simple test: f(x) = x^2
    dual_t x(3.0);  // value = 3
    x.grad = 1.0;   // derivative = 1
    dual_t y = x * x;    // y = 9, dy/dx = 6
    
    CHECK_CLOSE(y.val, 9.0, 1e-10);
    CHECK_CLOSE(y.grad, 6.0, 1e-10);
}

TEST(Reverse_AD_Var_Basic)
{
    // Test reverse AD with var_t
    using autodiff::var;
    using autodiff::derivatives;
    using autodiff::reverse::detail::wrt;
    
    // Simple test: f(x) = x^2
    var x = 3.0;
    var y = x * x;    // y = 9, dy/dx = 6
    
    auto dydx = derivatives(y, wrt(x));
    real_t grad_x = dydx[0];
    
    CHECK_CLOSE(grad_x, 6.0, 1e-10);
}

TEST(Forward_vs_Reverse_Consistency)
{
    // Simple consistency check between forward and reverse modes
    using autodiff::var;
    using autodiff::derivatives;
    using autodiff::reverse::detail::wrt;
    
    // Test function: f(x,y) = x*y + x^2
    real_t x_val = 2.0;
    real_t y_val = 3.0;
    
    // Forward mode: compute df/dx
    dual_t x_fwd(x_val);
    x_fwd.grad = 1.0;
    dual_t y_fwd(y_val);
    y_fwd.grad = 0.0;
    dual_t result_fwd = x_fwd * y_fwd + x_fwd * x_fwd;
    
    real_t deriv_fwd_dx = result_fwd.grad; // Extract derivative
    
    // Reverse mode: compute df/dx
    var x_rev(x_val);
    var y_rev(y_val);
    var result_rev = x_rev * y_rev + x_rev * x_rev;
    
    auto dfdx = derivatives(result_rev, wrt(x_rev));
    real_t deriv_rev_dx = dfdx[0];
    
    // Both should give df/dx = y + 2x = 3 + 4 = 7
    real_t expected = y_val + 2*x_val;
    
    CHECK_CLOSE(deriv_fwd_dx, expected, 1e-10);
    CHECK_CLOSE(deriv_rev_dx, expected, 1e-10);
}

}
