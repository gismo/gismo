#include "gismo_unittest.h"
#include <gsAutoDiff/gsAutoDiff2.h>

// Include implementation files to force instantiation for custom types
#include <gsCore/gsBasis.hpp>
#include <gsCore/gsFunctionSet.hpp>
#include <gsCore/gsFunction.hpp>
#include <gsCore/gsFunctionExpr.hpp>
#include <gsCore/gsConstantFunction.hpp>
#include <gsCore/gsGeometry.hpp>
#include <gsCore/gsCurve.hpp>
#include <gsCore/gsSurface.hpp>
#include <gsTensor/gsTensorBasis.hpp>
#include <gsNurbs/gsTensorBSplineBasis.hpp>
#include <gsNurbs/gsBSplineBasis.hpp>
#include <gsNurbs/gsBSpline.hpp>
#include <gsNurbs/gsTensorBSpline.hpp>
#include <gsNurbs/gsKnotVector.hpp>
#include <gsNurbs/gsBoehm.hpp>
#include <gsNurbs/gsBSplineSolver.hpp>
#include <gsUtils/gsMesh/gsMesh.hpp>

using namespace gismo;

SUITE(gsAutoDiff_test)
{

TEST(VarBSplineTest)
{
    gismo::var_t a = 0;
    gismo::var_t b = 1;
    index_t interior = 4;
    index_t multEnd = 3;

    gsKnotVector<gismo::var_t> kv(a, b, interior, multEnd);
    gsBSplineBasis<gismo::var_t> bsb(kv);
    
    CHECK_EQUAL(7, bsb.size());
    CHECK_EQUAL(2, bsb.degree());
}

TEST(VarArithmeticTest)
{
    // Test basic arithmetic
    var_t x = 2.0;
    var_t y = x * x + 1.0;
    
    // Compute derivatives
    auto [dydx] = autodiff::reverse::detail::derivatives(y.var, autodiff::reverse::detail::wrt(x.var));
    
    CHECK_CLOSE(5.0, y.val(), 1e-10);
    CHECK_CLOSE(4.0, dydx, 1e-10);

    // Test math functions
    var_t z = sin(x);
    auto [dzdx] = autodiff::reverse::detail::derivatives(z.var, autodiff::reverse::detail::wrt(x.var));
    
    CHECK_CLOSE(std::sin(2.0), z.val(), 1e-10);
    CHECK_CLOSE(std::cos(2.0), dzdx, 1e-10);

    // Test integration with G+Smo types
    gsMatrix<var_t> M(2, 2);
    M(0,0) = x; M(0,1) = 1.0;
    M(1,0) = 0.0; M(1,1) = x;
    
    gsMatrix<var_t> M2 = M * M;
    
    // Compute derivatives for matrix elements
    auto [dM2_00_dx] = autodiff::reverse::detail::derivatives(M2(0,0).var, autodiff::reverse::detail::wrt(x.var));
    auto [dM2_01_dx] = autodiff::reverse::detail::derivatives(M2(0,1).var, autodiff::reverse::detail::wrt(x.var));
    
    CHECK_CLOSE(4.0, M2(0,0).val(), 1e-10);
    CHECK_CLOSE(4.0, dM2_00_dx, 1e-10);
    
    CHECK_CLOSE(4.0, M2(0,1).val(), 1e-10);
    CHECK_CLOSE(2.0, dM2_01_dx, 1e-10);
}

TEST(DualArithmeticTest)
{
    // Test basic arithmetic
    gismo::dual_t x = 2.0;
    x.get().grad = 1.0; // Seed gradient
    gismo::dual_t y = x * x + 1.0;
    
    CHECK_CLOSE(5.0, y.val(), 1e-10);
    CHECK_CLOSE(4.0, y.grad(), 1e-10);
    
    // Test math functions
    gismo::dual_t z = sin(x);
    CHECK_CLOSE(std::sin(2.0), z.val(), 1e-10);
    CHECK_CLOSE(std::cos(2.0), z.grad(), 1e-10);

    // Test integration with G+Smo types
    gsMatrix<dual_t> M(2, 2);
    M(0,0) = x; M(0,1) = 1.0;
    M(1,0) = 0.0; M(1,1) = x;
    
    gsMatrix<dual_t> M2 = M * M;
    
    CHECK_CLOSE(4.0, M2(0,0).val(), 1e-10);
    CHECK_CLOSE(4.0, M2(0,0).grad(), 1e-10);
    
    CHECK_CLOSE(4.0, M2(0,1).val(), 1e-10);
    CHECK_CLOSE(2.0, M2(0,1).grad(), 1e-10);
}

TEST(BSplineSurface_Forward_Dual)
{
    // Case 1: Derivative of multiple points w.r.t. one coefficient change
    // Forward mode is efficient here because we have 1 input (coeff) and many outputs (points)
    
    using T = gismo::dual_t;
    
    index_t n = 5;
    index_t m = 5;
    index_t degree = 2;
    index_t num_points = 10; // Reduced for unittest
    
    gsKnotVector<T> kv_u(0, 1, n - degree - 1, degree + 1);
    gsKnotVector<T> kv_v(0, 1, m - degree - 1, degree + 1);
    
    gsTensorBSplineBasis<2, T> basis(kv_u, kv_v);
    gsMatrix<T> coefs(basis.size(), 3);
    coefs.setRandom();
    
    // Target coefficient to differentiate with respect to
    index_t target_coeff = basis.size() / 2;
    
    // Seed the gradient for the target coefficient
    coefs(target_coeff, 2).get().grad = 1.0;
    
    gsTensorBSpline<2, T> surface(basis, coefs);
    
    gsMatrix<T> eval_points(2, num_points);
    eval_points.setRandom();
    
    gsMatrix<T> result;
    surface.eval_into(eval_points, result);
    
    // Validation against exact derivatives
    gsMatrix<T> basis_values = basis.evalSingle(target_coeff, eval_points);
    
    for (index_t i = 0; i < result.cols(); ++i) {
        double ad_deriv = result(2, i).grad();
        double exact_deriv = (double)basis_values(0, i);
        CHECK_CLOSE(exact_deriv, ad_deriv, 1e-10);
    }
}

TEST(BSplineSurface_Reverse_Var)
{
    // Case 2: Derivative of one point w.r.t. all coefficient changes
    // Reverse mode is efficient here because we have many inputs (coeffs) and 1 output (point)
    
    using T = gismo::var_t;
    
    index_t n = 5;
    index_t m = 5;
    index_t degree = 2;
    
    gsKnotVector<T> kv_u(0, 1, n - degree - 1, degree + 1);
    gsKnotVector<T> kv_v(0, 1, m - degree - 1, degree + 1);
    
    gsTensorBSplineBasis<2, T> basis(kv_u, kv_v);
    gsMatrix<T> coefs(basis.size(), 3);
    coefs.setRandom();
    
    gsTensorBSpline<2, T> surface(basis, coefs);
    
    gsMatrix<T> point(2, 1);
    point.setRandom();
    
    gsMatrix<T> result;
    surface.eval_into(point, result);
    
    // We compute derivatives of result(2,0) w.r.t. coefs(k, 2)
    
    index_t k = basis.size() / 2;
    auto [ad_deriv] = autodiff::reverse::detail::derivatives(result(2, 0).var, autodiff::reverse::detail::wrt(coefs(k, 2).var));
    
    gsMatrix<T> basis_val = basis.evalSingle(k, point);
    double exact_deriv = (double)basis_val(0, 0);
    CHECK_CLOSE(exact_deriv, ad_deriv, 1e-10);
    
    // Check another one
    k = 0;
    auto [ad_deriv_0] = autodiff::reverse::detail::derivatives(result(2, 0).var, autodiff::reverse::detail::wrt(coefs(k, 2).var));
    gsMatrix<T> basis_val_0 = basis.evalSingle(k, point);
    double exact_deriv_0 = (double)basis_val_0(0, 0);
    CHECK_CLOSE(exact_deriv_0, ad_deriv_0, 1e-10);
}

}
