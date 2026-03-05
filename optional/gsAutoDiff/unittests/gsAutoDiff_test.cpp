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
    // Test var_t with B-spline basis
    autodiff::var a = 0;
    autodiff::var b = 1;
    index_t interior = 4;
    index_t multEnd = 3;

    gsKnotVector<autodiff::var> kv(a, b, interior, multEnd);
    gsBSplineBasis<autodiff::var> bsb(kv);
    
    CHECK_EQUAL(7, bsb.size());
    CHECK_EQUAL(2, bsb.degree());
}

TEST(VarArithmeticTest)
{
    // Test basic arithmetic with var
    autodiff::var x = 2.0;
    autodiff::var y = x * x + 1.0;
    
    // Compute derivatives using autodiff API
    auto [dydx] = autodiff::derivatives(y, autodiff::reverse::detail::wrt(x));
    
    CHECK_CLOSE(static_cast<double>(y), 5.0, 1e-10);
    CHECK_CLOSE(dydx, 4.0, 1e-10);
    
    // Test math functions
    autodiff::var z = autodiff::reverse::detail::sin(x);
    auto [dzdx] = autodiff::derivatives(z, autodiff::reverse::detail::wrt(x));
    
    CHECK_CLOSE(static_cast<double>(z), std::sin(2.0), 1e-10);
    CHECK_CLOSE(dzdx, std::cos(2.0), 1e-10);
    
    // Test integration with G+Smo types
    gsMatrix<autodiff::var> M(2, 2);
    M(0,0) = x; M(0,1) = 1.0;
    M(1,0) = 0.0; M(1,1) = x;
    
    gsMatrix<autodiff::var> M2 = M * M;
    
    // Compute derivatives for matrix elements
    auto [dM2_00_dx] = autodiff::derivatives(M2(0,0), autodiff::reverse::detail::wrt(x));
    auto [dM2_01_dx] = autodiff::derivatives(M2(0,1), autodiff::reverse::detail::wrt(x));
    
    CHECK_CLOSE(static_cast<double>(M2(0,0)), 4.0, 1e-10);
    CHECK_CLOSE(dM2_00_dx, 4.0, 1e-10);
    
    CHECK_CLOSE(static_cast<double>(M2(0,1)), 4.0, 1e-10);
    CHECK_CLOSE(dM2_01_dx, 2.0, 1e-10);
}

TEST(DualArithmeticTest)
{
    // Test basic arithmetic with dual numbers
    autodiff::detail::Dual<double, double> x = 2.0;
    autodiff::detail::seed<1>(x, 1.0);  // Seed x for differentiation
    autodiff::detail::Dual<double, double> y = x * x + 1.0;
    
    CHECK_CLOSE(autodiff::detail::derivative<0>(y), 5.0, 1e-10);
    CHECK_CLOSE(autodiff::detail::derivative<1>(y), 4.0, 1e-10);
        
    // Test integration with gsMatrix
    /* 
        M = | x  0   |
            | 0  x^2 |

        M^T = | x  0   |
              | 0  x^2 |

        M^2 = M^T * M = | x^2       0         |
                        | 0         x^4     |
        d(M^2_00)/dx = 2x = 4
        d(M^2_01)/dx = 0
        d(M^2_10)/dx = 0
        d(M^2_11)/dx = x^2*2*x+x^2*2*x = 4*x^3 = 32
    */
    
    gsMatrix<autodiff::detail::Dual<double, double>> M(2, 2);
    M(0,0) = x; M(0,1) = 0.0;
    M(1,0) = 0.0; M(1,1) = x*x;
    gsMatrix<autodiff::detail::Dual<double, double>> MT = M.transpose();
    
    gsMatrix<autodiff::detail::Dual<double, double>> M2 = MT * M;
    CHECK_CLOSE(autodiff::detail::derivative<1>(M2(0,0)), 4.0, 1e-10);
    CHECK_CLOSE(autodiff::detail::derivative<1>(M2(1,0)), 0.0, 1e-10);        
    CHECK_CLOSE(autodiff::detail::derivative<1>(M2(0,1)), 0.0, 1e-10);
    CHECK_CLOSE(autodiff::detail::derivative<1>(M2(1,1)), 32.0, 1e-10);
}

TEST(BSplineSurface_Forward_Dual)
{
    // Case 1: Derivative of multiple points w.r.t. one coefficient change
    // Forward mode is efficient here because we have 1 input (coeff) and many outputs (points)
    
    using T = autodiff::detail::Dual<double, double>;
    
    index_t n = 5;
    index_t m = 5;
    index_t degree = 2;
    
    gsKnotVector<T> kv_u(0, 1, n - degree - 1, degree + 1);
    gsKnotVector<T> kv_v(0, 1, m - degree - 1, degree + 1);
    
    gsTensorBSplineBasis<2, T> basis(kv_u, kv_v);
    gsMatrix<T> coefs(basis.size(), 3);
    coefs.setRandom();
    
    // Target coefficient to differentiate with respect to
    index_t target_coeff = basis.size() / 2;
    
    // Seed the gradient for the target coefficient
    autodiff::detail::seed<1>(coefs(target_coeff, 2), 1.0);
    
    gsTensorBSpline<2, T> surface(basis, coefs);
    
    gsVector<unsigned> np(2);
    np << 10, 10;
    gsVector<double> a(2);
    a << 0, 0;
    gsVector<double> b(2);
    b << 1, 1;
    gsMatrix<T> eval_points = gsPointGrid<double>(a, b, np).template cast<T>();

    
    gsMatrix<T> result;
    surface.eval_into(eval_points, result);
    
    // Validation against exact derivatives
    gsMatrix<T> basis_values = basis.evalSingle(target_coeff, eval_points);

    for (index_t i = 0; i < result.cols(); ++i) {
        double ad_deriv = autodiff::detail::derivative<1>(result(2, i));
        double exact_deriv = static_cast<double>(basis_values(0, i));
        CHECK_CLOSE(exact_deriv, ad_deriv, 1e-10);
    }
}

TEST(BSplineSurface_Reverse_Var)
{
    // Case 2: Derivative of one point w.r.t. all coefficient changes
    // Reverse mode is efficient here because we have many inputs (coeffs) and 1 output (point)
    
    using T = autodiff::var;
    
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
    auto [ad_deriv] = autodiff::derivatives(result(2, 0), autodiff::reverse::detail::wrt(coefs(k, 2)));
    
    gsMatrix<T> basis_val = basis.evalSingle(k, point);
    double exact_deriv = static_cast<double>(basis_val(0, 0));
    CHECK_CLOSE(exact_deriv, ad_deriv, 1e-10);
    
    // Check another one
    k = 0;
    auto [ad_deriv_0] = autodiff::derivatives(result(2, 0), autodiff::reverse::detail::wrt(coefs(k, 2)));
    gsMatrix<T> basis_val_0 = basis.evalSingle(k, point);
    double exact_deriv_0 = static_cast<double>(basis_val_0(0, 0));
    CHECK_CLOSE(exact_deriv_0, ad_deriv_0, 1e-10);
    
}

}
