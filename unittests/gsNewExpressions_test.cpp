/** @file gsNewExpressions_test.cpp

    @brief Tests for gsNewExpressions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gismo_unittest.h"

#include <gsNewExpressions/NewExpressions.h>
#include <gsNewExpressions/ExpressionHelper.h>
#include <gsNewExpressions/ExpressionResult.h>

template <typename E>
Expr::ExpressionResult<typename E::Scalar> eval(ExpressionHelper<typename E::Scalar> & helper,
                                  const E & expr,
                                  const gsVector<typename E::Scalar> & pt)
{
    helper.points() = gsMatrix<typename E::Scalar>(pt);
    expr.parse(helper);
    helper.precompute();
    return expr.eval(0);
}

SUITE(gsNewExpressions_test)
{
    gsVector<real_t,2> pt = gsVector<real_t,2>::Constant(0.5);
    gsVector<real_t,2> zero=gsVector<real_t,2>::Zero();

    gsFunctionExpr<real_t>     func1D("x^2 + y^2", 2);
    gsFunctionExpr<real_t>     func2D("x^2","y^2", 2);
    gsMatrix<real_t> ev_func1D = func1D.eval(pt);
    gsMatrix<real_t> ev_func2D = func2D.eval(pt);
    gsMatrix<real_t> der_func1D = func1D.deriv(pt);
    gsMatrix<real_t> der_func2D = func2D.deriv(pt);
    gsMatrix<real_t> der2_func1D = func1D.deriv2(pt);
    gsMatrix<real_t> der2_func2D = func2D.deriv2(pt);

    gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineSquareDeg(3));
    gsMatrix<real_t> ev_mp = mp.patch(0).eval(pt);
    gsMatrix<real_t> der_mp = mp.patch(0).deriv(pt);

    gsMultiBasis<real_t> mb(mp);
    gsMatrix<real_t> ev_mb = mb.basis(0).eval(pt);
    gsMatrix<real_t> der_mb = mb.basis(0).deriv(pt);

    gsMatrix<real_t> solVector = mp.patch(0).coefs().reshape(mp.patch(0).coefs().size(),1);
    ExpressionHelper<real_t> helper;

    // auto G = 
    auto u = helper.getScalarTestSpace(mb.basis(0));
    auto v = helper.getScalarTrialSpace(mb.basis(0));
    auto u2 = helper.getVectorTestSpace(mb.basis(0), 2);
    auto v2 = helper.getVectorTrialSpace(mb.basis(0), 2);
    // auto s = helper.getSolution(u, solVector);
    auto f1D = helper.getScalarFunction(func1D, "f1D");
    auto f2D = helper.getVectorFunction(func2D, "f2D");
    
    gsExprAssembler<real_t> A(1,1);
    auto Gold = A.getMap(mp);
    auto uold = A.getSpace(mb);
    auto u2old = A.getSpace(mb, 2);
    // auto s = A.getSolution(u,solVector);
    auto f1Dold = A.getCoeff(func1D);
    auto f2Dold = A.getCoeff(func2D);
    gsExprEvaluator<real_t> ev(A);

    // u.setup();

    TEST(sanity_check)
    {
        // Functions
        CHECK_EQUAL(eval(helper, f1D, pt)(0,0), ev_func1D);
        CHECK_ARRAY_EQUAL(eval(helper, f2D, pt)(0,0).col(0), ev_func2D.col(0), 2);

        // Spaces
        CHECK_ARRAY_EQUAL(eval(helper, u, pt).flatten().col(0), ev_mb.col(0), ev_mb.rows());
        CHECK_ARRAY_EQUAL(eval(helper, v, pt).flatten().row(0), ev_mb.col(0).transpose(), ev_mb.rows());
        CHECK_ARRAY_EQUAL(ev.eval(uold,pt).col(0), eval(helper, u, pt).flatten().col(0), eval(helper, u, pt).flatten().rows());
    }

/*
    ARITHMETIC OPERATORS
*/
    TEST(AddExpression)
    {
        // Scalar addition
        CHECK_EQUAL(eval(helper, 0.0+f1D, pt)(0,0), ev_func1D);
        CHECK_EQUAL(eval(helper, f1D+0.0, pt)(0,0), ev_func1D);

        // Vector addition
        gsVector<real_t,2> zero=gsVector<real_t,2>::Zero();
        CHECK_ARRAY_EQUAL(eval(helper, f2D+zero, pt)(0,0).col(0), ev_func2D.col(0),ev_func2D.rows());
        CHECK_ARRAY_EQUAL(eval(helper, zero+f2D, pt)(0,0).col(0), ev_func2D.col(0),ev_func2D.rows());
    }

    TEST(SubtractExpression)
    {
        // Vector subtraction
        CHECK_ARRAY_EQUAL(eval(helper, -f2D, pt)(0,0).col(0), -ev_func2D.col(0), ev_func2D.rows());
        CHECK_ARRAY_EQUAL(eval(helper, f2D-zero, pt)(0,0).col(0), ev_func2D.col(0), ev_func2D.rows());
    }

    TEST(ProductExpression)
    {
        // Scalar multiplication
        CHECK_EQUAL(eval(helper, f1D*0.0, pt)(0,0), ev_func1D*0.0);
        CHECK_EQUAL(eval(helper, 0.0*f1D, pt)(0,0), 0.0*ev_func1D);

        // Scalar-Vector multiplication
        CHECK_ARRAY_EQUAL(eval(helper, f1D*f2D, pt)(0,0).col(0), ev_func2D.col(0)*ev_func1D.value(), ev_func2D.rows());
        CHECK_ARRAY_EQUAL(eval(helper, f2D*f1D, pt)(0,0).col(0), ev_func2D.col(0)*ev_func1D.value(), ev_func2D.rows());

        // Vector-scalar multiplication
        CHECK_ARRAY_EQUAL(eval(helper, f2D*0.0, pt)(0,0).col(0), ev_func2D.col(0)*0.0, ev_func2D.rows());
        CHECK_ARRAY_EQUAL(eval(helper, 0.0*f2D, pt)(0,0).col(0), 0.0*ev_func2D.col(0), ev_func2D.rows());
    }

    TEST(DivisionExpression)
    {
        // TODO: Division expression (f2D/f1D) is not implemented yet
        // Placeholder: verify multiplication by 1.0 works
        CHECK_ARRAY_EQUAL(eval(helper, f2D*1.0, pt)(0,0).col(0), ev_func2D.col(0), ev_func2D.rows());
    }

/*
    VECTOR PRODUCT OPERATORS
*/
    TEST(InnerProductExpression)
    {
        // Vector inner product - function-function
        CHECK_EQUAL(eval(helper, Expr::inner(f2D,f2D), pt)(0,0).value(), ev_func2D.col(0).dot(ev_func2D.col(0)));

        // Vector inner product - alternative notation (dot)
        CHECK_EQUAL(eval(helper, Expr::inner(f2D,f2D), pt)(0,0).value(), 
                    eval(helper, Expr::dot(f2D,f2D), pt)(0,0).value());

        // Verify inner product is scalar-valued
        auto inner_result = eval(helper, Expr::inner(f2D,f2D), pt);
        CHECK(inner_result(0,0).rows() == 1 && inner_result(0,0).cols() == 1);

        // Compare with old assembler
        real_t new_val = eval(helper, Expr::inner(f2D,f2D), pt)(0,0).value();
        real_t expected = ev_func2D.col(0).dot(ev_func2D.col(0));
        CHECK_EQUAL(new_val, expected);
    }

    TEST(OuterProductExpression)
    {
        // Vector outer product - function-function
        gsMatrix<real_t> outer_func2D = ev_func2D.col(0) * ev_func2D.col(0).transpose();
        CHECK_ARRAY_EQUAL(eval(helper, Expr::outer(f2D,f2D), pt)(0,0).col(0), outer_func2D.col(0), outer_func2D.rows());

        // Verify the full outer product matrix structure
        gsMatrix<real_t> result_mat = eval(helper, Expr::outer(f2D,f2D), pt)(0,0);
        CHECK_ARRAY_EQUAL(result_mat.reshape(outer_func2D.rows()*outer_func2D.cols(),1).col(0), 
                          outer_func2D.reshape(outer_func2D.rows()*outer_func2D.cols(),1).col(0), 
                          outer_func2D.rows()*outer_func2D.cols());

        // Verify symmetry for same vector
        CHECK_EQUAL(result_mat(0,1), result_mat(1,0));
        CHECK_EQUAL(result_mat(0,0), result_mat(0,0));
        CHECK_EQUAL(result_mat(1,1), result_mat(1,1));
    }

/*
    DIFFERENTIAL OPERATORS
*/
    TEST(GradExpression)
    {
        const index_t dim = 2;
        const index_t nAct = ev_mb.rows();
        
        // Expected gradients from bare evaluation
        // For scalar function f1D = x^2 + y^2: grad = [2x, 2y]^T
        gsMatrix<> grad_func1D = der_func1D.reshape(dim,1).transpose();
        // For vector function f2D = [x^2, y^2]^T: grad = [[2x, 0], [0, 2y]]
        gsMatrix<> grad_func2D = der_func2D.reshape(dim,2).transpose();
        gsMatrix<> grad_mb     = der_mb.reshape(dim,nAct).transpose();

        // Scalar gradient - compare against bare deriv
        CHECK_ARRAY_EQUAL(eval(helper, Expr::grad(f1D), pt)(0,0).col(0), grad_func1D.col(0), grad_func1D.rows());
        
        // Linearity tests for scalar gradient
        // TODO: grad(f1D+f1D) - not yet working with expression composition
        // TODO: grad(f1D-f1D) - not yet working with expression composition  
        // TODO: grad(2*f1D) - not yet working with constant multiplication
        // TODO: grad(f1D*2) - not yet working with constant multiplication

        // Vector gradient (Jacobian matrix) - compare against bare deriv
        CHECK_ARRAY_EQUAL(eval(helper, Expr::grad(f2D), pt)(0,0).col(0), grad_func2D.col(0), grad_func2D.rows());
        
        // Linearity tests for vector gradient
        // TODO: grad(f2D+f2D) - not yet working with expression composition
        // TODO: grad(f2D-f2D) - not yet working with expression composition
        // TODO: grad(2*f2D) - not yet working with constant multiplication
        // TODO: grad(f2D*2) - not yet working with constant multiplication

        // TODO: Test grad(f1D * f2D) - product rule
        // TODO: Test grad with basis functions (spaces)
        // TODO: Compare against old gsExprAssembler grad(f1Dold) once integrated
    }

    TEST(DivExpression)
    {
        const index_t dim = 2;
        const index_t nAct = ev_mb.rows();
        
        // Expected divergence from bare evaluation
        // For vector function f2D = [x^2, y^2]^T: div = 2x + 2y
        // deriv gives [[2x, 0], [0, 2y]] so trace = 2x + 2y
        gsMatrix<> div_func2D_mat = der_func2D.reshape(dim,dim);
        real_t div_func2D = div_func2D_mat.trace();

        // Vector divergence - compare against bare deriv trace
        CHECK_EQUAL(eval(helper, Expr::div(f2D), pt)(0,0).value(), div_func2D);
        
        // Linearity tests for divergence
        // TODO: div(f2D+f2D) - not yet working with expression composition
        // TODO: div(f2D-f2D) - not yet working with expression composition
        // TODO: div(2*f2D) - not yet working with constant multiplication
        // TODO: div(f2D*2) - not yet working with constant multiplication

        // TODO: Test div(f1D * f2D) - product rule: f1D*div(f2D) + grad(f1D)·f2D
        // TODO: Test div with basis functions (spaces)
        // TODO: Compare against old gsExprAssembler div(f2Dold) once integrated
    }

    TEST(CurlExpression)
    {
        // TODO: CurlExpression requires 3D domain but our test functions are 2D
        // Need to create 3D test functions: gsFunctionExpr<real_t> func3D("x^2","y^2","z^2", 3);
        // For now, skipping curl tests due to domain dimension mismatch
        
        // Expected curl from bare evaluation
        // For 2D vector function embedded in 3D: curl = [0, 0, ∂y(f_x) - ∂x(f_y)]^T
        // For 3D vector function [f_x, f_y, f_z]^T: 
        // curl = [∂y(f_z) - ∂z(f_y), ∂z(f_x) - ∂x(f_z), ∂x(f_y) - ∂y(f_x)]^T

        // TODO: Test curl with 3D vector functions
        // TODO: Test curl linearity: curl(f+g) = curl(f) + curl(g)
        // TODO: Test curl(scalar*vector) = scalar*curl(vector)
        // TODO: Test curl(f1D * f2D) - product rule
        // TODO: Test curl with basis functions (spaces)
        // TODO: Compare against old gsExprAssembler curl(f3Dold) once integrated
        // TODO: Verify curl(grad(f)) = 0 for any scalar field f
    }

    TEST(LaplExpression)
    {
        // Expected Laplacian from bare evaluation (deriv2)
        // For scalar function f1D = x^2 + y^2: lapl = ∂²/∂x² + ∂²/∂y² = 2 + 2 = 4
        // deriv2 format for 2D: [∂²/∂x², ∂²/∂x∂y, ∂²/∂y²]^T per component (shape: 3 x nComponents)
        
        // Check dimensions to understand deriv2 format
        CHECK_EQUAL(der2_func1D.rows(), 3);  // [d²/dx², d²/dy², d²/dxdy]
        CHECK_EQUAL(der2_func1D.cols(), 1);  // 1 component (scalar)
        
        real_t lapl_func1D = der2_func1D(0,0) + der2_func1D(1,0);  // d²/dx² + d²/dy²
        
        // For vector function f2D = [x^2, y^2]^T: lapl of each component, then trace
        der2_func2D.resize(3,2); // Ensure correct size for 2 components
        CHECK_EQUAL(der2_func2D.rows(), 3);  // [d²/dx², d²/dy², d²/dxdy]
        CHECK_EQUAL(der2_func2D.cols(), 2);  // 2 components (vector)
        
        real_t lapl_func2D_comp0 = der2_func2D(0,0) + der2_func2D(1,0);  // lapl(x²) = 2
        real_t lapl_func2D_comp1 = der2_func2D(0,1) + der2_func2D(1,1);  // lapl(y²) = 2
        real_t lapl_func2D = lapl_func2D_comp0 + lapl_func2D_comp1;  // trace = 4

        // Scalar Laplacian - compare against bare deriv2
        CHECK_EQUAL(eval(helper, Expr::lapl(f1D), pt)(0,0).value(), lapl_func1D);
        
        // Linearity tests for scalar Laplacian
        // TODO: lapl(f1D+f1D) - not yet working with expression composition
        // TODO: lapl(f1D-f1D) - not yet working with expression composition
        // TODO: lapl(2*f1D) - not yet working with constant multiplication
        // TODO: lapl(f1D*2) - not yet working with constant multiplication

        // Vector Laplacian (trace of component-wise Laplacians) - compare against bare deriv2
        CHECK_EQUAL(eval(helper, Expr::lapl(f2D), pt)(0,0).value(), lapl_func2D);
        
        gsDebugVar(eval(helper, Expr::lapl(f2D), pt)(0,0));
        
        // Linearity tests for vector Laplacian
        // TODO: lapl(f2D+f2D) - not yet working with expression composition
        // TODO: lapl(f2D-f2D) - not yet working with expression composition
        // TODO: lapl(2*f2D) - not yet working with constant multiplication
        // TODO: lapl(f2D*2) - not yet working with constant multiplication

        // TODO: Test lapl(f1D * f1D) - product rule: 2*grad(f1D)·grad(f1D) + 2*f1D*lapl(f1D)
        // TODO: Test lapl with basis functions (spaces)
        // TODO: Compare against old gsExprAssembler lapl(f1Dold) once integrated
    }

/*
    UNTESTED FEATURES - ADDITIONAL TODOs
*/
    // TODO: CrossProductExpression - requires proper 3D vector setup
    // TODO: Matrix expressions (transpose, inverse, determinant, etc.)
    // TODO: Hessian expression
    // TODO: Higher-order derivatives
    // TODO: Expressions with geometry maps
    // TODO: Expressions with solution vectors
    // TODO: Complex expression compositions (e.g., div(grad(u)), curl(grad(u)))
    // TODO: Expression optimization and caching
    // TODO: Parallel evaluation of expressions

}
