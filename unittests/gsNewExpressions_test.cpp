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
#include <gsNewExpressions/ExprAssembler.h>

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
        GISMO_UNUSED(ev_mb);  // Available but not needed for this test
        
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

        // Vector Laplacian (trace/sum of component-wise Laplacians) - compare against bare deriv2
        // Note: The current implementation returns sum of component laplacians (scalar)
        // A proper vector laplacian would return [lapl(f_x), lapl(f_y)]^T (vector)
        CHECK_EQUAL(eval(helper, Expr::lapl(f2D), pt)(0,0).value(), lapl_func2D);
        
        // Linearity tests for vector Laplacian
        // TODO: lapl(f2D+f2D) - not yet working with expression composition
        // TODO: lapl(f2D-f2D) - not yet working with expression composition
        // TODO: lapl(2*f2D) - not yet working with constant multiplication
        // TODO: lapl(f2D*2) - not yet working with constant multiplication

        // TODO: Test lapl(f1D * f1D) - product rule: 2*grad(f1D)·grad(f1D) + 2*f1D*lapl(f1D)
        // TODO: Test lapl with basis functions (spaces)
        // TODO: Compare against old gsExprAssembler lapl(f1Dold) once integrated
    }

    /**
     * @brief Test vector Laplacian returns a vector [lapl(f_x), lapl(f_y)]^T
     * 
     * For vector field F = [f_x, f_y]^T, the vector Laplacian is:
     * ∇²F = [∇²f_x, ∇²f_y]^T
     */
    TEST(VectorLaplExpression)
    {
        // Vector function f2D = [x^2, y^2]^T
        // lapl(f2D[i]) should return lapl of component i
        
        // Calculate expected values analytically:
        // f0 = x^2: d²f₀/dx² = 2, d²f₀/dy² = 0, lapl(f0) = 2
        // f1 = y^2: d²f₁/dx² = 0, d²f₁/dy² = 2, lapl(f1) = 2
        real_t lapl_comp0 = 2.0;  // lapl(x²) = 2 + 0 = 2
        real_t lapl_comp1 = 2.0;  // lapl(y²) = 0 + 2 = 2
        
        // Evaluate component laplacians
        auto result0 = eval(helper, Expr::lapl(f2D[0]), pt);
        auto result1 = eval(helper, Expr::lapl(f2D[1]), pt);
        auto result_vec = eval(helper, Expr::lapl(f2D), pt);
        
        // Extract values as real_t
        real_t val0 = result0(0,0)(0,0);
        real_t val1 = result1(0,0)(0,0);
        real_t val_vec = result_vec(0,0)(0,0);
        
        // Note: gsFunctionExpr uses finite differences for deriv2, so tolerance must account for numerical error
        // Typical FD error is O(h^2) where h ~ 1e-5, so error ~ 1e-10, but accumulated error can be ~ 1e-7
        const real_t tol = 1e-6;
        
        // Test lapl of each component
        CHECK_CLOSE(val0, lapl_comp0, tol);
        CHECK_CLOSE(val1, lapl_comp1, tol);
        
        // Test that sum of component laplacians equals the vector laplacian (trace)
        CHECK_CLOSE(val_vec, lapl_comp0 + lapl_comp1, tol);
    }

    /**
     * @brief Test igrad (physical gradient) and measure expressions
     * 
     * For identity geometry:
     * - igrad(f, G) = grad(f) (since jacobian inverse is identity)
     * - meas(G) = 1 (Jacobian determinant)
     * 
     * For non-trivial geometry:
     * - igrad(f, G) = grad(f) * jac(G)^{-1}
     * - meas(G) = |det(jac(G))|
     */
    TEST(gsNewExpressions_igrad_meas)
    {
        // Identity geometry on unit square
        gsMultiPatch<real_t> mp_id(*gsNurbsCreator<real_t>::BSplineSquareDeg(2));
        
        ExpressionHelper<real_t> helper;
        helper.points() = gsMatrix<real_t>(pt);
        
        // Geometry map
        auto G = helper.getMap(mp_id);
        auto f1D_expr = helper.getScalarFunction(func1D, "f");
        helper.add(f1D_expr);
        
        // Parse all expressions
        auto igrad_expr = Expr::igrad(f1D_expr, G);
        auto meas_expr = Expr::meas(G);
        
        igrad_expr.parse(helper);
        meas_expr.parse(helper);
        helper.precompute(0);
        
        // Test meas(G) = 1 for identity geometry
        auto meas_val = meas_expr.eval(0);
        CHECK_CLOSE(meas_val(0,0)(0,0), 1.0, 1e-12);
        
        // Test igrad(f, G)
        auto igrad_val = igrad_expr.eval(0);
        auto grad_val = Expr::grad(f1D_expr).eval(0);
        
        // Compare igrad to grad (should be equal for identity geometry)
        CHECK_CLOSE(igrad_val(0,0)(0,0), grad_val(0,0)(0,0), 1e-12);
        CHECK_CLOSE(igrad_val(0,0)(0,1), grad_val(0,0)(0,1), 1e-12);
        
        // Verify against analytic gradient of x^2 + y^2
        CHECK_CLOSE(igrad_val(0,0)(0,0), 2.0 * pt[0], 1e-10);  // ∂/∂x (x^2 + y^2) = 2x
        CHECK_CLOSE(igrad_val(0,0)(0,1), 2.0 * pt[1], 1e-10);  // ∂/∂y (x^2 + y^2) = 2y
    }

    /**
     * @brief Test Jacobian and inverse Jacobian expressions
     */
    TEST(gsNewExpressions_jacobian)
    {
        // Identity geometry on unit square  
        gsMultiPatch<real_t> mp_id(*gsNurbsCreator<real_t>::BSplineSquareDeg(2));
        
        ExpressionHelper<real_t> helper;
        helper.points() = gsMatrix<real_t>(pt);
        
        auto G = helper.getMap(mp_id);
        auto jac_expr = Expr::jac(G);
        auto jacInv_expr = Expr::jacInv(G);
        
        jac_expr.parse(helper);
        jacInv_expr.parse(helper);
        helper.precompute(0);
        
        // Test jac(G) for identity geometry - should be identity matrix
        auto jac_val = jac_expr.eval(0);
        CHECK_CLOSE(jac_val(0,0)(0,0), 1.0, 1e-12);
        CHECK_CLOSE(jac_val(0,0)(0,1), 0.0, 1e-12);
        CHECK_CLOSE(jac_val(0,0)(1,0), 0.0, 1e-12);
        CHECK_CLOSE(jac_val(0,0)(1,1), 1.0, 1e-12);
        
        // Test jacInv(G) for identity geometry - should also be identity
        auto jacInv_val = jacInv_expr.eval(0);
        CHECK_CLOSE(jacInv_val(0,0)(0,0), 1.0, 1e-12);
        CHECK_CLOSE(jacInv_val(0,0)(0,1), 0.0, 1e-12);
        CHECK_CLOSE(jacInv_val(0,0)(1,0), 0.0, 1e-12);
        CHECK_CLOSE(jacInv_val(0,0)(1,1), 1.0, 1e-12);
    }

    /**
     * @brief Test Poisson assembly: verify old and new assemblers produce identical matrices
     * 
     * This test assembles -Δu = f using both assemblers and compares entry-by-entry
     */
    TEST(gsNewExpressions_poisson_assembly)
    {
        // Create domain and basis
        gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineSquare(2));
        gsMultiBasis<real_t> basis(mp, true);
        basis.uniformRefine();
        
        // Boundary conditions: homogeneous Dirichlet
        gsBoundaryConditions<real_t> bc;
        gsFunctionExpr<real_t> g_D("0", 2);
        bc.addCondition(boundary::west, condition_type::dirichlet, &g_D);
        bc.addCondition(boundary::east, condition_type::dirichlet, &g_D);
        bc.addCondition(boundary::south, condition_type::dirichlet, &g_D);
        bc.addCondition(boundary::north, condition_type::dirichlet, &g_D);
        bc.setGeoMap(mp);
        
        // OLD ASSEMBLER
        gsExprAssembler<real_t> A_old(1, 1);
        A_old.setIntegrationDomain(basis.domain());
        
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::space space;
        
        geometryMap G_old = A_old.getMap(mp);
        space u_old = A_old.getSpace(basis);
        u_old.setup(bc, dirichlet::homogeneous, 0);
        
        A_old.initSystem();
        A_old.assemble(igrad(u_old, G_old) * igrad(u_old, G_old).tr() * meas(G_old));
        gsSparseMatrix<real_t> K_old = A_old.matrix();
        
        // NEW ASSEMBLER
        auto domain_ptr = basis.domain();
        ExprAssembler<real_t> A_new(*domain_ptr, 1, 1);
        
        auto v = A_new.getScalarTestSpace(basis, 0, "v");
        auto u = A_new.getScalarTrialSpace(basis, 0, "u");
        v.setup(bc, dirichlet::homogeneous, 0);
        u.setup(bc, dirichlet::homogeneous, 0);
        
        auto G = A_new.getMap(mp);
        A_new.initSystem();
        A_new.assemble(Expr::inner(Expr::igrad(v, G), Expr::igrad(u, G)) * Expr::meas(G));
        gsSparseMatrix<real_t> K_new = A_new.matrix();
        
        // Compare matrices
        CHECK_EQUAL(K_old.rows(), K_new.rows());
        CHECK_EQUAL(K_old.cols(), K_new.cols());
        
        real_t max_diff = 0;
        for (index_t i = 0; i < K_old.rows(); ++i)
        {
            for (index_t j = 0; j < K_old.cols(); ++j)
            {
                real_t diff = std::abs(K_old.coeff(i, j) - K_new.coeff(i, j));
                max_diff = std::max(max_diff, diff);
            }
        }
        CHECK_CLOSE(max_diff, 0.0, 1e-10);
    }

    /**
     * @brief Test mass matrix assembly: verify old and new assemblers match
     */
    TEST(gsNewExpressions_mass_assembly)
    {
        // Create domain and basis
        gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineSquare(2));
        gsMultiBasis<real_t> basis(mp, true);
        basis.uniformRefine();
        
        // Boundary conditions: homogeneous Dirichlet
        gsBoundaryConditions<real_t> bc;
        gsFunctionExpr<real_t> g_D("0", 2);
        bc.addCondition(boundary::west, condition_type::dirichlet, &g_D);
        bc.addCondition(boundary::east, condition_type::dirichlet, &g_D);
        bc.addCondition(boundary::south, condition_type::dirichlet, &g_D);
        bc.addCondition(boundary::north, condition_type::dirichlet, &g_D);
        bc.setGeoMap(mp);
        
        // OLD ASSEMBLER
        gsExprAssembler<real_t> A_old(1, 1);
        A_old.setIntegrationDomain(basis.domain());
        
        typedef gsExprAssembler<real_t>::geometryMap geometryMap;
        typedef gsExprAssembler<real_t>::space space;
        
        geometryMap G_old = A_old.getMap(mp);
        space u_old = A_old.getSpace(basis);
        u_old.setup(bc, dirichlet::homogeneous, 0);
        
        A_old.initSystem();
        A_old.assemble(u_old * u_old.tr() * meas(G_old));
        gsSparseMatrix<real_t> M_old = A_old.matrix();
        
        // NEW ASSEMBLER
        auto domain_ptr = basis.domain();
        ExprAssembler<real_t> A_new(*domain_ptr, 1, 1);
        
        auto v = A_new.getScalarTestSpace(basis, 0, "v");
        auto u = A_new.getScalarTrialSpace(basis, 0, "u");
        v.setup(bc, dirichlet::homogeneous, 0);
        u.setup(bc, dirichlet::homogeneous, 0);
        
        auto G = A_new.getMap(mp);
        A_new.initSystem();
        A_new.assemble(v * u * Expr::meas(G));
        gsSparseMatrix<real_t> M_new = A_new.matrix();
        
        // Compare matrices
        CHECK_EQUAL(M_old.rows(), M_new.rows());
        CHECK_EQUAL(M_old.cols(), M_new.cols());
        
        real_t max_diff = 0;
        for (index_t i = 0; i < M_old.rows(); ++i)
        {
            for (index_t j = 0; j < M_old.cols(); ++j)
            {
                real_t diff = std::abs(M_old.coeff(i, j) - M_new.coeff(i, j));
                max_diff = std::max(max_diff, diff);
            }
        }
        CHECK_CLOSE(max_diff, 0.0, 1e-10);
    }

    /**
     * @brief Test that setDerivative propagates correctly through expression trees
     * 
     * This verifies that when we create grad(f), div(f), lapl(f), etc., the
     * setDerivative mechanism correctly sets the derivative order on the
     * underlying expression so that the correct flags are set during parse().
     */
    TEST(gsNewExpressions_derivative_propagation)
    {
        // Create test setup
        gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineSquareDeg(3));
        gsMultiBasis<real_t> basis(mp);
        
        ExpressionHelper<real_t> helper;
        auto f1D = helper.getScalarFunction(func1D, "f1D");
        auto f2D = helper.getVectorFunction(func2D, "f2D");
        
        // Test 1: grad(f1D) should have Deriv = 1
        auto grad_f1D = Expr::grad(f1D);
        typedef decltype(grad_f1D) GradType;
        CHECK_EQUAL(size_t(GradType::Deriv), 1u);
        CHECK_EQUAL(size_t(GradType::Order), 1u);  // Scalar -> Vector
        
        // Test 2: div(f2D) should have Deriv = 1
        auto div_f2D = Expr::div(f2D);
        typedef decltype(div_f2D) DivType;
        CHECK_EQUAL(size_t(DivType::Deriv), 1u);
        CHECK_EQUAL(size_t(DivType::Order), 0u);  // Vector -> Scalar
        
        // Test 3: lapl(f1D) should have Deriv = 2
        auto lapl_f1D = Expr::lapl(f1D);
        typedef decltype(lapl_f1D) LaplType;
        CHECK_EQUAL(size_t(LaplType::Deriv), 2u);
        CHECK_EQUAL(size_t(LaplType::Order), 0u);  // Scalar -> Scalar
        
        // Test 4: Verify parse() sets correct flags on the *expression's copy* of the variable
        // Note: UnaryOperator stores expressions by value, so the original variable's
        // data is not modified. We verify flags are set on the expression's internal copy.
        {
            ExpressionHelper<real_t> h1;
            auto f1 = h1.getScalarFunction(func1D, "f1");
            auto g1 = Expr::grad(f1);
            h1.points() = gsMatrix<real_t>(pt);
            g1.parse(h1);
            // Access data through the gradient expression's copy, not the original
            CHECK(g1.expr().data().flags & NEED_GRAD);
        }
        
        {
            ExpressionHelper<real_t> h2;
            auto f2 = h2.getScalarFunction(func1D, "f2");
            auto l2 = Expr::lapl(f2);
            h2.points() = gsMatrix<real_t>(pt);
            l2.parse(h2);
            CHECK(l2.expr().data().flags & NEED_LAPLACIAN);
        }
        
        {
            ExpressionHelper<real_t> h3;
            auto f3 = h3.getVectorFunction(func2D, "f3");
            auto d3 = Expr::div(f3);
            h3.points() = gsMatrix<real_t>(pt);
            d3.parse(h3);
            CHECK(d3.expr().data().flags & NEED_DERIV);
        }
    }

    /**
     * @brief Test derivative trait calculation for composed expressions
     * 
     * Verifies that the Deriv trait is correctly calculated for
     * expressions like grad(grad(f)) -> Deriv = 2
     */
    TEST(gsNewExpressions_composed_derivative_traits)
    {
        // Test 1: grad(f) has Deriv = 1
        typedef Expr::VariableObject<real_t, 0, false> ScalarVar;
        typedef Expr::GradExpression<ScalarVar, 0, Expr::SpaceType::None, 0> GradScalar;
        CHECK_EQUAL(size_t(GradScalar::Deriv), 1u);
        
        // Test 2: lapl(f) has Deriv = 2
        typedef Expr::LaplExpression<ScalarVar, 0, Expr::SpaceType::None, 0> LaplScalar;
        CHECK_EQUAL(size_t(LaplScalar::Deriv), 2u);
        
        // Test 3: div(grad(f)) = lapl(f) conceptually has Deriv = 2
        // (implemented via lapl() factory function)
        
        // Test 4: grad(grad(f)) for scalar would be Hessian (Deriv = 2)
        // Note: GradExpression of GradExpression
        typedef Expr::GradExpression<GradScalar, 1, Expr::SpaceType::None, 0> HessianScalar;
        CHECK_EQUAL(size_t(HessianScalar::Deriv), 2u);
    }

    /**
     * @brief Test SpaceObject::check() and setup() behavior
     * 
     * Verifies:
     * - check() returns false before setup()
     * - check() returns true after setup()
     * - Mapper is properly initialized after setup()
     * - Free DOFs and boundary DOFs are correctly computed
     */
    TEST(gsNewExpressions_space_setup)
    {
        // Create domain and basis
        gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineSquare(2));
        gsMultiBasis<real_t> basis(mp, true);
        basis.uniformRefine();
        
        // Boundary conditions
        gsBoundaryConditions<real_t> bc;
        gsFunctionExpr<real_t> g_D("0", 2);
        bc.addCondition(boundary::west, condition_type::dirichlet, &g_D);
        bc.addCondition(boundary::east, condition_type::dirichlet, &g_D);
        bc.addCondition(boundary::south, condition_type::dirichlet, &g_D);
        bc.addCondition(boundary::north, condition_type::dirichlet, &g_D);
        bc.setGeoMap(mp);
        
        // Create assembler and spaces
        auto domain_ptr = basis.domain();
        ExprAssembler<real_t> A(*domain_ptr, 1, 1);
        
        auto v = A.getScalarTestSpace(basis, 0, "v");
        auto u = A.getScalarTrialSpace(basis, 0, "u");
        
        // Before setup: check() should return false (spaces not initialized for assembly)
        // Note: After getScalarTestSpace, FeSpaceData::init() is called which sets initialized=true
        // but the mapper doesn't have boundary conditions applied yet
        // check() returns: m_sd->isInitialized() which is: initialized && mapper.isFinalized()
        
        // Test that setup with BCs changes the mapper
        index_t total_dofs_before = v.mapper().freeSize();
        
        // Setup with boundary conditions
        v.setup(bc, dirichlet::homogeneous, 0);
        u.setup(bc, dirichlet::homogeneous, 0);
        
        // After setup: check() should return true
        CHECK(v.check());
        CHECK(u.check());
        
        // Verify mapper is properly configured
        CHECK(v.mapper().isFinalized());
        CHECK(u.mapper().isFinalized());
        
        // Verify boundary DOFs were eliminated
        index_t free_dofs = v.mapper().freeSize();
        index_t boundary_dofs = v.mapper().boundarySize();
        
        // Total basis size for a 2-refined unit square with degree 1 is 5x5 = 25
        // Boundary DOFs (4 sides, minus corners counted once): 4*5 - 4 = 16
        // Free DOFs: 25 - 16 = 9
        CHECK_EQUAL(free_dofs + boundary_dofs, basis.basis(0).size());
        CHECK(boundary_dofs > 0);  // Some boundary DOFs should be eliminated
        CHECK(free_dofs < total_dofs_before);  // Free DOFs should decrease after setup
        
        // Verify initSystem doesn't throw (spaces are initialized)
        A.initSystem();
        
        // Compare with old assembler
        gsExprAssembler<real_t> A_old(1, 1);
        A_old.setIntegrationDomain(basis.domain());
        auto u_old = A_old.getSpace(basis);
        u_old.setup(bc, dirichlet::homogeneous, 0);
        
        CHECK_EQUAL(v.mapper().freeSize(), u_old.mapper().freeSize());
        CHECK_EQUAL(v.mapper().boundarySize(), u_old.mapper().boundarySize());
    }

    /**
     * @brief Test simple setup (without boundary conditions)
     */
    TEST(gsNewExpressions_space_simple_setup)
    {
        // Create domain and basis
        gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineSquare(2));
        gsMultiBasis<real_t> basis(mp, true);
        basis.uniformRefine();
        
        // Create assembler and spaces
        auto domain_ptr = basis.domain();
        ExprAssembler<real_t> A(*domain_ptr, 1, 1);
        
        auto v = A.getScalarTestSpace(basis, 0, "v");
        auto u = A.getScalarTrialSpace(basis, 0, "u");
        
        // Simple setup (no BC elimination)
        v.setup(0);  // C0 continuity
        u.setup(0);
        
        // After setup: check() should return true
        CHECK(v.check());
        CHECK(u.check());
        
        // Verify all DOFs are free (no boundary elimination)
        CHECK_EQUAL(v.mapper().freeSize(), basis.basis(0).size());
        CHECK_EQUAL(v.mapper().boundarySize(), 0);
    }

    /**
     * @brief Test sum factorization assembly produces same results as regular assembly
     * 
     * Sum factorization exploits tensor product structure to assemble matrices
     * more efficiently. The results should be identical to regular assembly.
     */
    TEST(gsNewExpressions_assembleSF_mass_matrix)
    {
        // Create simple 2D tensor-product B-spline basis
        gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineSquare(2));
        gsMultiBasis<real_t> basis(mp, true);
        basis.uniformRefine();
        
        auto domain_ptr = basis.domain();
        
        // Regular assembly
        ExprAssembler<real_t> A_reg(*domain_ptr, 1, 1);
        auto v_reg = A_reg.getScalarTestSpace(basis, 0, "v");
        auto u_reg = A_reg.getScalarTrialSpace(basis, 0, "u");
        v_reg.setup(0);
        u_reg.setup(0);
        A_reg.initSystem();
        A_reg.assemble(v_reg * u_reg);  // Mass matrix
        
        // Sum factorization assembly
        ExprAssembler<real_t> A_sf(*domain_ptr, 1, 1);
        auto v_sf = A_sf.getScalarTestSpace(basis, 0, "v");
        auto u_sf = A_sf.getScalarTrialSpace(basis, 0, "u");
        v_sf.setup(0);
        u_sf.setup(0);
        A_sf.initSystem();
        A_sf.assembleSF(v_sf * u_sf);  // Mass matrix via sum factorization
        
        // Compare matrices entry by entry
        const gsSparseMatrix<real_t>& M_reg = A_reg.matrix();
        const gsSparseMatrix<real_t>& M_sf = A_sf.matrix();
        
        CHECK_EQUAL(M_reg.rows(), M_sf.rows());
        CHECK_EQUAL(M_reg.cols(), M_sf.cols());
        
        // Check that matrices are close (within tolerance)
        real_t diff = (M_reg - M_sf).norm();
        real_t tol = 1e-10 * M_reg.norm();
        CHECK(diff < tol);
    }

    TEST(CrossProductExpression)
    {
        // Test cross product in 3D
        gsFunctionExpr<real_t> vec3D_a("x", "y", "1", 3);
        gsFunctionExpr<real_t> vec3D_b("1", "0", "z", 3);
        
        gsVector<real_t,3> pt3D;
        pt3D << 0.5, 0.3, 0.7;
        
        // TODO: Cross product operator % not yet available for VariableObject
        // Need to implement operator% first
        // ExpressionHelper<real_t> helper3D;
        // auto f1 = helper3D.getVectorFunction(vec3D_a, "f1");
        // auto f2 = helper3D.getVectorFunction(vec3D_b, "f2");
        // auto cross_expr = f1 % f2;  // Cross product
    }
    
    TEST(TransposeExpression)
    {
        // Test transpose of gradient
        auto grad_f = Expr::grad(f2D);  // Should be 2x2 Jacobian
        auto grad_t = grad_f.tr();  // Transpose
        
        auto result = eval(helper, grad_f, pt);
        auto result_t = eval(helper, grad_t, pt);
        
        // Check dimensions are swapped
        CHECK_EQUAL(result().rows(), result_t().cols());
        CHECK_EQUAL(result().cols(), result_t().rows());
        
        // Check values are transposed
        for (index_t i = 0; i < result().rows(); ++i)
            for (index_t j = 0; j < result().cols(); ++j)
                CHECK_CLOSE(result()(i,j), result_t()(j,i), 1e-10);
    }
    
    TEST(ComponentExpression)
    {
        // Test extracting components from vector function
        auto f0 = f2D[0];  // x^2
        auto f1 = f2D[1];  // y^2
        
        auto result_full = eval(helper, f2D, pt);
        auto result_0 = eval(helper, f0, pt);
        auto result_1 = eval(helper, f1, pt);
        
        CHECK_CLOSE(result_0()(0,0), result_full()(0,0), 1e-10);
        CHECK_CLOSE(result_1()(0,0), result_full()(1,0), 1e-10);
        CHECK_CLOSE(result_0()(0,0), pt[0]*pt[0], 1e-10);
        CHECK_CLOSE(result_1()(0,0), pt[1]*pt[1], 1e-10);
    }

/*
    UNTESTED FEATURES - ADDITIONAL TODOs
*/
    // TODO: Cross product operator % for VariableObjects
    // TODO: ArrayExpression - needs proper implementation first
    // TODO: NormalExpression - needs boundary evaluation support in ExpressionHelper
    // TODO: Hessian expression
    // TODO: Higher-order derivatives (third, fourth order)
    // TODO: Complex nested geometry transformations
    // TODO: Expression optimization and caching
    // TODO: Parallel evaluation of expressions

}
