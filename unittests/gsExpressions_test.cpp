/** @file gsExpressions_test.cpp

    @brief Tests for gsExpressions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gismo_unittest.h"


/* TODO:
    -[ ] Add gsFeSolution
*/

SUITE(gsExpressions_test)
{
    gsVector<real_t,2> pt = gsVector<real_t,2>::Constant(0.5);
    gsVector<real_t,2> zero=gsVector<real_t,2>::Zero();

    gsFunctionExpr<real_t>     func1D("x^2 + y^2", 2);
    gsFunctionExpr<real_t>     func2D("x^2","y^2", 2);
    gsMatrix<real_t> ev_func1D = func1D.eval(pt);
    gsMatrix<real_t> ev_func2D = func2D.eval(pt);
    gsMatrix<real_t> der_func1D = func1D.deriv(pt);
    gsMatrix<real_t> der_func2D = func2D.deriv(pt);

    gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineSquare());
    gsMatrix<real_t> ev_mp = mp.patch(0).eval(pt);
    gsMatrix<real_t> der_mp = mp.patch(0).deriv(pt);

    gsMultiBasis<real_t> mb(mp);
    gsMatrix<real_t> ev_mb = mb.basis(0).eval(pt);
    gsMatrix<real_t> der_mb = mb.basis(0).deriv(pt);

    gsMatrix<real_t> solVector = mp.patch(0).coefs().reshape(mp.patch(0).coefs().size(),1);
    gsExprAssembler<real_t> A(1,1);
    auto G = A.getMap(mp);
    auto u = A.getSpace(mb);
    // auto s = A.getSolution(u,solVector);
    auto f1D = A.getCoeff(func1D);
    auto f2D = A.getCoeff(func2D);
    gsExprEvaluator<real_t> ev(A);

    // u.setup();

    TEST(sanity_check)
    {
        CHECK_EQUAL(ev.eval(G,pt), ev_mp);
        CHECK_EQUAL(ev.eval(f1D,pt), ev_func1D);
        CHECK_ARRAY_EQUAL(ev.eval(f2D,pt).col(0), ev_func2D.col(0), 2);
    }
/*
    OPERATORS
*/
    TEST(add_expr)
    {
        // Scalar addition
        CHECK_EQUAL(ev.eval(f1D.val()+0.0,pt), ev_func1D);
        CHECK_EQUAL(ev.eval(0.0+f1D.val(),pt), ev_func1D);

        // Vector addition
        CHECK_EQUAL(ev.eval(f2D+G,pt), ev_func2D+ev_mp);
        CHECK_EQUAL(ev.eval(G+f2D,pt), ev_func2D+ev_mp);

        // NOT WORKING:
        // CHECK_EQUAL(ev.eval(f1D+0.0,pt), pt);
        // CHECK_EQUAL(ev.eval(f2D+zero,pt), pt);
        // CHECK_EQUAL(ev.eval(u+0.0,pt), ev_func1D);
        // CHECK_EQUAL(ev.eval(0.0+u,pt), ev_func1D);
    }

    TEST(sub_expr)
    {
        // Scalar subtraction
        CHECK_EQUAL(ev.eval(f1D.val()-0.0,pt), ev_func1D);
        CHECK_EQUAL(ev.eval(0.0-f1D.val(),pt), -ev_func1D);

        // Vector subtraction
        CHECK_EQUAL(ev.eval(f2D-G,pt), ev_func2D-ev_mp);
        CHECK_EQUAL(ev.eval(G-f2D,pt), ev_mp-ev_func2D);

        // NOT WORKING:
        // CHECK_EQUAL(ev.eval(f1D-0.0,pt), pt);
        // CHECK_EQUAL(ev.eval(f2D-zero,pt), pt);
        // CHECK_EQUAL(ev.eval(u-0.0,pt), ev_func1D);
        // CHECK_EQUAL(ev.eval(0.0-u,pt), -ev_func1D);
    }

    TEST(mult_expr)
    {
        // Scalar multiplication
        CHECK_EQUAL(ev.eval(f1D.val()*0.0,pt), ev_func1D*0.0);
        CHECK_EQUAL(ev.eval(0.0*f1D.val(),pt), 0.0*ev_func1D);

        // Scalar-Vector multiplication
        CHECK_EQUAL((ev.eval(f1D.val()*f2D,pt)-ev_func2D*ev_func1D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(f2D*f1D.val(),pt)-ev_func2D*ev_func1D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(f1D.val()*u,pt)-ev_func1D.value()*ev_mb).norm(), 0.0);
        CHECK_EQUAL((ev.eval(u*f1D.val(),pt)-ev_mb*ev_func1D.value()).norm(), 0.0);

        CHECK_EQUAL((ev.eval(f2D*0.0,pt)-ev_func2D*0.0).norm(), 0.0);
        CHECK_EQUAL((ev.eval(0.0*f2D,pt)-0.0*ev_func2D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(u*0.0,pt)-ev_mb*0.0).norm(), 0.0);
        CHECK_EQUAL((ev.eval(0.0*u,pt)-0.0*ev_mb).norm(), 0.0);

        CHECK_EQUAL((ev.eval(G*0.0,pt)-ev_mp*0.0).norm(), 0.0);
        CHECK_EQUAL((ev.eval(0.0*G,pt)-0.0*ev_mp).norm(), 0.0);
    }

    TEST(div_expr)
    {
        // Scalar division
        CHECK_EQUAL(ev.eval(f1D.val()/1.0,pt).value(), ev_func1D.value()/1.0);
        CHECK_EQUAL(ev.eval(1.0/f1D.val(),pt).value(), 1.0/ev_func1D.value());

        CHECK_EQUAL((ev.eval(u/1.0,pt)-ev_mb/1.0).norm(), 0.0);
        // CHECK_EQUAL((ev.eval(1.0/u,pt)-1.0*ev_mb).norm(), 0.0);

        // Vector division
        CHECK_EQUAL((ev.eval(f2D/1.0,pt)-ev_func2D/1.0).norm(), 0.0);
        CHECK_EQUAL((ev.eval(f2D/f1D.val(),pt)-ev_func2D/ev_func1D.value()).norm(), 0.0);
    }

/*
    MATHEMATICAL FUNCTIONS
*/
    TEST(pow_expr)
    {
        // Scalar power
        CHECK_EQUAL(ev.eval(pow(f1D,2),pt).value(), math::pow(ev_func1D.value(),2.0));
        // CHECK_EQUAL((ev.eval(pow(u,2),pt)-ev_mb.cwiseProduct(ev_mb)).norm(), 0.0);
    }

/*
    DIFFERENTIAL OPERATORS
*/
    TEST(grad_expr)
    {
        const index_t dim = 2;
        const index_t nAct = ev_mb.rows();
        gsMatrix<> grad_func1D = der_func1D.reshape(dim,1).transpose();
        gsMatrix<> grad_func2D = der_func2D.reshape(dim,2).transpose();
        gsMatrix<> grad_mb     = der_mb.reshape(dim,nAct).transpose();

        // Scalar gradient
        CHECK_EQUAL((ev.eval(grad(f1D),pt)    -grad_func1D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(f1D+f1D),pt)-2*grad_func1D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(f1D-f1D),pt)-0*grad_func1D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(2*f1D),pt)  -2*grad_func1D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(f1D*2),pt)  -2*grad_func1D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(f1D/2),pt)  -grad_func1D/2).norm(), 0.0);

        // Vector gradient
        CHECK_EQUAL((ev.eval(grad(f2D),pt)    -grad_func2D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(f2D+f2D),pt)-2*grad_func2D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(f2D-f2D),pt)-0*grad_func2D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(2*f2D),pt)  -2*grad_func2D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(f2D*2),pt)  -2*grad_func2D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(f2D/2),pt)  -grad_func2D/2).norm(), 0.0);

        // Space gradient
        CHECK_EQUAL((ev.eval(grad(u),pt)-grad_mb).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(u+u),pt)-2*grad_mb).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(u-u),pt)-0*grad_mb).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(2*u),pt)-2*grad_mb).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(u*2),pt)-2*grad_mb).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(u/2),pt)-grad_mb/2).norm(), 0.0);
    }

    TEST(jac_expr)
    {
        const index_t dim = 2;
        const index_t nAct = ev_mb.rows();
        gsMatrix<> grad_func1D = der_func1D.reshape(dim,1).transpose();
        gsMatrix<> grad_func2D = der_func2D.reshape(dim,2).transpose();
        gsMatrix<> grad_mb     = der_mb.reshape(dim,nAct).transpose();

        // Scalar gradient
        // CHECK_EQUAL((ev.eval(jac(f1D),pt)    -grad_func1D).norm(), 0.0);

        // CHECK_EQUAL((ev.eval(jac(f1D+f1D),pt)-2*grad_func1D).norm(), 0.0);
        // CHECK_EQUAL((ev.eval(jac(f1D-f1D),pt)-0*grad_func1D).norm(), 0.0);

        // CHECK_EQUAL((ev.eval(jac(2*f1D),pt)  -2*grad_func1D).norm(), 0.0);
        // CHECK_EQUAL((ev.eval(jac(f1D*2),pt)  -2*grad_func1D).norm(), 0.0);

        // CHECK_EQUAL((ev.eval(jac(f1D/2),pt)  -grad_func1D/2).norm(), 0.0);

    }

    // Test for the positive part of an expression ( Macaulay bracket )
    TEST(ppart_expr)
    {

        gsMatrix<real_t> zero1D = gsMatrix<real_t>::Zero(1,1);
        CHECK_EQUAL( ev.eval(f2D.ppart(), zero) , zero);
        CHECK_EQUAL( ev.eval(f1D.ppart(), zero) , zero1D);

        CHECK_EQUAL( ev.eval((-f2D).ppart(), zero) , zero);
        CHECK_EQUAL( ev.eval((-f1D).ppart(), zero) , zero1D);

        // Test positive part with positive values
        gsMatrix<real_t> pos_result2D = ev_func2D.cwiseMax(0.0);
        gsMatrix<real_t> pos_result1D = gsMatrix<real_t>::Constant(1,1,std::max(ev_func1D.value(),0.0));
        CHECK_EQUAL((ev.eval(f2D.ppart(), pt) - pos_result2D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(f1D.ppart(), pt) - pos_result1D).norm(), 0.0);

        // Test for the negative part of an expression
        CHECK_EQUAL( ev.eval(f2D.npart(), zero) , zero);
        CHECK_EQUAL( ev.eval(f1D.npart(), zero) , zero1D);

        CHECK_EQUAL( ev.eval((-f2D).npart(), pt) , -ev_func2D);
        CHECK_EQUAL( ev.eval((-f1D).npart(), pt) , -ev_func1D);

        // Test negative part with positive values (should be zero)
        gsMatrix<real_t> neg_result2D = (-ev_func2D).cwiseMax(0.0);
        gsMatrix<real_t> neg_result1D = gsMatrix<real_t>::Constant(1,1,std::max(-ev_func1D.value(),0.0));
        CHECK_EQUAL((ev.eval(f2D.npart(), pt) - neg_result2D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(f1D.npart(), pt) - neg_result1D).norm(), 0.0);

    }

}
