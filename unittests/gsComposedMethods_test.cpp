/** @file gsComposedMethods_test.cpp

    @brief Testing gsComposedX classes

    == BASIC REFERENCE ==
         - TEST(NAME_OF_TEST) { body_of_test }
         - TEST_FIXTURE(NAME_OF_FIXTURE,NAME_OF_TEST){ body_of_test }

    == CHECK MACRO REFERENCE ==
         - CHECK(EXPR);
         - CHECK_EQUAL(EXPECTED,ACTUAL);
         - CHECK_CLOSE(EXPECTED,ACTUAL,EPSILON);
         - CHECK_ARRAY_EQUAL(EXPECTED,ACTUAL,LENGTH);
         - CHECK_ARRAY_CLOSE(EXPECTED,ACTUAL,LENGTH,EPSILON);
         - CHECK_ARRAY2D_EQUAL(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT);
         - CHECK_ARRAY2D_CLOSE(EXPECTED,ACTUAL,ROWCOUNT,COLCOUNT,EPSILON);
         - CHECK_THROW(EXPR,EXCEPTION_TYPE_EXPECTED);

    == TIME CONSTRAINTS ==
         - UNITTEST_TIME_CONSTRAINT(TIME_IN_MILLISECONDS);
         - UNITTEST_TIME_CONSTRAINT_EXEMPT();

    == MORE INFO ==
         See: https://unittest-cpp.github.io/

    Author(s): A. Mantzaflaris,  H. Weiner
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework

SUITE(gsComposedMethods_test)                 // The suite should have the same name as the file
{
   TEST(gsComposedFunction_test)
    {
        gsGeometry<>::uPtr com = gsNurbsCreator<>::BSplineSquareDeg(2);
        com->coefs().row(4).setConstant(0.8);

        gsVector<> low = com->support().col(0);
        gsVector<> upp = com->support().col(1);
        gsMatrix<> par = uniformPointGrid(low,upp,10);
        gsMatrix<> cpar = com->eval(par);

        gsFunctionExpr<> fun("sin(x)*tanh(y)",2);

        gsComposedFunction<> cfun(*com,fun);
        gsMatrix<> eval = fun.eval(cpar);
        gsMatrix<> ceval= cfun.eval(par);
        CHECK_CLOSE(0.0, (eval-ceval).norm(), 1e-10);
    }

    TEST(gsComposedGeometry_test)
    {
        gsGeometry<>::uPtr com = gsNurbsCreator<>::BSplineSquareDeg(2);
        com->coefs().row(4).setConstant(0.8);

        gsVector<> low = com->support().col(0);
        gsVector<> upp = com->support().col(1);
        gsMatrix<> par = uniformPointGrid(low,upp,10);
        gsMatrix<> cpar = com->eval(par);
        for (index_t k=0; k!=cpar.cols(); k++)
        {
            cpar.col(k) = cpar.col(k).cwiseMax(low);
            cpar.col(k) = cpar.col(k).cwiseMin(upp);
        }

        gsGeometry<>::uPtr geo = gsNurbsCreator<>::NurbsAnnulus();
        geo->degreeElevate();

        gsComposedGeometry<> cgeom(*com,*geo);
        gsMatrix<> eval = geo->eval(cpar);
        gsMatrix<> ceval= cgeom.eval(par);
        CHECK_CLOSE(0.0, (eval-ceval).norm(), 1e-10);
        /* \todo: test derivatives */
    }

    TEST(gsComposedBasis_test)
    {
        gsGeometry<>::uPtr com = gsNurbsCreator<>::BSplineSquareDeg(2);
        com->coefs().row(4).setConstant(0.8);

        gsVector<> low = com->support().col(0);
        gsVector<> upp = com->support().col(1);
        gsMatrix<> par = uniformPointGrid(low,upp,10);
        gsMatrix<> cpar = com->eval(par);
        for (index_t k=0; k!=cpar.cols(); k++)
        {
            cpar.col(k) = cpar.col(k).cwiseMax(low);
            cpar.col(k) = cpar.col(k).cwiseMin(upp);
        }

        gsKnotVector<> KV(0,1,10,3,1,2); // start, end, #interior, mult. end, mult. int, degree
        gsTensorBSplineBasis<2> tbasis(KV,KV);
        gsTHBSplineBasis<2> thbbasis(tbasis);
        gsMatrix<> box(2,2);
        box.col(0).setConstant(0.0);
        box.col(1).setConstant(0.5);
        thbbasis.refine(box);

        gsComposedBasis<> cbasis(*com,thbbasis);
        gsMatrix<> eval = thbbasis.eval(cpar);
        gsMatrix<> ceval= cbasis.eval(par);
        CHECK_CLOSE(0.0, (eval-ceval).norm(), 1e-10);
        gsMatrix<index_t> act = thbbasis.active(cpar);
        gsMatrix<index_t> cact= cbasis.active(par);
        CHECK_CLOSE(0.0, (act-cact).norm(), 1e-10);
    }
}
