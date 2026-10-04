/** @file gsGeometry_test.cpp

    @brief Tsting gsGeometry

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

SUITE(gsGeometry_test)                 // The suite should have the same name as the file
{
    gsGeometry<>::uPtr g = gsReadFile<>("surfaces/simple.xml");   
    gsMatrix<> uv = gsPointGrid<>(g->support(), 5);
    gsMatrix<> xyz = g->eval(uv);
    
    TEST(recoverPoints)
    {
        gsMatrix<> puv, pxy = xyz;
        g->recoverPoints(pxy, puv, 2, 1e-8);
        CHECK(  (puv-uv ).norm() < 1e-6 );
        CHECK(  (pxy-xyz).norm() < 1e-6 );
    }

    TEST(rotate2D)
    {
        gsGeometry<>::uPtr g = gsNurbsCreator<>::BSplineSquare();
        gsMatrix<> orig = g->coefs();
        g->rotate(EIGEN_PI / 2);
        // 90-degree CCW rotation around origin: (x, y) -> (-y, x)
        for (index_t i = 0; i < orig.rows(); ++i)
        {
            CHECK_CLOSE(g->coefs()(i, 0), -orig(i, 1), EPSILON);
            CHECK_CLOSE(g->coefs()(i, 1),  orig(i, 0), EPSILON);
        }
    }

    TEST(rotate3D)
    {
        gsGeometry<>::uPtr g = gsNurbsCreator<>::BSplineCube();
        gsMatrix<> orig = g->coefs();
        gsVector<real_t, 3> axis;
        axis << 1, 0, 0;
        g->rotate(EIGEN_PI / 2, axis);
        // 90-degree rotation around X-axis: (x, y, z) -> (x, -z, y)
        for (index_t i = 0; i < orig.rows(); ++i)
        {
            CHECK_CLOSE(g->coefs()(i, 0),  orig(i, 0), EPSILON);
            CHECK_CLOSE(g->coefs()(i, 1), -orig(i, 2), EPSILON);
            CHECK_CLOSE(g->coefs()(i, 2),  orig(i, 1), EPSILON);
        }
    }

}
