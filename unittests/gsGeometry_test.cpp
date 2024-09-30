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

}
