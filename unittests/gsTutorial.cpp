/** @file gsTutorial.cpp

    @brief Tutorial for setting up unittests

    This is an example unit test that doesn't really do anything useful.
    It is here as a reference for you when creating additional unit tests.
    For additional reference information, see the "test.h" header.

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

SUITE(gsTutorial)                 // The suite should have the same name as the file
{

    TEST(tutorial1)               // Declares a test named "gsTutorial:tutorial1"
    {
        CHECK(true);              // We certainly hope that true is true
        CHECK_EQUAL(2,1+1);       // The value 1+1 should equal 2

        //CHECK_EQUAL(3,1+1);     // The value 1+1 should NOT equal 3

        int x[] = {1,2,3};
        int y[] = {1,2,3};
        CHECK_ARRAY_EQUAL(x,y,3); // These arrays of length 3 are equal

        double a = 1.51;
        double b = 1.52;
        CHECK_CLOSE(a,b,0.1);     // These equal within 0.1
    }

}
