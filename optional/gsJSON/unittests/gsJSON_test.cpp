/** @file gsJSON_test.cpp

    @brief Provides unittests for the gsJSON class

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

    Author(s): H.M.Verhelst (2019 - ..., TU Delft, 2023 - ... UniFi)
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework
#include <gsJSON/gsJSON.h>

SUITE(gsJSON_test)                 // The suite should have the same name as the file
{

    TEST(JSON_input)
    {
        // Create some objects
        gsMatrix<> mat(2,2);
        mat << 1, 2, 3, 4;
        gsVector<index_t> vec(3);
        vec << 1, 2, 3;
        gsKnotVector<> kv(0,1,0,2);
        gsBSplineBasis<> basis(kv);
        gsTensorBSplineBasis<2> tbasis(kv,kv);

        gsJSON j;
        j["a"] = 2;
        j["b"] = 3.0;
        j["c"]["d"] = "four";
        j["e"] = {"e1", "e2", "e3"};
        j["A"] = mat;
        j["v"] = vec;
        j["kv"] = kv;
        j["basis"] = basis;
        j["tbasis"] = tbasis;

        CHECK(j["a"] == 2);
        CHECK(j["b"] == 3.0);
        CHECK(j["c"]["d"] == "four");
        CHECK(j["e"][0] == "e1");
        CHECK(j["e"][1] == "e2");
        CHECK(j["e"][2] == "e3");
        CHECK(j["A"].get<gsMatrix<>>() == mat);
        CHECK(j["v"].get<gsVector<index_t>>() == vec);
        CHECK(j["kv"].get<gsKnotVector<>>() == kv);
        CHECK(j["basis"].get<gsBSplineBasis<>>().knots() == basis.knots());
        CHECK(j["tbasis"].get<gsTensorBSplineBasis<2>>().knots(0) == tbasis.knots(0));
        CHECK(j["tbasis"].get<gsTensorBSplineBasis<2>>().knots(1) == tbasis.knots(1));
    }

    TEST(JSON_options)
    {
        gsOptionList opt;
        opt.addInt("a", "", 2);
        opt.addReal("b", "", 3.0);
        opt.addString("c", "", "four");
        gsJSON j(opt);
        CHECK(j["a"] == opt.getInt("a"));
        CHECK(j["b"] == opt.getReal("b"));
        CHECK(j["c"] == opt.getString("c"));
    }


}