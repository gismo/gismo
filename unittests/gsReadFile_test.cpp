/** @file gsReadFile_test.cpp

    @brief Tests gsReadFile class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

    Author(s): J. Vogl
**/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework

using namespace gismo;

const std::string fmatrix = "unittests/matrix.xml";
const std::string f1 = "curves3d/curve_boundary.xml";    // string to existing file with gsMultiPatch
const std::string f2 = "curves3d/curve_boundary2.xml";   // string to non existing file

SUITE(gsReadFile_test)                // The suite should have the same name as the file
{
// Tests for the constructors that give back an object as parameter
TEST(gsMultiPatch_ref)
{
    gsMultiPatch<> mp;
    gsReadFile<>(f1, mp);
    CHECK(!mp.empty());
}

TEST(gsMultiPatch_ref_no_file)
{
    gsMultiPatch<real_t> mp;
    CHECK_THROW(gsReadFile<real_t>(f2, mp), std::runtime_error);
}

TEST(gsMultiPatch_ref_not_in_file)
{
    gsMultiPatch<real_t> mp;
    gsReadFile<real_t>(fmatrix, mp);
    CHECK(mp.empty());
}

TEST(Obj_ref)
{
    gsMatrix<real_t, 2, 2> exp;
    exp << 1,2,3,4;
    gsMatrix<> basis;
    gsReadFile<>(fmatrix, basis);
    for (int i = 0; i < 4; ++i) {
        if (basis.at(i) != exp.at(i))
            CHECK(false);
    }
    CHECK(true);
}

TEST(Obj_ref_no_file)
{
    gsMatrix<> basis;
    CHECK_THROW(gsReadFile<>(f2, basis), std::runtime_error);
}

TEST(Obj_ref_not_in_file)
{
    gsMatrix<> basis;
    gsReadFile<>(f1, basis);
    CHECK((basis.size() == 0));
}

// Tests for the constructors that take use of casts
/*TEST(Obj_uPtr)
{
    gsMatrix<real_t>::uPtr mp;
    if(mp = gsReadFile<gsMatrix<real_t> >(f1))
        CHECK(true);
    else
        CHECK(false);
}
TEST(Obj_uPtr_fail)
{
    gsMatrix<real_t>::uPtr mp;
    if(mp = gsReadFile<gsMatrix<real_t> >(f1))
        CHECK(false);
    else
        CHECK(true);
}*/
}
