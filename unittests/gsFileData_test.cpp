/** @file gsFileData_test.cpp

    @brief Tests for gsFileData

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Claude Fable 5.1
**/

#include "gismo_unittest.h"
#include <cstdio>

SUITE(gsFileData_test)
{

// -------------------------------------------------------------------------
// getLabel
// -------------------------------------------------------------------------

TEST(getLabel_existing_and_missing)
{
    // Write a single labelled matrix to a temporary xml file, then
    // read it back through gsFileData::getLabel.
    gsMatrix<real_t> M(2, 2);
    M << 1, 2, 3, 4;

    gsFileData<> fdOut;
    fdOut.addWithLabel(M, "myMatrix");
    const std::string fn = "gsFileData_test_getLabel.xml";
    fdOut.save(fn);

    gsFileData<> fdIn(fn);

    gsMatrix<real_t> result;
    fdIn.getLabel("myMatrix", result);
    CHECK_EQUAL(M.rows(), result.rows());
    CHECK_EQUAL(M.cols(), result.cols());
    for (index_t i = 0; i < M.rows(); ++i)
        for (index_t j = 0; j < M.cols(); ++j)
            CHECK_EQUAL(M(i,j), result(i,j));

    gsMatrix<real_t> missing;
    CHECK_THROW(fdIn.getLabel("noSuchLabel", missing), std::runtime_error);

    std::remove(fn.c_str());
}

}
