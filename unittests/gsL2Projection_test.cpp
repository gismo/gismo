/** @file gsL2Projection_test.cpp

    @brief Tests the L2 projection

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
 **/

//#define TEST_INFO

#include "gismo_unittest.h"

SUITE(gsL2Projection_test)
{

    TEST(test) // Declares test
    {
        // Make a geometry
        gsMultiPatch<> mp;
        mp.addPatch(gsNurbsCreator<>::NurbsAnnulus());

        // Make a basis
        gsMultiBasis<> mb(mp);
        mb.degreeElevate(1,1); // elevate in x
        mb.degreeElevate(1); // elevate in x
        mb.degreeElevate(1); // elevate in x
        mb.uniformRefine();

        // Coefficients object
        gsMatrix<> coefs_SP, coefs_MP;

        real_t error;

        // Project the geometry on the finer basis
        // - Single patch version
        //   * same basis for integration
        error = gsL2Projection<real_t>::project(mb.basis(0),mp.patch(0),coefs_SP);
        gsTestInfo<<"Geometry projection 1, error = "<<error<<"\n";
        CHECK_CLOSE(error,0.0,1e-10);
        //   * basis for integration specified
        error = gsL2Projection<real_t>::project(mb.basis(0),mb.basis(0),mp.patch(0),coefs_SP);
        gsTestInfo<<"Geometry projection 2, error = "<<error<<"\n";
        CHECK_CLOSE(error,0.0,1e-10);

        // - Single patch version
        //   * same basis for integration
        error = gsL2Projection<real_t>::project(mb,mp,coefs_MP);
        gsTestInfo<<"Geometry projection 3, error = "<<error<<"\n";
        CHECK_CLOSE(error,0.0,1e-10);
        //   * basis for integration specified
        error = gsL2Projection<real_t>::project(mb,mb,mp,coefs_MP);
        gsTestInfo<<"Geometry projection 4, error = "<<error<<"\n";
        CHECK_CLOSE(error,0.0,1e-10);

        std::vector<std::string> expressions = {"x^2*y","y^2*x","x^3","y^4"};
        gsFunctionExpr<> f(expressions,2);

        // Project a function on a basis
        // - Single patch version
        //   * same basis for integration
        error = gsL2Projection<real_t>::project(mb.basis(0),mp.patch(0),f,coefs_SP);
        gsTestInfo<<"Function projection 1, error = "<<error<<"\n";
        CHECK_CLOSE(error,0.0,1e-10);
        //   * basis for integration specified
        error = gsL2Projection<real_t>::project(mb.basis(0),mb.basis(0),mp.patch(0),f,coefs_SP);
        gsTestInfo<<"Function projection 2, error = "<<error<<"\n";
        CHECK_CLOSE(error,0.0,1e-10);

        // - Single patch version
        //   * same basis for integration
        error = gsL2Projection<real_t>::project(mb,mp,f,coefs_MP);
        gsTestInfo<<"Function projection 3, error = "<<error<<"\n";
        CHECK_CLOSE(error,0.0,1e-10);
        //   * basis for integration specified
        error = gsL2Projection<real_t>::project(mb,mb,mp,f,coefs_MP);
        gsTestInfo<<"Function projection 4, error = "<<error<<"\n";
        CHECK_CLOSE(error,0.0,1e-10);
    }

}
