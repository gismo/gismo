/** @file gsKnotVectors.cpp

    @brief test move sematics

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/
 
#include "gismo_unittest.h"  // Brings in G+Smo and the UnitTest++ framework

SUITE(gsMoveSemantics)
{

#if __cplusplus >=201103 && EIGEN_HAS_RVALUE_REFERENCES
TEST(gsMatrix_mc) { CHECK(util::has_move_constructor<gismo::gsMatrix<> >::value); }
TEST(gsVector_mc) { CHECK(util::has_move_constructor<gismo::gsVector<> >::value); }
TEST(gsSparseMatrix_mc) { CHECK(util::has_move_constructor<gismo::gsSparseMatrix<> >::value); }
TEST(gsSparseVector_mc) { CHECK(util::has_move_constructor<gismo::gsSparseVector<> >::value); }
TEST(gsKnotVector_mc) { CHECK(util::has_move_constructor<gismo::gsKnotVector<> >::value); }
TEST(gsBSplineBasis_mc) { CHECK(util::has_move_constructor<gismo::gsBSplineBasis<> >::value); }
TEST(gsBSpline_mc) { CHECK(util::has_move_constructor<gismo::gsBSpline<> >::value); }
TEST(gsTensorBSplineBasis_mc) { CHECK(util::has_move_constructor<gismo::gsTensorBSplineBasis<2> >::value); }
TEST(gsTensorBSpline_mc) { CHECK(util::has_move_constructor<gismo::gsTensorBSpline<2> >::value); }
#endif

TEST(gsMatrix_ms)
{
    gsMatrix<> b(2,2);
    gsMatrix<> a = give(b);
    CHECK(0 == b.size());
    /*
      gsMatrix<> c(2,2);
      c = give(b);
      // fail: Eigen::Matrix move assignment is implemented as swap
      CHECK(0 == b.size());
    */
}

TEST(gsVector_ms)
{
    gsVector<> b(4);
    gsVector<> a = give(b);
    CHECK(0 == b.size());
}

TEST(gsSparseMatrix_ms)
{
    gsSparseMatrix<> b(2,2); b(0,0)=1;
    gsSparseMatrix<> a = give(b);
    CHECK(0 == b.nonZeros());
}

TEST(gsSparseVector_ms)
{
    gsSparseVector<> b(2); b(0)=1;
    gsSparseVector<> a = give(b);
    CHECK(0 == b.nonZeros());
}

TEST(gsKnotVector_ms)
{
    gsKnotVector<> b(0, 1, 2, 3);
    gsKnotVector<> a = give(b);
    CHECK(0 == b.size());
}

TEST(gsBSplineBasis_ms)
{
    gsBSplineBasis<> b(0., 1., 2, 3);
    gsBSplineBasis<> a = give(b);
    CHECK(b.knots().size() <= 2);
}

TEST(gsBSpline_ms)
{
    gsBSplineBasis<> x(0., 1., 2, 3);
    gsMatrix<> cc = gsMatrix<>::Zero(x.size(),1);
    gsBSpline<> b(x,cc);
    gsBSpline<> a = give(b);
    CHECK(b.coefs().rows()==0);
}

TEST(gsTensorBSplineBasis_ms)
{
    //gsTensorBSplineBasis<2> c(1); // error
    gsKnotVector<> x(0, 1, 2, 3), y(0, 1, 2, 3);
    gsTensorBSplineBasis<2> b(give(x), give(y));
    gsTensorBSplineBasis<2> a = give(b);
    CHECK( !b.isValid() || b.knots(0).size()<=2 );
}

TEST(gsTensorBSpline_ms)
{
    gsKnotVector<> x(0, 1, 2, 3);
    gsTensorBSplineBasis<2> basis(x, x);
    gsMatrix<> cc = gsMatrix<>::Zero(basis.size(),1);
    gsTensorBSpline<2> b(give(basis), give(cc));
    gsTensorBSpline<2> a = give(b);
    CHECK( b.coefs().rows()==0 );
}

}
