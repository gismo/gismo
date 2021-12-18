/** @file gsMoveSemantics_test.cpp

    @brief test move sematics

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include "gismo_unittest.h"
#if (__cplusplus >= 201103L || _MSC_VER >= 1600)
#include <type_traits>
#endif

SUITE(gsMoveSemantics_test)
{

#if (__cplusplus >= 201103L || _MSC_VER >= 1600) && EIGEN_HAS_RVALUE_REFERENCES

#ifndef _MSC_VER // util::has_move_constructor doesn't work with MSVC
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

// std::is_move_constructible satisfy also copy constructors with const T&
TEST(gsMatrix_std_mc) { CHECK(std::is_move_constructible<gismo::gsMatrix<> >::value); }
TEST(gsVector_std_mc) { CHECK(std::is_move_constructible<gismo::gsVector<> >::value); }
TEST(gsSparseMatrix_std_mc) { CHECK(std::is_move_constructible<gismo::gsSparseMatrix<> >::value); }
TEST(gsSparseVector_std_mc) { CHECK(std::is_move_constructible<gismo::gsSparseVector<> >::value); }
TEST(gsKnotVector_std_mc) { CHECK(std::is_move_constructible<gismo::gsKnotVector<> >::value); }
TEST(gsBSplineBasis_std_mc) { CHECK(std::is_move_constructible<gismo::gsBSplineBasis<> >::value); }
TEST(gsBSpline_std_mc) { CHECK(std::is_move_constructible<gismo::gsBSpline<> >::value); }
TEST(gsTensorBSplineBasis_std_mc) { CHECK(std::is_move_constructible<gismo::gsTensorBSplineBasis<2> >::value); }
TEST(gsTensorBSpline_std_mc) { CHECK(std::is_move_constructible<gismo::gsTensorBSpline<2> >::value); }

// T needs to be a complete type.
//TEST(gsMatrix_nothrow_mc) { CHECK(std::is_nothrow_move_constructible<gismo::gsMatrix<> >::value); }
//TEST(gsVector_nothrow_mc) { CHECK(std::is_nothrow_move_constructible<gismo::gsVector<> >::value); }
TEST(gsKnotVector_nothrow_mc) { CHECK(std::is_nothrow_move_constructible<gismo::gsKnotVector<> >::value); }
#endif

TEST(gsMatrix_swap)
{
    gsMatrix<> a(2, 2);
    gsMatrix<> b(2, 2);

    a << 1, 2, 3, 4;
    b << 5, 6, 7, 8;

    gsMatrix<> c = a;                   // make a deep copy
    gsMatrix<> d = b;                   // make a deep copy

    size_t add_a = (size_t) &a;         // gdb p &a
    size_t add_b = (size_t) &b;         // gdb p &b
    size_t add_m_a = (size_t) &a.at(0); // gdb p a.m_storage.m_data
    size_t add_m_b = (size_t) &b.at(0); // gdb p b.m_storage.m_data

    // assert c is deep copy
    CHECK(a == c);                      // same values
    CHECK(!(b == c));
    CHECK(&a != &c);                    // copy
    CHECK(&a.at(0) != &c.at(0));        // deep copy

    // assert d is deep copy
    CHECK(!(a == d));
    CHECK(b == d);                      // same values
    CHECK(&b != &d);                    // copy
    CHECK(&b.at(0) != &d.at(0));        // deep copy

    a = give(b);                        // give swaps with flat copies under C++11

    // a
    CHECK((size_t) &a == add_a);
    CHECK((size_t) &a != add_b);
    CHECK(!(a == c));
    CHECK(a == d);                      // same values
    CHECK((size_t) &a.at(0) != add_m_a);// a isn't a any more
#if !defined(_MSC_VER) || _MSC_VER >= 1900
    CHECK((size_t) &a.at(0) == add_m_b);// flat copy, move worked
#endif

    // b
    CHECK((size_t) &b != add_a);
    CHECK((size_t) &b == add_b);
#if __cplusplus >= 201103L || _MSC_VER >= 1900
    CHECK(b == c);                      // same values
    CHECK(!(b == d));
    CHECK((size_t) &b.at(0) == add_m_a);// flat copy
    CHECK(0 != b.size());
#else
    CHECK(0 == b.size());
#endif
    CHECK((size_t) &b.at(0) != add_m_b);// b isn't b any more
}

TEST(gsMatrix_ms)
{
    gsMatrix<> b(2,2);
    b << 1, 2, 3, 4;
    CHECK_EQUAL(1, b.at(0));
    CHECK_EQUAL(3, b.at(1));
    CHECK_EQUAL(2, b.at(2));
    CHECK_EQUAL(4, b.at(3));
    size_t add_m_b = (size_t)& b.at(0);

    gsMatrix<> a(give(b));
    CHECK(0 == b.size());
#if !defined(_MSC_VER) || _MSC_VER >= 1900
    CHECK((size_t)& a.at(0) == add_m_b);	// flat copy, move worked
#endif
    CHECK_EQUAL(1, a.at(0));
    CHECK_EQUAL(3, a.at(1));
    CHECK_EQUAL(2, a.at(2));
    CHECK_EQUAL(4, a.at(3));
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
