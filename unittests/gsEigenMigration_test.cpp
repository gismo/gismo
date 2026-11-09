/** @file gsEigenMigration_test.cpp

    @brief Pins Eigen-5 migration invariants: the reported Eigen version, namespace
    isolation of the gsEigen rename, move ASSIGNMENT into a non-empty target for
    gsSparseVector, gsTensorBSpline<2> and gsTensorBSplineBasis<2>, the
    gsIncompleteLUT factor/permutation accessors, and that static asserts stay
    active.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gismo_unittest.h"

#include <type_traits>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <limits>

// Eigen 5 compiles EIGEN_STATIC_ASSERT to nothing when EIGEN_NO_STATIC_ASSERT
// is defined (E5/Eigen/src/Core/util/StaticAssert.h); it must stay undefined.
#ifdef EIGEN_NO_STATIC_ASSERT
#error "EIGEN_NO_STATIC_ASSERT is defined: Eigen 5 then compiles EIGEN_STATIC_ASSERT to nothing"
#endif

namespace {

// n x n path-graph (1D SPD) Laplacian: 2 on the diagonal, -1 on the two
// off-diagonals. A minimum-degree ordering of this pattern produces no fill.
gsSparseMatrix<real_t> gsEigenMigration_laplacian(index_t n)
{
    gsSparseMatrix<real_t> A(n, n);
    for (index_t i = 0; i != n; ++i)
    {
        A(i, i) = 2;
        if (i > 0)     A(i, i - 1) = -1;
        if (i < n - 1) A(i, i + 1) = -1;
    }
    A.makeCompressed();
    return A;
}

// Cyclic shift by \a shift, as a permutation matrix: index i maps to (i+shift)%n.
gsEigen::PermutationMatrix<gsEigen::Dynamic, gsEigen::Dynamic, index_t>
gsEigenMigration_cyclicShift(index_t n, index_t shift)
{
    gsEigen::PermutationMatrix<gsEigen::Dynamic, gsEigen::Dynamic, index_t> Q(n);
    for (index_t i = 0; i != n; ++i)
        Q.indices()(i) = (i + shift) % n;
    return Q;
}

}

SUITE(gsEigenMigration_test)
{

// --- Version -----------------------------------------------------------
// gsSysInfo::getEigenVersion() must build "<major>.<minor>.<patch>" from
// EIGEN_MAJOR_VERSION/MINOR/PATCH. EIGEN_WORLD_VERSION is pinned at 3, so a
// string built from WORLD.MAJOR.MINOR would misreport Eigen 5.0.0 as "3.5.0".
// Parsed independently of the macros that build the string, so this does not
// just restate the implementation.

TEST(eigen_version_reported)
{
    CHECK(EIGEN_MAJOR_VERSION >= 5);

    const std::string v = gsSysInfo::getEigenVersion();

    CHECK_EQUAL(std::size_t(2), std::size_t(std::count(v.begin(), v.end(), '.')));

    std::vector<std::string> fields;
    std::string field;
    for (std::string::const_iterator it = v.begin(); it != v.end(); ++it)
    {
        if (*it == '.') { fields.push_back(field); field.clear(); }
        else field.push_back(*it);
    }
    fields.push_back(field);

    CHECK_EQUAL(std::size_t(3), fields.size());
    for (std::size_t i = 0; i != fields.size(); ++i)
    {
        CHECK(!fields[i].empty());
        for (std::string::const_iterator c = fields[i].begin(); c != fields[i].end(); ++c)
            CHECK(0 != std::isdigit(static_cast<unsigned char>(*c)));
    }

    const int major = std::atoi(fields[0].c_str());
    CHECK(major >= 5);
    CHECK_EQUAL(EIGEN_MAJOR_VERSION, major);
}

// --- Namespace isolation -------------------------------------------------
// src/gsCore/gsLinearAlgebra.h brackets '#define Eigen gsEigen' with a
// matching '#undef Eigen' before the header closes, so the plain gsEigen
// namespace -- not a macro -- carries the renamed Eigen sources.

TEST(gsMatrix_derives_from_gsEigen)
{
    CHECK((std::is_base_of<gsEigen::Matrix<real_t, gsEigen::Dynamic, gsEigen::Dynamic>,
                            gsMatrix<real_t> >::value));
    CHECK((std::is_base_of<gsEigen::SparseMatrix<real_t, 0, index_t>,
                            gsSparseMatrix<real_t> >::value));
}

TEST(Eigen_macro_does_not_leak)
{
#ifdef Eigen
    CHECK(false);   // the renaming macro leaked out of gismo.h
#else
    CHECK(true);
#endif
}

// --- Move assignment into a non-empty target ------------------------------

TEST(gsSparseVector_move_assign)
{
    // gsSparseVector's move assignment is swap+clear (src/gsMatrix/gsSparseVector.h),
    // unlike gsEigen::SparseVector's own move members, which are a bare swap
    // and leave the source holding the target's old contents. Only a move
    // into a NON-EMPTY target tells the two apart.
    gsSparseVector<> a(5); a(3) = 7;
    gsSparseVector<> b(2); b(0) = 1;
    const real_t * pb = b.valuePtr();

    a = give(b);

    CHECK_EQUAL(2, a.size());
    CHECK_EQUAL(real_t(1), a.coeff(0));
    CHECK(a.valuePtr() == pb);

    CHECK_EQUAL(0, b.size());
    CHECK_EQUAL(0, b.nonZeros());

    // Copy control: a copy-and-swap or a plain copy neither transfers the
    // buffer nor empties the source.
    gsSparseVector<> a2(5); a2(3) = 7;
    gsSparseVector<> b2(2); b2(0) = 1;
    a2 = b2;
    CHECK(a2.valuePtr() != b2.valuePtr());
    CHECK_EQUAL(1, b2.nonZeros());
}

TEST(gsTensorBSpline_move_assign)
{
    // gsTensorBSpline<2> declares no special members: its implicit move
    // assignment resolves to gsGeometry<T>::operator=(gsGeometry&&)
    // (src/gsCore/gsGeometry.h), which swaps the coefficient matrix and
    // transfers basis ownership, nulling the source's basis pointer. b.basis()
    // must never be called after the move: it would dereference that pointer.
    gsKnotVector<> kv(0, 1, 2, 3);
    gsTensorBSplineBasis<2> basis(kv, kv);
    gsMatrix<> cc = gsMatrix<>::Zero(basis.size(), 1);
    gsTensorBSpline<2> b(basis, cc);

    gsKnotVector<> kv2(0, 1, 1, 2);
    gsTensorBSplineBasis<2> basis2(kv2, kv2);
    gsMatrix<> cc2 = gsMatrix<>::Zero(basis2.size(), 2);
    gsTensorBSpline<2> a(basis2, cc2);

    const gsBasis<real_t> * pb = &b.basis();
    const real_t * pc = b.coefs().data();

    a = give(b);

    CHECK(&a.basis() == pb);
    CHECK(a.coefs().data() == pc);
    CHECK_EQUAL(1, a.coefs().cols());

    CHECK_EQUAL(0, b.coefs().rows());

    // Copy control: a copy-and-swap or a plain copy allocates a fresh basis
    // and coefficient buffer, leaving the source intact.
    gsTensorBSpline<2> b2(basis, cc);
    gsTensorBSpline<2> a2(basis2, cc2);
    a2 = b2;
    CHECK(&a2.basis() != &b2.basis());
    CHECK(b2.coefs().rows() != 0);
}

TEST(gsTensorBSplineBasis_move_assign)
{
    // gsTensorBSplineBasis<2> declares no special members either: its
    // implicit move assignment resolves to gsTensorBasis<d>::operator=
    // (src/gsTensor/gsTensorBasis.h), which frees the target's own component
    // bases, transfers the source's pointers, and nulls them in the source --
    // isValid() then reports the moved-from source as invalid. component(dir)
    // must never be called after the move.
    gsKnotVector<> kv(0, 1, 2, 3), kv2(0, 1, 1, 2);
    gsTensorBSplineBasis<2> b(kv, kv);
    gsTensorBSplineBasis<2> a(kv2, kv2);

    const gsBSplineBasis<real_t> * p0 = &b.component(0);
    const gsBSplineBasis<real_t> * p1 = &b.component(1);

    a = give(b);

    CHECK(&a.component(0) == p0);
    CHECK(&a.component(1) == p1);
    CHECK(!b.isValid());

    // Copy control: a copy-and-swap or a plain copy allocates its own
    // component bases, leaving the source valid and untouched.
    gsTensorBSplineBasis<2> b2(kv, kv);
    gsTensorBSplineBasis<2> a2(kv2, kv2);
    a2 = b2;
    CHECK(&a2.component(0) != &b2.component(0));
    CHECK(b2.isValid());
}

// --- gsIncompleteLUT -------------------------------------------------------
// A = Q . A0 . Q^T (A0 the SPD Laplacian, Q a cyclic shift) so that the
// factorization's AMD ordering is not the identity permutation.

TEST(gsIncompleteLUT_accessors)
{
    const index_t n = 12;
    const index_t shift = 5;

    gsSparseMatrix<real_t> A0 = gsEigenMigration_laplacian(n);
    gsEigen::PermutationMatrix<gsEigen::Dynamic, gsEigen::Dynamic, index_t> Q
        = gsEigenMigration_cyclicShift(n, shift);
    gsSparseMatrix<real_t> A = A0.twistedBy(Q);

    gsIncompleteLUT<real_t> ilu(A);

    CHECK_EQUAL(n, ilu.factors().rows());
    CHECK_EQUAL(n, ilu.factors().cols());
    CHECK_EQUAL(n, ilu.fillReducingPermutation().size());
    CHECK_EQUAL(n, ilu.inversePermutation().size());

    const gsIncompleteLUT<real_t>::PermutationType & P    = ilu.fillReducingPermutation();
    const gsIncompleteLUT<real_t>::PermutationType & Pinv = ilu.inversePermutation();

    for (index_t i = 0; i != n; ++i)
        CHECK_EQUAL(i, P.indices()(Pinv.indices()(i)));

    // Precondition, kept as a permanent check: without it, a bug that swaps
    // the P/Pinv accessors would go undetected by the comparisons below.
    CHECK(P.indices() != Pinv.indices());

    gsVector<real_t> b(n);
    for (index_t i = 0; i != n; ++i)
        b(i) = real_t(i + 1);

    gsMatrix<real_t> LU = ilu.factors().toDense();
    gsVector<real_t> y = Pinv * b;
    y = LU.triangularView<gsEigen::UnitLower>().solve(y);
    y = LU.triangularView<gsEigen::Upper>().solve(y);
    y = P * y;

    gsVector<real_t> x = ilu.solve(b);

    const real_t tol = 100 * std::numeric_limits<real_t>::epsilon() * x.lpNorm<gsEigen::Infinity>();
    CHECK((y - x).lpNorm<gsEigen::Infinity>() <= tol);

    // Same accessors, same solve() code path as gsEigen::IncompleteLUT itself.
    gsEigen::IncompleteLUT<real_t, index_t> eigenIlu(A);
    gsVector<real_t> xe = eigenIlu.solve(b);
    CHECK((x - xe).lpNorm<gsEigen::Infinity>() <= tol);
}

TEST(gsIncompleteLUT_exact_reconstruction)
{
    const index_t n = 12;
    const index_t shift = 5;

    gsSparseMatrix<real_t> A0 = gsEigenMigration_laplacian(n);
    gsEigen::PermutationMatrix<gsEigen::Dynamic, gsEigen::Dynamic, index_t> Q
        = gsEigenMigration_cyclicShift(n, shift);
    gsSparseMatrix<real_t> A = A0.twistedBy(Q);

    // droptol = 0, fillfactor = n: the path-graph pattern produces no fill under a
    // minimum-degree ordering and an SPD matrix has no zero pivot, so LU is exact.
    gsIncompleteLUT<real_t> ilu;
    ilu.setDroptol(0);
    ilu.setFillfactor(n);
    ilu.compute(A);
    CHECK(ilu.info() == gsEigen::Success);

    gsMatrix<real_t> LU = ilu.factors().toDense();
    gsMatrix<real_t> L = LU.triangularView<gsEigen::StrictlyLower>();
    L.diagonal().setOnes();
    gsMatrix<real_t> U = LU.triangularView<gsEigen::Upper>();

    const gsIncompleteLUT<real_t>::PermutationType & P    = ilu.fillReducingPermutation();
    const gsIncompleteLUT<real_t>::PermutationType & Pinv = ilu.inversePermutation();

    gsMatrix<real_t> R = P * (L * U) * Pinv;

    const real_t tol = 100 * std::numeric_limits<real_t>::epsilon() * A.toDense().norm();
    CHECK((R - A.toDense()).norm() <= tol);
}

}
