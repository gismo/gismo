/** @file gsFiberMatrix_test.cpp

    @brief Unit tests for gsFiberMatrix lazy column allocation.

    Verifies that lazy mode produces identical values, pattern, and
    toSparseMatrix output as eager mode for the same insertion sequence,
    and that untouched (null) fibers report zero nonzeros.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): opencode
*/

#include "gismo_unittest.h"

#include <gsMatrix/gsFiberMatrix.h>

using namespace gismo;

SUITE(gsFiberMatrix_test)
{

    // Lazy and eager construction must agree on rows/cols/size before any
    // insertion. Untouched fibers report zero nonzeros.
    TEST(LazyShapeAndNullFibers)
    {
        const index_t rows = 5, cols = 7;
        gsFiberMatrix<real_t, ColMajor> eager(rows, cols);
        gsFiberMatrix<real_t, ColMajor> lazy;
        lazy.resizeLazy(rows, cols);

        CHECK_EQUAL(eager.rows(), lazy.rows());
        CHECK_EQUAL(eager.cols(), lazy.cols());
        CHECK_EQUAL(rows,        lazy.rows());
        CHECK_EQUAL(cols,        lazy.cols());
        CHECK_EQUAL(eager.outerSize(), lazy.outerSize());
        CHECK_EQUAL(eager.innerSize(), lazy.innerSize());

        // Untouched lazy fibers: nonZeros==0, nonZerosPerFiber all zero
        CHECK_EQUAL(0, lazy.nonZeros());
        gsVector<index_t> perFiber = lazy.nonZerosPerFiber();
        CHECK_EQUAL(lazy.outerSize(), perFiber.size());
        for (index_t i = 0; i < lazy.outerSize(); ++i)
            CHECK_EQUAL(0, perFiber[i]);

        // coeff on a null fiber returns 0 (default-constructed T)
        CHECK_EQUAL(real_t(0), lazy.coeff(0, 0));
        CHECK_EQUAL(real_t(0), lazy.coeff(rows - 1, cols - 1));

        // isExplicitZero on a null fiber returns true (no slot present)
        CHECK(lazy.isExplicitZero(0, 0));
    }

    // The same insertion sequence into lazy and eager matrices must yield
    // identical stored values and identical toSparseMatrix output.
    TEST(LazyMatchesEager)
    {
        const index_t rows = 6, cols = 6;
        gsFiberMatrix<real_t, ColMajor> eager(rows, cols);
        gsFiberMatrix<real_t, ColMajor> lazy;
        lazy.resizeLazy(rows, cols);

        // Reserve the same hint on both — lazy applies it at first touch
        eager.reservePerColumn(4);
        lazy.reservePerColumn(4);

        // A small non-trivial pattern: insert a few entries, leave some
        // columns entirely untouched.
        const index_t touchedCols[] = {0, 2, 3, 5};
        const index_t nTouched = sizeof(touchedCols) / sizeof(touchedCols[0]);
        for (index_t c = 0; c < nTouched; ++c)
        {
            const index_t j = touchedCols[c];
            for (index_t i = 0; i < rows; ++i)
            {
                const real_t v = static_cast<real_t>((i + 1) * (j + 1));
                eager.coeffRef(i, j) += v;
                lazy.coeffRef (i, j) += v;
            }
        }

        CHECK_EQUAL(eager.nonZeros(), lazy.nonZeros());

        // Compare stored values entry by entry through coeff
        for (index_t i = 0; i < rows; ++i)
            for (index_t j = 0; j < cols; ++j)
                CHECK_CLOSE(eager.coeff(i, j), lazy.coeff(i, j), 1e-12);

        // Compare the full toSparseMatrix output
        gsSparseMatrix<real_t> Me, Ml;
        eager.toSparseMatrix_into(Me);
        lazy.toSparseMatrix_into (Ml);
        CHECK_EQUAL(Me.rows(), Ml.rows());
        CHECK_EQUAL(Me.cols(), Ml.cols());
        CHECK_EQUAL(Me.nonZeros(), Ml.nonZeros());
        // Eigen SparseMatrix - SparseMatrix keeps structural zeros from
        // cancellation, so diff.nonZeros() is misleading. Value-compare via
        // the dense norm instead.
        CHECK((Me - Ml).toDense().norm() < 1e-12);

        // Untouched columns must report zero nonzeros in lazy
        gsVector<index_t> lazyPerFiber = lazy.nonZerosPerFiber();
        for (index_t j = 0; j < cols; ++j)
        {
            const bool touched = (j == 0 || j == 2 || j == 3 || j == 5);
            CHECK_EQUAL(touched ? rows : 0, lazyPerFiber[j]);
        }
    }

    // insertExplicitZero must allocate a lazy fiber (explicit zeros need a
    // slot). isExplicitZero follows the gismo/gsSparseMatrix contract: it
    // returns true when the (i,j) entry is NOT in the sparsity (i.e. is a
    // structural / "explicit" zero) and false when a slot exists.
    TEST(LazyInsertExplicitZero)
    {
        const index_t rows = 4, cols = 4;
        gsFiberMatrix<real_t, ColMajor> lazy;
        lazy.resizeLazy(rows, cols);

        // Null fiber → no slot → isExplicitZero true
        CHECK(lazy.isExplicitZero(1, 1));
        CHECK(lazy.isExplicitZero(1, 2));

        lazy.insertExplicitZero(1, 2);
        // Slot now exists at (1,2) → isExplicitZero false (it has a real
        // slot, NOT a structural/absent zero)
        CHECK(!lazy.isExplicitZero(1, 2));
        // (1,1) still has no slot
        CHECK(lazy.isExplicitZero(1, 1));
        // nonZeros counts the explicit-zero slot that was inserted
        CHECK_EQUAL(1, lazy.nonZeros());
    }

    // RowMajor lazy matrix: same property — untouched rows null, matches
    // eager on the same insertion sequence.
    TEST(LazyRowMajor)
    {
        const index_t rows = 5, cols = 4;
        gsFiberMatrix<real_t, RowMajor> eager(rows, cols);
        gsFiberMatrix<real_t, RowMajor> lazy;
        lazy.resizeLazy(rows, cols);

        const index_t touchedRows[] = {0, 2, 4};
        for (index_t r = 0; r < 3; ++r)
        {
            const index_t i = touchedRows[r];
            for (index_t j = 0; j < cols; ++j)
            {
                const real_t v = static_cast<real_t>((i + 1) * (j + 1));
                eager.row(i).coeffRef(j) += v;
                lazy.row (i).coeffRef(j) += v;
            }
        }

        CHECK_EQUAL(eager.nonZeros(), lazy.nonZeros());
        gsSparseMatrix<real_t> Me, Ml;
        eager.toSparseMatrix_into(Me);
        lazy.toSparseMatrix_into (Ml);
        CHECK_EQUAL(Me.rows(), Ml.rows());
        CHECK_EQUAL(Me.cols(), Ml.cols());
        CHECK_EQUAL(Me.nonZeros(), Ml.nonZeros());
        CHECK((Me - Ml).toDense().norm() < 1e-12);
    }

    // setZero on a lazy matrix with some null and some allocated fibers.
    // Eigen SparseVector::setZero() clears all slots (nonZeros → 0), so
    // values become 0 AND the sparsity pattern is gone. assignZero() is the
    // keep-pattern-zero-values alternative (used by clearMatrix's
    // save_sparsety_pattern branch in gsExprAssembler).
    TEST(LazySetZero)
    {
        const index_t rows = 4, cols = 4;
        gsFiberMatrix<real_t, ColMajor> lazy;
        lazy.resizeLazy(rows, cols);
        lazy.coeffRef(0, 0) = 1.0;
        lazy.coeffRef(1, 1) = 2.0;
        CHECK_EQUAL(2, lazy.nonZeros());

        // setZero() wipes the pattern
        lazy.setZero();
        CHECK_EQUAL(0, lazy.nonZeros());
        CHECK_EQUAL(real_t(0), lazy.coeff(0, 0));
        CHECK_EQUAL(real_t(0), lazy.coeff(1, 1));
        CHECK_EQUAL(real_t(0), lazy.coeff(2, 2));
    }

    // assignZero() zeroes values but keeps the sparsity pattern.
    TEST(LazyAssignZero)
    {
        const index_t rows = 4, cols = 4;
        gsFiberMatrix<real_t, ColMajor> lazy;
        lazy.resizeLazy(rows, cols);
        lazy.coeffRef(0, 0) = 1.0;
        lazy.coeffRef(1, 1) = 2.0;
        CHECK_EQUAL(2, lazy.nonZeros());

        lazy.assignZero();
        CHECK_EQUAL(2, lazy.nonZeros()); // slots preserved
        CHECK_EQUAL(real_t(0), lazy.coeff(0, 0));
        CHECK_EQUAL(real_t(0), lazy.coeff(1, 1));
        // Untouched columns still null — coeff returns 0
        CHECK_EQUAL(real_t(0), lazy.coeff(2, 2));
    }

    // Copy ctor of a lazy matrix must preserve null-vs-allocated per fiber
    // and the values of allocated fibers.
    TEST(LazyCopy)
    {
        const index_t rows = 4, cols = 4;
        gsFiberMatrix<real_t, ColMajor> lazy;
        lazy.resizeLazy(rows, cols);
        lazy.coeffRef(0, 0) = 7.0;
        lazy.coeffRef(2, 1) = 3.0;

        gsFiberMatrix<real_t, ColMajor> copy(lazy);
        CHECK_EQUAL(rows, copy.rows());
        CHECK_EQUAL(cols, copy.cols());
        CHECK_EQUAL(lazy.nonZeros(), copy.nonZeros());
        CHECK_EQUAL(7.0, copy.coeff(0, 0));
        CHECK_EQUAL(3.0, copy.coeff(2, 1));
        CHECK_EQUAL(real_t(0), copy.coeff(1, 1));
        // Per-fiber nnz matches
        gsVector<index_t> a = lazy.nonZerosPerFiber();
        gsVector<index_t> b = copy.nonZerosPerFiber();
        for (index_t i = 0; i < a.size(); ++i)
            CHECK_EQUAL(a[i], b[i]);
    }

    // fiberPointerBytes()/fiberDataBytes() must report a real eager-vs-lazy
    // split: identical pointer-array footprint (same outer size) but a
    // drastically smaller data footprint for lazy mode when only a handful
    // of fibers have been touched.
    TEST(ByteAccountingSplit)
    {
        const index_t rows = 4096, cols = 4096;

        gsFiberMatrix<real_t, ColMajor> eager;
        eager.resize(rows, cols);
        eager.reservePerColumn(8);

        gsFiberMatrix<real_t, ColMajor> lazy;
        lazy.resizeLazy(rows, cols);
        lazy.reservePerColumn(8);
        lazy.coeffRef(0, 0) = 1.0;
        lazy.coeffRef(1, 1) = 2.0;

        CHECK_EQUAL(eager.fiberPointerBytes(), lazy.fiberPointerBytes());
        CHECK(lazy.fiberDataBytes() * 100 < eager.fiberDataBytes());
        CHECK(lazy.fiberDataBytes() > 0);
    }

} // SUITE(gsFiberMatrix)