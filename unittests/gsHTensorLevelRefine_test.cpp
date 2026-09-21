/** @file gsHTensorLevelRefine_test.cpp

    @brief Tests level-based refine/unrefine for hierarchical tensor bases.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
**/

#include "gismo_unittest.h"

using namespace gismo;

namespace {

gsTHBSplineBasis<2> makeBasis()
{
    gsKnotVector<> kv(0, 1, 3, 3); // degree 2, 3 interior knots
    gsTensorBSplineBasis<2> tens(kv, kv);
    return gsTHBSplineBasis<2>(tens);
}

index_t minLeafLevel(const gsHTensorBasis<2> & basis)
{
    index_t minLvl = basis.maxLevel();
    for (auto leafIt = basis.tree().beginLeafIterator(); leafIt.good(); leafIt.next())
        minLvl = math::min(minLvl, static_cast<index_t>(leafIt.level()));
    return minLvl;
}

index_t maxLeafLevel(const gsHTensorBasis<2> & basis)
{
    index_t maxLvl = 0;
    for (auto leafIt = basis.tree().beginLeafIterator(); leafIt.good(); leafIt.next())
        maxLvl = math::max(maxLvl, static_cast<index_t>(leafIt.level()));
    return maxLvl;
}

} // namespace

SUITE(gsHTensorLevelRefine_test)
{
    TEST(refineToLevel_raises_coarse_leaves)
    {
        gsTHBSplineBasis<2> basis = makeBasis();
        CHECK_EQUAL(0, maxLeafLevel(basis));

        basis.refineToLevel(2);
        CHECK(minLeafLevel(basis) >= 2);
        CHECK(maxLeafLevel(basis) >= 2);
        CHECK(basis.size() > 0);
    }

    TEST(unrefineToLevel_lowers_fine_leaves)
    {
        gsTHBSplineBasis<2> basis = makeBasis();
        basis.refineToLevel(2);
        CHECK(maxLeafLevel(basis) >= 2);

        basis.unrefineToLevel(1);
        CHECK(maxLeafLevel(basis) <= 1);
        CHECK(basis.size() > 0);
    }

    TEST(refineCoarsest_and_unrefineFinest)
    {
        gsTHBSplineBasis<2> basis = makeBasis();
        const index_t size0 = basis.size();

        basis.refineCoarsestLevel();
        CHECK(maxLeafLevel(basis) >= 1);
        CHECK(basis.size() > size0);

        basis.unrefineFinestLevel();
        CHECK(maxLeafLevel(basis) == 0);
    }

    TEST(refineToLevel_withTransfer_dimensions)
    {
        gsTHBSplineBasis<2> basis = makeBasis();
        const index_t nCoarse = basis.size();

        gsSparseMatrix<real_t, RowMajor> transfer;
        basis.refineToLevel_withTransfer(1, transfer);

        CHECK_EQUAL(nCoarse, transfer.cols());
        CHECK_EQUAL(basis.size(), transfer.rows());
        CHECK(transfer.nonZeros() > 0);
    }

    TEST(unrefineToLevel_withTransfer_dimensions)
    {
        gsTHBSplineBasis<2> basis = makeBasis();
        basis.refineToLevel(2);
        const index_t nFine = basis.size();

        gsSparseMatrix<real_t, RowMajor> transfer;
        basis.unrefineToLevel_withTransfer(1, transfer);

        CHECK_EQUAL(basis.size(), transfer.cols());
        CHECK_EQUAL(nFine, transfer.rows());
        CHECK(transfer.nonZeros() > 0);
    }

    TEST(refineToLevel_withCoefs_preserves_geometry)
    {
        gsTHBSplineBasis<2> basis = makeBasis();
        gsMatrix<> coefs = gsMatrix<>::Random(basis.size(), 2);

        gsTHBSpline<2> geom(basis, coefs);
        gsMatrix<> pts(2, 5);
        pts << 0.1, 0.3, 0.5, 0.7, 0.9,
               0.2, 0.4, 0.6, 0.8, 0.1;
        gsMatrix<> vals0;
        geom.eval_into(pts, vals0);

        gsTHBSplineBasis<2> basisRef = basis;
        gsMatrix<> coefsRef = coefs;
        basisRef.refineToLevel_withCoefs(1, coefsRef);
        gsTHBSpline<2> geomRef(basisRef, coefsRef);
        gsMatrix<> vals1;
        geomRef.eval_into(pts, vals1);

        CHECK((vals0 - vals1).array().abs().maxCoeff() < 1e-10);
    }
}
