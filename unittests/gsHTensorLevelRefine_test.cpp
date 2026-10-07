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


// Two THB patches side by side; patch 1 is refined along its west side.
gsMultiBasis<real_t> twoPatchTHB()
{
    gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
    gsMultiBasis<real_t> mbT(mp);
    mbT.degreeElevate(1);
    mbT.uniformRefine(1);
    gsMultiBasis<real_t> mb;
    for (size_t p = 0; p != mbT.nBases(); ++p)
        mb.addBasis(new gsTHBSplineBasis<2, real_t>(
            static_cast<const gsTensorBSplineBasis<2, real_t> &>(mbT.basis(p))));
    mb.setTopology(mp.topology());
    std::vector<index_t> box = {1, 0, 0, 2, 4};
    mb.basis(1).refineElements(box);
    return mb;
}

// Lower corners and patch indices from it to end.
std::vector<std::pair<gsVector<real_t>, index_t> >
visit(gsDomainIteratorWrapper<real_t> & it, const gsDomainIteratorWrapper<real_t> & end)
{
    std::vector<std::pair<gsVector<real_t>, index_t> > out;
    for (; it != end; ++it)
        out.emplace_back(it.lowerCorner(), it.patchIndex());
    return out;
}

// A copy made partway through must visit the same remaining elements after the
// original has moved on and been destroyed.
void checkCopy(gsDomainIteratorWrapper<real_t> it, const gsDomainIteratorWrapper<real_t> & end,
               index_t patch)
{
    ++it;
    gsDomainIteratorWrapper<real_t> copy(it);
    std::vector<std::pair<gsVector<real_t>, index_t> > ref;
    {
        gsDomainIteratorWrapper<real_t> orig(give(it));
        ref = visit(orig, end);
    }
    const std::vector<std::pair<gsVector<real_t>, index_t> > got = visit(copy, end);
    CHECK(ref.size() > 1);
    CHECK_EQUAL(ref.size(), got.size());
    for (size_t i = 0; i != std::min(ref.size(), got.size()); ++i)
    {
        CHECK((ref[i].first - got[i].first).norm() < 1e-14);
        CHECK_EQUAL(patch, got[i].second);
    }
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

    TEST(hdomain_iterator_copies)
    {
        gsMultiBasis<real_t> mb = twoPatchTHB();
        typename gsDomain<real_t>::Ptr sub = mb.domain()->subdomain(1);
        checkCopy(sub->beginAll(), sub->endAll(), 1);
        checkCopy(sub->beginBdr(boundary::west), sub->endBdr(boundary::west), 1);
    }
}
