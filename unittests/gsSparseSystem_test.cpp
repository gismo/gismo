/** @file gsSparseSystem_test.cpp

    @brief Tests for gsSparseSystem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M.M. Hodzelmans
**/

#include "gismo_unittest.h"

SUITE(gsSparseSystem_test)
{
    // push() of several blocks at once checks the dofs of a column block
    // against the mapper of that column block, and adds the local rhs once
    // per row block: the same as pushing the blocks one at a time.
    TEST(PushBlocksWithOwnMappers)
    {
        gsVector<index_t> sz(1);
        sz << 3;
        gsDofMapper m0(sz), m1(sz);
        m1.eliminateDof(0, 0);
        m0.finalize();
        m1.finalize();

        std::vector<gsMatrix<index_t> > actives(2, gsMatrix<index_t>(3, 1));
        for (index_t i = 0; i != 3; ++i)
        {
            actives[0](i, 0) = m0.index(i);
            actives[1](i, 0) = m1.index(i);
        }
        std::vector<gsMatrix<real_t> > fixed(2);
        fixed[0].setZero(0, 1);
        fixed[1].setConstant(1, 1, 2.0);

        gsMatrix<real_t> localMat(6, 6), localRhs(6, 1);
        for (index_t i = 0; i != 6; ++i)
        {
            localRhs(i, 0) = i + 1;
            for (index_t j = 0; j != 6; ++j)
                localMat(i, j) = 1 + i + 10 * j;
        }

        std::vector<gsDofMapper> mappers;
        mappers.push_back(m0);
        mappers.push_back(m1);
        gsSparseSystem<real_t> sys(mappers, 2, 2);
        sys.reserve(6, 1);
        sys.push(localMat, localRhs, actives, fixed);

        mappers.clear();
        mappers.push_back(m0);
        mappers.push_back(m1);
        gsSparseSystem<real_t> ref(mappers, 2, 2);
        ref.reserve(6, 1);
        const gsMatrix<real_t> noRhs = gsMatrix<real_t>::Zero(3, 1);
        for (index_t r = 0; r != 2; ++r)
            for (index_t c = 0; c != 2; ++c)
            {
                const gsMatrix<real_t> block = localMat.block(3 * r, 3 * c, 3, 3);
                const gsMatrix<real_t> rhs = (0 == c ? gsMatrix<real_t>(localRhs.middleRows(3 * r, 3)) : noRhs);
                ref.push(block, rhs, actives[r], actives[c], fixed[c], r, c);
            }

        CHECK_EQUAL(5, sys.matrix().rows());
        CHECK_EQUAL(5, sys.matrix().cols());
        const gsMatrix<real_t> A = sys.matrix(), Aref = ref.matrix();
        CHECK(Aref.norm() > 0);
        CHECK((A - Aref).norm() <= 1e-12 * Aref.norm());
        CHECK((sys.rhs() - ref.rhs()).norm() <= 1e-12 * ref.rhs().norm());
    }
}
