/** @file gsGreville_test.cpp

    @brief Greville abscissae at the domain end and hierarchical bases at the upper corner

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#include "gismo_unittest.h"

SUITE(gsGreville_test)
{
    // Last Greville point of a knot vector on [0,1] must be exactly the domain end
    void checkLast(const gsKnotVector<real_t> & kv, short_t p)
    {
        const index_t i = kv.size() - p - 2;
        gsMatrix<real_t> g;
        kv.greville_into(g);
        CHECK_EQUAL((index_t)(kv.size() - p - 1), g.cols());
        CHECK_EQUAL(1.0, g(0, i));
        CHECK_EQUAL(1.0, kv.grevilleOf(i));
        CHECK(kv.grevilleOf(i) == g(0, i));
    }

    TEST(grevLastOpen)
    {
        for (short_t p = 1; p <= 3; ++p)
        {
            gsKnotVector<real_t> kv(0, 1, 3, p + 1, 1, p);
            checkLast(kv, p);
        }
    }

    TEST(grevLastCm1)
    {
        for (short_t p = 1; p <= 3; ++p)
        {
            gsKnotVector<real_t> kv(0, 1, 3, p + 1, p + 1, p);
            checkLast(kv, p);
        }
    }

    // Hierarchical basis refined twice towards the upper corner:
    // level 1 on [0.5,1]^d and level 2 on [0.75,1]^d
    template <bool Trunc>
    gsTHBSplineBasis<2, real_t, Trunc> makeBasis2(short_t p)
    {
        gsKnotVector<real_t> kv(0, 1, 3, p + 1);
        gsTensorBSplineBasis<2, real_t> tb(kv, kv);
        gsTHBSplineBasis<2, real_t, Trunc> hb(tb);
        std::vector<index_t> b1(5), b2(5);
        b1[0] = 1; b1[1] = 4; b1[2] = 4; b1[3] = 8; b1[4] = 8;
        b2[0] = 2; b2[1] = 12; b2[2] = 12; b2[3] = 16; b2[4] = 16;
        hb.refineElements(b1);
        hb.refineElements(b2);
        return hb;
    }

    template <bool Trunc>
    gsTHBSplineBasis<3, real_t, Trunc> makeBasis3(short_t p)
    {
        gsKnotVector<real_t> kv(0, 1, 3, p + 1);
        gsTensorBSplineBasis<3, real_t> tb(kv, kv, kv);
        gsTHBSplineBasis<3, real_t, Trunc> hb(tb);
        std::vector<index_t> b1(7), b2(7);
        b1[0] = 1; b2[0] = 2;
        for (int k = 0; k != 3; ++k)
        {
            b1[1 + k] = 4;  b1[4 + k] = 8;
            b2[1 + k] = 12; b2[4 + k] = 16;
        }
        hb.refineElements(b1);
        hb.refineElements(b2);
        return hb;
    }

    template <short_t d, bool Trunc>
    void checkCorner(const gsTHBSplineBasis<d, real_t, Trunc> & hb)
    {
        CHECK(hb.maxLevel() >= 2);
        gsMatrix<real_t> x(d, 1), val;
        gsMatrix<index_t> act;
        x.setOnes();
        hb.active_into(x, act);
        CHECK(act.rows() > 0);
        hb.eval_into(x, val);
        CHECK_EQUAL(act.rows(), val.rows());
        CHECK_EQUAL(1, val.cols());
        if (val.size() == 0) return;
        // The finest-level open-knot corner B-spline interpolates 1 there
        CHECK_CLOSE(1.0, val.maxCoeff(), 1e-12);
        if (Trunc)
            CHECK_CLOSE(1.0, val.sum(), 1e-12);
    }

    // Every anchor lies in the support of its own function (and the THB values
    // form a partition of unity there). A truncated function may vanish at its
    // anchor, so positivity is only asserted for the untruncated HB case.
    template <short_t d, bool Trunc>
    void checkAnchors(const gsTHBSplineBasis<d, real_t, Trunc> & hb)
    {
        CHECK(hb.maxLevel() >= 2);
        gsMatrix<real_t> A;
        hb.anchors_into(A);
        CHECK_EQUAL(hb.size(), A.cols());
        for (index_t j = 0; j != A.cols(); ++j)
        {
            gsMatrix<real_t> x = A.col(j), val;
            gsMatrix<index_t> act;
            hb.active_into(x, act);
            index_t r = -1;
            for (index_t k = 0; k != act.rows(); ++k)
                if (act(k, 0) == j) r = k;
            CHECK(r >= 0);
            hb.eval_into(x, val);
            CHECK_EQUAL(act.rows(), val.rows());
            if (r < 0 || val.rows() != act.rows()) continue;
            if (!Trunc)
                CHECK(val(r, 0) > 0);
            else
                CHECK_CLOSE(1.0, val.sum(), 1e-12);
        }
    }

    TEST(grevHierCorner2d)
    {
        for (short_t p = 1; p <= 3; ++p)
        {
            checkCorner<2, false>(makeBasis2<false>(p));
            checkCorner<2, true>(makeBasis2<true>(p));
        }
    }

    TEST(grevHierCorner3d)
    {
        checkCorner<3, false>(makeBasis3<false>(2));
        checkCorner<3, true>(makeBasis3<true>(2));
    }

    TEST(grevHierAnchors2d)
    {
        for (short_t p = 1; p <= 3; ++p)
        {
            checkAnchors<2, false>(makeBasis2<false>(p));
            checkAnchors<2, true>(makeBasis2<true>(p));
        }
    }

    TEST(grevHierAnchors3d)
    {
        checkAnchors<3, false>(makeBasis3<false>(2));
        checkAnchors<3, true>(makeBasis3<true>(2));
    }
}
