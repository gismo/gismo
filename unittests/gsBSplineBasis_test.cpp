/** @file gsBSplineBasis_test.cpp

    @brief Degree-0 derivative evaluation and derivative consistency of univariate B-spline bases

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#include "gismo_unittest.h"

SUITE(gsBSplineBasis_test)
{
    // Degree-0 basis with 4 elements of width 1/4 and a point in each element
    struct P0
    {
        P0() : kv(0, 1, 3, 1), b(kv), u(1, 4)
        {
            u << 0.125, 0.375, 0.625, 0.875;
        }
        gsKnotVector<real_t> kv;
        gsBSplineBasis<real_t> b;
        gsMatrix<real_t> u;
    };

    // Derivatives of a piecewise constant vanish inside every element
    void checkAllZero(const gsMatrix<real_t> & r)
    {
        CHECK_EQUAL(1, r.rows());
        CHECK_EQUAL(4, r.cols());
        for (index_t i = 0; i != r.rows(); ++i)
            for (index_t j = 0; j != r.cols(); ++j)
                CHECK(r(i, j) == 0);
    }

    TEST(bsplP0DerivInto)
    {
        P0 d;
        gsMatrix<real_t> r;
        d.b.deriv_into(d.u, r);
        checkAllZero(r);
    }

    TEST(bsplP0DerivSingleInto)
    {
        P0 d;
        for (index_t i = 0; i != 4; ++i)
        {
            gsMatrix<real_t> r;
            d.b.derivSingle_into(i, d.u, r);
            checkAllZero(r);
        }
    }

    TEST(bsplP0Deriv2Into)
    {
        P0 d;
        gsMatrix<real_t> r;
        d.b.deriv2_into(d.u, r);
        checkAllZero(r);
    }

    TEST(bsplP0Deriv2SingleInto)
    {
        P0 d;
        for (index_t i = 0; i != 4; ++i)
        {
            gsMatrix<real_t> r;
            d.b.deriv2Single_into(i, d.u, r);
            checkAllZero(r);
        }
    }

    TEST(bsplP0EvalAllDersInto)
    {
        P0 d;
        std::vector<gsMatrix<real_t> > r;
        d.b.evalAllDers_into(d.u, 2, r);
        CHECK_EQUAL((size_t)3, r.size());
        if (3 != r.size()) return;
        for (size_t k = 0; k != 3; ++k)
        {
            CHECK_EQUAL(1, r[k].rows());
            CHECK_EQUAL(4, r[k].cols());
        }
        for (index_t j = 0; j != r[0].cols(); ++j)
            CHECK(r[0](0, j) == 1);
        for (size_t k = 1; k != 3; ++k)
            for (index_t j = 0; j != r[k].cols(); ++j)
                CHECK(r[k](0, j) == 0);
    }

    TEST(bsplFdDerivSingle)
    {
        // Points lie 0.125 >> 2h from any knot, so x +- h stays in one element
        // where B_i is a polynomial of degree p. The central difference is exact
        // for p <= 2 and has error h^2/6 |B_i'''| ~ 1e-10 for p = 3; rounding is
        // ~eps/h ~ 2e-10. The tolerance 1e-7 is far below the O(1) size of a
        // wrong derivative.
        const real_t h = 1e-6;
        const real_t xs[4] = {0.125, 0.375, 0.625, 0.875};
        for (short_t p = 1; p <= 3; ++p)
        {
            gsKnotVector<real_t> kv(0, 1, 3, p + 1);
            gsBSplineBasis<real_t> b(kv);
            for (index_t i = 0; i != b.size(); ++i)
                for (int k = 0; k != 4; ++k)
                {
                    gsMatrix<real_t> x(1, 1), xp(1, 1), xm(1, 1), d, fp, fm;
                    x << xs[k];
                    xp << xs[k] + h;
                    xm << xs[k] - h;
                    b.derivSingle_into(i, x, d);
                    b.evalSingle_into(i, xp, fp);
                    b.evalSingle_into(i, xm, fm);
                    CHECK_EQUAL(1, d.rows());
                    CHECK_EQUAL(1, d.cols());
                    CHECK_CLOSE((fp(0, 0) - fm(0, 0)) / (2 * h), d(0, 0), 1e-7);
                }
        }
    }
}
