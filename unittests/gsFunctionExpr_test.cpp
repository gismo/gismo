/** @file gsFunctionExpr_test.cpp

    @brief Tests derivative evaluation of gsFunctionExpr.

    The focus is deriv2_into(), whose result is a *packed* column holding the
    upper triangle of the Hessian in the order documented for
    gsFunction::deriv2_into(),

        ( d_xx, d_yy, d_zz, d_xy, d_xz, d_yz )

    i.e. the d diagonal entries first, then the d(d-1)/2 off-diagonal ones in
    lexicographic (k,l) order. That layout is only exercised non-trivially for
    d >= 3: in 2D there is a single off-diagonal entry, so any indexing scheme
    happens to be correct and 2D tests cannot see a packing defect.

    gsFunctionExpr has three independent ways to reach the same second
    derivatives -- deriv2_into(), hess() and mderiv() -- plus laplacian() for
    the trace. The cross-checks below play them against each other and against
    a finite-difference Hessian built here from eval() alone, so a defect in
    any single one of them cannot hide.

    Run only these tests with:  ./unittests gsFunctionExpr_test

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
 **/

#include "gismo_unittest.h"

SUITE(gsFunctionExpr_test)
{
    // Without the optional autodiff module, gsFunctionExpr differentiates by
    // central finite differences with h = 1e-5, so second derivatives carry
    // roughly 1e-5 of absolute roundoff at the magnitudes used here (the
    // O(h^-2) amplification of cancellation, not truncation -- the polynomial
    // fixtures are differentiated exactly). The fixtures are chosen so that
    // every entry of the packed column differs from every other by at least
    // 1.0, which leaves this tolerance ~200x away from ever accepting a
    // mis-packed or unwritten entry, and ~200x away from a spurious failure.
    const real_t tol = 5e-3;

    // Hessian of component \a comp by central differences on eval() alone.
    // Independent of deriv2_into/hess/mderiv, hence usable as their oracle.
    gsMatrix<real_t> hessianByFD(const gsFunction<real_t> & f,
                                 const gsMatrix<real_t> & pt,
                                 index_t comp = 0)
    {
        const index_t d = pt.rows();
        const real_t h = 1e-4;
        gsMatrix<real_t> H(d,d), q(d,1);

        for (index_t k = 0; k!=d; ++k)
        {
            for (index_t l = k; l!=d; ++l)
            {
                real_t acc = 0;
                // 4-point stencil for k!=l; collapses to the 3-point second
                // difference for k==l, where the (+h,+h)/(-h,-h) samples merge.
                for (int sk = -1; sk <= 1; sk += 2)
                    for (int sl = -1; sl <= 1; sl += 2)
                    {
                        q = pt;
                        q(k,0) += sk*h;
                        q(l,0) += sl*h;
                        acc += sk*sl*f.eval(q)(comp,0);
                    }
                H(k,l) = H(l,k) = acc / (4*h*h);
            }
        }
        return H;
    }

    // Repack a dense Hessian into the layout deriv2_into() must produce.
    gsVector<real_t> packHessian(const gsMatrix<real_t> & H)
    {
        const index_t d = H.rows();
        gsVector<real_t> v(d + d*(d-1)/2);
        index_t m = d;
        for (index_t k = 0; k!=d; ++k)
        {
            v[k] = H(k,k);
            for (index_t l = k+1; l!=d; ++l)
                v[m++] = H(k,l);
        }
        return v;
    }

    gsMatrix<real_t> pointAt(real_t x, real_t y, real_t z)
    {
        gsMatrix<real_t> p(3,1);
        p << x, y, z;
        return p;
    }

    // f = x^2 y + 3 y^2 z + 5 z^2 x + 7 x y z, whose Hessian is
    //   d_xx = 2y,  d_yy = 6z,  d_zz = 10x,
    //   d_xy = 2x + 7z,  d_xz = 10z + 7y,  d_yz = 6y + 7x.
    const std::string cubic3d = "x^2*y + 3*y^2*z + 5*z^2*x + 7*x*y*z";

    TEST(deriv2_packing_3d)
    {
        gsFunctionExpr<real_t> f(cubic3d, 3);
        gsMatrix<real_t> u(3,2);
        u.col(0) = pointAt(1.0, 2.0, 3.0);
        u.col(1) = pointAt(0.5,-1.5, 2.0);

        gsMatrix<real_t> D;
        f.deriv2_into(u, D);
        CHECK_EQUAL(6, D.rows());   // stride = 3 + 3 = 6 for one component
        CHECK_EQUAL(2, D.cols());

        // (1,2,3): entries 4, 18, 10, 23, 44, 19 -- all distinct
        const real_t e0[6] = {4, 18, 10, 23, 44, 19};
        CHECK_ARRAY_CLOSE(e0, D.col(0).data(), 6, tol);

        // (0.5,-1.5,2): entries -3, 12, 5, 15, 9.5, -5.5 -- all distinct
        const real_t e1[6] = {-3, 12, 5, 15, 9.5, -5.5};
        CHECK_ARRAY_CLOSE(e1, D.col(1).data(), 6, tol);
    }

    TEST(deriv2_packing_2d)
    {
        // f = x^3 y + 2 x y^2:  d_xx = 6xy, d_yy = 4x, d_xy = 3x^2 + 4y
        gsFunctionExpr<real_t> f("x^3*y + 2*x*y^2", 2);
        gsMatrix<real_t> u(2,1);
        u << 1.0, 2.0;

        gsMatrix<real_t> D;
        f.deriv2_into(u, D);
        CHECK_EQUAL(3, D.rows());   // stride = 2 + 1

        const real_t e[3] = {12, 4, 11};
        CHECK_ARRAY_CLOSE(e, D.col(0).data(), 3, tol);
    }

    TEST(deriv2_multicomponent_3d)
    {
        // Two components must land at offsets 0 and stride, each packed
        // independently. g = 2xyz has a hollow diagonal, so a component
        // offset that slipped would show up as zeros in the wrong block.
        gsFunctionExpr<real_t> f(cubic3d, "2*x*y*z", 3);
        gsMatrix<real_t> u = pointAt(1.0, 2.0, 3.0);

        gsMatrix<real_t> D;
        f.deriv2_into(u, D);
        CHECK_EQUAL(12, D.rows());  // 2 components x stride 6

        const real_t ef[6] = {4, 18, 10, 23, 44, 19};
        const real_t eg[6] = {0, 0, 0, 6, 4, 2};   // d_xy=2z, d_xz=2y, d_yz=2x
        CHECK_ARRAY_CLOSE(ef, D.col(0).data(),     6, tol);
        CHECK_ARRAY_CLOSE(eg, D.col(0).data() + 6, 6, tol);
    }

    TEST(deriv2_matches_hess_3d)
    {
        // hess() reaches the same second derivatives through the same
        // differentiation but without a packing counter, so a disagreement
        // here isolates the packing; agreement while deriv2_matches_fd_3d
        // fails isolates the differentiation.
        gsFunctionExpr<real_t> f("exp(x)*sin(y) + x*cos(z) + sin(y*z)", 3);
        gsMatrix<real_t> u = pointAt(0.3, 0.7, 1.1);

        gsMatrix<real_t> D;
        f.deriv2_into(u, D);
        const gsVector<real_t> expected = packHessian(f.hess(u));

        CHECK_ARRAY_CLOSE(expected.data(), D.col(0).data(), 6, tol);
    }

    TEST(deriv2_matches_fd_3d)
    {
        // Oracle built from eval() only, on a wider stencil than either
        // internal route uses. Without the autodiff module this bounds the
        // packing and the differencing step against an independent stencil;
        // with it enabled the oracle becomes genuinely independent of the
        // differentiation itself, and this is the check that would catch an
        // exprtk adaptor dropping a cross term in a nested dual -- a defect
        // that leaves first derivatives correct and silently corrupts
        // second ones.
        gsFunctionExpr<real_t> f("exp(x)*sin(y) + x*cos(z) + sin(y*z)", 3);
        gsMatrix<real_t> u = pointAt(0.3, 0.7, 1.1);

        gsMatrix<real_t> D;
        f.deriv2_into(u, D);
        const gsVector<real_t> expected = packHessian(hessianByFD(f, u));

        CHECK_ARRAY_CLOSE(expected.data(), D.col(0).data(), 6, tol);
    }

    TEST(deriv_first_3d)
    {
        // Separates "first derivatives broken" from "packing broken": the
        // gradient uses the same seeding machinery but no packed layout.
        gsFunctionExpr<real_t> f(cubic3d, 3);
        gsMatrix<real_t> u = pointAt(1.0, 2.0, 3.0);

        gsMatrix<real_t> D;
        f.deriv_into(u, D);
        CHECK_EQUAL(3, D.rows());

        // d_x = 2xy + 5z^2 + 7yz, d_y = x^2 + 6yz + 7xz, d_z = 3y^2 + 10zx + 7xy
        const real_t e[3] = {4 + 45 + 42, 1 + 36 + 21, 12 + 30 + 14};
        CHECK_ARRAY_CLOSE(e, D.col(0).data(), 3, tol);
    }

    TEST(laplacian_matches_deriv2_trace_3d)
    {
        // The Laplacian must equal the sum of the d leading packed entries.
        // Reading past them would pull in an off-diagonal.
        gsFunctionExpr<real_t> f(cubic3d, 3);
        gsMatrix<real_t> u = pointAt(1.0, 2.0, 3.0);

        gsMatrix<real_t> D;
        f.deriv2_into(u, D);

        const gsMatrix<real_t> lap = f.laplacian(u);
        CHECK_EQUAL(1, lap.rows());
        CHECK_CLOSE(D.col(0).head(3).sum(), lap(0,0), tol);
        CHECK_CLOSE(32.0, lap(0,0), tol);   // 2y + 6z + 10x at (1,2,3)
    }

    TEST(mderiv_matches_deriv2_3d)
    {
        // mderiv(u,k,l) is a fourth route to the mixed entries; it must agree
        // with the off-diagonal block of the packed column.
        gsFunctionExpr<real_t> f(cubic3d, 3);
        gsMatrix<real_t> u = pointAt(1.0, 2.0, 3.0);

        gsMatrix<real_t> D;
        f.deriv2_into(u, D);

        index_t m = 3;
        for (index_t k = 0; k!=3; ++k)
            for (index_t l = k+1; l!=3; ++l)
            {
                memory::unique_ptr< gsMatrix<real_t> > md(f.mderiv(u,k,l));
                CHECK_CLOSE(D(m++,0), (*md)(0,0), tol);
            }
    }
}
