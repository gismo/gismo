/** @file gsQuasiInterpolate_test.cpp

    @brief Tests the quasi-interpolation schemes in gsQuasiInterpolate.

    Covered methods (each identified by a test-name prefix):
      * taylor_*      : localTaylor / Taylor / Taylor2D (Taylor QI), defined
                        for tensor B-spline bases only. The
                        tensor-product Taylor QI Q_{p,p} reproduces every
                        spline in its tensor B-spline space (Lyche-Morken
                        Thm 8.5). Verified by exact polynomial reproduction
                        (order-<=2 regime, gsFunctionExpr) and by exact
                        reproduction of an arbitrary tensor B-spline geometry
                        in 2D/3D (needs mixed partials to order d*p, which
                        spline geometries provide but gsFunctionExpr does not).
      * intpl_*       : localIntpl  (local interpolation QI); projector tests
                        on THB, HB, rational THB and NURBS bases
      * l2_*          : localL2     (local L2-projection QI); projector tests
                        on HB, THB and rational THB bases
      * schoenberg_*  : Schoenberg  (variation-diminishing, reproduces affine)
      * evalbased_*   : EvalBased   (evaluation-based, 1D)

    Run only the Taylor tests with:  ./unittests taylor_

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
 **/

#include "gismo_unittest.h"

SUITE(gsQuasiInterpolate_test)
{
    // Reproduction is exact in exact arithmetic; the sampled-L2 metric
    // accumulates per-point roundoff (~1e-8 over the sample grid), so the
    // tolerance is set well below any genuine defect (stride/aliasing bugs
    // produce O(1e-1) errors).
    const real_t tol = 1e-6;

    // Coefficients are compared directly (no sampling) and every local solve
    // is a small dense (p+1)^d LU, so the projector checks can be far tighter.
    const real_t projTol = 1e-12;

    // Level-to-level roundoff accumulates in the residual sweeps; the
    // tolerance stays far below the O(1e-1) error of a hierarchy-unaware QI.
    const real_t hbTol = 1e-11;

    // The local L2 functional solves a Gram system B W B^T whose condition
    // number is the square of that of the local collocation matrix; coefficient
    // round-off reaches ~3e-11 (THB, p=3) and ~3e-10 (rational THB, p=3), far
    // below the O(1e-1) error of a hierarchy-unaware scheme.
    const real_t l2Tol = 1e-8;

    // L2-type sampled error between two functions over the spline's domain.
    real_t sampleError(const gsFunction<real_t> & fun, const gsGeometry<real_t> & spl)
    {
        gsGridIterator<real_t,CUBE> pt(spl.support(), 15);
        real_t error = 0;
        for (; pt; ++pt)
            error += (fun.eval(*pt) - spl.eval(*pt)).squaredNorm();
        return math::sqrt(error);
    }

    // Build a geometry from coefs and measure its distance to fun.
    real_t errVsFun(const gsBasis<real_t> & basis, const gsFunction<real_t> & fun,
                    gsMatrix<real_t> & coefs)
    {
        gsGeometry<real_t>::uPtr geo = basis.makeGeometry(give(coefs));
        return sampleError(fun, *geo);
    }

    // Deterministic, varied control points (n x targetDim).
    gsMatrix<real_t> makeCoefs(index_t n, index_t targetDim)
    {
        gsMatrix<real_t> c(n, targetDim);
        for (index_t i = 0; i!=n; ++i)
            for (index_t j = 0; j!=targetDim; ++j)
                c(i,j) = std::sin(0.7*i + 1.3*j) + 0.15*i - 0.4*j;
        return c;
    }

    // Reproducible pseudo-random matrix with entries in [lo,hi]: 64-bit LCG
    // (Knuth MMIX constants), top 53 bits mapped to [0,1). Independent of the
    // global std::rand state, so results do not depend on suite order.
    gsMatrix<real_t> lcgMatrix(index_t rows, index_t cols, real_t lo, real_t hi,
                               uint64_t seed)
    {
        gsMatrix<real_t> m(rows, cols);
        for (index_t j = 0; j != cols; ++j)
            for (index_t i = 0; i != rows; ++i)
            {
                seed = 6364136223846793005ULL * seed + 1442695040888963407ULL;
                m(i,j) = lo + (hi - lo) * static_cast<real_t>(seed >> 11)
                                        / static_cast<real_t>(1ULL << 53);
            }
        return m;
    }

    // Quasi-interpolate the spline with coefficients c (in basis b) in b itself.
    gsMatrix<real_t> qiOfSpline(const gsBasis<real_t> & b, const gsMatrix<real_t> & c)
    {
        gsGeometry<real_t>::uPtr f = b.makeGeometry(c);
        gsMatrix<real_t> res;
        gsQuasiInterpolate<real_t>::localIntpl(b, *f, res);
        return res;
    }

    // ================================ Taylor ================================

    // 1D deg 2: reproduce constants, linears, quadratics (a gsBSplineBasis is a
    // gsTensorBSplineBasis<1>, exercising the d=1 path of localTaylor).
    TEST(taylor_reproduction_1D)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 3); // deg 2
        gsBSplineBasis<real_t> basis(kv);

        gsFunctionExpr<real_t> f1("1",         1);
        gsFunctionExpr<real_t> f2("2*x - 1",   1);
        gsFunctionExpr<real_t> f3("x^2 - 3*x", 1);

        for (const gsFunctionExpr<real_t>* f : {&f1,&f2,&f3})
        {
            gsMatrix<real_t> coefs;
            gsQuasiInterpolate<real_t>::localTaylor(basis, *f, basis.degree(0), coefs);
            CHECK_CLOSE(0.0, errVsFun(basis, *f, coefs), tol);
        }
    }

    // 2D deg 1: reproduce the full bilinear space {1,x,y,x*y}. The mixed term
    // x*y (order-2 derivative) is exactly what the old aliased code and the
    // targetDim stride bug got wrong; here total order = 2, so an analytic
    // gsFunctionExpr can supply the needed derivatives. Scalar and vector.
    TEST(taylor_reproduction_2D_bilinear)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 2); // deg 1
        gsTensorBSplineBasis<2,real_t> basis(kv, kv);

        gsFunctionExpr<real_t> fxy  ("x*y", 2);
        gsFunctionExpr<real_t> ffull("1 - 2*x + 3*y + 4*x*y", 2);
        gsFunctionExpr<real_t> fvec ("x*y", "1 - x + 2*y - 3*x*y", 2);

        for (const gsFunctionExpr<real_t>* f : {&fxy,&ffull,&fvec})
        {
            gsMatrix<real_t> coefs;
            gsQuasiInterpolate<real_t>::localTaylor(basis, *f, basis.degree(0), coefs);
            CHECK_CLOSE(0.0, errVsFun(basis, *f, coefs), tol);
        }
    }

    // 2D deg 2: reproduce an arbitrary tensor B-spline (needs mixed partials up
    // to order 4, incl. the terms the old code aliased). Scalar and vector.
    TEST(taylor_reproduction_2D_spline)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 4, 3); // deg 2
        gsTensorBSplineBasis<2,real_t> basis(kv, kv);

        for (index_t td : {1, 3})
        {
            gsMatrix<real_t> coefs = makeCoefs(basis.size(), td);
            gsTensorBSpline<2,real_t> g(basis, coefs);
            gsMatrix<real_t> qcoefs;
            gsQuasiInterpolate<real_t>::localTaylor(basis, g, basis.degree(0), qcoefs);
            CHECK_CLOSE(0.0, (qcoefs - coefs).cwiseAbs().maxCoeff(), tol);
        }
    }

    // 3D deg 2: reproduce an arbitrary tensor B-spline. Needs mixed partials up
    // to order 6 and exercises the order-3+ composition-order extraction in
    // derivRow (the (1,1,1)-type partials).
    TEST(taylor_reproduction_3D_spline)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 2, 3); // deg 2
        gsTensorBSplineBasis<3,real_t> basis(kv, kv, kv);

        for (index_t td : {1, 3})
        {
            gsMatrix<real_t> coefs = makeCoefs(basis.size(), td);
            gsTensorBSpline<3,real_t> g(basis, coefs);
            gsMatrix<real_t> qcoefs;
            gsQuasiInterpolate<real_t>::localTaylor(basis, g, basis.degree(0), qcoefs);
            CHECK_CLOSE(0.0, (qcoefs - coefs).cwiseAbs().maxCoeff(), tol);
        }
    }

    // Whole-basis 1D Taylor (the dedicated Taylor(...) overload): reproduce a
    // quadratic on a deg-2 basis.
    TEST(taylor_1D_wholebasis)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 3); // deg 2
        gsBSplineBasis<real_t> basis(kv);

        gsFunctionExpr<real_t> f("x^2 - 3*x + 1", 1);
        gsMatrix<real_t> coefs;
        gsQuasiInterpolate<real_t>::Taylor(basis, f, basis.degree(0), coefs);
        CHECK_CLOSE(0.0, errVsFun(basis, f, coefs), tol);
    }

    // Taylor2D (now delegating to the dimension-independent localTaylor):
    // reproduce an arbitrary 2D tensor B-spline.
    TEST(taylor2d_reproduction_spline)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 4, 3); // deg 2
        gsTensorBSplineBasis<2,real_t> basis(kv, kv);

        gsMatrix<real_t> coefs = makeCoefs(basis.size(), 1);
        gsTensorBSpline<2,real_t> g(basis, coefs);
        gsMatrix<real_t> qcoefs;
        gsQuasiInterpolate<real_t>::Taylor2D(basis, g, basis.degree(0), qcoefs);
        CHECK_CLOSE(0.0, (qcoefs - coefs).cwiseAbs().maxCoeff(), tol);
    }

    // ============================= localIntpl ==============================

    // Local interpolation reproduces polynomials up to the basis degree.
    TEST(intpl_reproduction_2D)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 3); // deg 2
        gsTensorBSplineBasis<2,real_t> basis(kv, kv);

        gsFunctionExpr<real_t> f("1 + x - 2*y + 3*x*y + x^2 - y^2", 2);
        gsMatrix<real_t> coefs;
        gsQuasiInterpolate<real_t>::localIntpl(basis, f, coefs);
        CHECK_CLOSE(0.0, errVsFun(basis, f, coefs), tol);
    }

    TEST(intpl_reproduction_THB)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 3); // deg 2
        gsTensorBSplineBasis<2,real_t> tbasis(kv, kv);
        gsTHBSplineBasis<2,real_t> basis(tbasis);
        gsMatrix<real_t> box(2,2);
        box.col(0) << 0.0, 0.0;
        box.col(1) << 0.5, 0.5;
        basis.refine(box);

        gsFunctionExpr<real_t> f("1 + x - 2*y + 3*x*y + x^2 - y^2", 2);
        gsMatrix<real_t> coefs;
        gsQuasiInterpolate<real_t>::localIntpl(basis, f, coefs);
        CHECK_CLOSE(0.0, errVsFun(basis, f, coefs), tol);
    }

    // localIntpl is a projector on truncated hierarchical spaces: quasi-
    // interpolating a spline of the basis returns its coefficients. With a
    // local QI on each tensor level whose functionals for level-l functions
    // are supported in Omega^l \ Omega^{l+1}, the hierarchical QI over the
    // truncated basis reproduces the whole THB space (H. Speleers, C. Manni,
    // "Effortless quasi-interpolation in hierarchical spaces", Numer. Math.
    // 132 (2016) 155-184).
    TEST(intpl_projector_THB_2D)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsTHBSplineBasis<2,real_t> basis(tb);
        gsMatrix<real_t> box(2,2);
        box.col(0) << 0.0, 0.0;
        box.col(1) << 0.5, 0.5;
        basis.refine(box);
        CHECK(basis.size() > tb.size());

        const gsMatrix<real_t> c = lcgMatrix(basis.size(), 1, -1.0, 1.0, 11);
        const gsMatrix<real_t> res = qiOfSpline(basis, c);
        CHECK_EQUAL(c.rows(), res.rows());
        CHECK_EQUAL(c.cols(), res.cols());
        if (res.rows() == c.rows() && res.cols() == c.cols())
            CHECK_CLOSE(0.0, (res - c).cwiseAbs().maxCoeff(), projTol);
    }

    // Projector property on the boundary basis (1D THB of mixed level).
    TEST(intpl_projector_THB_2D_boundary)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsTHBSplineBasis<2,real_t> thb(tb);
        gsMatrix<real_t> box(2,2);
        box.col(0) << 0.0, 0.0;
        box.col(1) << 0.5, 0.5;
        thb.refine(box);
        gsBasis<real_t>::uPtr basis = thb.boundaryBasis(boundary::west);
        CHECK(basis->size() > tb.boundaryBasis(boundary::west)->size());

        const gsMatrix<real_t> c = lcgMatrix(basis->size(), 1, -1.0, 1.0, 12);
        const gsMatrix<real_t> res = qiOfSpline(*basis, c);
        CHECK_EQUAL(c.rows(), res.rows());
        CHECK_EQUAL(c.cols(), res.cols());
        if (res.rows() == c.rows() && res.cols() == c.cols())
            CHECK_CLOSE(0.0, (res - c).cwiseAbs().maxCoeff(), projTol);
    }

    // Projector property on the 2D boundary basis of a 3D THB basis.
    TEST(intpl_projector_THB_3D_boundary)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
        gsTensorBSplineBasis<3,real_t> tb(kv, kv, kv);
        gsTHBSplineBasis<3,real_t> thb(tb);
        gsMatrix<real_t> box(3,2);
        box.col(0) << 0.0, 0.0, 0.0;
        box.col(1) << 0.5, 0.5, 0.5;
        thb.refine(box);
        gsBasis<real_t>::uPtr basis = thb.boundaryBasis(boundary::west);
        CHECK(basis->size() > tb.boundaryBasis(boundary::west)->size());

        const gsMatrix<real_t> c = lcgMatrix(basis->size(), 1, -1.0, 1.0, 13);
        const gsMatrix<real_t> res = qiOfSpline(*basis, c);
        CHECK_EQUAL(c.rows(), res.rows());
        CHECK_EQUAL(c.cols(), res.cols());
        if (res.rows() == c.rows() && res.cols() == c.cols())
            CHECK_CLOSE(0.0, (res - c).cwiseAbs().maxCoeff(), projTol);
    }

    // Projector property when the refined region is a single level-1 element
    // thick strip along the west side.
    TEST(intpl_projector_THB_2D_strip)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsTHBSplineBasis<2,real_t> basis(tb);
        gsMatrix<real_t> box(2,2);
        box.col(0) << 0.0, 0.0;
        box.col(1) << 0.125, 1.0;
        basis.refine(box);
        CHECK(basis.size() > tb.size());

        const gsMatrix<real_t> c = lcgMatrix(basis.size(), 1, -1.0, 1.0, 15);
        const gsMatrix<real_t> res = qiOfSpline(basis, c);
        CHECK_EQUAL(c.rows(), res.rows());
        CHECK_EQUAL(c.cols(), res.cols());
        if (res.rows() == c.rows() && res.cols() == c.cols())
            CHECK_CLOSE(0.0, (res - c).cwiseAbs().maxCoeff(), projTol);
    }

    // Rational THB spline f = sum c_i w_i T_i / sum w_k T_k with non-constant
    // weights: quasi-interpolation in the rational basis must return c.
    TEST(intpl_rationalTHB_2D)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsTHBSplineBasis<2,real_t> thb(tb);
        gsMatrix<real_t> box(2,2);
        box.col(0) << 0.0, 0.0;
        box.col(1) << 0.5, 0.5;
        thb.refine(box);
        CHECK(thb.size() > tb.size());

        const gsMatrix<real_t> w = lcgMatrix(thb.size(), 1, 0.5, 2.0, 21);
        CHECK(w.maxCoeff() - w.minCoeff() > 0.5);
        gsRationalTHBSplineBasis<2,real_t> rb(thb.clone().release(), w);

        const gsMatrix<real_t> c = lcgMatrix(rb.size(), 1, -1.0, 1.0, 22);
        const gsMatrix<real_t> res = qiOfSpline(rb, c);
        CHECK_EQUAL(c.rows(), res.rows());
        CHECK_EQUAL(c.cols(), res.cols());
        if (res.rows() == c.rows() && res.cols() == c.cols())
            CHECK_CLOSE(0.0, (res - c).cwiseAbs().maxCoeff(), projTol);
    }

    // Same on the 1D rational boundary basis with the restricted weights.
    TEST(intpl_rationalTHB_2D_boundary)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsTHBSplineBasis<2,real_t> thb(tb);
        gsMatrix<real_t> box(2,2);
        box.col(0) << 0.0, 0.0;
        box.col(1) << 0.5, 0.5;
        thb.refine(box);

        const gsMatrix<real_t> w = lcgMatrix(thb.size(), 1, 0.5, 2.0, 23);
        CHECK(w.maxCoeff() - w.minCoeff() > 0.5);
        gsRationalTHBSplineBasis<2,real_t> rb(thb.clone().release(), w);
        gsBasis<real_t>::uPtr basis = rb.boundaryBasis(boundary::west);
        CHECK(basis->size() > tb.boundaryBasis(boundary::west)->size());

        const gsMatrix<real_t> c = lcgMatrix(basis->size(), 1, -1.0, 1.0, 24);
        const gsMatrix<real_t> res = qiOfSpline(*basis, c);
        CHECK_EQUAL(c.rows(), res.rows());
        CHECK_EQUAL(c.cols(), res.cols());
        if (res.rows() == c.rows() && res.cols() == c.cols())
            CHECK_CLOSE(0.0, (res - c).cwiseAbs().maxCoeff(), projTol);
    }

    // Non-truncated hierarchical (HB) spaces. The scheme is level-by-level
    // quasi-interpolation: the coefficient of a level-l function is the
    // level-l local interpolant of f - s_{<l}, evaluated on a level-l cell of
    // Omega^l \ Omega^{l+1}, where s_{<l} is the spline built from the
    // coefficients of the coarser levels. On such a cell all HB functions of
    // level > l vanish, so every spline of the HB space is reproduced.
    // HB functions do not sum to one, so reproduction of the constant 1 is
    // measured on function values, not on coefficients.

    // d x 2 refinement box, column 0 = lower and column 1 = upper corner.
    gsMatrix<real_t> hierBox(index_t d, real_t lo, real_t hi, real_t hiFirst = -1.0)
    {
        gsMatrix<real_t> box(d, 2);
        box.col(0).setConstant(lo);
        box.col(1).setConstant(hi);
        if (hiFirst >= 0.0)
            box(0,1) = hiFirst;
        return box;
    }

    // Hierarchical basis over a tensor B-spline basis of degree p with 4
    // elements per direction, refined successively on the given boxes.
    template<short_t d, bool Trunc>
    gsTHBSplineBasis<d,real_t,Trunc> makeHier(int p, const std::vector<gsMatrix<real_t> > & boxes)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, p+1);
        std::vector<gsKnotVector<real_t> > kvs(d, kv);
        gsTensorBSplineBasis<d,real_t> tb(kvs);
        gsTHBSplineBasis<d,real_t,Trunc> b(tb);
        for (size_t k = 0; k != boxes.size(); ++k)
            b.refine(boxes[k]);
        return b;
    }

    // Number of functions of the unrefined tensor basis used by makeHier.
    template<short_t d>
    index_t tensorSize(int p)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, p+1);
        std::vector<gsKnotVector<real_t> > kvs(d, kv);
        return gsTensorBSplineBasis<d,real_t>(kvs).size();
    }

    // Maximum pointwise difference of two functions on a grid over support.
    real_t maxDiffOnGrid(const gsFunction<real_t> & a, const gsFunction<real_t> & b,
                         const gsMatrix<real_t> & support, index_t n)
    {
        gsGridIterator<real_t,CUBE> pt(support, n);
        real_t err = 0;
        for (; pt; ++pt)
        {
            const gsMatrix<real_t> d = a.eval(*pt) - b.eval(*pt);
            if (!d.allFinite())
                return std::numeric_limits<real_t>::infinity();
            err = (std::max)(err, d.cwiseAbs().maxCoeff());
        }
        return err;
    }

    // Quasi-interpolating a random spline of the basis returns its coefficients.
    void checkProjector(const gsBasis<real_t> & basis, index_t tbSize, uint64_t seed, real_t eps)
    {
        CHECK(basis.size() > tbSize);
        const gsMatrix<real_t> c = lcgMatrix(basis.size(), 1, -1.0, 1.0, seed);
        const gsMatrix<real_t> res = qiOfSpline(basis, c);
        CHECK_EQUAL(c.rows(), res.rows());
        CHECK_EQUAL(c.cols(), res.cols());
        if (res.rows() == c.rows() && res.cols() == c.cols())
        {
            CHECK(res.allFinite());
            CHECK_CLOSE(0.0, (res - c).cwiseAbs().maxCoeff(), eps);
        }
    }

    // Level-0 refinement patterns shared by the HB tests.
    std::vector<gsMatrix<real_t> > hbCorner()
    { return std::vector<gsMatrix<real_t> >(1, hierBox(2, 0.0, 0.5)); }

    std::vector<gsMatrix<real_t> > hbNested()
    {
        std::vector<gsMatrix<real_t> > b;
        b.push_back(hierBox(2, 0.0, 0.5));
        b.push_back(hierBox(2, 0.0, 0.25));
        return b;
    }

    // One level-0 element wide (the knot vector has 4 elements).
    std::vector<gsMatrix<real_t> > hbWestStrip()
    {
        gsMatrix<real_t> box = hierBox(2, 0.0, 1.0, 0.25);
        return std::vector<gsMatrix<real_t> >(1, box);
    }

    TEST(intpl_projector_HB_2D_corner)
    {
        for (int p = 2; p <= 3; ++p)
        {
            gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(p, hbCorner());
            checkProjector(hb, tensorSize<2>(p), 41 + p, hbTol);
        }
    }

    TEST(intpl_projector_HB_2D_nested)
    {
        for (int p = 2; p <= 3; ++p)
        {
            gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(p, hbNested());
            CHECK_EQUAL(2u, hb.maxLevel());
            checkProjector(hb, tensorSize<2>(p), 44 + p, hbTol);
        }
    }

    TEST(intpl_projector_HB_2D_westStrip)
    {
        for (int p = 2; p <= 3; ++p)
        {
            gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(p, hbWestStrip());
            checkProjector(hb, tensorSize<2>(p), 47 + p, hbTol);
        }
    }

    // The boundary basis of an HB basis is again HB.
    TEST(intpl_projector_HB_2D_boundary)
    {
        const int p = 3;
        gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(p, hbCorner());
        gsKnotVector<real_t> kv(0.0, 1.0, 3, p+1);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        const boundary::side sides[2] = { boundary::west, boundary::south };
        for (int s = 0; s != 2; ++s)
        {
            gsBasis<real_t>::uPtr bb = hb.boundaryBasis(sides[s]);
            const bool isHB = (dynamic_cast<const gsHBSplineBasis<1,real_t>*>(bb.get()) != nullptr);
            CHECK(isHB);
            checkProjector(*bb, tb.boundaryBasis(sides[s])->size(), 52 + s, hbTol);
        }
    }

    TEST(intpl_projector_HB_3D_boundary)
    {
        const int p = 2;
        std::vector<gsMatrix<real_t> > boxes(1, hierBox(3, 0.0, 0.5));
        gsHBSplineBasis<3,real_t> hb = makeHier<3,false>(p, boxes);
        gsKnotVector<real_t> kv(0.0, 1.0, 3, p+1);
        gsTensorBSplineBasis<3,real_t> tb(kv, kv, kv);
        const boundary::side sides[3] = { boundary::west, boundary::south, boundary::front };
        for (int s = 0; s != 3; ++s)
        {
            gsBasis<real_t>::uPtr bb = hb.boundaryBasis(sides[s]);
            const bool isHB = (dynamic_cast<const gsHBSplineBasis<2,real_t>*>(bb.get()) != nullptr);
            CHECK(isHB);
            checkProjector(*bb, tb.boundaryBasis(sides[s])->size(), 54 + s, hbTol);
        }
    }

    // The constant 1 lies in the HB space; its coefficients are not all 1.
    TEST(intpl_HB_constant)
    {
        gsFunctionExpr<real_t> one("1", 2);
        for (int k = 0; k != 2; ++k)
        {
            gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(3, k == 0 ? hbCorner() : hbNested());
            gsMatrix<real_t> coefs;
            gsQuasiInterpolate<real_t>::localIntpl(hb, one, coefs);
            gsGeometry<real_t>::uPtr geo = hb.makeGeometry(coefs);
            CHECK_CLOSE(0.0, maxDiffOnGrid(*geo, one, hb.support(), 21), hbTol);
        }
    }

    // For out-of-space data the HB scheme and the THB scheme on the same
    // hierarchy produce the same spline (truncation does not change the span).
    TEST(intpl_HB_matchesTHB)
    {
        gsFunctionExpr<real_t> f("sin(3*x)*cos(2*y) + exp(x*y)", 2);
        for (int p = 2; p <= 3; ++p)
        {
            gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(p, hbNested());
            gsTHBSplineBasis<2,real_t> thb = makeHier<2,true>(p, hbNested());
            CHECK_EQUAL(hb.size(), thb.size());

            gsMatrix<real_t> cHB, cTHB;
            gsQuasiInterpolate<real_t>::localIntpl(hb, f, cHB);
            gsQuasiInterpolate<real_t>::localIntpl(thb, f, cTHB);
            gsGeometry<real_t>::uPtr gHB  = hb.makeGeometry(cHB);
            gsGeometry<real_t>::uPtr gTHB = thb.makeGeometry(cTHB);
            CHECK_CLOSE(0.0, maxDiffOnGrid(*gHB, *gTHB, hb.support(), 21), hbTol);
        }
    }

    // A single functional of a level-l function needs the coarser-level
    // spline, so the per-index entry points are rejected on HB bases while
    // they stay valid on THB bases.
    TEST(intpl_HB_perIndexThrows)
    {
        gsFunctionExpr<real_t> f("sin(3*x)*cos(2*y) + exp(x*y)", 2);
        gsHBSplineBasis<2,real_t>  hb  = makeHier<2,false>(3, hbCorner());
        gsTHBSplineBasis<2,real_t> thb = makeHier<2,true>(3, hbCorner());
        CHECK_EQUAL(hb.size(), thb.size());

        const index_t idx[2] = { 0, hb.size() - 1 };
        CHECK_EQUAL(0, hb.levelOf(idx[0]));
        CHECK_EQUAL(1, hb.levelOf(idx[1]));
        for (int k = 0; k != 2; ++k)
        {
            const index_t i = idx[k];
            CHECK_THROW(gsQuasiInterpolate<real_t>::localIntpl(
                            static_cast<const gsBasis<real_t>&>(hb), f, i), std::exception);

            const gsMatrix<real_t> r = gsQuasiInterpolate<real_t>::localIntpl(
                            static_cast<const gsBasis<real_t>&>(thb), f, i);
            CHECK_EQUAL(1, r.size());
        }
    }

    // Rational basis over a non-truncated hierarchical source.
    class RationalHB2 : public gsRationalBasis<gsHBSplineBasis<2,real_t> >
    {
    public:
        typedef gsRationalBasis<gsHBSplineBasis<2,real_t> > Base;
        RationalHB2(gsHBSplineBasis<2,real_t> * src, gsMatrix<real_t> w) : Base(src, give(w)) { }
        GISMO_CLONE_FUNCTION(RationalHB2)
        gsGeometry<real_t>::uPtr makeGeometry(gsMatrix<real_t>) const override { GISMO_NO_IMPLEMENTATION }
        std::ostream & print(std::ostream & os) const override { return os << "RationalHB2\n"; }
    };

    TEST(intpl_rationalOverHB_bulkThrows)
    {
        gsFunctionExpr<real_t> f("sin(3*x)*cos(2*y) + exp(x*y)", 2);
        gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(3, hbCorner());
        const gsMatrix<real_t> w = lcgMatrix(hb.size(), 1, 0.5, 2.0, 56);
        CHECK(w.maxCoeff() - w.minCoeff() > 0.5);
        RationalHB2 rb(hb.clone().release(), w);
        gsMatrix<real_t> coefs;
        CHECK_THROW(gsQuasiInterpolate<real_t>::localIntpl(rb, f, coefs), std::exception);
    }

    TEST(l2_rationalOverHB_bulkThrows)
    {
        gsFunctionExpr<real_t> f("sin(3*x)*cos(2*y) + exp(x*y)", 2);
        gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(3, hbCorner());
        const gsMatrix<real_t> w = lcgMatrix(hb.size(), 1, 0.5, 2.0, 57);
        CHECK(w.maxCoeff() - w.minCoeff() > 0.5);
        RationalHB2 rb(hb.clone().release(), w);
        gsMatrix<real_t> coefs;
        CHECK_THROW(gsQuasiInterpolate<real_t>::localL2(rb, f, coefs), std::exception);
    }

    // Rational B-splines (NURBS) with non-constant weights: projector property.
    TEST(intpl_projector_NURBS_2D_full)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        const gsMatrix<real_t> w = lcgMatrix(tb.size(), 1, 0.5, 2.0, 57);
        CHECK(w.maxCoeff() - w.minCoeff() > 0.5);
        gsTensorNurbsBasis<2,real_t> nb(tb.clone().release(), w);

        const gsMatrix<real_t> c = lcgMatrix(nb.size(), 1, -1.0, 1.0, 58);
        const gsMatrix<real_t> res = qiOfSpline(nb, c);
        CHECK_EQUAL(c.rows(), res.rows());
        CHECK_EQUAL(c.cols(), res.cols());
        if (res.rows() == c.rows() && res.cols() == c.cols())
        {
            CHECK(res.allFinite());
            CHECK_CLOSE(0.0, (res - c).cwiseAbs().maxCoeff(), projTol);
        }
    }

    TEST(intpl_projector_NURBS_2D_boundary)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        const gsMatrix<real_t> w = lcgMatrix(tb.size(), 1, 0.5, 2.0, 59);
        CHECK(w.maxCoeff() - w.minCoeff() > 0.5);
        gsTensorNurbsBasis<2,real_t> nb(tb.clone().release(), w);
        gsBasis<real_t>::uPtr bb = nb.boundaryBasis(boundary::west);
        CHECK(dynamic_cast<const gsNurbsBasis<real_t>*>(bb.get()) != nullptr);

        const gsMatrix<real_t> c = lcgMatrix(bb->size(), 1, -1.0, 1.0, 60);
        const gsMatrix<real_t> res = qiOfSpline(*bb, c);
        CHECK_EQUAL(c.rows(), res.rows());
        CHECK_EQUAL(c.cols(), res.cols());
        if (res.rows() == c.rows() && res.cols() == c.cols())
        {
            CHECK(res.allFinite());
            CHECK_CLOSE(0.0, (res - c).cwiseAbs().maxCoeff(), projTol);
        }
    }

    // ============================== localL2 ================================

    // Local L2 projection reproduces polynomials up to the basis degree.
    TEST(l2_reproduction_2D)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 3); // deg 2
        gsTensorBSplineBasis<2,real_t> basis(kv, kv);

        gsFunctionExpr<real_t> f("1 + x - 2*y + 3*x*y + x^2 - y^2", 2);
        gsMatrix<real_t> coefs;
        gsQuasiInterpolate<real_t>::localL2(basis, f, coefs);
        CHECK_CLOSE(0.0, errVsFun(basis, f, coefs), tol);
    }

    TEST(l2_reproduction_THB)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 3); // deg 2
        gsTensorBSplineBasis<2,real_t> tbasis(kv, kv);
        gsTHBSplineBasis<2,real_t> basis(tbasis);
        gsMatrix<real_t> box(2,2);
        box.col(0) << 0.0, 0.0;
        box.col(1) << 0.5, 0.5;
        basis.refine(box);

        gsFunctionExpr<real_t> f("1 + x - 2*y + 3*x*y + x^2 - y^2", 2);
        gsMatrix<real_t> coefs;
        gsQuasiInterpolate<real_t>::localL2(basis, f, coefs);
        CHECK_CLOSE(0.0, errVsFun(basis, f, coefs), tol);
    }

    // localL2 of the spline with coefficients c (in basis b), in b itself.
    gsMatrix<real_t> l2OfSpline(const gsBasis<real_t> & b, const gsMatrix<real_t> & c)
    {
        gsGeometry<real_t>::uPtr f = b.makeGeometry(c);
        gsMatrix<real_t> res;
        gsQuasiInterpolate<real_t>::localL2(b, *f, res);
        return res;
    }

    // L2-projecting a random spline of the basis returns its coefficients.
    void checkL2Projector(const gsBasis<real_t> & basis, index_t tbSize, uint64_t seed, real_t eps)
    {
        CHECK(basis.size() > tbSize);
        const gsMatrix<real_t> c = lcgMatrix(basis.size(), 1, -1.0, 1.0, seed);
        const gsMatrix<real_t> res = l2OfSpline(basis, c);
        CHECK_EQUAL(c.rows(), res.rows());
        CHECK_EQUAL(c.cols(), res.cols());
        if (res.rows() == c.rows() && res.cols() == c.cols())
        {
            CHECK(res.allFinite());
            CHECK_CLOSE(0.0, (res - c).cwiseAbs().maxCoeff(), eps);
        }
    }

    // Runs call(); sets threw if it throws a std::exception and returns what was
    // written to std::cerr meanwhile (GISMO_ENSURE / GISMO_ERROR throw a literal
    // and write their message there).
    template<class Call>
    std::string cerrOfCall(Call call, bool & threw)
    {
        std::ostringstream err;
        std::streambuf * old = std::cerr.rdbuf(err.rdbuf());
        threw = false;
        try { call(); } catch (const std::exception &) { threw = true; }
        std::cerr.rdbuf(old);
        return err.str();
    }

    // Local L2 on an HB basis: the coefficient of a level-l function
    // is the Gram solve on one level-l cell of the tensor basis of level l,
    // applied to the residual f - s_{<l} of the coarser levels (same
    // level-by-level argument as for localIntpl above). HB and THB spaces on
    // the same hierarchy coincide, so both schemes project onto the same space.

    TEST(l2_projector_THB_2D_nested)
    {
        for (int p = 2; p <= 3; ++p)
        {
            gsTHBSplineBasis<2,real_t> thb = makeHier<2,true>(p, hbNested());
            checkL2Projector(thb, tensorSize<2>(p), 61 + p, l2Tol);
        }
    }

    TEST(l2_projector_HB_2D_corner)
    {
        for (int p = 2; p <= 3; ++p)
        {
            gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(p, hbCorner());
            checkL2Projector(hb, tensorSize<2>(p), 63 + p, l2Tol);
        }
    }

    TEST(l2_projector_HB_2D_nested)
    {
        for (int p = 2; p <= 3; ++p)
        {
            gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(p, hbNested());
            CHECK_EQUAL(2u, hb.maxLevel());
            checkL2Projector(hb, tensorSize<2>(p), 66 + p, l2Tol);
        }
    }

    TEST(l2_projector_HB_2D_westStrip)
    {
        for (int p = 2; p <= 3; ++p)
        {
            gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(p, hbWestStrip());
            checkL2Projector(hb, tensorSize<2>(p), 69 + p, l2Tol);
        }
    }

    // The boundary basis of an HB basis is again HB.
    TEST(l2_projector_HB_2D_boundary)
    {
        const int p = 3;
        gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(p, hbCorner());
        gsKnotVector<real_t> kv(0.0, 1.0, 3, p+1);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        const boundary::side sides[2] = { boundary::west, boundary::south };
        for (int s = 0; s != 2; ++s)
        {
            gsBasis<real_t>::uPtr bb = hb.boundaryBasis(sides[s]);
            const bool isHB = (dynamic_cast<const gsHBSplineBasis<1,real_t>*>(bb.get()) != nullptr);
            CHECK(isHB);
            checkL2Projector(*bb, tb.boundaryBasis(sides[s])->size(), 72 + s, l2Tol);
        }
    }

    TEST(l2_projector_HB_3D_boundary)
    {
        const int p = 2;
        std::vector<gsMatrix<real_t> > boxes(1, hierBox(3, 0.0, 0.5));
        gsHBSplineBasis<3,real_t> hb = makeHier<3,false>(p, boxes);
        gsKnotVector<real_t> kv(0.0, 1.0, 3, p+1);
        gsTensorBSplineBasis<3,real_t> tb(kv, kv, kv);
        const boundary::side sides[3] = { boundary::west, boundary::south, boundary::front };
        for (int s = 0; s != 3; ++s)
        {
            gsBasis<real_t>::uPtr bb = hb.boundaryBasis(sides[s]);
            const bool isHB = (dynamic_cast<const gsHBSplineBasis<2,real_t>*>(bb.get()) != nullptr);
            CHECK(isHB);
            checkL2Projector(*bb, tb.boundaryBasis(sides[s])->size(), 74 + s, l2Tol);
        }
    }

    // The constant 1 lies in the HB space; its coefficients are not all 1.
    TEST(l2_HB_constant)
    {
        gsFunctionExpr<real_t> one("1", 2);
        for (int k = 0; k != 2; ++k)
        {
            gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(3, k == 0 ? hbCorner() : hbNested());
            gsMatrix<real_t> coefs;
            gsQuasiInterpolate<real_t>::localL2(hb, one, coefs);
            gsGeometry<real_t>::uPtr geo = hb.makeGeometry(coefs);
            CHECK_CLOSE(0.0, maxDiffOnGrid(*geo, one, hb.support(), 21), l2Tol);
        }
    }

    // For out-of-space data the HB and THB schemes on the same hierarchy give
    // the same spline. Let lambda_i^l be the dual functional of the level-l
    // tensor basis supported on one level-l cell (lambda_i^l(B_k^l) = delta_ik).
    // Every coarser B_j^m restricted to that cell lies in the level-l tensor
    // space, so lambda_i^l(B_j^m) is exactly the coefficient that truncation
    // removes. Hence the HB residual coefficient lambda_i^l(f - s_{<l}) equals
    // the level-l HB coefficient of the THB scheme, and induction over the
    // levels gives the same spline for any f. This holds when both schemes use
    // the same functional: the same cell and the same Gram solve with the
    // (p+1)^d Gauss-Lobatto points.
    TEST(l2_HB_matchesTHB)
    {
        gsFunctionExpr<real_t> f("sin(3*x)*cos(2*y) + exp(x*y)", 2);
        for (int p = 2; p <= 3; ++p)
        {
            gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(p, hbNested());
            gsTHBSplineBasis<2,real_t> thb = makeHier<2,true>(p, hbNested());
            CHECK_EQUAL(hb.size(), thb.size());

            gsMatrix<real_t> cHB, cTHB;
            gsQuasiInterpolate<real_t>::localL2(hb, f, cHB);
            gsQuasiInterpolate<real_t>::localL2(thb, f, cTHB);
            gsGeometry<real_t>::uPtr gHB  = hb.makeGeometry(cHB);
            gsGeometry<real_t>::uPtr gTHB = thb.makeGeometry(cTHB);
            CHECK_CLOSE(0.0, maxDiffOnGrid(*gHB, *gTHB, hb.support(), 21), l2Tol);
        }
    }

    // A single functional of a level-l function needs the coarser-level
    // spline, so the per-index entry point is rejected on HB bases while it
    // stays valid on THB bases.
    TEST(l2_HB_perIndexThrows)
    {
        gsFunctionExpr<real_t> f("sin(3*x)*cos(2*y) + exp(x*y)", 2);
        gsHBSplineBasis<2,real_t>  hb  = makeHier<2,false>(3, hbCorner());
        gsTHBSplineBasis<2,real_t> thb = makeHier<2,true>(3, hbCorner());
        CHECK_EQUAL(hb.size(), thb.size());

        const index_t idx[2] = { 0, hb.size() - 1 };
        CHECK_EQUAL(0, hb.levelOf(idx[0]));
        CHECK_EQUAL(1, hb.levelOf(idx[1]));
        for (int k = 0; k != 2; ++k)
        {
            const index_t i = idx[k];
            const gsBasis<real_t> & hbRef = hb;
            bool threw = false;
            const std::string msg = cerrOfCall([&]() {
                gsQuasiInterpolate<real_t>::localL2(hbRef, f, i); }, threw);
            CHECK(threw);
            CHECK(msg.find("gsQuasiInterpolate::localL2") != std::string::npos);

            const gsMatrix<real_t> r = gsQuasiInterpolate<real_t>::localL2(
                            static_cast<const gsBasis<real_t>&>(thb), f, i);
            CHECK_EQUAL(1, r.size());
        }
    }

    // Rational THB bases with genuine weights: projector property.
    TEST(l2_projector_rationalTHB_2D_corner)
    {
        for (int p = 2; p <= 3; ++p)
        {
            gsTHBSplineBasis<2,real_t> thb = makeHier<2,true>(p, hbCorner());
            const gsMatrix<real_t> w = lcgMatrix(thb.size(), 1, 0.5, 2.0, 77 + p);
            CHECK(w.maxCoeff() - w.minCoeff() > 0.5);
            gsRationalTHBSplineBasis<2,real_t> rb(thb.clone().release(), w);
            checkL2Projector(rb, tensorSize<2>(p), 79 + p, l2Tol);
        }
    }

    TEST(l2_projector_rationalTHB_2D_nested)
    {
        for (int p = 2; p <= 3; ++p)
        {
            gsTHBSplineBasis<2,real_t> thb = makeHier<2,true>(p, hbNested());
            const gsMatrix<real_t> w = lcgMatrix(thb.size(), 1, 0.5, 2.0, 82 + p);
            CHECK(w.maxCoeff() - w.minCoeff() > 0.5);
            gsRationalTHBSplineBasis<2,real_t> rb(thb.clone().release(), w);
            checkL2Projector(rb, tensorSize<2>(p), 84 + p, l2Tol);
        }
    }

    TEST(l2_projector_rationalTHB_2D_boundary)
    {
        const int p = 3;
        gsTHBSplineBasis<2,real_t> thb = makeHier<2,true>(p, hbCorner());
        const gsMatrix<real_t> w = lcgMatrix(thb.size(), 1, 0.5, 2.0, 87);
        CHECK(w.maxCoeff() - w.minCoeff() > 0.5);
        gsRationalTHBSplineBasis<2,real_t> rb(thb.clone().release(), w);
        gsKnotVector<real_t> kv(0.0, 1.0, 3, p+1);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsBasis<real_t>::uPtr bb = rb.boundaryBasis(boundary::west);
        checkL2Projector(*bb, tb.boundaryBasis(boundary::west)->size(), 88, l2Tol);
    }

    // 3-D rational THB, p=3, non-constant weights. More rational THB functions
    // can be active on a level-l cell than the local polynomial dimension, so a
    // local L2 solve with the active rational functions is not well posed; the
    // scheme must instead project f*W with the level-l tensor basis and divide by
    // w_i. The 1e-6 tolerance sits above the 1e-8 round-off of the 3-D p=3 solves.
    TEST(l2_projector_rationalTHB_3D_p3)
    {
        const int p = 3;
        std::vector<gsMatrix<real_t> > boxes(1, hierBox(3, 0.0, 0.5));
        gsTHBSplineBasis<3,real_t> thb = makeHier<3,true>(p, boxes);
        const gsMatrix<real_t> w = lcgMatrix(thb.size(), 1, 0.5, 2.0, 211);
        CHECK(w.maxCoeff() - w.minCoeff() > 0.5);
        gsRationalTHBSplineBasis<3,real_t> rb(thb.clone().release(), w);
        checkL2Projector(rb, tensorSize<3>(p), 212, 1e-6);
    }

    // ============================= Schoenberg ==============================

    // Variation-diminishing spline approximation reproduces affine functions.
    TEST(schoenberg_reproduction_affine)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 3); // deg 2
        gsTensorBSplineBasis<2,real_t> basis(kv, kv);

        gsFunctionExpr<real_t> f("1 + 2*x - 3*y", 2);
        gsMatrix<real_t> coefs;
        gsQuasiInterpolate<real_t>::Schoenberg(basis, f, coefs);
        CHECK_CLOSE(0.0, errVsFun(basis, f, coefs), tol);
    }

    // ============================== EvalBased ==============================

    // Evaluation-based QI (1D). Special-case (deg 1/2/3) and general formulas
    // both reproduce polynomials up to the basis degree.
    TEST(evalbased_reproduction_1D_special)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 4, 3); // deg 2
        gsBSplineBasis<real_t> basis(kv);

        gsFunctionExpr<real_t> f("x^2 - 3*x + 1", 1);
        gsMatrix<real_t> coefs;
        gsQuasiInterpolate<real_t>::EvalBased(basis, f, true, coefs);
        CHECK_CLOSE(0.0, errVsFun(basis, f, coefs), tol);
    }

    TEST(evalbased_reproduction_1D_general)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 4, 3); // deg 2
        gsBSplineBasis<real_t> basis(kv);

        gsFunctionExpr<real_t> f("x^2 - 3*x + 1", 1);
        gsMatrix<real_t> coefs;
        gsQuasiInterpolate<real_t>::EvalBased(basis, f, false, coefs);
        CHECK_CLOSE(0.0, errVsFun(basis, f, coefs), tol);
    }

    // ============ localTaylor: supported and rejected bases ============

    // Taylor QI needs the derivatives of order d*p at the anchor of a single
    // tensor B-spline and is defined for tensor B-spline bases only. The
    // sources below are spline geometries because gsFunctionExpr provides
    // derivatives up to order 2 only.

    void checkTaylorProjector(const gsBasis<real_t> & basis, uint64_t seed)
    {
        const index_t r = basis.degree(0);
        const gsMatrix<real_t> c = lcgMatrix(basis.size(), 1, -1.0, 1.0, seed);
        gsGeometry<real_t>::uPtr g = basis.makeGeometry(c);
        gsMatrix<real_t> res;
        gsQuasiInterpolate<real_t>::localTaylor(basis, *g, r, res);
        CHECK_EQUAL(c.rows(), res.rows());
        CHECK_EQUAL(c.cols(), res.cols());
        if (res.rows() == c.rows() && res.cols() == c.cols())
        {
            CHECK(res.allFinite());
            CHECK_CLOSE(0.0, (res - c).cwiseAbs().maxCoeff(), projTol);
        }
    }

    TEST(taylor_projector_tensor_1D)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
        gsBSplineBasis<real_t> basis(kv);
        checkTaylorProjector(basis, 90);
    }

    TEST(taylor_projector_tensor_2D)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
        gsTensorBSplineBasis<2,real_t> basis(kv, kv);
        checkTaylorProjector(basis, 91);
    }

    TEST(taylor_projector_tensor_3D)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 4);
        gsTensorBSplineBasis<3,real_t> basis(kv, kv, kv);
        checkTaylorProjector(basis, 92);
    }

    TEST(taylor_reject_THB_bulk)
    {
        const int p = 2;
        gsTHBSplineBasis<2,real_t> thb = makeHier<2,true>(p, hbCorner());
        gsKnotVector<real_t> kv(0.0, 1.0, 3, p+1);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsGeometry<real_t>::uPtr g = tb.makeGeometry(lcgMatrix(tb.size(), 1, -1.0, 1.0, 93));
        const index_t r = p;
        gsMatrix<real_t> coefs;
        bool threw = false;
        const std::string msg = cerrOfCall([&]() {
            gsQuasiInterpolate<real_t>::localTaylor(thb, *g, r, coefs); }, threw);
        CHECK(threw);
        CHECK(msg.find("gsQuasiInterpolate::localTaylor") != std::string::npos);
        CHECK(msg.find("gsQuasiInterpolate::localIntpl") != std::string::npos);
    }

    TEST(taylor_reject_HB_bulk)
    {
        const int p = 2;
        gsHBSplineBasis<2,real_t> hb = makeHier<2,false>(p, hbCorner());
        gsKnotVector<real_t> kv(0.0, 1.0, 3, p+1);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsGeometry<real_t>::uPtr g = tb.makeGeometry(lcgMatrix(tb.size(), 1, -1.0, 1.0, 94));
        const index_t r = p;
        gsMatrix<real_t> coefs;
        bool threw = false;
        const std::string msg = cerrOfCall([&]() {
            gsQuasiInterpolate<real_t>::localTaylor(hb, *g, r, coefs); }, threw);
        CHECK(threw);
        CHECK(msg.find("gsQuasiInterpolate::localTaylor") != std::string::npos);
        CHECK(msg.find("gsQuasiInterpolate::localIntpl") != std::string::npos);
    }

    // A tensor B-spline basis keeps the per-index entry point; hierarchical
    // bases are rejected.
    TEST(taylor_reject_hier_perIndex)
    {
        const int p = 2;
        gsKnotVector<real_t> kv(0.0, 1.0, 3, p+1);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsGeometry<real_t>::uPtr g = tb.makeGeometry(lcgMatrix(tb.size(), 1, -1.0, 1.0, 95));
        const index_t r = p;

        gsTHBSplineBasis<2,real_t> thb = makeHier<2,true>(p, hbCorner());
        gsHBSplineBasis<2,real_t>  hb  = makeHier<2,false>(p, hbCorner());
        const gsBasis<real_t> * hier[2] = { &thb, &hb };
        for (int b = 0; b != 2; ++b)
        {
            const gsBasis<real_t> & h = *hier[b];
            const index_t idx[2] = { 0, h.size() - 1 };
            for (int k = 0; k != 2; ++k)
            {
                const index_t i = idx[k];
                bool threw = false;
                const std::string msg = cerrOfCall([&]() {
                    gsQuasiInterpolate<real_t>::localTaylor(h, *g, r, i); }, threw);
                CHECK(threw);
                CHECK(msg.find("gsQuasiInterpolate::localTaylor") != std::string::npos);
                CHECK(msg.find("gsQuasiInterpolate::localIntpl") != std::string::npos);
            }
        }

        const gsMatrix<real_t> x = gsQuasiInterpolate<real_t>::localTaylor(
                        static_cast<const gsBasis<real_t>&>(tb), *g, r, index_t(0));
        CHECK_EQUAL(1, x.size());
    }

    // Rational bases are not tensor B-spline bases: the bulk call throws.
    TEST(taylor_rational_NURBS_bulkThrows)
    {
        const int p = 2;
        gsKnotVector<real_t> kv(0.0, 1.0, 3, p+1);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsGeometry<real_t>::uPtr g = tb.makeGeometry(lcgMatrix(tb.size(), 1, -1.0, 1.0, 96));
        const gsMatrix<real_t> w = lcgMatrix(tb.size(), 1, 0.5, 2.0, 97);
        CHECK(w.maxCoeff() - w.minCoeff() > 0.5);
        gsTensorNurbsBasis<2,real_t> nb(tb.clone().release(), w);
        const index_t r = p;
        gsMatrix<real_t> coefs;
        CHECK_THROW(gsQuasiInterpolate<real_t>::localTaylor(nb, *g, r, coefs), std::exception);
    }

    TEST(taylor_rational_THB_bulkThrows)
    {
        const int p = 2;
        gsKnotVector<real_t> kv(0.0, 1.0, 3, p+1);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsGeometry<real_t>::uPtr g = tb.makeGeometry(lcgMatrix(tb.size(), 1, -1.0, 1.0, 98));
        gsTHBSplineBasis<2,real_t> thb = makeHier<2,true>(p, hbCorner());
        const gsMatrix<real_t> w = lcgMatrix(thb.size(), 1, 0.5, 2.0, 99);
        CHECK(w.maxCoeff() - w.minCoeff() > 0.5);
        gsRationalTHBSplineBasis<2,real_t> rb(thb.clone().release(), w);
        const index_t r = p;
        gsMatrix<real_t> coefs;
        CHECK_THROW(gsQuasiInterpolate<real_t>::localTaylor(rb, *g, r, coefs), std::exception);
    }

    // The 1-D wrapper rejects a basis that is not a univariate B-spline basis.
    TEST(taylor_whole1D_rejects2D)
    {
        const int p = 2;
        gsKnotVector<real_t> kv(0.0, 1.0, 3, p+1);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsGeometry<real_t>::uPtr g = tb.makeGeometry(lcgMatrix(tb.size(), 1, -1.0, 1.0, 100));
        const index_t r = p;
        gsMatrix<real_t> coefs;
        CHECK_THROW(gsQuasiInterpolate<real_t>::Taylor(tb, *g, r, coefs), std::exception);
    }
}
