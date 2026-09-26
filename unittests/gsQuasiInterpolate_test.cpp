/** @file gsQuasiInterpolate_test.cpp

    @brief Tests the quasi-interpolation schemes in gsQuasiInterpolate.

    Covered methods (each identified by a test-name prefix):
      * taylor_*      : localTaylor / Taylor / Taylor2D (Taylor QI). The
                        tensor-product Taylor QI Q_{p,p} reproduces every
                        spline in its tensor B-spline space (Lyche-Morken
                        Thm 8.5). Verified by exact polynomial reproduction
                        (order-<=2 regime, gsFunctionExpr) and by exact
                        reproduction of an arbitrary tensor B-spline geometry
                        in 2D/3D (needs mixed partials to order d*p, which
                        spline geometries provide but gsFunctionExpr does not).
      * intpl_*       : localIntpl  (local interpolation QI)
      * l2_*          : localL2     (local L2-projection QI)
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

    // Hierarchical (THB) dispatch path: gsHTensorBasis + tensorLevel(lvl) +
    // flatTensorIndexOf. Whole-domain refinement promotes every function to
    // level 1 (a pure tensor level), exercising the lvl>0 dispatch. Bilinear
    // source keeps derivatives analytically available.
    TEST(taylor_reproduction_THB_2D)
    {
        gsKnotVector<real_t> kv(0.0, 1.0, 3, 2); // deg 1
        gsTensorBSplineBasis<2,real_t> tbasis(kv, kv);
        gsTHBSplineBasis<2,real_t> basis(tbasis);
        gsMatrix<real_t> box(2,2);
        box.col(0) << 0.0, 0.0;
        box.col(1) << 1.0, 1.0;
        basis.refine(box);

        gsFunctionExpr<real_t> ffull("1 - 2*x + 3*y + 4*x*y", 2);
        gsMatrix<real_t> coefs;
        gsQuasiInterpolate<real_t>::localTaylor(basis, ffull, basis.degree(0), coefs);
        CHECK_CLOSE(0.0, errVsFun(basis, ffull, coefs), tol);
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
}
