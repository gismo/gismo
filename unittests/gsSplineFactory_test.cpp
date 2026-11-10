/** @file gsSplineFactory_test.cpp

    @brief Unit tests for spline factory methods (toBezier, squared, cubed, grad, lapl)
    compared against L2 projection, quasi-interpolation, and interpolateAtAnchors.

    For each derived quantity the factory result serves as the reference.
    All three alternative methods project / interpolate the same derived function
    onto the factory basis and should reproduce the factory coefficients to
    floating-point precision.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#include "gismo_unittest.h"

SUITE(gsSplineFactory_test)
{

// =============================================================================
// Function wrappers: pointwise evaluation of derived quantities from a 2D
// scalar B-spline c. Required because gsExprAssembler::getCoeff calls clone().
// =============================================================================

template <class T>
class SquaredSplineFunction final : public gsFunction<T>
{
    const gsTensorBSpline<2,T>& m_c;
public:
    explicit SquaredSplineFunction(const gsTensorBSpline<2,T>& c) : m_c(c) {}

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        m_c.eval_into(u, result);                         // 1 × npts
        result = result.cwiseProduct(result);
    }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    GISMO_CLONE_FUNCTION(SquaredSplineFunction)
};

template <class T>
class CubedSplineFunction final : public gsFunction<T>
{
    const gsTensorBSpline<2,T>& m_c;
public:
    explicit CubedSplineFunction(const gsTensorBSpline<2,T>& c) : m_c(c) {}

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        gsMatrix<T> val;
        m_c.eval_into(u, val);
        result = val.cwiseProduct(val.cwiseProduct(val));  // c³ element-wise
    }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    GISMO_CLONE_FUNCTION(CubedSplineFunction)
};

// ∂c/∂x_{dir} wrapper; dir ∈ {0,1}
template <class T>
class GradSplineFunction final : public gsFunction<T>
{
    const gsTensorBSpline<2,T>& m_c;
    short_t m_dir;
public:
    GradSplineFunction(const gsTensorBSpline<2,T>& c, short_t dir)
        : m_c(c), m_dir(dir) {}

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        gsMatrix<T> deriv;
        m_c.deriv_into(u, deriv);   // (domainDim*targetDim) × npts = 2 × npts
        result = deriv.row(m_dir);  // 1 × npts
    }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    GISMO_CLONE_FUNCTION(GradSplineFunction)
};

// Δc = ∂²c/∂x² + ∂²c/∂y² (second-order Voigt rows 0 and 1)
template <class T>
class LaplSplineFunction final : public gsFunction<T>
{
    const gsTensorBSpline<2,T>& m_c;
public:
    explicit LaplSplineFunction(const gsTensorBSpline<2,T>& c) : m_c(c) {}

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        gsMatrix<T> d2;
        m_c.deriv2_into(u, d2);     // 3 × npts for scalar 2D: [∂²/∂x², ∂²/∂y², ∂²/∂x∂y]
        result = d2.row(0) + d2.row(1);
    }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    GISMO_CLONE_FUNCTION(LaplSplineFunction)
};

// =============================================================================
// Helpers
// =============================================================================

static gsTensorBSpline<2,real_t> makeIdentityGeo()
{
    gsKnotVector<real_t> kv(0.0, 1.0, 0, 2); // [0 0 1 1]
    gsMatrix<real_t> coefs(4, 2);
    coefs << 0,0,  1,0,  0,1,  1,1;
    return gsTensorBSpline<2,real_t>(kv, kv, coefs);
}

template <class T>
void printRow(const std::string& qty, const std::string& method,
              const gsMatrix<T>& ref, const gsMatrix<T>& other, double secs)
{
#ifdef TEST_INFO
    gsMatrix<T> diff = ref - other;
    gsInfo << "  " << std::left
           << std::setw(10) << qty
           << std::setw(22) << method
           << std::scientific << std::setprecision(3)
           << std::setw(14) << diff.norm()
           << std::setw(14) << diff.template lpNorm<gsEigen::Infinity>()
           << secs << "\n";
#else
    (void)qty; (void)method; (void)ref; (void)other; (void)secs;
#endif
}

template <class T>
void printRef(const std::string& qty, double secs)
{
#ifdef TEST_INFO
    gsInfo << "  " << std::left
           << std::setw(10) << qty
           << std::setw(22) << "Factory (ref)"
           << std::setw(14) << "(exact)"
           << std::setw(14) << "(exact)"
           << secs << "\n";
#else
    (void)qty; (void)secs;
#endif
}

// =============================================================================
// Tests
// =============================================================================

// -----------------------------------------------------------------------------
// toBezier tests
// -----------------------------------------------------------------------------
TEST(toBezier_comparison_factory)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.toBezier(); (void)r; }
    double tFact = timer.stop() / numReps;

    auto bz = c.toBezier();

    gsMatrix<real_t> pts(2, 100);
    pts.setRandom();
    pts = (pts.array() * 0.5 + 0.5).matrix();
    gsMatrix<real_t> vc, vbz;
    c.eval_into(pts, vc);
    bz.eval_into(pts, vbz);
    real_t funcErr = (vc - vbz).norm();

#ifdef TEST_INFO
    gsInfo << std::left
           << "  " << std::setw(10) << "toBezier"
           << std::setw(22) << "Factory (ref)"
           << std::setw(14) << "(func check)"
           << std::scientific << std::setprecision(3)
           << std::setw(14) << funcErr
           << tFact << "  (func err at 100 pts)\n";
#endif

    CHECK_CLOSE(0.0, funcErr, 1e-10);
}

// -----------------------------------------------------------------------------
// squared tests
// -----------------------------------------------------------------------------
TEST(squared_comparison_factory)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.squared(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto c2 = c.squared();
    printRef<real_t>("squared", tFact);

    const gsTensorBSplineBasis<2,real_t>& bSq = c2.basis();
    const gsMatrix<real_t>& refCoefs = c2.coefs();

    CHECK_EQUAL(bSq.size(), refCoefs.size());
}

TEST(squared_comparison_L2)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.squared(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto c2 = c.squared();
    printRef<real_t>("squared", tFact);

    const gsTensorBSplineBasis<2,real_t>& bSq = c2.basis();
    const gsMatrix<real_t>& refCoefs = c2.coefs();

    SquaredSplineFunction<real_t> sq(c);

    gsMatrix<real_t> coefs_l2;
    timer.restart();
    gsL2Projection<real_t>::project(bSq, idGeo, sq, coefs_l2);
    double tL2 = timer.stop();
    printRow<real_t>("squared", "L2 projection", refCoefs, coefs_l2, tL2);
    CHECK_ARRAY_CLOSE(refCoefs.data(), coefs_l2.data(), refCoefs.size(), 1e-10);
}

TEST(squared_comparison_quasi)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.squared(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto c2 = c.squared();
    printRef<real_t>("squared", tFact);

    const gsTensorBSplineBasis<2,real_t>& bSq = c2.basis();
    const gsMatrix<real_t>& refCoefs = c2.coefs();

    SquaredSplineFunction<real_t> sq(c);

    gsMatrix<real_t> coefs_qi;
    timer.restart();
    gsQuasiInterpolate<real_t>::localIntpl(bSq, sq, coefs_qi);
    double tQI = timer.stop();
    printRow<real_t>("squared", "Quasi-interp", refCoefs, coefs_qi, tQI);
    CHECK_ARRAY_CLOSE(refCoefs.data(), coefs_qi.data(), refCoefs.size(), 1e-10);
}

TEST(squared_comparison_interpolate)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.squared(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto c2 = c.squared();
    printRef<real_t>("squared", tFact);

    const gsTensorBSplineBasis<2,real_t>& bSq = c2.basis();
    const gsMatrix<real_t>& refCoefs = c2.coefs();

    SquaredSplineFunction<real_t> sq(c);

    timer.restart();
    gsMatrix<real_t> anch = bSq.anchors();
    gsMatrix<real_t> vals;
    sq.eval_into(anch, vals);
    auto geom_ia = bSq.interpolateAtAnchors(vals);
    double tIA = timer.stop();
    printRow<real_t>("squared", "InterpAtAnchors", refCoefs, geom_ia->coefs(), tIA);
    CHECK_ARRAY_CLOSE(refCoefs.data(), geom_ia->coefs().data(), refCoefs.size(), 1e-10);
}

// -----------------------------------------------------------------------------
// cubed tests
// -----------------------------------------------------------------------------
TEST(cubed_comparison_factory)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.cubed(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto c3 = c.cubed();
    printRef<real_t>("cubed", tFact);

    const gsTensorBSplineBasis<2,real_t>& bCu = c3.basis();
    const gsMatrix<real_t>& refCoefs = c3.coefs();

    CHECK_EQUAL(bCu.size(), refCoefs.size());
}

TEST(cubed_comparison_quasi)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.cubed(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto c3 = c.cubed();
    printRef<real_t>("cubed", tFact);

    const gsTensorBSplineBasis<2,real_t>& bCu = c3.basis();
    const gsMatrix<real_t>& refCoefs = c3.coefs();

    CubedSplineFunction<real_t> cu(c);

    gsMatrix<real_t> coefs_qi;
    timer.restart();
    gsQuasiInterpolate<real_t>::localIntpl(bCu, cu, coefs_qi);
    printRow<real_t>("cubed", "Quasi-interp", refCoefs, coefs_qi, timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), coefs_qi.data(), refCoefs.size(), 1e-10);
}

TEST(cubed_comparison_interpolate)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.cubed(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto c3 = c.cubed();
    printRef<real_t>("cubed", tFact);

    const gsTensorBSplineBasis<2,real_t>& bCu = c3.basis();
    const gsMatrix<real_t>& refCoefs = c3.coefs();

    CubedSplineFunction<real_t> cu(c);

    timer.restart();
    gsMatrix<real_t> anch = bCu.anchors();
    gsMatrix<real_t> vals;
    cu.eval_into(anch, vals);
    auto geom_ia = bCu.interpolateAtAnchors(vals);
    printRow<real_t>("cubed", "InterpAtAnchors", refCoefs, geom_ia->coefs(), timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), geom_ia->coefs().data(), refCoefs.size(), 1e-10);
}

// -----------------------------------------------------------------------------
// grad_x tests
// -----------------------------------------------------------------------------
TEST(grad_x_comparison_factory)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.grad(0); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cdx = c.grad(0);
    printRef<real_t>("grad_x", tFact);

    const gsTensorBSplineBasis<2,real_t>& bGx = cdx.basis();
    const gsMatrix<real_t>& refCoefs = cdx.coefs();

    CHECK_EQUAL(bGx.size(), refCoefs.size());
}

TEST(grad_x_comparison_L2)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.grad(0); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cdx = c.grad(0);
    printRef<real_t>("grad_x", tFact);

    const gsTensorBSplineBasis<2,real_t>& bGx = cdx.basis();
    const gsMatrix<real_t>& refCoefs = cdx.coefs();

    GradSplineFunction<real_t> gx(c, 0);

    gsMatrix<real_t> coefs_l2;
    timer.restart();
    gsL2Projection<real_t>::project(bGx, idGeo, gx, coefs_l2);
    printRow<real_t>("grad_x", "L2 projection", refCoefs, coefs_l2, timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), coefs_l2.data(), refCoefs.size(), 1e-10);
}

TEST(grad_x_comparison_quasi)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.grad(0); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cdx = c.grad(0);
    printRef<real_t>("grad_x", tFact);

    const gsTensorBSplineBasis<2,real_t>& bGx = cdx.basis();
    const gsMatrix<real_t>& refCoefs = cdx.coefs();

    GradSplineFunction<real_t> gx(c, 0);

    gsMatrix<real_t> coefs_qi;
    timer.restart();
    gsQuasiInterpolate<real_t>::localIntpl(bGx, gx, coefs_qi);
    printRow<real_t>("grad_x", "Quasi-interp", refCoefs, coefs_qi, timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), coefs_qi.data(), refCoefs.size(), 1e-10);
}

TEST(grad_x_comparison_interpolate)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.grad(0); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cdx = c.grad(0);
    printRef<real_t>("grad_x", tFact);

    const gsTensorBSplineBasis<2,real_t>& bGx = cdx.basis();
    const gsMatrix<real_t>& refCoefs = cdx.coefs();

    GradSplineFunction<real_t> gx(c, 0);

    timer.restart();
    gsMatrix<real_t> anch = bGx.anchors();
    gsMatrix<real_t> vals;
    gx.eval_into(anch, vals);
    auto geom_ia = bGx.interpolateAtAnchors(vals);
    printRow<real_t>("grad_x", "InterpAtAnchors", refCoefs, geom_ia->coefs(), timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), geom_ia->coefs().data(), refCoefs.size(), 1e-10);
}

// -----------------------------------------------------------------------------
// grad_y tests
// -----------------------------------------------------------------------------
TEST(grad_y_comparison_factory)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.grad(1); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cdy = c.grad(1);
    printRef<real_t>("grad_y", tFact);

    const gsTensorBSplineBasis<2,real_t>& bGy = cdy.basis();
    const gsMatrix<real_t>& refCoefs = cdy.coefs();

    CHECK_EQUAL(bGy.size(), refCoefs.size());
}

TEST(grad_y_comparison_L2)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.grad(1); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cdy = c.grad(1);
    printRef<real_t>("grad_y", tFact);

    const gsTensorBSplineBasis<2,real_t>& bGy = cdy.basis();
    const gsMatrix<real_t>& refCoefs = cdy.coefs();

    GradSplineFunction<real_t> gy(c, 1);

    gsMatrix<real_t> coefs_l2;
    timer.restart();
    gsL2Projection<real_t>::project(bGy, idGeo, gy, coefs_l2);
    printRow<real_t>("grad_y", "L2 projection", refCoefs, coefs_l2, timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), coefs_l2.data(), refCoefs.size(), 1e-10);
}

TEST(grad_y_comparison_quasi)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.grad(1); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cdy = c.grad(1);
    printRef<real_t>("grad_y", tFact);

    const gsTensorBSplineBasis<2,real_t>& bGy = cdy.basis();
    const gsMatrix<real_t>& refCoefs = cdy.coefs();

    GradSplineFunction<real_t> gy(c, 1);

    gsMatrix<real_t> coefs_qi;
    timer.restart();
    gsQuasiInterpolate<real_t>::localIntpl(bGy, gy, coefs_qi);
    printRow<real_t>("grad_y", "Quasi-interp", refCoefs, coefs_qi, timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), coefs_qi.data(), refCoefs.size(), 1e-10);
}

TEST(grad_y_comparison_interpolate)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.grad(1); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cdy = c.grad(1);
    printRef<real_t>("grad_y", tFact);

    const gsTensorBSplineBasis<2,real_t>& bGy = cdy.basis();
    const gsMatrix<real_t>& refCoefs = cdy.coefs();

    GradSplineFunction<real_t> gy(c, 1);

    timer.restart();
    gsMatrix<real_t> anch = bGy.anchors();
    gsMatrix<real_t> vals;
    gy.eval_into(anch, vals);
    auto geom_ia = bGy.interpolateAtAnchors(vals);
    printRow<real_t>("grad_y", "InterpAtAnchors", refCoefs, geom_ia->coefs(), timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), geom_ia->coefs().data(), refCoefs.size(), 1e-10);
}

// -----------------------------------------------------------------------------
// lapl tests (requires GISMO_WITH_ADIFF for accuracy)
// -----------------------------------------------------------------------------
#ifdef GISMO_WITH_ADIFF
TEST(lapl_comparison_degree2_factory)
{
    const index_t degree   = 2;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.lapl(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cl = c.lapl();
    printRef<real_t>("lapl", tFact);

    const gsTensorBSplineBasis<2,real_t>& bL = cl.basis();
    const gsMatrix<real_t>& refCoefs = cl.coefs();

    CHECK_EQUAL(bL.size(), refCoefs.size());
}

TEST(lapl_comparison_degree2_L2)
{
    const index_t degree   = 2;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.lapl(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cl = c.lapl();
    printRef<real_t>("lapl", tFact);

    const gsTensorBSplineBasis<2,real_t>& bL = cl.basis();
    const gsMatrix<real_t>& refCoefs = cl.coefs();

    LaplSplineFunction<real_t> lp(c);

    gsMatrix<real_t> coefs_l2;
    timer.restart();
    gsL2Projection<real_t>::project(bL, idGeo, lp, coefs_l2);
    printRow<real_t>("lapl", "L2 projection", refCoefs, coefs_l2, timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), coefs_l2.data(), refCoefs.size(), 1e-10);
}

TEST(lapl_comparison_degree2_quasi)
{
    const index_t degree   = 2;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.lapl(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cl = c.lapl();
    printRef<real_t>("lapl", tFact);

    const gsTensorBSplineBasis<2,real_t>& bL = cl.basis();
    const gsMatrix<real_t>& refCoefs = cl.coefs();

    LaplSplineFunction<real_t> lp(c);

    gsMatrix<real_t> coefs_qi;
    timer.restart();
    gsQuasiInterpolate<real_t>::localIntpl(bL, lp, coefs_qi);
    printRow<real_t>("lapl", "Quasi-interp", refCoefs, coefs_qi, timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), coefs_qi.data(), refCoefs.size(), 1e-10);
}

TEST(lapl_comparison_degree2_interpolate)
{
    const index_t degree   = 2;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.lapl(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cl = c.lapl();
    printRef<real_t>("lapl", tFact);

    const gsTensorBSplineBasis<2,real_t>& bL = cl.basis();
    const gsMatrix<real_t>& refCoefs = cl.coefs();

    LaplSplineFunction<real_t> lp(c);

    timer.restart();
    gsMatrix<real_t> anch = bL.anchors();
    gsMatrix<real_t> vals;
    lp.eval_into(anch, vals);
    auto geom_ia = bL.interpolateAtAnchors(vals);
    printRow<real_t>("lapl", "InterpAtAnchors", refCoefs, geom_ia->coefs(), timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), geom_ia->coefs().data(), refCoefs.size(), 1e-10);
}

TEST(lapl_comparison_degree3_factory)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.lapl(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cl = c.lapl();
    printRef<real_t>("lapl", tFact);

    const gsTensorBSplineBasis<2,real_t>& bL = cl.basis();
    const gsMatrix<real_t>& refCoefs = cl.coefs();

    CHECK_EQUAL(bL.size(), refCoefs.size());
}

TEST(lapl_comparison_degree3_L2)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.lapl(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cl = c.lapl();
    printRef<real_t>("lapl", tFact);

    const gsTensorBSplineBasis<2,real_t>& bL = cl.basis();
    const gsMatrix<real_t>& refCoefs = cl.coefs();

    LaplSplineFunction<real_t> lp(c);

    gsMatrix<real_t> coefs_l2;
    timer.restart();
    gsL2Projection<real_t>::project(bL, idGeo, lp, coefs_l2);
    printRow<real_t>("lapl", "L2 projection", refCoefs, coefs_l2, timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), coefs_l2.data(), refCoefs.size(), 1e-10);
}

TEST(lapl_comparison_degree3_quasi)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.lapl(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cl = c.lapl();
    printRef<real_t>("lapl", tFact);

    const gsTensorBSplineBasis<2,real_t>& bL = cl.basis();
    const gsMatrix<real_t>& refCoefs = cl.coefs();

    LaplSplineFunction<real_t> lp(c);

    gsMatrix<real_t> coefs_qi;
    timer.restart();
    gsQuasiInterpolate<real_t>::localIntpl(bL, lp, coefs_qi);
    printRow<real_t>("lapl", "Quasi-interp", refCoefs, coefs_qi, timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), coefs_qi.data(), refCoefs.size(), 1e-10);
}

TEST(lapl_comparison_degree3_interpolate)
{
    const index_t degree   = 3;
    const index_t nRefine  = 2;
    const index_t numReps  = 5;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> anchors = basis.anchors();
    gsMatrix<real_t> fvals;
    f.eval_into(anchors, fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c = dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;

    timer.restart();
    for (index_t k = 0; k < numReps; ++k) { auto r = c.lapl(); (void)r; }
    double tFact = timer.stop() / numReps;
    auto cl = c.lapl();
    printRef<real_t>("lapl", tFact);

    const gsTensorBSplineBasis<2,real_t>& bL = cl.basis();
    const gsMatrix<real_t>& refCoefs = cl.coefs();

    LaplSplineFunction<real_t> lp(c);

    timer.restart();
    gsMatrix<real_t> anch = bL.anchors();
    gsMatrix<real_t> vals;
    lp.eval_into(anch, vals);
    auto geom_ia = bL.interpolateAtAnchors(vals);
    printRow<real_t>("lapl", "InterpAtAnchors", refCoefs, geom_ia->coefs(), timer.stop());
    CHECK_ARRAY_CLOSE(refCoefs.data(), geom_ia->coefs().data(), refCoefs.size(), 1e-10);
}
#endif

// =============================================================================
// Coverage for the general product, the common-space arithmetic and the
// Hessian.  These check VALUES pointwise rather than coefficient counts: the
// product of two splines is defined pointwise, so a size assertion cannot
// distinguish a correct convolution from a wrong one.
// =============================================================================

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
static gsTensorBSpline<2,real_t> makeSpline2D(index_t degree, index_t nRefine,
                                              const std::string & expr,
                                              index_t interior = 1)
{
    gsKnotVector<real_t> kv(0.0, 1.0, interior, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f(expr, 2);
    gsMatrix<real_t> fvals;
    f.eval_into(basis.anchors(), fvals);
    auto c = basis.interpolateAtAnchors(fvals);
    return gsTensorBSpline<2,real_t>(
        dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c));
}

static gsTensorBSpline<3,real_t> makeSpline3D(index_t degree, index_t nRefine,
                                              const std::string & expr,
                                              index_t interior = 1)
{
    gsKnotVector<real_t> kv(0.0, 1.0, interior, degree+1);
    gsTensorBSplineBasis<3,real_t> basis(kv, kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f(expr, 3);
    gsMatrix<real_t> fvals;
    f.eval_into(basis.anchors(), fvals);
    auto c = basis.interpolateAtAnchors(fvals);
    return gsTensorBSpline<3,real_t>(
        dynamic_cast<const gsTensorBSpline<3,real_t>&>(*c));
}

static gsMatrix<real_t> unitCubePoints(short_t dim, index_t n)
{
    gsMatrix<real_t> pts(dim, n);
    pts.setRandom();
    return (pts.array() * 0.5 + 0.5).matrix();
}

/// Max-norm of prod - A*B over \a pts.
template <short_t d>
static real_t productError(const gsTensorBSpline<d,real_t> & prod,
                           const gsTensorBSpline<d,real_t> & A,
                           const gsTensorBSpline<d,real_t> & B,
                           const gsMatrix<real_t> & pts)
{
    gsMatrix<real_t> vp, va, vb;
    prod.eval_into(pts, vp);
    A.eval_into(pts, va);
    B.eval_into(pts, vb);
    return (vp - va.cwiseProduct(vb)).template lpNorm<gsEigen::Infinity>();
}

/// Max-norm of pw - f(pts), where f is given as an expression.
template <short_t d>
static real_t exprError(const gsTensorBSpline<d,real_t> & spl, index_t col,
                        const std::string & expr, const gsMatrix<real_t> & pts)
{
    gsFunctionExpr<real_t> f(expr, d);
    gsMatrix<real_t> vs, vf;
    spl.eval_into(pts, vs);
    f.eval_into(pts, vf);
    return (vs.row(col) - vf.row(0)).template lpNorm<gsEigen::Infinity>();
}

// -----------------------------------------------------------------------------
// multiply
// -----------------------------------------------------------------------------
TEST(multiply_mixedDegrees_pointwise)
{
    // Different degree AND different knot vectors: the union of breakpoints
    // and the degree lift both have to be handled.
    gsTensorBSpline<2,real_t> A = makeSpline2D(3, 2, "sin(pi*x)*sin(pi*y)");
    gsTensorBSpline<2,real_t> B = makeSpline2D(2, 1, "1 + x*y + 0.3*y^2");

    gsTensorBSpline<2,real_t> P = gsTensorBSpline<2,real_t>::multiply(A, B);

    CHECK_EQUAL(5, P.degree(0));
    CHECK_EQUAL(5, P.degree(1));

    gsMatrix<real_t> pts = unitCubePoints(2, 200);
    CHECK_CLOSE(0.0, productError<2>(P, A, B, pts), 1e-12);
}

TEST(multiply_compresses_below_bezier)
{
    // The default path must actually leave Bézier form; if the knot-removal
    // tolerance silently rejects every removal, this is the check that fails.
    gsTensorBSpline<2,real_t> A = makeSpline2D(3, 2, "sin(pi*x)*sin(pi*y)");
    gsTensorBSpline<2,real_t> B = makeSpline2D(2, 1, "1 + x*y + 0.3*y^2");

    gsTensorBSpline<2,real_t> P   = gsTensorBSpline<2,real_t>::multiply(A, B, false);
    gsTensorBSpline<2,real_t> Pbz = gsTensorBSpline<2,real_t>::multiply(A, B, true);

    CHECK(P.coefsSize() < Pbz.coefsSize());

    // Minimal space: A is C^2 and B is C^1 at their common breakpoints, so the
    // product is C^1, i.e. interior multiplicity 5-1 = 4 at a knot of B.
    CHECK_EQUAL(4, (int)P.knots(0).multiplicity(0.5));

    gsMatrix<real_t> pts = unitCubePoints(2, 200);
    CHECK_CLOSE(0.0, productError<2>(P,   A, B, pts), 1e-12);
    CHECK_CLOSE(0.0, productError<2>(Pbz, A, B, pts), 1e-12);
}

TEST(multiply_matches_squared)
{
    // squared() must be the same kernel, not a second implementation.
    gsTensorBSpline<2,real_t> c = makeSpline2D(3, 2, "sin(pi*x)*sin(pi*y)");

    gsTensorBSpline<2,real_t> viaMultiply = gsTensorBSpline<2,real_t>::multiply(c, c);
    gsTensorBSpline<2,real_t> viaSquared  = c.squared();

    CHECK_EQUAL(viaSquared.coefsSize(), viaMultiply.coefsSize());
    CHECK_ARRAY_CLOSE(viaSquared.coefs().data(), viaMultiply.coefs().data(),
                      viaSquared.coefs().size(), 1e-14);
}

TEST(multiply_matches_cubed)
{
    gsTensorBSpline<2,real_t> c = makeSpline2D(3, 2, "sin(pi*x)*sin(pi*y)");

    gsTensorBSpline<2,real_t> viaMultiply =
        gsTensorBSpline<2,real_t>::multiply(c.squared(), c);
    gsTensorBSpline<2,real_t> viaCubed = c.cubed();

    CHECK_EQUAL(viaCubed.coefsSize(), viaMultiply.coefsSize());
    CHECK_ARRAY_CLOSE(viaCubed.coefs().data(), viaMultiply.coefs().data(),
                      viaCubed.coefs().size(), 1e-14);
}

// -----------------------------------------------------------------------------
// Low degrees — the Bernstein convolution has no p>=2 restriction
// -----------------------------------------------------------------------------
TEST(product_degree1)
{
    gsTensorBSpline<2,real_t> A = makeSpline2D(1, 2, "x + 2*y");
    gsTensorBSpline<2,real_t> B = makeSpline2D(1, 1, "3 - x");

    gsTensorBSpline<2,real_t> P  = gsTensorBSpline<2,real_t>::multiply(A, B);
    gsTensorBSpline<2,real_t> A2 = A.squared();
    gsTensorBSpline<2,real_t> A3 = A.cubed();

    CHECK_EQUAL(2, P.degree(0));
    CHECK_EQUAL(2, A2.degree(0));
    CHECK_EQUAL(3, A3.degree(0));

    gsMatrix<real_t> pts = unitCubePoints(2, 200);
    CHECK_CLOSE(0.0, productError<2>(P,  A, B,  pts), 1e-12);
    CHECK_CLOSE(0.0, productError<2>(A2, A, A,  pts), 1e-12);
    CHECK_CLOSE(0.0, productError<2>(A3, A2, A, pts), 1e-12);
}

TEST(product_degree2)
{
    gsTensorBSpline<2,real_t> A = makeSpline2D(2, 2, "cos(pi*x)*(1+y)");
    gsTensorBSpline<2,real_t> B = makeSpline2D(2, 2, "1 + x*y");

    gsTensorBSpline<2,real_t> P  = gsTensorBSpline<2,real_t>::multiply(A, B);
    gsTensorBSpline<2,real_t> A2 = A.squared();
    gsTensorBSpline<2,real_t> A3 = A.cubed();

    CHECK_EQUAL(4, P.degree(0));
    CHECK_EQUAL(6, A3.degree(1));

    gsMatrix<real_t> pts = unitCubePoints(2, 200);
    CHECK_CLOSE(0.0, productError<2>(P,  A, B,  pts), 1e-12);
    CHECK_CLOSE(0.0, productError<2>(A2, A, A,  pts), 1e-12);
    CHECK_CLOSE(0.0, productError<2>(A3, A2, A, pts), 1e-12);
}

// -----------------------------------------------------------------------------
// Repeated interior knots
// -----------------------------------------------------------------------------
static gsTensorBSpline<2,real_t> makeSplineWithDoubleKnot()
{
    // degree 3, interior knots 0.25, 0.5 (twice), 0.75  →  C^1 at 0.5
    real_t kn[] = {0,0,0,0, 0.25, 0.5, 0.5, 0.75, 1,1,1,1};
    std::vector<real_t> k(kn, kn + 12);
    gsKnotVector<real_t> kv(k, 3);

    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    gsFunctionExpr<real_t> f("sin(pi*x)*cos(pi*y) + x*y", 2);
    gsMatrix<real_t> fvals;
    f.eval_into(basis.anchors(), fvals);
    auto c = basis.interpolateAtAnchors(fvals);
    return gsTensorBSpline<2,real_t>(
        dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c));
}

TEST(toBezier_repeatedInteriorKnot)
{
    gsTensorBSpline<2,real_t> c = makeSplineWithDoubleKnot();
    CHECK_EQUAL(2, (int)c.knots(0).multiplicity(0.5));

    gsTensorBSpline<2,real_t> bz = c.toBezier();

    // Normalisation: every interior knot at multiplicity p+1, including the
    // one that already had multiplicity 2 — a "insert p+1 times" loop that
    // ignores the existing multiplicity would overshoot here.
    for (short_t k = 0; k != 2; ++k)
    {
        CHECK_EQUAL(4, (int)bz.knots(k).multiplicity(0.25));
        CHECK_EQUAL(4, (int)bz.knots(k).multiplicity(0.5));
        CHECK_EQUAL(4, (int)bz.knots(k).multiplicity(0.75));
    }

    gsMatrix<real_t> pts = unitCubePoints(2, 200);
    gsMatrix<real_t> vc, vb;
    c.eval_into(pts, vc);
    bz.eval_into(pts, vb);
    CHECK_CLOSE(0.0, (vc - vb).lpNorm<gsEigen::Infinity>(), 1e-12);
}

TEST(squared_repeatedInteriorKnot)
{
    gsTensorBSpline<2,real_t> c  = makeSplineWithDoubleKnot();
    gsTensorBSpline<2,real_t> c2 = c.squared();

    // c is C^2 at 0.25/0.75 and C^1 at 0.5, so c² inherits exactly that:
    // interior multiplicity 2p-cont = 6-2 = 4 and 6-1 = 5 respectively.
    // A hardcoded target of p+1 = 4 everywhere would leave 0.5 uncompressed.
    for (short_t k = 0; k != 2; ++k)
    {
        CHECK_EQUAL(4, (int)c2.knots(k).multiplicity(0.25));
        CHECK_EQUAL(5, (int)c2.knots(k).multiplicity(0.5));
        CHECK_EQUAL(4, (int)c2.knots(k).multiplicity(0.75));
    }

    gsMatrix<real_t> pts = unitCubePoints(2, 200);
    CHECK_CLOSE(0.0, productError<2>(c2, c, c, pts), 1e-12);
}

// -----------------------------------------------------------------------------
// 3D
// -----------------------------------------------------------------------------
TEST(multiply_3D_pointwise)
{
    gsTensorBSpline<3,real_t> A = makeSpline3D(2, 1, "sin(pi*x)*y*(1+z)");
    gsTensorBSpline<3,real_t> B = makeSpline3D(2, 1, "1 + x*z");

    gsTensorBSpline<3,real_t> P = gsTensorBSpline<3,real_t>::multiply(A, B);

    for (short_t k = 0; k != 3; ++k) CHECK_EQUAL(4, P.degree(k));

    gsMatrix<real_t> pts = unitCubePoints(3, 200);
    CHECK_CLOSE(0.0, productError<3>(P, A, B, pts), 1e-12);
}

TEST(squared_cubed_3D_pointwise)
{
    gsTensorBSpline<3,real_t> c = makeSpline3D(2, 1, "sin(pi*x)*y*(1+z)");

    gsTensorBSpline<3,real_t> c2 = c.squared();
    gsTensorBSpline<3,real_t> c3 = c.cubed();

    for (short_t k = 0; k != 3; ++k)
    {
        CHECK_EQUAL(4, c2.degree(k));
        CHECK_EQUAL(6, c3.degree(k));
    }

    gsMatrix<real_t> pts = unitCubePoints(3, 200);
    CHECK_CLOSE(0.0, productError<3>(c2, c,  c, pts), 1e-12);
    CHECK_CLOSE(0.0, productError<3>(c3, c2, c, pts), 1e-12);
}

TEST(grad_3D_pointwise)
{
    // An exact polynomial, so the derivative is known analytically and the
    // check does not lean on gismo's own derivative evaluation.
    gsTensorBSpline<3,real_t> c = makeSpline3D(2, 1, "x^2*y + z^2 - x*z");

    gsMatrix<real_t> pts = unitCubePoints(3, 200);
    CHECK_CLOSE(0.0, exprError<3>(c.grad(0), 0, "2*x*y - z",  pts), 1e-12);
    CHECK_CLOSE(0.0, exprError<3>(c.grad(1), 0, "x^2",        pts), 1e-12);
    CHECK_CLOSE(0.0, exprError<3>(c.grad(2), 0, "2*z - x",    pts), 1e-12);

    std::vector< gsTensorBSpline<3,real_t> > g = c.grad();
    CHECK_EQUAL(3u, (unsigned)g.size());
    for (short_t k = 0; k != 3; ++k)
    {
        gsMatrix<real_t> va, vb;
        g[k].eval_into(pts, va);
        c.grad(k).eval_into(pts, vb);
        CHECK_CLOSE(0.0, (va - vb).lpNorm<gsEigen::Infinity>(), 1e-14);
    }
}

// -----------------------------------------------------------------------------
// cubed vs L2 projection (2D) — the L2 counterpart the suite was missing
// -----------------------------------------------------------------------------
TEST(cubed_comparison_L2)
{
    // Lower degree than the neighbouring cubed tests on purpose: the cube of a
    // degree-3 spline lives in a degree-9 space, whose mass matrix is far too
    // ill-conditioned for the L2 solve to resolve the factory coefficients to
    // 1e-10 — the gate would then measure gsL2Projection, not cubed().  At
    // degree 2 the product space is degree 6, matching squared_comparison_L2.
    const index_t degree  = 2;
    const index_t nRefine = 1;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree+1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t i = 0; i < nRefine; ++i) basis.uniformRefine();

    gsFunctionExpr<real_t> f("sin(pi*x)*sin(pi*y)", 2);
    gsMatrix<real_t> fvals;
    f.eval_into(basis.anchors(), fvals);
    auto c_ptr = basis.interpolateAtAnchors(fvals);
    const gsTensorBSpline<2,real_t>& c =
        dynamic_cast<const gsTensorBSpline<2,real_t>&>(*c_ptr);

    gsTensorBSpline<2,real_t> idGeo = makeIdentityGeo();

    gsStopwatch timer;
    timer.restart();
    gsTensorBSpline<2,real_t> c3 = c.cubed();
    printRef<real_t>("cubed", timer.stop());

    const gsTensorBSplineBasis<2,real_t>& bCu = c3.basis();
    const gsMatrix<real_t>& refCoefs = c3.coefs();

    CubedSplineFunction<real_t> cu(c);

    gsMatrix<real_t> coefs_l2;
    timer.restart();
    gsL2Projection<real_t>::project(bCu, idGeo, cu, coefs_l2);
    double tL2 = timer.stop();
    printRow<real_t>("cubed", "L2 projection", refCoefs, coefs_l2, tL2);
    CHECK_ARRAY_CLOSE(refCoefs.data(), coefs_l2.data(), refCoefs.size(), 1e-10);
}

// -----------------------------------------------------------------------------
// makeCompatible / linearCombination
// -----------------------------------------------------------------------------
TEST(makeCompatible_commonSpace)
{
    std::vector< gsTensorBSpline<2,real_t> > v;
    v.push_back(makeSpline2D(3, 2, "sin(pi*x)*sin(pi*y)"));
    v.push_back(makeSpline2D(2, 1, "1 + x*y"));
    v.push_back(makeSpline2D(1, 3, "x - y"));

    std::vector< gsTensorBSpline<2,real_t> > before = v;
    gsTensorBSpline<2,real_t>::makeCompatible(v);

    // One space for all three...
    for (size_t i = 1; i != v.size(); ++i)
    {
        CHECK_EQUAL(v[0].coefsSize(), v[i].coefsSize());
        for (short_t k = 0; k != 2; ++k)
        {
            CHECK_EQUAL(v[0].degree(k), v[i].degree(k));
            CHECK_EQUAL((int)v[0].knots(k).size(), (int)v[i].knots(k).size());
        }
    }
    // ...and degree elevation plus knot insertion are both value-preserving.
    gsMatrix<real_t> pts = unitCubePoints(2, 200);
    for (size_t i = 0; i != v.size(); ++i)
    {
        gsMatrix<real_t> va, vb;
        v[i].eval_into(pts, va);
        before[i].eval_into(pts, vb);
        CHECK_CLOSE(0.0, (va - vb).lpNorm<gsEigen::Infinity>(), 1e-12);
    }
}

TEST(linearCombination_pointwise)
{
    gsTensorBSpline<2,real_t> A = makeSpline2D(3, 2, "sin(pi*x)*sin(pi*y)");
    gsTensorBSpline<2,real_t> B = makeSpline2D(2, 1, "1 + x*y + 0.3*y^2");

    gsTensorBSpline<2,real_t> R =
        gsTensorBSpline<2,real_t>::linearCombination(2.0, A, -3.0, B);

    gsMatrix<real_t> pts = unitCubePoints(2, 200);
    gsMatrix<real_t> vr, va, vb;
    R.eval_into(pts, vr);
    A.eval_into(pts, va);
    B.eval_into(pts, vb);
    CHECK_CLOSE(0.0, (vr - (2.0*va - 3.0*vb)).lpNorm<gsEigen::Infinity>(), 1e-12);
}

// The motivating case: det(J) of a planar map, built exactly as a spline from
// two products and one difference.
TEST(determinantOfJacobian_exact)
{
    typedef gsTensorBSpline<2,real_t> Spline;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, 4);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    basis.uniformRefine();

    gsFunctionExpr<real_t> fx("x + 0.2*sin(pi*x)*sin(pi*y)", 2);
    gsFunctionExpr<real_t> fy("y - 0.1*x*y", 2);
    gsMatrix<real_t> ax = basis.anchors(), vx, vy;
    fx.eval_into(ax, vx);
    fy.eval_into(ax, vy);

    gsMatrix<real_t> coefs(basis.size(), 2);
    coefs.col(0) = vx.transpose();
    coefs.col(1) = vy.transpose();
    Spline geo(basis, coefs);

    gsMatrix<real_t> c0 = coefs.col(0), c1 = coefs.col(1);
    Spline x(basis, c0), y(basis, c1);

    Spline det = Spline::linearCombination(
                     1.0, Spline::multiply(x.grad(0), y.grad(1)),
                    -1.0, Spline::multiply(x.grad(1), y.grad(0)));

    gsMatrix<real_t> pts = unitCubePoints(2, 100);
    gsMatrix<real_t> vdet, J;
    det.eval_into(pts, vdet);
    geo.jacobian_into(pts, J);

    real_t err = 0.0;
    for (index_t i = 0; i != pts.cols(); ++i)
        err = math::max(err,
              math::abs(vdet(0,i) - J.block(0, 2*i, 2, 2).determinant()));
    CHECK_CLOSE(0.0, err, 1e-12);
}

// -----------------------------------------------------------------------------
// hess
// -----------------------------------------------------------------------------
TEST(hess_2D_exactPolynomial)
{
    // c is exactly representable in the degree-2 space, so the Hessian is
    // known in closed form and the test does not depend on gismo's deriv2.
    gsTensorBSpline<2,real_t> c = makeSpline2D(2, 1, "x^2*y^2 + 3*x^2 - y^2 + x*y");

    gsTensorBSpline<2,real_t> H = c.hess();
    CHECK_EQUAL(4, H.targetDim());

    gsMatrix<real_t> pts = unitCubePoints(2, 200);
    CHECK_CLOSE(0.0, exprError<2>(H, 0, "2*y^2 + 6",  pts), 1e-11); // xx
    CHECK_CLOSE(0.0, exprError<2>(H, 1, "4*x*y + 1",  pts), 1e-11); // xy
    CHECK_CLOSE(0.0, exprError<2>(H, 2, "4*x*y + 1",  pts), 1e-11); // yx
    CHECK_CLOSE(0.0, exprError<2>(H, 3, "2*x^2 - 2",  pts), 1e-11); // yy
}

TEST(hess_symmetry_and_trace)
{
    gsTensorBSpline<3,real_t> c = makeSpline3D(3, 1, "sin(pi*x)*cos(pi*y)*(1+z^2)");

    gsTensorBSpline<3,real_t> H = c.hess();
    CHECK_EQUAL(9, H.targetDim());

    // Mixed partials commute, and hess() puts them all on one basis, so the
    // symmetric entries must agree in the coefficients, not just pointwise.
    for (short_t i = 0; i != 3; ++i)
        for (short_t j = i+1; j != 3; ++j)
            CHECK_CLOSE(0.0,
                (H.coefs().col(i*3+j) - H.coefs().col(j*3+i))
                    .lpNorm<gsEigen::Infinity>(), 1e-12);

    // trace(H) == lapl()
    gsTensorBSpline<3,real_t> L = c.lapl();
    gsMatrix<real_t> pts = unitCubePoints(3, 200);
    gsMatrix<real_t> vh, vl;
    H.eval_into(pts, vh);
    L.eval_into(pts, vl);
    gsMatrix<real_t> tr = vh.row(0) + vh.row(4) + vh.row(8);
    CHECK_CLOSE(0.0, (tr - vl).lpNorm<gsEigen::Infinity>(), 1e-10);
}

// What lapl(keepBezier) actually does, as opposed to what its doc claimed.
// The keepBezier branch inserts each interior knot exactly TWICE in the
// directions that were not differentiated, while degree elevation has already
// raised the differentiated direction by two.  The union of the two is then
// the same space in both cases whenever the input knots are simple, so the
// flag buys nothing and never reaches Bezier for p >= 3.
TEST(lapl_keepBezier_versus_minimal)
{
    gsTensorBSpline<2,real_t> c = makeSpline2D(3, 2, "sin(pi*x)*sin(pi*y)");

    gsTensorBSpline<2,real_t> Lmin = c.lapl(false);
    gsTensorBSpline<2,real_t> Lbez = c.lapl(true);

    const real_t xi = 0.5;
    const int pPlus1 = 4;                       // Bezier multiplicity at p = 3
    for (short_t k = 0; k != 2; ++k)
    {
        CHECK_EQUAL(3, (int)Lmin.knots(k).multiplicity(xi));   // m + 2, i.e. C^0
        CHECK_EQUAL(3, (int)Lbez.knots(k).multiplicity(xi));   // NOT pPlus1
        CHECK(Lbez.knots(k).multiplicity(xi) < pPlus1);
    }
    CHECK_EQUAL(Lmin.coefsSize(), Lbez.coefsSize());

    // Whatever the space, the values must agree with each other and with
    // the trace of the Hessian.
    gsMatrix<real_t> pts = unitCubePoints(2, 200);
    gsMatrix<real_t> vmin, vbez, vh;
    Lmin.eval_into(pts, vmin);
    Lbez.eval_into(pts, vbez);
    c.hess().eval_into(pts, vh);
    CHECK_CLOSE(0.0, (vmin - vbez).lpNorm<gsEigen::Infinity>(), 1e-10);
    CHECK_CLOSE(0.0, (vmin - (vh.row(0) + vh.row(3))).lpNorm<gsEigen::Infinity>(), 1e-10);
}

// The 3-D counterpart of determinantOfJacobian_exact: a 3x3 determinant is six
// TRIPLE products and a six-term sum, so it exercises chained multiply() and
// repeated linearCombination() in a way the 2x2 case does not.
TEST(determinantOfJacobian_exact_3D)
{
    typedef gsTensorBSpline<3,real_t> Spline;

    gsKnotVector<real_t> kv(0.0, 1.0, 1, 3);   // degree 2, one interior knot
    gsTensorBSplineBasis<3,real_t> basis(kv, kv, kv);

    const char * expr[3] = { "x + 0.15*x*y*z", "y - 0.10*x*z", "z + 0.05*x*y" };
    gsMatrix<real_t> ax = basis.anchors();
    gsMatrix<real_t> coefs(basis.size(), 3), v;
    for (short_t i = 0; i != 3; ++i)
    {
        gsFunctionExpr<real_t>(expr[i], 3).eval_into(ax, v);
        coefs.col(i) = v.transpose();
    }
    Spline geo(basis, coefs);

    // J(i,j) = d(component i) / d(x_j), each a scalar spline
    std::vector< std::vector<Spline> > J(3);
    for (short_t i = 0; i != 3; ++i)
    {
        gsMatrix<real_t> ci = coefs.col(i);
        Spline comp(basis, ci);
        for (short_t j = 0; j != 3; ++j)
            J[i].push_back(comp.grad(j));
    }

    // Leibniz expansion of det(J) over the six permutations of (0,1,2).
    const int perm[6][3] = {{0,1,2},{1,2,0},{2,0,1},{0,2,1},{2,1,0},{1,0,2}};
    const real_t sign[6] = {1,1,1,-1,-1,-1};

    Spline det;
    for (int t = 0; t != 6; ++t)
    {
        Spline term = Spline::multiply(
                          Spline::multiply(J[0][perm[t][0]], J[1][perm[t][1]]),
                          J[2][perm[t][2]]);
        det = (t == 0) ? term
                       : Spline::linearCombination(1.0, det, sign[t], term);
    }

    for (short_t k = 0; k != 3; ++k)
        CHECK_EQUAL(3*2 - 1, det.degree(k));   // three factors: (p-1) + p + p

    gsMatrix<real_t> pts = unitCubePoints(3, 100);
    gsMatrix<real_t> vdet, Jm;
    det.eval_into(pts, vdet);
    geo.jacobian_into(pts, Jm);

    real_t err = 0.0;
    for (index_t i = 0; i != pts.cols(); ++i)
        err = math::max(err,
              math::abs(vdet(0,i) - Jm.block(0, 3*i, 3, 3).determinant()));
    CHECK_CLOSE(0.0, err, 1e-11);

    // Each term differentiates exactly one factor in each direction, so every
    // term is C^0 at the interior knot and so is their sum: multiplicity
    // 5 - 0 = 5, giving 17 knots and 11 functions per direction.  The Bezier
    // form would be 12 per direction, so the compression really happened.
    for (short_t k = 0; k != 3; ++k)
        CHECK_EQUAL(5, (int)det.knots(k).multiplicity(0.5));
    CHECK_EQUAL(11*11*11, det.coefsSize());
}

} // SUITE
