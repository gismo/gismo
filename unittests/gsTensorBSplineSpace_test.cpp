/** @file gsTensorBSplineSpace_test.cpp

    @brief Unit tests for tensor B-spline space specs and transfer maps.
*/

#include "gismo_unittest.h"

SUITE(gsTensorBSplineSpace_test)
{

namespace
{

template <short_t d, class T>
class SplineViewFunction final : public gsFunction<T>
{
    const gsTensorBSpline<d,T>& m_spline;

public:
    explicit SplineViewFunction(const gsTensorBSpline<d,T>& spline)
        : m_spline(spline)
    { }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        m_spline.eval_into(u, result);
    }

    short_t domainDim() const override { return d; }
    short_t targetDim() const override { return m_spline.targetDim(); }

    GISMO_CLONE_FUNCTION(SplineViewFunction)
};

template <short_t d, class T>
class ProductSplineFunction final : public gsFunction<T>
{
    const gsTensorBSpline<d,T>& m_lhs;
    const gsTensorBSpline<d,T>& m_rhs;

public:
    ProductSplineFunction(const gsTensorBSpline<d,T>& lhs,
                          const gsTensorBSpline<d,T>& rhs)
        : m_lhs(lhs), m_rhs(rhs)
    { }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        gsMatrix<T> lhs, rhs;
        m_lhs.eval_into(u, lhs);
        m_rhs.eval_into(u, rhs);
        result = lhs.cwiseProduct(rhs);
    }

    short_t domainDim() const override { return d; }
    short_t targetDim() const override { return 1; }

    GISMO_CLONE_FUNCTION(ProductSplineFunction)
};

template <class T>
class GradSplineFunction final : public gsFunction<T>
{
    const gsTensorBSpline<2,T>& m_spline;
    short_t m_dir;

public:
    GradSplineFunction(const gsTensorBSpline<2,T>& spline, short_t dir)
        : m_spline(spline), m_dir(dir)
    { }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        gsMatrix<T> deriv;
        m_spline.deriv_into(u, deriv);
        result = deriv.row(m_dir);
    }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    GISMO_CLONE_FUNCTION(GradSplineFunction)
};

template <class T>
class LaplSplineFunction final : public gsFunction<T>
{
    const gsTensorBSpline<2,T>& m_spline;

public:
    explicit LaplSplineFunction(const gsTensorBSpline<2,T>& spline)
        : m_spline(spline)
    { }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        gsMatrix<T> deriv2;
        m_spline.deriv2_into(u, deriv2);
        result = deriv2.row(0) + deriv2.row(1);
    }

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    GISMO_CLONE_FUNCTION(LaplSplineFunction)
};

template <short_t d>
void checkSameBasis(const gsTensorBSplineBasis<d,real_t>& a,
                    const gsTensorBSplineBasis<d,real_t>& b)
{
    CHECK_EQUAL(a.size(), b.size());
    for (short_t dir = 0; dir < d; ++dir)
    {
        CHECK_EQUAL(a.degree(dir), b.degree(dir));
        const gsKnotVector<real_t>& akv = a.knots(dir);
        const gsKnotVector<real_t>& bkv = b.knots(dir);
        CHECK_EQUAL(akv.uSize(), bkv.uSize());
        for (index_t i = 0; i < static_cast<index_t>(akv.uSize()); ++i)
        {
            const real_t xi = akv.uValue(i);
            CHECK_EQUAL(xi, bkv.uValue(i));
            CHECK_EQUAL(akv.multiplicity(xi), bkv.multiplicity(xi));
        }
    }
}

template <short_t d>
gsMatrix<real_t> evalPts(index_t n)
{
    index_t total = 1;
    for (short_t dir = 0; dir < d; ++dir)
        total *= n;

    gsMatrix<real_t> pts(d, total);
    for (index_t k = 0; k < total; ++k)
    {
        index_t q = k;
        for (short_t dir = 0; dir < d; ++dir)
        {
            pts(dir, k) = (real_t(q % n) + real_t(0.5)) / real_t(n);
            q /= n;
        }
    }
    return pts;
}

template <short_t d>
void checkSplineMatchesFunction(const gsTensorBSpline<d,real_t>& spline,
                                const gsFunction<real_t>& func,
                                real_t tol)
{
    CHECK_EQUAL(static_cast<int>(d), static_cast<int>(func.domainDim()));
    CHECK_EQUAL(static_cast<int>(spline.targetDim()),
                static_cast<int>(func.targetDim()));

    gsMatrix<real_t> vs, vf;
    const gsMatrix<real_t> pts = evalPts<d>(d == 3 ? 3 : 5);
    spline.eval_into(pts, vs);
    func.eval_into(pts, vf);
    const real_t rel = (vs - vf).norm() / (vf.norm() + real_t(1e-300));
    CHECK(rel < tol);
}

template <short_t d>
void checkAnchorInterpolationMatches(const gsTensorBSpline<d,real_t>& spline,
                                     const gsFunction<real_t>& func,
                                     real_t tol)
{
    gsMatrix<real_t> values;
    func.eval_into(spline.basis().anchors(), values);
    const gsGeometry<real_t>::uPtr interpolated =
        spline.basis().interpolateAtAnchors(values);
    CHECK_EQUAL(spline.coefs().rows(), interpolated->coefs().rows());
    CHECK_EQUAL(spline.coefs().cols(), interpolated->coefs().cols());
    CHECK_ARRAY_CLOSE(spline.coefs().data(), interpolated->coefs().data(),
                      spline.coefs().size(), tol);
}

template <short_t d>
gsTensorBSpline<d,real_t> makeSpline(
    const gsTensorBSplineBasis<d,real_t>& basis,
    real_t scale = real_t(1))
{
    gsMatrix<real_t> coefs(basis.size(), 1);
    for (index_t i = 0; i < coefs.rows(); ++i)
        coefs(i, 0) = scale * (real_t(0.2)
                     + real_t(0.03) * math::sin(real_t(i + 1))
                     - real_t(0.01) * math::cos(real_t(2 * i + 1)));
    return gsTensorBSpline<d,real_t>(basis, coefs);
}

gsTensorBSplineBasis<2,real_t> makeBasis2(index_t degree, index_t nRefine)
{
    gsKnotVector<real_t> kv(0.0, 1.0, 1, degree + 1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    for (index_t r = 0; r < nRefine; ++r)
        basis.uniformRefine();
    return basis;
}

gsTensorBSplineBasis<2,real_t> makeNonUniformRepeatedBasis2()
{
    const real_t xKnots[] = {
        0.0, 0.0, 0.0, 0.0,
        0.2, 0.5, 0.5, 0.8,
        1.0, 1.0, 1.0, 1.0
    };
    const real_t yKnots[] = {
        0.0, 0.0, 0.0,
        0.3, 0.6,
        1.0, 1.0, 1.0
    };
    gsKnotVector<real_t> kvx(3, xKnots, xKnots + sizeof(xKnots) / sizeof(real_t));
    gsKnotVector<real_t> kvy(2, yKnots, yKnots + sizeof(yKnots) / sizeof(real_t));
    return gsTensorBSplineBasis<2,real_t>(kvx, kvy);
}

} // anonymous namespace

TEST(space_specs_match_existing_factory_bases_2d)
{
    const gsTensorBSplineBasis<2,real_t> basis = makeBasis2(3, 2);
    const gsTensorBSpline<2,real_t> c = makeSpline<2>(basis);

    typedef gsTensorBSplineSpaceSpec<2,real_t> Spec;
    const Spec s = Spec::fromBasis(basis);

    checkSameBasis(bezierSpace(s).toBasis(), c.toBezier().basis());
    checkSameBasis(powerSpace(s, 2).toBasis(), c.squared(false).basis());
    checkSameBasis(powerSpace(s, 2, gsTensorBSplineSpacePolicy::Bezier).toBasis(),
                   c.squared(true).basis());
    checkSameBasis(powerSpace(s, 3).toBasis(), c.cubed(false).basis());

    gsVector<index_t,2> dx;
    dx.setZero();
    dx[0] = 1;
    checkSameBasis(derivativeSpace(s, dx).toBasis(), c.grad(0).basis());

    checkSameBasis(laplacianSpace(s).toBasis(), c.lapl(false).basis());
    checkSameBasis(laplacianSpace(s, gsTensorBSplineSpacePolicy::Bezier).toBasis(),
                   c.lapl(true).basis());
}

TEST(transfer_matrix_reproduces_wrapped_source_function_in_superspace)
{
    const gsTensorBSplineBasis<2,real_t> basis = makeBasis2(2, 2);
    const gsTensorBSpline<2,real_t> c = makeSpline<2>(basis);

    typedef gsTensorBSplineSpaceSpec<2,real_t> Spec;
    const Spec source = Spec::fromBasis(basis);
    const Spec target = powerSpace(source, 2);

    const gsSparseMatrix<real_t,RowMajor> L = transferMatrix(source, target);
    const gsTensorBSplineBasis<2,real_t> targetBasis = target.toBasis();
    const gsTensorBSpline<2,real_t> lifted(targetBasis, gsMatrix<real_t>(L * c.coefs()));

    const SplineViewFunction<2,real_t> sourceFunction(c);
    checkSplineMatchesFunction(lifted, sourceFunction, real_t(1e-12));
}

TEST(bezier_transfer_matches_toBezier_coefficients)
{
    const gsTensorBSplineBasis<2,real_t> basis = makeBasis2(3, 1);
    const gsTensorBSpline<2,real_t> c = makeSpline<2>(basis);

    typedef gsTensorBSplineSpaceSpec<2,real_t> Spec;
    const Spec source = Spec::fromBasis(basis);
    const Spec target = bezierSpace(source);
    const gsSparseMatrix<real_t,RowMajor> L = transferMatrix(source, target);

    const gsTensorBSpline<2,real_t> bez = c.toBezier();
    CHECK((gsMatrix<real_t>(L * c.coefs()) - bez.coefs()).norm()
          / (bez.coefs().norm() + real_t(1e-300)) < real_t(1e-12));

    const SplineViewFunction<2,real_t> sourceFunction(c);
    checkSplineMatchesFunction(bez, sourceFunction, real_t(1e-12));
}

TEST(multiply_matches_function_interpolation)
{
    const gsTensorBSplineBasis<2,real_t> basis = makeBasis2(3, 1);
    const gsTensorBSpline<2,real_t> f = makeSpline<2>(basis, real_t(0.7));
    const gsTensorBSpline<2,real_t> g = makeSpline<2>(basis, real_t(-0.4));
    const gsTensorBSpline<2,real_t> fg = gsTensorBSpline<2,real_t>::multiply(f, g, false);

    const ProductSplineFunction<2,real_t> productFunction(f, g);
    checkSplineMatchesFunction(fg, productFunction, real_t(1e-12));
    checkAnchorInterpolationMatches(fg, productFunction, real_t(1e-10));
}

TEST(derivative_space_matches_grad_function_interpolation)
{
    const gsTensorBSplineBasis<2,real_t> basis = makeNonUniformRepeatedBasis2();
    const gsTensorBSpline<2,real_t> c = makeSpline<2>(basis);

    typedef gsTensorBSplineSpaceSpec<2,real_t> Spec;
    const Spec s = Spec::fromBasis(basis);

    for (short_t dir = 0; dir < 2; ++dir)
    {
        gsVector<index_t,2> orders;
        orders.setZero();
        orders[dir] = 1;

        const gsTensorBSpline<2,real_t> grad = c.grad(dir);
        checkSameBasis(derivativeSpace(s, orders).toBasis(), grad.basis());

        const GradSplineFunction<real_t> gradFunction(c, dir);
        checkSplineMatchesFunction(grad, gradFunction, real_t(1e-12));
        checkAnchorInterpolationMatches(grad, gradFunction, real_t(1e-10));
    }
}

TEST(laplacian_space_matches_lapl_function_interpolation)
{
    const gsTensorBSplineBasis<2,real_t> basis = makeNonUniformRepeatedBasis2();
    const gsTensorBSpline<2,real_t> c = makeSpline<2>(basis);

    typedef gsTensorBSplineSpaceSpec<2,real_t> Spec;
    const Spec s = Spec::fromBasis(basis);

    const LaplSplineFunction<real_t> laplFunction(c);

    const gsTensorBSpline<2,real_t> laplMin = c.lapl(false);
    checkSameBasis(laplacianSpace(s).toBasis(), laplMin.basis());
    checkSplineMatchesFunction(laplMin, laplFunction, real_t(1e-12));
    checkAnchorInterpolationMatches(laplMin, laplFunction, real_t(1e-10));

    const gsTensorBSpline<2,real_t> laplBez = c.lapl(true);
    checkSameBasis(laplacianSpace(s, gsTensorBSplineSpacePolicy::Bezier).toBasis(),
                   laplBez.basis());
    checkSplineMatchesFunction(laplBez, laplFunction, real_t(1e-12));
    checkAnchorInterpolationMatches(laplBez, laplFunction, real_t(1e-10));
}

TEST(nonuniform_repeated_asymmetric_product_matches_wrapped_function)
{
    const gsTensorBSplineBasis<2,real_t> basis = makeNonUniformRepeatedBasis2();
    const gsTensorBSpline<2,real_t> f = makeSpline<2>(basis, real_t(1.3));
    const gsTensorBSpline<2,real_t> g = makeSpline<2>(basis, real_t(-0.2));

    typedef gsTensorBSplineSpaceSpec<2,real_t> Spec;
    const Spec s = Spec::fromBasis(basis);
    checkSameBasis(productSpace(s, s).toBasis(), gsTensorBSpline<2,real_t>::multiply(f, g, false).basis());
    checkSameBasis(productSpace(s, s, gsTensorBSplineSpacePolicy::Bezier).toBasis(),
                   gsTensorBSpline<2,real_t>::multiply(f, g, true).basis());

    const ProductSplineFunction<2,real_t> productFunction(f, g);
    checkSplineMatchesFunction(gsTensorBSpline<2,real_t>::multiply(f, g, false), productFunction,
                               real_t(1e-12));
    checkSplineMatchesFunction(gsTensorBSpline<2,real_t>::multiply(f, g, true), productFunction,
                               real_t(1e-12));
}

TEST(product_space_merges_breaks_and_so_does_the_product)
{
    const real_t aKnots[] = {0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0};
    const real_t bKnots[] = {0.0, 0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0, 1.0};
    gsKnotVector<real_t> kva(2, aKnots, aKnots + sizeof(aKnots) / sizeof(real_t));
    gsKnotVector<real_t> kvb(2, bKnots, bKnots + sizeof(bKnots) / sizeof(real_t));
    gsTensorBSplineBasis<2,real_t> basisA(kva, kva);
    gsTensorBSplineBasis<2,real_t> basisB(kvb, kvb);

    typedef gsTensorBSplineSpaceSpec<2,real_t> Spec;
    const Spec a = Spec::fromBasis(basisA);
    const Spec b = Spec::fromBasis(basisB);
    const Spec product = productSpace(a, b);

    CHECK(product.knots(0).has(real_t(0.25)));
    CHECK(product.knots(0).has(real_t(0.5)));
    CHECK(product.knots(0).has(real_t(0.75)));

    const gsTensorBSpline<2,real_t> f = makeSpline<2>(basisA);
    const gsTensorBSpline<2,real_t> g = makeSpline<2>(basisB);

    // multiply() unions the breakpoints of the two factors, so it lands in
    // exactly the space productSpace() predicts and still reproduces f*g.
    const gsTensorBSpline<2,real_t> fg = gsTensorBSpline<2,real_t>::multiply(f, g);
    checkSameBasis(product.toBasis(), fg.basis());

    const ProductSplineFunction<2,real_t> productFunction(f, g);
    checkSplineMatchesFunction(fg, productFunction, real_t(1e-12));
}

TEST(product_space_treats_absent_breakpoint_as_smooth_factor)
{
    const real_t lowKnots[] = {0.0, 0.0, 1.0, 1.0};
    const real_t highKnots[] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.25,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0
    };
    gsKnotVector<real_t> low(1, lowKnots,
                             lowKnots + sizeof(lowKnots) / sizeof(real_t));
    gsKnotVector<real_t> high(5, highKnots,
                              highKnots + sizeof(highKnots) / sizeof(real_t));

    typedef gsTensorBSplineSpaceSpec<1,real_t> Spec;
    const Spec product = productSpace(
        Spec::fromBasis(gsTensorBSplineBasis<1,real_t>(low)),
        Spec::fromBasis(gsTensorBSplineBasis<1,real_t>(high)));

    CHECK_EQUAL(6, product.degree(0));
    CHECK(product.knots(0).has(real_t(0.25)));
    CHECK_EQUAL(2, product.knots(0).multiplicity(real_t(0.25)));
}

TEST(space_specs_reject_invalid_inputs)
{
    gsTensorBSplineBasis<2,real_t> basis = makeBasis2(2, 2);
    typedef gsTensorBSplineSpaceSpec<2,real_t> Spec;
    const Spec s = Spec::fromBasis(basis);

    CHECK_THROW(powerSpace(s, 0), std::runtime_error);

    gsVector<index_t,2> tooHigh;
    tooHigh.setZero();
    tooHigh[0] = static_cast<index_t>(basis.degree(0)) + 1;
    CHECK_THROW(derivativeSpace(s, tooHigh), std::runtime_error);

    const real_t shiftedKnots[] = {0.0, 0.0, 0.0, 0.4, 1.2, 1.2, 1.2};
    gsKnotVector<real_t> shifted(2, shiftedKnots,
                                 shiftedKnots + sizeof(shiftedKnots) / sizeof(real_t));
    const Spec shiftedSpec =
        Spec::fromBasis(gsTensorBSplineBasis<2,real_t>(shifted, shifted));
    CHECK_THROW(productSpace(s, shiftedSpec), std::runtime_error);

    CHECK_THROW(transferMatrix(powerSpace(s, 2), s), std::runtime_error);

    const real_t sourceKnots[] = {0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0};
    const real_t targetKnots[] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    gsKnotVector<real_t> sourceKv(2, sourceKnots,
                                  sourceKnots + sizeof(sourceKnots) / sizeof(real_t));
    gsKnotVector<real_t> targetKv(2, targetKnots,
                                  targetKnots + sizeof(targetKnots) / sizeof(real_t));
    typedef gsTensorBSplineSpaceSpec<1,real_t> Spec1;
    CHECK_THROW(transferMatrix(
        Spec1::fromBasis(gsTensorBSplineBasis<1,real_t>(sourceKv)),
        Spec1::fromBasis(gsTensorBSplineBasis<1,real_t>(targetKv))),
        std::runtime_error);

    gsTensorBSplineBasis<2,real_t> periodicBasis = makeBasis2(2, 3);
    periodicBasis.setPeriodic(0);
    CHECK(periodicBasis.isPeriodic());
    CHECK_THROW(Spec::fromBasis(periodicBasis), std::runtime_error);
}

TEST(space_specs_and_transfers_instantiate_for_1d_and_3d)
{
    gsKnotVector<real_t> kv1(0.0, 1.0, 2, 3);
    gsTensorBSplineBasis<1,real_t> b1(kv1);
    typedef gsTensorBSplineSpaceSpec<1,real_t> Spec1;
    const Spec1 s1 = Spec1::fromBasis(b1);
    CHECK_EQUAL(4, powerSpace(s1, 2).degree(0));

    const gsTensorBSpline<1,real_t> c1 = makeSpline<1>(b1);
    const Spec1 target1 = bezierSpace(s1);
    const gsSparseMatrix<real_t,RowMajor> L1 = transferMatrix(s1, target1);
    const gsTensorBSpline<1,real_t> lifted1(target1.toBasis(),
                                            gsMatrix<real_t>(L1 * c1.coefs()));
    const SplineViewFunction<1,real_t> f1(c1);
    checkSplineMatchesFunction(lifted1, f1, real_t(1e-12));

    gsKnotVector<real_t> kv3(0.0, 1.0, 1, 3);
    gsTensorBSplineBasis<3,real_t> b3(kv3, kv3, kv3);
    b3.uniformRefine();
    typedef gsTensorBSplineSpaceSpec<3,real_t> Spec3;
    const Spec3 s3 = Spec3::fromBasis(b3);
    gsVector<index_t,3> orders;
    orders.setZero();
    orders[2] = 1;
    CHECK_EQUAL(1, derivativeSpace(s3, orders).degree(2));
    CHECK_EQUAL(2, derivativeSpace(s3, orders).degree(0));

    const gsTensorBSpline<3,real_t> c3 = makeSpline<3>(b3);
    const Spec3 target3 = bezierSpace(s3);
    const gsSparseMatrix<real_t,RowMajor> L3 = transferMatrix(s3, target3);
    const gsTensorBSpline<3,real_t> lifted3(target3.toBasis(),
                                            gsMatrix<real_t>(L3 * c3.coefs()));
    const SplineViewFunction<3,real_t> f3(c3);
    checkSplineMatchesFunction(lifted3, f3, real_t(1e-12));
}

} // SUITE
