/** @file gsAdaptiveParametrization_test.cpp

    @brief Tests for gsOptMesh gradient computation (ValueBased & GradientBased modes)

    Compares analytical gradients from gradObj_into against finite-difference
    gradients from gradObj_FD_into. Tolerance depends on whether the build
    has ADIFF enabled: with ADIFF, gsFunctionExpr derivatives are exact
    (~1e-9 accuracy); without ADIFF, they use FD internally (~1e-4 accuracy).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
 **/

#include "gismo_unittest.h"
#include <gsNurbs/gsSquareDomain.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsOptimizer/gsGradientDescent.h>
#include <gsUtils/gsStopwatch.h>

#include <algorithm>
#include <array>

#ifdef GISMO_WITH_ADIFF
static const real_t GRAD_TOL = 1e-6;
#else
static const real_t GRAD_TOL = 1e-2;
#endif

namespace {

template <MonitorMode MODE>
real_t gradientRelError(gsOptMesh<real_t,MODE> & opt,
                        const gsVector<real_t> & controls)
{
    const index_t nc = controls.size();

    gsVector<real_t> grad(nc);
    gsAsVector<real_t> asgrad(grad.data(), grad.rows());
    opt.gradObj_into(gsAsConstVector<real_t>(controls.data(), controls.size()), asgrad);

    gsVector<real_t> grad_fd(nc);
    gsAsVector<real_t> asgrad_fd(grad_fd.data(), grad_fd.rows());
    opt.gradObj_FD_into(gsAsConstVector<real_t>(controls.data(), controls.size()), asgrad_fd);

    real_t absErr = (grad - grad_fd).norm();
    real_t relErr = absErr / (grad_fd.norm() > 1e-15 ? grad_fd.norm() : 1.0);
    return relErr;
}

struct TestFixture
{
    gsGeometry<>::uPtr composition;
    gsSquareDomain<real_t> domain;
    gsVector<real_t> controls;
    gsKnotVector<> kv1;
    gsKnotVector<> kv2;
    gsTensorBSplineBasis<2> tbasis;
    gsComposedBasis<> cbasis;
    gsMatrix<> anchors;

    TestFixture()
        : composition(gsNurbsCreator<>::BSplineSquareDeg(2))
        , domain(*composition)
        , controls(domain.nControls())
        , kv1({0,0,0,1,1,1}, 2)
        , kv2({0,0,0,0.50,1,1,1}, 2)
        , tbasis(kv1, kv2)
        , cbasis(*composition, tbasis)
        , anchors(cbasis.anchors().transpose())
    {
        controls << 0.95, 0.95;
    }

    gsGeometry<>::uPtr makeSurfaceGeom()
    {
        gsMatrix<> coefs3(cbasis.size(), 3);
        coefs3.leftCols(2) = anchors;
        for (index_t i = 0; i < coefs3.rows(); ++i)
            coefs3(i,2) = math::sin(coefs3(i,0)*EIGEN_PI) * math::cos(coefs3(i,1)*EIGEN_PI);
        return tbasis.makeGeometry(coefs3);
    }

    gsGeometry<>::uPtr makePlanarGeom()
    {
        gsMatrix<> coefs2(cbasis.size(), 2);
        coefs2 = anchors;
        for (index_t i = 0; i < coefs2.rows(); ++i)
        {
            coefs2(i,0) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,1));
            coefs2(i,1) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,0));
        }
        return tbasis.makeGeometry(coefs2);
    }

    gsGeometry<>::uPtr makeIdentityGeom()
    {
        gsMatrix<> coefs2(cbasis.size(), 2);
        coefs2 = anchors;
        return tbasis.makeGeometry(coefs2);
    }
};

// Element corners {lo_x, lo_y, up_x, up_y} of every element of \a b, sorted
// lexicographically. Works uniformly for a hierarchical (gsHTensorBasis) or
// tensor basis, since both expose the same domain-iterator interface -- the
// comparison below is at the element-partition level, not the tree level, so
// it can compare a hierarchical result against a tensor one.
std::vector<std::array<real_t,4>> elementSignatures(const gsBasis<real_t> & b)
{
    std::vector<std::array<real_t,4>> result;
    gsBasis<real_t>::domainIter domIt    = b.domain()->beginAll();
    gsBasis<real_t>::domainIter domItEnd = b.domain()->endAll();
    for (; domIt != domItEnd; ++domIt)
    {
        std::array<real_t,4> sig = { domIt.lowerCorner()(0), domIt.lowerCorner()(1),
                                      domIt.upperCorner()(0), domIt.upperCorner()(1) };
        result.push_back(sig);
    }
    std::sort(result.begin(), result.end());
    return result;
}

// True iff \a a and \a b have the same number of elements and the same
// element corners up to \a tol, order-independent.
bool sameElements(const gsBasis<real_t> & a, const gsBasis<real_t> & b, real_t tol = 1e-14)
{
    std::vector<std::array<real_t,4>> sa = elementSignatures(a);
    std::vector<std::array<real_t,4>> sb = elementSignatures(b);
    if (sa.size() != sb.size())
        return false;
    for (size_t i = 0; i != sa.size(); ++i)
        for (short_t j = 0; j != 4; ++j)
            if (math::abs(sa[i][j] - sb[i][j]) > tol)
                return false;
    return true;
}

// True iff some element of \a b has a lower or upper corner within \a tol of
// \a x in direction \a dir -- the release-safe way to check that a knot line
// survived into the integration mesh, since gsHTensorBasis::merge's own
// nesting check is compiled out under -DNDEBUG.
bool hasBreak(const gsBasis<real_t> & b, short_t dir, real_t x, real_t tol = 1e-12)
{
    gsBasis<real_t>::domainIter domIt    = b.domain()->beginAll();
    gsBasis<real_t>::domainIter domItEnd = b.domain()->endAll();
    for (; domIt != domItEnd; ++domIt)
    {
        if (math::abs(domIt.lowerCorner()(dir) - x) <= tol) return true;
        if (math::abs(domIt.upperCorner()(dir) - x) <= tol) return true;
    }
    return false;
}

// Geometry S on an arbitrary tensor B-spline basis \a tb, following
// TestFixture::makePlanarGeom but usable with a fixture-local basis.
gsGeometry<real_t>::uPtr makePlanarS(const gsTensorBSplineBasis<2,real_t> & tb)
{
    gsMatrix<real_t> coefs = tb.anchors().transpose();
    for (index_t i = 0; i < coefs.rows(); ++i)
    {
        const real_t x = coefs(i,0), y = coefs(i,1);
        coefs(i,0) = x + 0.05 * math::sin(2.0 * EIGEN_PI * y);
        coefs(i,1) = y + 0.05 * math::sin(2.0 * EIGEN_PI * x);
    }
    return tb.makeGeometry(coefs);
}

// Identity sigma basis at degree \a sigmaDeg with sigmaKnots = 2^sigmaRef - 1
// interior knots, in the exact knot-vector convention
// gsSquareDomain(const gsBasis<T>&, bool) expects.
gsTensorBSplineBasis<2,real_t> makeSigmaBasis(short_t sigmaDeg, index_t sigmaRef)
{
    gsKnotVector<real_t> ks(0, 1, (1 << sigmaRef) - 1, sigmaDeg + 1);
    return gsTensorBSplineBasis<2,real_t>(ks, ks);
}

} // anonymous namespace

SUITE(gsAdaptiveParametrizatin_test)
{

#ifndef GISMO_WITH_ADIFF
    TEST(ADIFF_warning)
    {
        gsWarn << "GISMO_WITH_ADIFF is not enabled. "
               << "GradientBased gradient tests use relaxed tolerance ("
               << GRAD_TOL << ") because gsFunctionExpr derivatives "
               << "fall back to finite differences.\n";
        CHECK(true);
    }
#endif

    TEST(B1_Surface_NoMonitor)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makeSurfaceGeom();
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &f.tbasis);
        real_t relErr = gradientRelError(opt, f.controls);
        gsTestInfo << "B1 relErr = " << relErr << "\n";
        CHECK(relErr < GRAD_TOL);
    }

    TEST(A1_Planar_NoMonitor)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &f.tbasis);
        real_t relErr = gradientRelError(opt, f.controls);
        gsTestInfo << "A1 relErr = " << relErr << "\n";
        CHECK(relErr < GRAD_TOL);
    }

    TEST(B2_Surface_ValueBased)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makeSurfaceGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        real_t relErr = gradientRelError(opt, f.controls);
        gsTestInfo << "B2 relErr = " << relErr << "\n";
        CHECK(relErr < GRAD_TOL);
    }

    TEST(A2_Planar_ValueBased)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        real_t relErr = gradientRelError(opt, f.controls);
        gsTestInfo << "A2 relErr = " << relErr << "\n";
        CHECK(relErr < GRAD_TOL);
    }

    TEST(C1_Surface_GradientBased_Parametric)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makeSurfaceGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        real_t relErr = gradientRelError<MonitorMode::GradientBased>(opt, f.controls);
        gsTestInfo << "C1 relErr = " << relErr << "\n";
        CHECK(relErr < GRAD_TOL);
    }

    TEST(C2_Planar_GradientBased_Parametric)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        real_t relErr = gradientRelError<MonitorMode::GradientBased>(opt, f.controls);
        gsTestInfo << "C2 relErr = " << relErr << "\n";
        CHECK(relErr < GRAD_TOL);
    }

    TEST(C3_Surface_GradientBased_Physical)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makeSurfaceGeom();
        gsFunctionExpr<> monFun("1 + 0.3*x + 0.2*y + 0.1*z", 3);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, false);
        real_t relErr = gradientRelError<MonitorMode::GradientBased>(opt, f.controls);
        CHECK(relErr < GRAD_TOL);
    }

    TEST(C4_Planar_GradientBased_Physical)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*x*x + 0.3*y*y", 2);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, false);
        real_t relErr = gradientRelError<MonitorMode::GradientBased>(opt, f.controls);
        gsTestInfo << "C4 relErr = " << relErr << "\n";
        CHECK(relErr < GRAD_TOL);
    }

    TEST(C5_Planar_GradientBased_LinearMonitor)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("0.3*x + 0.7*y", 2);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        real_t relErr = gradientRelError<MonitorMode::GradientBased>(opt, f.controls);
        gsTestInfo << "C5 relErr = " << relErr << "\n";
        CHECK(relErr < GRAD_TOL);
    }

    TEST(C6_IdentityGeom_GradientBased_LinearMonitor)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makeIdentityGeom();
        gsFunctionExpr<> monFun("0.3*x + 0.7*y", 2);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        real_t relErr = gradientRelError<MonitorMode::GradientBased>(opt, f.controls);
        gsTestInfo << "C6 relErr = " << relErr << "\n";
        CHECK(relErr < GRAD_TOL);
    }

    TEST(C7_GradientBased_ThetaZero)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        opt.options().setReal("Smoothing", 0.0);
        real_t relErr = gradientRelError<MonitorMode::GradientBased>(opt, f.controls);
        gsTestInfo << "C7 relErr = " << relErr << "\n";
        CHECK(relErr < GRAD_TOL);
    }

    // -----------------------------------------------------------------
    // gsOptFit / gsOptL2 analytic gradients agree with gsCheckSigmaGradient()'s
    // central-difference self-check. Both are built away from the
    // identity/stationary point so the check is not vacuous
    // (gsCheckSigmaGradient() skips components whose gradient AND FD are
    // both below 1e-10).
    // -----------------------------------------------------------------

    TEST(D1_OptFit_GradientMatchesFD_BarrierInactive)
    {
        TestFixture f;
        gsGeometry<>::uPtr S = f.makeSurfaceGeom();

        // Displace the sigma controls off the identity by a deterministic
        // offset (not gsSquareDomain::perturb(), which uses unseeded rand()).
        gsVector<real_t> c = f.domain.getControls();
        c.array() += 0.05;
        f.domain.setControls(c);

        gsVector<real_t> a(2), b(2);
        a.setZero(); b.setOnes();
        gsVector<unsigned> np(2); np << 6, 6;
        gsMatrix<real_t> uv = gsPointGrid(a, b, np);

        gsMatrix<real_t> xyz;
        S->eval_into(uv, xyz);
        // Perturb the target points by a fixed smooth offset so the residual
        // does not vanish.
        for (index_t i = 0; i != xyz.cols(); i++)
            xyz(2, i) += 0.1 * math::sin(3.0 * EIGEN_PI * uv(0, i));

        gsOptFit<real_t> opt(f.domain, *S, uv, xyz, /*mu*/0.0, /*eps*/0.1,
                             gsFoldBarrierMode::Sampled, /*quB*/2);
        real_t err = gsCheckSigmaGradient(opt, f.domain);
        gsTestInfo << "D1 OptFit (barrier inactive) checkGradient = " << err << "\n";
        CHECK(err < 1e-5);
    }

    TEST(D2_OptFit_GradientMatchesFD_BarrierActive)
    {
        TestFixture f;
        gsGeometry<>::uPtr S = f.makeSurfaceGeom();

        gsVector<real_t> c = f.domain.getControls();
        c.array() += 0.05;
        f.domain.setControls(c);

        gsVector<real_t> a(2), b(2);
        a.setZero(); b.setOnes();
        gsVector<unsigned> np(2); np << 6, 6;
        gsMatrix<real_t> uv = gsPointGrid(a, b, np);

        gsMatrix<real_t> xyz;
        S->eval_into(uv, xyz);
        for (index_t i = 0; i != xyz.cols(); i++)
            xyz(2, i) += 0.1 * math::sin(3.0 * EIGEN_PI * uv(0, i));

        // mu=0 baseline objective at the same controls, to verify the
        // barrier is genuinely active in the mu>0 case below.
        gsOptFit<real_t> opt0(f.domain, *S, uv, xyz, /*mu*/0.0, /*eps*/2.0,
                             gsFoldBarrierMode::Sampled, /*quB*/2);
        gsVector<real_t> u = f.domain.getControls();
        real_t f0 = opt0.evalObj(gsAsConstVector<real_t>(u.data(), u.rows()));

        // eps=2.0 is well above det(J_sigma) (~1 near this mild perturbation),
        // so the barrier term is active at every quadrature point.
        gsOptFit<real_t> optMu(f.domain, *S, uv, xyz, /*mu*/10.0, /*eps*/2.0,
                               gsFoldBarrierMode::Sampled, /*quB*/2);
        real_t fMu = optMu.evalObj(gsAsConstVector<real_t>(u.data(), u.rows()));

        gsTestInfo << "D2 f0=" << f0 << " fMu=" << fMu << " |diff|=" << math::abs(fMu - f0) << "\n";
        CHECK(math::abs(fMu - f0) > 1e-12);

        real_t err = gsCheckSigmaGradient(optMu, f.domain);
        gsTestInfo << "D2 OptFit (barrier active) checkGradient = " << err << "\n";
        CHECK(err < 1e-5);
    }

    TEST(D3_OptL2_GradientMatchesFD)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom(); // gsOptL2 requires a
                                                       // planar (domainDim==
                                                       // targetDim==2) S

        gsVector<real_t> c = f.domain.getControls();
        c.array() += 0.05;
        f.domain.setControls(c);

        gsFunctionExpr<> solution("sin(pi*x)*cos(pi*y)", 2); // frozen u_h(v)
        gsFunctionExpr<> target("1 + 0.5*sin(pi*x)*sin(pi*y)", 2); // target f

        gsOptL2<real_t> opt(f.domain, *geom, solution, target, &f.tbasis,
                            /*mu*/0.0, /*eps*/0.1, gsFoldBarrierMode::Sampled, /*quB*/2,
                            /*parametric*/false);
        real_t err = gsCheckSigmaGradient(opt, f.domain);
        gsTestInfo << "D3 OptL2 checkGradient = " << err << "\n";
        CHECK(err < 1e-5);
    }

    TEST(D4_OptL2_Surface_GradientMatchesFD)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makeSurfaceGeom();   // 2 -> 3

        gsVector<real_t> c = f.domain.getControls();
        c.array() += 0.05;
        f.domain.setControls(c);

        gsFunctionExpr<> solution("sin(pi*x)*cos(pi*y)", 2);          // u_h(v), v in [0,1]^2
        gsFunctionExpr<> target("1 + 0.5*sin(pi*x)*sin(pi*y)*z", 3);  // f in AMBIENT x,y,z

        gsOptL2<real_t> opt(f.domain, *geom, solution, target, &f.tbasis,
                            /*mu*/0.0, /*eps*/0.1, gsFoldBarrierMode::Sampled, /*quB*/2,
                            /*parametric*/false);
        real_t err = gsCheckSigmaGradient(opt, f.domain);
        gsTestInfo << "D4 OptL2 surface checkGradient = " << err << "\n";
        CHECK(err < 1e-5);
    }

    // -----------------------------------------------------------------
    // gsFoldBarrierMode::Coefficient. D5 pins the oracle (the per-element
    // Bezier convolution reproduces the factory's multiply(...,keepBezier=true)
    // coefficients); D6 is the analytic-gradient gate for the new mode; D7 is
    // the certificate ordering that must always hold. D5/D7 use an ASYMMETRIC
    // bidegree (degree 3 x degree 2, different element counts per direction)
    // instead of TestFixture's symmetric degree-2x2 sigma: term 1 (a*b') and
    // term 2 (a'*b) of the Leibniz product use per-direction weight tables
    // built from DIFFERENT operand bidegrees, and with degU==degV those
    // tables happen to line up in a way that would hide a "reused the wrong
    // table" bug.
    // -----------------------------------------------------------------

    TEST(D5_FoldBarrier_CoefficientBoundMatchesFactory)
    {
        typedef gsTensorBSpline<2,real_t> Spline;

        gsKnotVector<real_t> kvU(0, 1, 2, 4); // degree 3, 3 elements
        gsKnotVector<real_t> kvV(0, 1, 1, 3); // degree 2, 2 elements
        gsTensorBSplineBasis<2,real_t> basis(kvU, kvV);
        gsSquareDomain<real_t> domain(basis, true);

        gsVector<real_t> c = domain.getControls();
        c.array() += 0.05;
        domain.setControls(c);

        gsFoldBarrier<real_t> barrier(domain, /*mu*/1.0, /*eps*/0.0,
                                      gsFoldBarrierMode::Coefficient, /*quB*/0);
        const real_t barrierMin = barrier.minDetJacobian();
        // This mild perturbation's Bezier minimum lands on a u-edge element,
        // where d_v(sigma_0) is IDENTICALLY zero (Slide=true fixes sigma_0
        // there, and the fixed value is x=const along that edge) -- so this
        // test pins term 1's convolution and the Bezier layout identity, not
        // term 2's sign; term 2's sign is what D6's FD gradient check pins
        // (a sign error there desyncs the objective from its own gradient
        // everywhere, not just at one boundary-adjacent coefficient).

        const gsMatrix<real_t> & coefs = domain.domain().coefs();
        gsMatrix<real_t> c0 = coefs.col(0), c1 = coefs.col(1);
        Spline x(basis, c0), y(basis, c1);
        Spline det = Spline::linearCombination(
                         1.0, Spline::multiply(x.grad(0), y.grad(1), true),
                        -1.0, Spline::multiply(x.grad(1), y.grad(0), true));
        const real_t factoryMin = det.coefs().minCoeff();

        gsTestInfo << "D5 barrierMin=" << barrierMin << " factoryMin=" << factoryMin << "\n";
        CHECK(math::abs(factoryMin) > 1e-8);
        CHECK(math::abs(barrierMin - factoryMin) / math::abs(factoryMin) < 1e-9);

        // The check above pins one scalar of m_c (its min) against the
        // factory oracle -- not the Kronecker-factored extraction as a
        // whole. Extend it to every coefficient at
        // once via addObj()'s public Coefficient-mode formula (gsFoldBarrier.hpp
        // addObj()), E = mu*sum_i max(0,eps-c_i)^2 + box(controls): two
        // evaluations of the SAME domain/mu (so the box term is identical and
        // cancels exactly) isolate sum_i(eps-c_i)^2 --
        //   E0: eps=0 -- every c_i > 0 here (barrierMin>0, checked below), so
        //       the clamp is zero everywhere and E0 IS the box term alone.
        //   E1: eps = max(det.coefs())+1 -- every c_i < eps, so the clamp
        //       never fires and E1-E0 = sum_i(eps-c_i)^2 exactly.
        // A wrong 1-D factor (a role swap between m_Du/m_Ev/m_Cu/m_Fv, a
        // transposed reshape, or a scale error) perturbs nearly every c_i, so
        // this is sensitive to the WHOLE extraction rather than one entry --
        // and the reference (det.coefs()) is gsTensorBSpline's own 2-D
        // grad()+toBezier()+multiply(), never gsFoldBarrier's internals, so
        // the check does not validate against its own output.
        CHECK(barrierMin > 0);
        const real_t epsFull = det.coefs().maxCoeff() + 1.0;

        gsFoldBarrier<real_t> barrier0(domain, /*mu*/1.0, /*eps*/0.0,
                                       gsFoldBarrierMode::Coefficient, /*quB*/0);
        real_t E0 = 0; barrier0.addObj(E0);

        gsFoldBarrier<real_t> barrier1(domain, /*mu*/1.0, /*eps*/epsFull,
                                       gsFoldBarrierMode::Coefficient, /*quB*/0);
        real_t E1 = 0; barrier1.addObj(E1);

        const real_t sumSqRef     = (epsFull - det.coefs().array()).square().sum();
        const real_t sumSqBarrier = E1 - E0;

        gsTestInfo << "D5 nPoints=" << barrier.nPoints() << " detCoefs.rows()=" << det.coefs().rows()
                   << " sumSqBarrier=" << sumSqBarrier << " sumSqRef=" << sumSqRef << "\n";
        CHECK(barrier.nPoints() == det.coefs().rows());
        CHECK(math::abs(sumSqRef) > 1e-8);
        CHECK(math::abs(sumSqBarrier - sumSqRef) / math::abs(sumSqRef) < 1e-12);
    }

    // Same "pin the full Bezier coefficient array" check as D5's second
    // half, on a basis with the degree/element-count asymmetry
    // REVERSED (v now has the larger degree AND more elements) -- catches a
    // m_Du<->m_Cu or m_Ev<->m_Fv role swap that D5's basis alone could hide,
    // since both are non-square but D5 alone always has u as the "larger"
    // direction.
    TEST(D5b_FoldBarrier_CoefficientExtractionMatchesFactory_ReversedAsymmetry)
    {
        typedef gsTensorBSpline<2,real_t> Spline;

        gsKnotVector<real_t> kvU(0, 1, 1, 3); // degree 2, 2 elements
        gsKnotVector<real_t> kvV(0, 1, 2, 4); // degree 3, 3 elements
        gsTensorBSplineBasis<2,real_t> basis(kvU, kvV);
        gsSquareDomain<real_t> domain(basis, true);

        gsVector<real_t> c = domain.getControls();
        c.array() += 0.05;
        domain.setControls(c);

        const gsMatrix<real_t> & coefs = domain.domain().coefs();
        gsMatrix<real_t> c0 = coefs.col(0), c1 = coefs.col(1);
        Spline x(basis, c0), y(basis, c1);
        Spline det = Spline::linearCombination(
                         1.0, Spline::multiply(x.grad(0), y.grad(1), true),
                        -1.0, Spline::multiply(x.grad(1), y.grad(0), true));

        gsFoldBarrier<real_t> barrier0(domain, /*mu*/1.0, /*eps*/0.0,
                                       gsFoldBarrierMode::Coefficient, /*quB*/0);
        const real_t barrierMin0 = barrier0.minDetJacobian();
        CHECK(barrierMin0 > 0); // needed for E0 below to equal the box term alone
        real_t E0 = 0; barrier0.addObj(E0);

        const real_t epsFull = det.coefs().maxCoeff() + 1.0;
        gsFoldBarrier<real_t> barrier1(domain, /*mu*/1.0, /*eps*/epsFull,
                                       gsFoldBarrierMode::Coefficient, /*quB*/0);
        real_t E1 = 0; barrier1.addObj(E1);

        const real_t sumSqRef     = (epsFull - det.coefs().array()).square().sum();
        const real_t sumSqBarrier = E1 - E0;

        gsTestInfo << "D5b nPoints=" << barrier0.nPoints() << " detCoefs.rows()=" << det.coefs().rows()
                   << " sumSqBarrier=" << sumSqBarrier << " sumSqRef=" << sumSqRef << "\n";
        CHECK(barrier0.nPoints() == det.coefs().rows());
        CHECK(math::abs(sumSqRef) > 1e-8);
        CHECK(math::abs(sumSqBarrier - sumSqRef) / math::abs(sumSqRef) < 1e-12);
    }

    TEST(D6_OptFit_CoefficientBarrier_GradientMatchesFD)
    {
        TestFixture f;
        gsGeometry<>::uPtr S = f.makeSurfaceGeom();

        gsVector<real_t> c = f.domain.getControls();
        c.array() += 0.05;
        f.domain.setControls(c);

        gsVector<real_t> a(2), b(2);
        a.setZero(); b.setOnes();
        gsVector<unsigned> np(2); np << 6, 6;
        gsMatrix<real_t> uv = gsPointGrid(a, b, np);

        gsMatrix<real_t> xyz;
        S->eval_into(uv, xyz);
        for (index_t i = 0; i != xyz.cols(); i++)
            xyz(2, i) += 0.1 * math::sin(3.0 * EIGEN_PI * uv(0, i));

        // mu=0 baseline, to verify the coefficient barrier is genuinely
        // active in the mu>0 case below (same pattern as D2).
        gsOptFit<real_t> opt0(f.domain, *S, uv, xyz, /*mu*/0.0, /*eps*/2.0,
                              gsFoldBarrierMode::Coefficient, /*quB*/0);
        gsVector<real_t> u = f.domain.getControls();
        real_t f0 = opt0.evalObj(gsAsConstVector<real_t>(u.data(), u.rows()));

        gsOptFit<real_t> optMu(f.domain, *S, uv, xyz, /*mu*/10.0, /*eps*/2.0,
                               gsFoldBarrierMode::Coefficient, /*quB*/0);
        real_t fMu = optMu.evalObj(gsAsConstVector<real_t>(u.data(), u.rows()));

        gsTestInfo << "D6 f0=" << f0 << " fMu=" << fMu << " |diff|=" << math::abs(fMu - f0) << "\n";
        CHECK(math::abs(fMu - f0) > 1e-12);

        real_t err = gsCheckSigmaGradient(optMu, f.domain);
        gsTestInfo << "D6 OptFit (coefficient barrier active) checkGradient = " << err << "\n";
        CHECK(err < 1e-5);

        // Second fixture, same TEST: an asymmetric MULTI-ELEMENT sigma (3
        // elements in u, 2 in v; degree 3x2). TestFixture's sigma above is a
        // SINGLE Bezier element in both directions (BSplineSquareDeg(2),
        // knot vector {0,0,0,1,1,1}) -- every element-offset term in
        // addCoefficientGrad()'s closed-form index formulas
        // ((eu*m_degU+i0)+nA0*(ev*(m_degV+1)+i1), etc.) is then zero for
        // eu=ev=0, so a misindexed block extraction or scatter in the
        // adjoint would pass D6's first half, D5 and D7 (which only exercise
        // the forward path) undetected. This repeats D6's mu=0-vs-mu>0 +
        // gsCheckSigmaGradient pattern on the same asymmetric basis D5/D7
        // use, which does have eu,ev > 0 element offsets.
        gsKnotVector<real_t> kvU2(0, 1, 2, 4); // degree 3, 3 elements
        gsKnotVector<real_t> kvV2(0, 1, 1, 3); // degree 2, 2 elements
        gsTensorBSplineBasis<2,real_t> basis2(kvU2, kvV2);
        gsSquareDomain<real_t> domain2(basis2, true);
        gsVector<real_t> c2 = domain2.getControls();
        c2.array() += 0.05;
        domain2.setControls(c2);

        gsGeometry<real_t>::uPtr S2 = makePlanarS(basis2);
        gsMatrix<real_t> uv2 = gsPointGrid(a, b, np);
        gsMatrix<real_t> xyz2;
        S2->eval_into(uv2, xyz2);
        for (index_t i = 0; i != xyz2.cols(); i++)
            xyz2(1, i) += 0.1 * math::sin(3.0 * EIGEN_PI * uv2(0, i));

        gsOptFit<real_t> opt0b(domain2, *S2, uv2, xyz2, /*mu*/0.0, /*eps*/2.0,
                               gsFoldBarrierMode::Coefficient, /*quB*/0);
        gsVector<real_t> u2 = domain2.getControls();
        real_t f0b = opt0b.evalObj(gsAsConstVector<real_t>(u2.data(), u2.rows()));

        gsOptFit<real_t> optMub(domain2, *S2, uv2, xyz2, /*mu*/10.0, /*eps*/2.0,
                                gsFoldBarrierMode::Coefficient, /*quB*/0);
        real_t fMub = optMub.evalObj(gsAsConstVector<real_t>(u2.data(), u2.rows()));

        gsTestInfo << "D6b (multi-element) f0=" << f0b << " fMu=" << fMub
                   << " |diff|=" << math::abs(fMub - f0b) << "\n";
        CHECK(math::abs(fMub - f0b) > 1e-12);

        real_t errB = gsCheckSigmaGradient(optMub, domain2);
        gsTestInfo << "D6b OptFit (coefficient barrier active, multi-element) "
                      "checkGradient = " << errB << "\n";
        CHECK(errB < 1e-5);
    }

    // Seam accumulation. matchDof (gsSquareDomain::_initMapper, triggered by
    // addInterface()+applyOptions() -- see gsSquareDomain_test.cpp's
    // CopyCtorPreservesInterfacesAcrossApplyOptions, whose west/east pairing
    // this mirrors) can map TWO active sigma coefficients to the SAME free
    // index; addCoefficientGrad()'s final scatter (gsFoldBarrier.hpp) must
    // accumulate (+=) there, not overwrite (=), or one coefficient's
    // contribution is silently dropped -- the same hazard
    // gsSquareDomain::detJacobianDeriv_into's own scatter documents.
    // m_gX/m_gY (what addCoefficientGrad() scatters) are produced via a
    // reshaped contraction indexed by the flat tensor index k, not a direct
    // nb-row dense product, so the scatter's correctness at a seam is worth
    // pinning directly rather than assumed. Uses D5/D6b's asymmetric basis
    // (degU=3,nelU=3 / degV=2,nelV=2) so the interface case is non-square too.
    TEST(D8_FoldBarrier_SeamAccumulation_GradientMatchesFD)
    {
        gsKnotVector<real_t> kvU(0, 1, 2, 4); // degree 3, 3 elements
        gsKnotVector<real_t> kvV(0, 1, 1, 3); // degree 2, 2 elements
        gsTensorBSplineBasis<2,real_t> basis(kvU, kvV);
        gsSquareDomain<real_t> domain(basis, true);

        const size_t nBefore = domain.nControls();
        domain.addInterface(boundaryInterface(patchSide(0, boundary::west),
                                               patchSide(0, boundary::east), true));
        domain.applyOptions();
        const size_t nWithIfc = domain.nControls();
        gsTestInfo << "D8 nBefore=" << nBefore << " nWithIfc=" << nWithIfc << "\n";
        CHECK(nWithIfc < nBefore); // confirms matchDof actually coupled two dofs

        gsVector<real_t> c = domain.getControls();
        c.array() += 0.05;
        domain.setControls(c);

        gsGeometry<real_t>::uPtr S = makePlanarS(basis);
        gsVector<real_t> a(2), b(2); a.setZero(); b.setOnes();
        gsVector<unsigned> np(2); np << 6, 6;
        gsMatrix<real_t> uv = gsPointGrid(a, b, np);
        gsMatrix<real_t> xyz;
        S->eval_into(uv, xyz);
        for (index_t i = 0; i != xyz.cols(); i++)
            xyz(1, i) += 0.1 * math::sin(3.0 * EIGEN_PI * uv(0, i));

        gsOptFit<real_t> opt(domain, *S, uv, xyz, /*mu*/10.0, /*eps*/2.0,
                             gsFoldBarrierMode::Coefficient, /*quB*/0);
        real_t err = gsCheckSigmaGradient(opt, domain);
        gsTestInfo << "D8 seam checkGradient = " << err << "\n";
        CHECK(err < 1e-5);
    }

    TEST(D7_CoefficientBoundIsBelowSampledMinimum)
    {
        gsKnotVector<real_t> kvU(0, 1, 2, 4); // degree 3, 3 elements
        gsKnotVector<real_t> kvV(0, 1, 1, 3); // degree 2, 2 elements
        gsTensorBSplineBasis<2,real_t> basis(kvU, kvV);
        gsSquareDomain<real_t> domain(basis, true);

        gsVector<real_t> c = domain.getControls();
        c.array() += 0.05;
        domain.setControls(c);

        gsFoldBarrier<real_t> barrier(domain, /*mu*/1.0, /*eps*/0.0,
                                      gsFoldBarrierMode::Coefficient, /*quB*/0);
        const real_t coeffMin   = barrier.minDetJacobian();
        const real_t sampledMin = domain.minJacobian(7);

        gsTestInfo << "D7 coeffMin=" << coeffMin << " sampledMin=" << sampledMin << "\n";
        CHECK(coeffMin <= sampledMin + 1e-12);
    }

    // -----------------------------------------------------------------
    // The hierarchical arm of makeIntegrationBasis (built on
    // gsHTensorBasis::merge) and its hard nested-only precondition: sigma's
    // mesh must be a dyadic level of the analysis hierarchy, sharing its
    // level 0. Tests pin: (F1) on an unrefined THB basis the hierarchical
    // and tensor arms agree exactly (same partition, same degree, same
    // gsOptMesh objective) -- the load-bearing check that the new arm did
    // not silently change the quadrature; (F2) on a locally refined THB
    // basis the hierarchical arm strictly undercounts the tensor arm by a
    // derived amount; (F4) a non-nested sigma mesh throws (GISMO_ENSURE,
    // release-active), while the tensor overload on the same pair keeps
    // working, with sigma's knot line present in its result.
    // -----------------------------------------------------------------

    typedef gsAdaptiveParametrization<real_t,MonitorMode::ValueBased> APV;

    TEST(F1a_IntegrationBasis_HierarchicalMatchesTensor_Uniform_Sigma3)
    {
        gsKnotVector<real_t> kv(0, 1, 3, 3);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsTHBSplineBasis<2,real_t> thbU(tb);
        CHECK_EQUAL(index_t(0), thbU.maxLevel());
        CHECK_EQUAL(size_t(16), thbU.numElements());

        gsSquareDomain<real_t> sigma3(makeSigmaBasis(2, 3), true);
        gsTensorBSplineBasis<2,real_t> sigmaBasis3 =
            static_cast<const gsTensorBSplineBasis<2,real_t>&>(sigma3.domain().basis());

        gsBasis<real_t>::uPtr hier = APV::makeIntegrationBasis<2>(thbU, sigmaBasis3);
        gsTensorBSplineBasis<2,real_t> tens = APV::makeIntegrationBasis<2>(tb, sigmaBasis3);

        gsTestInfo << "F1a hier elements=" << hier->numElements()
                   << " tensor elements=" << tens.numElements() << "\n";
        CHECK_EQUAL(size_t(64), hier->numElements());
        CHECK_EQUAL(size_t(64), tens.numElements());
        CHECK_EQUAL(4, hier->degree(0)); CHECK_EQUAL(4, hier->degree(1));
        CHECK_EQUAL(4, tens.degree(0));  CHECK_EQUAL(4, tens.degree(1));
        CHECK(sameElements(*hier, tens));
    }

    TEST(F1b_IntegrationBasis_HierarchicalMatchesTensor_Uniform_Sigma2)
    {
        gsKnotVector<real_t> kv(0, 1, 3, 3);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsTHBSplineBasis<2,real_t> thbU(tb);

        gsSquareDomain<real_t> sigma2(makeSigmaBasis(2, 2), true);
        gsTensorBSplineBasis<2,real_t> sigmaBasis2 =
            static_cast<const gsTensorBSplineBasis<2,real_t>&>(sigma2.domain().basis());

        gsBasis<real_t>::uPtr hier = APV::makeIntegrationBasis<2>(thbU, sigmaBasis2);
        gsTensorBSplineBasis<2,real_t> tens = APV::makeIntegrationBasis<2>(tb, sigmaBasis2);

        gsTestInfo << "F1b hier elements=" << hier->numElements()
                   << " tensor elements=" << tens.numElements() << "\n";
        CHECK_EQUAL(size_t(16), hier->numElements());
        CHECK_EQUAL(size_t(16), tens.numElements());
        CHECK_EQUAL(4, hier->degree(0)); CHECK_EQUAL(4, hier->degree(1));
        CHECK_EQUAL(4, tens.degree(0));  CHECK_EQUAL(4, tens.degree(1));
        CHECK(sameElements(*hier, tens));
    }

    TEST(F1c_IntegrationBasis_HierarchicalAndTensorArmsAgree)
    {
        gsKnotVector<real_t> kv(0, 1, 3, 3);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);

        gsSquareDomain<real_t> sigma3(makeSigmaBasis(2, 3), true);
        gsTensorBSplineBasis<2,real_t> sigmaBasis3 =
            static_cast<const gsTensorBSplineBasis<2,real_t>&>(sigma3.domain().basis());
        gsGeometry<real_t>::uPtr S = makePlanarS(tb);

        const gsVector<real_t> ctrlsIdentity = sigma3.getControls();
        gsVector<real_t> ctrlsPerturbed = ctrlsIdentity;
        ctrlsPerturbed.array() += 0.02;

        // --- Fixture A: uniform THB (no refinement). Both arms produce the
        // SAME 64-element partition here, so this only proves identical
        // meshes give identical energies -- fixture B below is where the
        // objective-agreement claim actually gets exercised, under a
        // genuine mesh mismatch between the two arms.
        gsTHBSplineBasis<2,real_t> thbU(tb);
        gsBasis<real_t>::uPtr ibHierarchical = APV::makeIntegrationBasis<2>(thbU, sigmaBasis3);
        gsTensorBSplineBasis<2,real_t> ibTensor = APV::makeIntegrationBasis<2>(tb, sigmaBasis3);
        CHECK(sameElements(*ibHierarchical, ibTensor));
        CHECK_EQUAL(ibHierarchical->degree(0), ibTensor.degree(0));
        CHECK_EQUAL(ibHierarchical->degree(1), ibTensor.degree(1));

        gsOptMesh<real_t,MonitorMode::ValueBased> objHier(sigma3, *S, ibHierarchical.get());
        gsOptMesh<real_t,MonitorMode::ValueBased> objTensor(sigma3, *S, &ibTensor);
        objHier.options().setReal("quA", 1.0); objHier.options().setInt("quB", 1);
        objTensor.options().setReal("quA", 1.0); objTensor.options().setInt("quB", 1);

        sigma3.setControls(ctrlsIdentity);
        real_t energyHierIdentity   = objHier.evalObj(gsAsConstVector<real_t>(ctrlsIdentity.data(), ctrlsIdentity.rows()));
        real_t energyTensorIdentity = objTensor.evalObj(gsAsConstVector<real_t>(ctrlsIdentity.data(), ctrlsIdentity.rows()));
        gsTestInfo << "F1c uniform identity energyHier=" << energyHierIdentity
                   << " energyTensor=" << energyTensorIdentity << "\n";
        CHECK(math::abs(energyTensorIdentity) > 1e-12);
        CHECK(math::abs(energyHierIdentity - energyTensorIdentity) / math::abs(energyTensorIdentity) < 1e-12);

        sigma3.setControls(ctrlsPerturbed);
        real_t energyHierPerturbed   = objHier.evalObj(gsAsConstVector<real_t>(ctrlsPerturbed.data(), ctrlsPerturbed.rows()));
        real_t energyTensorPerturbed = objTensor.evalObj(gsAsConstVector<real_t>(ctrlsPerturbed.data(), ctrlsPerturbed.rows()));
        gsTestInfo << "F1c uniform perturbed energyHier=" << energyHierPerturbed
                   << " energyTensor=" << energyTensorPerturbed << "\n";
        CHECK(math::abs(energyTensorPerturbed) > 1e-12);
        CHECK(math::abs(energyHierPerturbed - energyTensorPerturbed) / math::abs(energyTensorPerturbed) < 1e-12);

        // --- Fixture B: locally refined THB, same construction as F2 below.
        // The hierarchical and tensor arms now produce DIFFERENT element
        // counts (different quadrature meshes) -- this is the case that
        // actually defends the nested-only assert in makeIntegrationBasis.
        //
        // S must be a SINGLE-ELEMENT (Bezier) geometry here, not the
        // multi-knot \a S used in fixture A above: under a non-identity
        // sigma, S's own knot lines pull back to curves that are element
        // boundaries of neither integration mesh, so the composed
        // integrand has a genuine C^1 kink there and Gauss quadrature over
        // it converges only algebraically -- no attainable quA/quB brings
        // the two meshes' quadrature error below roundoff (verified: the
        // multi-knot S plateaus at a ~1e-5 relative mismatch from quB=1 up
        // to quB=8). A single global polynomial patch has no such kink --
        // its only breaklines are sigma's own knots, which both
        // integration meshes honor by construction -- so integration
        // error, not partition choice, is the only source of mismatch,
        // and it vanishes well inside quadrature order margin.
        gsKnotVector<real_t> kvBezier(0, 1, 0, 3);
        gsTensorBSplineBasis<2,real_t> tbBezier(kvBezier, kvBezier);
        gsGeometry<real_t>::uPtr SBezier = makePlanarS(tbBezier);

        gsTHBSplineBasis<2,real_t> thbL(tb);
        thbL.refineElements({{1, 0, 0, 2, 2}});
        thbL.refineElements({{2, 0, 0, 2, 2}});

        gsBasis<real_t>::uPtr ibHierarchicalRefined = APV::makeIntegrationBasis<2>(thbL, sigmaBasis3);
        gsTensorBSplineBasis<2,real_t> ibTensorRefined = APV::makeIntegrationBasis<2>(thbL.tensorLevel(2), sigmaBasis3);
        gsTestInfo << "F1c refined element counts: hierarchical=" << ibHierarchicalRefined->numElements()
                   << " tensor=" << ibTensorRefined.numElements() << "\n";
        CHECK(ibHierarchicalRefined->numElements() != ibTensorRefined.numElements());

        gsOptMesh<real_t,MonitorMode::ValueBased> objHierRefined(sigma3, *SBezier, ibHierarchicalRefined.get());
        gsOptMesh<real_t,MonitorMode::ValueBased> objTensorRefined(sigma3, *SBezier, &ibTensorRefined);
        objHierRefined.options().setReal("quA", 1.0); objHierRefined.options().setInt("quB", 4);
        objTensorRefined.options().setReal("quA", 1.0); objTensorRefined.options().setInt("quB", 4);

        sigma3.setControls(ctrlsIdentity);
        real_t energyHierRefinedIdentity   = objHierRefined.evalObj(gsAsConstVector<real_t>(ctrlsIdentity.data(), ctrlsIdentity.rows()));
        real_t energyTensorRefinedIdentity = objTensorRefined.evalObj(gsAsConstVector<real_t>(ctrlsIdentity.data(), ctrlsIdentity.rows()));
        gsTestInfo << "F1c refined identity energyHier=" << energyHierRefinedIdentity
                   << " energyTensor=" << energyTensorRefinedIdentity << "\n";
        CHECK(math::abs(energyTensorRefinedIdentity) > 1e-12);
        CHECK(math::abs(energyHierRefinedIdentity - energyTensorRefinedIdentity) / math::abs(energyTensorRefinedIdentity) < 1e-12);

        sigma3.setControls(ctrlsPerturbed);
        real_t energyHierRefinedPerturbed   = objHierRefined.evalObj(gsAsConstVector<real_t>(ctrlsPerturbed.data(), ctrlsPerturbed.rows()));
        real_t energyTensorRefinedPerturbed = objTensorRefined.evalObj(gsAsConstVector<real_t>(ctrlsPerturbed.data(), ctrlsPerturbed.rows()));
        gsTestInfo << "F1c refined perturbed energyHier=" << energyHierRefinedPerturbed
                   << " energyTensor=" << energyTensorRefinedPerturbed << "\n";
        CHECK(math::abs(energyTensorRefinedPerturbed) > 1e-12);
        CHECK(math::abs(energyHierRefinedPerturbed - energyTensorRefinedPerturbed) / math::abs(energyTensorRefinedPerturbed) < 1e-12);
    }

    TEST(F2_IntegrationBasis_HierarchicalUndercountsTensor_LocallyRefined)
    {
        gsKnotVector<real_t> kv(0, 1, 3, 3);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsTHBSplineBasis<2,real_t> thbL(tb);
        thbL.refineElements({{1, 0, 0, 2, 2}});  // level-1 box [0,2]^2 == parameter [0,0.25]^2
        thbL.refineElements({{2, 0, 0, 2, 2}});  // level-2 box [0,2]^2 == parameter [0,0.125]^2
        CHECK_EQUAL(index_t(2), thbL.maxLevel());
        CHECK_EQUAL(size_t(22), thbL.numElements());
        CHECK_EQUAL(size_t(256), thbL.tensorLevel(2).numElements());

        gsSquareDomain<real_t> sigma3(makeSigmaBasis(2, 3), true);
        gsTensorBSplineBasis<2,real_t> sigmaBasis3 =
            static_cast<const gsTensorBSplineBasis<2,real_t>&>(sigma3.domain().basis());

        gsBasis<real_t>::uPtr hier = APV::makeIntegrationBasis<2>(thbL, sigmaBasis3);
        gsTensorBSplineBasis<2,real_t> tens = APV::makeIntegrationBasis<2>(thbL.tensorLevel(2), sigmaBasis3);

        gsTestInfo << "F2 nHier=" << hier->numElements()
                   << " nTensor=" << tens.numElements() << "\n";
        CHECK_EQUAL(size_t(67), hier->numElements());
        CHECK_EQUAL(size_t(256), tens.numElements());
        CHECK(hier->numElements() < tens.numElements());
    }

    TEST(F4_IntegrationBasis_NonNestedSigma_Throws)
    {
        gsKnotVector<real_t> kv(0, 1, 3, 3);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsTHBSplineBasis<2,real_t> thbU(tb);
        CHECK_EQUAL(index_t(0), thbU.maxLevel());

        gsKnotVector<real_t> kvS({0,0,0,0.3,1,1,1}, 2);   // degree 2, break at 0.3
        gsTensorBSplineBasis<2,real_t> sigmaBasisNN(kvS, kvS);
        gsSquareDomain<real_t> sigmaNN(sigmaBasisNN);

        CHECK_THROW(APV::makeIntegrationBasis<2>(thbU, sigmaBasisNN), std::runtime_error);
        gsTensorBSplineBasis<2,real_t> tens = APV::makeIntegrationBasis<2>(thbU.tensorLevel(0), sigmaBasisNN);

        gsTestInfo << "F4 tens elements=" << tens.numElements()
                   << " degree=" << tens.degree(0) << "x" << tens.degree(1) << "\n";
        CHECK_EQUAL(size_t(25), tens.numElements());
        CHECK(hasBreak(tens, 0, 0.3));
        CHECK(hasBreak(tens, 1, 0.3));
        CHECK_EQUAL(4, tens.degree(0));
        CHECK_EQUAL(4, tens.degree(1));
    }

    // The pass-through (integrationBasisIsFinal_t) constructor must store the
    // handed-in basis VERBATIM -- no knot union, no degree
    // raise. The ordinary constructors raise the degree to
    // p_analysis*p_sigma (see makeIntegrationBasis above); feeding an
    // ALREADY-raised basis through one of those would raise it a second
    // time, which is exactly what this constructor exists to avoid.
    TEST(G1_PassThroughCtor_NoDoubleElevation)
    {
        TestFixture f;

        gsTensorBSplineBasis<2,real_t> sigmaBasis =
            static_cast<const gsTensorBSplineBasis<2,real_t>&>(f.domain.domain().basis());
        // Union basis built EXACTLY as a driver would (knots inserted, degree
        // raised) BEFORE handing it to the pass-through constructor.
        gsTensorBSplineBasis<2,real_t> ibasis = APV::makeIntegrationBasis(f.tbasis, sigmaBasis);

        // Local probe: exposes the protected m_integrationBasis so the test
        // can inspect what the constructor actually stored (mirrors the
        // gsAdaptiveParametrizationProbe pattern in monitor_gradcheck.cpp).
        struct Probe : public APV
        {
            using APV::APV;
            const gsBasis<real_t> & integrationBasis() const { return *this->m_integrationBasis; }
        };

        gsGradientDescent<real_t> optimizer; // required by the ctor only; never solve()d
        Probe probe(f.domain, *f.composition, nullptr, ibasis, optimizer,
                    /*parametric*/true, integrationBasisIsFinal);

        gsTestInfo << "G1 handed-in degree=" << ibasis.degree(0) << "x" << ibasis.degree(1)
                   << " elements=" << ibasis.numElements()
                   << " | stored degree=" << probe.integrationBasis().degree(0) << "x"
                   << probe.integrationBasis().degree(1)
                   << " elements=" << probe.integrationBasis().numElements() << "\n";

        CHECK_EQUAL(ibasis.degree(0), probe.integrationBasis().degree(0));
        CHECK_EQUAL(ibasis.degree(1), probe.integrationBasis().degree(1));
        CHECK_EQUAL(ibasis.numElements(), probe.integrationBasis().numElements());
    }

    // -----------------------------------------------------------------
    // The pass-through (integrationBasisIsFinal_t) constructor's domainDim +
    // knot-containment GISMO_ENSUREs
    // (gsAdaptiveParametrization.hpp:2181-2271), including the widening
    // from a narrower gsTHBSplineBasis<d,T> cast to gsHTensorBasis<d,T> that
    // makes gsHBSplineBasis (== gsTHBSplineBasis<d,T,false>, the SAME class
    // template as gsTHBSplineBasis<d,T,true> with a different non-type
    // parameter, NOT a derived class -- gsForwardDeclarations.h:172)
    // actually checked.
    //
    // Every reject/accept case below builds its own sigma via
    // makeSigmaBasis(2,1) (one interior knot at 0.5) rather than using
    // TestFixture::domain: the fixture's BSplineSquareDeg(2) knot vector is
    // {0,0,0,1,1,1} -- NO interior knots -- against which the containment
    // check is vacuously true, proving nothing for either the reject or the
    // accept case. sigmaBasis.numElements() > 1 pins that this sigma is
    // never accidentally re-vacuated.
    // -----------------------------------------------------------------

    TEST(G2a_PassThroughCtor_RejectsWrongDomainDim)
    {
        TestFixture f;
        gsSquareDomain<real_t> sigma(makeSigmaBasis(2, 1), true);
        gsTensorBSplineBasis<2,real_t> sigmaBasis =
            static_cast<const gsTensorBSplineBasis<2,real_t>&>(sigma.domain().basis());
        CHECK(sigmaBasis.numElements() > 1);

        gsKnotVector<real_t> kv1D(0, 1, 0, 2);
        gsBSplineBasis<real_t> wrongDim(kv1D);   // domainDim() == 1 != composition's 2

        gsGradientDescent<real_t> optimizer; // required by the ctor only; never solve()d
        CHECK_THROW(APV(sigma, *f.composition, nullptr, wrongDim, optimizer,
                        /*parametric*/true, integrationBasisIsFinal),
                    std::runtime_error);
    }

    TEST(G2b_PassThroughCtor_RejectsCoarseTensorBasis)
    {
        TestFixture f;
        gsSquareDomain<real_t> sigma(makeSigmaBasis(2, 1), true);
        gsTensorBSplineBasis<2,real_t> sigmaBasis =
            static_cast<const gsTensorBSplineBasis<2,real_t>&>(sigma.domain().basis());
        CHECK(sigmaBasis.numElements() > 1);

        gsKnotVector<real_t> kvCoarse({0,0,0,1,1,1}, 2);   // no interior knots
        gsTensorBSplineBasis<2,real_t> coarseTensor(kvCoarse, kvCoarse);

        gsGradientDescent<real_t> optimizer;
        CHECK_THROW(APV(sigma, *f.composition, nullptr, coarseTensor, optimizer,
                        /*parametric*/true, integrationBasisIsFinal),
                    std::runtime_error);
    }

    TEST(G2c_PassThroughCtor_RejectsCoarseHBSplineBasis)
    {
        // The regression for the gsHTensorBasis widening: a narrower
        // gsTHBSplineBasis<2,T> cast at the constructor's call site would
        // silently miss gsHBSplineBasis and let this through unchecked.
        TestFixture f;
        gsSquareDomain<real_t> sigma(makeSigmaBasis(2, 1), true);
        gsTensorBSplineBasis<2,real_t> sigmaBasis =
            static_cast<const gsTensorBSplineBasis<2,real_t>&>(sigma.domain().basis());
        CHECK(sigmaBasis.numElements() > 1);

        gsKnotVector<real_t> kvCoarse({0,0,0,1,1,1}, 2);
        gsTensorBSplineBasis<2,real_t> coarseTensor(kvCoarse, kvCoarse);
        gsHBSplineBasis<2,real_t> coarseHB(coarseTensor);   // unrefined

        gsGradientDescent<real_t> optimizer;
        CHECK_THROW(APV(sigma, *f.composition, nullptr, coarseHB, optimizer,
                        /*parametric*/true, integrationBasisIsFinal),
                    std::runtime_error);
    }

    TEST(G2d_PassThroughCtor_RejectsCoarseTHBSplineBasis)
    {
        TestFixture f;
        gsSquareDomain<real_t> sigma(makeSigmaBasis(2, 1), true);
        gsTensorBSplineBasis<2,real_t> sigmaBasis =
            static_cast<const gsTensorBSplineBasis<2,real_t>&>(sigma.domain().basis());
        CHECK(sigmaBasis.numElements() > 1);

        gsKnotVector<real_t> kvCoarse({0,0,0,1,1,1}, 2);
        gsTensorBSplineBasis<2,real_t> coarseTensor(kvCoarse, kvCoarse);
        gsTHBSplineBasis<2,real_t> coarseTHB(coarseTensor);   // unrefined

        gsGradientDescent<real_t> optimizer;
        CHECK_THROW(APV(sigma, *f.composition, nullptr, coarseTHB, optimizer,
                        /*parametric*/true, integrationBasisIsFinal),
                    std::runtime_error);
    }

    TEST(G2e_PassThroughCtor_AcceptsMatchingBasis)
    {
        TestFixture f;
        gsSquareDomain<real_t> sigma(makeSigmaBasis(2, 1), true);
        gsTensorBSplineBasis<2,real_t> sigmaBasis =
            static_cast<const gsTensorBSplineBasis<2,real_t>&>(sigma.domain().basis());
        CHECK(sigmaBasis.numElements() > 1);

        // The same route a real driver uses to build a final integration
        // basis, so this exercises the real accepted configuration rather
        // than a contrivance.
        gsTensorBSplineBasis<2,real_t> matching = APV::makeIntegrationBasis(f.tbasis, sigmaBasis);

        gsGradientDescent<real_t> optimizer;
        APV probe(sigma, *f.composition, nullptr, matching, optimizer,
                  /*parametric*/true, integrationBasisIsFinal);
        CHECK(true);   // reaching here means construction did not throw
    }

    // -----------------------------------------------------------------
    // Guard tests for gsOptMesh: quadrature-node options (checkQuadOptions,
    // gsAdaptiveParametrization.hpp), Penalty>0, the ValueBased
    // 1+theta*f>0 GISMO_ENSURE (post-sweep), and the constructor's
    // det(J_g)>0 orientation sample.
    // -----------------------------------------------------------------

    TEST(Q1_QuadOptions_NegativeQuBThrows)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &f.tbasis);
        opt.options().setInt("quB", -1);
        CHECK_THROW(opt.evalObj(gsAsConstVector<real_t>(f.controls.data(), f.controls.rows())),
                    std::runtime_error);
    }

    TEST(Q2_QuadOptions_ZeroNodesThrows)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &f.tbasis);
        opt.options().setReal("quA", 0.1);
        opt.options().setInt("quB", 0);
        CHECK_THROW(opt.evalObj(gsAsConstVector<real_t>(f.controls.data(), f.controls.rows())),
                    std::runtime_error);
    }

    TEST(Q3_QuadOptions_SinglePointRuleAccepted)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &f.tbasis);
        opt.options().setReal("quA", 0.0);
        opt.options().setInt("quB", 1);
        real_t val = opt.evalObj(gsAsConstVector<real_t>(f.controls.data(), f.controls.rows()));
        gsTestInfo << "Q3 val=" << val << "\n";
        CHECK(std::isfinite(val));
    }

    TEST(P1_PenaltyZeroThrows)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &f.tbasis);
        opt.options().setReal("Penalty", 0.0);
        CHECK_THROW(opt.evalObj(gsAsConstVector<real_t>(f.controls.data(), f.controls.rows())),
                    std::runtime_error);
    }

    TEST(T1_SmoothingNegativeThrows)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        opt.options().setReal("Smoothing", -1.0);
        CHECK_THROW(opt.evalObj(gsAsConstVector<real_t>(f.controls.data(), f.controls.rows())),
                    std::runtime_error);
    }

    TEST(V1_ValueBased_NegativeDenominatorThrows)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("0*x - 5", 2);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        opt.options().setReal("Smoothing", 1.0); // 1+theta*f = 1+1*(-5) = -4 < 0
        CHECK_THROW(opt.evalObj(gsAsConstVector<real_t>(f.controls.data(), f.controls.rows())),
                    std::runtime_error);
    }

    TEST(V2_GradientBased_ConstantMonitorNoThrow)
    {
        TestFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("0*x - 5", 2);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        opt.options().setReal("Smoothing", 1.0); // grad(const) = 0, so 1+theta*||grad f||^2 = 1 > 0
        real_t val = opt.evalObj(gsAsConstVector<real_t>(f.controls.data(), f.controls.rows()));
        gsTestInfo << "V2 val=" << val << "\n";
        CHECK(std::isfinite(val));
    }

    TEST(O1_MirroredSquare_NegativeOrientationThrows)
    {
        TestFixture f;
        gsGeometry<>::uPtr sq = gsNurbsCreator<>::BSplineSquare(real_t(1));
        gsMatrix<real_t> coefs = sq->coefs();
        coefs.col(0).swap(coefs.col(1)); // mirror: reflection, det(J_g) < 0 everywhere
        gsGeometry<>::uPtr mirrored = sq->basis().makeGeometry(coefs);
        CHECK_THROW(
            (gsOptMesh<real_t,MonitorMode::ValueBased>(f.domain, *mirrored, &f.tbasis)),
            std::runtime_error);
    }

    TEST(O2_UnmirroredSquare_NoThrow)
    {
        TestFixture f;
        gsGeometry<>::uPtr sq = gsNurbsCreator<>::BSplineSquare(real_t(1));
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *sq, &f.tbasis);
        real_t val = opt.evalObj(gsAsConstVector<real_t>(f.controls.data(), f.controls.rows()));
        gsTestInfo << "O2 val=" << val << "\n";
        CHECK(std::isfinite(val));
    }

    TEST(G3_TensorBasisAdmissible_ContainedAndMissing)
    {
        // Exercises tensorBasisAdmissible() (now an internal helper of
        // gsAdaptiveParametrization.hpp) only through the public pass-through
        // (integrationBasisIsFinal_t) constructor, which is its sole caller.
        TestFixture f;
        gsSquareDomain<real_t> sigma(makeSigmaBasis(2, 1), true);
        gsTensorBSplineBasis<2,real_t> sigmaBasis =
            static_cast<const gsTensorBSplineBasis<2,real_t>&>(sigma.domain().basis());
        CHECK(sigmaBasis.numElements() > 1);

        gsGradientDescent<real_t> optimizer;

        // Contains all of sigma's interior knots: accepted.
        gsTensorBSplineBasis<2,real_t> containing = APV::makeIntegrationBasis(f.tbasis, sigmaBasis);
        APV probe(sigma, *f.composition, nullptr, containing, optimizer,
                  /*parametric*/true, integrationBasisIsFinal);
        CHECK(true);   // reaching here means construction did not throw

        // Missing sigma's interior knot in direction 0 only: rejected.
        gsKnotVector<real_t> kvNoInterior({0,0,0,1,1,1}, 2);
        gsTensorBSplineBasis<2,real_t> missingU(kvNoInterior, sigmaBasis.knots(1));
        CHECK_THROW(APV(sigma, *f.composition, nullptr, missingU, optimizer,
                        /*parametric*/true, integrationBasisIsFinal),
                    std::runtime_error);

        // Fine in u, coarse in v: missing sigma's interior knot in
        // direction 1 only: rejected.
        gsTensorBSplineBasis<2,real_t> missingV(sigmaBasis.knots(0), kvNoInterior);
        CHECK_THROW(APV(sigma, *f.composition, nullptr, missingV, optimizer,
                        /*parametric*/true, integrationBasisIsFinal),
                    std::runtime_error);
    }

}
