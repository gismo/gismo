/** @file gsAdaptiveParametrizationNewton_test.cpp

    @brief Tests for gsAdaptiveParametrizationNewton: independent residual
    assembly (vs gradObj_into), analytic Hessian (vs FD of the analytic
    gradient + symmetry), Picard SPD property, and Newton convergence order.

    See doc/derivation_hessian.md in the sibling gismo_opt checkout
    (../gismo_opt/doc/derivation_hessian.md) for the assembled formulas --
    that document is not part of this repository.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
 **/

#include "gismo_unittest.h"
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsAssembler/gsAdaptiveParametrizationNewton.h>

#include <sstream>

#ifdef GISMO_WITH_ADIFF
static const real_t RES_TOL  = 1e-9;   // residual vs gradObj_into (same quadrature => machine precision)
static const real_t HESS_TOL = 5e-6;   // analytic K vs central-FD of analytic gradient (h=1e-7)
#else
static const real_t RES_TOL  = 1e-9;   // independent of f-derivative quality (same f calls)
static const real_t HESS_TOL = 1e-2;   // gsFunctionExpr derivs are FD without ADIFF
#endif
// H7/H8/H9 use AnalyticCubicMonitor (closed-form deriv2_into), so their FD
// comparison is genuinely accurate in BOTH the ADIFF and non-ADIFF build --
// unlike HESS_TOL above, this tolerance does not need an #ifdef. Measured
// achieved err ~ O(1e-9); 1e-6 leaves ~3 orders of margin.
static const real_t HESS_TOL_ANALYTIC = 1e-6;
// Relative tolerance for the E1-E4 energy-consistency regression:
// _evalObjAndMinJ (hand-rolled sweep) vs gsOptMesh::evalObj (gsExprEvaluator),
// both fed the SAME controls and the SAME pinned quadrature (quA=1, quB=1).
// A stale Tikhonov shift in this energy assembly produced a relative energy
// mismatch of ~5e-2; 1e-12 catches that by ten orders while still tolerating
// OpenMP sum-reduction order differences between the two code paths.
static const real_t ENERGY_TOL = 1e-12;
// minJ is a min-reduction over the same elements with the same rule -- no
// summation-order ambiguity, so it is safe to pin far tighter than ENERGY_TOL.
static const real_t MINJ_TOL = 1e-14;

namespace {

struct NewtonFixture
{
    gsGeometry<>::uPtr composition;
    gsSquareDomain<real_t> domain;
    gsVector<real_t> controls;
    gsKnotVector<> kv1;
    gsKnotVector<> kv2;
    gsTensorBSplineBasis<2> tbasis;
    gsComposedBasis<> cbasis;
    gsMatrix<> anchors;

    NewtonFixture()
        : composition(gsNurbsCreator<>::BSplineSquareDeg(2))
        , domain(*composition)
        , controls(domain.nControls())
        , kv1({0,0,0,1,1,1}, 2)
        , kv2({0,0,0,0.50,1,1,1}, 2)
        , tbasis(kv1, kv2)
        , cbasis(*composition, tbasis)
        , anchors(cbasis.anchors().transpose())
    {
        controls << 0.95, 0.95;   // distorted sigma (away from a stationary point)
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

// Syncs the gsOptMesh oracle's options onto the Newton's, for a fair
// comparison between the two independent code paths.
//
// Penalty/Smoothing are genuinely shared options on both sides (gsOptMesh's
// Smoothing default is 0.1 against the Newton's 1e-2, hence the sync).
// quA/quB are NOT: gsAdaptiveParametrizationNewton<T,MODE>::defaultOptions()
// (gsAdaptiveParametrizationNewton.hpp:86-100) declares no quA/quB at all --
// its energy sweeps hardcode the rule instead, in _evalObjAndMinJ
// (gsAdaptiveParametrizationNewton.hpp:175-177):
//   quadOptions.addReal("quA","",1.0);
//   quadOptions.addInt ("quB","",1);
// gsOptMesh::evalObj/gradObj_into/computeMinJacobian, in contrast, read
// quA/quB from their own option list (gsAdaptiveParametrization.hpp:215-223,
// 419-420, 827-828, 1427-1429), whose defaults (1.0, 1) happen to coincide
// with the Newton's hardcode today. Pin them explicitly here rather than
// relying on that coincidence, and guard the assumption with GISMO_ENSURE
// (not GISMO_ASSERT: build_rel is -DNDEBUG and GISMO_ASSERT compiles out)
// so a future divergence is loud rather than silent.
template <MonitorMode MODE>
void syncEnergyOptions(gsAdaptiveParametrizationNewton<real_t,MODE> & newton,
                       gsOptMesh<real_t,MODE> & opt)
{
    opt.options().setReal("Penalty",   newton.options().getReal("Penalty"));
    opt.options().setReal("Smoothing", newton.options().getReal("Smoothing"));

    // The Newton's list has no quA/quB (defaultOptions,
    // gsAdaptiveParametrizationNewton.hpp:86-100); its sweeps hardcode
    // (1.0, 1). If that ever becomes an option, sync it here instead.
    // gsOptionList::exists()/isReal()/isInt() are private (gsOptionList.h:
    // 238-249), so the public "ask with an impossible sentinel default"
    // idiom is used to detect presence: quA/quB are always non-negative
    // quadrature-rule parameters, so a returned negative value means the
    // label was absent and askReal/askInt fell back to the sentinel.
    GISMO_ENSURE(newton.options().askReal("quA", -1.0) < 0.0 &&
                 newton.options().askInt ("quB", -1)   < 0,
                 "gsAdaptiveParametrizationNewton now exposes quA/quB: "
                 "sync them onto gsOptMesh instead of pinning (1.0, 1).");
    opt.options().setReal("quA", 1.0);
    opt.options().setInt ("quB", 1);
}

// Relative error of assembleResidual vs gradObj_into (same options)
template <MonitorMode MODE>
real_t residualRelError(gsAdaptiveParametrizationNewton<real_t,MODE> & newton,
                        gsOptMesh<real_t,MODE> & opt,
                        const gsVector<real_t> & u)
{
    syncEnergyOptions(newton, opt);

    gsVector<real_t> R;
    newton.assembleResidual(u, R);

    gsVector<real_t> g(u.size());
    gsAsVector<real_t> asg(g.data(), g.rows());
    opt.gradObj_into(gsAsConstVector<real_t>(u.data(), u.size()), asg);

    const real_t den = g.norm() > 1e-15 ? g.norm() : 1.0;
    return (R - g).norm() / den;
}

// Relative error of analytic Hessian vs central FD of the analytic gradient
template <MonitorMode MODE>
real_t hessianRelError(gsAdaptiveParametrizationNewton<real_t,MODE> & newton,
                       gsOptMesh<real_t,MODE> & opt,
                       const gsVector<real_t> & u,
                       real_t & symErr)
{
    syncEnergyOptions(newton, opt);

    gsSparseMatrix<real_t> K;
    newton.assembleHessian(u, K);

    gsMatrix<real_t> Kd = K.toDense();
    symErr = (Kd - Kd.transpose()).norm() / std::max(Kd.norm(), (real_t)1e-15);

    // FD of analytic gradient
    const index_t n = u.size();
    const real_t h = 1e-7;
    gsMatrix<real_t> Kfd(n, n);
    gsVector<real_t> up(n), um(n), uu = u;
    gsAsVector<real_t> aup(up.data(), n), aum(um.data(), n);
    for (index_t j = 0; j != n; ++j)
    {
        uu[j] = u[j] + h;
        opt.gradObj_into(gsAsConstVector<real_t>(uu.data(), n), aup);
        uu[j] = u[j] - h;
        opt.gradObj_into(gsAsConstVector<real_t>(uu.data(), n), aum);
        uu[j] = u[j];
        Kfd.col(j) = (up - um) / (2.0 * h);
    }
    Kfd = 0.5 * (Kfd + Kfd.transpose());

    return (Kd - Kfd).norm() / std::max(Kfd.norm(), (real_t)1e-15);
}

// Analytic-derivative monitor for H7/H8/H9.
//
// Why this exists: without GISMO_WITH_ADIFF, gsFunctionExpr<T>::deriv2_into
// (gsFunctionExpr.hpp:620-634) is itself numerical -- it evaluates exprtk's
// second_derivative/mixed_derivative at h=1e-5. gsAdaptiveParametrizationNewton's
// _funDeriv3_into (gsAdaptiveParametrizationNewton.hpp:473-506) then
// central-differences THAT result again at h=1e-5 to build d^3 f, so the third
// derivative carries round-off of order eps/h^3. That noise enters the
// assembled Hessian only through the theta * 2 v . dgf term of dms(d,m), which
// is why the observed H7-H9 asymmetry scales with the Smoothing option theta
// and would vanish only at theta = 0 -- it is exprtk differentiation noise,
// not a defect in assembleHessian. Supplying a monitor with CLOSED-FORM
// eval_into/deriv_into/deriv2_into removes that noise source entirely, so
// CHECK(symErr < 1e-8) becomes a genuine test of the Hessian FORMULA rather
// than of exprtk's finite-difference quality.
//
// The polynomial is capped at degree 3 (a cubic term per coordinate) so that
// deriv2_into is affine in the coordinates: _funDeriv3_into's outer central
// difference of deriv2 is then EXACT (no truncation error, only ~eps/h
// rounding), which is what pushes symErr down to ~1e-13..1e-17 instead of
// merely to exprtk-analytic (~1e-10) accuracy.
//
// f(x)   = c0 + a.x + 0.5 x^T B x + sum_i g_i x_i^3      (B symmetric)
// grad_i = a_i + (B x)_i + 3 g_i x_i^2
// H_ii   = B_ii + 6 g_i x_i
// H_kl   = B_kl   (k != l)
//
// deriv2_into row layout (verified against gsFunctionExpr.hpp:595,620-635):
// rows 0..d-1 are the pure second derivatives H_00..H_{d-1,d-1}, followed by
// the mixed ones in lexicographic (k,l), k<l. For d=2: f_xx,f_yy,f_xy.
// For d=3: f_xx,f_yy,f_zz,f_xy,f_xz,f_yz.
class AnalyticCubicMonitor : public gsFunction<real_t>
{
    short_t m_d;
    real_t m_c0;
    gsVector<real_t> m_a;
    gsMatrix<real_t> m_B;
    gsVector<real_t> m_g;
public:
    AnalyticCubicMonitor(real_t c0, const gsVector<real_t> & a,
                          const gsMatrix<real_t> & B, const gsVector<real_t> & g)
    : m_d(static_cast<short_t>(a.size())), m_c0(c0), m_a(a), m_B(B), m_g(g)
    { GISMO_ENSURE(B.rows()==m_d && B.cols()==m_d && g.size()==m_d,
                    "AnalyticCubicMonitor: inconsistent dimensions"); }

    short_t domainDim() const override { return m_d; }
    short_t targetDim() const override { return 1; }

    void eval_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(1, u.cols());
        for (index_t p = 0; p != u.cols(); ++p)
        {
            const gsVector<real_t> x = u.col(p);
            real_t val = m_c0 + m_a.dot(x) + 0.5 * x.dot(m_B * x);
            for (short_t i = 0; i != m_d; ++i) val += m_g[i] * x[i]*x[i]*x[i];
            result(0,p) = val;
        }
    }

    void deriv_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        result.resize(m_d, u.cols());
        for (index_t p = 0; p != u.cols(); ++p)
        {
            const gsVector<real_t> x = u.col(p);
            gsVector<real_t> grad = m_a + m_B * x;
            for (short_t i = 0; i != m_d; ++i) grad[i] += 3.0 * m_g[i] * x[i]*x[i];
            result.col(p) = grad;
        }
    }

    void deriv2_into(const gsMatrix<real_t> & u, gsMatrix<real_t> & result) const override
    {
        const short_t stride = m_d + m_d*(m_d-1)/2;
        result.resize(stride, u.cols());
        for (index_t p = 0; p != u.cols(); ++p)
        {
            const gsVector<real_t> x = u.col(p);
            for (short_t i = 0; i != m_d; ++i)
                result(i,p) = m_B(i,i) + 6.0 * m_g[i] * x[i];
            short_t row = m_d;
            for (short_t k = 0; k != m_d; ++k)
                for (short_t l = k+1; l != m_d; ++l)
                    result(row++,p) = m_B(k,l);
        }
    }
};

// Identity sigma basis at degree \a sigmaDeg with sigmaKnots = 2^sigmaRef - 1
// interior knots, in the exact knot-vector convention
// gsSquareDomain(const gsBasis<T>&, bool) expects. Copied from
// unittests/gsAdaptiveParametrizatin_test.cpp / gsSquareDomain_test.cpp
// (this file has no shared header with either).
gsTensorBSplineBasis<2,real_t> makeSigmaBasis(short_t sigmaDeg, index_t sigmaRef)
{
    gsKnotVector<real_t> ks(0, 1, (1 << sigmaRef) - 1, sigmaDeg + 1);
    return gsTensorBSplineBasis<2,real_t>(ks, ks);
}

// gsWarn (gsCore/gsDebug.h:50-55) is a macro expanding to std::cout<<"Warning: ",
// so a warning is only observable by capturing std::cout itself. The RAII
// restore in the destructor keeps the swap exception-safe: a throw out of
// solveNewton() still leaves std::cout usable afterward.
struct CoutCapture
{
    std::ostringstream oss;
    std::streambuf * saved;
    CoutCapture() : saved(std::cout.rdbuf(oss.rdbuf())) {}
    ~CoutCapture() { std::cout.rdbuf(saved); }
    std::string str() const { return oss.str(); }
};

} // anonymous namespace

SUITE(gsAdaptiveParametrizationNewton_test)
{

    // =====================================================================
    // Residual: independent weak-form assembly == gradObj_into  [H-3.2]
    // =====================================================================

    TEST(R1_Residual_Planar_NoMonitor)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, f.tbasis);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &f.tbasis);
        real_t err = residualRelError(newton, opt, f.controls);
        gsTestInfo << "R1 relErr = " << err << "\n";
        CHECK(err < RES_TOL);
    }

    TEST(R2_Residual_Surface_NoMonitor)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makeSurfaceGeom();
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, f.tbasis);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &f.tbasis);
        real_t err = residualRelError(newton, opt, f.controls);
        gsTestInfo << "R2 relErr = " << err << "\n";
        CHECK(err < RES_TOL);
    }

    TEST(R3_Residual_Planar_ValueBased)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, monFun, f.tbasis, true);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        real_t err = residualRelError(newton, opt, f.controls);
        gsTestInfo << "R3 relErr = " << err << "\n";
        CHECK(err < RES_TOL);
    }

    TEST(R4_Residual_Surface_ValueBased)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makeSurfaceGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, monFun, f.tbasis, true);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        real_t err = residualRelError(newton, opt, f.controls);
        gsTestInfo << "R4 relErr = " << err << "\n";
        CHECK(err < RES_TOL);
    }

    TEST(R5_Residual_Planar_GradientBased_Parametric)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::GradientBased> newton(f.domain, *geom, monFun, f.tbasis, true);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        real_t err = residualRelError(newton, opt, f.controls);
        gsTestInfo << "R5 relErr = " << err << "\n";
        CHECK(err < RES_TOL);
    }

    TEST(R6_Residual_Planar_GradientBased_Physical)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*x*x + 0.3*y*y", 2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::GradientBased> newton(f.domain, *geom, monFun, f.tbasis, false);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, false);
        real_t err = residualRelError(newton, opt, f.controls);
        gsTestInfo << "R6 relErr = " << err << "\n";
        CHECK(err < RES_TOL);
    }

    TEST(R7_Residual_Surface_GradientBased_Physical)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makeSurfaceGeom();
        gsFunctionExpr<> monFun("1 + 0.3*x + 0.2*y + 0.1*z", 3);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::GradientBased> newton(f.domain, *geom, monFun, f.tbasis, false);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, false);
        real_t err = residualRelError(newton, opt, f.controls);
        gsTestInfo << "R7 relErr = " << err << "\n";
        CHECK(err < RES_TOL);
    }

    // =====================================================================
    // Hessian: analytic K vs FD of analytic gradient + symmetry  [H-4]
    // =====================================================================

    TEST(H1_Hessian_Planar_NoMonitor)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, f.tbasis);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &f.tbasis);
        real_t symErr;
        real_t err = hessianRelError(newton, opt, f.controls, symErr);
        gsTestInfo << "H1 relErr = " << err << ", symErr = " << symErr << "\n";
        CHECK(symErr < 1e-12);
        CHECK(err < HESS_TOL);
    }

    TEST(H2_Hessian_Identity_NoMonitor)
    {
        // Identity geometry: geometric (third-derivative) block vanishes
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makeIdentityGeom();
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, f.tbasis);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &f.tbasis);
        real_t symErr;
        real_t err = hessianRelError(newton, opt, f.controls, symErr);
        gsTestInfo << "H2 relErr = " << err << ", symErr = " << symErr << "\n";
        CHECK(symErr < 1e-12);
        CHECK(err < HESS_TOL);
    }

    TEST(H3_Hessian_Surface_NoMonitor)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makeSurfaceGeom();
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, f.tbasis);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &f.tbasis);
        real_t symErr;
        real_t err = hessianRelError(newton, opt, f.controls, symErr);
        gsTestInfo << "H3 relErr = " << err << ", symErr = " << symErr << "\n";
        CHECK(symErr < 1e-12);
        CHECK(err < HESS_TOL);
    }

    TEST(H4_Hessian_Planar_ValueBased)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, monFun, f.tbasis, true);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        real_t symErr;
        real_t err = hessianRelError(newton, opt, f.controls, symErr);
        gsTestInfo << "H4 relErr = " << err << ", symErr = " << symErr << "\n";
        CHECK(symErr < 1e-10);
        CHECK(err < HESS_TOL);
    }

    TEST(H5_Hessian_Surface_ValueBased)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makeSurfaceGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, monFun, f.tbasis, true);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        real_t symErr;
        real_t err = hessianRelError(newton, opt, f.controls, symErr);
        gsTestInfo << "H5 relErr = " << err << ", symErr = " << symErr << "\n";
        CHECK(symErr < 1e-10);
        CHECK(err < HESS_TOL);
    }

    TEST(H6_Hessian_Planar_ValueBased_Physical)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*x*x + 0.3*y*y", 2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, monFun, f.tbasis, false);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &monFun, &f.tbasis, false);
        real_t symErr;
        real_t err = hessianRelError(newton, opt, f.controls, symErr);
        gsTestInfo << "H6 relErr = " << err << ", symErr = " << symErr << "\n";
        CHECK(symErr < 1e-10);
        CHECK(err < HESS_TOL);
    }

    TEST(H7_Hessian_Planar_GradientBased_Parametric)
    {
        // Analytic-derivative monitor (see AnalyticCubicMonitor above): removes
        // the exprtk FD-of-FD noise in d^3 f, so symErr tests the Hessian FORMULA.
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsVector<real_t> a2(2); a2 << 0.3, 0.2;
        gsMatrix<real_t> B2(2,2); B2 << 2,1, 1,4;
        gsVector<real_t> g2(2); g2 << 1.0, 1.0;
        AnalyticCubicMonitor monFun(1.0, a2, B2, g2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::GradientBased> newton(f.domain, *geom, monFun, f.tbasis, true);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        real_t symErr;
        real_t err = hessianRelError(newton, opt, f.controls, symErr);
        gsTestInfo << "H7 relErr = " << err << ", symErr = " << symErr << "\n";
        CHECK(symErr < 1e-8);
        CHECK(err < HESS_TOL_ANALYTIC);
    }

    TEST(H8_Hessian_Planar_GradientBased_Physical)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsVector<real_t> a2(2); a2 << 0.3, 0.2;
        gsMatrix<real_t> B2(2,2); B2 << 2,1, 1,4;
        gsVector<real_t> g2(2); g2 << 1.0, 1.0;
        AnalyticCubicMonitor monFun(1.0, a2, B2, g2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::GradientBased> newton(f.domain, *geom, monFun, f.tbasis, false);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, false);
        real_t symErr;
        real_t err = hessianRelError(newton, opt, f.controls, symErr);
        gsTestInfo << "H8 relErr = " << err << ", symErr = " << symErr << "\n";
        CHECK(symErr < 1e-8);
        CHECK(err < HESS_TOL_ANALYTIC);
    }

    TEST(H9_Hessian_Surface_GradientBased_Physical)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makeSurfaceGeom();
        gsVector<real_t> a3(3); a3 << 0.3, 0.2, 0.1;
        gsMatrix<real_t> B3(3,3); B3 << 2,1,0.5, 1,4,0.3, 0.5,0.3,1;
        gsVector<real_t> g3(3); g3 << 1.0, 1.0, 0.5;
        AnalyticCubicMonitor monFun(1.0, a3, B3, g3);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::GradientBased> newton(f.domain, *geom, monFun, f.tbasis, false);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, false);
        real_t symErr;
        real_t err = hessianRelError(newton, opt, f.controls, symErr);
        gsTestInfo << "H9 relErr = " << err << ", symErr = " << symErr << "\n";
        CHECK(symErr < 1e-8);
        CHECK(err < HESS_TOL_ANALYTIC);
    }

    // =====================================================================
    // Picard surrogate: symmetric positive definite  [H-5]
    // =====================================================================

    TEST(P1_Picard_SPD_NoMonitor)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, f.tbasis);
        gsSparseMatrix<real_t> B;
        newton.assemblePicard(f.controls, B);
        gsMatrix<real_t> Bd = B.toDense();
        const real_t symErr = (Bd - Bd.transpose()).norm() / Bd.norm();
        gsEigen::SelfAdjointEigenSolver<gsMatrix<real_t>::Base> eig(Bd);
        gsTestInfo << "P1 symErr = " << symErr
                   << ", min eig = " << eig.eigenvalues().minCoeff() << "\n";
        CHECK(symErr < 1e-12);
        CHECK(eig.eigenvalues().minCoeff() > 0.0);
    }

    TEST(P2_Picard_SPD_ValueBased)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, monFun, f.tbasis, true);
        gsSparseMatrix<real_t> B;
        newton.assemblePicard(f.controls, B);
        gsMatrix<real_t> Bd = B.toDense();
        const real_t symErr = (Bd - Bd.transpose()).norm() / Bd.norm();
        gsEigen::SelfAdjointEigenSolver<gsMatrix<real_t>::Base> eig(Bd);
        gsTestInfo << "P2 symErr = " << symErr
                   << ", min eig = " << eig.eigenvalues().minCoeff() << "\n";
        CHECK(symErr < 1e-12);
        CHECK(eig.eigenvalues().minCoeff() > 0.0);
    }

    // =====================================================================
    // Newton convergence: superlinear contraction (exact-Hessian signature)
    // =====================================================================

    TEST(N1_Newton_Converges_Planar_NoMonitor)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, f.tbasis);
        f.domain.setControls(f.controls);
        newton.options().setSwitch("Verbose", false);
        newton.options().setReal("Tolerance", 1e-10);
        newton.options().setInt("MaxIter", 30);

        const index_t iters = newton.solveNewton();
        const gsMatrix<real_t> & log = newton.iterLog();
        const real_t finalR = log(1, log.cols()-1);
        gsTestInfo << "N1 iters = " << iters << ", final ||R|| = " << finalR << "\n";
        CHECK(finalR < 1e-8);

        // Convergence-order estimate q = log(r_k/r_{k-1}) / log(r_{k-1}/r_{k-2})
        // over the last iterations above the roundoff floor; quadratic => q ~ 2
        if (log.cols() >= 4)
        {
            std::vector<real_t> rs;
            for (index_t i = 0; i != log.cols(); ++i)
                if (log(1, i) > 1e-13) rs.push_back(log(1, i));
            if (rs.size() >= 3)
            {
                const size_t k = rs.size() - 1;
                const real_t q = math::log(rs[k]/rs[k-1]) / math::log(rs[k-1]/rs[k-2]);
                gsTestInfo << "N1 convergence order q = " << q << "\n";
                CHECK(q > 1.5);
            }
        }
    }

    TEST(N2_Newton_vs_Picard_SameEnergy)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();

        // Newton
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, f.tbasis);
        f.domain.setControls(f.controls);
        newton.options().setSwitch("Verbose", false);
        newton.options().setReal("Tolerance", 1e-10);
        newton.solveNewton();
        const real_t E_newton = newton.evalObj();
        const real_t minJ_newton = newton.computeMinJacobian();

        // Picard (fresh object, same start)
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> picard(f.domain, *geom, f.tbasis);
        f.domain.setControls(f.controls);
        picard.options().setSwitch("Verbose", false);
        picard.options().setReal("Tolerance", 1e-8);
        picard.options().setInt("MaxIter", 200);
        picard.solvePicard();
        const real_t E_picard = picard.evalObj();
        const real_t minJ_picard = picard.computeMinJacobian();

        gsTestInfo << "N2 E_newton = " << E_newton << ", E_picard = " << E_picard
                   << ", minJ = " << minJ_newton << " / " << minJ_picard << "\n";
        CHECK(minJ_newton > 0.0);
        CHECK(minJ_picard > 0.0);
        CHECK(math::abs(E_newton - E_picard) < 1e-5 * math::abs(E_newton));
    }

    // =====================================================================
    // Energy consistency: private _evalObjAndMinJ (the sweep the Armijo
    // line search descends) vs the independent gsOptMesh::evalObj oracle,
    // at the DEFAULT Penalty = 1e-2 (not the Penalty -> 0 limit).
    //
    // _evalObjAndMinJ is private (gsAdaptiveParametrizationNewton.h:194) with
    // no public accessor; evalObj() delegates to m_optMesh.evalObj and so
    // does NOT exercise the independent assembly (it IS the oracle side).
    // Reached instead via iterLog(): _solveLoop (.hpp:1871-1891) computes
    // _evalObjAndMinJ(u, E, minJ) and writes it into columns (0,iter)/(2,iter)
    // BEFORE testing the stopping criterion, so with MaxIter=1 and a
    // Tolerance the residual can never beat, iterLog()(0,0)/(2,0) are exactly
    // _evalObjAndMinJ's outputs at the controls the solver started from --
    // no assembly, no line search in between.
    //
    // Note the asymmetry with the H*-suite checks above: symErr there
    // validates the Hessian assembly's internal consistency, but a
    // self-consistently-wrong formula is still symmetric -- corrupting
    // deriv2_into's row layout, for instance, does not move symErr. Only a
    // comparison against an INDEPENDENT oracle -- here gsOptMesh::evalObj, a
    // different code path (hand-rolled sweep vs gsExprEvaluator) -- can
    // catch a defect that a self-consistency check like symErr is blind to.
    // =====================================================================

    TEST(E1_Energy_FusedSweep_vs_EvalObj_Planar_NoMonitor)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, f.tbasis);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &f.tbasis);

        f.domain.setControls(f.controls);
        newton.options().setSwitch("Verbose", false);
        newton.options().setInt ("MaxIter", 1);
        newton.options().setReal("Tolerance", 1e30);   // log row 0, then stop:
                                                        // no assembly, no line search
        newton.solveNewton();
        const real_t E_newton = newton.iterLog()(0, 0);   // == _evalObjAndMinJ
        const real_t minJ_newton = newton.iterLog()(2, 0);

        // Oracle at the SAME controls (evalObj/computeMinJacobian call
        // setControls themselves, so order does not matter).
        syncEnergyOptions(newton, opt);
        gsAsConstVector<real_t> uc(f.controls.data(), f.controls.size());
        const real_t E_ref   = opt.evalObj(uc);
        const real_t minJ_ref = opt.computeMinJacobian(uc);

        gsTestInfo << "E1 E_newton = " << E_newton << ", E_ref = " << E_ref
                   << ", relErr(E) = " << math::abs(E_newton - E_ref) / math::abs(E_ref)
                   << ", minJ_newton = " << minJ_newton << ", minJ_ref = " << minJ_ref << "\n";
        CHECK(math::abs(E_newton - E_ref) < ENERGY_TOL * math::abs(E_ref));
        CHECK(math::abs(minJ_newton - minJ_ref) < MINJ_TOL * math::abs(minJ_ref));
    }

    TEST(E2_Energy_FusedSweep_vs_EvalObj_Planar_ValueBased)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, monFun, f.tbasis, true);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);

        f.domain.setControls(f.controls);
        newton.options().setSwitch("Verbose", false);
        newton.options().setInt ("MaxIter", 1);
        newton.options().setReal("Tolerance", 1e30);
        newton.solveNewton();
        const real_t E_newton = newton.iterLog()(0, 0);

        syncEnergyOptions(newton, opt);
        gsAsConstVector<real_t> uc(f.controls.data(), f.controls.size());
        const real_t E_ref = opt.evalObj(uc);

        gsTestInfo << "E2 E_newton = " << E_newton << ", E_ref = " << E_ref
                   << ", relErr(E) = " << math::abs(E_newton - E_ref) / math::abs(E_ref) << "\n";
        CHECK(math::abs(E_newton - E_ref) < ENERGY_TOL * math::abs(E_ref));
    }

    TEST(E3_Energy_FusedSweep_vs_EvalObj_Planar_GradientBased)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::GradientBased> newton(f.domain, *geom, monFun, f.tbasis, true);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);

        f.domain.setControls(f.controls);
        newton.options().setSwitch("Verbose", false);
        newton.options().setInt ("MaxIter", 1);
        newton.options().setReal("Tolerance", 1e30);
        newton.solveNewton();
        const real_t E_newton = newton.iterLog()(0, 0);

        syncEnergyOptions(newton, opt);
        gsAsConstVector<real_t> uc(f.controls.data(), f.controls.size());
        const real_t E_ref = opt.evalObj(uc);

        gsTestInfo << "E3 E_newton = " << E_newton << ", E_ref = " << E_ref
                   << ", relErr(E) = " << math::abs(E_newton - E_ref) / math::abs(E_ref) << "\n";
        CHECK(math::abs(E_newton - E_ref) < ENERGY_TOL * math::abs(E_ref));
    }

    TEST(E4_Energy_FusedSweep_vs_EvalObj_Surface_ValueBased)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makeSurfaceGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased> newton(f.domain, *geom, monFun, f.tbasis, true);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);

        f.domain.setControls(f.controls);
        newton.options().setSwitch("Verbose", false);
        newton.options().setInt ("MaxIter", 1);
        newton.options().setReal("Tolerance", 1e30);
        newton.solveNewton();
        const real_t E_newton = newton.iterLog()(0, 0);

        syncEnergyOptions(newton, opt);
        gsAsConstVector<real_t> uc(f.controls.data(), f.controls.size());
        const real_t E_ref = opt.evalObj(uc);

        gsTestInfo << "E4 E_newton = " << E_newton << ", E_ref = " << E_ref
                   << ", relErr(E) = " << math::abs(E_newton - E_ref) / math::abs(E_ref) << "\n";
        CHECK(math::abs(E_newton - E_ref) < ENERGY_TOL * math::abs(E_ref));
    }

    // =====================================================================
    // The iter==0 && minJ<=0 gsWarn in _solveLoop
    // (gsAdaptiveParametrizationNewton.hpp:2015-2018), which exists because
    // the GISMO_ENSURE(R.allFinite(), ...) immediately above it
    // (:1996-1999) is too weak to catch a fold -- a folded start gives
    // det C = (det J)^2 > 0, so the normal-equation residual stays finite
    // even though the line search can never make progress from there.
    // =====================================================================

    TEST(W1_FoldedStart_WarnsOnFirstIteration)
    {
        // gsSquareDomain_test.cpp:433-477 (TEST(DetJCertificate_DetectsFold))
        // shows control 0 of a makeSigmaBasis(2,3) sigma folds under a
        // displacement of order 0.13, as seen on a 401x401 dense grid. But
        // _evalObjAndMinJ instead sweeps the ANALYSIS mesh at hardcoded
        // quA=1/quB=1 -- a much coarser sample -- and a standalone probe
        // against computeMinJacobian() (same fixture, same tbFine below)
        // showed that 0.13125 sliver is invisible to it (minJ=+0.37), while
        // grossly folding with delta=1.0 gives minJ=-3.8: comfortably
        // negative and visible to the coarse sweep.
        const real_t delta = 1.0;

        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();

        // Finer analysis basis than NewtonFixture::tbasis (which has only 2
        // elements): 16 elements gives the quA=1/quB=1 sweep enough
        // quadrature points to see the fold.
        gsKnotVector<real_t> kv(0, 1, 3, 3);
        gsTensorBSplineBasis<2,real_t> tbFine(kv, kv);

        gsSquareDomain<real_t> sigma(makeSigmaBasis(2, 3), true);
        gsVector<real_t> c = sigma.getControls();
        c[0] += delta;
        sigma.setControls(c);   // setControls does NOT clamp (gsSquareDomain.hpp:156-168)

        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased>
            newton(sigma, *geom, tbFine);
        newton.options().setSwitch("Verbose", false);
        newton.options().setInt ("MaxIter", 1);
        newton.options().setReal("Tolerance", 1e30);   // log row 0, then break:
                                                        // no assembly, no line search

        const real_t minJ = newton.computeMinJacobian();
        gsTestInfo << "W1 computeMinJacobian() = " << minJ << "\n";
        CHECK(minJ <= 0.0);   // calibration: the fold is visible to the sweep

        std::string out;
        {
            CoutCapture cap;
            newton.solveNewton();
            out = cap.str();
        }
        // The fragment must be exactly this (not "fold"/"folded"/"sigma is
        // folded"): _solveLoop has a SECOND, pre-existing gsWarn at :1937-1940
        // (the det(J_sigma) coefficient-bound warning) that also fires on this
        // same folded input and also contains "fold"/"folded". Matching on a
        // shorter fragment would pass even with the guard at :2015-2018
        // deleted. This exact string occurs exactly once in
        // gsAdaptiveParametrizationNewton.hpp, so the match is unambiguous.
        CHECK(out.find("sigma is folded at the starting mesh") != std::string::npos);

        // Negative control: an unfolded start (plain identity sigma, no
        // displacement) must not produce this warning.
        gsSquareDomain<real_t> sigmaOK(makeSigmaBasis(2, 3), true);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased>
            newtonOK(sigmaOK, *geom, tbFine);
        newtonOK.options().setSwitch("Verbose", false);
        newtonOK.options().setInt ("MaxIter", 1);
        newtonOK.options().setReal("Tolerance", 1e30);

        std::string outOK;
        {
            CoutCapture cap;
            newtonOK.solveNewton();
            outOK = cap.str();
        }
        CHECK(outOK.find("sigma is folded at the starting mesh") == std::string::npos);
    }

    TEST(W2_FoldedStart_SurvivesAllFiniteAndLogsNonPositiveMinJ)
    {
        // Companion documenting WHY the guard exists: on the same folded
        // start, solveNewton() must return without throwing (the
        // GISMO_ENSURE(R.allFinite(), ...) at :1996-1999 does NOT catch a
        // fold) while iterLog()(2,0) -- the SAME minJ local the warning
        // branches on -- is non-positive. Same calibrated delta as W1 (see
        // its comment): 1.0 gives computeMinJacobian()=-3.8 against tbFine.
        const real_t delta = 1.0;

        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();

        gsKnotVector<real_t> kv(0, 1, 3, 3);
        gsTensorBSplineBasis<2,real_t> tbFine(kv, kv);

        gsSquareDomain<real_t> sigma(makeSigmaBasis(2, 3), true);
        gsVector<real_t> c = sigma.getControls();
        c[0] += delta;
        sigma.setControls(c);

        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased>
            newton(sigma, *geom, tbFine);
        newton.options().setSwitch("Verbose", false);
        newton.options().setInt ("MaxIter", 1);
        newton.options().setReal("Tolerance", 1e30);

        {
            // Warning text is swallowed here (W1 already covers its
            // content); this test only asserts on iterLog()/no-throw.
            CoutCapture cap;
            newton.solveNewton();   // must not throw despite the fold
        }
        const real_t minJ = newton.iterLog()(2, 0);
        gsTestInfo << "W2 folded iterLog()(2,0) = " << minJ << "\n";
        CHECK(minJ <= 0.0);

        // Contrasting unfolded start: minJ must be strictly positive, so
        // the check above is not vacuously true for every input.
        gsSquareDomain<real_t> sigmaOK(makeSigmaBasis(2, 3), true);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::ValueBased>
            newtonOK(sigmaOK, *geom, tbFine);
        newtonOK.options().setSwitch("Verbose", false);
        newtonOK.options().setInt ("MaxIter", 1);
        newtonOK.options().setReal("Tolerance", 1e30);

        {
            CoutCapture cap;
            newtonOK.solveNewton();
        }
        const real_t minJOK = newtonOK.iterLog()(2, 0);
        gsTestInfo << "W2 unfolded iterLog()(2,0) = " << minJOK << "\n";
        CHECK(minJOK > 0.0);
    }

}
