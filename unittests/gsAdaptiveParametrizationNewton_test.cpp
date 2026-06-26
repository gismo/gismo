/** @file gsAdaptiveParametrizationNewton_test.cpp

    @brief Tests for gsAdaptiveParametrizationNewton: independent residual
    assembly (vs gradObj_into), analytic Hessian (vs FD of the analytic
    gradient + symmetry), Picard SPD property, and Newton convergence order.

    See doc/derivation_hessian.md for the assembled formulas.

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

#ifdef GISMO_WITH_ADIFF
static const real_t RES_TOL  = 1e-9;   // residual vs gradObj_into (same quadrature => machine precision)
static const real_t HESS_TOL = 5e-6;   // analytic K vs central-FD of analytic gradient (h=1e-7)
#else
static const real_t RES_TOL  = 1e-9;   // independent of f-derivative quality (same f calls)
static const real_t HESS_TOL = 1e-2;   // gsFunctionExpr derivs are FD without ADIFF
#endif

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

// Relative error of assembleResidual vs gradObj_into (same options)
template <MonitorMode MODE>
real_t residualRelError(gsAdaptiveParametrizationNewton<real_t,MODE> & newton,
                        gsOptMesh<real_t,MODE> & opt,
                        const gsVector<real_t> & u)
{
    opt.options().setReal("Penalty",   newton.options().getReal("Penalty"));
    opt.options().setReal("Smoothing", newton.options().getReal("Smoothing"));

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
    opt.options().setReal("Penalty",   newton.options().getReal("Penalty"));
    opt.options().setReal("Smoothing", newton.options().getReal("Smoothing"));

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
        // d^3 f via FD-of-deriv2 (h=1e-5): slightly relaxed tolerance
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::GradientBased> newton(f.domain, *geom, monFun, f.tbasis, true);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, true);
        real_t symErr;
        real_t err = hessianRelError(newton, opt, f.controls, symErr);
        gsTestInfo << "H7 relErr = " << err << ", symErr = " << symErr << "\n";
        CHECK(symErr < 1e-8);
        CHECK(err < 100*HESS_TOL);
    }

    TEST(H8_Hessian_Planar_GradientBased_Physical)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makePlanarGeom();
        gsFunctionExpr<> monFun("1 + 0.5*x*x + 0.3*y*y", 2);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::GradientBased> newton(f.domain, *geom, monFun, f.tbasis, false);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, false);
        real_t symErr;
        real_t err = hessianRelError(newton, opt, f.controls, symErr);
        gsTestInfo << "H8 relErr = " << err << ", symErr = " << symErr << "\n";
        CHECK(symErr < 1e-8);
        CHECK(err < 100*HESS_TOL);
    }

    TEST(H9_Hessian_Surface_GradientBased_Physical)
    {
        NewtonFixture f;
        gsGeometry<>::uPtr geom = f.makeSurfaceGeom();
        gsFunctionExpr<> monFun("1 + 0.3*x + 0.2*y + 0.1*z", 3);
        gsAdaptiveParametrizationNewton<real_t,MonitorMode::GradientBased> newton(f.domain, *geom, monFun, f.tbasis, false);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(f.domain, *geom, &monFun, &f.tbasis, false);
        real_t symErr;
        real_t err = hessianRelError(newton, opt, f.controls, symErr);
        gsTestInfo << "H9 relErr = " << err << ", symErr = " << symErr << "\n";
        CHECK(symErr < 1e-8);
        CHECK(err < 100*HESS_TOL);
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

}
