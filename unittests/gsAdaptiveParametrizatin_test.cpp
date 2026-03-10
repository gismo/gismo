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
#include <gsAssembler/gsAdaptiveParametrization.h>

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

}
