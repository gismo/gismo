#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>

using namespace gismo;
using namespace gismo::expr;

template <MonitorMode MODE>
bool testGradient(const std::string & label,
                  gsOptMesh<real_t,MODE> & opt,
                  const gsVector<real_t> & controls,
                  real_t tol = 1e-6)
{
    const index_t nc = controls.size();

    gsStopwatch timer;

    gsVector<real_t> grad(nc);
    gsAsVector<real_t> asgrad(grad.data(), grad.rows());
    timer.restart();
    opt.gradObj_into(gsAsConstVector<real_t>(controls.data(), controls.size()), asgrad);
    real_t grad_time = timer.stop();

    gsVector<real_t> grad_fd(nc);
    gsAsVector<real_t> asgrad_fd(grad_fd.data(), grad_fd.rows());
    timer.restart();
    opt.gradObj_FD_into(gsAsConstVector<real_t>(controls.data(), controls.size()), asgrad_fd);
    real_t grad_FD_time = timer.stop();

    real_t absErr = (grad - grad_fd).norm();
    real_t relErr = absErr / (grad_fd.norm() > 1e-15 ? grad_fd.norm() : 1.0);

    bool pass = relErr < tol;
    gsInfo << "[" << (pass ? "PASS" : "FAIL") << "] " << label
           << "  (relErr=" << relErr << ")\n";
    gsInfo << "  Time:        " << grad_time << "s (analytical), " << grad_FD_time << "s (FD)\n";

    if (!pass || relErr > 1e-8)
    {
        gsInfo << "  Analytical: " << grad.transpose() << "\n";
        gsInfo << "  FD:         " << grad_fd.transpose() << "\n";
        gsInfo << "  Diff:       " << (grad - grad_fd).transpose() << "\n";
        gsInfo << "  Time:        " << grad_time << "s (analytical), " << grad_FD_time << "s (FD)\n";
    }
    return pass;
}

int main(int argc, char *argv[])
{
    int nFail = 0;

    gsGeometry<>::uPtr composition = gsNurbsCreator<>::BSplineSquareDeg(2);
    gsSquareDomain<real_t> domain(*composition);

    gsVector<real_t> controls(domain.nControls());
    controls << 0.95, 0.95;

    gsKnotVector<> kv1({0,0,0,1,1,1}, 2);
    gsKnotVector<> kv2({0,0,0,0.50,1,1,1}, 2);
    gsTensorBSplineBasis<2> tbasis(kv1, kv2);
    gsComposedBasis<> cbasis(*composition, tbasis);

    gsMatrix<> anchors = cbasis.anchors().transpose();

    gsInfo << "===== B1: Surface (2D->3D), no monitor =====\n";
    {
        gsMatrix<> coefs3(cbasis.size(), 3);
        coefs3.leftCols(2) = anchors;
        for (index_t i = 0; i < coefs3.rows(); ++i)
            coefs3(i,2) = math::sin(coefs3(i,0)*EIGEN_PI) * math::cos(coefs3(i,1)*EIGEN_PI);
        gsGeometry<>::uPtr geom3 = tbasis.makeGeometry(coefs3);

        gsOptMesh<real_t,MonitorMode::ValueBased> opt(domain, *geom3, &tbasis);
        if (!testGradient("B1: tr(Cinv)/detG, analytical", opt, controls))
            nFail++;
    }

    gsInfo << "\n===== A1: Planar (2D->2D), no monitor =====\n";
    {
        gsMatrix<> coefs2(cbasis.size(), 2);
        coefs2 = anchors;
        for (index_t i = 0; i < coefs2.rows(); ++i)
        {
            coefs2(i,0) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,1));
            coefs2(i,1) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,0));
        }
        gsGeometry<>::uPtr geom2 = tbasis.makeGeometry(coefs2);

        gsOptMesh<real_t,MonitorMode::ValueBased> opt(domain, *geom2, &tbasis);
        if (!testGradient("A1: tr(Cinv)/detG, analytical", opt, controls))
            nFail++;
    }

    gsInfo << "\n===== B2: Surface (2D->3D), with monitor =====\n";
    {
        gsMatrix<> coefs3(cbasis.size(), 3);
        coefs3.leftCols(2) = anchors;
        for (index_t i = 0; i < coefs3.rows(); ++i)
            coefs3(i,2) = math::sin(coefs3(i,0)*EIGEN_PI) * math::cos(coefs3(i,1)*EIGEN_PI);
        gsGeometry<>::uPtr geom3 = tbasis.makeGeometry(coefs3);

        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(domain, *geom3, &monFun, &tbasis, true);
        if (!testGradient("B2: m2*tr(Cinv)*detG, analytical", opt, controls))
            nFail++;
    }

    gsInfo << "\n===== A2: Planar (2D->2D), with monitor =====\n";
    {
        gsMatrix<> coefs2(cbasis.size(), 2);
        coefs2 = anchors;
        for (index_t i = 0; i < coefs2.rows(); ++i)
        {
            coefs2(i,0) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,1));
            coefs2(i,1) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,0));
        }
        gsGeometry<>::uPtr geom2 = tbasis.makeGeometry(coefs2);

        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(domain, *geom2, &monFun, &tbasis, true);
        if (!testGradient("A2: m2*tr(Cinv)*detG, analytical", opt, controls))
            nFail++;
    }

    gsInfo << "\n===== C1: Surface (2D->3D), GradientBased monitor, parametric =====\n";
    {
        gsMatrix<> coefs3(cbasis.size(), 3);
        coefs3.leftCols(2) = anchors;
        for (index_t i = 0; i < coefs3.rows(); ++i)
            coefs3(i,2) = math::sin(coefs3(i,0)*EIGEN_PI) * math::cos(coefs3(i,1)*EIGEN_PI);
        gsGeometry<>::uPtr geom3 = tbasis.makeGeometry(coefs3);

        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(domain, *geom3, &monFun, &tbasis, true);
        if (!testGradient<MonitorMode::GradientBased>("C1: GradBased surface, parametric", opt, controls))
            nFail++;
    }

    gsInfo << "\n===== C2: Planar (2D->2D), GradientBased monitor, parametric =====\n";
    {
        gsMatrix<> coefs2(cbasis.size(), 2);
        coefs2 = anchors;
        for (index_t i = 0; i < coefs2.rows(); ++i)
        {
            coefs2(i,0) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,1));
            coefs2(i,1) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,0));
        }
        gsGeometry<>::uPtr geom2 = tbasis.makeGeometry(coefs2);

        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(domain, *geom2, &monFun, &tbasis, true);
        if (!testGradient<MonitorMode::GradientBased>("C2: GradBased planar, parametric", opt, controls))
            nFail++;
    }

    gsInfo << "\n===== C3: Surface (2D->3D), GradientBased monitor, physical =====\n";
    {
        gsMatrix<> coefs3(cbasis.size(), 3);
        coefs3.leftCols(2) = anchors;
        for (index_t i = 0; i < coefs3.rows(); ++i)
            coefs3(i,2) = math::sin(coefs3(i,0)*EIGEN_PI) * math::cos(coefs3(i,1)*EIGEN_PI);
        gsGeometry<>::uPtr geom3 = tbasis.makeGeometry(coefs3);

        // Monitor function defined on physical space (3D)
        gsFunctionExpr<> monFun("1 + 0.3*x + 0.2*y + 0.1*z", 3);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(domain, *geom3, &monFun, &tbasis, false);
        if (!testGradient<MonitorMode::GradientBased>("C3: GradBased surface, physical", opt, controls))
            nFail++;
    }

    gsInfo << "\n===== C4: Planar (2D->2D), GradientBased monitor, physical =====\n";
    {
        gsMatrix<> coefs2(cbasis.size(), 2);
        coefs2 = anchors;
        for (index_t i = 0; i < coefs2.rows(); ++i)
        {
            coefs2(i,0) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,1));
            coefs2(i,1) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,0));
        }
        gsGeometry<>::uPtr geom2 = tbasis.makeGeometry(coefs2);

        // Monitor function defined on physical space (2D)
        gsFunctionExpr<> monFun("1 + 0.5*x*x + 0.3*y*y", 2);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(domain, *geom2, &monFun, &tbasis, false);
        if (!testGradient<MonitorMode::GradientBased>("C4: GradBased planar, physical", opt, controls))
            nFail++;
    }

    gsInfo << "\n===== C5: Planar (2D->2D), GradientBased, parametric, LINEAR monitor (term2=0) =====\n";
    {
        gsMatrix<> coefs2(cbasis.size(), 2);
        coefs2 = anchors;
        for (index_t i = 0; i < coefs2.rows(); ++i)
        {
            coefs2(i,0) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,1));
            coefs2(i,1) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,0));
        }
        gsGeometry<>::uPtr geom2 = tbasis.makeGeometry(coefs2);
        // Linear monitor: Hessian is zero, so term2=0. Only term1 is non-trivial.
        gsFunctionExpr<> monFun("0.3*x + 0.7*y", 2);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(domain, *geom2, &monFun, &tbasis, true);
        if (!testGradient<MonitorMode::GradientBased>("C5: GradBased planar, linear monitor (term2=0)", opt, controls))
            nFail++;
    }

    gsInfo << "\n===== C6: Planar (2D->2D), GradientBased, parametric, CONSTANT gradient (term1 check) =====\n";
    {
        // Use a geometry that is the identity (Jg = I), so Cg = I, Cg_inv = I.
        // Then eta2 = ||grad_xi_f||^2 (simple Euclidean norm).
        // dCg/dalpha = 0 (Jg is constant), so term1 = 0. Only term2 remains.
        // With f=0.3*x+0.7*y and identity geometry: eta2 = 0.09+0.49 = 0.58 (constant), dm2/dalpha = 0.
        // But with Hess=0: dm2_dalpha = 0, so gradient reduces to m2 * (standard geometric gradient).
        // This tests: is m2 correct? Is the C-based part of dE correct?
        gsMatrix<> coefs2_id(cbasis.size(), 2);
        coefs2_id = anchors; // identity-like geometry (no deformation)
        gsGeometry<>::uPtr geom2_id = tbasis.makeGeometry(coefs2_id);
        gsFunctionExpr<> monFun("0.3*x + 0.7*y", 2);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(domain, *geom2_id, &monFun, &tbasis, true);
        if (!testGradient<MonitorMode::GradientBased>("C6: GradBased identity geom, linear monitor", opt, controls))
            nFail++;
    }

    gsInfo << "\n===== C7: GradBased, theta=0 (m2=1 const, dm2=0; tests geometric part only) =====\n";
    {
        gsMatrix<> coefs2(cbasis.size(), 2);
        coefs2 = anchors;
        for (index_t i = 0; i < coefs2.rows(); ++i)
        {
            coefs2(i,0) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,1));
            coefs2(i,1) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs2(i,0));
        }
        gsGeometry<>::uPtr geom2 = tbasis.makeGeometry(coefs2);
        gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(domain, *geom2, &monFun, &tbasis, true);
        opt.options().setReal("Smoothing", 0.0); // theta=0: m2=1 always, dm2/dalpha=0
        if (!testGradient<MonitorMode::GradientBased>("C7: GradBased theta=0 (geometric part only)", opt, controls))
            nFail++;
    }

    // ===========================================================================
    // 3D volumetric (planar) tests: dd=3, td=3
    // Composition: unit cube sigma (gsSquareDomain built from BSplineCube(1))
    // Geometry:    slightly deformed unit cube (td=dd=3)
    // ===========================================================================

    gsGeometry<>::uPtr composition3d = gsNurbsCreator<real_t>::BSplineCube(2);
    gsSquareDomain<real_t> domain3d(*composition3d);

    // Free controls for 3D domain (all interior + edge dofs, boundary locked).
    gsVector<real_t> controls3d = domain3d.getControls();

    // Build a 3D tensor B-spline basis for the geometry (deg 2 in each direction)
    gsKnotVector<> kv3_1({0,0,0,1,1,1}, 2);
    gsKnotVector<> kv3_2({0,0,0,0.5,1,1,1}, 2);
    gsKnotVector<> kv3_3({0,0,0,1,1,1}, 2);
    gsTensorBSplineBasis<3> tbasis3(kv3_1, kv3_2, kv3_3);

    gsComposedBasis<> cbasis3(*composition3d, tbasis3);
    gsMatrix<> anchors3 = cbasis3.anchors().transpose(); // nDof x 3

    // Perturb the 3D controls away from the symmetric center to get a non-trivial gradient.
    // Without perturbation, controls3d = (0.5, 0.5, 0.5) is a symmetry point where the
    // gradient is exactly zero (fp - fm = O(machine_eps) with h=1e-7 => FD noise only).
    controls3d(0) += 0.04;
    controls3d(1) -= 0.03;
    controls3d(2) += 0.02;
    domain3d.setControls(controls3d);

    gsInfo << "\n===== D1: Planar (3D->3D), no monitor =====\n";
    {
        // Slightly deformed cube geometry (td=dd=3)
        gsMatrix<> coefs3d(cbasis3.size(), 3);
        coefs3d = anchors3;
        for (index_t i = 0; i < coefs3d.rows(); ++i)
        {
            coefs3d(i,0) += 0.03 * math::sin(2*EIGEN_PI * coefs3d(i,1)) * math::sin(2*EIGEN_PI * coefs3d(i,2));
            coefs3d(i,1) += 0.03 * math::sin(2*EIGEN_PI * coefs3d(i,2)) * math::sin(2*EIGEN_PI * coefs3d(i,0));
            coefs3d(i,2) += 0.03 * math::sin(2*EIGEN_PI * coefs3d(i,0)) * math::sin(2*EIGEN_PI * coefs3d(i,1));
        }
        gsGeometry<>::uPtr geom3d = tbasis3.makeGeometry(coefs3d);

        gsOptMesh<real_t,MonitorMode::ValueBased> opt(domain3d, *geom3d, &tbasis3);
        if (!testGradient("D1: 3D planar no monitor", opt, controls3d))
            nFail++;
    }

    gsInfo << "\n===== D2: Planar (3D->3D), sphere monitor (ValueBased) =====\n";
    {
        gsMatrix<> coefs3d(cbasis3.size(), 3);
        coefs3d = anchors3;
        for (index_t i = 0; i < coefs3d.rows(); ++i)
        {
            coefs3d(i,0) += 0.03 * math::sin(2*EIGEN_PI * coefs3d(i,1)) * math::sin(2*EIGEN_PI * coefs3d(i,2));
            coefs3d(i,1) += 0.03 * math::sin(2*EIGEN_PI * coefs3d(i,2)) * math::sin(2*EIGEN_PI * coefs3d(i,0));
            coefs3d(i,2) += 0.03 * math::sin(2*EIGEN_PI * coefs3d(i,0)) * math::sin(2*EIGEN_PI * coefs3d(i,1));
        }
        gsGeometry<>::uPtr geom3d = tbasis3.makeGeometry(coefs3d);

        // Sphere-shaped monitor: concentrated at (0.5, 0.5, 0.5)
        gsFunctionExpr<> monFun("exp(-50*((x-0.5)^2+(y-0.5)^2+(z-0.5)^2))", 3);
        gsOptMesh<real_t,MonitorMode::ValueBased> opt(domain3d, *geom3d, &monFun, &tbasis3, true);
        if (!testGradient("D2: 3D planar sphere monitor (ValueBased)", opt, controls3d))
            nFail++;
    }

    gsInfo << "\n===== D3: Planar (3D->3D), sphere monitor (GradientBased, parametric) =====\n";
    {
        gsMatrix<> coefs3d(cbasis3.size(), 3);
        coefs3d = anchors3;
        for (index_t i = 0; i < coefs3d.rows(); ++i)
        {
            coefs3d(i,0) += 0.03 * math::sin(2*EIGEN_PI * coefs3d(i,1)) * math::sin(2*EIGEN_PI * coefs3d(i,2));
            coefs3d(i,1) += 0.03 * math::sin(2*EIGEN_PI * coefs3d(i,2)) * math::sin(2*EIGEN_PI * coefs3d(i,0));
            coefs3d(i,2) += 0.03 * math::sin(2*EIGEN_PI * coefs3d(i,0)) * math::sin(2*EIGEN_PI * coefs3d(i,1));
        }
        gsGeometry<>::uPtr geom3d = tbasis3.makeGeometry(coefs3d);

        // Sphere-shaped monitor: concentrated at (0.5, 0.5, 0.5), parametric
        gsFunctionExpr<> monFun("exp(-50*((x-0.5)^2+(y-0.5)^2+(z-0.5)^2))", 3);
        gsOptMesh<real_t,MonitorMode::GradientBased> opt(domain3d, *geom3d, &monFun, &tbasis3, true);
        if (!testGradient<MonitorMode::GradientBased>("D3: 3D planar sphere monitor (GradBased, parametric)", opt, controls3d))
            nFail++;
    }

    gsInfo << "\n===== Summary =====\n";
    if (nFail == 0)
        gsInfo << "All tests passed!\n";
    else
        gsInfo << nFail << " test(s) failed!\n";

    return nFail;
}