#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>

using namespace gismo;
using namespace gismo::expr;

bool testGradient(const std::string & label,
                  gsOptMesh<real_t,MonitorMode::ValueBased> & opt,
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

    gsInfo << "\n===== Summary =====\n";
    if (nFail == 0)
        gsInfo << "All tests passed!\n";
    else
        gsInfo << nFail << " test(s) failed!\n";

    return nFail;
}
