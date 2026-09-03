#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    index_t numRefine     = 3;
    index_t compRefine    = 3;
    index_t deg           = 3;
    index_t numIter       = 5;
    bool runFD            = false;

    gsCmdLine cmd("Profile evalObj / gradObj_into at scale.");
    cmd.addInt("r", "refine",     "Uniform refinements of integration basis", numRefine);
    cmd.addInt("R", "compRefine", "Uniform refinements of composition (controls ~4*(2^R-1)^2)", compRefine);
    cmd.addInt("p", "degree",     "Polynomial degree of composition", deg);
    cmd.addInt("n", "niter",      "Number of repeated evaluations for timing", numIter);
    cmd.addSwitch("fd", "Also run FD gradient (slow!)", runFD);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsGeometry<>::uPtr composition = gsNurbsCreator<>::BSplineSquareDeg(deg);
    for (index_t i = 0; i < compRefine; ++i)
        composition->uniformRefine();
    gsSquareDomain<real_t> domain(*composition);

    gsVector<real_t> controls = domain.getControls();

    gsKnotVector<> kv0(0, 1, 0, deg + 1);
    for (index_t i = 0; i < numRefine; ++i)
        kv0.uniformRefine();

    gsTensorBSplineBasis<2> tbasis(kv0, kv0);
    gsComposedBasis<> cbasis(*composition, tbasis);

    gsMatrix<> anchors = cbasis.anchors().transpose();
    gsMatrix<> coefs(cbasis.size(), 2);
    coefs = anchors;
    for (index_t i = 0; i < coefs.rows(); ++i)
    {
        coefs(i,0) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs(i,1));
        coefs(i,1) += 0.05 * math::sin(2.0 * EIGEN_PI * coefs(i,0));
    }
    gsGeometry<>::uPtr geom = tbasis.makeGeometry(coefs);

    gsFunctionExpr<> monFun("1 + 0.5*sin(pi*x)*sin(pi*y)", 2);

    gsOptMesh<real_t, MonitorMode::GradientBased> opt(domain, *geom, &monFun, &tbasis, true);

    const index_t nc = opt.numDesignVars();
    gsInfo << "Degree:           " << deg << "\n";
    gsInfo << "Comp refinements: " << compRefine << "\n";
    gsInfo << "Int refinements:  " << numRefine << "\n";
    gsInfo << "Basis size:       " << tbasis.size() << "\n";
    gsInfo << "Controls:         " << nc << "\n";
    gsInfo << "Iterations:       " << numIter << "\n\n";

    gsAsConstVector<real_t> asu(controls.data(), controls.size());

    gsStopwatch timer;

    timer.restart();
    real_t objVal = 0;
    for (index_t it = 0; it < numIter; ++it)
        objVal = opt.evalObj(asu);
    real_t t_eval = timer.stop();
    gsInfo << "evalObj:       " << t_eval << " s  (" << numIter << " calls, "
           << t_eval / numIter << " s/call)  val=" << objVal << "\n";

    gsVector<real_t> grad(nc);
    gsAsVector<real_t> asgrad(grad.data(), grad.rows());
    timer.restart();
    for (index_t it = 0; it < numIter; ++it)
        opt.gradObj_into(asu, asgrad);
    real_t t_grad = timer.stop();
    gsInfo << "gradObj_into:  " << t_grad << " s  (" << numIter << " calls, "
           << t_grad / numIter << " s/call)  |grad|=" << grad.norm() << "\n";

    if (runFD)
    {
        gsVector<real_t> grad_fd(nc);
        gsAsVector<real_t> asgrad_fd(grad_fd.data(), grad_fd.rows());
        timer.restart();
        opt.gradObj_FD_into(asu, asgrad_fd);
        real_t t_fd = timer.stop();
        gsInfo << "gradObj_FD:    " << t_fd << " s  (1 call)\n";

        real_t relErr = (grad - grad_fd).norm() / (grad_fd.norm() > 1e-15 ? grad_fd.norm() : 1.0);
        gsInfo << "FD relErr:     " << relErr << "\n";
        gsInfo << "Speedup (analytical/FD): " << (t_fd / (t_grad / numIter)) << "x\n";
    }

    gsInfo << "\nGrad/Eval ratio: " << (t_grad / t_eval) << "x\n";

    return 0;
}
