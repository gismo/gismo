#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // 1. Geometry (Identity in 2D)
    gsKnotVector<> kv({0,0,0,1,1,1}, 2);
    gsTensorBSplineBasis<2> geomBasis(kv, kv);
    gsMatrix<> geomCoefs = geomBasis.anchors().transpose();
    gsGeometry<>::uPtr geom = geomBasis.makeGeometry(geomCoefs);

    // 2. Composition (Identity in 2D)
    gsGeometry<>::uPtr composition = gsNurbsCreator<>::BSplineSquareDeg(2);
    gsSquareDomain<real_t> domain(*composition);
    domain.setAllActive();
    
    // 3. Integration Basis
    gsKnotVector<> kv_int({0,0,0,0.5,1,1,1}, 2);
    gsTensorBSplineBasis<2> intBasis(kv, kv_int);

    // 4. OptMesh
    gsOptMesh<real_t,MonitorMode::ValueBased> opt(domain,*geom,&intBasis);
    
    gsInfo << "Domain nControls: " << domain.nControls() << "\n";
    gsInfo << "Opt numDesignVars: " << opt.numDesignVars() << "\n";

    // Identity controls
    gsVector<real_t> controls(domain.nControls());
    gsMatrix<> anchors = composition->basis().anchors();
    for (index_t i=0; i<anchors.cols(); ++i) {
        controls(i) = anchors(0,i);
        controls(i + anchors.cols()) = anchors(1,i);
    }
    
    gsInfo << "Evaluating objective at identity...\n";
    real_t val = opt.evalObj(gsAsConstVector<real_t>(controls.data(), controls.size()));
    gsInfo << "Objective value: " << val << " (expected 2.0 for 2D identity Winslow)\n";

    gsVector<real_t> grad(controls.size());
    gsAsVector<real_t> asgrad(grad.data(), grad.rows());
    opt.gradObj_into(gsAsConstVector<real_t>(controls.data(), controls.size()), asgrad);
    
    gsInfo << "Analytical gradient norm: " << grad.norm() << "\n";
    
    // Finite difference check
    real_t eps = 1e-6;
    gsVector<real_t> grad_fd(controls.size());
    for (index_t i = 0; i < controls.size(); ++i)
    {
        gsVector<real_t> controls_plus = controls;
        gsVector<real_t> controls_minus = controls;
        controls_plus(i) += eps;
        controls_minus(i) -= eps;
        real_t f_plus = opt.evalObj(gsAsConstVector<real_t>(controls_plus.data(), controls_plus.size()));
        real_t f_minus = opt.evalObj(gsAsConstVector<real_t>(controls_minus.data(), controls_minus.size()));
        grad_fd(i) = (f_plus - f_minus) / (2 * eps);
    }
    
    gsInfo << "FD gradient norm: " << grad_fd.norm() << "\n";
    gsInfo << "Gradient difference norm: " << (grad - grad_fd).norm() << "\n";

    return EXIT_SUCCESS;
}
