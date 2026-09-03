#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>

using namespace gismo;
using namespace gismo::expr;

int main(int argc, char *argv[])
{
    gsGeometry<>::uPtr composition = gsNurbsCreator<>::BSplineSquareDeg(2);
    gsSquareDomain<real_t> domain(*composition);

    gsKnotVector<> kv1({0,0,0,1,1,1}, 2);
    gsKnotVector<> kv2({0,0,0,0.50,1,1,1}, 2);
    gsTensorBSplineBasis<2> tbasis(kv1, kv2);

    gsComposedBasis<> cbasis(*composition, tbasis);

    gsMatrix<> coefs(cbasis.size(),3);
    coefs.leftCols(2) = cbasis.anchors().transpose();
    for (index_t i=0; i<coefs.rows(); ++i)
        coefs(i,2) = math::sin(coefs(i,0)*M_PI) * math::cos(coefs(i,1)*M_PI);

    gsGeometry<>::uPtr geom = tbasis.makeGeometry(coefs);

    gsVector<real_t> controls(domain.nControls());
    controls << 0.95, 0.95;
    domain.setControls(controls);

    gsMatrix<real_t> pt(2,1);
    pt << 0.3, 0.7;

    gsInfo << "=== STEP 1: Verify J_sigma layout ===\n";
    gsMatrix<real_t> sigma_val;
    domain.eval_into(pt, sigma_val);
    gsInfo << "sigma(pt) = " << sigma_val.transpose() << "\n";

    gsMatrix<real_t> Jsigma_flat;
    domain.deriv_into(pt, Jsigma_flat);
    gsInfo << "deriv_into flat = " << Jsigma_flat.transpose() << "\n";

    gsMatrix<real_t> JsT = Jsigma_flat.col(0).reshaped(2, 2);
    gsInfo << "Js^T (reshaped(dd,dd)):\n" << JsT << "\n";
    gsMatrix<real_t> Js = JsT.transpose();
    gsInfo << "Js = Js^T.transpose():\n" << Js << "\n";

    const real_t eps = 1e-7;
    gsMatrix<real_t> ptx = pt; ptx(0) += eps;
    gsMatrix<real_t> pty = pt; pty(1) += eps;
    gsMatrix<real_t> sx, sy;
    domain.eval_into(ptx, sx);
    domain.eval_into(pty, sy);
    gsMatrix<real_t> Js_fd(2,2);
    Js_fd.col(0) = (sx - sigma_val) / eps;
    Js_fd.col(1) = (sy - sigma_val) / eps;
    gsInfo << "Js (FD):\n" << Js_fd << "\n";
    gsInfo << "Js error = " << (Js - Js_fd).norm() << "\n\n";

    gsInfo << "=== STEP 2: Verify J_geom layout ===\n";
    gsMatrix<real_t> Jgeom_flat;
    geom->deriv_into(sigma_val, Jgeom_flat);
    gsInfo << "J_geom deriv_into flat = " << Jgeom_flat.transpose() << "\n";
    gsMatrix<real_t> JgT = Jgeom_flat.col(0).reshaped(2, 3);
    gsMatrix<real_t> Jg = JgT.transpose();
    gsInfo << "Jg (3x2):\n" << Jg << "\n";

    gsMatrix<real_t> gval, gx, gy;
    geom->eval_into(sigma_val, gval);
    gsMatrix<real_t> sx2 = sigma_val; sx2(0) += eps;
    gsMatrix<real_t> sy2 = sigma_val; sy2(1) += eps;
    geom->eval_into(sx2, gx);
    geom->eval_into(sy2, gy);
    gsMatrix<real_t> Jg_fd(3,2);
    Jg_fd.col(0) = (gx - gval) / eps;
    Jg_fd.col(1) = (gy - gval) / eps;
    gsInfo << "Jg (FD):\n" << Jg_fd << "\n";
    gsInfo << "Jg error = " << (Jg - Jg_fd).norm() << "\n\n";

    gsInfo << "=== STEP 3: Verify J_composed ===\n";
    gsMatrix<real_t> Jc = Jg * Js;
    gsInfo << "Jc = Jg * Js:\n" << Jc << "\n";

    gsComposedBasis<> cbasis2(*composition, tbasis);
    gsGeometry<>::uPtr cgeom = cbasis2.makeGeometry(coefs);
    gsMatrix<real_t> cval;
    cgeom->eval_into(pt, cval);
    gsMatrix<real_t> cx, cy;
    gsMatrix<real_t> ptxc = pt; ptxc(0) += eps;
    gsMatrix<real_t> ptyc = pt; ptyc(1) += eps;
    cgeom->eval_into(ptxc, cx);
    cgeom->eval_into(ptyc, cy);
    gsMatrix<real_t> Jc_fd(3,2);
    Jc_fd.col(0) = (cx - cval) / eps;
    Jc_fd.col(1) = (cy - cval) / eps;
    gsInfo << "Jc (FD from composed geom):\n" << Jc_fd << "\n";
    gsInfo << "Jc error = " << (Jc - Jc_fd).norm() << "\n\n";

    gsInfo << "=== STEP 4: Verify integrand E ===\n";
    gsMatrix<real_t> A = Jc.transpose() * Jc;
    real_t trA = A.trace();
    real_t detA = A.determinant();
    real_t meas_c = math::sqrt(detA);
    real_t E = trA / meas_c;
    gsInfo << "E = trA/meas_c = " << E << "\n";
    gsInfo << "trA = " << trA << ", detA = " << detA << ", meas_c = " << meas_c << "\n\n";

    gsInfo << "=== STEP 5: FD of integrand dE/dalpha ===\n";
    gsMatrix<real_t> Ainv = A.inverse();

    const index_t nc = domain.nControls();
    const index_t dd = 2;
    const index_t td = 3;
    const gsDofMapper & mapper = domain.mapper();
    const gsBasis<real_t> & sigmaBasis = domain.domain().basis();
    const index_t S = dd * (dd + 1) / 2;

    gsMatrix<real_t> deriv2_geom;
    geom->deriv2_into(sigma_val, deriv2_geom);

    gsMatrix<index_t> actives;
    sigmaBasis.active_into(pt, actives);

    gsMatrix<real_t> basisVals;
    sigmaBasis.eval_into(pt, basisVals);

    gsMatrix<real_t> basisDerivs;
    sigmaBasis.deriv_into(pt, basisDerivs);

    gsVector<real_t> dE_an(nc);
    dE_an.setZero();

    gsInfo << "deriv2_geom at sigma_val:\n" << deriv2_geom.col(0).transpose() << "\n";

    const index_t nActive = actives.rows();
    for (index_t loc = 0; loc != nActive; ++loc)
    {
        const index_t k = actives(loc, 0);
        real_t Nk = basisVals(loc, 0);

        gsMatrix<real_t> gradNk(dd, 1);
        for (index_t j = 0; j != dd; ++j)
            gradNk(j) = basisDerivs(loc * dd + j, 0);

        for (index_t d = 0; d != dd; ++d)
        {
            if (!mapper.is_free(k, 0, d))
                continue;
            index_t ii = mapper.index(k, 0, d);

            gsMatrix<real_t> dJs = gsMatrix<real_t>::Zero(dd, dd);
            for (index_t j = 0; j != dd; ++j)
                dJs(d, j) = gradNk(j);

            gsMatrix<real_t> dJg = gsMatrix<real_t>::Zero(td, dd);
            for (index_t a = 0; a != td; ++a)
            {
                for (index_t j = 0; j != dd; ++j)
                {
                    index_t lo = math::min(d, j);
                    index_t hi = math::max(d, j);
                    index_t hess_idx = (lo == hi) ? lo : dd + lo * (2 * dd - lo - 3) / 2 + hi - 1;
                    dJg(a, j) = Nk * deriv2_geom(a * S + hess_idx, 0);
                }
            }

            gsMatrix<real_t> dJc = dJg * Js + Jg * dJs;
            gsMatrix<real_t> dA = dJc.transpose() * Jc + Jc.transpose() * dJc;
            real_t dE_val = (dA.trace() - trA / 2.0 * (Ainv * dA).trace()) / meas_c;
            dE_an(ii) = dE_val;

            gsInfo << "  control " << ii << " (basis=" << k << ", comp=" << d << "): "
                   << "Nk=" << Nk << ", gradNk=" << gradNk.transpose()
                   << "\n    dJs:\n" << dJs
                   << "\n    dJg:\n" << dJg
                   << "\n    dJc:\n" << dJc
                   << "\n    dA:\n" << dA
                   << "\n    dE=" << dE_val << "\n";
        }
    }

    gsInfo << "\n=== STEP 5a: Verify A, Ainv, E ===\n";
    gsInfo << "A:\n" << A << "\n";
    gsInfo << "Ainv:\n" << Ainv << "\n";
    gsInfo << "detA = " << detA << ", meas_c = " << meas_c << ", E = " << E << "\n";

    gsInfo << "\n=== STEP 5b: Verify dJc/dalpha and dA/dalpha by FD ===\n";
    for (index_t i = 0; i < nc; ++i)
    {
        gsVector<real_t> cplus = controls;
        cplus(i) += eps;
        domain.setControls(cplus);
        gsMatrix<real_t> sigma_plus_5b;
        domain.eval_into(pt, sigma_plus_5b);
        gsMatrix<real_t> Jsigma_plus_flat_5b;
        domain.deriv_into(pt, Jsigma_plus_flat_5b);
        gsMatrix<real_t> Js_plus_5b = Jsigma_plus_flat_5b.col(0).reshaped(dd, dd).transpose();
        gsMatrix<real_t> Jg_plus_flat_5b;
        geom->deriv_into(sigma_plus_5b, Jg_plus_flat_5b);
        gsMatrix<real_t> Jg_plus_5b = Jg_plus_flat_5b.col(0).reshaped(dd, td).transpose();
        gsMatrix<real_t> Jc_plus_5b = Jg_plus_5b * Js_plus_5b;

        gsMatrix<real_t> dJc_fd = (Jc_plus_5b - Jc) / eps;
        gsMatrix<real_t> A_plus_5b = Jc_plus_5b.transpose() * Jc_plus_5b;
        gsMatrix<real_t> dA_fd = (A_plus_5b - A) / eps;
        real_t E_plus_5b = A_plus_5b.trace() / math::sqrt(A_plus_5b.determinant());
        real_t dE_fd_5b = (E_plus_5b - E) / eps;

        gsInfo << "  control " << i << ":\n"
               << "    dJc (FD):\n" << dJc_fd
               << "\n    dA (FD):\n" << dA_fd
               << "\n    dE (FD) = " << dE_fd_5b << "\n";
    }
    domain.setControls(controls);

    gsVector<real_t> dE_fd(nc);
    dE_fd.setZero();
    for (index_t i = 0; i < nc; ++i)
    {
        gsVector<real_t> cplus = controls;
        gsVector<real_t> cminus = controls;
        cplus(i) += eps;
        cminus(i) -= eps;

        domain.setControls(cplus);
        gsMatrix<real_t> sigma_plus;
        domain.eval_into(pt, sigma_plus);
        gsMatrix<real_t> Jsigma_plus_flat;
        domain.deriv_into(pt, Jsigma_plus_flat);
        gsMatrix<real_t> Js_plus = Jsigma_plus_flat.col(0).reshaped(dd, dd).transpose();
        gsMatrix<real_t> Jg_plus_flat;
        geom->deriv_into(sigma_plus, Jg_plus_flat);
        gsMatrix<real_t> Jg_plus = Jg_plus_flat.col(0).reshaped(dd, td).transpose();
        gsMatrix<real_t> Jc_plus = Jg_plus * Js_plus;
        gsMatrix<real_t> A_plus = Jc_plus.transpose() * Jc_plus;
        real_t E_plus = A_plus.trace() / math::sqrt(A_plus.determinant());

        domain.setControls(cminus);
        gsMatrix<real_t> sigma_minus;
        domain.eval_into(pt, sigma_minus);
        gsMatrix<real_t> Jsigma_minus_flat;
        domain.deriv_into(pt, Jsigma_minus_flat);
        gsMatrix<real_t> Js_minus = Jsigma_minus_flat.col(0).reshaped(dd, dd).transpose();
        gsMatrix<real_t> Jg_minus_flat;
        geom->deriv_into(sigma_minus, Jg_minus_flat);
        gsMatrix<real_t> Jg_minus = Jg_minus_flat.col(0).reshaped(dd, td).transpose();
        gsMatrix<real_t> Jc_minus = Jg_minus * Js_minus;
        gsMatrix<real_t> A_minus = Jc_minus.transpose() * Jc_minus;
        real_t E_minus = A_minus.trace() / math::sqrt(A_minus.determinant());

        dE_fd(i) = (E_plus - E_minus) / (2 * eps);

        if (i == 0) {
            gsInfo << "  FD control 0: E_plus=" << E_plus << ", E_minus=" << E_minus
                   << ", A_plus:\n" << A_plus << "\n"
                   << "  trA_plus=" << A_plus.trace() << ", detA_plus=" << A_plus.determinant()
                   << ", meas_plus=" << math::sqrt(A_plus.determinant()) << "\n";
        }
    }
    domain.setControls(controls);

    gsInfo << "dE/dalpha (analytical): " << dE_an.transpose() << "\n";
    gsInfo << "dE/dalpha (FD):         " << dE_fd.transpose() << "\n";
    gsInfo << "Difference:             " << (dE_an - dE_fd).transpose() << "\n";
    gsInfo << "Relative error:         " << (dE_an - dE_fd).norm() / dE_fd.norm() << "\n\n";

    gsInfo << "=== STEP 6: Full gradient (integral) ===\n";
    gsOptMesh<real_t,MonitorMode::ValueBased> opt(domain,*geom,&tbasis);

    gsVector<real_t> grad(nc);
    gsAsVector<real_t> asgrad(grad.data(), grad.rows());
    opt.gradObj_into(gsAsConstVector<real_t>(controls.data(), controls.size()), asgrad);

    gsVector<real_t> grad_fd(nc);
    for (index_t i = 0; i < nc; ++i)
    {
        gsVector<real_t> cplus = controls;
        gsVector<real_t> cminus = controls;
        cplus(i) += eps;
        cminus(i) -= eps;
        real_t f_plus = opt.evalObj(gsAsConstVector<real_t>(cplus.data(), cplus.size()));
        real_t f_minus = opt.evalObj(gsAsConstVector<real_t>(cminus.data(), cminus.size()));
        grad_fd(i) = (f_plus - f_minus) / (2 * eps);
    }

    gsInfo << "Analytical gradient: " << grad.transpose() << "\n";
    gsInfo << "FD gradient:         " << grad_fd.transpose() << "\n";
    gsInfo << "Difference:          " << (grad - grad_fd).transpose() << "\n";
    gsInfo << "Relative error:      " << (grad - grad_fd).norm() / grad_fd.norm() << "\n";

    return EXIT_SUCCESS;
}
