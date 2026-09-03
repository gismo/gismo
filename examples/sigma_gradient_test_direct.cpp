#include <gismo.h>
#include <gsAssembler/gsExpressions.h>
#include <gsAssembler/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>

using namespace gismo;

template <class T>
class DirectOpt
{
public:
    DirectOpt(gsGeometry<T> & comp, const gsGeometry<T> & geom, const gsBasis<T> & ib)
    : m_comp(comp), m_geom(geom), m_ib(ib), m_mb(ib, true)
    {
        m_cgeom = memory::make_unique(new gsComposedGeometry<T>(m_comp, m_geom));
        m_mp.addPatch(give(m_cgeom));
    }

    T eval(const gsVector<T> & u)
    {
        index_t nb = m_comp.basis().size();
        gsMatrix<T> coefs(nb, 2);
        for (index_t d=0; d<2; ++d)
            for (index_t i=0; i<nb; ++i)
                coefs(i, d) = u(i + d*nb);
        m_comp.setCoefs(coefs);

        gsExprEvaluator<T> evaluator;
        evaluator.setIntegrationElements(m_mb);
        auto G = evaluator.getMap(m_mp);
        auto fform = jac(G).tr() * jac(G);
        return evaluator.integral( fform.trace() / meas(G) );
    }

    void grad(const gsVector<T> & u, gsVector<T> & result)
    {
        index_t nb = m_comp.basis().size();
        gsMatrix<T> coefs(nb, 2);
        for (index_t d=0; d<2; ++d)
            for (index_t i=0; i<nb; ++i)
                coefs(i, d) = u(i + d*nb);
        m_comp.setCoefs(coefs);

        result.resize(u.size());
        result.setZero();

        gsExprEvaluator<T> evaluator;
        evaluator.setIntegrationElements(m_mb);
        auto G = evaluator.getMap(m_mp);
        auto ExprJ = jac(G);
        auto ExprMeas = meas(G);

        const gsBasis<T> & compBasis = m_comp.basis();

        typename gsBasis<T>::domainIter domIt;
        for (unsigned patchInd=0; patchInd < m_mb.nBases(); ++patchInd)
        {
            gsQuadRule<T> QuRule = gsQuadrature::get(m_mb.basis(patchInd), evaluator.options());
            domIt = m_mb.piece(patchInd).makeDomainIterator();
            for (; domIt->good(); domIt->next() )
            {
                gsMatrix<T> uvPoints; gsVector<T> tmpWeights;
                QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), uvPoints, tmpWeights);
                
                gsVector<T> diag = (domIt->upperCorner() - domIt->lowerCorner())/2.0;
                T detJacElement = diag.prod();

                for (index_t p = 0; p!=uvPoints.cols(); ++p)
                {
                    gsVector<T> uv = uvPoints.col(p);
                    gsMatrix<T> J = evaluator.eval(ExprJ, uv);
                    T m = evaluator.eval(ExprMeas, uv)(0);
                    T TrJJ = (J.transpose() * J).trace();
                    T E = TrJJ / m;
                    
                    gsMatrix<T> H(J.rows(), J.cols());
                    T eps_h = 1e-7;
                    for (index_t row=0; row<J.rows(); ++row) {
                        for (index_t col=0; col<J.cols(); ++col) {
                            gsMatrix<T> Jplus = J; Jplus(row, col) += eps_h;
                            T m_plus = std::sqrt((Jplus.transpose()*Jplus).determinant());
                            T E_plus = (Jplus.transpose()*Jplus).trace() / m_plus;
                            gsMatrix<T> Jminus = J; Jminus(row, col) -= eps_h;
                            T m_minus = std::sqrt((Jminus.transpose()*Jminus).determinant());
                            T E_minus = (Jminus.transpose()*Jminus).trace() / m_minus;
                            H(row, col) = (E_plus - E_minus) / (2 * eps_h);
                        }
                    }
                    
                    gsVector<T> sigma = m_comp.eval(uv);
                    gsMatrix<T> A; m_geom.deriv_into(sigma, A);
                    A.resize(m_geom.targetDim(), m_geom.domainDim());
                    gsMatrix<T> d2G; m_geom.deriv2_into(sigma, d2G);
                    gsMatrix<T> Jsigma = m_comp.deriv(uv);
                    Jsigma.resize(2,2);
                    
                    gsMatrix<T> B_vals; compBasis.eval_into(uv, B_vals);
                    gsMatrix<T> B_derivs_flat; compBasis.deriv_into(uv, B_derivs_flat);
                    gsMatrix<T> B_derivs(nb, 2);
                    for (index_t ii=0; ii<nb; ++ii) {
                        B_derivs(ii, 0) = B_derivs_flat(ii*2 + 0, 0);
                        B_derivs(ii, 1) = B_derivs_flat(ii*2 + 1, 0);
                    }
                    
                    for (index_t i = 0; i < nb; ++i)
                    {
                        for (index_t k = 0; k < 2; ++k)
                        {
                            // Numerical deltaJ
                            gsMatrix<T> deltaJ_num(J.rows(), J.cols());
                            T eps_dj = 1e-7;
                            gsVector<T> u_plus = u; u_plus(i + k*nb) += eps_dj;
                            gsVector<T> u_minus = u; u_minus(i + k*nb) -= eps_dj;
                            
                            DirectOpt<T> opt_temp(m_comp, m_geom, m_ib); // wait, expensive
                            // We just need the Jacobian of the composed map at uv
                            // J = (dG/dsigma) * (dsigma/dxi)
                            auto getJ_at = [&](const gsVector<T>& u_vec) {
                                m_comp.setCoefs(gsAsMatrix<T>(const_cast<T*>(u_vec.data()), nb, 2));
                                gsVector<T> sig = m_comp.eval(uv);
                                gsMatrix<T> JacG; m_geom.deriv_into(sig, JacG);
                                JacG.resize(m_geom.targetDim(), m_geom.domainDim());
                                gsMatrix<T> JacS = m_comp.deriv(uv);
                                JacS.resize(2,2);
                                return gsMatrix<T>(JacG * JacS);
                            };
                            
                            deltaJ_num = (getJ_at(u_plus) - getJ_at(u_minus)) / (2 * eps_dj);
                            m_comp.setCoefs(coefs); // reset
                            
                            gsMatrix<T> HessG_k(m_geom.targetDim(), 2);
                            for (index_t j = 0; j < m_geom.targetDim(); ++j) {
                                if (k == 0) {
                                    HessG_k(j, 0) = d2G(j*3 + 0, 0);
                                    HessG_k(j, 1) = d2G(j*3 + 2, 0);
                                } else {
                                    HessG_k(j, 0) = d2G(j*3 + 2, 0);
                                    HessG_k(j, 1) = d2G(j*3 + 1, 0);
                                }
                            }
                            gsMatrix<T> ak = A.col(k);
                            gsMatrix<T> b_grad_i = B_derivs.row(i);
                            b_grad_i.resize(1, 2);
                            gsMatrix<T> deltaJ = B_vals(i, 0) * (HessG_k * Jsigma) + ak * b_grad_i;
                            
                            if (p==0 && i==4 && k==0) {
                                gsInfo << "deltaJ analytical:\n" << deltaJ << "\n";
                                gsInfo << "deltaJ numerical: \n" << deltaJ_num << "\n";
                            }

                            result(i + k * nb) += (H.transpose() * deltaJ).trace() * tmpWeights(p) * detJacElement;
                        }
                    }
                }
            }
        }
    }

private:
    gsGeometry<T> & m_comp;
    const gsGeometry<T> & m_geom;
    const gsBasis<T> & m_ib;
    gsMultiBasis<T> m_mb;
    typename gsGeometry<T>::uPtr m_cgeom;
    gsMultiPatch<T> m_mp;
};

int main()
{
    gsKnotVector<> kv({0,0,0,1,1,1}, 2);
    gsTensorBSplineBasis<2> geomBasis(kv, kv);
    gsMatrix<> geomCoefs(geomBasis.size(), 3);
    geomCoefs.leftCols(2) = geomBasis.anchors().transpose();
    for (index_t i=0; i<geomCoefs.rows(); ++i)
        geomCoefs(i,2) = math::sin(geomCoefs(i,0)*M_PI) * math::cos(geomCoefs(i,1)*M_PI);
    gsGeometry<>::uPtr geom = geomBasis.makeGeometry(geomCoefs);

    gsGeometry<>::uPtr comp = geomBasis.makeGeometry(geomBasis.anchors().transpose());

    gsKnotVector<> kv_int({0,0,0,0.5,1,1,1}, 2);
    gsTensorBSplineBasis<2> intBasis(kv, kv_int);

    DirectOpt<double> opt(*comp, *geom, intBasis);

    index_t nb = comp->basis().size();
    gsVector<double> u(nb * 2);
    for (index_t d=0; d<2; ++d)
        for (index_t i=0; i<nb; ++i)
            u(i + d*nb) = comp->coefs()(i, d);

    u(4) += 0.1;

    gsInfo << "Evaluating Winslow objective (DeltaJ check)...";
    double val = opt.eval(u);
    gsInfo << "Value: " << val << "\n";

    gsVector<double> analytical_grad;
    opt.grad(u, analytical_grad);

    double eps = 1e-6;
    gsVector<double> grad_fd(u.size());
    for (index_t i = 0; i < u.size(); ++i)
    {
        gsVector<double> u_plus = u;
        gsVector<double> u_minus = u;
        u_plus(i) += eps;
        u_minus(i) -= eps;
        grad_fd(i) = (opt.eval(u_plus) - opt.eval(u_minus)) / (2 * eps);
    }

    gsInfo << "Analytical gradient norm: " << analytical_grad.norm() << "\n";
    gsInfo << "FD gradient norm: " << grad_fd.norm() << "\n";
    if (grad_fd.norm() > 1e-12)
        gsInfo << "Relative difference: " << (analytical_grad - grad_fd).norm() / grad_fd.norm() << "\n";

    return 0;
}
