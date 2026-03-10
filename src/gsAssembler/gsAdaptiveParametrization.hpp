/** @file gsSquareDomain.hpp

    @brief Implementation of the gsSquareDomain class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once


#include <gsCore/gsComposedFunction.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsExpressions/gsExpressions.h>
#include <gsExpressions/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsNurbs/gsTensorNurbsBasis.h>
#include <gsUtils/gsStopwatch.h>

namespace gismo
{

template<class T, enum MonitorMode MODE>
gsOptMesh<T,MODE>::gsOptMesh()
{}

template<class T, enum MonitorMode MODE>
gsOptMesh<T,MODE>::gsOptMesh(         gsSquareDomain<T> & composition,
                                const gsGeometry<T> & geometry,
                                const gsBasis<T>    * integrationBasis)
:
gsOptMesh(composition,geometry,nullptr,integrationBasis,false)
{}

template<class T, enum MonitorMode MODE>
gsOptMesh<T,MODE>::gsOptMesh(         gsSquareDomain<T> & composition,
                                const gsGeometry<T> & geometry,
                                const gsFunction<T> * fun,
                                const gsBasis<T>    * integrationBasis,
                                const bool            parametric)
:
m_comp(&composition),
m_geom(&geometry),
m_fun(fun),
m_ib(integrationBasis),
m_mb(*m_ib,true),
m_cgeom(*m_comp,geometry),
// THIS DOES NOT WORK FOR PARAMETRIC=FALSE. POINTER IS LOST WHEN ARRIVING IN evalObj
// m_cfun(parametric ? gsComposedFunction<T>(*m_comp,fun) : gsComposedFunction<T>(m_cgeom,fun)),
m_mp(m_cgeom),
m_parametric(parametric)
{
    m_numDesignVars = m_comp->nControls();
    m_curDesign.resize(m_numDesignVars,1);
    m_controls.resize(m_numDesignVars,1);
    m_controls.col(0) = m_comp->getControls();

    m_options.addReal("Smoothing","Smoothing parameter for the monitor function",0.1);
    m_options.addReal("Penalty","Penalty parameter for the monitor function",1e-2);
}

template<class T, enum MonitorMode MODE>
gsSquareDomain<T> & gsOptMesh<T,MODE>::composition()
{
    return *m_comp;
}

template<class T, enum MonitorMode MODE>
gsOptionList & gsOptMesh<T,MODE>::options()
{
    return m_options;
}


// // ORIGINAL
// template<class T, enum MonitorMode MODE>
// T gsOptMesh<T,MODE>::evalObj(const gsAsConstVector<T> &u) const
// {
//     typedef typename gsExprHelper<T>::geometryMap geometryMap; ///< Geometry map type
//     gsExprEvaluator<T> evaluator;

//     evaluator.setIntegrationElements(m_mb); // does not work when in constructor
//     m_comp->setControls(u);

//     // Penalty constant
//     gsConstantFunction<T> pen(m_options.getReal("Penalty"), m_cgeom.domainDim());
//     geometryMap G = evaluator.getMap(m_mp);
//     auto eps = evaluator.getVariable(pen);

//     // auto chi = 0.5 * (jac(G).det() + pow(eps.val() + pow(jac(G).det(), 2.0), 0.5));
//     // auto invJacMat = jac(G).adj()/chi;
//     // auto eta = evaluator.getVariable(m_fun);
//     // return evaluator.integral( (monitor(eta,G).asDiag()*invJacMat).sqNorm()*meas(G));

//     // gsComposedFunction<T> fun(*m_comp,m_fun.function(0));
//     // gsComposedFunction<T> fun(mpG.patch(0),m_fun.function(0));


//     // gsVector<unsigned> np(m_fun.domainDim());
//     // np.setConstant(m_options.getInt("nSamplingPoints"));
//     // gsMatrix<T> grid = gsPointGrid<T>(m_cgeom.support().col(0),m_cgeom.support().col(1),np);
//     if (m_cgeom.domainDim()==m_cgeom.targetDim())
//     {
//         auto detG = jac(G).det();
//         auto chi = 0.5 * (detG + pow(pow(eps.val(),2.0) + pow(detG, 2.0), 0.5));
//         auto invJacMat = jac(G).adj()/chi; // inverse of jacobian matrix with 'determinant' replaced

//         if (m_fun==nullptr)
//         {
//             auto M = 1/detG;
//             return evaluator.integral( (M*invJacMat).sqNorm()*meas(G));
//         }
//         else
//         {
//             // TEMPORARY FIX: SEE COMMENT IN CONSTRUCTOR
//             m_cfun = m_parametric ? gsComposedFunction<T>(*m_comp,*m_fun) : gsComposedFunction<T>(m_cgeom,*m_fun);
//             auto eta = evaluator.getVariable(m_cfun);
//             auto M   = gismo::expr::monitor<MODE>(eta,G,m_options.getReal("Smoothing"));
//             return evaluator.integral( (M*invJacMat).sqNorm()*meas(G));
//         }
//     }
//     else if (m_cgeom.domainDim()<m_cgeom.targetDim())
//     {
//         auto fform = jac(G).tr()*jac(G);
//         // SQUARE ROOT:
//         auto detG = pow(fform.det().val(),0.5).val(); //jacobian determinant for a surface, i.e. the measure
//         // SQUARED:
//         // auto detG = fform.det().val(); //jacobian determinant for a surface, i.e. the measure

//         // OPTION 1
//         // Compute the chi part
//         // auto chiPPart = eps * ((detG.val() - eps.val()).exp());
//         // Ternary operation to compute chi and chip
//         // auto chi = ternary(eps.val() - detG, chiPPart.val(), detG.val());
//         // OPTION 2
//         // auto chi = 0.5 * (detG + pow(pow(eps.val(),2.0) + pow(detG, 2.0), 0.5));
//         // OPTION 3
//         auto chi = detG;
//         // auto chiPPart = -(-2.0 + pow(eps.val(),2) + eps.val())*pow(detG,3)/pow(eps.val(),4) + (- 3.0 + 2.0*pow(eps.val(),2) + 2.0*eps.val())*pow(detG,2)/pow(eps.val(),3) - detG/eps.val() + 1.0/eps.val();
//         // auto chi = ternary(eps.val() - detG, chiPPart.val(), detG.val());
//         // SQUARED
//         //
//         // auto chi = 0.5 * (pow(detG,0.5) + pow(pow(eps.val(),2.0) + detG, 0.5));

//         // SQUARE ROOT:
//         auto invJacMat = fform.sqrt().adj()/chi; // inverse of jacobian matrix with 'determinant' replaced
//         // SQUARED:
//         // auto invJacMat = fform.adj()/chi; // inverse of jacobian matrix with 'determinant' replaced

//         if (m_fun==nullptr)
//         {
//             auto M = 1/detG;
//             return evaluator.integral( M*(invJacMat).sqNorm()*meas(G));
//         }
//         else
//         {
//             // TEMPORARY FIX: SEE COMMENT IN CONSTRUCTOR
//             m_cfun = m_parametric ? gsComposedFunction<T>(*m_comp,*m_fun) : gsComposedFunction<T>(m_cgeom,*m_fun);
//             auto eta = evaluator.getVariable(m_cfun);
//             auto M   = gismo::expr::monitor<MODE>(eta,G,m_options.getReal("Smoothing"));
//             return evaluator.integral( M*(invJacMat).sqNorm()*meas(G));
//         }
//     }
//     else
//         GISMO_ERROR("Domain dimension must be smaller than or equal to the target dimension, but domainDim = "<<m_cgeom.domainDim()<<" and targetDim = "<<m_cgeom.targetDim());
// }


/// @cond
//
// Objective function evaluation.
//
// Map chain:   hat{Omega} --sigma(u;alpha)--> tilde{Omega} --G--> Omega
//
// Jacobians:   J_sigma (dd x dd),  J_g (td x dd),  J_c = J_g * J_sigma (td x dd)
// Metric:      C = J_c^T J_c  (dd x dd)
// Geometry metric: C_g = J_g^T J_g  (dd x dd)
//
// WITHOUT monitor:
//   E(alpha) = int  tr(C^{-1}) / sqrt(det C)  d hat{Omega}
//
// WITH monitor (weight m^2 = 1/(1 + theta * eta^2)):
//   E(alpha) = int  m^2 * tr(C^{-1}) * sqrt(det C)  d hat{Omega}
//
// where eta^2 depends on MonitorMode:
//
// ValueBased:
//   eta = f(xi(alpha))          (monitor value at moving parametric point)
//   m^2 = 1 / (1 + theta * eta^2)
//
// GradientBased:
//   eta^2 = (nabla_xi f)^T C_g^{-1} (nabla_xi f)
//         = || nabla_x f ||^2   (squared physical gradient norm, for td == dd)
//   m^2 = 1 / (1 + theta * eta^2)
//
//   For m_parametric=true:   nabla_xi f = deriv of f w.r.t. xi (direct)
//   For m_parametric=false:  nabla_xi f = J_g^T * nabla_x f  (chain rule)
//
//   The C_g^{-1} metric accounts for the coordinate distortion so that
//   eta^2 measures the physical gradient norm regardless of the
//   parametric representation.
//
/// @endcond
template<class T, enum MonitorMode MODE>
T gsOptMesh<T,MODE>::evalObj(const gsAsConstVector<T> &u) const
{
    const index_t dd = m_comp->domainDim();
    const index_t td = m_geom->targetDim();
    GISMO_ASSERT(dd <= td, "domainDim must be <= targetDim");

    m_comp->setControls(u);

    const T theta = m_options.getReal("Smoothing");
    const T penalty = m_options.getReal("Penalty");
    const bool hasMonitor = (m_fun != nullptr);

    T result = T(0);

    gsOptionList quadOptions;
    quadOptions.addReal("quA","",1.0);
    quadOptions.addInt("quB","",1);

    auto dom = m_mb.domain();

#   pragma omp parallel
    {
        T thResult = T(0);

        gsFuncData<T> compData, geomData, funData;
        compData.flags = NEED_VALUE | NEED_DERIV;
        geomData.flags = NEED_VALUE | NEED_DERIV;
        if (hasMonitor)
        {
            if (MODE == ValueBased)
                funData.flags = NEED_VALUE;
            else if (MODE == GradientBased)
                funData.flags = NEED_VALUE | NEED_DERIV;
        }

        gsMatrix<T> Js(dd, dd), Jg(td, dd), Jc(td, dd);
        gsMatrix<T> C(dd, dd), Cinv(dd, dd);
        gsMatrix<T> Cg(dd, dd), Cg_inv(dd, dd);
        gsMatrix<T> grad_xi_f(dd, 1);

        gsQuadRule<T> QuRule;
        index_t QuPatch = -1;

        gsMatrix<T> uvPoints;
        gsVector<T> tmpWeights;
        gsMatrix<T> monVals, monDerivs_eval;
        gsMatrix<T> Jgeom_flat, Jsigma_flat;

        for (auto & elem : dom->allElements())
        {
            if (QuPatch != elem.patch())
            {
                QuPatch = elem.patch();
                QuRule = gsQuadrature::get(m_ib->basis(QuPatch), quadOptions);
            }

            QuRule.mapTo(elem.lowerCorner(), elem.upperCorner(),
                         uvPoints, tmpWeights);

            m_comp->compute(uvPoints, compData);
            Jsigma_flat = compData.values[1];
            m_geom->compute(compData.values[0], geomData);
            Jgeom_flat = geomData.values[1];

            if (hasMonitor)
            {
                m_fun->compute((m_parametric ? compData.values[0] : geomData.values[0]), funData);
                monVals = funData.values[0];
                if (MODE == GradientBased)
                    monDerivs_eval = funData.values[1];
            }

            const index_t nPts = uvPoints.cols();
            for (index_t p = 0; p != nPts; ++p)
            {
                Js.noalias() = Jsigma_flat.col(p).reshaped(dd, dd).transpose();
                Jg.noalias() = Jgeom_flat.col(p).reshaped(dd, td).transpose();
                Jc.noalias() = Jg * Js;
                C.noalias() = Jc.transpose() * Jc;
                // Tikhonov shift keeps C^{-1} numerically stable (C = J_c^T J_c
                // is always PSD, but can be ill-conditioned near singularity).
                C.diagonal().array() += penalty;
                Cinv.noalias() = C.inverse();
                // Chi barrier: sign-aware regularization of the Jacobian determinant.
                //
                // Planar case (td==dd): use det(J_c) = det(J_g * J_s) — the composed Jacobian.
                //   chi_c = 0.5*(det_c + sqrt(penalty^2 + det_c^2)) > 0 always.
                //   When det_c < 0 (folded mesh), chi_c -> 0+, so 1/chi_c^2 -> +inf (fold barrier).
                //   No-monitor:   tr(C^{-1}) * |det_c|   / chi_c^2
                //   With-monitor: tr(C^{-1}) * |det_c|^3 / chi_c^2
                //
                // Surface case (td>dd): J_c is rectangular (td x dd), no scalar determinant.
                //   Use det_s = det(J_s) (sigma Jacobian, dd x dd) as the area element.
                //   chi_s = det_s (OPTION 3, no regularization — matches old expression evaluator).
                //   No-monitor:   tr(C^{-1}) / chi_s
                //   With-monitor: m2 * tr(C^{-1}) * chi_s  (old surface behavior)

                T integrand;
                if (td == dd)
                {
                    // Planar case: fold-barrier via composed Jacobian det
                    T det_c = Jc.determinant();
                    T chi_c = T(0.5) * (det_c + math::sqrt(penalty * penalty + det_c * det_c));
                    T abs_det_c = math::abs(det_c);

                    T m2 = T(0);
                    if (hasMonitor)
                    {
                        if (MODE == ValueBased)
                        {
                            T eta = monVals(0, p);
                            m2 = T(1) / (T(1) + theta * eta * eta);
                        }
                        else
                        {
                            Cg.noalias() = Jg.transpose() * Jg;
                            Cg_inv.noalias() = Cg.inverse();

                            if (m_parametric)
                                grad_xi_f.noalias() = monDerivs_eval.col(p);
                            else
                                grad_xi_f.noalias() = Jg.transpose() * monDerivs_eval.col(p).reshaped(td, 1);

                            T eta2 = (grad_xi_f.transpose() * Cg_inv * grad_xi_f)(0, 0);
                            m2 = T(1) / (T(1) + theta * eta2);
                        }
                        // With-monitor integrand (p=3): m2 * tr(C^{-1}) * |det_c|^3 / chi_c^2
                        integrand = m2 * Cinv.trace() * abs_det_c * abs_det_c * abs_det_c / (chi_c * chi_c);
                    }
                    else
                    {
                        // No-monitor integrand (p=1): tr(C^{-1}) * |det_c| / chi_c^2
                        integrand = Cinv.trace() * abs_det_c / (chi_c * chi_c);
                    }
                }
                else
                {
                    // Surface case (td>dd): J_c is rectangular (td x dd), so no scalar det(J_c).
                    // Instead use det_s = det(J_s), the sigma Jacobian determinant, as the area
                    // element in the integral.  Since sigma is a parametrization of the parameter
                    // domain (a bijective map from hat{Omega} to tilde{Omega}), det(J_s) > 0 is
                    // a precondition — not something that needs a fold-barrier.  Therefore chi_s
                    // is just det_s (no regularization, no chi barrier).
                    // No-monitor:   tr(C^{-1}) / chi_s          (old surface behavior)
                    // With-monitor: m2 * tr(C^{-1}) * chi_s     (old surface behavior)
                    T det_s = Js.determinant();
                    T chi_s = det_s;

                    T m2 = T(0);
                    if (hasMonitor)
                    {
                        if (MODE == ValueBased)
                        {
                            T eta = monVals(0, p);
                            m2 = T(1) / (T(1) + theta * eta * eta);
                        }
                        else
                        {
                            Cg.noalias() = Jg.transpose() * Jg;
                            Cg_inv.noalias() = Cg.inverse();

                            if (m_parametric)
                                grad_xi_f.noalias() = monDerivs_eval.col(p);
                            else
                                grad_xi_f.noalias() = Jg.transpose() * monDerivs_eval.col(p).reshaped(td, 1);

                            T eta2 = (grad_xi_f.transpose() * Cg_inv * grad_xi_f)(0, 0);
                            m2 = T(1) / (T(1) + theta * eta2);
                        }
                        // Old surface with-monitor behavior
                        integrand = m2 * Cinv.trace() * chi_s;
                    }
                    else
                    {
                        // Old surface no-monitor behavior
                        integrand = Cinv.trace() / chi_s;
                    }
                }

                thResult += tmpWeights[p] * integrand;
            }
        }
#       pragma omp atomic update
        result += thResult;
    } // omp parallel
    return result;
}


/// @cond
//
// Analytical gradient of the objective function w.r.t. the control
// variables alpha of the composition sigma.
//
// ---------------------------------------------------------------------------
// Notation
// ---------------------------------------------------------------------------
//
// alpha_{k,d}  : d-th coordinate of the k-th basis function coefficient
// N_k          : k-th basis function of sigma, evaluated at quadrature point
// nabla_hat N_k: gradient of N_k w.r.t. hat{u} (the reference coordinates)
// xi = sigma(hat{u}; alpha)  : parametric point in tilde{Omega}
// G(xi)        : geometry map  tilde{Omega} -> Omega
// J_sigma      : Jacobian of sigma  (dd x dd)
// J_g          : Jacobian of G      (td x dd)
// J_c = J_g J_sigma : composed Jacobian (td x dd)
// C = J_c^T J_c     : composed metric   (dd x dd) + penalty*I (Tikhonov shift)
// C_g = J_g^T J_g   : geometry metric   (dd x dd)
// det_c = det(J_c)  : composed Jacobian determinant (square case td==dd)
// chi_c = 0.5*(det_c + sqrt(penalty^2 + det_c^2))  : fold-barrier regularizer
//
// ---------------------------------------------------------------------------
// Kinematic derivatives  (d(...)/d alpha_{k,d})
// ---------------------------------------------------------------------------
//
// d xi_i / d alpha_{k,d}  =  N_k * delta_{id}
//
// d J_sigma / d alpha_{k,d}  =:  dJ_s
//   (dJ_s)_{ij} = delta_{id} * (nabla_hat N_k)_j
//   i.e. only row d is nonzero, equal to nabla_hat N_k transposed
//
// d J_g / d alpha_{k,d}  =:  dJ_g       (td x dd)
//   (dJ_g)_{a,j} = N_k * d^2 G_a / (d xi_j d xi_d)
//   This is the second derivative of the geometry map, contracted with
//   d xi/d alpha = N_k e_d.
//
// d J_c / d alpha_{k,d} = dJ_g * J_s + J_g * dJ_s
// d C   / d alpha_{k,d} = dJ_c^T J_c + J_c^T dJ_c
//
// d(det_c)/d(alpha_{k,d}) = tr(adj(J_c)^T * d(J_c)/d(alpha))
//                         = adj(J_c)^T.col(d) . gradNk
//                         = (J_g^T * adj(J_c)).row(d) . gradNk  [td==dd]
//
// ---------------------------------------------------------------------------
// Objective and integrand
// ---------------------------------------------------------------------------
//
// phi1 = |det_c|   / chi_c^2   (no-monitor exponent p=1)
// phi3 = |det_c|^3 / chi_c^2   (with-monitor exponent p=3)
//
// No-monitor:   E_nm = tr(C^{-1}) * phi1
// With-monitor: E_wm = m^2 * tr(C^{-1}) * phi3
//
// Both integrands blow up as det_c -> 0 (fold), providing a fold barrier.
// This matches the old expression-evaluator code (commit 8c42f2d05).
//
// ---------------------------------------------------------------------------
// Integrand derivative  (no monitor)
// ---------------------------------------------------------------------------
//
// dE_nm = d(tr(C^{-1}))/dalpha * phi1 + tr(C^{-1}) * d(phi1)/dalpha
//       = -trCinvdCCinv * phi1 + tr(C^{-1}) * dphi1_dalpha
//
// where trCinvdCCinv = tr(C^{-1} dC C^{-1})  (= -d(tr(C^{-1}))/dalpha)
//
// d(phi1)/dalpha = phi1 * (1/det_c - 2*dchi_ddet_c/chi_c) * ddet_c/dalpha
//
// ---------------------------------------------------------------------------
// Integrand derivative  (with monitor, both ValueBased & GradientBased)
// ---------------------------------------------------------------------------
//
// dE_wm = dm^2/dalpha * tr(C^{-1}) * phi3
//       + m^2 * (-trCinvdCCinv * phi3 + tr(C^{-1}) * dphi3_dalpha)
//
// d(phi3)/dalpha = phi3 * (3/det_c - 2*dchi_ddet_c/chi_c) * ddet_c/dalpha
//
// dm^2/dalpha = (dm^2/d eta^2) * (d eta^2 / d alpha):
//
// ValueBased:
//   eta = f(xi(alpha))
//   m^2 = 1/(1 + theta eta^2)
//   dm^2/deta = -2 theta eta / (1 + theta eta^2)^2
//   deta/dalpha_{k,d} = (nabla_xi f)_d * N_k
//     (chain rule: d f(xi(alpha))/dalpha = nabla_xi f . d xi/dalpha)
//
// GradientBased:
//   eta^2 = (nabla_xi f)^T C_g^{-1} (nabla_xi f)
//   m^2 = 1/(1 + theta eta^2)
//   dm^2/d(eta^2) = -theta / (1 + theta eta^2)^2
//
//   d(eta^2)/dalpha_{k,d} = term1 + term2, where:
//
//   Term 1 (from C_g changing):
//     Using d(M^{-1}) = -M^{-1} dM M^{-1}, let v = C_g^{-1} nabla_xi f:
//     term1 = -v^T dC_g v
//     with dC_g = dJ_g^T J_g + J_g^T dJ_g
//
//   Term 2 (from nabla_xi f changing at the moving point):
//     term2 = 2 v^T d(nabla_xi f)/dalpha
//
//     For m_parametric=true:
//       d(nabla_xi f)/dalpha_{k,d} = N_k * H_f * e_d
//         (H_f is the dd x dd parametric Hessian of f)
//
//     For m_parametric=false:
//       nabla_xi f = J_g^T nabla_x f,  so differentiating:
//       d(nabla_xi f)/dalpha_{k,d}
//         = dJ_g^T * nabla_x f
//         + N_k * J_g^T * H_f * J_g e_d
//       where H_f is the td x td physical-space Hessian
//       and J_g e_d is the d-th column of J_g (= d G/d xi_d).
//
/// @endcond
template<class T, enum MonitorMode MODE>
void gsOptMesh<T,MODE>::gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
{
    // this->gradObj_FD_into(u, result);
    // return;

    const index_t nc = m_comp->nControls();
    const index_t dd = m_comp->domainDim();
    const index_t td = m_geom->targetDim();
    GISMO_ASSERT(dd <= td, "domainDim must be <= targetDim");

    result.resize(nc);
    result.setZero();
    m_comp->setControls(u);

    const T theta = m_options.getReal("Smoothing");
    const T penalty = m_options.getReal("Penalty");
    const bool hasMonitor = (m_fun != nullptr);

    const gsBasis<T> & sigmaBasis = m_comp->domain().basis();
    const gsDofMapper & mapper = m_comp->mapper();
    const index_t S = dd * (dd + 1) / 2;

    gsOptionList quadOptions;
    quadOptions.addReal("quA","",1.0);
    quadOptions.addInt("quB","",1);

    auto dom = m_mb.domain();

#   pragma omp parallel
    {
        gsVector<T> thResult(nc);
        thResult.setZero();

        gsFuncData<T> geomData, funData, sigmaData, sigmaBasisData;
        geomData.flags = NEED_VALUE | NEED_DERIV | NEED_DERIV2;
        funData.flags = NEED_VALUE;
        if (hasMonitor)
        {
            if (MODE == ValueBased)
                funData.flags |= NEED_DERIV;
            else
                funData.flags |= NEED_DERIV | NEED_DERIV2;
        }
        sigmaData.flags = NEED_VALUE | NEED_DERIV;
        sigmaBasisData.flags = NEED_ACTIVE | NEED_VALUE | NEED_DERIV;

        gsMatrix<T> Js(dd, dd), Jg(td, dd), Jc(td, dd), JcT(dd, td);
        gsMatrix<T> C(dd, dd), Cinv(dd, dd);
        gsMatrix<T> Cg(dd, dd), Cg_inv(dd, dd);
        gsMatrix<T> gradMon(dd, 1), grad_xi_f(dd, 1), grad_x_f(td, 1);
        gsMatrix<T> Hess_f(m_parametric ? dd : td, m_parametric ? dd : td);
        gsVector<T> gradNk(dd), Cinv_gradNk(dd);
        gsVector<T> v(dd);

        gsMatrix<T> D_d(td, dd), A_d(dd, dd), CinvA_d(dd, dd);
        gsVector<T> b_d(dd);
        gsMatrix<T> Cinv_b_all(dd, dd);
        gsVector<T> trCAC(dd);
        gsVector<T> mon_scalar_d(dd);

        gsQuadRule<T> QuRule;
        index_t QuPatch = -1;

        gsMatrix<T> uvPoints;
        gsVector<T> tmpWeights;
        gsMatrix<T> monVals, monDerivs, monDeriv2;
        gsMatrix<T> Jsigma_flat, Jgeom_flat, deriv2_geom;
        gsMatrix<index_t> actives;
        gsMatrix<T> basisVals, basisDerivs;

        for (auto & elem : dom->allElements())
        {
            if (QuPatch != elem.patch())
            {
                QuPatch = elem.patch();
                QuRule = gsQuadrature::get(m_ib->basis(QuPatch), quadOptions);
            }

            QuRule.mapTo(elem.lowerCorner(), elem.upperCorner(),
                         uvPoints, tmpWeights);

            m_comp->compute(uvPoints, sigmaData);
            Jsigma_flat = sigmaData.values[1];

            sigmaBasis.compute(uvPoints, sigmaBasisData);
            actives = sigmaBasisData.actives;
            basisVals = sigmaBasisData.values[0];
            basisDerivs = sigmaBasisData.values[1];

            m_geom->compute(sigmaData.values[0], geomData);
            Jgeom_flat = geomData.values[1];
            deriv2_geom = geomData.values[2];

            if (hasMonitor)
            {
                m_fun->compute((m_parametric ? sigmaData.values[0] : geomData.values[0]), funData);
                monVals = funData.values[0];
                monDerivs = funData.values[1];
                if (MODE == GradientBased)
                    monDeriv2 = funData.values[2];
            }
            const index_t nPts = uvPoints.cols();
            for (index_t p = 0; p != nPts; ++p)
            {
                Js.noalias() = Jsigma_flat.col(p).reshaped(dd, dd).transpose();
                Jg.noalias() = Jgeom_flat.col(p).reshaped(dd, td).transpose();
                Jc.noalias() = Jg * Js;

                C.noalias() = Jc.transpose() * Jc;
                // Tikhonov shift keeps C^{-1} numerically stable (same shift as evalObj).
                C.diagonal().array() += penalty;
                Cinv.noalias() = C.inverse();
                T trCinv = Cinv.trace();

                // Chi barrier: depends on planar (td==dd) vs surface (td>dd).
                //
                // Planar (td==dd): use det(J_c) — the composed Jacobian determinant.
                //   chi_c = 0.5*(det_c + sqrt(penalty^2 + det_c^2))
                //   phi_noMon = |det_c|   / chi_c^2   (no-monitor, exponent p=1)
                //   phi_mon   = |det_c|^3 / chi_c^2   (with-monitor, exponent p=3)
                //   d(phi_p)/d(alpha) uses ddet_c/d(alpha) via adj(J_c)
                //
                // Surface (td>dd): J_c is rectangular — no scalar determinant.
                //   Use det_s = det(J_s) as area element.  Since sigma is a bijective
                //   parametrization of the parameter domain, det_s > 0 is a precondition
                //   (no folding possible) — so chi_s = det_s with no regularization/barrier.
                //   ddet_s/d(alpha) via adj(J_s) and sigma basis derivatives.
                T det_c = T(0), sqrt_reg = T(0), chi_c = T(0), dchi_ddet_c = T(0);
                // phi_noMon = |det_c|^1 / chi_c^2  (no-monitor integrand factor, exponent p=1)
                // phi_mon   = |det_c|^3 / chi_c^2  (with-monitor integrand factor, exponent p=3)
                T abs_det_c = T(0), phi_noMon = T(0), phi_mon = T(0);
                T det_s_val = T(0), chi_s = T(0);
                if (td == dd)
                {
                    det_c     = Jc.determinant();
                    sqrt_reg  = math::sqrt(penalty * penalty + det_c * det_c);
                    chi_c     = T(0.5) * (det_c + sqrt_reg);
                    // d(chi_c)/d(det_c) = 0.5*(1 + det_c/sqrt(penalty^2+det_c^2))
                    dchi_ddet_c = T(0.5) * (T(1) + det_c / sqrt_reg);
                    abs_det_c = math::abs(det_c);
                    phi_noMon = abs_det_c / (chi_c * chi_c);
                    phi_mon   = abs_det_c * abs_det_c * abs_det_c / (chi_c * chi_c);
                }
                else
                {
                    det_s_val = Js.determinant();
                    chi_s     = det_s_val;
                }

                T m2 = T(0), dm2_deta2 = T(0);

                if (hasMonitor)
                {
                    if (MODE == ValueBased)
                    {
                        T eta = monVals(0, p);
                        T denom = T(1) + theta * eta * eta;
                        m2 = T(1) / denom;
                        T dm2_deta = -T(2) * theta * eta / (denom * denom);

                        if (m_parametric)
                            gradMon.noalias() = monDerivs.col(p);
                        else
                            gradMon.noalias() = Jg.transpose() * monDerivs.col(p).reshaped(td, 1);
                        dm2_deta2 = dm2_deta;
                    }
                    else
                    {
                        Cg.noalias() = Jg.transpose() * Jg;
                        Cg_inv.noalias() = Cg.inverse();

                        if (m_parametric)
                        {
                            grad_xi_f.noalias() = monDerivs.col(p);

                            for (index_t i = 0; i != dd; ++i)
                                for (index_t j = i; j != dd; ++j)
                                {
                                    index_t hidx = (i==j) ? i : dd + i*(2*dd-i-3)/2+j-1;
                                    Hess_f(i,j) = monDeriv2(hidx, p);
                                    Hess_f(j,i) = Hess_f(i,j);
                                }
                        }
                        else
                        {
                            grad_x_f.noalias() = monDerivs.col(p).reshaped(td, 1);
                            grad_xi_f.noalias() = Jg.transpose() * grad_x_f;

                            for (index_t i = 0; i != td; ++i)
                                for (index_t j = i; j != td; ++j)
                                {
                                    index_t hidx = (i==j) ? i : td + i*(2*td-i-3)/2+j-1;
                                    Hess_f(i,j) = monDeriv2(hidx, p);
                                    Hess_f(j,i) = Hess_f(i,j);
                                }
                        }

                        T eta2 = (grad_xi_f.transpose() * Cg_inv * grad_xi_f)(0,0);
                        T denom = T(1) + theta * eta2;
                        m2          = T(1) / denom;
                        dm2_deta2   = -theta / (denom * denom);
                    }
                }

                if (hasMonitor && MODE == GradientBased)
                    v.noalias() = Cg_inv * grad_xi_f;

                JcT.noalias() = Jc.transpose();

                // Term1 of d(det_c)/d(alpha) = tr(adj(J_c)^T * D_d * J_s) * N_k  [planar only]
                // Precomputed per d in outer loop; used in inner k loop.
                gsVector<T> trAdjJcDdJs;
                // adj(J_c)^T [pre-computed once per quad point for dd==3]
                // avoids O(nActive * dd) redundant matrix inversions in the inner loops.
                // For dd==3 we compute adj(J_c)^T directly via cofactors (18 mults, cheaper than LU).
                // adj(A)^T(i,j) = (-1)^{i+j} * minor(A,i,j)
                gsMatrix<T> adjJcT_precomp;
                // adj(J_s)^T = det_s * J_s^{-T}  [pre-computed once per quad point for dd==3, surface]
                gsMatrix<T> adjJs_precomp;
                if (td == dd)
                {
                    trAdjJcDdJs.resize(dd);
                    if (dd == 3)
                    {
                        // adj(A)^T(i,j) = (-1)^{i+j} * minor(A, i, j)
                        // minor(A, i, j) = det of A with row i and col j removed.
                        adjJcT_precomp.resize(3, 3);
                        adjJcT_precomp(0,0) =  Jc(1,1)*Jc(2,2) - Jc(1,2)*Jc(2,1); // +M(0,0)
                        adjJcT_precomp(0,1) = -(Jc(1,0)*Jc(2,2) - Jc(1,2)*Jc(2,0)); // -M(0,1)
                        adjJcT_precomp(0,2) =  Jc(1,0)*Jc(2,1) - Jc(1,1)*Jc(2,0); // +M(0,2)
                        adjJcT_precomp(1,0) = -(Jc(0,1)*Jc(2,2) - Jc(0,2)*Jc(2,1)); // -M(1,0)
                        adjJcT_precomp(1,1) =  Jc(0,0)*Jc(2,2) - Jc(0,2)*Jc(2,0); // +M(1,1)
                        adjJcT_precomp(1,2) = -(Jc(0,0)*Jc(2,1) - Jc(0,1)*Jc(2,0)); // -M(1,2)
                        adjJcT_precomp(2,0) =  Jc(0,1)*Jc(1,2) - Jc(0,2)*Jc(1,1); // +M(2,0)
                        adjJcT_precomp(2,1) = -(Jc(0,0)*Jc(1,2) - Jc(0,2)*Jc(1,0)); // -M(2,1)
                        adjJcT_precomp(2,2) =  Jc(0,0)*Jc(1,1) - Jc(0,1)*Jc(1,0); // +M(2,2)
                    }
                    else if (dd > 3)
                        adjJcT_precomp.noalias() = det_c * Jc.inverse().transpose();
                }
                else if (dd == 3)
                {
                    adjJs_precomp.resize(3, 3);
                    adjJs_precomp(0,0) =  Js(1,1)*Js(2,2) - Js(1,2)*Js(2,1); // +M(0,0)
                    adjJs_precomp(0,1) = -(Js(1,0)*Js(2,2) - Js(1,2)*Js(2,0)); // -M(0,1)
                    adjJs_precomp(0,2) =  Js(1,0)*Js(2,1) - Js(1,1)*Js(2,0); // +M(0,2)
                    adjJs_precomp(1,0) = -(Js(0,1)*Js(2,2) - Js(0,2)*Js(2,1)); // -M(1,0)
                    adjJs_precomp(1,1) =  Js(0,0)*Js(2,2) - Js(0,2)*Js(2,0); // +M(1,1)
                    adjJs_precomp(1,2) = -(Js(0,0)*Js(2,1) - Js(0,1)*Js(2,0)); // -M(1,2)
                    adjJs_precomp(2,0) =  Js(0,1)*Js(1,2) - Js(0,2)*Js(1,1); // +M(2,0)
                    adjJs_precomp(2,1) = -(Js(0,0)*Js(1,2) - Js(0,2)*Js(1,0)); // -M(2,1)
                    adjJs_precomp(2,2) =  Js(0,0)*Js(1,1) - Js(0,1)*Js(1,0); // +M(2,2)
                }
                else if (dd > 3)
                    adjJs_precomp.noalias() = det_s_val * Js.inverse().transpose();

                for (index_t d = 0; d != dd; ++d)
                {
                    for (index_t a = 0; a != td; ++a)
                        for (index_t j = 0; j != dd; ++j)
                        {
                            index_t lo = math::min(d, j);
                            index_t hi = math::max(d, j);
                            index_t hess_idx = (lo == hi) ? lo : dd + lo * (2 * dd - lo - 3) / 2 + hi - 1;
                            D_d(a, j) = deriv2_geom(a * S + hess_idx, p);
                        }

                    A_d.noalias() = Js.transpose() * D_d.transpose() * Jc + JcT * D_d * Js;
                    b_d.noalias() = JcT * Jg.col(d);
                    CinvA_d.noalias() = Cinv * A_d;
                    Cinv_b_all.col(d).noalias() = Cinv * b_d;
                    trCAC(d) = (CinvA_d.array() * Cinv.transpose().array()).sum();

                         if (td == dd)
                         {
                             // trAdjJcDdJs(d) = tr(adj(J_c) * D_d * J_s)  [using non-transposed adj]
                             if (dd == 2)
                             {
                                 // 2x2 special case: avoids Jc.inverse() near singularity.
                                 // adj(J_c) = [[Jc(1,1), -Jc(0,1)], [-Jc(1,0), Jc(0,0)]]
                                 // For dd>2 the general branch below is used instead.
                                 gsMatrix<T> adjJc(2, 2);
                                 adjJc << Jc(1,1), -Jc(0,1), -Jc(1,0), Jc(0,0);
                                 trAdjJcDdJs(d) = (adjJc * D_d * Js).trace();
                             }
                             else
                             {
                                 // General: adj(J_c) = det_c * J_c^{-1} = adjJcT_precomp^T
                                 // adjJcT_precomp pre-computed once per quad point (avoids O(dd) inversions).
                                 trAdjJcDdJs(d) = (adjJcT_precomp.transpose() * D_d * Js).trace();
                             }
                         }

                    if (hasMonitor && MODE == GradientBased)
                    {
                        gsMatrix<T> E_d = D_d.transpose() * Jg + Jg.transpose() * D_d;
                        T vEv_d = (v.transpose() * E_d * v)(0,0);

                        if (m_parametric)
                        {
                            mon_scalar_d(d) = -vEv_d + T(2) * (v.transpose() * Hess_f.col(d))(0,0);
                        }
                        else
                        {
                            gsMatrix<T> HfJgd = Hess_f * Jg.col(d);
                            gsMatrix<T> Dg_d = D_d.transpose() * grad_x_f;
                            gsMatrix<T> JtHJd = Jg.transpose() * HfJgd;
                            mon_scalar_d(d) = -vEv_d + T(2) * (v.transpose() * (Dg_d + JtHJd))(0,0);
                        }
                    }
                }

                const index_t nActive = actives.rows();
                for (index_t loc = 0; loc != nActive; ++loc)
                {
                    const index_t k = actives(loc, p);
                    const T Nk = basisVals(loc, p);

                    for (index_t j = 0; j != dd; ++j)
                        gradNk(j) = basisDerivs(loc * dd + j, p);

                    Cinv_gradNk.noalias() = Cinv * gradNk;

                    for (index_t d = 0; d != dd; ++d)
                    {
                        if (!mapper.is_free(k, 0, d))
                            continue;
                        const index_t ii = mapper.index(k, 0, d);

                        T cb_cgN = Cinv_b_all.col(d).dot(Cinv_gradNk);

                        // trCinvdCCinv = d(tr(C^{-1}))/d(alpha_{k,d})
                        //              = -tr(C^{-1} dC C^{-1}) where dC = A_d Nk + 2 b_d gradNk^T (sym)
                        //              = -Nk*tr(C^{-1}A_dC^{-1}) - 2*(C^{-1}b_d).(C^{-1}gradNk)
                        T trCinvdCCinv = Nk * trCAC(d) + T(2) * cb_cgN;

                        // d(det)/d(alpha_{k,d}) and gradient contribution dE
                        T dE;
                        if (td == dd)
                        {
                            // Planar case: det_c = det(J_c), fold-barrier phi_noMon/phi_mon
                            //
                            // d(det_c)/d(alpha_{k,d}): J_c = J_g * J_s, d(J_c)/d(alpha) = J_g.col(d)*gradNk^T
                            //   d(det_c)/d(alpha) = (J_g^T * adj(J_c)).row(d) . gradNk
                            T ddet_c_dalpha;
                            if (dd == 2)
                            {
                                // 2x2 special case: avoids Jc.inverse() near singularity.
                                // Term1: trAdjJcDdJs(d) * Nk
                                // Term2: Jg.col(d) . (adj(J_c)^T * gradNk)
                                //   adj(J_c)^T = [[Jc(1,1), -Jc(1,0)], [-Jc(0,1), Jc(0,0)]]
                                // For dd>2 the general branch below is used instead.
                                gsVector<T> adjJcT_gN(2);
                                adjJcT_gN(0) =  Jc(1,1)*gradNk(0) - Jc(1,0)*gradNk(1);
                                adjJcT_gN(1) = -Jc(0,1)*gradNk(0) + Jc(0,0)*gradNk(1);
                                ddet_c_dalpha = trAdjJcDdJs(d) * Nk + Jg.col(d).dot(adjJcT_gN);
                            }
                            else
                             {
                                 // General square case: adj(J_c)^T = det_c * J_c^{-T}
                                 // Use adjJcT_precomp pre-computed once per quad point (avoids O(nActive*dd) inversions).
                                 T term2 = Jg.col(d).dot(adjJcT_precomp * gradNk);
                                 ddet_c_dalpha = trAdjJcDdJs(d) * Nk + term2;
                             }

                            // d(phi_p)/d(alpha) = phi_p * (p/det_c - 2*dchi/chi) * ddet_c_dalpha
                            T inv_det_c = (det_c != T(0)) ? T(1) / det_c : T(0);
                            T common_factor = -T(2) * dchi_ddet_c / chi_c;
                            T dphi_noMon_dalpha = phi_noMon * (inv_det_c + common_factor) * ddet_c_dalpha;
                            T dphi_mon_dalpha   = phi_mon   * (T(3)*inv_det_c + common_factor) * ddet_c_dalpha;

                            if (hasMonitor)
                            {
                                T dm2_dalpha;
                                if (MODE == ValueBased)
                                    dm2_dalpha = dm2_deta2 * gradMon(d) * Nk;
                                else
                                    dm2_dalpha = dm2_deta2 * Nk * mon_scalar_d(d);

                                // d/d(alpha) [ m2 * tr(C^{-1}) * |det_c|^3 / chi_c^2 ]
                                dE = m2 * (-trCinvdCCinv * phi_mon + trCinv * dphi_mon_dalpha)
                                   + dm2_dalpha * trCinv * phi_mon;
                            }
                            else
                            {
                                // d/d(alpha) [ tr(C^{-1}) * |det_c| / chi_c^2 ]
                                dE = -trCinvdCCinv * phi_noMon + trCinv * dphi_noMon_dalpha;
                            }
                        }
                        else
                        {
                            // Surface case (td>dd): det_s = det(J_s), chi_s = det_s (OPTION 3)
                            //
                            // d(det_s)/d(alpha_{k,d}): J_s is dd x dd, only row d changes.
                            //   d(J_s)_{id,j} = (nabla_hat N_k)_j  (row d of dJ_s = gradNk^T)
                            //   d(det_s)/d(alpha) = adj(J_s)^T.row(d) . gradNk
                            //                    = (adj(J_s)).col(d) . gradNk
                            T ddet_s_dalpha;
                            if (dd == 2)
                            {
                                // 2x2 special case: avoids Js.inverse() near singularity.
                                // adj(J_s) col d for 2x2.
                                // For dd>2 the general branch below is used instead.
                                if (d == 0)
                                    ddet_s_dalpha =  Js(1,1)*gradNk(0) - Js(1,0)*gradNk(1);
                                else
                                    ddet_s_dalpha = -Js(0,1)*gradNk(0) + Js(0,0)*gradNk(1);
                            }
                            else
                             {
                                 // General: adj(J_s)^T = det_s * J_s^{-T}
                                 // Use adjJs_precomp pre-computed once per quad point (avoids O(nActive*dd) inversions).
                                 ddet_s_dalpha = adjJs_precomp.col(d).dot(gradNk);
                             }

                            // phi_s = 1/chi_s (no-monitor) or chi_s (with-monitor)
                            // d/d(alpha) [tr(C^{-1}) / chi_s]
                            //   = -trCinvdCCinv / chi_s - tr(C^{-1}) / chi_s^2 * ddet_s_dalpha
                            // d/d(alpha) [m2 * tr(C^{-1}) * chi_s]
                            //   = m2 * (-trCinvdCCinv * chi_s + trCinv * ddet_s_dalpha)
                            //   + dm2_dalpha * trCinv * chi_s
                            if (hasMonitor)
                            {
                                T dm2_dalpha;
                                if (MODE == ValueBased)
                                    dm2_dalpha = dm2_deta2 * gradMon(d) * Nk;
                                else
                                    dm2_dalpha = dm2_deta2 * Nk * mon_scalar_d(d);

                                dE = m2 * (-trCinvdCCinv * chi_s + trCinv * ddet_s_dalpha)
                                   + dm2_dalpha * trCinv * chi_s;
                            }
                            else
                            {
                                T inv_chi_s = (chi_s != T(0)) ? T(1) / chi_s : T(0);
                                dE = -trCinvdCCinv * inv_chi_s
                                     - trCinv * inv_chi_s * inv_chi_s * ddet_s_dalpha;
                            }
                        }

                        thResult(ii) += tmpWeights[p] * dE;
                    }
                }
            }
        }

#       pragma omp critical
        result += thResult;
    } // omp parallel
}

template<class T, enum MonitorMode MODE>
void gsOptMesh<T,MODE>::gradObj_FD_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
{
    const index_t n = u.rows();
    const T h = T(1e-7);
    result.resize(n);

    gsVector<T> uu = u;
    gsAsVector<T> tmp(uu.data(), n);
    gsAsConstVector<T> ctmp(uu.data(), n);

    for (index_t i = 0; i < n; ++i)
    {
        tmp[i] = u[i] + h;
        const T fp = this->evalObj(ctmp);
        tmp[i] = u[i] - h;
        const T fm = this->evalObj(ctmp);
        tmp[i] = u[i];
        result[i] = (fp - fm) / (T(2) * h);
    }
}

template<class T, enum MonitorMode MODE>
T gsOptMesh<T,MODE>::computeMinJacobian(const gsAsConstVector<T> & u) const
{
    const index_t dd = m_comp->domainDim();
    m_comp->setControls(u);

    T minDet = std::numeric_limits<T>::max();

    gsOptionList quadOptions;
    quadOptions.addReal("quA","",1.0);
    quadOptions.addInt("quB","",1);

    auto dom = m_mb.domain();

    gsFuncData<T> compData;
    compData.flags = NEED_VALUE | NEED_DERIV;

    gsMatrix<T> Js(dd, dd);
    gsQuadRule<T> QuRule;
    index_t QuPatch = -1;

    gsMatrix<T> uvPoints;
    gsVector<T> tmpWeights;

    for (auto & elem : dom->allElements())
    {
        if (QuPatch != elem.patch())
        {
            QuPatch = elem.patch();
            QuRule = gsQuadrature::get(m_ib->basis(QuPatch), quadOptions);
        }

        QuRule.mapTo(elem.lowerCorner(), elem.upperCorner(), uvPoints, tmpWeights);
        m_comp->compute(uvPoints, compData);

        const index_t nPts = uvPoints.cols();
        for (index_t p = 0; p != nPts; ++p)
        {
            Js.noalias() = compData.values[1].col(p).reshaped(dd, dd).transpose();
            T det = Js.determinant();
            if (det < minDet) minDet = det;
        }
    }

    return minDet;
}

/*
TODO:
ADD constructor in .h file
delegate construction to a separate (DIM) templated function
 */
template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsSquareDomain<T> & composition,
                                                                const gsGeometry<T>     & geometry,
                                                                const gsBasis<T>        & integrationBasis,
                                                                      gsOptimizer<T>    & optimizer,
                                                                const bool                parametric)
:
gsAdaptiveParametrization(composition,geometry,nullptr,integrationBasis,optimizer,parametric)
{
}

template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsSquareDomain<T> & composition,
                                                                const gsGeometry<T>     & geometry,
                                                                const gsFunction<T>     & function,
                                                                const gsBasis<T>        & integrationBasis,
                                                                      gsOptimizer<T>    & optimizer,
                                                                const bool                parametric)
:
gsAdaptiveParametrization(composition,geometry,&function,integrationBasis,optimizer,parametric)
{
}

template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsSquareDomain<T> & composition,
                                                                const gsGeometry<T>     & geometry,
                                                                      gsOptimizer<T>    & optimizer,
                                                                const bool                parametric)
:
gsAdaptiveParametrization(composition,geometry,nullptr,geometry.basis(),optimizer,parametric)
{
}

template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsSquareDomain<T> & composition,
                                                                const gsGeometry<T>     & geometry,
                                                                const gsFunction<T>     & function,
                                                                      gsOptimizer<T>    & optimizer,
                                                                const bool                parametric)
:
gsAdaptiveParametrization(composition,geometry,&function,geometry.basis(),optimizer,parametric)
{
}

template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsSquareDomain<T> & composition,
                                                                const gsGeometry<T>     & geometry,
                                                                const gsFunction<T>     * function,
                                                                const gsBasis<T>        & integrationBasis,
                                                                      gsOptimizer<T>    & optimizer,
                                                                const bool                parametric)
:
m_comp(composition),
m_geom(geometry),
m_fun(function),
m_optimizer(optimizer)
{
    const index_t dd = composition.domain().domainDim();
    if (dd == 2)
    {
        const gsTensorBSplineBasis<2,T> * comp_tbasis_ptr = dynamic_cast<const gsTensorBSplineBasis<2,T> *>(&composition.domain().basis());
        GISMO_ASSERT(comp_tbasis_ptr,"The composition must be a tensor B-spline or tensor NURBS basis (2D)");
        const gsTensorBSplineBasis<2,T> & comp_tbasis = *comp_tbasis_ptr;
        if (const gsTensorBSplineBasis<2,T> * tbasis = dynamic_cast<const gsTensorBSplineBasis<2,T> *>(&integrationBasis))
        {
            gsTensorBSplineBasis<2,T> ibasis = makeIntegrationBasis(*tbasis,comp_tbasis);
            m_integrationBasis = memory::make_unique(new gsTensorBSplineBasis<2,T>(ibasis));
        }
        else if (const gsTensorNurbsBasis<2,T> * nbasis = dynamic_cast<const gsTensorNurbsBasis<2,T> *>(&integrationBasis))
        {
            gsTensorNurbsBasis<2,T> ibasis = makeIntegrationBasis(nbasis->source(),comp_tbasis);
            m_integrationBasis = memory::make_unique(new gsTensorNurbsBasis<2,T>(ibasis));
        }
        else
            GISMO_ERROR("The integration basis must be either a tensor B-spline or a tensor NURBS basis (2D)");
    }
    else if (dd == 3)
    {
        const gsTensorBSplineBasis<3,T> * comp_tbasis_ptr = dynamic_cast<const gsTensorBSplineBasis<3,T> *>(&composition.domain().basis());
        GISMO_ASSERT(comp_tbasis_ptr,"The composition must be a tensor B-spline basis (3D)");
        const gsTensorBSplineBasis<3,T> & comp_tbasis = *comp_tbasis_ptr;
        if (const gsTensorBSplineBasis<3,T> * tbasis = dynamic_cast<const gsTensorBSplineBasis<3,T> *>(&integrationBasis))
        {
            gsTensorBSplineBasis<3,T> ibasis = makeIntegrationBasis(*tbasis,comp_tbasis);
            m_integrationBasis = memory::make_unique(new gsTensorBSplineBasis<3,T>(ibasis));
        }
        else if (const gsTensorNurbsBasis<3,T> * nbasis = dynamic_cast<const gsTensorNurbsBasis<3,T> *>(&integrationBasis))
        {
            gsTensorNurbsBasis<3,T> ibasis = makeIntegrationBasis(nbasis->source(),comp_tbasis);
            m_integrationBasis = memory::make_unique(new gsTensorNurbsBasis<3,T>(ibasis));
        }
        else
            GISMO_ERROR("The integration basis must be either a tensor B-spline or a tensor NURBS basis (3D)");
    }
    else
        GISMO_ERROR("Only 2D and 3D composition domains are supported");
    m_optProblem = gsOptMesh<T,MODE>(composition,m_geom,m_fun,m_integrationBasis.get(),parametric);
    this->defaultOptions();
}

template <class T, enum MonitorMode MODE>
template <short_t _d>
gsTensorBSplineBasis<_d,T> gsAdaptiveParametrization<T,MODE>::makeIntegrationBasis(const gsTensorBSplineBasis<_d,T> & basis1,
                                                                                  const gsTensorBSplineBasis<_d,T> & basis2)
{
    gsTensorBSplineBasis<_d,T> ibasis(basis1);
    // Integration basis: parent basis with knots of composition basis inserted.
    // The integration degree is the product of the two degrees, since G(sigma(u))
    // has (Bezier) polynomial degree p_G * p_sigma in u.  This is conservative
    // (over-integrates) but ensures sufficient quadrature accuracy.
    index_t targetDegree;
    for (size_t d = 0; d!=_d; d++)
    {
        // 1. Insert interior knots of composition basis
        for (typename gsKnotVector<T>::uiterator it = std::next(basis2.knots(d).ubegin());
                                                    it!= std::prev(basis2.knots(d).uend());
                                                    ++it)
            {
                if (ibasis.knots(d).has(*it))
                    continue;
                ibasis.insertKnot(*it,d);
            }
        // 2. Increase the degree
        targetDegree = ibasis.degree(d) * basis2.degree(d);
        ibasis.degreeIncrease(targetDegree-ibasis.degree(d),d);

    }
    return ibasis;
}

template <class T, enum MonitorMode MODE>
gsOptionList & gsAdaptiveParametrization<T,MODE>::options()
{
    return m_options;
}

template <class T, enum MonitorMode MODE>
void gsAdaptiveParametrization<T,MODE>::defaultOptions()
{
    m_options.addReal("Penalty","Penalization coefficient for Jacobian determinant [default=0.01]",1e-2);
    m_options.addReal("Smoothing","Smoothing parameter in the monitor function [default=0.1]",0.1);
    m_options.addInt("Mode","0: Relocate based on f [default]; 1: Relocate based on grad(f)",0);
}

template <class T, enum MonitorMode MODE>
void gsAdaptiveParametrization<T,MODE>::solve()
{
    // Set the optimization problem
    m_optProblem.options().setReal("Smoothing",m_options.getReal("Smoothing"));
    m_optProblem.options().setReal("Penalty",m_options.getReal("Penalty"));
    m_optimizer.setProblem(&m_optProblem);
    // Solve the optimization problem
    gsVector<T> controls = m_optProblem.composition().getControls();
    m_optimizer.solve(controls);
    controls = m_optimizer.currentDesign();
    m_optProblem.composition().setControls(controls);
    gsInfo<<"Finished with objective value: "<<m_optimizer.objective()<<std::endl;
}

template<class T, enum MonitorMode MODE>
T gsAdaptiveParametrization<T,MODE>::computeMinJacobian() const
{
    gsVector<T> controls = m_comp.getControls();
    gsAsConstVector<T> u(controls.data(), controls.rows());
    return m_optProblem.computeMinJacobian(u);
}


} // namespace gismo

