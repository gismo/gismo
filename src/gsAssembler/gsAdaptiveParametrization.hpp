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
#include <gsExpressions/gsExpressions.h>
#include <gsExpressions/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsNurbs/gsTensorNurbsBasis.h>
#include <gsUtils/gsPointGrid.h>
#include <iomanip>
#include <sstream>

namespace gismo
{

/// Validates the "quA"/"quB" element-sweep quadrature options shared by
/// gsOptMesh, gsOptFit and gsOptL2: per direction i, gsQuadrature builds
/// nnodes[i] = truncate(quA*deg_i + quB + 0.5) (see gsQuadrature::numNodes()),
/// a rounding, not a ceiling, cast. A negative value, an all-zero pair, or a
/// small-enough positive pair (e.g. quA=0.1, quB=0 on a degree-2 basis) can
/// therefore still truncate to zero nodes in some direction, which silently
/// makes the element sweep visit zero quadrature points and the integrand
/// evaluate to exactly 0 rather than raising a diagnosable error. Checked
/// against the actual integration basis \a ib, so this must run once per
/// gsBasis in play; a multi-patch \a ib whose patches carry different
/// degrees is only checked against patch 0.
template <class T>
inline void checkQuadOptions(const gsOptionList & options, const gsBasis<T> & ib)
{
    const real_t quA = options.getReal("quA");
    const index_t quB = options.getInt("quB");
    GISMO_ENSURE(quA >= 0, "quA must be >= 0, got " << quA);
    GISMO_ENSURE(quB >= 0, "quB must be >= 0, got " << quB);
    GISMO_ENSURE(quA > 0 || quB > 0,
        "quA and quB cannot both be zero (quA=" << quA << ", quB=" << quB
        << "); the element sweep would use zero quadrature points per direction");

    const gsVector<index_t> nn = gsQuadrature::numNodes(ib, quA, quB);
    index_t worstDir = 0;
    for (index_t i = 1; i < nn.size(); ++i)
        if (nn[i] < nn[worstDir]) worstDir = i;
    GISMO_ENSURE(nn[worstDir] >= 1,
        "quA=" << quA << ", quB=" << quB << " give a zero-node quadrature rule "
        "in direction " << worstDir << " on this integration basis (nodes per "
        "direction: " << nn.transpose() << ")");
}

/// True iff \a candidate is admissible as an integration basis for the
/// composition basis \a comp, per direction: same parameter range,
/// degree at least \a comp's, and every interior knot of \a comp
/// present up to tolerance. Used by the pass-through
/// (integrationBasisIsFinal_t) constructor. Coverage note: \ref adaptparam_supermesh.
template <short_t d, class T>
inline bool tensorBasisAdmissible(const gsTensorBSplineBasis<d,T>  & comp,
                          const gsTensorBSplineBasis<d,T>  & candidate)
{
    // Absolute knot tolerance (comp is always a gsSquareDomain's basis, fixed
    // domain [0,1]); shared rationale: \ref adaptparam_supermesh.
    const T tol = 1e-12;
    for (short_t dir = 0; dir != d; ++dir)
    {
        // 1. Parameter range: equality, not containment. A wider integration
        // domain evaluates sigma outside its own parameter domain
        // (extrapolation); a narrower one leaves part of the domain
        // unintegrated. Both are silently wrong.
        const T compFirst = comp.knots(dir).first(), compLast = comp.knots(dir).last();
        const T candFirst = candidate.knots(dir).first(), candLast = candidate.knots(dir).last();
        if (math::abs(candFirst - compFirst) > tol || math::abs(candLast - compLast) > tol)
        {
            // Full precision, not the stream's default (6 significant
            // digits): a mismatch below default precision but above \a tol
            // would otherwise print as two identical intervals.
            std::ostringstream oss;
            oss << std::setprecision(REAL_DIG+1);
            oss << "tensorBasisAdmissible: parameter range mismatch in direction " << dir
                << ": sigma is on [" << compFirst << "," << compLast << "], the supplied "
                << "integration basis is on [" << candFirst << "," << candLast << "]\n";
            gsWarn << oss.str();
            return false;
        }

        // 2. Degree: >=, deliberately not ==, and deliberately not the product
        // degree p_analysis*p_sigma that makeIntegrationBasis() produces --
        // over-integration is safe, and the product degree is unknowable here
        // (this constructor only holds the caller-supplied candidate, not the
        // analysis basis it may have been built from).
        if (candidate.degree(dir) < comp.degree(dir))
        {
            gsWarn << "tensorBasisAdmissible: degree mismatch in direction " << dir
                   << ": sigma has degree " << comp.degree(dir) << ", the supplied "
                   << "integration basis has degree " << candidate.degree(dir) << "\n";
            return false;
        }

        // 3. Interior-knot containment.
        const gsKnotVector<T> & candidateKnots = candidate.knots(dir);
        for (typename gsKnotVector<T>::uiterator it = std::next(comp.knots(dir).ubegin());
                                                    it!= std::prev(comp.knots(dir).uend());
                                                    ++it)
        {
            bool found = false;
            for (typename gsKnotVector<T>::uiterator jt = candidateKnots.ubegin();
                                                        jt!= candidateKnots.uend();
                                                        ++jt)
            {
                if (math::abs(*it - *jt) <= tol)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                // Same precision rationale as the parameter-range warning
                // above: a near-tolerance mismatch is invisible at the
                // stream's default precision.
                std::ostringstream oss;
                oss << std::setprecision(REAL_DIG+1);
                oss << "tensorBasisAdmissible: sigma's interior knot " << *it
                    << " in direction " << dir
                    << " is not present in the supplied integration basis\n";
                gsWarn << oss.str();
                return false;
            }
        }
    }
    return true;
}

/// True iff sigma's knot mesh (\a comp) is exactly a refinement LEVEL of
/// the dyadic hierarchy of \a analysis; then \a sigmaLevel is set to that
/// level, else \a reason receives a human-readable failure cause.
/// Preconditions and the pass-through-constructor relationship:
/// \ref adaptparam_supermesh.
template <short_t d, class T>
inline bool sigmaLevelInHierarchy(const gsTHBSplineBasis<d,T>      & analysis,
                          const gsTensorBSplineBasis<d,T>  & comp,
                          index_t                          & sigmaLevel,
                          std::string                      & reason)
{
    // merge() and the dyadic level-index arithmetic used below both require
    // automatic (non-manual) levels.
    if (analysis.manualLevels())
    {
        reason = "analysis basis has manual levels (merge and level-index arithmetic require automatic dyadic levels)";
        return false;
    }

    const typename gsHTensorBasis<d,T>::tensorBasis & level0 = analysis.tensorLevel(0);

    // Candidate level: per direction, the smallest dyadic refinement of level 0
    // whose element count is >= comp's element count; clamp at 0 so a coarser
    // comp does not request a negative level.
    sigmaLevel = 0;
    for (short_t dir = 0; dir != d; ++dir)
    {
        const index_t N0 = static_cast<index_t>(level0.knots(dir).numElements());
        const index_t Ns = static_cast<index_t>(comp.knots(dir).numElements());
        // N <<= 1 does not advance past 0, so N0 == 0 would spin the loop
        // forever for any Ns > 0. Not reachable for a valid basis
        // (numElements() >= 1 always), but guarded explicitly rather than
        // relied upon.
        if (0 == N0)
        {
            std::ostringstream oss;
            oss << "analysis basis' level-0 tensor mesh has zero elements in direction " << dir;
            reason = oss.str();
            return false;
        }
        index_t k = 0;
        for (index_t N = N0; N < Ns; N <<= 1)
            ++k;
        sigmaLevel = math::max(sigmaLevel, k);
    }

    // Containment: every unique interior knot of comp must appear, up to
    // tolerance, among the unique knots of analysis' tensor level sigmaLevel.
    // A tolerance (not gsKnotVector::has's exact std::binary_search) is
    // required because a knot read from an XML file need not be the exactly
    // representable dyadic value the hierarchy carries internally. An
    // absolute tolerance is scale-appropriate here (rather than a relative
    // one) because \a comp is always a gsSquareDomain's basis, on the fixed
    // interval [0,1]^d -- it would not be on a general parameter domain.
    // Consequence: a knot that matches within \a tol is treated as present,
    // and the union mesh then carries the analysis hierarchy's knot value,
    // not \a comp's -- the element boundary is displaced by O(tol), harmless
    // at 1e-12 but a real substitution.
    const T tol = 1e-12;
    const typename gsHTensorBasis<d,T>::tensorBasis & Lsigma = analysis.tensorLevel(sigmaLevel);
    for (short_t dir = 0; dir != d; ++dir)
    {
        const gsKnotVector<T> & analysisKnots = Lsigma.knots(dir);
        for (typename gsKnotVector<T>::uiterator it = std::next(comp.knots(dir).ubegin());
                                                    it!= std::prev(comp.knots(dir).uend());
                                                    ++it)
        {
            bool found = false;
            for (typename gsKnotVector<T>::uiterator jt = analysisKnots.ubegin();
                                                        jt!= analysisKnots.uend();
                                                        ++jt)
            {
                if (math::abs(*it - *jt) <= tol)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                // Full precision, not the stream's default (6 significant
                // digits): a near-tolerance mismatch would otherwise print an
                // indistinguishable knot value. This text reaches the user
                // directly, via the caller's GISMO_ENSURE appending \a reason.
                std::ostringstream oss;
                oss << std::setprecision(REAL_DIG+1);
                oss << "sigma's interior knot " << *it << " in direction " << dir
                    << " is not present in the analysis hierarchy at level " << sigmaLevel;
                reason = oss.str();
                return false;
            }
        }
    }
    return true;
}

template<class T>
gsOptSigma<T>::gsOptSigma()
: m_domain(nullptr)
{}

template<class T>
gsOptSigma<T>::gsOptSigma(gsSquareDomain<T> & domain)
: m_domain(&domain)
{
    m_numDesignVars = m_domain->nControls();
    m_curDesign.resize(m_numDesignVars,1);
}

template<class T>
gsSquareDomain<T> & gsOptSigma<T>::composition()
{
    return *m_domain;
}

template<class T>
void gsOptSigma<T>::fillCurDesign()
{
    m_curDesign.col(0) = m_domain->getControls();
}

template<class T>
void gsOptSigma<T>::setCtrls(const gsAsConstVector<T> & u) const
{ const gsVector<T> c = u; m_domain->setControls(c); }

template<class T>
T gsCheckSigmaGradient(gsOptProblem<T> & problem, gsSquareDomain<T> & domain,
                       T h, T floor)
{
    gsVector<T> u = domain.getControls();
    gsAsConstVector<T> uMap(u.data(), u.rows());
    gsVector<T> g(u.rows());
    gsAsVector<T> gMap(g.data(), g.rows());
    problem.gradObj_into(uMap, gMap);
    T maxRel = 0;
    for (index_t j = 0; j != u.rows(); j++)
    {
        // Scale the step with the magnitude of the control -- see gradObj_FD_into's
        // matching convention.
        const T hj = h * math::max(T(1), math::abs(u[j]));
        gsVector<T> up = u, um = u;
        up[j] += hj; um[j] -= hj;
        const T fd =
            (problem.evalObj(gsAsConstVector<T>(up.data(), up.rows())) -
             problem.evalObj(gsAsConstVector<T>(um.data(), um.rows()))) / (T(2) * hj);
        const T scale = math::max(math::abs(fd), math::abs(g[j]));
        if (scale > floor)
            maxRel = math::max(maxRel, math::abs(fd - g[j]) / scale);
    }
    domain.setControls(u); // restore
    return maxRel;
}

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
Base(composition),
m_geom(&geometry),
m_fun(fun),
m_ib(integrationBasis),
m_mb(*m_ib,true),
m_cgeom(*m_domain,geometry),
m_mp(m_cgeom),
m_parametric(parametric)
{
    m_controls.resize(m_numDesignVars,1);
    m_controls.col(0) = m_domain->getControls();

    m_options.addReal("Smoothing","Smoothing parameter for the monitor function",0.1);
    m_options.addReal("Penalty","Penalty parameter for the monitor function",1e-2);
    // Quadrature of the element sweeps: nodes per direction =
    // ceil(deg*quA) + quB on each element of the integration basis.  These used
    // to be hardcoded to (1.0, 1), which left the integration MESH as the only
    // way to control quadrature accuracy -- and hence tempted callers to hand in
    // a finely refined analysis basis purely to buy quadrature points.
    m_options.addReal("quA","Quadrature nodes per direction: deg*quA + quB",1.0);
    m_options.addInt ("quB","Quadrature nodes per direction: deg*quA + quB",1);

    // Orientation guard (planar case only, targetDim()==domainDim()): the
    // evalObj()/gradObj_into() fold barrier regularises det(J_c)=det(J_g)*det(J_sigma)
    // via chi(x)=0.5*(x+sqrt(pen^2+x^2)), which SILENTLY treats a globally
    // negative det(J_g) (a negatively oriented input geometry) as one giant
    // fold: chi stays finite and small, so the optimiser sees a consistent but
    // wrong energy landscape instead of an error. Sample det(J_g) on a coarse
    // grid of the geometry's own parameter domain and reject a non-positive
    // sample up front instead.
    if (m_geom->targetDim() == m_geom->domainDim())
    {
        const index_t gdd = m_geom->domainDim();
        const gsMatrix<T> support = m_geom->support();
        gsVector<unsigned> np(gdd);
        np.setConstant(5);
        const gsMatrix<T> pts = gsPointGrid<T>(support.col(0), support.col(1), np);

        for (index_t p = 0; p != pts.cols(); ++p)
        {
            const T detJg = m_geom->jacobian(pts.col(p)).determinant();
            GISMO_ENSURE(detJg > 0,
                "gsOptMesh requires a positively oriented planar geometry: "
                "det(J_g) = " << detJg << " <= 0 at parameter point "
                << pts.col(p).transpose());
        }
    }
}

template<class T, enum MonitorMode MODE>
gsOptionList & gsOptMesh<T,MODE>::options()
{
    return m_options;
}


// Size the persistent per-thread scratch.  Matrices that are written
// ELEMENT-WISE (D_d, gradNk, Hess_f, trA, ...) must be pre-sized here; those
// assigned wholesale via .noalias() resize themselves, but are sized anyway so
// that the very first sweep does not allocate.
template<class T, enum MonitorMode MODE>
void gsOptMesh<T,MODE>::EvalScratch::init(index_t dd, index_t td)
{
    Js.setZero(dd, dd);  Jg.setZero(td, dd);  Jc.setZero(td, dd);
    Cg.setZero(dd, dd);  Cg_inv.setZero(dd, dd);
    grad_xi_f.setZero(dd, 1);
    ready = true;
}

template<class T, enum MonitorMode MODE>
void gsOptMesh<T,MODE>::GradScratch::init(index_t dd, index_t td,
                                          bool parametric, index_t nc)
{
    Js.setZero(dd, dd);  Jg.setZero(td, dd);
    Jc.setZero(td, dd);  JcT.setZero(dd, td);
    Cg.setZero(dd, dd);  Cg_inv.setZero(dd, dd);
    gradMon.setZero(dd, 1);  grad_xi_f.setZero(dd, 1);  grad_x_f.setZero(td, 1);
    Hess_f.setZero(parametric ? dd : td, parametric ? dd : td);
    gradNk.setZero(dd);  v.setZero(dd);
    D_d.setZero(td, dd); DdJs.setZero(td, dd);
    b_d.setZero(dd);     trA.setZero(dd);
    b_all.setZero(dd, dd);
    mon_scalar_d.setZero(dd);
    trAdjJcDdJs.setZero(dd);  trCginvE_d.setZero(dd);
    adjJcT_precomp.setZero(td == dd ? dd : 0, td == dd ? dd : 0);
    adjJc.setZero(dd, dd);
    E_d.setZero(dd, dd);
    HfJgd.setZero(td, 1); Dg_d.setZero(dd, 1); JtHJd.setZero(dd, 1);
    thResult.setZero(nc);
    ready = true;
}

// Energy formulation (no-monitor / with-monitor, per-mode weight): \ref adaptparam_rstep.
template<class T, enum MonitorMode MODE>
T gsOptMesh<T,MODE>::evalObj(const gsAsConstVector<T> &u) const
{
    const index_t dd = m_domain->domainDim();
    const index_t td = m_geom->targetDim();
    GISMO_ASSERT(dd <= td, "domainDim must be <= targetDim");

    m_domain->setControls(u);

    const T theta = m_options.getReal("Smoothing");
    const T penalty = m_options.getReal("Penalty");
    GISMO_ENSURE(penalty > 0, "Penalty must be > 0, got " << penalty);
    GISMO_ENSURE(theta >= 0, "Smoothing must be >= 0, got " << theta);
    const bool hasMonitor = (m_fun != nullptr);

    T result = T(0);

    checkQuadOptions(m_options, m_ib->basis(0));
    gsOptionList quadOptions;
    quadOptions.addReal("quA","",m_options.getReal("quA"));
    quadOptions.addInt ("quB","",m_options.getInt ("quB"));

    auto dom = m_mb.domain();

    // Per-thread partials summed in THREAD-ID order (not as they finish):
    // makes evalObj bit-reproducible for a given thread count, hence the
    // whole R step reproducible.
#   ifdef _OPENMP
    const int nThreads = omp_get_max_threads();
#   else
    const int nThreads = 1;
#   endif
    std::vector<T> partial(nThreads, T(0));
    // Per-thread minimum of (1+theta*f) over the ValueBased monitor: an
    // exception cannot cross the parallel region below (UB / std::terminate),
    // so the violation is recorded here and raised once, after the join.
    std::vector<T> worst(nThreads, std::numeric_limits<T>::infinity());

#   pragma omp parallel
    {
        T thResult = T(0);
        T thWorst = std::numeric_limits<T>::infinity();

        // Persistent per-thread scratch (see gsOptMesh::EvalScratch): buffers,
        // gsFuncData and the Gauss rule are allocated on the first call and
        // reused by every later objective evaluation of this solve.
        EvalScratch & sc = m_evalScratch.mine();
        if (!sc.ready) sc.init(dd, td);

        gsFuncData<T> & compData = sc.compData;
        gsFuncData<T> & geomData = sc.geomData;
        gsFuncData<T> & funData  = sc.funData;
        compData.flags = NEED_VALUE | NEED_DERIV;
        geomData.flags = NEED_VALUE | NEED_DERIV;
        if (hasMonitor)
        {
            if (MODE == ValueBased)
                funData.flags = NEED_VALUE;
            else if (MODE == GradientBased)
                // eta^2 = (grad f)^T Cg^{-1} (grad f) -- the VALUE is not used.
                funData.flags = NEED_DERIV;
        }

        gsMatrix<T> & Js = sc.Js, & Jg = sc.Jg, & Jc = sc.Jc;
        gsMatrix<T> & Cg = sc.Cg, & Cg_inv = sc.Cg_inv;
        gsMatrix<T> & grad_xi_f = sc.grad_xi_f;
        gsMatrix<T> & uvPoints  = sc.uvPoints;
        gsVector<T> & tmpWeights = sc.tmpWeights;

        for (auto & elem : dom->allElements())
        {
            if (sc.QuPatch != elem.patch())
            {
                sc.QuPatch = elem.patch();
                sc.QuRule = gsQuadrature::get(m_ib->basis(sc.QuPatch), quadOptions);
            }

            sc.QuRule.mapTo(elem.lowerCorner(), elem.upperCorner(),
                            uvPoints, tmpWeights);

            m_domain->compute(uvPoints, compData);
            m_geom->compute(compData.values[0], geomData);

            if (hasMonitor)
                m_fun->compute((m_parametric ? compData.values[0] : geomData.values[0]), funData);

            // Alias, do not copy: assigning these into local gsMatrix members
            // would copy every Jacobian block for every element of every
            // objective evaluation.
            const gsMatrix<T> & Jsigma_flat = compData.values[1];
            const gsMatrix<T> & Jgeom_flat  = geomData.values[1];
            // (funData.values is left empty when there is no monitor, so it
            // cannot be indexed; s_noMon keeps the aliases well-formed.)
            static const gsMatrix<T> s_noMon;
            const gsMatrix<T> & monVals =
                (hasMonitor && MODE == ValueBased)    ? funData.values[0] : s_noMon;
            const gsMatrix<T> & monDerivs_eval =
                (hasMonitor && MODE == GradientBased) ? funData.values[1] : s_noMon;

            const index_t nPts = uvPoints.cols();
            for (index_t p = 0; p != nPts; ++p)
            {
                Js.noalias() = Jsigma_flat.col(p).reshaped(dd, dd).transpose();
                Jg.noalias() = Jgeom_flat.col(p).reshaped(dd, td).transpose();
                Jc.noalias() = Jg * Js;

                // T = tr(J_c^T J_c) = ||J_c||_F^2, computed DIRECTLY (not via
                // tr(C^{-1}) with a Tikhonov-shifted C -- that shift breaks
                // the identity tr(C^{-1})*det^2 == T and rewards a collapsing
                // mesh instead of penalising it). Fold barrier via the
                // Garanzha regulariser chi(x) = 0.5*(x + sqrt(eps^2+x^2)) on
                // the regularised composite area element g_reg (planar:
                // chi(det J_c); surface: chi(det J_s)*g_S). Full derivation
                // and case analysis: \ref adaptparam_rgrad.
                const T Tval = Jc.squaredNorm();

                T integrand;
                if (td == dd)
                {
                    // Planar case: fold-barrier via composed Jacobian det
                    T det_c = Jc.determinant();
                    T chi_c = T(0.5) * (det_c + math::sqrt(penalty * penalty + det_c * det_c));

                    T m2 = T(0);
                    if (hasMonitor)
                    {
                        if (MODE == ValueBased)
                        {
                            T eta = monVals(0, p);
                            // value-based weight (\ref adaptparam_rstep): m2 = 1/(1+theta*f).
                            // The precondition 1+theta*f>0 cannot be enforced by GISMO_ENSURE
                            // here (inside #pragma omp parallel, an escaping exception is UB);
                            // instead thWorst records the minimum seen and the aggregate is
                            // checked once after the region joins. A violation here makes m2
                            // negative and omega non-finite below, which the exit isfinite
                            // check at the end of evalObj already turns into a large finite
                            // value for the line search to reject.
                            thWorst = math::min(thWorst, T(1) + theta * eta);
                            m2 = T(1) / (T(1) + theta * eta);
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
                            // gradient-based weight (\ref adaptparam_rstep): m2 = 1/(1+theta*||grad f||^2)
                            m2 = T(1) / (T(1) + theta * eta2);
                        }
                        // With-monitor: omega * T / g_reg,  g_reg = chi_c   (paper E = int omega*T/g)
                        T omega = math::sqrt(m2);
                        integrand = omega * Tval / chi_c;
                    }
                    else
                    {
                        // No-monitor (omega = 1/g): T / g_reg^2,  g_reg = chi_c
                        integrand = Tval / (chi_c * chi_c);
                    }
                }
                else
                {
                    // Surface case (td>dd): J_c is rectangular (td x dd), so no scalar det(J_c).
                    // The composite area element factorises as g_c = |det J_s| * g_S, with
                    //   g_S = sqrt(det(Cg)),  Cg = J_g^T J_g   (the SOURCE surface area element).
                    // Only det_s = det(J_s) can change sign (a fold of sigma), so the barrier
                    // acts on it, and the regularised composite area element is
                    //   g_reg = chi(det_s) * g_S  > 0   (always).
                    T det_s      = Js.determinant();
                    T sqrt_reg_s = math::sqrt(penalty * penalty + det_s * det_s);
                    T chi_s      = T(0.5) * (det_s + sqrt_reg_s);

                    // g_S is needed in BOTH branches now (it enters g_reg), not just with a monitor.
                    // A degenerate geometry point (Cg singular, gS==0, or Cg already NaN from an
                    // upstream fold) is NOT rejected here: it drives g_reg -> 0, the integrand to
                    // overflow, and is caught once by the non-finite exit at the end of evalObj
                    // (see its comment) -- the same policy the no-monitor branch above already
                    // relies on for a degenerate chi_c.
                    Cg.noalias() = Jg.transpose() * Jg;
                    const T gS    = math::sqrt(math::max(Cg.determinant(), T(0)));
                    const T g_reg = chi_s * gS;

                    T m2 = T(0);
                    if (hasMonitor)
                    {
                        if (MODE == ValueBased)
                        {
                            T eta = monVals(0, p);
                            // value-based weight (\ref adaptparam_rstep): m2 = 1/(1+theta*f);
                            // see the planar branch above for why the guard is deferred.
                            thWorst = math::min(thWorst, T(1) + theta * eta);
                            m2 = T(1) / (T(1) + theta * eta);
                        }
                        else
                        {
                            Cg_inv.noalias() = Cg.inverse();

                            if (m_parametric)
                                grad_xi_f.noalias() = monDerivs_eval.col(p);
                            else
                                grad_xi_f.noalias() = Jg.transpose() * monDerivs_eval.col(p).reshaped(td, 1);

                            T eta2 = (grad_xi_f.transpose() * Cg_inv * grad_xi_f)(0, 0);
                            // gradient-based weight (\ref adaptparam_rstep): m2 = 1/(1+theta*||grad f||^2)
                            m2 = T(1) / (T(1) + theta * eta2);
                        }
                        // With-monitor: omega * T / g_reg   (paper E = int omega*T/g)
                        T omega = math::sqrt(m2);
                        integrand = omega * Tval / g_reg;
                    }
                    else
                    {
                        // No-monitor (omega = 1/g): T / g_reg^2
                        integrand = Tval / (g_reg * g_reg);
                    }
                }

                thResult += tmpWeights[p] * integrand;
            }
        }
#       ifdef _OPENMP
        partial[omp_get_thread_num()] = thResult;
        worst[omp_get_thread_num()] = thWorst;
#       else
        partial[0] = thResult;
        worst[0] = thWorst;
#       endif
    } // omp parallel

    for (int i = 0; i != nThreads; ++i)
        result += partial[i];

    T worstAgg = worst[0];
    for (int i = 1; i != nThreads; ++i)
        worstAgg = math::min(worstAgg, worst[i]);
    GISMO_ENSURE(worstAgg > T(0), "ValueBased monitor must satisfy 1+theta*f>0 at every "
        "quadrature point, min(1+theta*f) = " << worstAgg);

    // The line search probes trial designs far from the current one (Moore-Thuente starts
    // at a unit step). A badly folded trial can drive chi -> 0 or C -> singular, so the
    // integrand overflows to Inf/NaN. Handing NaN back to the optimiser is fatal: the
    // cubic interpolation of the line search becomes NaN, it can no longer backtrack, and
    // it returns a ZERO step -- i.e. the "optimisation" silently leaves sigma at the
    // identity even though a perfectly good descent direction exists at a smaller step.
    // Returning a large FINITE value instead lets the line search reject the trial and
    // backtrack normally.
    if (!math::isfinite(result))
        return std::numeric_limits<T>::max() / T(1e6);

    return result;
}


// Full analytic-gradient derivation (kinematic derivatives, integrand
// derivative per mode, chain rule through the monitor weight): \ref adaptparam_rgrad.
template<class T, enum MonitorMode MODE>
void gsOptMesh<T,MODE>::gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
{
    // this->gradObj_FD_into(u, result);
    // return;

    const index_t nc = m_domain->nControls();
    const index_t dd = m_domain->domainDim();
    const index_t td = m_geom->targetDim();
    GISMO_ASSERT(dd <= td, "domainDim must be <= targetDim");

    result.resize(nc);
    result.setZero();
    m_domain->setControls(u);

    const T theta = m_options.getReal("Smoothing");
    const T penalty = m_options.getReal("Penalty");
    GISMO_ENSURE(penalty > 0, "Penalty must be > 0, got " << penalty);
    GISMO_ENSURE(theta >= 0, "Smoothing must be >= 0, got " << theta);
    const bool hasMonitor = (m_fun != nullptr);

    const gsBasis<T> & sigmaBasis = m_domain->domain().basis();
    const gsDofMapper & mapper = m_domain->mapper();
    const index_t S = dd * (dd + 1) / 2;

    checkQuadOptions(m_options, m_ib->basis(0));
    gsOptionList quadOptions;
    quadOptions.addReal("quA","",m_options.getReal("quA"));
    quadOptions.addInt ("quB","",m_options.getInt ("quB"));

    auto dom = m_mb.domain();

    // DETERMINISTIC reduction -- see the comment in evalObj().  "#pragma omp
    // critical result += thResult" was order-dependent in exactly the same way.
#   ifdef _OPENMP
    const int nThreads = omp_get_max_threads();
#   else
    const int nThreads = 1;
#   endif
    std::vector<gsVector<T> > partial(nThreads, gsVector<T>::Zero(nc));
    // Per-thread minimum of (1+theta*f) over the ValueBased monitor -- see the
    // matching comment in evalObj() for why this cannot be a GISMO_ENSURE
    // inside the parallel region.
    std::vector<T> worst(nThreads, std::numeric_limits<T>::infinity());

#   pragma omp parallel
    {
        // Persistent per-thread scratch (see gsOptMesh::GradScratch).
        GradScratch & sc = m_gradScratch.mine();
        if (!sc.ready) sc.init(dd, td, m_parametric, nc);

        gsVector<T> & thResult = sc.thResult;
        thResult.setZero();
        T thWorst = std::numeric_limits<T>::infinity();

        gsFuncData<T> & geomData       = sc.geomData;
        gsFuncData<T> & funData        = sc.funData;
        gsFuncData<T> & sigmaData      = sc.sigmaData;
        gsFuncData<T> & sigmaBasisData = sc.sigmaBasisData;
        geomData.flags = NEED_VALUE | NEED_DERIV | NEED_DERIV2;
        // GradientBased never reads the monitor VALUE (only grad f and Hess f
        // enter eta^2 and its derivative), so do not ask for it -- for an
        // analytic monitor that is a whole extra function evaluation per
        // quadrature point.
        funData.flags = 0;
        if (hasMonitor)
        {
            if (MODE == ValueBased)
                funData.flags = NEED_VALUE | NEED_DERIV;
            else
                funData.flags = NEED_DERIV | NEED_DERIV2;
        }
        sigmaData.flags = NEED_VALUE | NEED_DERIV;
        sigmaBasisData.flags = NEED_ACTIVE | NEED_VALUE | NEED_DERIV;

        gsMatrix<T> & Js = sc.Js, & Jg = sc.Jg, & Jc = sc.Jc, & JcT = sc.JcT;
        gsMatrix<T> & Cg = sc.Cg, & Cg_inv = sc.Cg_inv;
        gsMatrix<T> & gradMon = sc.gradMon, & grad_xi_f = sc.grad_xi_f;
        gsMatrix<T> & grad_x_f = sc.grad_x_f, & Hess_f = sc.Hess_f;
        gsVector<T> & gradNk = sc.gradNk, & v = sc.v;

        // dT/d(alpha_{k,d}) = Nk*tr(A_d) + 2*(b_d.gradNk), where T = tr(C) =
        // ||J_c||_F^2 (see evalObj) and A_d = Js^T D_d^T J_c + J_c^T D_d Js.
        // T is available directly as the squared Frobenius norm of J_c, so this
        // gradient needs no C^{-1} (a dd x dd inversion PER QUADRATURE POINT).
        //
        // A_d's two summands are transposes of each other, so
        //     tr(A_d) = 2*tr(J_c^T D_d Js) = 2 * sum(J_c .* (D_d Js)),
        // i.e. one (td x dd) product and a Frobenius dot instead of forming
        // A_d explicitly and reducing it via two (dd x dd) triple products.
        gsMatrix<T> & D_d = sc.D_d, & DdJs = sc.DdJs;
        gsVector<T> & b_d = sc.b_d, & trA = sc.trA;
        gsMatrix<T> & b_all = sc.b_all;
        gsVector<T> & mon_scalar_d = sc.mon_scalar_d;

        // Member-owned scratch, sized once and reused across the quadrature-point
        // / active-function loops. Allocating them per iteration instead: for
        // the default planar dd=2 case, adjJc and adjJcT_gN alone would cost
        // O(nActive*dd) mallocs per quadrature point -- ~29k allocations per
        // gradient sweep at -r 2.
        gsVector<T> & trAdjJcDdJs = sc.trAdjJcDdJs, & trCginvE_d = sc.trCginvE_d;
        gsMatrix<T> & adjJcT_precomp = sc.adjJcT_precomp, & adjJc = sc.adjJc;
        gsMatrix<T> & E_d = sc.E_d;
        gsMatrix<T> & HfJgd = sc.HfJgd, & Dg_d = sc.Dg_d, & JtHJd = sc.JtHJd;

        gsMatrix<T> & uvPoints   = sc.uvPoints;
        gsVector<T> & tmpWeights = sc.tmpWeights;

        for (auto & elem : dom->allElements())
        {
            if (sc.QuPatch != elem.patch())
            {
                sc.QuPatch = elem.patch();
                sc.QuRule = gsQuadrature::get(m_ib->basis(sc.QuPatch), quadOptions);
            }

            sc.QuRule.mapTo(elem.lowerCorner(), elem.upperCorner(),
                            uvPoints, tmpWeights);

            m_domain->compute(uvPoints, sigmaData);
            sigmaBasis.compute(uvPoints, sigmaBasisData);
            m_geom->compute(sigmaData.values[0], geomData);

            if (hasMonitor)
                m_fun->compute((m_parametric ? sigmaData.values[0] : geomData.values[0]), funData);

            // Alias, do not copy: copying these nine blocks (Jacobian /
            // second-derivative / basis-value) would happen per element, on
            // every gradient evaluation.
            const gsMatrix<T> & Jsigma_flat = sigmaData.values[1];
            const gsMatrix<T> & Jgeom_flat  = geomData.values[1];
            const gsMatrix<T> & deriv2_geom = geomData.values[2];
            const gsMatrix<index_t> & actives = sigmaBasisData.actives;
            const gsMatrix<T> & basisVals   = sigmaBasisData.values[0];
            const gsMatrix<T> & basisDerivs = sigmaBasisData.values[1];
            static const gsMatrix<T> s_noMon;
            const gsMatrix<T> & monVals   = hasMonitor ? funData.values[0] : s_noMon;
            const gsMatrix<T> & monDerivs = hasMonitor ? funData.values[1] : s_noMon;
            const gsMatrix<T> & monDeriv2 =
                (hasMonitor && MODE == GradientBased) ? funData.values[2] : s_noMon;
            const index_t nPts = uvPoints.cols();
            for (index_t p = 0; p != nPts; ++p)
            {
                Js.noalias() = Jsigma_flat.col(p).reshaped(dd, dd).transpose();
                Jg.noalias() = Jgeom_flat.col(p).reshaped(dd, td).transpose();
                Jc.noalias() = Jg * Js;

                // T = tr(C) = ||J_c||_F^2 -- the quantity the CORRECTED integrand uses
                // directly (see evalObj). C/Cinv below are retained only because A_d and
                // b_d are built alongside them; the corrected gradient does not use C^{-1}.
                const T Tval = Jc.squaredNorm();

                // Chi barrier: depends on planar (td==dd) vs surface (td>dd).
                //
                // Planar vs surface fold-barrier case analysis: \ref adaptparam_rgrad.
                // The integrand actually assembled (see evalObj) is omega*T/chi_c
                // resp. T/chi_c^2, so only chi_c and dchi_ddet_c are needed below.
                T det_c = T(0), sqrt_reg = T(0), chi_c = T(0), dchi_ddet_c = T(0);
                if (td == dd)
                {
                    det_c     = Jc.determinant();
                    sqrt_reg  = math::sqrt(penalty * penalty + det_c * det_c);
                    chi_c     = T(0.5) * (det_c + sqrt_reg);
                    // d(chi_c)/d(det_c) = 0.5*(1 + det_c/sqrt(penalty^2+det_c^2))
                    dchi_ddet_c = T(0.5) * (T(1) + det_c / sqrt_reg);
                }
                // Surface case (td>dd): det_s, chi_s, phi_s are computed per-alpha in the
                // inner loop (same Js for all active functions, just wastefully recomputed).

                T m2 = T(0), dm2_deta2 = T(0);

                if (hasMonitor)
                {
                    if (MODE == ValueBased)
                    {
                        T eta = monVals(0, p);
                        // value-based weight (\ref adaptparam_rstep): m2 = 1/(1+theta*f); see
                        // evalObj()'s matching comment for why the guard is deferred to a
                        // per-thread minimum checked after the parallel region joins.
                        T denom = T(1) + theta * eta;
                        thWorst = math::min(thWorst, denom);
                        m2 = T(1) / denom;
                        // d(m2)/d(f) = d(1/(1+theta*f))/d(f) = -theta/(1+theta*f)^2
                        T dm2_deta = -theta / (denom * denom);

                        if (m_parametric)
                            gradMon.noalias() = monDerivs.col(p);
                        else
                            gradMon.noalias() = Jg.transpose() * monDerivs.col(p).reshaped(td, 1);
                        dm2_deta2 = dm2_deta;
                    }
                    else
                    {
                        // No guard on a singular Cg here: a degenerate/near-singular
                        // geometry point yields a non-finite eta2/gradient contribution
                        // that the exit isfinite pass at the end of gradObj_into zeroes.
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
                        // eta2 is the quadratic form (grad f)^T Cg^{-1} (grad f) >= 0, so
                        // denom = 1+theta*eta2 > 0 for theta > 0 -- no guard needed here
                        // (unlike ValueBased, where the raw monitor value f can be negative).
                    }
                }

                if (hasMonitor && MODE == GradientBased)
                    v.noalias() = Cg_inv * grad_xi_f;

                // Surface area element g_S = sqrt(det Cg) (Cg = Jg^T Jg) and Cg^{-1},
                // needed for the paper-exact surface with-monitor area weight (Eq.18).
                // For GradientBased these are already set above; ensure they exist for
                // ValueBased too. Planar (td==dd) does not need g_S.
                // g_S (source surface area element) now enters the REGULARISED area element
                // g_reg = chi(det J_s) * g_S for BOTH the monitor and the no-monitor surface
                // integrand, so it must be available whenever td > dd -- not only with a monitor.
                T gS = T(1);
                // trCginvE_d = tr(Cg^{-1} E_d), E_d = D_d^T Jg + Jg^T D_d (surface only)
                if (td > dd)
                {
                    if (!hasMonitor || MODE == ValueBased)
                    {
                        Cg.noalias()     = Jg.transpose() * Jg;
                        Cg_inv.noalias() = Cg.inverse();
                    }
                    // A degenerate/singular Cg here (gS==0 or NaN) is not rejected: it
                    // propagates to a non-finite gradient component, zeroed by the
                    // exit isfinite pass at the end of gradObj_into (same policy as
                    // evalObj's non-finite exit for the matching degenerate point).
                    gS = math::sqrt(math::max(Cg.determinant(), T(0)));
                    trCginvE_d.setZero(dd);
                }

                JcT.noalias() = Jc.transpose();

                // Term1 of d(det_c)/d(alpha) = tr(adj(J_c)^T * D_d * J_s) * N_k  [planar only]
                // Precomputed per d in outer loop; used in inner k loop.
                // adj(J_c)^T [pre-computed once per quad point for dd>=3]
                // avoids O(nActive * dd) redundant matrix inversions in the inner loops.
                // For dd==3 we compute adj(J_c)^T directly via cofactors (18 mults, cheaper than LU).
                // adj(A)^T(i,j) = (-1)^{i+j} * minor(A,i,j)
                if (td == dd)
                {
                    if (dd == 2)
                    {
                        // adj(J_c) for the 2x2 case, reused by every direction d
                        // (and, transposed, by every active function below).
                        // Avoids Jc.inverse() near singularity.
                        adjJc(0,0) =  Jc(1,1); adjJc(0,1) = -Jc(0,1);
                        adjJc(1,0) = -Jc(1,0); adjJc(1,1) =  Jc(0,0);
                    }
                    else if (dd == 3)
                    {
                        // adj(A)^T(i,j) = (-1)^{i+j} * minor(A, i, j)
                        // minor(A, i, j) = det of A with row i and col j removed.
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
                // Surface case (td>dd): det_s, chi_s, phi_s computed per-(k,d) in the inner loop.

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

                    DdJs.noalias() = D_d * Js;
                    b_d.noalias()  = JcT * Jg.col(d);
                    // tr(A_d) = 2*tr(J_c^T D_d J_s) = 2 * sum(J_c .* (D_d J_s)),
                    // for dT/d(alpha) = Nk*tr(A_d) + 2*(b_d.gradNk); A_d itself is
                    // never formed explicitly (see the note at the temporaries above).
                    trA(d)        = T(2) * (Jc.array() * DdJs.array()).sum();
                    b_all.col(d)  = b_d;

                    if (td == dd)
                    {
                        // trAdjJcDdJs(d) = tr(adj(J_c) * D_d * J_s)  [using non-transposed adj]
                        // adjJc / adjJcT_precomp are pre-computed once per quadrature
                        // point above, so this is one (dd x dd) product per direction.
                        if (dd == 2)
                            trAdjJcDdJs(d) = (adjJc * DdJs).trace();
                        else
                            // General: adj(J_c) = det_c * J_c^{-1} = adjJcT_precomp^T
                            trAdjJcDdJs(d) = (adjJcT_precomp.transpose() * DdJs).trace();
                    }

                    // E_d = D_d^T Jg + Jg^T D_d  (= d(Cg)/d(xi_d) contracted) is needed
                    // by the GradientBased monitor term AND by the surface area-element
                    // derivative; form it ONCE when either wants it.
                    const bool needE_d = (hasMonitor && MODE == GradientBased) || (td > dd);
                    if (needE_d)
                        E_d.noalias() = D_d.transpose() * Jg + Jg.transpose() * D_d;

                    if (hasMonitor && MODE == GradientBased)
                    {
                        T vEv_d = (v.transpose() * E_d * v)(0,0);

                        if (m_parametric)
                        {
                            mon_scalar_d(d) = -vEv_d + T(2) * (v.transpose() * Hess_f.col(d))(0,0);
                        }
                        else
                        {
                            HfJgd.noalias() = Hess_f * Jg.col(d);
                            Dg_d.noalias()  = D_d.transpose() * grad_x_f;
                            JtHJd.noalias() = Jg.transpose() * HfJgd;
                            mon_scalar_d(d) = -vEv_d + T(2) * (v.transpose() * (Dg_d + JtHJd))(0,0);
                        }
                    }

                    // Surface area-element derivative factor (paper-exact Eq.18):
                    //   d(g_S)/d(alpha_{k,d}) = 0.5 * g_S * tr(Cg^{-1} E_d) * N_k
                    if (td > dd)   // needed for d(g_S)/d(alpha) in BOTH monitor and no-monitor
                        trCginvE_d(d) = (Cg_inv * E_d).trace();
                }

                const index_t nActive = actives.rows();
                for (index_t loc = 0; loc != nActive; ++loc)
                {
                    const index_t k = actives(loc, p);
                    const T Nk = basisVals(loc, p);

                    for (index_t j = 0; j != dd; ++j)
                        gradNk(j) = basisDerivs(loc * dd + j, p);

                    for (index_t d = 0; d != dd; ++d)
                    {
                        if (!mapper.is_free(k, 0, d))
                            continue;
                        const index_t ii = mapper.index(k, 0, d);

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
                                // Kept as two scalars rather than a gsVector: this is
                                // the innermost loop (nActive * dd per quadrature
                                // point), where a heap-allocating gsVector would cost
                                // one malloc per iteration.
                                const T adjJcT_gN0 =  Jc(1,1)*gradNk(0) - Jc(1,0)*gradNk(1);
                                const T adjJcT_gN1 = -Jc(0,1)*gradNk(0) + Jc(0,0)*gradNk(1);
                                ddet_c_dalpha = trAdjJcDdJs(d) * Nk
                                              + Jg(0,d)*adjJcT_gN0 + Jg(1,d)*adjJcT_gN1;
                            }
                            else
                             {
                                 // General square case: adj(J_c)^T = det_c * J_c^{-T}
                                 // Use adjJcT_precomp pre-computed once per quad point (avoids O(nActive*dd) inversions).
                                 T term2 = Jg.col(d).dot(adjJcT_precomp * gradNk);
                                 ddet_c_dalpha = trAdjJcDdJs(d) * Nk + term2;
                             }

                            // ---- CORRECTED planar gradient (g_reg = chi_c) -------------------
                            // dT/d(alpha) = tr(dC) = Nk*tr(A_d) + 2*(b_d . gradNk)
                            const T dTval        = Nk * trA(d) + T(2) * b_all.col(d).dot(gradNk);
                            // d(chi_c)/d(alpha) = chi'(det_c) * d(det_c)/d(alpha)
                            const T dchi_c_dalpha = dchi_ddet_c * ddet_c_dalpha;

                            if (hasMonitor)
                            {
                                T dm2_dalpha;
                                if (MODE == ValueBased)
                                    dm2_dalpha = dm2_deta2 * gradMon(d) * Nk;
                                else
                                    dm2_dalpha = dm2_deta2 * Nk * mon_scalar_d(d);

                                // P = omega * T / chi_c
                                // dP = domega*T/chi + omega*dT/chi - omega*T*dchi/chi^2
                                const T omega         = math::sqrt(m2);
                                const T domega_dalpha = dm2_dalpha / (T(2) * omega);
                                dE = domega_dalpha * Tval / chi_c
                                   + omega * dTval / chi_c
                                   - omega * Tval * dchi_c_dalpha / (chi_c * chi_c);
                            }
                            else
                            {
                                // Q = T / chi_c^2
                                // dQ = dT/chi^2 - 2*T*dchi/chi^3
                                dE = dTval / (chi_c * chi_c)
                                   - T(2) * Tval * dchi_c_dalpha / (chi_c * chi_c * chi_c);
                            }
                        }
                        else
                        {
                            // Surface case (td>dd): use det_s = det(J_s) as the signed area element.
                            // Mirrors the planar case with det_s playing the role of det_c.
                            //
                            // E_noMon   = tr(C^{-1}) * |det_s|^2 / chi_s^2   (p=2, E≈T/g^2)
                            // E_withMon = m2 * tr(C^{-1}) * |det_s|^3 / chi_s^2   (p=3, E≈m2*T/g)
                            //   chi_s = 0.5*(det_s + sqrt(pen^2 + det_s^2))  -- fold barrier
                            //   When det_s < 0 (fold), chi_s -> 0+ so energy -> +inf.
                            //
                            // d(det_s)/d(alpha_{k,d}):  J_s changes as dJ_s/dalpha = e_d * gradNk^T
                            //   d(det_s)/dalpha = adj(J_s).row(d) . gradNk
                            //   For dd=2: adj(J_s) = [[Js(1,1),-Js(0,1)],[-Js(1,0),Js(0,0)]]
                            //
                            // phi_noMon = |det_s|^2 / chi_s^2   (p=2, E ≈ T/g^2)
                            // phi_mon   = |det_s|^3 / chi_s^2   (p=3, E ≈ m^2*T/g)
                            // d(phi_p)/d(alpha) = phi_p * (p/det_s - 2*dchi_s/chi_s) * ddet_s_dalpha

                            T det_s     = Js.determinant();
                            T sqrt_reg_s = math::sqrt(penalty * penalty + det_s * det_s);
                            T chi_s     = T(0.5) * (det_s + sqrt_reg_s);
                            T dchi_s_ddet_s = T(0.5) * (T(1) + det_s / sqrt_reg_s);

                            // d(det_s)/d(alpha_{k,d}) = adj(Js)^T.row(d) . gradNk
                            //   adj(Js)^T_{di} = adj(Js)_{id} = cofactor of Js at (i,d)
                            //   For dd=2: adj(Js)^T = [[Js(1,1), -Js(1,0)], [-Js(0,1), Js(0,0)]]
                            //   General: adj(Js)^T = det_s * Js^{-T}, so
                            //     adj(Js)^T.row(d) = det_s * Js^{-T}.row(d) = det_s * Js^{-1}.col(d)
                            T ddet_s_dalpha;
                            if (dd == 2)
                            {
                                // adj(Js)^T row 0: [ Js(1,1), -Js(1,0)]
                                // adj(Js)^T row 1: [-Js(0,1),  Js(0,0)]
                                if (d == 0)
                                    ddet_s_dalpha =  Js(1,1)*gradNk(0) - Js(1,0)*gradNk(1);
                                else
                                    ddet_s_dalpha = -Js(0,1)*gradNk(0) + Js(0,0)*gradNk(1);
                            }
                            else
                            {
                                // General: adj(Js)^T = det_s * Js^{-T}
                                // d(det_s)/dalpha = det_s * Js^{-T}.row(d) . gradNk
                                //                 = det_s * Js^{-1}.col(d) . gradNk
                                ddet_s_dalpha = det_s * Js.inverse().col(d).dot(gradNk);
                            }

                            // ---- CORRECTED surface gradient (g_reg = chi_s * g_S) ------------
                            // dT/d(alpha) = tr(dC) = Nk*tr(A_d) + 2*(b_d . gradNk)
                            const T dTval = Nk * trA(d) + T(2) * b_all.col(d).dot(gradNk);
                            // d(g_S)/d(alpha_{k,d}) = 0.5 * g_S * tr(Cg^{-1} E_d) * N_k
                            const T dgS_dalpha = T(0.5) * gS * trCginvE_d(d) * Nk;
                            // g_reg = chi_s * g_S  =>  d(g_reg) = chi'(det_s)*ddet_s*g_S + chi_s*dg_S
                            const T g_reg      = chi_s * gS;
                            const T dg_reg     = dchi_s_ddet_s * ddet_s_dalpha * gS
                                               + chi_s * dgS_dalpha;

                            if (hasMonitor)
                            {
                                T dm2_dalpha;
                                if (MODE == ValueBased)
                                    dm2_dalpha = dm2_deta2 * gradMon(d) * Nk;
                                else
                                    dm2_dalpha = dm2_deta2 * Nk * mon_scalar_d(d);

                                // P = omega * T / g_reg
                                // dP = domega*T/g_reg + omega*dT/g_reg - omega*T*dg_reg/g_reg^2
                                const T omega         = math::sqrt(m2);
                                const T domega_dalpha = dm2_dalpha / (T(2) * omega);
                                dE = domega_dalpha * Tval / g_reg
                                   + omega * dTval / g_reg
                                   - omega * Tval * dg_reg / (g_reg * g_reg);
                            }
                            else
                            {
                                // Q = T / g_reg^2  =>  dQ = dT/g_reg^2 - 2*T*dg_reg/g_reg^3
                                dE = dTval / (g_reg * g_reg)
                                   - T(2) * Tval * dg_reg / (g_reg * g_reg * g_reg);
                            }
                        }

                        thResult(ii) += tmpWeights[p] * dE;
                    }
                }
            }
        }

#       ifdef _OPENMP
        partial[omp_get_thread_num()] = thResult;
        worst[omp_get_thread_num()] = thWorst;
#       else
        partial[0] = thResult;
        worst[0] = thWorst;
#       endif
    } // omp parallel

    for (int i = 0; i != nThreads; ++i)
        result += partial[i];

    T worstAgg = worst[0];
    for (int i = 1; i != nThreads; ++i)
        worstAgg = math::min(worstAgg, worst[i]);
    GISMO_ENSURE(worstAgg > T(0), "ValueBased monitor must satisfy 1+theta*f>0 at every "
        "quadrature point, min(1+theta*f) = " << worstAgg);

    // Mirrors evalObj()'s non-finite policy: a badly folded trial design can
    // drive a component to Inf/NaN. Zeroing it lets the objective's own
    // isfinite check (which returns a large finite value) make the line
    // search back off, instead of poisoning the whole descent direction with
    // a single non-finite entry.
    for (index_t i = 0; i != result.rows(); ++i)
        if (!math::isfinite(result[i]))
            result[i] = T(0);
}

template<class T, enum MonitorMode MODE>
void gsOptMesh<T,MODE>::gradObj_FD_into( const gsAsConstVector<T> & u, gsAsVector<T> & result, T h0) const
{
    const index_t n = u.rows();
    result.resize(n);

    gsVector<T> uu = u;
    gsAsVector<T> tmp(uu.data(), n);
    gsAsConstVector<T> ctmp(uu.data(), n);

    for (index_t i = 0; i < n; ++i)
    {
        // Scale the step with the magnitude of the control, as gsCheckSigmaGradient does.
        const T h = h0 * math::max(T(1), math::abs(u[i]));
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
    const index_t dd = m_domain->domainDim();
    m_domain->setControls(u);

    T minDet = std::numeric_limits<T>::max();

    checkQuadOptions(m_options, m_ib->basis(0));
    gsOptionList quadOptions;
    quadOptions.addReal("quA","",m_options.getReal("quA"));
    quadOptions.addInt ("quB","",m_options.getInt ("quB"));

    auto dom = m_mb.domain();

    // OpenMP element sweep (mirrors evalObj): gsDomain::allElements() hands each
    // thread its own element chunk, m_domain->compute is read-only after
    // setControls (already called concurrently in evalObj). The per-element
    // minima are combined by a min-reduction, which — unlike a sum — is
    // order-independent and exact, so the result is bit-identical to a serial
    // sweep.
#   pragma omp parallel
    {
        T thMin = std::numeric_limits<T>::max();

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
            m_domain->compute(uvPoints, compData);

            const index_t nPts = uvPoints.cols();
            for (index_t p = 0; p != nPts; ++p)
            {
                Js.noalias() = compData.values[1].col(p).reshaped(dd, dd).transpose();
                const T det = Js.determinant();
                if (det < thMin) thMin = det;
            }
        }
#       pragma omp critical (gsOptMesh_computeMinJacobian)
        if (thMin < minDet) minDet = thMin;
    } // omp parallel

    return minDet;
}

// ---------------------------------------------------------------------------
// gsOptFit: rebased on gsOptSigma. See the gsOptFit doxygen block in
// gsAdaptiveParametrization.h for the derivation; det J_sigma is piecewise
// polynomial on sigma's own knot mesh, so a Gauss rule there resolves it
// exactly. The fold-barrier/box terms themselves live in m_barrier
// (gsFoldBarrier::addObj/addGrad); this class adds only the point-cloud data
// term.
// ---------------------------------------------------------------------------

template<class T>
gsOptFit<T>::gsOptFit(      gsSquareDomain<T> & domain,
                       const gsGeometry<T>     & S,
                       const gsMatrix<T>       & uv,
                       const gsMatrix<T>       & xyz,
                             T mu, T eps, gsFoldBarrierMode mode, index_t quB)
: Base(domain), m_S(&S), m_uv(uv), m_xyz(xyz)
{
    GISMO_ENSURE(domain.domainDim()==2 && domain.targetDim()==2,
                 "gsOptFit is implemented for a planar (domainDim==targetDim==2) "
                 "sigma domain only");
    // The gsFoldBarrier constructor forwards quB straight into
    // gsQuadrature::get() in Sampled mode, which segfaults (does not throw)
    // for a negative order -- validate here too, BEFORE constructing
    // m_barrier, with GISMO_ENSURE (not GISMO_ASSERT) so the check survives
    // -DNDEBUG release builds, which is exactly where this argument is
    // reachable from unchecked user input (--dirQuB on the fitting driver).
    GISMO_ENSURE(quB >= 0, "gsOptFit: barrier quB must be >= 0 (got " << quB << ")");
    m_barrier = gsFoldBarrier<T>(domain, mu, eps, mode, quB);

    this->fillCurDesign();
    m_colloc = domain.domain().basis().collocationMatrix(m_uv); // N x nb

    // Cache sigma's own basis evaluation (active functions + basis values)
    // at the FIXED fitting points m_uv. Only sigma's control coefficients
    // change per iteration, not these points, so the expensive part of
    // eval_into -- element lookup and basis evaluation, redone from scratch
    // on every call otherwise -- needs to happen only once here. Per-iteration
    // evaluation then reduces to gsBasis::linearCombination_into, a sparse
    // weighted sum with the CURRENT coefficients (see evalObj/gradObj_into).
    {
        const gsBasis<T> & sb = domain.domain().basis();
        gsFuncData<T> uvData(NEED_VALUE | NEED_ACTIVE);
        sb.compute(m_uv, uvData);
        m_uvActives = uvData.actives;
        m_uvVals    = uvData.values[0];
    }
}

template<class T>
T gsOptFit<T>::evalObj(const gsAsConstVector<T> & u) const
{
    // m_barrier.addObj() reads m_domain->domain().coefs() and getControls()
    // directly (no design-vector argument -- see the gsFoldBarrier class
    // doc), so it needs sigma's controls in sync with u just like the data
    // term below does.
    this->setCtrls(u);

    gsMatrix<T> xi, vals;
    gsBasis<T>::linearCombination_into(m_domain->domain().coefs(), m_uvActives, m_uvVals, xi);
    xi = xi.cwiseMax(T(0)).cwiseMin(T(1));
    m_S->eval_into(xi, vals);

    const index_t N = m_uv.cols();
    T E = (vals - m_xyz).squaredNorm() / T(N);

    m_barrier.addObj(E);

    if (!math::isfinite(E))
        return std::numeric_limits<T>::max() / T(1e6);

    return E;
}

template<class T>
void gsOptFit<T>::gradObj_into(const gsAsConstVector<T> & u, gsAsVector<T> & result) const
{
    this->setCtrls(u);

    const index_t N = m_uv.cols();
    gsMatrix<T> xi, vals, dS;
    gsBasis<T>::linearCombination_into(m_domain->domain().coefs(), m_uvActives, m_uvVals, xi);
    xi = xi.cwiseMax(T(0)).cwiseMin(T(1));
    // ONE compute() call gets value AND derivative together (shared active
    // function lookup), instead of eval_into then deriv_into separately
    // redoing that lookup at the same N points.
    gsFuncData<T> sData(NEED_VALUE | NEED_DERIV | NEED_ACTIVE);
    m_S->compute(xi, sData);
    vals = sData.values[0];
    dS   = sData.values[1]; // (targetDim*2) x N, row c*2+j: dS_c/dxi_j

    const short_t d = m_S->targetDim();
    gsMatrix<T> W(N, 2); // W(i,j) = 2/N (r_i^T J_S)_j
    for (index_t i = 0; i != N; i++)
        for (short_t j = 0; j != 2; j++)
        {
            T s = 0;
            for (short_t c = 0; c != d; c++)
                s += (vals(c, i) - m_xyz(c, i)) * dS(c * 2 + j, i);
            W(i, j) = T(2) * s / T(N);
        }
    gsMatrix<T> g = m_colloc.transpose() * W; // nb x 2, g(k,j) per coefficient
    result.setZero();
    const gsDofMapper & mapper = m_domain->mapper();
    for (index_t k = 0; k != g.rows(); k++)
        for (short_t j = 0; j != 2; j++)
            if (mapper.is_free(k, 0, j))
                result[mapper.index(k, 0, j)] += g(k, j);

    m_barrier.addGrad(result);

    // Mirrors evalObj()'s non-finite policy (see gsOptMesh::evalObj's comment).
    for (index_t i = 0; i != result.rows(); ++i)
        if (!math::isfinite(result[i]))
            result[i] = T(0);
}

// ---------------------------------------------------------------------------
// gsOptL2: L2-projection-error objective with frozen analysis-space
// coefficients. See the class doc in gsAdaptiveParametrization.h for the
// derivation; this is new (no examples-level source to port).
// ---------------------------------------------------------------------------

template<class T>
gsOptL2<T>::gsOptL2(      gsSquareDomain<T> & domain,
                     const gsGeometry<T>     & geometry,
                     const gsFunction<T>     & solution,
                     const gsFunction<T>     & fun,
                     const gsBasis<T>        * integrationBasis,
                           T mu, T eps, gsFoldBarrierMode mode, index_t quB,
                     const bool                parametric)
: Base(domain), m_geom(&geometry), m_solution(&solution), m_fun(&fun),
  m_ib(integrationBasis), m_mb(*m_ib,true), m_parametric(parametric)
{
    // GISMO_ASSERT is INERT in build_rel; the dimension contract below must
    // hold at run time (fun's dimension mismatch is otherwise a silent
    // garbage read of funData.values[1] in gradObj_into's chain-rule term),
    // so it is GISMO_ENSURE -- same reasoning as the quB check right below.
    GISMO_ENSURE(domain.domainDim()==2 && geometry.domainDim()==2,
                 "gsOptL2: sigma domain and geometry must both have "
                 "domainDim()==2 (sigma is the unit square)");
    GISMO_ENSURE(geometry.targetDim() >= geometry.domainDim(),
                 "gsOptL2: geometry.targetDim() must be >= geometry.domainDim() "
                 "(planar ==2, surface >2), got targetDim()="
                 << geometry.targetDim());
    GISMO_ENSURE(parametric ? (fun.domainDim() == 2)
                             : (fun.domainDim() == geometry.targetDim()),
                 "gsOptL2: fun.domainDim() (" << fun.domainDim() << ") must "
                 "equal 2 when parametric==true, or geometry.targetDim() ("
                 << geometry.targetDim() << ") when parametric==false");
    // See gsOptFit::gsOptFit() for why this is GISMO_ENSURE, not
    // GISMO_ASSERT, and why it must run BEFORE m_barrier is constructed:
    // the gsFoldBarrier constructor forwards quB straight into
    // gsQuadrature::get() in Sampled mode, which segfaults for a negative
    // value, in release builds too.
    GISMO_ENSURE(quB >= 0, "gsOptL2: barrier quB must be >= 0 (got " << quB << ")");
    m_barrier = gsFoldBarrier<T>(domain, mu, eps, mode, quB);

    this->fillCurDesign();

    // Quadrature of the element sweep -- deliberately separate option names
    // from the barrier's own quB (constructor argument above): this quA/quB
    // controls the L2-error integral over m_ib, the barrier's controls the
    // fold quadrature on sigma's own knot mesh. Mirrors gsOptMesh's
    // constructor (gsAdaptiveParametrization.hpp, gsOptMesh ctor).
    m_options.addReal("quA","Quadrature nodes per direction: deg*quA + quB (element sweep)",1.0);
    m_options.addInt ("quB","Quadrature nodes per direction: deg*quA + quB (element sweep)",1);
}

template<class T>
void gsOptL2<T>::EvalScratch::init(index_t dd, index_t tdS)
{
    Js.setZero(dd, dd);
    JS.setZero(tdS, dd);
    Cg.setZero(dd, dd);
    ready = true;
}

template<class T>
void gsOptL2<T>::GradScratch::init(index_t dd, index_t tdS, index_t td_sol, index_t nc)
{
    Js.setZero(dd, dd);
    JS.setZero(tdS, dd);
    Cg.setZero(dd, dd);      Cg_inv.setZero(dd, dd);
    adjJS.setZero(dd, dd);
    D_d.setZero(tdS, dd);
    E_d.setZero(dd, dd);
    r.setZero(td_sol);
    dFdxi.setZero(td_sol, dd);
    thResult.setZero(nc);
    ready = true;
}

template<class T>
T gsOptL2<T>::evalObj(const gsAsConstVector<T> & u) const
{
    // m_barrier.addObj() reads m_domain's CURRENT coefficients/controls
    // directly (see the gsFoldBarrier class doc), so sync first, same as
    // the data term below needs.
    this->setCtrls(u);

    const short_t dd = 2, tdS = m_geom->targetDim();
    const short_t td_sol = m_solution->targetDim();

    checkQuadOptions(m_options, m_ib->basis(0));
    gsOptionList quadOptions;
    quadOptions.addReal("quA","",m_options.getReal("quA"));
    quadOptions.addInt ("quB","",m_options.getInt ("quB"));

    auto dom = m_mb.domain();

    // DETERMINISTIC reduction -- see the comment in gsOptMesh::evalObj():
    // per-thread partials are summed in THREAD-ID order, not completion
    // order, so the objective (and hence anything built on top of it, such
    // as an L-BFGS iteration count) stays reproducible for a given thread
    // count.
#   ifdef _OPENMP
    const int nThreads = omp_get_max_threads();
#   else
    const int nThreads = 1;
#   endif
    std::vector<T> partial(nThreads, T(0));

#   pragma omp parallel
    {
        T thResult = T(0);

        // Persistent per-thread scratch (see gsOptL2::EvalScratch).
        EvalScratch & sc = m_evalScratch.mine();
        if (!sc.ready) sc.init(dd, tdS);

        gsFuncData<T> & compData = sc.compData;
        gsFuncData<T> & geomData = sc.geomData;
        gsFuncData<T> & funData  = sc.funData;
        compData.flags = NEED_VALUE | NEED_DERIV;
        geomData.flags = NEED_VALUE | NEED_DERIV;
        funData.flags  = NEED_VALUE;

        gsMatrix<T> & solVals = sc.solVals;
        gsMatrix<T> & Js = sc.Js, & JS = sc.JS, & Cg = sc.Cg;
        gsMatrix<T> & uvPoints   = sc.uvPoints;
        gsVector<T> & tmpWeights = sc.tmpWeights;

        for (auto & elem : dom->allElements())
        {
            if (sc.QuPatch != elem.patch())
            {
                sc.QuPatch = elem.patch();
                sc.QuRule = gsQuadrature::get(m_ib->basis(sc.QuPatch), quadOptions);
            }
            sc.QuRule.mapTo(elem.lowerCorner(), elem.upperCorner(), uvPoints, tmpWeights);

            m_domain->compute(uvPoints, compData);           // xi = sigma(v), J_sigma
            m_geom->compute(compData.values[0], geomData);   // S(xi), J_S(xi)
            m_solution->eval_into(uvPoints, solVals);         // u_h(v) -- does NOT depend on alpha
            m_fun->compute(m_parametric ? compData.values[0] : geomData.values[0], funData);

            const index_t nPts = uvPoints.cols();
            for (index_t p = 0; p != nPts; ++p)
            {
                Js.noalias() = compData.values[1].col(p).reshaped(dd, dd).transpose();
                JS.noalias() = geomData.values[1].col(p).reshaped(dd, tdS).transpose();
                const T detJs = Js.determinant();

                // g_S: planar keeps the SIGNED det(J_S); surface uses the
                // non-negative area element sqrt(det(J_S^T J_S)) (Cg = J_S^T J_S).
                // det(J_sigma) stays signed in both cases -- the fold barrier is
                // what keeps it positive, not this objective.
                T gS;
                if (tdS == dd)
                    gS = JS.determinant();
                else
                {
                    Cg.noalias() = JS.transpose() * JS;
                    gS = math::sqrt(math::max(Cg.determinant(), T(0)));
                }
                const T w = gS * detJs;

                T r2 = T(0);
                for (short_t c = 0; c != td_sol; ++c)
                {
                    const T rc = solVals(c, p) - funData.values[0](c, p);
                    r2 += rc * rc;
                }

                thResult += tmpWeights[p] * r2 * w;
            }
        }

#       ifdef _OPENMP
        partial[omp_get_thread_num()] = thResult;
#       else
        partial[0] = thResult;
#       endif
    } // omp parallel

    T E = T(0);
    for (int i = 0; i != nThreads; ++i)
        E += partial[i];

    m_barrier.addObj(E);

    if (!math::isfinite(E))
        return std::numeric_limits<T>::max() / T(1e6);

    return E;
}

template<class T>
void gsOptL2<T>::gradObj_into(const gsAsConstVector<T> & u, gsAsVector<T> & result) const
{
    const index_t nc = m_domain->nControls();
    result.resize(nc);
    result.setZero();
    this->setCtrls(u);

    const short_t dd = 2, tdS = m_geom->targetDim();
    const short_t S2 = dd * (dd + 1) / 2; // = 3, Hessian-triplet stride
    const short_t td_sol = m_solution->targetDim();

    const gsBasis<T> & sigmaBasis = m_domain->domain().basis();
    const gsDofMapper & mapper = m_domain->mapper();

    checkQuadOptions(m_options, m_ib->basis(0));
    gsOptionList quadOptions;
    quadOptions.addReal("quA","",m_options.getReal("quA"));
    quadOptions.addInt ("quB","",m_options.getInt ("quB"));

    auto dom = m_mb.domain();

    // DETERMINISTIC reduction -- see the comment in evalObj(). The only
    // shared write is the per-control result[]; each thread accumulates it
    // into its own scratch vector, and the partials are summed in
    // THREAD-ID order afterwards.
#   ifdef _OPENMP
    const int nThreads = omp_get_max_threads();
#   else
    const int nThreads = 1;
#   endif
    std::vector<gsVector<T> > partial(nThreads, gsVector<T>::Zero(nc));

#   pragma omp parallel
    {
        // Persistent per-thread scratch (see gsOptL2::GradScratch).
        GradScratch & sc = m_gradScratch.mine();
        if (!sc.ready) sc.init(dd, tdS, td_sol, nc);

        gsVector<T> & thResult = sc.thResult;
        thResult.setZero();

        gsFuncData<T> & compData       = sc.compData;
        gsFuncData<T> & geomData       = sc.geomData;
        gsFuncData<T> & sigmaBasisData = sc.sigmaBasisData;
        gsFuncData<T> & funData        = sc.funData;
        compData.flags       = NEED_VALUE | NEED_DERIV;
        geomData.flags       = NEED_VALUE | NEED_DERIV | NEED_DERIV2;
        sigmaBasisData.flags = NEED_ACTIVE | NEED_VALUE;
        funData.flags        = NEED_VALUE | NEED_DERIV;

        gsMatrix<T> & solVals = sc.solVals;
        gsMatrix<T> & dDetJs  = sc.dDetJs; // nc x nPts, d(det J_sigma)/d(alpha_i), free-control indexed
        gsMatrix<T> & Js = sc.Js, & JS = sc.JS;
        gsMatrix<T> & Cg = sc.Cg, & Cg_inv = sc.Cg_inv;
        gsMatrix<T> & adjJS = sc.adjJS, & D_d = sc.D_d, & E_d = sc.E_d;
        gsVector<T> & r = sc.r;
        gsMatrix<T> & dFdxi = sc.dFdxi;
        gsMatrix<T> & uvPoints   = sc.uvPoints;
        gsVector<T> & tmpWeights = sc.tmpWeights;

        for (auto & elem : dom->allElements())
        {
            if (sc.QuPatch != elem.patch())
            {
                sc.QuPatch = elem.patch();
                sc.QuRule = gsQuadrature::get(m_ib->basis(sc.QuPatch), quadOptions);
            }
            sc.QuRule.mapTo(elem.lowerCorner(), elem.upperCorner(), uvPoints, tmpWeights);

            m_domain->compute(uvPoints, compData);
            // d(det J_sigma)/d(alpha_i), already indexed by FREE control i -- the
            // same call gsFoldBarrier's Sampled mode uses.
            m_domain->detJacobianDeriv_into(uvPoints, dDetJs);
            m_geom->compute(compData.values[0], geomData);
            sigmaBasis.compute(uvPoints, sigmaBasisData);
            m_solution->eval_into(uvPoints, solVals);
            m_fun->compute(m_parametric ? compData.values[0] : geomData.values[0], funData);

            const gsMatrix<index_t> & actives = sigmaBasisData.actives;
            const gsMatrix<T> & basisVals     = sigmaBasisData.values[0];
            const gsMatrix<T> & deriv2_geom   = geomData.values[2];

            const index_t nPts = uvPoints.cols();
            for (index_t p = 0; p != nPts; ++p)
            {
                Js.noalias() = compData.values[1].col(p).reshaped(dd, dd).transpose();
                JS.noalias() = geomData.values[1].col(p).reshaped(dd, tdS).transpose();
                const T detJs = Js.determinant();

                // g_S: planar keeps the SIGNED det(J_S); surface uses the
                // non-negative area element sqrt(det(Cg)), Cg = J_S^T J_S.
                // det(J_sigma) stays signed in both cases -- the fold barrier
                // is what keeps it positive, not this objective.
                T gS;
                if (tdS == dd)
                    gS = JS.determinant();
                else
                {
                    Cg.noalias() = JS.transpose() * JS;
                    Cg_inv.noalias() = Cg.inverse();
                    gS = math::sqrt(math::max(Cg.determinant(), T(0)));
                }
                const T w = gS * detJs;

                // r(v) = u_h(v) - F(xi), r2 = ||r||^2 (sum over solution components)
                T r2 = T(0);
                for (short_t c = 0; c != td_sol; ++c)
                {
                    r(c) = solVals(c, p) - funData.values[0](c, p);
                    r2 += r(c) * r(c);
                }

                // dF/dxi (td_sol x dd): parametric==true -> fun's own gradient at
                // xi directly; parametric==false -> chain rule through S, J_f*J_S.
                if (m_parametric)
                {
                    for (short_t c = 0; c != td_sol; ++c)
                        for (short_t d = 0; d != dd; ++d)
                            dFdxi(c, d) = funData.values[1](c * dd + d, p);
                }
                else
                {
                    for (short_t c = 0; c != td_sol; ++c)
                        for (short_t d = 0; d != dd; ++d)
                        {
                            T s = 0;
                            for (short_t a = 0; a != tdS; ++a)
                                s += funData.values[1](c * tdS + a, p) * JS(a, d);
                            dFdxi(c, d) = s;
                        }
                }

                // d_d(g_S): planar -- tr(adj(J_S) * D_d), D_d(a,j) = d^2 S_a/(dxi_j
                // dxi_d) (the same Hessian-slice indexing gsOptMesh::gradObj_into
                // uses); surface -- 0.5*g_S*tr(Cg^{-1} E_d), E_d = D_d^T J_S +
                // J_S^T D_d (gsOptMesh's own Eq.18 formula, imitated verbatim --
                // see gsAdaptiveParametrization.hpp:1061-1072, :1143-1170, :1299).
                if (tdS == dd)
                {
                    adjJS(0,0) =  JS(1,1); adjJS(0,1) = -JS(0,1);
                    adjJS(1,0) = -JS(1,0); adjJS(1,1) =  JS(0,0);
                }

                T rDotDF[2] = {T(0), T(0)};
                T d_d_gS[2] = {T(0), T(0)};
                for (short_t d = 0; d != dd; ++d)
                {
                    for (short_t a = 0; a != tdS; ++a)
                        for (short_t j = 0; j != dd; ++j)
                        {
                            index_t lo = math::min((index_t)d, (index_t)j);
                            index_t hi = math::max((index_t)d, (index_t)j);
                            index_t hess_idx = (lo == hi) ? lo : dd + lo * (2 * dd - lo - 3) / 2 + hi - 1;
                            D_d(a, j) = deriv2_geom(a * S2 + hess_idx, p);
                        }
                    if (tdS == dd)
                        d_d_gS[d] = (adjJS * D_d).trace();
                    else
                    {
                        E_d.noalias() = D_d.transpose() * JS + JS.transpose() * D_d;
                        d_d_gS[d] = T(0.5) * gS * (Cg_inv * E_d).trace();
                    }

                    T s = 0;
                    for (short_t c = 0; c != td_sol; ++c)
                        s += r(c) * dFdxi(c, d);
                    rDotDF[d] = s;
                }

                const T coefA0 = -T(2) * rDotDF[0] * w;
                const T coefA1 = -T(2) * rDotDF[1] * w;
                const T coefB0 = r2 * d_d_gS[0] * detJs;
                const T coefB1 = r2 * d_d_gS[1] * detJs;

                const index_t nActive = actives.rows();
                for (index_t loc = 0; loc != nActive; ++loc)
                {
                    const index_t k = actives(loc, p);
                    const T Nk = basisVals(loc, p);

                    if (mapper.is_free(k, 0, 0))
                        thResult(mapper.index(k, 0, 0)) += tmpWeights[p] * Nk * (coefA0 + coefB0);
                    if (mapper.is_free(k, 0, 1))
                        thResult(mapper.index(k, 0, 1)) += tmpWeights[p] * Nk * (coefA1 + coefB1);
                }

                // g_S * d(det J_sigma)/d(alpha_i) -- already indexed by free
                // control i by detJacobianDeriv_into, so no N_k scatter here.
                for (index_t i = 0; i != nc; ++i)
                    thResult(i) += tmpWeights[p] * r2 * gS * dDetJs(i, p);
            }
        }

#       ifdef _OPENMP
        partial[omp_get_thread_num()] = thResult;
#       else
        partial[0] = thResult;
#       endif
    } // omp parallel

    for (int i = 0; i != nThreads; ++i)
        result += partial[i];

    m_barrier.addGrad(result);

    // Mirrors evalObj()'s non-finite policy (see gsOptMesh::evalObj's comment).
    for (index_t i = 0; i != result.rows(); ++i)
        if (!math::isfinite(result[i]))
            result[i] = T(0);
}

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
        GISMO_ENSURE(comp_tbasis_ptr,"The composition must be a tensor B-spline or tensor NURBS basis (2D)");
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
        // gsTensorNurbsBasis derives from gsTensorBSplineBasis (checked above first),
        // gsTHBSplineBasis does not derive from either, so this arm is disjoint.
        else if (const gsTHBSplineBasis<2,T> * hbasis = dynamic_cast<const gsTHBSplineBasis<2,T> *>(&integrationBasis))
        {
            m_integrationBasis = makeIntegrationBasis<2>(*hbasis,comp_tbasis);
        }
        else
            GISMO_ERROR("The integration basis must be either a tensor B-spline or a tensor NURBS basis (2D)");
    }
    else if (dd == 3)
    {
        const gsTensorBSplineBasis<3,T> * comp_tbasis_ptr = dynamic_cast<const gsTensorBSplineBasis<3,T> *>(&composition.domain().basis());
        GISMO_ENSURE(comp_tbasis_ptr,"The composition must be a tensor B-spline basis (3D)");
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

// Pass-through constructor (integrationBasisIsFinal_t tag): unlike the
// dispatching constructor above, this one does NOT insert sigma's knots or
// raise the degree -- \a integrationBasis IS the integration mesh already,
// verbatim. No dimension dispatch for BUILDING the basis either: the caller
// committed to a concrete basis, so there is nothing left to infer. What IS
// still checked is that the caller's basis is admissible: same domainDim as
// the composition, and at least as fine as sigma's own knot mesh -- silently
// accepting a coarser basis would under-integrate without any indication.
template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsSquareDomain<T> & composition,
                                                                const gsGeometry<T>     & geometry,
                                                                const gsFunction<T>     * function,
                                                                const gsBasis<T>        & integrationBasis,
                                                                      gsOptimizer<T>    & optimizer,
                                                                const bool                parametric,
                                                                      integrationBasisIsFinal_t)
:
m_comp(composition),
m_geom(geometry),
m_fun(function),
m_optimizer(optimizer)
{
    GISMO_ENSURE(integrationBasis.domainDim() == composition.domainDim(),
        "The integration basis has domainDim()="<<integrationBasis.domainDim()
        <<", but the composition has domainDim()="<<composition.domainDim()<<".");

    // Admissibility check: the supplied basis must have the same parameter
    // range as sigma, a degree at least sigma's, and all of sigma's interior
    // knots present. Checked for tensor B-spline/NURBS bases directly and for
    // gsHTensorBasis-derived hierarchical bases via their finest tensor level;
    // gsRationalTHBSplineBasis and any other concrete basis type are left
    // unchecked beyond domainDim. Coverage rationale: \ref adaptparam_supermesh.
    const index_t dd = composition.domainDim();
    if (dd == 2)
    {
        const gsTensorBSplineBasis<2,T> * comp_tbasis_ptr = dynamic_cast<const gsTensorBSplineBasis<2,T> *>(&composition.domain().basis());
        GISMO_ENSURE(comp_tbasis_ptr,"The composition must be a tensor B-spline or tensor NURBS basis (2D)");
        const gsTensorBSplineBasis<2,T> & comp_tbasis = *comp_tbasis_ptr;

        if (const gsTensorBSplineBasis<2,T> * tbasis = dynamic_cast<const gsTensorBSplineBasis<2,T> *>(&integrationBasis))
            GISMO_ENSURE(tensorBasisAdmissible<2>(comp_tbasis,*tbasis),
                "The supplied integration basis is not admissible for sigma (see the preceding gsWarn for the offending property).");
        else if (const gsTensorNurbsBasis<2,T> * nbasis = dynamic_cast<const gsTensorNurbsBasis<2,T> *>(&integrationBasis))
            GISMO_ENSURE(tensorBasisAdmissible<2>(comp_tbasis,nbasis->source()),
                "The supplied integration basis is not admissible for sigma (see the preceding gsWarn for the offending property).");
        else if (const gsHTensorBasis<2,T> * hbasis = dynamic_cast<const gsHTensorBasis<2,T> *>(&integrationBasis))
            // gsHTensorBasis is the shared base of every hierarchical basis
            // (gsTHBSplineBasis<2,T,true> and its alias gsHBSplineBasis<2,T>);
            // checked against the finest tensor level, deliberately not via
            // sigmaLevelInHierarchy(). Rationale: \ref adaptparam_supermesh.
            GISMO_ENSURE(tensorBasisAdmissible<2>(comp_tbasis,hbasis->tensorLevel(hbasis->maxLevel())),
                "The supplied integration basis is not admissible for sigma (see the preceding gsWarn for the offending property).");
        // else: basis type not covered by a cheap admissibility check -- left
        // unchecked, see the comment above the domainDim test.
    }
    else if (dd == 3)
    {
        const gsTensorBSplineBasis<3,T> * comp_tbasis_ptr = dynamic_cast<const gsTensorBSplineBasis<3,T> *>(&composition.domain().basis());
        GISMO_ENSURE(comp_tbasis_ptr,"The composition must be a tensor B-spline basis (3D)");
        const gsTensorBSplineBasis<3,T> & comp_tbasis = *comp_tbasis_ptr;

        if (const gsTensorBSplineBasis<3,T> * tbasis = dynamic_cast<const gsTensorBSplineBasis<3,T> *>(&integrationBasis))
            GISMO_ENSURE(tensorBasisAdmissible<3>(comp_tbasis,*tbasis),
                "The supplied integration basis is not admissible for sigma (see the preceding gsWarn for the offending property).");
        else if (const gsTensorNurbsBasis<3,T> * nbasis = dynamic_cast<const gsTensorNurbsBasis<3,T> *>(&integrationBasis))
            GISMO_ENSURE(tensorBasisAdmissible<3>(comp_tbasis,nbasis->source()),
                "The supplied integration basis is not admissible for sigma (see the preceding gsWarn for the offending property).");
        else if (const gsHTensorBasis<3,T> * hbasis = dynamic_cast<const gsHTensorBasis<3,T> *>(&integrationBasis))
            // Dimension-agnostic check against the finest tensor level, as
            // the 2D case above; this constructor only CHECKS an
            // already-built basis, it does not build a 3D hierarchical union.
            GISMO_ENSURE(tensorBasisAdmissible<3>(comp_tbasis,hbasis->tensorLevel(hbasis->maxLevel())),
                "The supplied integration basis is not admissible for sigma (see the preceding gsWarn for the offending property).");
        // else: basis type not covered by a cheap admissibility check (e.g.
        // gsRationalTHBSplineBasis) -- left unchecked.
    }
    else
        GISMO_ERROR("Only 2D and 3D composition domains are supported");

    m_integrationBasis = integrationBasis.clone();
    m_optProblem = gsOptMesh<T,MODE>(composition,m_geom,m_fun,m_integrationBasis.get(),parametric);
    this->defaultOptions();
}

template <class T, enum MonitorMode MODE>
template <short_t _d>
gsTensorBSplineBasis<_d,T> gsAdaptiveParametrization<T,MODE>::makeIntegrationBasis(const gsTensorBSplineBasis<_d,T> & basis1,
                                                                                  const gsTensorBSplineBasis<_d,T> & basis2)
{
    gsTensorBSplineBasis<_d,T> ibasis(basis1);
    // Integration basis: parent basis with knots of composition basis
    // inserted, then degree-raised per direction to
    // targetDegree = p_analysis * p_sigma. This is the uncoupled bound; the
    // true per-direction degree of a coupled sigma can reach d*p_G*p_sigma,
    // so the quadrature rule sized from targetDegree over-integrates rather
    // than integrates exactly. Full derivation, the N(x,y)=x*y counterexample
    // and the resulting quadrature-order caveats: \ref adaptparam_supermesh.
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
template <short_t d>
typename gsBasis<T>::uPtr
gsAdaptiveParametrization<T,MODE>::makeIntegrationBasis(const gsTHBSplineBasis<d,T>     & basis1,
                                                         const gsTensorBSplineBasis<d,T> & basis2)
{
    // Precondition of the gsHTensorBasis::merge() below: sigma's mesh must be a
    // dyadic refinement level of basis1's hierarchy, sharing its level 0. The
    // union mesh is built by transplanting both partitions onto that shared
    // level-0 basis; a non-nested pair has no such common root.
    index_t sigmaLevel;
    std::string reason;
    GISMO_ENSURE(sigmaLevelInHierarchy<d>(basis1,basis2,sigmaLevel,reason),
        "makeIntegrationBasis: sigma's mesh is not a dyadic level of the analysis "
        "hierarchy (" << reason << "). Use the tensor overload "
        "makeIntegrationBasis<d>(basis1.tensorLevel(basis1.maxLevel()),basis2) for a "
        "non-nested pair.");

    // B0: the shared level-0 basis of both hierarchies, degree-raised per
    // direction to the product degree p_analysis * p_sigma -- same target
    // degree and non-exactness argument as the tensor overload above (\ref
    // adaptparam_supermesh). degreeIncrease raises the degree but leaves the
    // breakpoints -- and hence every dyadic level built on B0 -- unchanged,
    // so B0 can be shared by the analysis-derived and sigma-derived hierarchies.
    gsTensorBSplineBasis<d,T> B0(basis1.tensorLevel(0));
    for (short_t dir = 0; dir != d; ++dir)
    {
        const index_t target = basis1.degree(dir) * basis2.degree(dir);
        B0.degreeIncrease(target - B0.degree(dir), dir);
    }

    // H_analysis: the analysis partition transplanted onto B0. Boxes are read
    // off in each box's own level-index space (getBoxesInLevelIndex).
    gsMatrix<index_t> b1, b2;
    gsVector<index_t> lvl;
    basis1.tree().getBoxesInLevelIndex(b1,b2,lvl);

    std::vector<index_t> boxes;
    boxes.reserve(b1.rows()*(2*d+1));
    for (index_t i = 0; i != b1.rows(); ++i)
    {
        boxes.push_back(lvl[i]);
        for (short_t dir = 0; dir != d; ++dir)
            boxes.push_back(b1(i,dir));
        for (short_t dir = 0; dir != d; ++dir)
            boxes.push_back(b2(i,dir));
    }
    gsTHBSplineBasis<d,T> Hanalysis(B0,boxes);

    // H_sigma: sigma's own mesh as one level (sigmaLevel) of the same B0
    // hierarchy. tensorLevel(sigmaLevel) creates the level on demand, so
    // sigmaLevel may exceed basis1's max level.  uniformRefine is NOT usable
    // here: it drops level 0 from the hierarchy, destroying the shared root
    // basis that merge() below requires.
    gsTHBSplineBasis<d,T> Hsigma(B0);
    std::vector<index_t> sbox;
    sbox.push_back(sigmaLevel);
    for (short_t dir = 0; dir != d; ++dir)
        sbox.push_back(0);
    for (short_t dir = 0; dir != d; ++dir)
        sbox.push_back(static_cast<index_t>(Hsigma.tensorLevel(sigmaLevel).knots(dir).numElements()));
    Hsigma.refineElements(sbox);

    // Element-wise union of the two partitions == max(analysisLevel,L_sigma).
    Hanalysis.merge(Hsigma);

    typename gsBasis<T>::uPtr result = memory::make_unique(new gsTHBSplineBasis<d,T>(Hanalysis));
    gsInfo << "makeIntegrationBasis: hierarchical, L_sigma=" << sigmaLevel
           << ", maxLevel=" << Hanalysis.maxLevel()
           << ", elements=" << result->numElements() << ", degree=";
    for (short_t dir = 0; dir != d; ++dir)
        gsInfo << Hanalysis.degree(dir) << (dir+1!=d ? "x" : "");
    gsInfo << "\n";
    return result;
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
    m_options.addReal("quA","Quadrature nodes per direction: deg*quA + quB",1.0);
    m_options.addInt ("quB","Quadrature nodes per direction: deg*quA + quB",1);
}

template <class T, enum MonitorMode MODE>
void gsAdaptiveParametrization<T,MODE>::solve()
{
    // Set the optimization problem
    m_optProblem.options().setReal("Smoothing",m_options.getReal("Smoothing"));
    m_optProblem.options().setReal("Penalty",m_options.getReal("Penalty"));
    m_optProblem.options().setReal("quA",m_options.getReal("quA"));
    m_optProblem.options().setInt ("quB",m_options.getInt ("quB"));
    m_optimizer.setProblem(&m_optProblem);

    // One-off sanity sweep (per solve, not per evaluation): the "Penalty"
    // option is the Garanzha fold-regulariser radius eps in
    // chi(x)=0.5*(x+sqrt(eps^2+x^2)); it only acts as a barrier when
    // eps << |det J_sigma| = O(1). Warn if it is comparable to or larger
    // than the current min det J_sigma (chi then degenerates towards a
    // constant and folds go unpenalised), and if it underflows chi(0).
    {
        const T penalty = m_options.getReal("Penalty");
        gsVector<T> controls0 = m_optProblem.composition().getControls();
        gsAsConstVector<T> u0(controls0.data(), controls0.rows());
        const T minDet0 = m_optProblem.computeMinJacobian(u0);
        if (minDet0 > T(0) && penalty >= T(0.1) * minDet0)
            gsWarn << "gsAdaptiveParametrization: Penalty (Garanzha radius eps="
                   << penalty << ") is not small compared to min det J_sigma = "
                   << minDet0 << "; the fold regulariser is then nearly constant "
                      "and folds are not penalised. Use eps in [1e-3, 1e-2] for a "
                      "unit-square sigma.\n";
        if (penalty < T(1e-12))
            gsWarn << "gsAdaptiveParametrization: Penalty (Garanzha radius eps="
                   << penalty << ") is smaller than 1e-12; chi(0) underflows and "
                      "the fold regulariser is effectively disabled.\n";
    }

    // Solve the optimization problem
    gsVector<T> controls = m_optProblem.composition().getControls();
    m_optimizer.solve(controls);
    controls = m_optimizer.currentDesign();
    m_optProblem.composition().setControls(controls);
    gsInfo<<"Finished with objective value: "<<m_optimizer.objective()<<std::endl;

    // Certified (sampling-free) lower bound on det J_sigma of the accepted
    // design; conservative, so this warning is not proof of a fold, but a
    // clean pass is a genuine unfolded-ness certificate.
    const T certJ = m_comp.minDetJCoefficient();
    if (certJ <= T(0))
        gsWarn << "gsAdaptiveParametrization: certified lower bound on det J_sigma "
                  "is <= 0 after relocation (bound = " << certJ << "); the map may "
                  "be folded (the bound is conservative).\n";
}

template<class T, enum MonitorMode MODE>
T gsAdaptiveParametrization<T,MODE>::computeMinJacobian() const
{
    gsVector<T> controls = m_comp.getControls();
    gsAsConstVector<T> u(controls.data(), controls.rows());
    return m_optProblem.computeMinJacobian(u);
}


} // namespace gismo

