/** @file monitor_gradcheck.cpp

    @brief Diagnostic driver: finite-difference check of the analytic gradient
    of gsOptMesh (the objective behind gsAdaptiveParametrization), for both
    MonitorMode::ValueBased and MonitorMode::GradientBased.

    This is a pure DIAGNOSTIC tool: it never calls solve(), it never modifies
    the library.  It answers three questions:

      1. Is f(x0) finite at the initial (identity) design vector, or already nan?
      2. Is the analytic gradient g(x0) finite, and does it match a central
         finite-difference gradient over a sweep of step sizes h?
      3. Are the monitor-function data (value, gradient, Hessian) that
         GradientBased needs actually finite on the domain?

    The relative-error measure matches gsOptFit::checkGradient in
    src/gsAssembler/gsAdaptiveParametrization.hpp:1608, so the two FD checks
    are directly comparable:
        rel = |g_fd - g_an| / max(|g_fd|,|g_an|),  skipped if the scale < 1e-10.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

//! [Include namespace]
#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsOptimizer/gsGradientDescent.h>

#include <iomanip>
#include <limits>

using namespace gismo;

/// Portable finite check (real_t may be float/double/long double)
static inline bool isFin(const real_t v)
{
    return (v == v) && (v <  std::numeric_limits<real_t>::infinity())
                    && (v > -std::numeric_limits<real_t>::infinity());
}

/// Human-readable classification of a scalar
static inline std::string classify(const real_t v)
{
    if (v != v)                                          return "nan";
    if (v ==  std::numeric_limits<real_t>::infinity())    return "+inf";
    if (v == -std::numeric_limits<real_t>::infinity())    return "-inf";
    return "finite";
}

/**
 * @brief Probe that exposes the protected gsOptMesh of gsAdaptiveParametrization.
 *
 * Construction is *inherited*, so the objective built here is bit-for-bit the
 * one the example builds (same union integration basis via
 * gsAdaptiveParametrization::makeIntegrationBasis, same quadrature).  Nothing
 * is re-derived in this driver.
 */
template<class T, enum MonitorMode MODE>
class gsAdaptiveParametrizationProbe : public gsAdaptiveParametrization<T,MODE>
{
    typedef gsAdaptiveParametrization<T,MODE> Base;
public:
    using Base::Base;   // inherit all constructors

    gsOptMesh<T,MODE> & problem() { return this->m_optProblem; }

    /// The union integration basis that gsOptMesh integrates on: needed to
    /// replicate its quadrature exactly in the term-by-term integrand scan.
    const gsBasis<T> & integrationBasis() const { return *this->m_integrationBasis; }
};

/// min/max/non-finite accumulator for one intermediate quantity
struct termStat
{
    termStat() : mn(std::numeric_limits<real_t>::max()),
                 mx(-std::numeric_limits<real_t>::max()), nBad(0) {}
    real_t mn, mx; index_t nBad;
    void add(const real_t v)
    {
        if (!isFin(v)) { ++nBad; return; }
        if (v < mn) mn = v;
        if (v > mx) mx = v;
    }
    void print(const std::string & name) const
    {
        gsInfo << "      " << std::left << std::setw(14) << name << std::right
               << " min = " << std::setw(14) << mn
               << "   max = " << std::setw(14) << mx
               << "   non-finite = " << nBad
               << (nBad ? "   <<<<<< " : "") << "\n";
    }
};

/**
 * @brief Term-by-term scan of the evalObj integrand at a prescribed design.
 *
 * Replicates gsOptMesh<T,MODE>::evalObj (gsAdaptiveParametrization.hpp:211-417)
 * exactly -- same union integration basis, same quadrature (quA=1,quB=1), same
 * formulas -- but reports the range and the non-finite count of every
 * intermediate quantity instead of only the integral.  This is what identifies
 * WHICH term manufactures a nan.
 */
template<enum MonitorMode MODE>
void scanIntegrand(gsAdaptiveParametrizationProbe<real_t,MODE> & probe,
                         gsSquareDomain<real_t>  & domain,
                   const gsGeometry<real_t>      & geom,
                   const gsFunction<real_t>      * fun,
                   const bool                      parametric,
                   const gsVector<real_t>        & x)
{
    const index_t dd = domain.domainDim();
    const index_t td = geom.targetDim();
    const real_t theta   = probe.problem().options().getReal("Smoothing");
    const real_t penalty = probe.problem().options().getReal("Penalty");
    const bool hasMonitor = (fun != nullptr);

    domain.setControls(x);

    const gsBasis<real_t> & ib = probe.integrationBasis();
    gsMultiBasis<real_t> mb(ib,true);

    gsOptionList quadOptions;
    quadOptions.addReal("quA","",1.0);
    quadOptions.addInt ("quB","",1);

    termStat sSigma, sGeomJ, sDet, sChi, sTrCinv, sEta, sM2, sOmega, sInt, sFunV, sFunD;
    real_t   integral = 0;
    index_t  nPtsTot = 0, nOutside = 0, nBadInt = 0, nBadOutside = 0;
    bool     firstBadPrinted = false;
    const gsMatrix<real_t> sup = geom.support();

    gsFuncData<real_t> compData, geomData, funData;
    compData.flags = NEED_VALUE | NEED_DERIV;
    geomData.flags = NEED_VALUE | NEED_DERIV;
    if (hasMonitor)
        funData.flags = (MODE==ValueBased) ? NEED_VALUE : (NEED_VALUE|NEED_DERIV);

    gsMatrix<real_t> Js(dd,dd), Jg(td,dd), Jc(td,dd), C(dd,dd), Cinv(dd,dd);
    gsMatrix<real_t> Cg(dd,dd), Cg_inv(dd,dd), gxi(dd,1);
    gsQuadRule<real_t> QuRule;
    index_t QuPatch = -1;
    gsMatrix<real_t> uvPoints;
    gsVector<real_t> tmpWeights;

    auto dom = mb.domain();
    for (auto & elem : dom->allElements())
    {
        if (QuPatch != elem.patch())
        {
            QuPatch = elem.patch();
            QuRule = gsQuadrature::get(mb.basis(QuPatch), quadOptions);
        }
        QuRule.mapTo(elem.lowerCorner(), elem.upperCorner(), uvPoints, tmpWeights);

        domain.compute(uvPoints, compData);
        geom.compute(compData.values[0], geomData);
        if (hasMonitor)
            fun->compute(parametric ? compData.values[0] : geomData.values[0], funData);

        for (index_t p = 0; p != uvPoints.cols(); ++p)
        {
            ++nPtsTot;
            Js.noalias() = compData.values[1].col(p).reshaped(dd,dd).transpose();
            Jg.noalias() = geomData.values[1].col(p).reshaped(dd,td).transpose();
            Jc.noalias() = Jg*Js;
            C.noalias()  = Jc.transpose()*Jc;
            C.diagonal().array() += penalty;
            Cinv.noalias() = C.inverse();
            const real_t trCinv = Cinv.trace();
            sSigma.add(Js.determinant());
            sGeomJ.add(Jg.norm());
            sTrCinv.add(trCinv);

            real_t detq, chi, absdet;
            if (td == dd) { detq = Jc.determinant(); }
            else          { detq = Js.determinant(); }
            chi    = real_t(0.5)*(detq + math::sqrt(penalty*penalty + detq*detq));
            absdet = math::abs(detq);
            sDet.add(detq);
            sChi.add(chi);

            real_t m2 = 0, omega = 0, eta = 0, integrand = 0, gS = 1;
            if (hasMonitor)
            {
                sFunV.add(funData.values[0](0,p));
                if (MODE == ValueBased)
                {
                    eta = funData.values[0](0,p);
                    m2  = real_t(1)/(real_t(1)+theta*eta);
                }
                else
                {
                    Cg.noalias() = Jg.transpose()*Jg;
                    Cg_inv.noalias() = Cg.inverse();
                    if (parametric) gxi.noalias() = funData.values[1].col(p);
                    else            gxi.noalias() = Jg.transpose()*funData.values[1].col(p).reshaped(td,1);
                    sFunD.add(gxi.norm());
                    eta = (gxi.transpose()*Cg_inv*gxi)(0,0);
                    m2  = real_t(1)/(real_t(1)+theta*eta);
                }
                if (td > dd) gS = math::sqrt(Cg.determinant());
                omega = math::sqrt(m2);
                integrand = omega*trCinv*absdet*absdet*absdet/(chi*chi)*(td>dd ? gS : real_t(1));
                sEta.add(eta); sM2.add(m2); sOmega.add(omega);
            }
            else
                integrand = trCinv*absdet*absdet/(chi*chi);

            sInt.add(integrand);
            integral += tmpWeights[p]*integrand;

            // Is sigma(uhat) inside the parameter domain of the geometry?
            // Outside it, G is evaluated by extrapolation/clamping and J_g can
            // degenerate -- which is what makes Cg = Jg^T Jg singular.
            bool inSup = true;
            for (index_t d = 0; d != dd; ++d)
                if (compData.values[0](d,p) < sup(d,0) - 1e-12 ||
                    compData.values[0](d,p) > sup(d,1) + 1e-12) inSup = false;
            if (!inSup) ++nOutside;
            if (!isFin(integrand)) { ++nBadInt; if (!inSup) ++nBadOutside; }

            if (!isFin(integrand) && !firstBadPrinted)
            {
                firstBadPrinted = true;
                gsMatrix<real_t> CgLoc = Jg.transpose()*Jg;
                gsInfo << "      FIRST non-finite integrand at uhat = ("
                       << uvPoints(0,p) << "," << uvPoints(1,p) << "):\n"
                       << "        sigma(uhat) = " << compData.values[0].col(p).transpose()
                       << "   inside geometry support [" << sup(0,0) << "," << sup(0,1)
                       << "]x[" << sup(1,0) << "," << sup(1,1) << "] ? "
                       << (inSup ? "YES" : "NO") << "\n"
                       << "        G(sigma)    = " << geomData.values[0].col(p).transpose() << "\n"
                       << "        J_g         = [" << Jg.row(0) << " ; " << Jg.row(1) << "]"
                       << "   |J_g|_F = " << Jg.norm() << "\n"
                       << "        Cg = Jg^T Jg, det(Cg) = " << CgLoc.determinant()
                       << ", cond-proxy |Cg^-1|_F = " << CgLoc.inverse().norm() << "\n"
                       << "        det(Js) = " << Js.determinant()
                       << ", detq = " << detq << ", chi = " << chi << "\n"
                       << "        tr(Cinv) = " << trCinv
                       << ", eta = " << eta << ", m2 = " << m2
                       << ", omega = " << omega << ", gS = " << gS << "\n";
            }
        }
    }

    gsInfo << "      integral = " << integral << " (" << classify(integral) << "), "
           << nPtsTot << " quadrature points\n";
    gsInfo << "      quadrature points with sigma(uhat) OUTSIDE the geometry "
              "parameter domain: " << nOutside << " / " << nPtsTot << "\n";
    gsInfo << "      non-finite integrands: " << nBadInt
           << " (of which " << nBadOutside << " at points with sigma outside)\n";
    sSigma.print("det(J_sigma)");
    sDet  .print(td==dd ? "det_c" : "det_s");
    sChi  .print("chi");
    sTrCinv.print("tr(C^-1)");
    sGeomJ.print("|J_g|_F");
    if (hasMonitor)
    {
        sFunV.print("f");
        if (MODE!=ValueBased) sFunD.print("|grad_xi f|");
        sEta  .print(MODE==ValueBased ? "eta = f" : "eta2");
        sM2   .print("m2");
        sOmega.print("omega=sqrt");
    }
    sInt  .print("integrand");
}

/**
 * @brief Write-coverage test of gsFunctionExpr::deriv2_into.
 *
 * Pre-fills the result matrix (already at the exact final size, so that the
 * internal resize() cannot reallocate) with a sentinel and checks which rows
 * deriv2_into actually writes.  A surviving sentinel proves that the row is
 * never written, i.e. that the caller reads uninitialised/stale memory.
 */
void deriv2WriteCoverage(const gsFunctionExpr<real_t> & fun)
{
    const index_t d = fun.domainDim();
    const index_t n = fun.targetDim();
    const index_t stride = d + d*(d-1)/2;

    gsMatrix<real_t> pts(d,3);
    pts.setConstant(real_t(0.3));
    pts(0,1) = 0.4; pts(0,2) = 0.5;   // three distinct points

    const real_t sentinel = -123456.75;
    gsMatrix<real_t> res(stride*n, pts.cols());
    res.setConstant(sentinel);
    fun.deriv2_into(pts, res);

    gsInfo << "\n[0] gsFunctionExpr::deriv2_into write-coverage test"
           << "   (d = " << d << ", targetDim = " << n
           << ", stride = " << stride << ")\n";
    gsInfo << "    expected row layout: [";
    for (index_t i = 0; i != d; ++i) gsInfo << (i?" ":"") << "d" << i << i;
    for (index_t i = 0; i != d; ++i)
        for (index_t j = i+1; j != d; ++j) gsInfo << " d" << i << j;
    gsInfo << "]\n";
    index_t nUnwritten = 0;
    for (index_t r = 0; r != res.rows(); ++r)
    {
        bool untouched = true;
        for (index_t p = 0; p != res.cols(); ++p)
            if (res(r,p) != sentinel) untouched = false;
        if (untouched) ++nUnwritten;
        gsInfo << "      row " << std::setw(2) << r << " : "
               << (untouched ? "NEVER WRITTEN  <<<<<< reads stale memory"
                             : "written") << "\n";
    }
    gsInfo << "    rows never written: " << nUnwritten << " / " << res.rows()
           << (nUnwritten ? "   *** BUG ***" : "   (ok)") << "\n";
}

/**
 * @brief Scan of the monitor function data on a uniform grid of the domain.
 *
 * GradientBased needs f, grad f (evalObj) and additionally the Hessian of f
 * (gradObj_into, gsAdaptiveParametrization.hpp:579,637).  ValueBased needs only
 * f (evalObj).  A non-finite entry anywhere in these arrays is a candidate
 * cause for a nan objective/gradient, so we report them separately.
 *
 * Also reports the range of eta^2 = (grad_xi f)^T Cg^{-1} (grad_xi f) and of
 * m2 = 1/(1+theta*eta^2) (GradientBased) resp. m2 = 1/(1+theta*f) (ValueBased),
 * because omega = sqrt(m2) sits in a division at
 * gsAdaptiveParametrization.hpp:930 (domega_dalpha = dm2_dalpha/(2*omega)):
 * m2 -> 0 makes that 0/0 and m2 < 0 makes sqrt(m2) nan.
 */
void scanMonitor(const gsFunction<real_t> & fun,
                 const gsGeometry<real_t> & geom,
                 const bool                 parametric,
                 const real_t               theta,
                 const index_t              n1d)
{
    const index_t dd = geom.domainDim();
    const index_t td = geom.targetDim();

    // Uniform grid on the parameter domain of the geometry.  At the identity
    // composition sigma = id, these ARE the points where the monitor is
    // evaluated (up to the quadrature-point positions inside each element).
    gsMatrix<real_t> sup = geom.support();
    gsVector<unsigned> np(dd);
    np.setConstant(n1d);
    gsMatrix<real_t> xi = gsPointGrid<real_t>(sup.col(0), sup.col(1), np);
    const index_t nPts = xi.cols();

    gsMatrix<real_t> phys, Jflat;
    geom.eval_into (xi, phys);
    geom.deriv_into(xi, Jflat);

    const gsMatrix<real_t> & evalPts = parametric ? xi : phys;

    gsMatrix<real_t> vals, ders, ders2;
    fun.eval_into  (evalPts, vals );
    fun.deriv_into (evalPts, ders );
    fun.deriv2_into(evalPts, ders2);

    index_t nNanVal = 0, nNanDer = 0, nNanDer2 = 0;
    real_t maxAbsVal = 0, maxAbsDer = 0, maxAbsDer2 = 0;
    real_t minVal =  std::numeric_limits<real_t>::max();
    real_t maxVal = -std::numeric_limits<real_t>::max();

    for (index_t p = 0; p != nPts; ++p)
    {
        for (index_t r = 0; r != vals.rows(); ++r)
        {
            const real_t v = vals(r,p);
            if (!isFin(v)) ++nNanVal;
            else { maxAbsVal = math::max(maxAbsVal, math::abs(v));
                   minVal = math::min(minVal,v); maxVal = math::max(maxVal,v); }
        }
        for (index_t r = 0; r != ders.rows(); ++r)
        {
            const real_t v = ders(r,p);
            if (!isFin(v)) ++nNanDer;
            else maxAbsDer = math::max(maxAbsDer, math::abs(v));
        }
        for (index_t r = 0; r != ders2.rows(); ++r)
        {
            const real_t v = ders2(r,p);
            if (!isFin(v)) ++nNanDer2;
            else maxAbsDer2 = math::max(maxAbsDer2, math::abs(v));
        }
    }

    gsInfo << "  grid points                 : " << nPts
           << " (" << n1d << "^" << dd << ")\n";
    gsInfo << "  f      rows=" << vals.rows()
           << "  non-finite=" << nNanVal
           << "  range=[" << minVal << ", " << maxVal << "]\n";
    gsInfo << "  grad f rows=" << ders.rows()
           << "  non-finite=" << nNanDer  << "  max|.|=" << maxAbsDer  << "\n";
    gsInfo << "  hess f rows=" << ders2.rows()
           << "  non-finite=" << nNanDer2 << "  max|.|=" << maxAbsDer2 << "\n";

    // eta^2 / m2 ranges (GradientBased weight), mirroring
    // gsAdaptiveParametrization.hpp:328-338 (planar) and 377-396 (surface)
    real_t minEta2 =  std::numeric_limits<real_t>::max();
    real_t maxEta2 = -std::numeric_limits<real_t>::max();
    index_t nNegEta2 = 0, nNanEta2 = 0;
    gsMatrix<real_t> Jg(td,dd), Cg(dd,dd), Cginv(dd,dd), gxi(dd,1);
    for (index_t p = 0; p != nPts; ++p)
    {
        Jg.noalias() = Jflat.col(p).reshaped(dd,td).transpose();
        Cg.noalias() = Jg.transpose()*Jg;
        Cginv.noalias() = Cg.inverse();
        if (parametric) gxi.noalias() = ders.col(p);
        else            gxi.noalias() = Jg.transpose()*ders.col(p).reshaped(td,1);
        const real_t eta2 = (gxi.transpose()*Cginv*gxi)(0,0);
        if (!isFin(eta2)) { ++nNanEta2; continue; }
        if (eta2 < 0) ++nNegEta2;
        minEta2 = math::min(minEta2,eta2);
        maxEta2 = math::max(maxEta2,eta2);
    }
    gsInfo << "  eta^2 (GradientBased) range = [" << minEta2 << ", " << maxEta2
           << "]  negative=" << nNegEta2 << "  non-finite=" << nNanEta2 << "\n";
    if (isFin(maxEta2))
    {
        const real_t m2min = real_t(1)/(real_t(1)+theta*maxEta2);
        const real_t m2max = real_t(1)/(real_t(1)+theta*minEta2);
        gsInfo << "  m2 = 1/(1+theta*eta^2), theta=" << theta
               << " -> range=[" << m2min << ", " << m2max
               << "], omega=sqrt(m2) in [" << math::sqrt(math::abs(m2min))
               << ", " << math::sqrt(math::abs(m2max)) << "]\n";
    }
    gsInfo << "  m2 = 1/(1+theta*f)  (ValueBased) -> range=["
           << real_t(1)/(real_t(1)+theta*maxVal) << ", "
           << real_t(1)/(real_t(1)+theta*minVal) << "]"
           << ((real_t(1)+theta*minVal) <= 0 ?
               "   *** 1+theta*f <= 0 : sqrt(m2) IS nan ***" : "") << "\n";
}

/**
 * @brief Term-level FD test of d(eta^2)/d(alpha_{k,d}), the ONLY quantity that
 *   is exclusive to MonitorMode::GradientBased.
 *
 * gradObj_into computes (gsAdaptiveParametrization.hpp:832-848)
 *   mon_scalar_d(d) = -v^T E_d v + 2 v^T ( D_d^T grad_x f + Jg^T H_f Jg e_d )
 *   d(eta^2)/d(alpha_{k,d}) = N_k * mon_scalar_d(d)
 * with v = Cg^{-1} grad_xi f, E_d = D_d^T Jg + Jg^T D_d, D_d(a,j) = d2G_a/dxi_j dxi_d.
 *
 * NOTE for td==dd:  Jg Cg^{-1} Jg^T = I, hence eta^2 = ||grad_x f||^2 does not
 * depend on Jg at all and the two terms -v^T E_d v and 2 v^T D_d^T grad_x f
 * CANCEL IDENTICALLY.  A planar test therefore validates only the Hessian term;
 * the E_d / D_d^T grad_x f pair is exercised only on a surface (td>dd).  That is
 * why this routine reports several variants: whichever one matches the FD value
 * identifies the faulty term.
 */
void checkEta2Gradient(      gsSquareDomain<real_t> & domain,
                       const gsGeometry<real_t>     & geom,
                       const gsFunction<real_t>     & fun,
                       const bool                     parametric,
                       const gsVector<real_t>       & x0,
                       const index_t                  nSample)
{
    const index_t dd = domain.domainDim();
    const index_t td = geom.targetDim();
    const index_t n  = x0.rows();

    gsInfo << "\n[8] term-level FD test of d(eta^2)/d(alpha_{k,d})"
           << "   (GradientBased-exclusive term)\n";
    if (parametric)
    {
        gsInfo << "    (Parametric=1: the monitor is evaluated in the parametric "
                  "domain; this test covers the Parametric=0 path only)\n";
        return;
    }

    // interior sample points, deliberately away from element boundaries
    gsMatrix<real_t> U(dd, nSample*nSample);
    for (index_t i = 0; i != nSample; ++i)
        for (index_t j = 0; j != nSample; ++j)
        {
            U(0, i*nSample+j) = (real_t(i)+real_t(0.5))/real_t(nSample);
            U(1, i*nSample+j) = (real_t(j)+real_t(0.5))/real_t(nSample);
        }
    const index_t nPts = U.cols();

    // eta^2 at all sample points for a given control vector
    auto eta2All = [&](const gsVector<real_t> & c, gsVector<real_t> & out)
    {
        domain.setControls(c);
        gsFuncData<real_t> cD, gD, fD;
        cD.flags = NEED_VALUE | NEED_DERIV;
        gD.flags = NEED_VALUE | NEED_DERIV;
        fD.flags = NEED_VALUE | NEED_DERIV;
        domain.compute(U, cD);
        geom.compute(cD.values[0], gD);
        fun.compute(gD.values[0], fD);
        out.resize(nPts);
        gsMatrix<real_t> Jg(td,dd), Cg(dd,dd), w(dd,1);
        for (index_t p = 0; p != nPts; ++p)
        {
            Jg.noalias() = gD.values[1].col(p).reshaped(dd,td).transpose();
            Cg.noalias() = Jg.transpose()*Jg;
            w.noalias()  = Jg.transpose()*fD.values[1].col(p).reshaped(td,1);
            out(p) = (w.transpose()*Cg.inverse()*w)(0,0);
        }
    };

    // ---- base data at x0 ----
    domain.setControls(x0);
    const gsBasis<real_t> & sb = domain.domain().basis();
    const gsDofMapper & mapper = domain.mapper();

    gsFuncData<real_t> cD, gD, fD, sbD;
    cD.flags = NEED_VALUE | NEED_DERIV;
    gD.flags = NEED_VALUE | NEED_DERIV | NEED_DERIV2;
    fD.flags = NEED_VALUE | NEED_DERIV | NEED_DERIV2;
    sbD.flags = NEED_ACTIVE | NEED_VALUE | NEED_DERIV;
    domain.compute(U, cD);
    geom.compute(cD.values[0], gD);
    fun.compute(gD.values[0], fD);
    sb.compute(U, sbD);

    const index_t S = dd*(dd+1)/2;

    // Hessian of f by central FD of fun.deriv (INDEPENDENT of deriv2_into):
    // lets us separate a wrong formula from a wrong/mis-indexed deriv2.
    const real_t hh = 1e-5;
    std::vector<gsMatrix<real_t> > HfFD(nPts, gsMatrix<real_t>(td,td));
    {
        gsMatrix<real_t> xp = gD.values[0], xm = gD.values[0], dp, dm;
        for (index_t b = 0; b != td; ++b)
        {
            xp = gD.values[0]; xp.row(b).array() += hh;
            xm = gD.values[0]; xm.row(b).array() -= hh;
            fun.deriv_into(xp, dp);
            fun.deriv_into(xm, dm);
            for (index_t p = 0; p != nPts; ++p)
                for (index_t a = 0; a != td; ++a)
                    HfFD[p](a,b) = (dp(a,p)-dm(a,p))/(real_t(2)*hh);
        }
        for (index_t p = 0; p != nPts; ++p)
            HfFD[p] = real_t(0.5)*(HfFD[p]+HfFD[p].transpose()).eval();
    }

    // ---- Hess_f as the library assembles it, vs. the FD reference ----
    {
        gsMatrix<real_t> Hlib(td,td), maxDiff(td,td);
        maxDiff.setZero();
        for (index_t p = 0; p != nPts; ++p)
        {
            for (index_t i = 0; i != td; ++i)
                for (index_t j = i; j != td; ++j)
                {
                    const index_t hidx = (i==j) ? i : td + i*(2*td-i-3)/2+j-1;
                    Hlib(i,j) = fD.values[2](hidx,p);
                    Hlib(j,i) = Hlib(i,j);
                }
            for (index_t i = 0; i != td; ++i)
                for (index_t j = 0; j != td; ++j)
                    maxDiff(i,j) = math::max(maxDiff(i,j),
                                             math::abs(Hlib(i,j)-HfFD[p](i,j)));
            if (p == 0)
            {
                gsInfo << "    monitor Hessian at the first sample point x = "
                       << gD.values[0].col(0).transpose() << "\n";
                gsInfo << "      raw deriv2 rows (gsFunction layout, "
                       << fD.values[2].rows() << " rows) : "
                       << fD.values[2].col(0).transpose() << "\n";
                gsInfo << "      H_f as assembled by the library (hpp:725-731):\n";
                for (index_t i = 0; i != td; ++i)
                    gsInfo << "        " << Hlib.row(i) << "\n";
                gsInfo << "      H_f by central FD of fun.deriv (reference):\n";
                for (index_t i = 0; i != td; ++i)
                    gsInfo << "        " << HfFD[p].row(i) << "\n";
            }
        }
        gsInfo << "    max |H_lib - H_FD| over " << nPts << " points, per entry:\n";
        for (index_t i = 0; i != td; ++i)
            gsInfo << "        " << maxDiff.row(i) << "\n";
    }

    // ---- FD reference: d(eta^2)/d(design_i) for every design component ----
    const real_t h = 1e-6;
    gsMatrix<real_t> dEta2FD(nPts, n);
    for (index_t i = 0; i != n; ++i)
    {
        gsVector<real_t> xp = x0, xm = x0, ep, em;
        xp(i) += h; xm(i) -= h;
        eta2All(xp, ep);
        eta2All(xm, em);
        for (index_t p = 0; p != nPts; ++p)
            dEta2FD(p,i) = (ep(p)-em(p))/(real_t(2)*h);
    }
    domain.setControls(x0);

    // ---- analytic variants ----
    const index_t nVar = 6;
    const char * varName[nVar] = {"full (as coded)",
                                  "hessian term only",
                                  "no  -v'E_d v",
                                  "no  2v'D_d'gradf",
                                  "sign flip -> +v'E_d v",
                                  "full, H_f from FD"};
    gsVector<real_t> maxRel(nVar), sumAbs(nVar);
    maxRel.setZero(); sumAbs.setZero();
    gsVector<index_t> nCmp(nVar); nCmp.setZero();

    gsMatrix<real_t> Jg(td,dd), Cg(dd,dd), Cginv(dd,dd), D_d(td,dd);
    gsMatrix<real_t> Hess_f(td,td), gx(td,1), w(dd,1);
    gsVector<real_t> v(dd), gradNk(dd);

    for (index_t p = 0; p != nPts; ++p)
    {
        Jg.noalias() = gD.values[1].col(p).reshaped(dd,td).transpose();
        Cg.noalias() = Jg.transpose()*Jg;
        Cginv.noalias() = Cg.inverse();
        gx.noalias() = fD.values[1].col(p).reshaped(td,1);
        w.noalias() = Jg.transpose()*gx;
        v.noalias() = Cginv*w;

        // Hess_f exactly as the library fills it (hpp:725-731)
        for (index_t i = 0; i != td; ++i)
            for (index_t j = i; j != td; ++j)
            {
                const index_t hidx = (i==j) ? i : td + i*(2*td-i-3)/2+j-1;
                Hess_f(i,j) = fD.values[2](hidx,p);
                Hess_f(j,i) = Hess_f(i,j);
            }

        for (index_t d = 0; d != dd; ++d)
        {
            // D_d exactly as the library builds it (hpp:797-804)
            for (index_t a = 0; a != td; ++a)
                for (index_t j = 0; j != dd; ++j)
                {
                    const index_t lo = math::min(d,j), hi = math::max(d,j);
                    const index_t hess_idx = (lo==hi) ? lo : dd + lo*(2*dd-lo-3)/2 + hi - 1;
                    D_d(a,j) = gD.values[2](a*S + hess_idx, p);
                }
            const gsMatrix<real_t> E_d = D_d.transpose()*Jg + Jg.transpose()*D_d;
            const real_t vEv  = (v.transpose()*E_d*v)(0,0);
            const gsMatrix<real_t> Dg_d  = D_d.transpose()*gx;
            const gsMatrix<real_t> JtHJd = Jg.transpose()*(Hess_f*Jg.col(d));
            const gsMatrix<real_t> JtHJdF= Jg.transpose()*(HfFD[p]*Jg.col(d));

            const real_t term[nVar] = {
                -vEv + real_t(2)*(v.transpose()*(Dg_d+JtHJd))(0,0),
                       real_t(2)*(v.transpose()*JtHJd)(0,0),
                       real_t(2)*(v.transpose()*(Dg_d+JtHJd))(0,0),
                -vEv + real_t(2)*(v.transpose()*JtHJd)(0,0),
                 vEv + real_t(2)*(v.transpose()*(Dg_d+JtHJd))(0,0),
                -vEv + real_t(2)*(v.transpose()*(Dg_d+JtHJdF))(0,0) };

            for (index_t loc = 0; loc != sbD.actives.rows(); ++loc)
            {
                const index_t k  = sbD.actives(loc,p);
                if (!mapper.is_free(k,0,d)) continue;
                const index_t ii = mapper.index(k,0,d);
                const real_t Nk = sbD.values[0](loc,p);
                const real_t fd = dEta2FD(p,ii);
                for (index_t vv = 0; vv != nVar; ++vv)
                {
                    const real_t an = Nk*term[vv];
                    const real_t scale = math::max(math::abs(an),math::abs(fd));
                    if (scale <= 1e-8) continue;
                    maxRel(vv) = math::max(maxRel(vv), math::abs(an-fd)/scale);
                    sumAbs(vv) += math::abs(an-fd);
                    ++nCmp(vv);
                }
            }
        }
    }

    gsInfo << "    " << nPts << " sample points, comparing N_k*mon_scalar_d(d) "
           << "against central FD (h=" << h << ")\n";
    gsInfo << "    " << std::left << std::setw(26) << "variant" << std::right
           << std::setw(14) << "maxRel" << std::setw(16) << "sum|diff|"
           << std::setw(10) << "nCmp" << "\n";
    for (index_t vv = 0; vv != nVar; ++vv)
        gsInfo << "    " << std::left << std::setw(26) << varName[vv] << std::right
               << std::scientific << std::setprecision(4)
               << std::setw(14) << maxRel(vv) << std::setw(16) << sumAbs(vv)
               << std::defaultfloat << std::setprecision(6)
               << std::setw(10) << nCmp(vv) << "\n";

    domain.setControls(x0);
}

/**
 * @brief FD-vs-analytic gradient check of gsOptMesh<real_t,MODE>.
 *
 * Complexity: the h-sweep costs 2*nControls*nH objective evaluations, each
 * O(nElements * nQuadPts) -- linear in the sweep length, so keep nH modest.
 */
template<enum MonitorMode MODE>
void runCheck(const std::string          & modeName,
                    gsSquareDomain<real_t> & domain,
              const gsGeometry<real_t>     & geom,
              const gsFunction<real_t>     * fun,
              const gsBasis<real_t>        & ibasis,
                    gsOptimizer<real_t>    & opt,
              const bool                     parametric,
              const real_t                   penalty,
              const real_t                   smoothing,
              const gsVector<real_t>       & x0,
              const index_t                  nRepeat,
              const index_t                  nPert,
              const real_t                   pertMag,
              const index_t                  nShow)
{
    gsInfo << "\n===============================================================\n";
    gsInfo << " MODE = " << modeName << "   (monitor "
           << (fun==nullptr ? "ABSENT" : "present") << ")\n";
    gsInfo << "===============================================================\n";

    domain.setControls(x0);   // always start from the same design
    gsAdaptiveParametrizationProbe<real_t,MODE> probe(domain,geom,fun,ibasis,opt,parametric);

    // gsAdaptiveParametrization only forwards these two options to the
    // internal gsOptMesh inside solve() (gsAdaptiveParametrization.hpp:1271-1272).
    // We never call solve(), so set them on BOTH lists explicitly, otherwise the
    // gsOptMesh defaults (Penalty=1e-2, Smoothing=0.1) would silently be used.
    probe.options().setReal("Penalty",penalty);
    probe.options().setReal("Smoothing",smoothing);
    probe.problem().options().setReal("Penalty",penalty);
    probe.problem().options().setReal("Smoothing",smoothing);
    gsInfo << "Penalty = " << probe.problem().options().getReal("Penalty")
           << ", Smoothing = " << probe.problem().options().getReal("Smoothing")
           << ", Parametric = " << parametric << "\n";

    const index_t n = x0.rows();

    // ---- objective evaluation helper (restores nothing: setControls is done
    //      by evalObj itself; we always pass the full design vector) ----
    index_t nNonFiniteEvals = 0;
    auto fEval = [&](const gsVector<real_t> & x) -> real_t
    {
        gsAsConstVector<real_t> xv(x.data(), x.rows());
        const real_t f = probe.problem().evalObj(xv);
        if (!isFin(f)) ++nNonFiniteEvals;
        return f;
    };

    // -----------------------------------------------------------------
    // 1. f(x0) and the reproducibility (noise floor) of the objective
    // -----------------------------------------------------------------
    const real_t f0 = fEval(x0);
    gsInfo << "\n[1] f(x0) = " << std::setprecision(17) << f0
           << std::setprecision(6) << "   -> " << classify(f0) << "\n";

    real_t spread = 0;
    for (index_t i = 0; i != nRepeat; ++i)
    {
        const real_t fi = fEval(x0);
        if (isFin(fi) && isFin(f0)) spread = math::max(spread, math::abs(fi-f0));
    }
    gsInfo << "    reproducibility over " << nRepeat << " repeats: max|f_i-f_0| = "
           << spread;
    if (isFin(f0) && f0 != real_t(0))
        gsInfo << "  (relative " << spread/math::abs(f0) << ")";
    gsInfo << "\n    -> FD noise floor at step h is about "
           << "spread/(2h); h below that is pure noise.\n";

    // -----------------------------------------------------------------
    // 2. analytic gradient at x0
    // -----------------------------------------------------------------
    gsVector<real_t> ga(n);
    {
        gsAsConstVector<real_t> x0v(x0.data(), n);
        gsAsVector<real_t> gav(ga.data(), n);
        probe.problem().gradObj_into(x0v, gav);
    }
    index_t nNanG = 0;
    real_t maxAbsG = 0;
    for (index_t i = 0; i != n; ++i)
    {
        if (!isFin(ga(i))) ++nNanG;
        else maxAbsG = math::max(maxAbsG, math::abs(ga(i)));
    }
    gsInfo << "\n[2] analytic gradient at x0: n = " << n
           << ", non-finite components = " << nNanG
           << ", max|g_i| (finite ones) = " << maxAbsG << "\n";
    if (nNanG == 0)
        gsInfo << "    ||g|| = " << ga.norm() << "\n";
    else
        gsInfo << "    ||g|| = nan (non-finite components present)\n";

    if (!isFin(f0))
    {
        gsInfo << "\n    *** f(x0) IS NOT FINITE: the FD gradient is meaningless.\n"
               << "    *** The bug is in evalObj, not in gradObj_into.\n";
    }

    // -----------------------------------------------------------------
    // 3. h-sweep of the central-difference gradient
    // -----------------------------------------------------------------
    const std::vector<real_t> hs = {1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8};

    gsInfo << "\n[3] central-difference sweep  g_fd[i] = (f(x+h e_i)-f(x-h e_i))/(2h)\n";
    gsInfo << "    rel = |g_fd-g_an|/max(|g_fd|,|g_an|), skipped when scale<1e-10\n";
    gsInfo << "         (same measure as gsOptFit::checkGradient)\n\n";
    gsInfo << std::setw(10) << "h"
           << std::setw(14) << "maxRel"
           << std::setw(8)  << "argmax"
           << std::setw(14) << "medianRel"
           << std::setw(10) << "nBadFD"
           << std::setw(10) << "nSkip"
           << std::setw(16) << "||g_fd||" << "\n";

    real_t bestMaxRel = std::numeric_limits<real_t>::max();
    real_t bestH = 0;
    gsVector<real_t> bestFD(n); bestFD.setZero();

    for (size_t ih = 0; ih != hs.size(); ++ih)
    {
        const real_t h = hs[ih];
        gsVector<real_t> gfd(n);
        gsVector<real_t> xp = x0, xm = x0;
        index_t nBadFD = 0, nSkip = 0;
        std::vector<real_t> rels;
        rels.reserve(n);
        real_t maxRel = 0; index_t argmax = -1;

        for (index_t i = 0; i != n; ++i)
        {
            xp = x0; xp(i) += h;
            xm = x0; xm(i) -= h;
            const real_t fp = fEval(xp);
            const real_t fm = fEval(xm);
            if (!isFin(fp) || !isFin(fm)) { ++nBadFD; gfd(i) = std::numeric_limits<real_t>::quiet_NaN(); continue; }
            gfd(i) = (fp-fm)/(real_t(2)*h);

            const real_t scale = math::max(math::abs(gfd(i)), math::abs(ga(i)));
            if (!isFin(ga(i))) { ++nBadFD; continue; }
            if (scale <= 1e-10) { ++nSkip; continue; }
            const real_t rel = math::abs(gfd(i)-ga(i))/scale;
            rels.push_back(rel);
            if (rel > maxRel) { maxRel = rel; argmax = i; }
        }
        real_t med = 0;
        if (!rels.empty())
        {
            std::sort(rels.begin(),rels.end());
            med = rels[rels.size()/2];
        }
        bool gfdFinite = true;
        for (index_t i = 0; i != n; ++i) if (!isFin(gfd(i))) gfdFinite = false;

        gsInfo << std::scientific << std::setprecision(1) << std::setw(10) << h
               << std::setprecision(4) << std::setw(14) << maxRel
               << std::setw(8) << argmax
               << std::setw(14) << med
               << std::setw(10) << nBadFD
               << std::setw(10) << nSkip
               << std::setw(16) << (gfdFinite ? gfd.norm()
                                    : std::numeric_limits<real_t>::quiet_NaN())
               << std::defaultfloat << std::setprecision(6) << "\n";

        if (nBadFD == 0 && maxRel < bestMaxRel) { bestMaxRel = maxRel; bestH = h; bestFD = gfd; }
    }

    // -----------------------------------------------------------------
    // 4. per-component table at the best h
    // -----------------------------------------------------------------
    if (bestH != real_t(0))
    {
        gsInfo << "\n[4] per-component table at the best h = " << bestH
               << " (max rel err " << bestMaxRel << ")\n";

        // invert the mapper so components can be attributed to a (basis
        // function, coordinate) pair and to a location in the domain: this is
        // what distinguishes a global gradient error from a boundary/corner one.
        const gsBasis<real_t> & sb = domain.domain().basis();
        const gsDofMapper & mapper = domain.mapper();
        const index_t dd = domain.domainDim();
        gsMatrix<real_t> anchors;
        sb.anchors_into(anchors);
        std::vector<index_t> ck(n,-1), cd(n,-1);
        for (index_t k = 0; k != sb.size(); ++k)
            for (index_t d = 0; d != dd; ++d)
                if (mapper.is_free(k,0,d))
                {
                    const index_t ii = mapper.index(k,0,d);
                    if (ii >= 0 && ii < n) { ck[ii]=k; cd[ii]=d; }
                }

        // order the components by decreasing relative error, show the worst nShow
        std::vector<std::pair<real_t,index_t> > order;
        for (index_t i = 0; i != n; ++i)
        {
            const real_t scale = math::max(math::abs(bestFD(i)), math::abs(ga(i)));
            const real_t rel = (scale>1e-10) ? math::abs(bestFD(i)-ga(i))/scale : real_t(0);
            order.push_back(std::make_pair(rel,i));
        }
        std::sort(order.begin(),order.end());
        std::reverse(order.begin(),order.end());

        gsInfo << std::setw(6) << "i" << std::setw(6) << "k" << std::setw(4) << "d"
               << std::setw(20) << "anchor"
               << std::setw(16) << "g_analytic"
               << std::setw(16) << "g_FD"
               << std::setw(14) << "abs diff"
               << std::setw(12) << "rel diff" << "\n";
        const index_t nP = math::min(nShow,n);
        for (index_t j = 0; j != nP; ++j)
        {
            const index_t i = order[j].second;
            std::ostringstream anch;
            if (ck[i] >= 0)
            {
                anch << "(";
                for (index_t d = 0; d != dd; ++d)
                    anch << (d?",":"") << std::setprecision(3) << anchors(d,ck[i]);
                anch << ")";
            }
            else anch << "(?)";
            gsInfo << std::setw(6) << i << std::setw(6) << ck[i] << std::setw(4) << cd[i]
                   << std::setw(20) << anch.str()
                   << std::scientific << std::setprecision(6)
                   << std::setw(16) << ga(i)
                   << std::setw(16) << bestFD(i)
                   << std::setw(14) << math::abs(bestFD(i)-ga(i))
                   << std::setprecision(3) << std::setw(12) << order[j].first
                   << std::defaultfloat << std::setprecision(6) << "\n";
        }
        gsInfo << "    (" << nP << " worst of " << n << " components shown)\n";
    }
    else
        gsInfo << "\n[4] no h produced a fully finite FD gradient: per-component "
                  "table skipped.\n";

    // -----------------------------------------------------------------
    // 5. the same check at a few perturbed designs (x0 is the identity, a
    //    special point: a gradient bug can vanish there by symmetry)
    // -----------------------------------------------------------------
    for (index_t ip = 0; ip != nPert; ++ip)
    {
        // deterministic pseudo-random perturbation
        gsVector<real_t> xp(n);
        for (index_t i = 0; i != n; ++i)
        {
            const real_t s = math::sin(real_t(1000*(ip+1) + 7*i));
            xp(i) = x0(i) + pertMag*s;
        }
        const real_t fp = fEval(xp);
        gsVector<real_t> gap(n);
        {
            gsAsConstVector<real_t> xv(xp.data(), n);
            gsAsVector<real_t> gv(gap.data(), n);
            probe.problem().gradObj_into(xv, gv);
        }
        index_t nNanGp = 0;
        for (index_t i = 0; i != n; ++i) if (!isFin(gap(i))) ++nNanGp;

        // FD at two representative steps
        gsInfo << "\n[5." << (ip+1) << "] perturbed design, |dx|_inf = " << pertMag
               << ":  f = " << fp << " (" << classify(fp) << "), "
               << "non-finite g comps = " << nNanGp << "\n";
        const real_t hh[2] = {1e-5,1e-6};
        for (index_t k = 0; k != 2; ++k)
        {
            const real_t h = hh[k];
            real_t maxRel = 0; index_t argmax=-1, nBad=0;
            for (index_t i = 0; i != n; ++i)
            {
                gsVector<real_t> a = xp, b = xp;
                a(i) += h; b(i) -= h;
                const real_t fa = fEval(a), fb = fEval(b);
                if (!isFin(fa) || !isFin(fb) || !isFin(gap(i))) { ++nBad; continue; }
                const real_t fd = (fa-fb)/(real_t(2)*h);
                const real_t scale = math::max(math::abs(fd), math::abs(gap(i)));
                if (scale <= 1e-10) continue;
                const real_t rel = math::abs(fd-gap(i))/scale;
                if (rel > maxRel) { maxRel = rel; argmax = i; }
            }
            gsInfo << "        h = " << h << ": maxRel = " << maxRel
                   << " at i = " << argmax << ", nBad = " << nBad << "\n";
        }
    }

    // -----------------------------------------------------------------
    // 6. descent-line probe: f(x0 - t*g/||g||) for a logarithmic sweep of t.
    //    This is what a line search walks along.  If f is finite at x0 and the
    //    gradient is correct, a stall can only come from the objective going
    //    non-finite (or non-descending) somewhere along this line.
    // -----------------------------------------------------------------
    gsInfo << "\n[6] descent-line probe  x(t) = x0 - t * g/||g||\n";
    gsInfo << std::setw(12) << "t" << std::setw(22) << "f(x(t))"
           << std::setw(12) << "state"
           << std::setw(16) << "min det J_s"
           << std::setw(16) << "max|x-x0|_inf" << "\n";
    real_t firstBadT = -1;
    gsVector<real_t> xBad;
    if (nNanG == 0 && ga.norm() > 0)
    {
        const gsVector<real_t> dir = ga/ga.norm();
        const real_t ts[] = {1e-8,1e-6,1e-4,1e-3,1e-2,1e-1,0.25,0.5,1.0,2.0,5.0,1e1,1e2};
        for (index_t k = 0; k != 13; ++k)
        {
            const gsVector<real_t> xt = x0 - ts[k]*dir;
            const real_t ft = fEval(xt);
            real_t minDet;
            {
                gsAsConstVector<real_t> xv(xt.data(), n);
                minDet = probe.problem().computeMinJacobian(xv);
            }
            gsInfo << std::scientific << std::setprecision(2) << std::setw(12) << ts[k]
                   << std::setprecision(10) << std::setw(22) << ft
                   << std::setw(12) << classify(ft)
                   << std::setprecision(4) << std::setw(16) << minDet
                   << std::setw(16) << (xt-x0).cwiseAbs().maxCoeff()
                   << std::defaultfloat << std::setprecision(6) << "\n";
            if (!isFin(ft) && firstBadT < 0) { firstBadT = ts[k]; xBad = xt; }
        }
    }
    else
        gsInfo << "    (skipped: analytic gradient is not usable)\n";

    // -----------------------------------------------------------------
    // 7. term-by-term integrand scan, at x0 and (if any) at the first
    //    design on the descent line where the objective is non-finite.
    // -----------------------------------------------------------------
    gsInfo << "\n[7] term-by-term integrand scan at x0:\n";
    scanIntegrand<MODE>(probe,domain,geom,fun,parametric,x0);
    if (firstBadT > 0)
    {
        gsInfo << "\n[7b] term-by-term integrand scan at the FIRST non-finite point"
               << " on the descent line (t = " << firstBadT << "):\n";
        scanIntegrand<MODE>(probe,domain,geom,fun,parametric,xBad);
    }

    gsInfo << "\n    total non-finite objective evaluations in this mode: "
           << nNonFiniteEvals << "\n";

    domain.setControls(x0);   // leave the domain as we found it
}

int main(int arg, char *argv[])
{
    //! [Parse command line]
    index_t numRefineG=0, numRefineD=0, numElevateG=0, numElevateD=0;
    index_t mode = -1;              // -1: both modes
    index_t nRepeat = 10;
    index_t nPert = 2;
    index_t nShow = 15;
    index_t nGrid = 21;
    real_t  pertMag = 1e-2;
    real_t  penalty = -1;
    real_t  smoothing = -1;
    std::string input = "parametrization/monitor_example_planar.xml";

    gsCmdLine cmd("FD gradient check of the gsAdaptiveParametrization objective.");
    cmd.addReal("P","penalty","Override penalty parameter",penalty);
    cmd.addReal("S","smoothing","Override smoothing parameter",smoothing);
    cmd.addInt("r","numRefG","Uniform h-refinement loops for the geometry",numRefineG);
    cmd.addInt("e","numElevG","Degree elevation steps for the geometry",numElevateG);
    cmd.addInt("R","numRefD","Uniform h-refinement loops for the composition",numRefineD);
    cmd.addInt("E","numElevD","Degree elevation steps for the composition",numElevateD);
    cmd.addString("i","input","Input file",input);
    cmd.addInt("m","mode","Monitor mode: 0 ValueBased, 1 GradientBased, -1 both",mode);
    cmd.addInt("","repeat","Number of repeated f(x0) evaluations (noise floor)",nRepeat);
    cmd.addInt("","nPert","Number of perturbed designs to check",nPert);
    cmd.addInt("","nShow","Number of worst components to tabulate",nShow);
    cmd.addInt("","nGrid","1D grid size for the monitor-function scan",nGrid);
    cmd.addReal("","pert","Magnitude of the design perturbation",pertMag);
    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // ---------------------------------------------------------------
    // Read input -- identical to monitor_parameterization_example.cpp:66-121
    // ---------------------------------------------------------------
    gsFileData<> fd(input);
    gsMultiPatch<>      geom;
    gsMultiPatch<>      composition;
    gsFunctionExpr<>    function;
    gsOptionList        OPToptions;
    gsOptionList        PARoptions;

    fd.getLabel("geometry",geom);
    if (geom.nInterfaces()==0)
        geom.computeTopology();
    fd.getLabel("composition",composition);
    if (fd.hasLabel("function"))
        fd.getLabel("function",function);
    fd.getLabel("optimizer_options",OPToptions);
    fd.getLabel("parametrization_options",PARoptions);

    GISMO_ENSURE(geom.nPatches()==composition.nPatches(),
                 "The number of patches in the geometry and the composition must be the same");
    GISMO_ENSURE(geom.nPatches()==1,"Only one patch is supported");

    for (size_t p=0; p!=geom.nPatches(); p++)
    {
        geom.patch(p).degreeElevate(numElevateG);
        for (index_t i = 0; i < numRefineG; i++)
            geom.patch(p).uniformRefine();

        composition.patch(p).degreeElevate(numElevateD);
        for (index_t i = 0; i < numRefineD; i++)
            composition.patch(p).uniformRefine();
    }

    const gsBasis<>    & cbasis0 = composition.basis(0);
    const gsGeometry<> & gpatch0 = geom.patch(0);
    const gsBasis<>    & gbasis0 = gpatch0.basis();

    gsSquareDomain<real_t> domain(cbasis0);
    for (auto ifc = geom.iBegin(); ifc != geom.iEnd(); ++ifc)
        domain.addInterface(*ifc);
    domain.options().addSwitch("Slide","",PARoptions.askSwitch("Slide",false));
    domain.applyOptions();

    const bool   parametric = PARoptions.askSwitch("Parametric",false);
    const real_t pen  = (penalty  ==-1) ? PARoptions.askReal("Penalty",1e-2)   : penalty;
    const real_t smoo = (smoothing==-1) ? PARoptions.askReal("Smoothing",1e-2) : smoothing;

    gsInfo << "===============================================================\n";
    gsInfo << " monitor_gradcheck : " << input << "\n";
    gsInfo << "===============================================================\n";
    gsInfo << "geometry     : dd = " << gpatch0.domainDim()
           << ", td = " << gpatch0.targetDim()
           << ", deg = " << gbasis0.degree(0) << "x" << gbasis0.degree(1)
           << ", size = " << gbasis0.size() << "\n";
    gsInfo << "composition  : deg = " << cbasis0.degree(0) << "x" << cbasis0.degree(1)
           << ", size = " << cbasis0.size() << "\n";
    gsInfo << "case         : " << (gpatch0.domainDim()==gpatch0.targetDim()
                                    ? "PLANAR (td==dd)" : "SURFACE (td>dd)") << "\n";
    gsInfo << "nControls    : " << domain.nControls() << "\n";
    gsInfo << "Slide        : " << PARoptions.askSwitch("Slide",false) << "\n";
    gsInfo << "Parametric   : " << parametric << "\n";
    gsInfo << "Penalty      : " << pen << "\n";
    gsInfo << "Smoothing    : " << smoo << "\n";
    gsInfo << "monitor      : " << (function.domainDim()==0 ? "NONE"
                                    : "present, dim=") ;
    if (function.domainDim()!=0) gsInfo << function.domainDim();
    gsInfo << "\n";
#ifdef _OPENMP
    gsInfo << "OpenMP max threads : " << omp_get_max_threads() << "\n";
#else
    gsInfo << "OpenMP : disabled\n";
#endif

    GISMO_ENSURE(domain.nControls()>0,
                 "No free controls: refine the composition (-R/-E) or enable Slide.");

    const gsFunction<real_t> * funPtr =
        (function.domainDim()==0) ? nullptr : static_cast<const gsFunction<real_t>*>(&function);

    if (funPtr!=nullptr)
    {
        gsInfo << "\n--- monitor-function scan (value / gradient / Hessian) ---\n";
        deriv2WriteCoverage(function);
        scanMonitor(function,gpatch0,parametric,smoo,nGrid);
        checkEta2Gradient(domain,gpatch0,function,parametric,domain.getControls(),4);
    }

    // The optimizer is required by the constructor only; solve() is NEVER called.
    gsGradientDescent<real_t> optimizer;

    const gsVector<real_t> x0 = domain.getControls();

    if (mode==-1 || mode==MonitorMode::ValueBased)
        runCheck<MonitorMode::ValueBased>("ValueBased (control)",domain,gpatch0,funPtr,
                                          gbasis0,optimizer,parametric,pen,smoo,x0,
                                          nRepeat,nPert,pertMag,nShow);
    if (mode==-1 || mode==MonitorMode::GradientBased)
        runCheck<MonitorMode::GradientBased>("GradientBased",domain,gpatch0,funPtr,
                                             gbasis0,optimizer,parametric,pen,smoo,x0,
                                             nRepeat,nPert,pertMag,nShow);

    return 0;
}// end main
