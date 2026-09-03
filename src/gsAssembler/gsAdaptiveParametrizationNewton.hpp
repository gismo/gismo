/** @file gsAdaptiveParametrizationNewton.hpp

    @brief Newton and Picard solvers for composite spline relocation.

    Implementation. See gsAdaptiveParametrizationNewton.h for the interface
    and doc/derivation_hessian.md for the mathematics.

    Assembly strategy (hand-rolled, mirroring gsExprAssembler):
    - OpenMP element parallelism via gsDomain::allElements() (per-thread
      element chunks, see gsDomain.h).
    - Two-phase sparse storage: a cached gsFiberMatrix sparsity pattern
      (built once from the sigma-basis support graph), then numeric
      assembly accumulates each thread's contributions into a private
      value buffer (_fpatternPos maps a fixed pattern entry to its flat
      offset), merged into m_fpattern in THREAD-ID order after the
      parallel region (_mergePartialK) for bit-reproducibility, finally
      gsFiberMatrix::toSparseMatrix_into. Identical in spirit to the
      _pattern/_eval idiom of gsExprAssembler.

    Complexity per assembly: O(n_qp * (dd*k)^2) local work + scatter,
    k = number of active sigma basis functions per element.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-)
*/

#pragma once

#include <gsAssembler/gsAdaptiveParametrizationNewton.h>
#include <gsUtils/gsStopwatch.h>
#include <gsMatrix/gsSparseSolver.h>
#include <iomanip>
#include <fstream>

namespace gismo
{

// ---------------------------------------------------------------------------
// Constructors
// ---------------------------------------------------------------------------

template<class T, enum MonitorMode MODE>
gsAdaptiveParametrizationNewton<T,MODE>::gsAdaptiveParametrizationNewton(
        gsSquareDomain<T>  & composition,
    const gsGeometry<T>    & geometry,
    const gsBasis<T>       & integrationBasis)
:
    m_comp(&composition),
    m_geom(&geometry),
    m_fun(nullptr),
    m_ib(&integrationBasis),
    m_parametric(false),
    m_mb(integrationBasis, true),
    m_optMesh(composition, geometry, &integrationBasis),
    m_options(defaultOptions()),
    m_patternReady(false),
    m_geomDeriv3Mode(-1)
{
}

template<class T, enum MonitorMode MODE>
gsAdaptiveParametrizationNewton<T,MODE>::gsAdaptiveParametrizationNewton(
        gsSquareDomain<T>  & composition,
    const gsGeometry<T>    & geometry,
    const gsFunction<T>    & fun,
    const gsBasis<T>       & integrationBasis,
          bool               parametric)
:
    m_comp(&composition),
    m_geom(&geometry),
    m_fun(&fun),
    m_ib(&integrationBasis),
    m_parametric(parametric),
    m_mb(integrationBasis, true),
    m_optMesh(composition, geometry, &fun, &integrationBasis, parametric),
    m_options(defaultOptions()),
    m_patternReady(false),
    m_geomDeriv3Mode(-1)
{
}

template<class T, enum MonitorMode MODE>
gsOptionList gsAdaptiveParametrizationNewton<T,MODE>::defaultOptions()
{
    gsOptionList opt;
    opt.addReal("Penalty",    "Fold-barrier / Tikhonov regularisation", 1e-2);
    opt.addReal("Smoothing",  "Monitor smoothing parameter theta",      1e-2);
    opt.addReal("Tolerance",  "Newton/Picard residual convergence tol", 1e-8);
    opt.addInt ("MaxIter",    "Maximum Newton/Picard iterations",        50);
    opt.addReal("ArmijoC1",   "Armijo sufficient decrease constant",    1e-4);
    opt.addReal("ArmijoRho",  "Armijo step reduction factor",           0.5);
    opt.addInt ("ArmijoMax",  "Armijo maximum backtracks",               20);
    opt.addReal("LevMarTau0", "Initial Levenberg-Marquardt shift",      1e-4);
    opt.addSwitch("Verbose",  "Print per-iteration info",               true);
    opt.addString("LogFile",  "CSV log file path (empty = no file)",    "");
    return opt;
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::_getControls(gsVector<T> & u) const
{
    u = m_comp->getControls();
}

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::_setControls(const gsVector<T> & u) const
{
    m_comp->setControls(u);
}

template<class T, enum MonitorMode MODE>
T gsAdaptiveParametrizationNewton<T,MODE>::evalObj() const
{
    // Sync regularisation options (m_optMesh has its own defaults)
    m_optMesh.options().setReal("Penalty",   m_options.getReal("Penalty"));
    m_optMesh.options().setReal("Smoothing", m_options.getReal("Smoothing"));
    gsVector<T> u;
    _getControls(u);
    return m_optMesh.evalObj(gsAsConstVector<T>(u.data(), u.size()));
}

template<class T, enum MonitorMode MODE>
T gsAdaptiveParametrizationNewton<T,MODE>::computeMinJacobian() const
{
    gsVector<T> u;
    _getControls(u);
    return m_optMesh.computeMinJacobian(gsAsConstVector<T>(u.data(), u.size()));
}

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::_residual(
    const gsVector<T> & u, gsVector<T> & R) const
{
    const index_t n = u.size();
    R.resize(n);
    gsAsVector<T> aR(R.data(), n);
    m_optMesh.gradObj_into(gsAsConstVector<T>(u.data(), n), aR);
}

// ---------------------------------------------------------------------------
// Persistent per-thread scratch: sizing matrices that are written wholesale
// via .noalias() self-size on first use, but are pre-sized here anyway so the
// very first sweep does not allocate (mirrors gsOptMesh::EvalScratch::init).
template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::EvalScratch::init(index_t dd, index_t td)
{
    Js.setZero(dd,dd);  Jg.setZero(td,dd);  Jc.setZero(td,dd);
    C.setZero(dd,dd);   Cinv.setZero(dd,dd);
    Cg.setZero(dd,dd);  Cg_inv.setZero(dd,dd);
    grad_xi_f.setZero(dd,1);
    ready = true;
}

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::ResidualScratch::init(index_t dd, index_t td)
{
    Js.setZero(dd,dd);  Jg.setZero(td,dd);  Jc.setZero(td,dd);
    C.setZero(dd,dd);   Q.setZero(dd,dd);
    D_d.setZero(td,dd); Md.setZero(td,dd); A_d.setZero(dd,dd); QA.setZero(dd,dd);
    adj.setZero(dd,dd); AdjJg.setZero(dd,dd);
    b_d.setZero(dd);    Qb.setZero(dd);    Q2b.setZero(dd);
    Pi.setZero(dd,dd);
    r.setZero(dd);
    gA.setZero(dd);     v.setZero(dd);
    Cg.setZero(dd,dd);  Cg_inv.setZero(dd,dd);
    gradMon.setZero(dd,1); grad_xi_f.setZero(dd,1); grad_x_f.setZero(td,1);
    mon_scalar_d.setZero(dd);
    E_d.setZero(dd,dd); HfJgd.setZero(td,1); Dg_d.setZero(dd,1); JtHJd.setZero(dd,1);
    vs.setZero(2);
    dm2coef.setZero(dd);
    ready = true;
}

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::HessianScratch::init(index_t dd, index_t td)
{
    Js.setZero(dd,dd);  Jg.setZero(td,dd);  Jc.setZero(td,dd);
    C.setZero(dd,dd);   Q.setZero(dd,dd);
    Cg.setZero(dd,dd);  Cg_inv.setZero(dd,dd);

    D.resize(dd); Md.resize(dd); A.resize(dd); QAQ.resize(dd);
    b.resize(dd); Qb.resize(dd); Q2b.resize(dd);
    trCAC.setZero(dd); ad.setZero(dd);
    AdjJg.setZero(dd,dd);
    T3.resize(dd); E3.resize(dd);
    Ddm.resize(dd); cdm.resize(dd); Q2c.resize(dd);
    edm.resize(dd); amg.resize(dd); adjD.resize(dd);
    E1.resize(dd); gfm.resize(dd);
    for (index_t d = 0; d != dd; ++d)
    {
        D[d].setZero(td,dd); Md[d].setZero(td,dd); A[d].setZero(dd,dd);
        QAQ[d].setZero(dd,dd);
        b[d].setZero(dd); Qb[d].setZero(dd); Q2b[d].setZero(dd);
        T3[d].resize(dd); E3[d].resize(dd);
        Ddm[d].resize(dd); cdm[d].resize(dd); Q2c[d].resize(dd);
        edm[d].resize(dd); amg[d].resize(dd); adjD[d].resize(dd);
        for (index_t m = 0; m != dd; ++m)
        {
            T3[d][m].setZero(td,dd); E3[d][m].setZero(dd,dd);
            Ddm[d][m].setZero(td); cdm[d][m].setZero(dd); Q2c[d][m].setZero(dd);
            edm[d][m].setZero(dd); amg[d][m].setZero(dd); adjD[d][m].setZero(dd);
        }
        E1[d].setZero(dd,dd);
        gfm[d].setZero(dd);
    }
    M2.setZero(dd,dd); trQHQ.setZero(dd,dd); S1.setZero(dd,dd); Ss2.setZero(dd,dd);
    cJg.setZero(dd,dd); tT.setZero(dd,dd); bdMM.setZero(dd,dd);
    adj.setZero(dd,dd); adjMd.setZero(dd,dd); tmpM.setZero(td,dd);
    Htmp.setZero(dd,dd); QQ.setZero(dd,dd);
    gA.setZero(dd); gB.setZero(dd); qA.setZero(dd); qB.setZero(dd);
    v.setZero(dd); v_m.setZero(dd);
    gradMon.setZero(dd,1); grad_xi_f.setZero(dd,1); grad_x_f.setZero(td,1);
    ms.setZero(dd);
    dms.setZero(dd,dd); hEta.setZero(dd,dd);
    dgf.setZero(dd);
    dHxf.setZero(td,td);
    trCginvEd_vec.setZero(dd);
    trCginvEdm_mat.setZero(dd,dd); trCginvEmCginvEd_mat.setZero(dd,dd);
    dm2coef.setZero(dd); d2m2coef.setZero(dd,dd);
    ready = true;
}

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::PicardScratch::init(index_t dd, index_t td)
{
    Js.setZero(dd,dd);  Jg.setZero(td,dd);  Jc.setZero(td,dd);
    C.setZero(dd,dd);   Q.setZero(dd,dd);
    Cg.setZero(dd,dd);  Cg_inv.setZero(dd,dd);  Bc.setZero(dd,dd);
    grad_xi_f.setZero(dd,1);
    gA.setZero(dd); gB.setZero(dd); BcgB.setZero(dd);
    ready = true;
}

// _evalObjAndMinJ: fused single-sweep energy + min sigma-Jacobian.
//
// One OpenMP element pass computing simultaneously
//   energy = sum_q w_q * integrand(q)      (== gsOptMesh::evalObj)
//   minJ   = min_q det(J_sigma(q))         (== gsOptMesh::computeMinJacobian)
// The integrand mirrors evalObj exactly (planar/surface x monitor mode);
// det(J_sigma) is already needed there (surface) or otherwise nearly free, so
// the fold guard adds no extra sweep. Energy summation order matches evalObj's
// (same dom->allElements() order), so the single-thread result is bit-identical
// to the two separate routines.
// ---------------------------------------------------------------------------
template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::_evalObjAndMinJ(
    const gsVector<T> & u, T & energy, T & minJ) const
{
    const index_t dd = m_comp->domainDim();
    const index_t td = m_geom->targetDim();

    _setControls(u);

    const T penalty = m_options.getReal("Penalty");
    const T theta   = m_options.getReal("Smoothing");
    const bool hasMonitor = (m_fun != nullptr);

    energy = T(0);
    minJ   = std::numeric_limits<T>::max();

    gsOptionList quadOptions;
    quadOptions.addReal("quA","",1.0);
    quadOptions.addInt ("quB","",1);

    auto dom = m_mb.domain();

    // DETERMINISTIC reduction -- see gsOptMesh::evalObj (gsAdaptiveParametrization.hpp:425-432).
    // Per-thread partials are summed in THREAD-ID order, making this routine
    // bit-reproducible for a fixed thread count. minJ stays a plain
    // comparison: min is associative/commutative exactly in IEEE arithmetic,
    // so its reduction order never mattered; it is folded in the same loop
    // only for uniformity.
#   ifdef _OPENMP
    const int nThreads = omp_get_max_threads();
#   else
    const int nThreads = 1;
#   endif
    std::vector<T> partial(nThreads, T(0));
    std::vector<T> partialMin(nThreads, std::numeric_limits<T>::max());

#   pragma omp parallel
    {
        T thEnergy = T(0);
        T thMin    = std::numeric_limits<T>::max();

        // Persistent per-thread scratch (see EvalScratch in the header):
        // buffers, gsFuncData and the Gauss rule are allocated on the first
        // call and reused by every later evaluation of this solve.
        EvalScratch & sc = m_evalScratch.mine();
        if (!sc.ready) sc.init(dd, td);

        gsFuncData<T> & compData = sc.compData;
        gsFuncData<T> & geomData = sc.geomData;
        gsFuncData<T> & funData  = sc.funData;
        compData.flags = NEED_VALUE | NEED_DERIV;
        geomData.flags = NEED_VALUE | NEED_DERIV;
        if (hasMonitor)
        {
            if (MODE == ValueBased)         funData.flags = NEED_VALUE;
            else /* GradientBased */        funData.flags = NEED_VALUE | NEED_DERIV;
        }

        gsMatrix<T> & Js = sc.Js, & Jg = sc.Jg, & Jc = sc.Jc;
        gsMatrix<T> & C = sc.C, & Cinv = sc.Cinv, & Cg = sc.Cg, & Cg_inv = sc.Cg_inv;
        gsMatrix<T> & grad_xi_f = sc.grad_xi_f;

        gsMatrix<T> & uvPoints = sc.uvPoints;
        gsMatrix<T> & Jsigma_flat = sc.Jsigma_flat, & Jgeom_flat = sc.Jgeom_flat;
        gsMatrix<T> & monVals = sc.monVals, & monDerivs = sc.monDerivs;
        gsVector<T> & tmpWeights = sc.tmpWeights;

        for (auto & elem : dom->allElements())
        {
            if (sc.QuPatch != elem.patch())
            {
                sc.QuPatch = elem.patch();
                sc.QuRule = gsQuadrature::get(m_ib->basis(sc.QuPatch), quadOptions);
            }

            sc.QuRule.mapTo(elem.lowerCorner(), elem.upperCorner(), uvPoints, tmpWeights);

            m_comp->compute(uvPoints, compData);
            Jsigma_flat = compData.values[1];
            m_geom->compute(compData.values[0], geomData);
            Jgeom_flat = geomData.values[1];

            if (hasMonitor)
            {
                m_fun->compute((m_parametric ? compData.values[0] : geomData.values[0]), funData);
                monVals = funData.values[0];
                if (MODE == GradientBased)
                    monDerivs = funData.values[1];
            }

            const index_t nPts = uvPoints.cols();
            for (index_t p = 0; p != nPts; ++p)
            {
                Js.noalias() = Jsigma_flat.col(p).reshaped(dd, dd).transpose();
                Jg.noalias() = Jgeom_flat.col(p).reshaped(dd, td).transpose();
                Jc.noalias() = Jg * Js;
                C.noalias() = Jc.transpose() * Jc;
                // NO Tikhonov shift on C (matches gsOptMesh's unregularised
                // energy); C is inverted unregularised, so a folded iterate
                // returns NaN here. Rationale and the fold hazard: \ref adaptparam_newton.
                Cinv.noalias() = C.inverse();

                // Fold guard: min sigma-Jacobian det (same quantity as
                // gsOptMesh::computeMinJacobian).
                const T detJs = Js.determinant();
                if (detJs < thMin) thMin = detJs;

                // Monitor weight m2 (only when a monitor is present)
                // The paper (see the relevant equation): omega = 1/sqrt(1+theta*f) or 1/sqrt(1+theta*||∇f||^2)
                // so m2 = omega^2 = 1/(1+theta*f) or 1/(1+theta*||∇f||^2)
                T m2 = T(0);
                if (hasMonitor)
                {
                    if (MODE == ValueBased)
                    {
                        const T eta = monVals(0, p);
                        m2 = T(1) / (T(1) + theta * eta);  // per the value-based weight formula, see \ref adaptparam_rstep
                    }
                    else
                    {
                        Cg.noalias() = Jg.transpose() * Jg;
                        Cg_inv.noalias() = Cg.inverse();
                        if (m_parametric)
                            grad_xi_f.noalias() = monDerivs.col(p);
                        else
                            grad_xi_f.noalias() = Jg.transpose() * monDerivs.col(p).reshaped(td, 1);
                        const T eta2 = (grad_xi_f.transpose() * Cg_inv * grad_xi_f)(0,0);
                        m2 = T(1) / (T(1) + theta * eta2);  // per the gradient-based weight formula, see \ref adaptparam_rstep
                    }
                }

                // Barrier: planar uses det(J_c); surface uses det(J_s)=detJs.
                T abs_det, chi;
                if (td == dd)
                {
                    const T det_c = Jc.determinant();
                    chi = T(0.5) * (det_c + math::sqrt(penalty*penalty + det_c*det_c));
                    abs_det = math::abs(det_c);
                }
                else
                {
                    chi = T(0.5) * (detJs + math::sqrt(penalty*penalty + detJs*detJs));
                    abs_det = math::abs(detJs);
                }

                T integrand;
                if (hasMonitor)   // omega = sqrt(m2); see q_exp note below
                {
                    const T omega = math::sqrt(m2);
                    // q_exp = 1 on chi (monitor branch): with the Tikhonov
                    // shift gone, tr(Cinv)*det^2 == T exactly, so the barrier
                    // exponent belongs on chi, not on |det|, to reproduce
                    // gsOptMesh's omega*T/chi.
                    integrand = omega * Cinv.trace() * abs_det*abs_det / chi;
                    // Surface (td>dd) paper-exact area weight g_S = sqrt(det Cg) (Eq.18).
                    // Cg already formed above for GradientBased; compute for ValueBased.
                    if (td > dd)
                    {
                        if (MODE == ValueBased)
                            Cg.noalias() = Jg.transpose() * Jg;
                        integrand *= math::sqrt(Cg.determinant());
                    }
                }
                else              // p=2 barrier exponent
                    integrand = Cinv.trace() * abs_det*abs_det / (chi*chi);

                thEnergy += tmpWeights[p] * integrand;
            }
        }
#       ifdef _OPENMP
        partial[omp_get_thread_num()]    = thEnergy;
        partialMin[omp_get_thread_num()] = thMin;
#       else
        partial[0]    = thEnergy;
        partialMin[0] = thMin;
#       endif
    } // omp parallel

    for (int i = 0; i != nThreads; ++i)
    {
        energy += partial[i];
        if (partialMin[i] < minJ) minJ = partialMin[i];
    }
}

// ---------------------------------------------------------------------------
// Sparsity pattern (phase 1 of the two-phase assembly, cf. gsExprAssembler
// ::_computePattern). All pairs of free DoFs whose sigma basis functions
// share an element, for all component pairs (d,m). Built once and cached:
// the support graph does not change between iterations.
// ---------------------------------------------------------------------------

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::_buildPattern() const
{
    if (m_patternReady) return;

    const index_t dd = m_comp->domainDim();
    const gsDofMapper & mapper = m_comp->mapper();
    const index_t N = (index_t)m_comp->nControls();
    const gsBasis<T> & sigmaBasis = m_comp->domain().basis();

    m_fpattern.resize(N, N);

    auto dom = m_mb.domain();
    gsMatrix<index_t> actives;
    gsMatrix<T> center;

    // Serial sweep: one active_into per element (cheap; done once).
    auto domIt    = dom->beginAll();
    auto domItEnd = dom->endAll();
    for (; domIt < domItEnd; ++domIt)
    {
        center = domIt.centerPoint();
        sigmaBasis.active_into(center, actives);
        const index_t nAct = actives.rows();

        for (index_t l = 0; l != nAct; ++l)
        {
            const index_t li = actives(l, 0);
            for (index_t m = 0; m != dd; ++m)
            {
                if (!mapper.is_free(li, 0, m)) continue;
                const index_t jj = mapper.index(li, 0, m);
                for (index_t k = 0; k != nAct; ++k)
                {
                    const index_t ki = actives(k, 0);
                    for (index_t d = 0; d != dd; ++d)
                    {
                        if (!mapper.is_free(ki, 0, d)) continue;
                        const index_t ii = mapper.index(ki, 0, d);
                        m_fpattern.insertExplicitZero(ii, jj);
                    }
                }
            }
        }
    }

    // Prefix sum of nonzeros per column, so a fixed pattern entry maps to a
    // unique offset in a flat per-thread value buffer (see _fpatternPos).
    const gsVector<index_t> nzPerFiber = m_fpattern.nonZerosPerFiber();
    m_fpatternColOffset.resize(N);
    index_t off = 0;
    for (index_t j = 0; j != N; ++j)
    {
        m_fpatternColOffset[j] = off;
        off += nzPerFiber[j];
    }

    m_patternReady = true;
}

template<class T, enum MonitorMode MODE>
index_t gsAdaptiveParametrizationNewton<T,MODE>::_fpatternPos(index_t ii, index_t jj) const
{
    // Column-major layout: m_fpattern's fiber jj stores the (sorted) row
    // indices of column jj; searchLowerIndex is the same lower_bound used
    // internally by gsFiberMatrix::isExplicitZero.
    const auto & vdata = m_fpattern.col(jj).data();
    const index_t local = vdata.searchLowerIndex(ii);
    GISMO_ASSERT(local < vdata.size() && vdata.index(local) == ii,
                 "Pattern entry ("<<ii<<","<<jj<<") missing from m_fpattern; "
                 "did _buildPattern() run?");
    if (local >= vdata.size() || vdata.index(local) != ii)
        return -1;
    return m_fpatternColOffset[jj] + local;
}

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::_mergePartialK(
    const std::vector<const gsVector<T> *> & partial, gsSparseMatrix<T> & out) const
{
    const index_t N = m_fpattern.cols();
    for (index_t j = 0; j != N; ++j)
    {
        auto & col = m_fpattern.col(j);
        const index_t nzj = col.nonZeros();
        if (nzj == 0) continue;
        gsAsVector<T> dst(col.valuePtr(), nzj);
        for (size_t t = 0; t != partial.size(); ++t)
            dst += partial[t]->segment(m_fpatternColOffset[j], nzj);
    }
    m_fpattern.toSparseMatrix_into(out);
    out.makeCompressed();
}

// ---------------------------------------------------------------------------
// Third derivatives of the geometry S at geometry-parametric points.
//
// Output layout ("derivative of deriv2"): with S2 = dd*(dd+1)/2,
//   out((m*td + a)*S2 + hidx(i,j), p) = d^3 S_a / dxi_i dxi_j dxi_m
// (hidx = the gismo deriv2 index of the unordered pair (i,j)).
// Redundant under permutations of (i,j,m) — harmless, simplifies indexing.
//
// Exact path: evalAllDers_into(u, 3, .) (tensor B-spline bases). The exact
// result[3] uses composition layout: for dd==2 row 4*a + c1 with c1 = number
// of derivatives w.r.t. xi_1; converted below.
// FD path (bases without third derivatives, e.g. NURBS): central differences
// of deriv2 with h=1e-5, clamped to the support box (one-sided near the
// boundary). Truncation error O(h^2) ~ 1e-10.
// ---------------------------------------------------------------------------

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::_geomDeriv3_into(
    const gsMatrix<T> & pts, gsMatrix<T> & out) const
{
    const index_t dd = m_geom->domainDim();
    const index_t td = m_geom->targetDim();
    const index_t S2 = dd * (dd + 1) / 2;
    const index_t nPts = pts.cols();

    out.resize(td * S2 * dd, nPts);

    if (m_geomDeriv3Mode == -1)
    {
        // Probe once whether the basis provides exact third derivatives.
        // GISMO_ERROR prints to std::cerr before throwing, so redirect cerr
        // to a discard buffer to suppress the spurious "order up to 2<3" output.
        std::ostringstream devnull;
        std::streambuf * saved_cerr = std::cerr.rdbuf(devnull.rdbuf());
        std::vector<gsMatrix<T>> ders;
        try
        {
            gsMatrix<T> probe = pts.col(0);
            m_geom->evalAllDers_into(probe, 3, ders);
            m_geomDeriv3Mode = 1;
        }
        catch (...)
        {
            m_geomDeriv3Mode = 0;
            std::cerr.rdbuf(saved_cerr);
            gsInfo << "gsAdaptiveParametrizationNewton: geometry has no exact "
                      "third derivatives; using FD-of-deriv2 (h=1e-5).\n";
        }
        std::cerr.rdbuf(saved_cerr);
    }

    if (m_geomDeriv3Mode == 1)
    {
        GISMO_ASSERT(dd == 2, "Exact third-derivative conversion implemented for dd==2 only");
        std::vector<gsMatrix<T>> ders;
        m_geom->evalAllDers_into(pts, 3, ders);
        const gsMatrix<T> & d3 = ders[3]; // rows: 4*a + c1, c1 = #ones in triple
        // hidx for dd==2: (0,0)->0, (1,1)->1, (0,1)->2; pair "ones" count:
        // hidx 0 -> 0 ones, hidx 1 -> 2 ones, hidx 2 -> 1 one.
        const index_t pairOnes[3] = {0, 2, 1};
        for (index_t m = 0; m != dd; ++m)
            for (index_t a = 0; a != td; ++a)
                for (index_t h = 0; h != S2; ++h)
                    out.row((m * td + a) * S2 + h) = d3.row(4 * a + pairOnes[h] + m);
    }
    else
    {
        // FD of deriv2 in each parametric direction m, clamped to support.
        gsMatrix<T> supp = m_geom->support();
        const T h = T(1e-5);
        gsMatrix<T> ptsP, ptsM, d2P, d2M;
        for (index_t m = 0; m != dd; ++m)
        {
            ptsP = pts; ptsM = pts;
            ptsP.row(m).array() += h;
            ptsM.row(m).array() -= h;
            if (supp.cols() >= 2)
                for (index_t p = 0; p != nPts; ++p)
                {
                    ptsP(m,p) = math::min(ptsP(m,p), supp(m,1));
                    ptsM(m,p) = math::max(ptsM(m,p), supp(m,0));
                }
            m_geom->deriv2_into(ptsP, d2P);
            m_geom->deriv2_into(ptsM, d2M);
            for (index_t p = 0; p != nPts; ++p)
            {
                const T den = ptsP(m,p) - ptsM(m,p);
                out.block(m * td * S2, p, td * S2, 1) =
                    (d2P.col(p) - d2M.col(p)) / den;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Third derivatives of the monitor f (same layout as _geomDeriv3_into, with
// td -> 1 and dd -> input dimension of f). Always FD of deriv2: generic
// gsFunction (e.g. gsFunctionExpr) provides at most n=2.
// ---------------------------------------------------------------------------

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::_funDeriv3_into(
    const gsMatrix<T> & pts, gsMatrix<T> & out) const
{
    GISMO_ASSERT(m_fun != nullptr, "No monitor function");
    const index_t df = pts.rows();             // input dimension of f
    const index_t S2 = df * (df + 1) / 2;
    const index_t nPts = pts.cols();

    out.resize(S2 * df, nPts);

    gsMatrix<T> supp = m_fun->support();
    const T h = T(1e-5);
    gsMatrix<T> ptsP, ptsM, d2P, d2M;
    for (index_t m = 0; m != df; ++m)
    {
        ptsP = pts; ptsM = pts;
        ptsP.row(m).array() += h;
        ptsM.row(m).array() -= h;
        if (supp.rows() >= df && supp.cols() >= 2)
            for (index_t p = 0; p != nPts; ++p)
            {
                ptsP(m,p) = math::min(ptsP(m,p), supp(m,1));
                ptsM(m,p) = math::max(ptsM(m,p), supp(m,0));
            }
        m_fun->deriv2_into(ptsP, d2P);
        m_fun->deriv2_into(ptsM, d2M);
        for (index_t p = 0; p != nPts; ++p)
        {
            const T den = ptsP(m,p) - ptsM(m,p);
            out.block(m * S2, p, S2, 1) = (d2P.col(p) - d2M.col(p)) / den;
        }
    }
}

// ---------------------------------------------------------------------------
// assembleResidual: independent weak-form assembly of the EL residual
//
//   R_A^(d) = ∫ Pi_d · ∇N_A + r_d N_A  dΩ̂                       [H-3.2]
//
// Flux Pi_d and reaction r_d are the gradNk / Nk coefficients of the
// per-quad-point gradient formula of gsOptMesh::gradObj_into:
//
//   Planar (td==dd), phi = |det_c|^p / chi^2, gamma = p/det - 2 chi'/chi:
//     Pi_d = m2 * ( -2 phi Q^2 b_d + trQ phi gamma (adj(Jc) Jg).col(d) )
//     r_d  = m2 * ( -phi tr(Q A_d Q) + trQ phi gamma tr(adj(Jc) D_d Js) )
//            + dm2_d * trQ * phi
//   Surface (td>dd), det_s = det(Js):
//     Pi_d = m2 * ( -2 phi Q^2 b_d + trQ phi gamma adj(Js).col(d) )
//     r_d  = m2 * ( -phi tr(Q A_d Q) ) + dm2_d * trQ * phi
//
// with Q = C_eps^{-1}, b_d = Jc^T Jg.col(d), A_d = Md^T Jc + Jc^T Md,
// Md = D_d Js, and dm2_d the monitor-weight derivative (reaction only).
// Verified == gradObj_into to machine precision in unittests.
// ---------------------------------------------------------------------------

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::assembleResidual(
    const gsVector<T> & u, gsVector<T> & R_out) const
{
    _setControls(u);
    m_optMesh.options().setReal("Penalty",   m_options.getReal("Penalty"));
    m_optMesh.options().setReal("Smoothing", m_options.getReal("Smoothing"));

    const index_t dd = m_comp->domainDim();
    const index_t td = m_geom->targetDim();
    GISMO_ASSERT(dd <= td, "domainDim must be <= targetDim");
    GISMO_ASSERT(dd == 2, "assembleResidual implemented for dd==2");

    const T penalty = m_options.getReal("Penalty");
    const T theta   = m_options.getReal("Smoothing");
    const bool hasMonitor = (m_fun != nullptr);
    const bool planar = (td == dd);
    // q_exp: barrier exponent on chi (not |det|), reproducing gsOptMesh's
    // T/chi^2 (no monitor) / omega*T/chi (with monitor). See \ref adaptparam_newton.
    const T q_exp = hasMonitor ? T(1) : T(2);   // [H-2]

    const gsDofMapper & mapper = m_comp->mapper();
    const index_t N = (index_t)m_comp->nControls();
    const gsBasis<T> & sigmaBasis = m_comp->domain().basis();
    const index_t S2 = dd * (dd + 1) / 2;

    R_out.resize(N);
    R_out.setZero();

    gsOptionList quadOptions;
    quadOptions.addReal("quA","",1.0);
    quadOptions.addInt ("quB","",1);

    auto dom = m_mb.domain();

    // DETERMINISTIC reduction -- see gsOptMesh::evalObj (gsAdaptiveParametrization.hpp:425-432).
    // Per-thread partials are summed in THREAD-ID order, making this routine
    // bit-reproducible for a fixed thread count.
#   ifdef _OPENMP
    const int nThreads = omp_get_max_threads();
#   else
    const int nThreads = 1;
#   endif
    std::vector<gsVector<T> > partial(nThreads, gsVector<T>::Zero(N));

#   pragma omp parallel
    {
        gsVector<T> thR(N);
        thR.setZero();

        // Persistent per-thread scratch (see ResidualScratch in the header).
        ResidualScratch & sc = m_residualScratch.mine();
        if (!sc.ready) sc.init(dd, td);

        gsFuncData<T> & geomData = sc.geomData;
        gsFuncData<T> & funData  = sc.funData;
        gsFuncData<T> & sigmaData = sc.sigmaData;
        gsFuncData<T> & sigmaBasisData = sc.sigmaBasisData;
        geomData.flags = NEED_VALUE | NEED_DERIV | NEED_DERIV2;
        funData.flags  = NEED_VALUE | NEED_DERIV;
        if (hasMonitor && MODE == GradientBased)
            funData.flags |= NEED_DERIV2;
        sigmaData.flags = NEED_VALUE | NEED_DERIV;
        // SAME_ELEMENT: all quad points of an element share the active set,
        // so actives are computed once per element (cf. gsExprAssembler).
        sigmaBasisData.flags = NEED_ACTIVE | NEED_VALUE | NEED_DERIV | SAME_ELEMENT;

        gsMatrix<T> & uvPoints = sc.uvPoints;
        gsVector<T> & tmpWeights = sc.tmpWeights;
        gsMatrix<T> & Jsigma_flat = sc.Jsigma_flat, & Jgeom_flat = sc.Jgeom_flat;
        gsMatrix<T> & deriv2_geom = sc.deriv2_geom;
        gsMatrix<T> & monVals = sc.monVals, & monDerivs = sc.monDerivs, & monDeriv2 = sc.monDeriv2;
        gsMatrix<index_t> & actives = sc.actives;
        gsMatrix<T> & basisVals = sc.basisVals, & basisDerivs = sc.basisDerivs;

        gsMatrix<T> & Js = sc.Js, & Jg = sc.Jg, & Jc = sc.Jc, & C = sc.C, & Q = sc.Q;
        gsMatrix<T> & D_d = sc.D_d, & Md = sc.Md, & A_d = sc.A_d, & QA = sc.QA;
        gsMatrix<T> & adj = sc.adj, & AdjJg = sc.AdjJg;
        gsVector<T> & b_d = sc.b_d, & Qb = sc.Qb, & Q2b = sc.Q2b;
        gsMatrix<T> & Pi = sc.Pi;          // flux: column d = Pi_d
        gsVector<T> & r = sc.r;            // reaction
        gsVector<T> & gA = sc.gA, & v = sc.v;
        gsMatrix<T> & Cg = sc.Cg, & Cg_inv = sc.Cg_inv;
        gsMatrix<T> & gradMon = sc.gradMon, & grad_xi_f = sc.grad_xi_f, & grad_x_f = sc.grad_x_f;
        gsMatrix<T> & Hess_f = sc.Hess_f;
        gsVector<T> & mon_scalar_d = sc.mon_scalar_d;
        // GradientBased / surface scratch (reused per quad point)
        gsMatrix<T> & E_d = sc.E_d, & HfJgd = sc.HfJgd, & Dg_d = sc.Dg_d, & JtHJd = sc.JtHJd;
        gsVector<T> & vs = sc.vs;

        for (auto & elem : dom->allElements())
        {
            if (sc.QuPatch != elem.patch())
            {
                sc.QuPatch = elem.patch();
                sc.QuRule = gsQuadrature::get(m_ib->basis(sc.QuPatch), quadOptions);
            }

            sc.QuRule.mapTo(elem.lowerCorner(), elem.upperCorner(), uvPoints, tmpWeights);

            m_comp->compute(uvPoints, sigmaData);
            Jsigma_flat = sigmaData.values[1];

            sigmaBasis.compute(uvPoints, sigmaBasisData);
            actives     = sigmaBasisData.actives;
            basisVals   = sigmaBasisData.values[0];
            basisDerivs = sigmaBasisData.values[1];
            const index_t nAct = actives.rows();

            m_geom->compute(sigmaData.values[0], geomData);
            Jgeom_flat  = geomData.values[1];
            deriv2_geom = geomData.values[2];

            if (hasMonitor)
            {
                m_fun->compute((m_parametric ? sigmaData.values[0] : geomData.values[0]), funData);
                monVals   = funData.values[0];
                monDerivs = funData.values[1];
                if (MODE == GradientBased)
                    monDeriv2 = funData.values[2];
            }

            const index_t nPts = uvPoints.cols();
            for (index_t p = 0; p != nPts; ++p)
            {
                const T w_q = tmpWeights[p];

                Js.noalias() = Jsigma_flat.col(p).reshaped(dd, dd).transpose();
                Jg.noalias() = Jgeom_flat.col(p).reshaped(dd, td).transpose();
                Jc.noalias() = Jg * Js;

                C.noalias() = Jc.transpose() * Jc;
                // NO Tikhonov shift on C (matches gsOptMesh's unregularised
                // energy); C is inverted unregularised, so a folded iterate
                // returns NaN here. Rationale and the fold hazard: \ref adaptparam_newton.
                Q.noalias() = C.inverse();
                const T trQ = Q.trace();

                // Barrier factor phi_p, gamma = d(ln phi)/d(det)   [H-4.1b]
                T det, chi, dchi, phi, gamma;
                if (planar) det = Jc.determinant();
                else        det = Js.determinant();
                {
                    const T sq = math::sqrt(penalty*penalty + det*det);
                    chi  = T(0.5) * (det + sq);
                    dchi = T(0.5) * (T(1) + det / sq);
                    // q_exp is always 1 or 2 (see its definition above); avoid a
                    // libm pow() call per quadrature point -- hasMonitor is
                    // loop-invariant, so this branch is predictable and free.
                    phi = hasMonitor ? (det * det / chi) : (det * det / (chi * chi));
                    const T inv_det = (det != T(0)) ? T(1)/det : T(0);
                    gamma = T(2) * inv_det - q_exp * dchi / chi;
                }

                if (planar)
                {
                    // adj(Jc), 2x2: adj*Jc = det*I
                    adj(0,0) =  Jc(1,1); adj(0,1) = -Jc(0,1);
                    adj(1,0) = -Jc(1,0); adj(1,1) =  Jc(0,0);
                    AdjJg.noalias() = adj * Jg;
                }

                // Monitor weight m2 = omega^2 and its alpha-derivative coefficient
                //  omega = 1/sqrt(1+theta*f) => m2 = 1/(1+theta*f)
                T m2 = T(1);
                gsVector<T> & dm2coef = sc.dm2coef;    // dm2_X = dm2coef(d) * N_A
                dm2coef.setZero();
                if (hasMonitor)
                {
                    if (MODE == ValueBased)
                    {
                        const T eta = monVals(0, p);
                        const T denom = T(1) + theta * eta;
                        m2 = T(1) / denom;
                        const T dm2_deta = -theta / (denom * denom);
                        if (m_parametric)
                            gradMon.noalias() = monDerivs.col(p);
                        else
                            gradMon.noalias() = Jg.transpose() * monDerivs.col(p).reshaped(td, 1);
                        for (index_t d = 0; d != dd; ++d)
                            dm2coef(d) = dm2_deta * gradMon(d, 0);
                    }
                    else
                    {
                        // GradientBased: eta^2 = grad_xi_f^T Cg^{-1} grad_xi_f
                        Cg.noalias() = Jg.transpose() * Jg;
                        Cg_inv.noalias() = Cg.inverse();

                        const index_t hd = m_parametric ? dd : td;
                        Hess_f.resize(hd, hd);
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

                        const T eta2 = (grad_xi_f.transpose() * Cg_inv * grad_xi_f)(0,0);
                        const T denom = T(1) + theta * eta2;
                        m2 = T(1) / denom;
                        const T dm2_deta2 = -theta / (denom * denom);

                        v.noalias() = Cg_inv * grad_xi_f;

                        // mon_scalar_d(d) = d(eta^2)/d(xi_d) (cf. gradObj_into)
                        for (index_t d = 0; d != dd; ++d)
                        {
                            for (index_t a = 0; a != td; ++a)
                                for (index_t j = 0; j != dd; ++j)
                                {
                                    const index_t lo = math::min(d, j), hi = math::max(d, j);
                                    const index_t hess_idx = (lo == hi) ? lo : dd + lo*(2*dd-lo-3)/2 + hi - 1;
                                    D_d(a, j) = deriv2_geom(a * S2 + hess_idx, p);
                                }
                            E_d.noalias() = D_d.transpose() * Jg;
                            E_d += Jg.transpose() * D_d;
                            const T vEv = (v.transpose() * E_d * v)(0,0);
                            if (m_parametric)
                                mon_scalar_d(d) = -vEv + T(2) * (v.transpose() * Hess_f.col(d))(0,0);
                            else
                            {
                                HfJgd.noalias() = Hess_f * Jg.col(d);
                                Dg_d.noalias()  = D_d.transpose() * grad_x_f;
                                JtHJd.noalias() = Jg.transpose() * HfJgd;
                                mon_scalar_d(d) = -vEv + T(2) * (v.transpose() * (Dg_d + JtHJd))(0,0);
                            }
                            dm2coef(d) = dm2_deta2 * mon_scalar_d(d);
                        }
                    }
                }

                // Surface paper-exact area weight g_S = sqrt(det Cg), Cg = Jg^T Jg (Eq.18).
                // Cg/Cg_inv already formed for GradientBased above; compute for ValueBased.
                // Planar (td==dd) does not need g_S.
                T gS = T(1);
                if (hasMonitor && !planar)
                {
                    if (MODE == ValueBased)
                    {
                        Cg.noalias()     = Jg.transpose() * Jg;
                        Cg_inv.noalias() = Cg.inverse();
                    }
                    gS = math::sqrt(Cg.determinant());
                }

                // Flux Pi and reaction r per component d
                for (index_t d = 0; d != dd; ++d)
                {
                    for (index_t a = 0; a != td; ++a)
                        for (index_t j = 0; j != dd; ++j)
                        {
                            const index_t lo = math::min(d, j), hi = math::max(d, j);
                            const index_t hess_idx = (lo == hi) ? lo : dd + lo*(2*dd-lo-3)/2 + hi - 1;
                            D_d(a, j) = deriv2_geom(a * S2 + hess_idx, p);
                        }
                    Md.noalias()  = D_d * Js;
                    A_d.noalias() = Md.transpose() * Jc + Jc.transpose() * Md;
                    b_d.noalias() = Jc.transpose() * Jg.col(d);
                    QA.noalias()  = Q * A_d;
                    const T trCAC = (QA * Q).trace();
                    Qb.noalias()  = Q * b_d;
                    Q2b.noalias() = Q * Qb;

                    Pi.col(d) = -T(2) * phi * Q2b;
                    r(d)      = -phi * trCAC;

                    if (planar)
                    {
                        // d(det_c) flux/reaction:  u_d = adj*Jg.col(d),
                        //                          a_d = tr(adj*Md)
                        Pi.col(d) += trQ * phi * gamma * AdjJg.col(d);
                        r(d)      += trQ * phi * gamma * (adj * Md).trace();
                    }
                    else
                    {
                        // d(det_s) is flux-only: adj(Js).col(d)
                        // adj(Js) 2x2:
                        if (d == 0) { vs(0) =  Js(1,1); vs(1) = -Js(1,0); }
                        else        { vs(0) = -Js(0,1); vs(1) =  Js(0,0); }
                        Pi.col(d) += trQ * phi * gamma * vs;
                    }

                    // integrand = omega * (T/g), where omega = sqrt(m2)
                    // So flux and reaction scaled by omega, not m2
                    const T omega = math::sqrt(m2);
                    Pi.col(d) *= omega;
                    // d(omega)/d(alpha) = dm2/d(alpha) / (2*omega)
                    const T domega_coef = (m2 > 0) ? (dm2coef(d) / (T(2) * omega)) : T(0);
                    r(d) = omega * r(d) + domega_coef * trQ * phi;

                    // Surface (td>dd): paper-exact area weight e = A * g_S, with the current
                    // (omega-scaled) Pi,r = the flux/reaction of A = omega*trQ*phi.
                    //   d(A*g_S) = dA*g_S + A*dg_S; dg_S/dalpha = 0.5*g_S*tr(Cg^{-1}E_d)*N_A
                    //   (a pure reaction term). E_d = D_d^T Jg + Jg^T D_d.
                    if (hasMonitor && !planar)
                    {
                        E_d.noalias() = D_d.transpose() * Jg;
                        E_d += Jg.transpose() * D_d;
                        const T trCgE = (Cg_inv * E_d).trace();
                        const T A = omega * trQ * phi;   // integrand value (before g_S)
                        Pi.col(d) *= gS;
                        r(d) = gS * r(d) + A * T(0.5) * gS * trCgE;
                    }
                }

                // Scatter: R(ii) += w * (Pi_d . gA + r_d * N_A)
                for (index_t loc = 0; loc != nAct; ++loc)
                {
                    const index_t k = actives(loc, 0);
                    const T NA = basisVals(loc, p);
                    for (index_t j = 0; j != dd; ++j)
                        gA(j) = basisDerivs(loc * dd + j, p);

                    for (index_t d = 0; d != dd; ++d)
                    {
                        if (!mapper.is_free(k, 0, d)) continue;
                        const index_t ii = mapper.index(k, 0, d);
                        thR(ii) += w_q * (Pi.col(d).dot(gA) + r(d) * NA);
                    }
                }
            }
        }

#       ifdef _OPENMP
        partial[omp_get_thread_num()] = thR;
#       else
        partial[0] = thR;
#       endif
    } // omp parallel

    for (int i = 0; i != nThreads; ++i)
        R_out += partial[i];
}

// ---------------------------------------------------------------------------
// _hessianFD: central differences of the analytic gradient (h=1e-7).
// Validation oracle for the analytic Hessian and fallback for dd==3.
// O(n) gradient evaluations — slow, use for testing only.
// ---------------------------------------------------------------------------

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::_hessianFD(
    const gsVector<T> & u, gsSparseMatrix<T> & K_out) const
{
    const index_t n = u.size();
    const T h = T(1e-7);

    gsVector<T> grad_p(n), grad_m(n);
    gsVector<T> uu = u;

    gsMatrix<T> K_dense(n, n);

    for (index_t j = 0; j < n; ++j)
    {
        uu[j] = u[j] + h;
        _residual(uu, grad_p);
        uu[j] = u[j] - h;
        _residual(uu, grad_m);
        uu[j] = u[j];
        K_dense.col(j) = (grad_p - grad_m) / (T(2) * h);
    }

    // Symmetrize
    K_dense = T(0.5) * (K_dense + K_dense.transpose());

    K_out = K_dense.sparseView(T(0), T(1e-20));
    K_out.makeCompressed();
}

// ---------------------------------------------------------------------------
// assembleHessian: exact analytic sparse Newton tangent       [H-4]
//
// Obtained by differentiating the *discrete* gradient of gradObj_into
// w.r.t. alpha_B^m. Per quad point, with X=(A,d), Y=(B,m), N=values,
// g=parametric gradients of the sigma basis, and
//   Q = C_eps^{-1}, A_d = Md^T Jc + Jc^T Md, Md = D_d Js, b_d = Jc^T Jg_d,
//   dC_X  = N_A A_d + b_d gA^T + gA b_d^T,
//   dtr_X = -tr(Q dC_X Q),    dphi_X = phi gamma ddet_X:
//
// no-monitor:  K_XY = phi d2tr_XY + dphi_Y dtr_X + dphi_X dtr_Y + trQ d2phi_XY
// monitor:     K_XY = m2 K^core_XY + dm2_X base_Y + dm2_Y base_X
//                     + d2m2_XY (trQ phi),   base_X = phi dtr_X + trQ dphi_X
//
// d2tr_XY = 2 tr(Q dC_Y Q dC_X Q) - tr(Q d2C_XY Q)             [H-4.1a]
// d2phi_XY = phi (gamma^2+gamma') ddet_X ddet_Y + phi gamma d2det_XY [H-4.1b]
//   gamma' = -p/det^2 - 2 (chi'' chi - chi'^2)/chi^2,
//   chi''  = eps^2 / (2 (eps^2+det^2)^{3/2})
// d2C_XY  = N_A dA_d_Y + db_d_Y gA^T + gA db_d_Y^T  with
//   db_d_Y = N_B c_dm + Cg(d,m) gB,   c_dm = Mm^T Jg_d + Jc^T D_dm,
//   dA_d_Y = N_B H_dm + (rank-1 terms in gB),  H_dm = Js^T E_dm Js,
//   E_dm = T_dm^T Jg + Jg^T T_dm + D_d^T D_m + D_m^T D_d   (T = third derivs of S)
// d2det (planar, dd==2; det quadratic in entries => exact bilinear):
//   d2det_XY = bd(dJc_X, dJc_Y) + tr(adj(Jc) d2Jc_XY),  bd(A,B)=tr(adj(A)B)
//   with bd(x a^T, y b^T) = cross(x,y) cross(a,b)
// d2det (surface): d2det_XY = eps_{dm} cross(gA,gB)
//
// Geometric/curvature contributions (third derivatives of S) appear in
// H_dm, tr(adj T_dm Js) and adj*D_dm — all vanish for affine S    [H-4.2].
// Monitor blocks [H-4.3]: see inline comments (ValueBased: Hess f;
// GradientBased: dms(d,m) = d^2(eta^2)/dxi_d dxi_m, needs third derivs of f).
// ---------------------------------------------------------------------------

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::assembleHessian(
    const gsVector<T> & u, gsSparseMatrix<T> & K_out) const
{
    const index_t dd = m_comp->domainDim();
    if (dd != 2)
    {
        // Analytic Hessian implemented for dd==2; FD fallback otherwise.
        _hessianFD(u, K_out);
        return;
    }

    _setControls(u);
    m_optMesh.options().setReal("Penalty",   m_options.getReal("Penalty"));
    m_optMesh.options().setReal("Smoothing", m_options.getReal("Smoothing"));

    const index_t td = m_geom->targetDim();
    const T penalty = m_options.getReal("Penalty");
    const T theta   = m_options.getReal("Smoothing");
    const bool hasMonitor = (m_fun != nullptr);
    const bool planar = (td == dd);
    // q_exp: barrier exponent on chi, see \ref adaptparam_newton.
    const T q_exp = hasMonitor ? T(1) : T(2);

    const gsDofMapper & mapper = m_comp->mapper();
    const gsBasis<T> & sigmaBasis = m_comp->domain().basis();
    const index_t S2 = dd * (dd + 1) / 2;       // = 3

    _buildPattern();
    m_fpattern.assignZero();

    gsOptionList quadOptions;
    quadOptions.addReal("quA","",1.0);
    quadOptions.addInt ("quB","",1);

    auto dom = m_mb.domain();

    // Trigger the geometry deriv3 probe serially (sets m_geomDeriv3Mode).
    {
        gsMatrix<T> c0 = dom->beginAll().centerPoint();
        gsFuncData<T> sd; sd.flags = NEED_VALUE;
        m_comp->compute(c0, sd);
        gsMatrix<T> dummy;
        _geomDeriv3_into(sd.values[0], dummy);
    }

    // DETERMINISTIC reduction (same invariant as the residual/energy sums
    // above): each thread accumulates the matrix entries it touches into
    // its own value buffer, merged in THREAD-ID order after the parallel
    // region, so K_out is bit-reproducible for a fixed thread count.
    const index_t nnzK = m_fpattern.nonZeros();
#   ifdef _OPENMP
    const int nThreads = omp_get_max_threads();
#   else
    const int nThreads = 1;
#   endif
    std::vector<const gsVector<T> *> partialK(nThreads, nullptr);
    std::vector<char> missFlags(nThreads, 0);

#   pragma omp parallel
    {
        // Persistent per-thread scratch (see HessianScratch in the header).
        // thPartialK accumulates this call's contribution and must be
        // re-zeroed here even though its allocation persists across calls.
        HessianScratch & sc = m_hessianScratch.mine();
        if (!sc.ready) sc.init(dd, td);
        if (sc.thPartialK.size() != nnzK) sc.thPartialK.resize(nnzK);
        sc.thPartialK.setZero();
        sc.patternMiss = false;
        gsVector<T> & thPartialK = sc.thPartialK;

        gsFuncData<T> & geomData = sc.geomData;
        gsFuncData<T> & funData  = sc.funData;
        gsFuncData<T> & sigmaData = sc.sigmaData;
        gsFuncData<T> & sigmaBasisData = sc.sigmaBasisData;
        geomData.flags = NEED_VALUE | NEED_DERIV | NEED_DERIV2;
        funData.flags  = NEED_VALUE | NEED_DERIV | NEED_DERIV2;
        sigmaData.flags = NEED_VALUE | NEED_DERIV;
        // SAME_ELEMENT: all quad points of an element share the active set,
        // so actives are computed once per element (cf. gsExprAssembler).
        sigmaBasisData.flags = NEED_ACTIVE | NEED_VALUE | NEED_DERIV | SAME_ELEMENT;

        gsMatrix<T> & uvPoints = sc.uvPoints;
        gsVector<T> & tmpWeights = sc.tmpWeights;
        gsMatrix<T> & Jsigma_flat = sc.Jsigma_flat, & Jgeom_flat = sc.Jgeom_flat;
        gsMatrix<T> & deriv2_geom = sc.deriv2_geom, & deriv3_geom = sc.deriv3_geom;
        gsMatrix<T> & monVals = sc.monVals, & monDerivs = sc.monDerivs,
                    & monDeriv2 = sc.monDeriv2, & monDeriv3 = sc.monDeriv3;
        gsMatrix<index_t> & actives = sc.actives;
        gsMatrix<T> & basisVals = sc.basisVals, & basisDerivs = sc.basisDerivs;

        // common per-point
        gsMatrix<T> & Js = sc.Js, & Jg = sc.Jg, & Jc = sc.Jc, & C = sc.C, & Q = sc.Q;
        gsMatrix<T> & Cg = sc.Cg, & Cg_inv = sc.Cg_inv;
        // per-direction d
        std::vector<gsMatrix<T> > & D = sc.D, & Md = sc.Md, & A = sc.A, & QAQ = sc.QAQ;
        std::vector<gsVector<T> > & b = sc.b, & Qb = sc.Qb, & Q2b = sc.Q2b;
        gsVector<T> & trCAC = sc.trCAC, & ad = sc.ad;
        gsMatrix<T> & AdjJg = sc.AdjJg;
        // per-(d,m)
        std::vector<std::vector<gsMatrix<T> > > & T3 = sc.T3, & E3 = sc.E3;
        std::vector<std::vector<gsVector<T> > > & Ddm = sc.Ddm, & cdm = sc.cdm, & Q2c = sc.Q2c,
                                                 & edm = sc.edm, & amg = sc.amg, & adjD = sc.adjD;
        gsMatrix<T> & M2 = sc.M2, & trQHQ = sc.trQHQ, & S1 = sc.S1, & Ss2 = sc.Ss2,
                    & cJg = sc.cJg, & tT = sc.tT, & bdMM = sc.bdMM;
        gsMatrix<T> & adj = sc.adj, & adjMd = sc.adjMd, & tmpM = sc.tmpM,
                    & Htmp = sc.Htmp, & QQ = sc.QQ;
        gsVector<T> & gA = sc.gA, & gB = sc.gB, & qA = sc.qA, & qB = sc.qB,
                    & v = sc.v, & v_m = sc.v_m;
        gsMatrix<T> & Hess_f = sc.Hess_f, & gradMon = sc.gradMon,
                    & grad_xi_f = sc.grad_xi_f, & grad_x_f = sc.grad_x_f;
        gsVector<T> & ms = sc.ms;
        gsMatrix<T> & dms = sc.dms, & hEta = sc.hEta;
        gsMatrix<T> & localK = sc.localK;
        // GradientBased scratch (reused per quad point)
        std::vector<gsMatrix<T> > & E1 = sc.E1;
        std::vector<gsVector<T> > & gfm = sc.gfm;
        gsVector<T> & dgf = sc.dgf;
        gsMatrix<T> & dHxf = sc.dHxf;
        // surface+monitor: gS product-rule scratch.  Initialised to 1 (the
        // planar no-op value) to match assembleResidual's `T gS = T(1)`: it is
        // only ever READ under (hasMonitor && !planar), which is exactly when
        // it is assigned, but gcc cannot prove those two flags invariant
        // between the assignment and the use and warns -Wmaybe-uninitialized.
        T gS = T(1);
        gsVector<T> & trCginvEd_vec = sc.trCginvEd_vec;
        gsMatrix<T> & trCginvEdm_mat = sc.trCginvEdm_mat, & trCginvEmCginvEd_mat = sc.trCginvEmCginvEd_mat;

        for (auto & elem : dom->allElements())
        {
            if (sc.QuPatch != elem.patch())
            {
                sc.QuPatch = elem.patch();
                sc.QuRule = gsQuadrature::get(m_ib->basis(sc.QuPatch), quadOptions);
            }

            sc.QuRule.mapTo(elem.lowerCorner(), elem.upperCorner(), uvPoints, tmpWeights);

            m_comp->compute(uvPoints, sigmaData);
            Jsigma_flat = sigmaData.values[1];

            sigmaBasis.compute(uvPoints, sigmaBasisData);
            actives     = sigmaBasisData.actives;
            basisVals   = sigmaBasisData.values[0];
            basisDerivs = sigmaBasisData.values[1];
            const index_t nAct = actives.rows();

            m_geom->compute(sigmaData.values[0], geomData);
            Jgeom_flat  = geomData.values[1];
            deriv2_geom = geomData.values[2];
            _geomDeriv3_into(sigmaData.values[0], deriv3_geom);

            if (hasMonitor)
            {
                const gsMatrix<T> & monPts =
                    (m_parametric ? sigmaData.values[0] : geomData.values[0]);
                m_fun->compute(monPts, funData);
                monVals   = funData.values[0];
                monDerivs = funData.values[1];
                monDeriv2 = funData.values[2];
                if (MODE == GradientBased)
                    _funDeriv3_into(monPts, monDeriv3);
            }

            const index_t ndof = dd * nAct;
            localK.setZero(ndof, ndof);

            const index_t nPts = uvPoints.cols();
            for (index_t p = 0; p != nPts; ++p)
            {
                const T w_q = tmpWeights[p];

                Js.noalias() = Jsigma_flat.col(p).reshaped(dd, dd).transpose();
                Jg.noalias() = Jgeom_flat.col(p).reshaped(dd, td).transpose();
                Jc.noalias() = Jg * Js;

                C.noalias() = Jc.transpose() * Jc;
                // NO Tikhonov shift on C (matches gsOptMesh's unregularised
                // energy); C is inverted unregularised, so a folded iterate
                // returns NaN here. Rationale and the fold hazard: \ref adaptparam_newton.
                Q.noalias() = C.inverse();
                QQ.noalias() = Q * Q;
                const T trQ = Q.trace();
                Cg.noalias() = Jg.transpose() * Jg;

                // ---- barrier scalars
                T det, chi, chip, chipp, phi, gamma, gammap;
                if (planar) det = Jc.determinant();
                else        det = Js.determinant();
                {
                    const T sq  = math::sqrt(penalty*penalty + det*det);
                    chi   = T(0.5) * (det + sq);
                    chip  = T(0.5) * (T(1) + det / sq);
                    chipp = T(0.5) * penalty * penalty / (sq * sq * sq);
                    // q_exp is always 1 or 2 (see its definition above); avoid a
                    // libm pow() call per quadrature point -- hasMonitor is
                    // loop-invariant, so this branch is predictable and free.
                    phi = hasMonitor ? (det * det / chi) : (det * det / (chi * chi));
                    const T inv_det = (det != T(0)) ? T(1)/det : T(0);
                    gamma  = T(2) * inv_det - q_exp * chip / chi;
                    // gamma' = d gamma / d det  (sign: see [H-4.1b])
                    gammap = -T(2) * inv_det * inv_det
                             - q_exp * (chipp * chi - chip * chip) / (chi * chi);
                }

                if (planar)
                {
                    adj(0,0) =  Jc(1,1); adj(0,1) = -Jc(0,1);
                    adj(1,0) = -Jc(1,0); adj(1,1) =  Jc(0,0);
                    AdjJg.noalias() = adj * Jg;
                }

                // ---- per-direction quantities
                for (index_t d = 0; d != dd; ++d)
                {
                    for (index_t a = 0; a != td; ++a)
                        for (index_t j = 0; j != dd; ++j)
                        {
                            const index_t lo = math::min(d, j), hi = math::max(d, j);
                            const index_t hidx = (lo == hi) ? lo : dd + lo*(2*dd-lo-3)/2 + hi - 1;
                            D[d](a, j) = deriv2_geom(a * S2 + hidx, p);
                        }
                    Md[d].noalias() = D[d] * Js;
                    A[d].noalias()  = Md[d].transpose() * Jc + Jc.transpose() * Md[d];
                    b[d].noalias()  = Jc.transpose() * Jg.col(d);
                    QAQ[d].noalias() = Q * A[d] * Q;
                    trCAC(d) = QAQ[d].trace();   // tr(Q A_d Q)
                    Qb[d].noalias()  = Q * b[d];
                    Q2b[d].noalias() = Q * Qb[d];

                    if (planar)
                        ad(d) = (adj * Md[d]).trace();
                }

                // ---- per-(d,m) quantities
                for (index_t d = 0; d != dd; ++d)
                    for (index_t m = 0; m != dd; ++m)
                    {
                        // third derivatives of S: T_dm(a,j) = S_{a,jdm}
                        for (index_t a = 0; a != td; ++a)
                            for (index_t j = 0; j != dd; ++j)
                            {
                                const index_t lo = math::min(d, j), hi = math::max(d, j);
                                const index_t hidx = (lo == hi) ? lo : dd + lo*(2*dd-lo-3)/2 + hi - 1;
                                T3[d][m](a, j) = deriv3_geom((m * td + a) * S2 + hidx, p);
                            }

                        Ddm[d][m] = D[d].col(m);

                        // E_dm = T^T Jg + Jg^T T + D_d^T D_m + D_m^T D_d
                        E3[d][m].noalias() = T3[d][m].transpose() * Jg;
                        E3[d][m] += E3[d][m].transpose().eval();
                        E3[d][m].noalias() += D[d].transpose() * D[m];
                        E3[d][m].noalias() += D[m].transpose() * D[d];

                        // H_dm = Js^T E_dm Js;  trQHQ = tr(Q H Q)
                        Htmp.noalias() = Js.transpose() * E3[d][m] * Js;
                        trQHQ(d,m) = (Q * Htmp * Q).trace();

                        // c_dm = Mm^T Jg_d + Jc^T D_dm
                        cdm[d][m].noalias() = Md[m].transpose() * Jg.col(d);
                        cdm[d][m].noalias() += Jc.transpose() * Ddm[d][m];
                        Q2c[d][m].noalias() = QQ * cdm[d][m];

                        // M2 = tr(Q A_m Q A_d Q)
                        M2(d,m) = (QAQ[m] * A[d] * Q).trace();
                        // mixed A/rank-1 terms of tr(Q dC_Y Q dC_X Q):
                        //   tr(Q b_m gB^T Q A_d Q) = gB . (QA_dQ) Q b_m
                        //   tr(Q gB b_m^T Q A_d Q) = gB . Q (QA_dQ) b_m
                        // (A_d and Q do not commute => two distinct vectors)
                        edm[d][m].noalias() = QAQ[d] * Qb[m];
                        edm[d][m].noalias() += Q * (QAQ[d] * b[m]);

                        S1(d,m)  = b[d].dot(Qb[m]);
                        Ss2(d,m) = Qb[d].dot(Qb[m]);

                        if (planar)
                        {
                            // adj(Md), bd(Md,Mm) = tr(adj(Md) Mm)
                            adjMd(0,0) =  Md[d](1,1); adjMd(0,1) = -Md[d](0,1);
                            adjMd(1,0) = -Md[d](1,0); adjMd(1,1) =  Md[d](0,0);
                            bdMM(d,m) = (adjMd * Md[m]).trace();
                            amg[d][m].noalias() = adjMd * Jg.col(m);
                            cJg(d,m) = Jg(0,d)*Jg(1,m) - Jg(1,d)*Jg(0,m);
                            tT(d,m)  = (adj * T3[d][m] * Js).trace();
                            adjD[d][m].noalias() = adj * Ddm[d][m];
                        }
                    }

                // ---- monitor scalars
                T m2 = T(1);
                gsVector<T> & dm2coef = sc.dm2coef;   // dm2_X = dm2coef(d) N_A
                gsMatrix<T> & d2m2coef = sc.d2m2coef; // d2m2_XY = d2m2coef(d,m) N_A N_B
                dm2coef.setZero(); d2m2coef.setZero();

                if (hasMonitor)
                {
                    if (MODE == ValueBased)
                    {
                        const T eta = monVals(0, p);
                        //  omega = 1/sqrt(1+theta*f) => m2 = omega^2 = 1/(1+theta*f)
                        const T denom = T(1) + theta * eta;
                        m2 = T(1) / denom;
                        const T dm2_deta  = -theta / (denom * denom);            // dm2/df
                        const T d2m2_deta2 = T(2) * theta * theta
                                           / (denom * denom * denom);           // d2m2/df2

                        if (m_parametric)
                        {
                            gradMon.noalias() = monDerivs.col(p);
                            // hEta(d,m) = Hess_xi f
                            for (index_t i = 0; i != dd; ++i)
                                for (index_t j = i; j != dd; ++j)
                                {
                                    const index_t hidx = (i==j) ? i : dd + i*(2*dd-i-3)/2+j-1;
                                    hEta(i,j) = monDeriv2(hidx, p);
                                    hEta(j,i) = hEta(i,j);
                                }
                        }
                        else
                        {
                            grad_x_f.noalias() = monDerivs.col(p).reshaped(td, 1);
                            gradMon.noalias() = Jg.transpose() * grad_x_f;
                            Hess_f.resize(td, td);
                            for (index_t i = 0; i != td; ++i)
                                for (index_t j = i; j != td; ++j)
                                {
                                    const index_t hidx = (i==j) ? i : td + i*(2*td-i-3)/2+j-1;
                                    Hess_f(i,j) = monDeriv2(hidx, p);
                                    Hess_f(j,i) = Hess_f(i,j);
                                }
                            // hEta(d,m) = D_dm . grad_x f + Jg_d^T Hxf Jg_m
                            for (index_t d = 0; d != dd; ++d)
                                for (index_t m = 0; m != dd; ++m)
                                    hEta(d,m) = Ddm[d][m].dot(grad_x_f.col(0))
                                              + (Jg.col(d).transpose() * Hess_f * Jg.col(m))(0,0);
                        }

                        for (index_t d = 0; d != dd; ++d)
                        {
                            dm2coef(d) = dm2_deta * gradMon(d,0);
                            for (index_t m = 0; m != dd; ++m)
                                d2m2coef(d,m) = d2m2_deta2 * gradMon(d,0) * gradMon(m,0)
                                              + dm2_deta * hEta(d,m);
                        }
                    }
                    else // GradientBased
                    {
                        Cg_inv.noalias() = Cg.inverse();

                        const index_t df = m_parametric ? dd : td;
                        const index_t Sf = df * (df + 1) / 2;
                        Hess_f.resize(df, df);
                        for (index_t i = 0; i != df; ++i)
                            for (index_t j = i; j != df; ++j)
                            {
                                const index_t hidx = (i==j) ? i : df + i*(2*df-i-3)/2+j-1;
                                Hess_f(i,j) = monDeriv2(hidx, p);
                                Hess_f(j,i) = Hess_f(i,j);
                            }

                        if (m_parametric)
                            grad_xi_f.noalias() = monDerivs.col(p);
                        else
                        {
                            grad_x_f.noalias()  = monDerivs.col(p).reshaped(td, 1);
                            grad_xi_f.noalias() = Jg.transpose() * grad_x_f;
                        }

                        const T eta2 = (grad_xi_f.transpose() * Cg_inv * grad_xi_f)(0,0);
                        const T denom = T(1) + theta * eta2;
                        m2 = T(1) / denom;
                        const T dm2_dq  = -theta / (denom * denom);
                        const T d2m2_dq2 = T(2) * theta * theta / (denom * denom * denom);

                        v.noalias() = Cg_inv * grad_xi_f;

                        // E_d = D_d^T Jg + Jg^T D_d (dCg per unit N in dir d)
                        // (E1, gfm: per-thread scratch reused each quad point)
                        for (index_t d = 0; d != dd; ++d)
                        {
                            E1[d].noalias() = D[d].transpose() * Jg;
                            E1[d] += E1[d].transpose().eval();

                            if (m_parametric)
                                gfm[d] = Hess_f.col(d);
                            else
                            {
                                gfm[d].noalias() = D[d].transpose() * grad_x_f;
                                gfm[d].noalias() += Jg.transpose() * (Hess_f * Jg.col(d));
                            }

                            // ms(d) = d(eta^2)/dxi_d = -v^T E_d v + 2 v^T gfm_d
                            ms(d) = -(v.transpose() * E1[d] * v)(0,0)
                                    + T(2) * v.dot(gfm[d]);
                        }

                        // dms(d,m) = d^2(eta^2)/dxi_d dxi_m
                        //  = -2 v_m^T E_d v - v^T E_dm v + 2 v_m^T gfm_d + 2 v^T d_m(gfm_d)
                        //  v_m = Cg^{-1}(gfm_m - E_m v)
                        for (index_t m = 0; m != dd; ++m)
                        {
                            v_m.noalias() = Cg_inv * (gfm[m] - E1[m] * v);
                            for (index_t d = 0; d != dd; ++d)
                            {
                                // d_m(gfm_d): third derivatives of f (dgf scratch)
                                if (m_parametric)
                                {
                                    // tf(j) = f_{,j d m} from monDeriv3 layout
                                    for (index_t j = 0; j != dd; ++j)
                                    {
                                        const index_t lo = math::min(d, j), hi = math::max(d, j);
                                        const index_t hidx = (lo == hi) ? lo : df + lo*(2*df-lo-3)/2 + hi - 1;
                                        dgf(j) = monDeriv3(m * Sf + hidx, p);
                                    }
                                }
                                else
                                {
                                    // d_m [ D_d^T gxf + Jg^T Hxf Jg_d ]
                                    // = T_dm^T gxf + D_d^T Hxf Jg_m
                                    // + D_m^T Hxf Jg_d + Jg^T (d_m Hxf) Jg_d + Jg^T Hxf D_dm
                                    // d_m Hxf = sum_c Txf(:,:,c) Jg(c,m) (dHxf scratch)
                                    dHxf.setZero();
                                    for (index_t c = 0; c != td; ++c)
                                    {
                                        const T jgcm = Jg(c, m);
                                        if (jgcm == T(0)) continue;
                                        for (index_t i = 0; i != td; ++i)
                                            for (index_t j = i; j != td; ++j)
                                            {
                                                const index_t hidx = (i==j) ? i : td + i*(2*td-i-3)/2+j-1;
                                                const T t3 = monDeriv3(c * Sf + hidx, p);
                                                dHxf(i,j) += t3 * jgcm;
                                                if (i != j) dHxf(j,i) += t3 * jgcm;
                                            }
                                    }
                                    dgf.noalias() = T3[d][m].transpose() * grad_x_f;
                                    dgf.noalias() += D[d].transpose() * (Hess_f * Jg.col(m));
                                    dgf.noalias() += D[m].transpose() * (Hess_f * Jg.col(d));
                                    dgf.noalias() += Jg.transpose() * (dHxf * Jg.col(d));
                                    dgf.noalias() += Jg.transpose() * (Hess_f * Ddm[d][m]);
                                }

                                dms(d,m) = -T(2) * (v_m.transpose() * E1[d] * v)(0,0)
                                           - (v.transpose() * E3[d][m] * v)(0,0)
                                           + T(2) * v_m.dot(gfm[d])
                                           + T(2) * v.dot(dgf);
                            }
                        }

                        for (index_t d = 0; d != dd; ++d)
                        {
                            dm2coef(d) = dm2_dq * ms(d);
                            for (index_t m = 0; m != dd; ++m)
                                d2m2coef(d,m) = d2m2_dq2 * ms(d) * ms(m)
                                              + dm2_dq * dms(d,m);
                        }
                    }
                }

                // gS product-rule precomputation for surface+monitor analytic Hessian
                if (!planar && hasMonitor)
                {
                    gS = math::sqrt(Cg.determinant());

                    if (MODE == ValueBased)
                    {
                        // Cg_inv and E1[d] not computed in the ValueBased monitor block above
                        Cg_inv.noalias() = Cg.inverse();
                        for (index_t d = 0; d != dd; ++d)
                        {
                            E1[d].noalias() = D[d].transpose() * Jg;
                            E1[d] += E1[d].transpose().eval();
                        }
                    }
                    // GradientBased: Cg_inv and E1[d] already set in the GradientBased
                    // monitor block above.

                    for (index_t d = 0; d != dd; ++d)
                        trCginvEd_vec(d) = (Cg_inv * E1[d]).trace();

                    for (index_t d = 0; d != dd; ++d)
                        for (index_t m = 0; m != dd; ++m)
                        {
                            trCginvEdm_mat(d,m)       = (Cg_inv * E3[d][m]).trace();
                            trCginvEmCginvEd_mat(d,m) = (Cg_inv * E1[m] * Cg_inv * E1[d]).trace();
                        }
                }

                const T g2gp = gamma * gamma + gammap;
                const T trQphi = trQ * phi;

                // ---- local kernel: double loop over (active, component)
                for (index_t lA = 0; lA != nAct; ++lA)
                {
                    const T NA = basisVals(lA, p);
                    for (index_t j = 0; j != dd; ++j)
                        gA(j) = basisDerivs(lA * dd + j, p);
                    qA.noalias() = Q * gA;

                    for (index_t lB = 0; lB != nAct; ++lB)
                    {
                        const T NB = basisVals(lB, p);
                        for (index_t j = 0; j != dd; ++j)
                            gB(j) = basisDerivs(lB * dd + j, p);
                        qB.noalias() = Q * gB;

                        const T qAqB = qA.dot(qB);
                        const T crossAB = gA(0)*gB(1) - gA(1)*gB(0);

                        for (index_t d = 0; d != dd; ++d)
                        {
                            const T dtrX = -(NA * trCAC(d) + T(2) * Qb[d].dot(qA));
                            T ddetX;
                            if (planar) ddetX = NA * ad(d) + gA.dot(AdjJg.col(d));
                            else        ddetX = (d == 0)
                                            ?  Js(1,1)*gA(0) - Js(1,0)*gA(1)
                                            : -Js(0,1)*gA(0) + Js(0,0)*gA(1);
                            const T dphiX = phi * gamma * ddetX;
                            const T baseX = phi * dtrX + trQ * dphiX;

                            for (index_t m = 0; m != dd; ++m)
                            {
                                const T dtrY = -(NB * trCAC(m) + T(2) * Qb[m].dot(qB));
                                T ddetY;
                                if (planar) ddetY = NB * ad(m) + gB.dot(AdjJg.col(m));
                                else        ddetY = (m == 0)
                                                ?  Js(1,1)*gB(0) - Js(1,0)*gB(1)
                                                : -Js(0,1)*gB(0) + Js(0,0)*gB(1);
                                const T dphiY = phi * gamma * ddetY;
                                const T baseY = phi * dtrY + trQ * dphiY;

                                // tr(Q dC_Y Q dC_X Q): 9 expanded terms
                                const T trQdCQdCQ =
                                      NA * NB * M2(d,m)
                                    + NA * gB.dot(edm[d][m])
                                    + NB * gA.dot(edm[m][d])
                                    + gB.dot(Qb[d]) * qA.dot(Qb[m])
                                    + qA.dot(gB) * Ss2(d,m)
                                    + S1(d,m) * qAqB
                                    + gA.dot(Qb[m]) * qB.dot(Qb[d]);

                                // tr(Q d2C_XY Q)
                                const T trQd2CQ =
                                      NA * NB * trQHQ(d,m)
                                    + T(2) * NA * gB.dot(Q2c[m][d])
                                    + T(2) * NB * gA.dot(Q2c[d][m])
                                    + T(2) * Cg(d,m) * qAqB;

                                const T d2tr = T(2) * trQdCQdCQ - trQd2CQ;

                                // d2det
                                T d2det;
                                if (planar)
                                {
                                    d2det = NA * NB * (bdMM(d,m) + tT(d,m))
                                          + NA * (gB.dot(amg[d][m]) + gB.dot(adjD[d][m]))
                                          + NB * (gA.dot(amg[m][d]) + gA.dot(adjD[d][m]))
                                          + cJg(d,m) * crossAB;
                                }
                                else
                                    d2det = (d == m) ? T(0)
                                          : ((d == 0) ? crossAB : -crossAB);

                                const T d2phi = phi * g2gp * ddetX * ddetY
                                              + phi * gamma * d2det;

                                T Kval = phi * d2tr + dphiY * dtrX + dphiX * dtrY
                                       + trQ * d2phi;

                                if (hasMonitor)
                                {
                                    // Integrand: A = omega*(tr*phi), omega=sqrt(m2).
                                    // d²A = omega*d²(tr*phi) + domega_d*NA*d(tr*phi)_m
                                    //     + domega_m*NB*d(tr*phi)_d + d2omega*NA*NB*(tr*phi)
                                    // (Kval = d²(tr*phi), baseX = d(tr*phi)_d, trQphi = tr*phi.)
                                    const T omega    = math::sqrt(m2);
                                    const T inv2w    = T(1) / (T(2) * omega);
                                    const T domega_d = dm2coef(d) * inv2w;
                                    const T domega_m = dm2coef(m) * inv2w;
                                    const T d2omega  = d2m2coef(d,m) * inv2w
                                        - dm2coef(d) * dm2coef(m) / (T(4) * omega * omega * omega);

                                    if (planar)
                                    {
                                        Kval = omega * Kval
                                             + (domega_d * NA) * baseY
                                             + (domega_m * NB) * baseX
                                             + (d2omega * NA * NB) * trQphi;
                                    }
                                    else
                                    {
                                        // Surface: integrand = A*gS; apply gS product rule.
                                        // d²A (assembled above):
                                        const T d2A = omega * Kval
                                                    + (domega_d * NA) * baseY
                                                    + (domega_m * NB) * baseX
                                                    + (d2omega * NA * NB) * trQphi;
                                        // dA/dα_{A,d} and dA/dβ_{B,m}:
                                        const T dA_alpha = omega * baseX + domega_d * NA * trQphi;
                                        const T dA_beta  = omega * baseY + domega_m * NB * trQphi;
                                        // A itself:
                                        const T A_val    = omega * trQphi;
                                        // d²gS curvature coefficient (NA*NB factored out):
                                        const T d2gS_coef = gS * (
                                            T(0.25) * trCginvEd_vec(d) * trCginvEd_vec(m)
                                            + T(0.5) * (trCginvEdm_mat(d,m) - trCginvEmCginvEd_mat(d,m))
                                        );
                                        // Product rule: d²(A*gS) = d²A*gS + dA_α*dgS_β + dA_β*dgS_α + A*d²gS
                                        Kval = d2A * gS
                                             + dA_alpha * (T(0.5) * gS * trCginvEd_vec(m) * NB)
                                             + dA_beta  * (T(0.5) * gS * trCginvEd_vec(d) * NA)
                                             + A_val * NA * NB * d2gS_coef;
                                    }
                                }

                                localK(d * nAct + lA, m * nAct + lB) += w_q * Kval;
                            }
                        }
                    }
                }
            } // quad points

            // ---- scatter element matrix into this thread's own value
            // buffer (accumulate, no cross-thread write) ----
            for (index_t lA = 0; lA != nAct; ++lA)
            {
                const index_t kA = actives(lA, 0);
                for (index_t d = 0; d != dd; ++d)
                {
                    if (!mapper.is_free(kA, 0, d)) continue;
                    const index_t ii = mapper.index(kA, 0, d);
                    for (index_t lB = 0; lB != nAct; ++lB)
                    {
                        const index_t kB = actives(lB, 0);
                        for (index_t m = 0; m != dd; ++m)
                        {
                            if (!mapper.is_free(kB, 0, m)) continue;
                            const index_t jj = mapper.index(kB, 0, m);
                            const T val = localK(d * nAct + lA, m * nAct + lB);
                            if (val == T(0)) continue;
                            const index_t pos = _fpatternPos(ii, jj);
                            if (pos < 0) { sc.patternMiss = true; continue; }
                            thPartialK[pos] += val;
                        }
                    }
                }
            }
        } // elements

        // Pointer (not copy): thPartialK aliases the persistent per-thread
        // scratch buffer; _mergePartialK reads it directly.
        const int tid =
#       ifdef _OPENMP
            omp_get_thread_num();
#       else
            0;
#       endif
        partialK[tid] = &thPartialK;
        missFlags[tid] = sc.patternMiss ? 1 : 0;
    } // omp parallel

    bool anyMiss = false;
    for (char m : missFlags) anyMiss = anyMiss || m;
    GISMO_ENSURE(!anyMiss, "gsAdaptiveParametrizationNewton: tangent assembly "
        "hit a (row,col) pair outside the precomputed sparsity pattern; the "
        "pattern is built from element centre points and the quadrature rule "
        "must not produce nodes outside that support graph.");

    _mergePartialK(partialK, K_out);
}

// ---------------------------------------------------------------------------
// assemblePicard: SPD frozen-coefficient surrogate B            [H-5]
//
//   B_{(A d)(B m)} = ∫ coeff * Cg(d,m) * (∇N_A)^T Bc ∇N_B  dΩ̂
//
//   no-monitor (q=2):  Bc = Q^2 (exact frozen EL operator, SPD)    [H-5.2]
//   with-monitor (q=1): Bc = (2 w / g) Q   (SPD surrogate;
//     the literal frozen q=1 operator is indefinite)               [H-5.3]
//
// w is the monitor weight m^2 frozen at the current iterate (both
// ValueBased and GradientBased computed exactly).
// Same two-phase OpenMP assembly as assembleHessian.
// ---------------------------------------------------------------------------

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::assemblePicard(
    const gsVector<T> & u, gsSparseMatrix<T> & B_out) const
{
    _setControls(u);
    m_optMesh.options().setReal("Penalty",   m_options.getReal("Penalty"));
    m_optMesh.options().setReal("Smoothing", m_options.getReal("Smoothing"));

    const index_t dd = m_comp->domainDim();
    const index_t td = m_geom->targetDim();
    GISMO_ASSERT(dd <= td, "domainDim must be <= targetDim");

    const T penalty = m_options.getReal("Penalty");
    // No chi/phi barrier is assembled in this routine (see the NO-Tikhonov-shift
    // comment below), so `penalty` has no remaining numerical use here; it is
    // still read for parity with the other three assembly routines' option handling.
    GISMO_UNUSED(penalty);
    const T theta   = m_options.getReal("Smoothing");
    const bool hasMonitor = (m_fun != nullptr);

    const gsDofMapper & mapper = m_comp->mapper();
    const gsBasis<T> & sigmaBasis = m_comp->domain().basis();

    _buildPattern();
    m_fpattern.assignZero();

    gsOptionList quadOptions;
    quadOptions.addReal("quA","",1.0);
    quadOptions.addInt ("quB","",1);

    auto dom = m_mb.domain();

    // DETERMINISTIC reduction (same invariant as assembleHessian): each
    // thread accumulates the matrix entries it touches into its own value
    // buffer, merged in THREAD-ID order after the parallel region, so
    // B_out is bit-reproducible for a fixed thread count.
    const index_t nnzB = m_fpattern.nonZeros();
#   ifdef _OPENMP
    const int nThreads = omp_get_max_threads();
#   else
    const int nThreads = 1;
#   endif
    std::vector<const gsVector<T> *> partialK(nThreads, nullptr);
    std::vector<char> missFlags(nThreads, 0);

#   pragma omp parallel
    {
        // Persistent per-thread scratch (see PicardScratch in the header).
        // thPartialK accumulates this call's contribution and must be
        // re-zeroed here even though its allocation persists across calls.
        PicardScratch & sc = m_picardScratch.mine();
        if (!sc.ready) sc.init(dd, td);
        if (sc.thPartialK.size() != nnzB) sc.thPartialK.resize(nnzB);
        sc.thPartialK.setZero();
        sc.patternMiss = false;
        gsVector<T> & thPartialK = sc.thPartialK;

        gsFuncData<T> & geomData = sc.geomData;
        gsFuncData<T> & funData  = sc.funData;
        gsFuncData<T> & sigmaData = sc.sigmaData;
        gsFuncData<T> & sigmaBasisData = sc.sigmaBasisData;
        geomData.flags = NEED_VALUE | NEED_DERIV;
        funData.flags  = NEED_VALUE;
        if (hasMonitor && MODE == GradientBased)
            funData.flags |= NEED_DERIV;
        sigmaData.flags = NEED_VALUE | NEED_DERIV;
        // SAME_ELEMENT: all quad points of an element share the active set,
        // so actives are computed once per element (cf. gsExprAssembler).
        sigmaBasisData.flags = NEED_ACTIVE | NEED_VALUE | NEED_DERIV | SAME_ELEMENT;

        gsMatrix<T> & uvPoints = sc.uvPoints;
        gsVector<T> & tmpWeights = sc.tmpWeights;
        gsMatrix<T> & Jsigma_flat = sc.Jsigma_flat, & Jgeom_flat = sc.Jgeom_flat;
        gsMatrix<T> & monVals = sc.monVals, & monDerivs = sc.monDerivs;
        gsMatrix<index_t> & actives = sc.actives;
        gsMatrix<T> & basisVals = sc.basisVals, & basisDerivs = sc.basisDerivs;

        gsMatrix<T> & Js = sc.Js, & Jg = sc.Jg, & Jc = sc.Jc, & C = sc.C, & Q = sc.Q;
        gsMatrix<T> & Cg = sc.Cg, & Cg_inv = sc.Cg_inv, & Bc = sc.Bc;
        gsMatrix<T> & grad_xi_f = sc.grad_xi_f;
        gsVector<T> & gA = sc.gA, & gB = sc.gB, & BcgB = sc.BcgB;
        gsMatrix<T> & localB = sc.localB;

        for (auto & elem : dom->allElements())
        {
            if (sc.QuPatch != elem.patch())
            {
                sc.QuPatch = elem.patch();
                sc.QuRule = gsQuadrature::get(m_ib->basis(sc.QuPatch), quadOptions);
            }

            sc.QuRule.mapTo(elem.lowerCorner(), elem.upperCorner(), uvPoints, tmpWeights);

            m_comp->compute(uvPoints, sigmaData);
            Jsigma_flat = sigmaData.values[1];

            sigmaBasis.compute(uvPoints, sigmaBasisData);
            actives     = sigmaBasisData.actives;
            basisVals   = sigmaBasisData.values[0];
            basisDerivs = sigmaBasisData.values[1];
            const index_t nAct = actives.rows();

            m_geom->compute(sigmaData.values[0], geomData);
            Jgeom_flat = geomData.values[1];

            if (hasMonitor)
            {
                m_fun->compute((m_parametric ? sigmaData.values[0] : geomData.values[0]), funData);
                monVals = funData.values[0];
                if (MODE == GradientBased)
                    monDerivs = funData.values[1];
            }

            const index_t ndof = dd * nAct;
            localB.setZero(ndof, ndof);

            const index_t nPts = uvPoints.cols();
            for (index_t p = 0; p != nPts; ++p)
            {
                const T w_q = tmpWeights[p];

                Js.noalias() = Jsigma_flat.col(p).reshaped(dd, dd).transpose();
                Jg.noalias() = Jgeom_flat.col(p).reshaped(dd, td).transpose();
                Jc.noalias() = Jg * Js;

                C.noalias() = Jc.transpose() * Jc;
                // NO Tikhonov shift on C (matches gsOptMesh's unregularised
                // energy); C is inverted unregularised, so a folded iterate
                // returns NaN here. Rationale and the fold hazard: \ref adaptparam_newton.
                Q.noalias() = C.inverse();
                Cg.noalias() = Jg.transpose() * Jg;

                if (!hasMonitor)
                {
                    // q=2: Bc = Q^2 (exact frozen EL operator)    [H-5.2]
                    Bc.noalias() = Q * Q;
                }
                else
                {
                    // q=1 surrogate: Bc = (2 w / g) Q             [H-5.3]
                    T g_val;
                    if (td == dd) g_val = math::abs(Jc.determinant());
                    else          g_val = math::abs(Js.determinant());
                    const T g_safe = math::max(g_val, T(1e-14));

                    // frozen monitor weight w = omega = sqrt(m2) at the current iterate
                    // m2 = 1/(1+theta*f) or 1/(1+theta*||∇f||^2)
                    T eta2 = T(0);
                    if (MODE == ValueBased)
                    {
                        const T eta = monVals(0, p);
                        eta2 = eta;   // value-based denom is linear in f (Eq.13)
                    }
                    else
                    {
                        Cg_inv.noalias() = Cg.inverse();
                        if (m_parametric)
                            grad_xi_f.noalias() = monDerivs.col(p);
                        else
                            grad_xi_f.noalias() = Jg.transpose() * monDerivs.col(p).reshaped(td, 1);
                        eta2 = (grad_xi_f.transpose() * Cg_inv * grad_xi_f)(0,0);
                    }
                    // omega = sqrt(m2) is the paper's monitor weight (integrand uses omega, not m2)
                    const T m2    = T(1) / (T(1) + theta * eta2);
                    const T w_mon = math::sqrt(m2);
                    Bc = (T(2) * w_mon / g_safe) * Q;
                }

                for (index_t lB = 0; lB != nAct; ++lB)
                {
                    for (index_t j = 0; j != dd; ++j)
                        gB(j) = basisDerivs(lB * dd + j, p);
                    BcgB.noalias() = Bc * gB;

                    for (index_t lA = 0; lA != nAct; ++lA)
                    {
                        for (index_t j = 0; j != dd; ++j)
                            gA(j) = basisDerivs(lA * dd + j, p);
                        const T quad_bc = gA.dot(BcgB);

                        for (index_t d = 0; d != dd; ++d)
                            for (index_t m = 0; m != dd; ++m)
                                localB(d * nAct + lA, m * nAct + lB)
                                    += w_q * Cg(d, m) * quad_bc;
                    }
                }
            }

            // scatter into this thread's own value buffer
            for (index_t lA = 0; lA != nAct; ++lA)
            {
                const index_t kA = actives(lA, 0);
                for (index_t d = 0; d != dd; ++d)
                {
                    if (!mapper.is_free(kA, 0, d)) continue;
                    const index_t ii = mapper.index(kA, 0, d);
                    for (index_t lB = 0; lB != nAct; ++lB)
                    {
                        const index_t kB = actives(lB, 0);
                        for (index_t m = 0; m != dd; ++m)
                        {
                            if (!mapper.is_free(kB, 0, m)) continue;
                            const index_t jj = mapper.index(kB, 0, m);
                            const T val = localB(d * nAct + lA, m * nAct + lB);
                            if (val == T(0)) continue;
                            const index_t pos = _fpatternPos(ii, jj);
                            if (pos < 0) { sc.patternMiss = true; continue; }
                            thPartialK[pos] += val;
                        }
                    }
                }
            }
        }

        // Pointer (not copy): thPartialK aliases the persistent per-thread
        // scratch buffer; _mergePartialK reads it directly.
        const int tid =
#       ifdef _OPENMP
            omp_get_thread_num();
#       else
            0;
#       endif
        partialK[tid] = &thPartialK;
        missFlags[tid] = sc.patternMiss ? 1 : 0;
    } // omp parallel

    bool anyMiss = false;
    for (char m : missFlags) anyMiss = anyMiss || m;
    GISMO_ENSURE(!anyMiss, "gsAdaptiveParametrizationNewton: tangent assembly "
        "hit a (row,col) pair outside the precomputed sparsity pattern; the "
        "pattern is built from element centre points and the quadrature rule "
        "must not produce nodes outside that support graph.");

    _mergePartialK(partialK, B_out);
}

// ---------------------------------------------------------------------------
// Make system matrix positive definite via Levenberg–Marquardt shift
// ---------------------------------------------------------------------------

template<class T, enum MonitorMode MODE>
void gsAdaptiveParametrizationNewton<T,MODE>::_makeDefinite(
    gsSparseMatrix<T> & K, T & tau) const
{
    // Descent is detected in _solveLoop after the solve; here just apply
    // the stored tau shift.
    if (tau > T(0))
    {
        const index_t n = K.rows();
        for (index_t i = 0; i < n; ++i)
            K.coeffRef(i, i) += tau;
        K.makeCompressed();
    }
}

// ---------------------------------------------------------------------------
// Armijo back-tracking line search
// ---------------------------------------------------------------------------

template<class T, enum MonitorMode MODE>
T gsAdaptiveParametrizationNewton<T,MODE>::_armijoStep(
    const gsVector<T> & u0,
    const gsVector<T> & R,
    const gsVector<T> & direction,
    T f0,
    index_t * nBacktracks) const
{
    const T c1  = m_options.getReal("ArmijoC1");
    const T rho = m_options.getReal("ArmijoRho");
    const int maxBT = m_options.getInt ("ArmijoMax");

    // Gradient · direction = R · direction (since ∇E = R)
    const T slope = R.dot(direction);

    T alpha = T(1);
    gsVector<T> u_new(u0.size());
    for (int bt = 0; bt < maxBT; ++bt)
    {
        if (nBacktracks) ++(*nBacktracks);

        u_new.noalias() = u0 + alpha * direction;

        // One fused sweep yields both the trial energy and its min Jacobian, so
        // the fold guard costs nothing beyond the energy evaluation.
        T f_new, minJ;
        _evalObjAndMinJ(u_new, f_new, minJ);

        // Accept iff sufficient decrease AND no fold.
        if (minJ > T(0) && f_new <= f0 + c1 * alpha * slope)
            return alpha;

        alpha *= rho;
    }
    // Line search failed: every trial step above left the controls set to the
    // last REJECTED trial point (_evalObjAndMinJ calls _setControls on each
    // candidate). Restore u0 so the caller and any subsequent query
    // (computeMinJacobian, evalObj, ...) see the last ACCEPTED iterate, not a
    // rejected/folded one.
    _setControls(u0);
    return T(0);   // line search failed: no acceptable step
}

// ---------------------------------------------------------------------------
// Shared solve loop: Newton (useHessian=true) or Picard (useHessian=false)
// ---------------------------------------------------------------------------

template<class T, enum MonitorMode MODE>
index_t gsAdaptiveParametrizationNewton<T,MODE>::_solveLoop(bool useHessian)
{
    const T tol       = m_options.getReal("Tolerance");
    const int maxIter = m_options.getInt ("MaxIter");
    const bool verbose = m_options.getSwitch("Verbose");
    const std::string logFilePath = m_options.getString("LogFile");

    m_optMesh.options().setReal("Penalty",   m_options.getReal("Penalty"));
    m_optMesh.options().setReal("Smoothing", m_options.getReal("Smoothing"));

    gsVector<T> u;
    _getControls(u);
    const index_t n = u.size();

    // A positive coefficient bound PROVES sigma is unfolded (non-negative basis,
    // partition of unity). It is conservative, so a negative value is not proof of
    // a fold - warn rather than abort, and let the finiteness check below catch an
    // actual folded start.
    const T certJ = m_comp->minDetJCoefficient();
    if (certJ <= T(0))
        gsWarn << "gsAdaptiveParametrizationNewton: det(J_sigma) coefficient bound is "
               << certJ << " <= 0. This bound is conservative, so it is not proof of a "
                  "fold; but if sigma IS folded the assembly below yields NaN.\n";

    gsVector<T> R(n);
    gsSparseMatrix<T> K;

    // Log matrix: rows = [energy, ||R||, minDet, alpha, wallTime].
    // Pre-allocate maxIter columns (trimmed to the used count on exit) to
    // avoid the O(maxIter^2) per-iteration conservativeResize copies.
    m_iterLog.resize(5, math::max(maxIter, 1));

    std::ofstream logStream;
    if (!logFilePath.empty())
    {
        logStream.open(logFilePath);
        logStream << "iter,energy,norm_R,minDet,alpha,wallTime\n";
    }

    // Per-phase profiling (gsStopwatch): accumulated wall time per phase.
    double t_residual = 0, t_eval = 0, t_assemble = 0, t_solve = 0,
           t_linesearch = 0;
    index_t n_backtracks = 0;   // trial steps = fused sweeps in the line search
    gsStopwatch sw;

    // Persistent factorisation: the K sparsity pattern is iteration-independent
    // (and the Levenberg shift only touches existing diagonal entries), so the
    // symbolic analyzePattern is done once and only factorize() repeats.
    typename gsSparseSolver<T>::LU solver;
    bool analyzed = false;

    gsStopwatch timer;
    timer.restart();

    T tau_lm = T(0);  // Levenberg–Marquardt shift (increases when K is indefinite)

    if (verbose)
    {
        gsInfo << "\n" << (useHessian ? "[Newton]" : "[Picard]") << " starting.\n";
        gsInfo << std::setw(5)  << "Iter"
               << std::setw(16) << "Energy"
               << std::setw(14) << "||R||"
               << std::setw(12) << "minDet"
               << std::setw(10) << "alpha"
               << std::setw(10) << "time[s]"
               << "\n";
        gsInfo << std::string(67, '-') << "\n";
    }

    index_t iter = 0;
    for (; iter < maxIter; ++iter)
    {
        sw.restart();
        _residual(u, R);
        t_residual += sw.stop();
        // Only checked on the starting mesh: the message below names it
        // explicitly, and a fold introduced by a later step is instead caught
        // by the line search's own minJ > 0 gate (see _armijoStep).
        if (iter == 0)
            GISMO_ENSURE(R.allFinite(),
                "gsAdaptiveParametrizationNewton: residual is not finite at the starting mesh "
                "- sigma is folded or the geometry is degenerate.");
        const T normR = R.norm();

        // One fused sweep yields both the energy (line-search base f0 + log) and
        // the min Jacobian (fold log) for the current iterate.
        sw.restart();
        T E, minJ;
        _evalObjAndMinJ(u, E, minJ);
        t_eval += sw.stop();

        // A folded start (det J_sigma < 0) gives det C = (det J)^2 > 0, so
        // C = J^T J stays invertible and R.allFinite() above passes even
        // though the line search can never make progress: every trial step
        // is rejected because it cannot improve on an already-inconsistent
        // gradient, not because of a NaN. Only an EXACTLY singular start
        // (det J == 0) produces the NaN the allFinite ensure above catches.
        if (iter == 0 && minJ <= T(0))
            gsWarn << "gsAdaptiveParametrizationNewton: sigma is folded at the starting "
                      "mesh (min det J_sigma = " << minJ << "); the line search cannot "
                      "make progress from a folded start.\n";

        const T t = timer.elapsed();

        // Record log row (columns pre-allocated)
        m_iterLog(0, iter) = E;
        m_iterLog(1, iter) = normR;
        m_iterLog(2, iter) = minJ;
        m_iterLog(3, iter) = T(1);  // alpha (filled after the line search)
        m_iterLog(4, iter) = t;

        if (normR < tol)
        {
            if (verbose)
                gsInfo << std::setw(5)  << iter
                       << std::setw(16) << E
                       << std::setw(14) << normR
                       << std::setw(12) << minJ
                       << std::setw(10) << "-"
                       << std::setw(10) << t
                       << "\n";
            break;
        }

        // Assemble tangent matrix
        sw.restart();
        if (useHessian)
            assembleHessian(u, K);
        else
            assemblePicard(u, K);
        t_assemble += sw.stop();

        // Solve: K Δ = -R. The pattern is constant, so analyzePattern() runs
        // once and only factorize() is repeated. _makeDefinite applies the
        // current Levenberg shift to the (existing) diagonal before factoring.
        sw.restart();
        if (!analyzed) { solver.analyzePattern(K); analyzed = true; }
        _makeDefinite(K, tau_lm);
        solver.factorize(K);
        // GISMO_ENSURE, not GISMO_ASSERT: solver.info() reports a RUNTIME
        // condition (the factorisation of a hard/singular K), not a programmer
        // error, so it must survive -DNDEBUG.  Compiled out, a failed
        // factorisation let solver.solve() return garbage or NaN, and the
        // descent test below (R.dot(delta) >= 0) then ACCEPTS a NaN step,
        // because every IEEE comparison with NaN is false.
        GISMO_ENSURE(solver.info() == 0,
            "gsAdaptiveParametrizationNewton::_solveLoop: linear solve "
            "factorisation failed (info = " << solver.info() << ") at "
            "Levenberg shift tau = " << tau_lm << "; the "
            << (useHessian ? "Newton" : "Picard") << " system matrix is "
            "singular or not factorisable.");
        gsVector<T> delta = solver.solve(-R);

        // Check descent: R·delta should be < 0 for a descent direction.
        // If not, add Levenberg shift and retry.
        if (R.dot(delta) >= T(0))
        {
            if (tau_lm == T(0))
                tau_lm = m_options.getReal("LevMarTau0");
            else
                tau_lm *= T(10);

            _makeDefinite(K, tau_lm);
            solver.factorize(K);
            delta = solver.solve(-R);
        }
        else
        {
            // Good descent direction: reduce shift (but not below min)
            tau_lm = math::max(tau_lm * T(0.1), T(0));
        }
        t_solve += sw.stop();

        // Armijo line search (returns 0 on failure). Each trial step is one
        // fused sweep (energy + fold guard).
        sw.restart();
        T alpha = _armijoStep(u, R, delta, E, &n_backtracks);
        t_linesearch += sw.stop();

        if (alpha > T(0))
        {
            u += alpha * delta;
            _setControls(u);
        }
        else
        {
            // Newton direction failed: retry with steepest descent (-R).
            // slope = R·(-R) = -||R||² < 0 → Armijo satisfiable from fold-free point.
            const gsVector<T> sdDir = -R;
            T alpha_sd = _armijoStep(u, R, sdDir, E, &n_backtracks);
            if (alpha_sd > T(0))
            {
                u += alpha_sd * sdDir;
                _setControls(u);
                alpha = alpha_sd;
            }
            else
            {
                // Even steepest descent failed (near fold boundary). Grow LM shift.
                tau_lm = (tau_lm == T(0)) ? m_options.getReal("LevMarTau0")
                                          : tau_lm * T(10);
            }
        }

        m_iterLog(3, iter) = alpha;

        if (verbose)
            gsInfo << std::setw(5)  << iter
                   << std::setw(16) << E
                   << std::setw(14) << normR
                   << std::setw(12) << minJ
                   << std::setw(10) << alpha
                   << std::setw(10) << t
                   << "\n";

        if (logStream.is_open())
            logStream << iter << "," << E << "," << normR << ","
                      << minJ << "," << alpha << "," << t << "\n";
    }

    // Every branch above that advances u also calls _setControls(u); a total
    // line-search failure (both Newton and steepest-descent directions
    // rejected) instead leaves the controls at u0 (restored by _armijoStep).
    // Re-assert here unconditionally so every return path leaves the domain
    // holding the accepted iterate u, never a rejected trial point.
    _setControls(u);

    // Certified (sampling-free) lower bound on det J_sigma of the accepted
    // iterate; conservative, so a warning here is not proof of a fold, but a
    // clean pass is a genuine unfolded-ness certificate (see minDetJCoefficient).
    {
        const T certJfinal = m_comp->minDetJCoefficient();
        if (certJfinal <= T(0))
            gsWarn << "gsAdaptiveParametrizationNewton: certified lower bound on "
                      "det J_sigma is <= 0 after relocation (bound = " << certJfinal
                   << "); the map may be folded (the bound is conservative).\n";
    }

    // Trim the pre-allocated log to the number of recorded iterations.
    const index_t usedCols = (iter < maxIter) ? (iter + 1) : maxIter;
    if (m_iterLog.cols() != usedCols)
        m_iterLog.conservativeResize(5, usedCols);

    if (verbose)
    {
        const double t_total = t_residual + t_eval + t_assemble + t_solve
                             + t_linesearch;
        const double inv = (t_total > 0) ? 100.0 / t_total : 0.0;
        gsInfo << std::string(67, '-') << "\n";
        gsInfo << (useHessian ? "[Newton]" : "[Picard]")
               << " profile  (N = " << n << " DoFs, " << iter << " iters)\n";
        // Deterministic work count: each line-search trial step is ONE fused
        // sweep (energy + fold guard fused), versus two separate sweeps before.
        gsInfo << "  line-search sweeps : " << n_backtracks
               << " (1 fused sweep each; was 2)\n";
        auto line = [&](const char * name, double tt)
        {
            gsInfo << std::setw(16) << name << std::setw(12) << tt << " s"
                   << std::setw(9) << (tt * inv) << " %\n";
        };
        line("residual",   t_residual);
        line("eval(E+minJ)",t_eval);
        line("assemble",   t_assemble);
        line("solve",      t_solve);
        line("line search",t_linesearch);
        line("total",      t_total);
        gsInfo << std::string(67, '-') << "\n";
    }

    return iter;
}

// ---------------------------------------------------------------------------
// Public solver entry points
// ---------------------------------------------------------------------------

template<class T, enum MonitorMode MODE>
index_t gsAdaptiveParametrizationNewton<T,MODE>::solveNewton()
{
    return _solveLoop(/*useHessian=*/true);
}

template<class T, enum MonitorMode MODE>
index_t gsAdaptiveParametrizationNewton<T,MODE>::solvePicard()
{
    return _solveLoop(/*useHessian=*/false);
}

} // namespace gismo
