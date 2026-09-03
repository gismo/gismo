/** @file gsFoldBarrier.hpp

    @brief Implementation of the gsFoldBarrier class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-)
*/

#pragma once

#include <gsAssembler/gsQuadrature.h>
#include <gsCore/gsFuncData.h>
#include <gsMatrix/gsAsMatrix.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <cstdint>
#include <limits>

namespace gismo
{

namespace
{
/// Binomial coefficient C(n,k), duplicated from the anonymous-namespace
/// helper in gsTensorBSpline.hpp (not exported there) under a distinguishing
/// name to avoid any risk of an unnamed-namespace collision in a
/// translation unit that includes both headers.
///
/// Accumulates in a fixed-width 64-bit integer regardless of index_t's own
/// width: after folding k = min(k,n-k), the loop invariant result == C(n,i+1)
/// holds exactly after iteration i (C(n,i)*(n-i) is always divisible by
/// (i+1)), but the pre-division intermediate transiently reaches
/// C(n,i+1)*(i+1), strictly larger than the final C(n,k). Since C(n,m) is
/// increasing for m <= n/2, that intermediate is maximised at the last
/// iteration, i.e. at C(n,k)*k with k = floor(n/2) -- which first exceeds the
/// signed 64-bit range at n = 62 (C(62,31)*31 ~= 1.44e19 > INT64_MAX ~=
/// 9.22e18; n = 61 is the largest safe value, C(61,30)*30 ~= 6.98e18).
inline std::int64_t gsFoldBarrierBinom(index_t n, index_t k)
{
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    GISMO_ENSURE(n <= 61, "gsFoldBarrierBinom: n = " << n << " exceeds 61, "
                 "the largest value for which the intermediate product "
                 "C(n,k)*k stays within the 64-bit accumulator");
    std::int64_t result = 1;
    for (index_t i = 0; i < k; ++i)
    {
        result *= static_cast<std::int64_t>(n - i);
        result /= static_cast<std::int64_t>(i + 1);
    }
    return result;
}
} // anonymous namespace

template<class T>
gsFoldBarrier<T>::gsFoldBarrier()
: m_domain(nullptr), m_mu(0), m_eps(0), m_mode(gsFoldBarrierMode::Sampled)
{}

template<class T>
gsFoldBarrier<T>::gsFoldBarrier(gsSquareDomain<T> & domain, T mu, T eps,
                                gsFoldBarrierMode mode, index_t quB)
: m_domain(&domain), m_mu(mu), m_eps(eps), m_mode(mode)
{
    // GISMO_ENSURE, not GISMO_ASSERT: asserts are inert in build_rel, and
    // the determinant formula (both modes) and the Coefficient mode's degree
    // arithmetic below are hard-coded for 2D, so a non-2D domain would
    // otherwise silently compute a wrong det J_sigma instead of failing.
    GISMO_ENSURE(m_domain->domainDim()==2, "gsFoldBarrier: the fold-barrier "
                 "determinant formula is hard-coded for a 2D sigma domain");
    // GISMO_ENSURE, not GISMO_ASSERT: a negative quB reaches
    // gsQuadrature::get() unguarded and segfaults, including in build_rel
    // (-DNDEBUG) release builds. Checked in BOTH modes -- Coefficient mode
    // never uses quB, but a caller must not be able to rely on the mode to
    // skip the check.
    GISMO_ENSURE(quB >= 0, "gsFoldBarrier: barrier quB must be >= 0 (got "
                 << quB << ")");

    if (m_mode == gsFoldBarrierMode::Sampled)
        initSampled(quB);
    else
        initCoefficient();
}

// ---------------------------------------------------------------------------
// Sampled mode: Gauss quadrature on sigma's own knot mesh. There is no
// uniform-grid alternative here: every caller always requests the Gauss
// rule, so a uniform-grid path would be dead code.
// ---------------------------------------------------------------------------
template<class T>
void gsFoldBarrier<T>::initSampled(index_t quB)
{
    // Gauss quadrature on sigma's own knot mesh. det J_sigma is piecewise
    // polynomial exactly on this mesh, so delegating the rule to the basis
    // makes the barrier resolution follow sigmaKnots automatically (same
    // delegation as the integration basis of gsAdaptiveParametrization; the
    // union of bases degenerates to sigma's own basis because the integrand
    // involves sigma only).
    const gsBasis<T> & sb = m_domain->domain().basis();
    gsOptionList quadOptions;
    quadOptions.addReal("quA", "", 1.0);
    quadOptions.addInt ("quB", "", quB);
    gsQuadRule<T> rule = gsQuadrature::get(sb, quadOptions);
    std::vector<gsMatrix<T>> pts;
    std::vector<gsVector<T>> wts;
    gsMatrix<T> p;
    gsVector<T> w;
    index_t M = 0;
    for (auto & elem : sb.domain()->allElements())
    {
        rule.mapTo(elem.lowerCorner(), elem.upperCorner(), p, w);
        pts.push_back(p);
        wts.push_back(w);
        M += p.cols();
    }
    m_barrierPts.resize(2, M);
    m_barrierWts.resize(M);
    index_t o = 0;
    for (size_t e = 0; e != pts.size(); e++)
    {
        m_barrierPts.middleCols(o, pts[e].cols()) = pts[e];
        m_barrierWts.segment(o, wts[e].rows()) = wts[e];
        o += pts[e].cols();
    }

    // Cache sigma's own basis evaluation (active functions + basis
    // derivatives) at the barrier points. Only sigma's CONTROL COEFFICIENTS
    // change per optimization iteration, not these points, so the expensive
    // part of deriv_into -- element lookup and basis evaluation, redone from
    // scratch on every call otherwise -- needs to happen only once here.
    gsFuncData<T> bpData(NEED_DERIV | NEED_ACTIVE);
    sb.compute(m_barrierPts, bpData);
    m_bpActives = bpData.actives;
    m_bpDerivs  = bpData.values[1];
}

// ---------------------------------------------------------------------------
// Coefficient mode: det(J_sigma) = d_u(sigma_0) d_v(sigma_1)
// - d_v(sigma_0) d_u(sigma_1), with each derivative+Bezier-extraction an
// exact linear map of sigma's control coefficients that factors as a
// Kronecker product of two 1-D operators (sigma's basis is a tensor
// product). Full derivation of the Kronecker factorisation, the adjoint used
// by addCoefficientGrad(), and the closed-form per-element Bezier indexing:
// \ref adaptparam_foldbarrier.
// ---------------------------------------------------------------------------

namespace
{
/// Weight table of ONE direction of ONE Leibniz term:
/// table(i,k) = C(degA,i) * C(degB,k-i) / C(degA+degB,k), zero where
/// k-i is outside [0,degB]. Two such tables (one per direction) combine by a
/// product into the full 2D binomial weight multiply() uses.
template<class T>
void buildFoldBarrierWeightTable(index_t degA, index_t degB, gsMatrix<T> & table)
{
    const index_t degR = degA + degB;
    table.setZero(degA + 1, degR + 1);
    for (index_t i = 0; i <= degA; ++i)
        for (index_t k = 0; k <= degR; ++k)
        {
            const index_t j = k - i;
            if (j < 0 || j > degB) continue;
            table(i, k) = static_cast<T>(gsFoldBarrierBinom(degA, i))
                        * static_cast<T>(gsFoldBarrierBinom(degB, j))
                        / static_cast<T>(gsFoldBarrierBinom(degR, k));
        }
}
} // anonymous namespace

namespace
{
/// 1-D analogue of gsTensorBSpline<d,T>::grad(dir), specialised to a single
/// direction. Differentiates the degree-p B-spline curve (basis, coefs) and
/// returns the degree-(p-1) result: d_i = p*(c_{i+1}-c_i)/(t_{i+p+1}-t_{i+1})
/// per coefficient, derivative knot vector KV[1:-1].
template<class T>
gsBSpline<T> gsFoldBarrierDeriv1D(const gsBSplineBasis<T> & basis, const gsMatrix<T> & coefs)
{
    const short_t p = basis.degree();
    const gsKnotVector<T> & kv = basis.knots();
    const index_t n = coefs.rows();
    gsMatrix<T> dcoefs(n - 1, coefs.cols());
    for (index_t i = 0; i < n - 1; ++i)
    {
        const T denom = kv[i + p + 1] - kv[i + 1];
        if (denom == T(0))
            dcoefs.row(i).setZero();
        else
            dcoefs.row(i) = (T(p) / denom) * (coefs.row(i + 1) - coefs.row(i));
    }
    const index_t szFull = static_cast<index_t>(kv.size());
    std::vector<T> kvNew;
    kvNew.reserve(static_cast<size_t>(szFull - 2));
    for (index_t i = 1; i < szFull - 1; ++i)
        kvNew.push_back(kv[i]);
    gsKnotVector<T> dkv(kvNew, static_cast<short_t>(p - 1));
    gsBSplineBasis<T> dbasis(dkv);
    return gsBSpline<T>(dbasis, give(dcoefs));
}

/// 1-D analogue of gsTensorBSpline<d,T>::toBezier(), specialised to a single
/// direction: raises every interior knot to full multiplicity degree+1 via
/// Boehm knot insertion, so consecutive (degree+1)-blocks of the returned
/// coefficient matrix are element e's Bezier control points. Deliberately
/// NOT gsBSpline<T>::toBezier(): that returns a gsMultiPatch of split
/// segments, with no single flat coefficient matrix to read off.
template<class T>
gsMatrix<T> gsFoldBarrierBezier1D(gsBSpline<T> B)
{
    const short_t p = B.degree();
    const gsKnotVector<T> & kv0 = B.knots();
    const T first = kv0.first(), last = kv0.last();
    std::vector<T> interior;
    interior.reserve(kv0.uSize());
    for (size_t i = 0; i != kv0.uSize(); ++i)
    {
        const T xi = kv0.uValue(i);
        if (xi > first && xi < last) interior.push_back(xi);
    }
    // gsBSpline<T>::insertKnot has a public 2-arg overload (knot, count) and
    // a PRIVATE virtual 3-arg override (knot, dir, count=1) -- with the
    // trailing default argument, a plain 2-argument call is ambiguous
    // between them (overload resolution ignores accessibility), so pin the
    // public overload's address explicitly by its exact 2-parameter type.
    void (gsBSpline<T>::*insertKnot2)(T, index_t) = &gsBSpline<T>::insertKnot;
    for (const T xi : interior)
    {
        const int mult = B.knots().multiplicity(xi);
        const int need = p + 1 - mult;
        if (need > 0) (B.*insertKnot2)(xi, static_cast<index_t>(need));
    }
    return B.coefs();
}

/// Builds ONE of the four 1-D Kronecker factors (m_Du/m_Ev/m_Cu/m_Fv): feeds
/// each of \a basis's n unit coefficient vectors through
/// gsFoldBarrierDeriv1D (iff \a differentiate) then gsFoldBarrierBezier1D,
/// and stores the resulting Bezier column. O(n) columns (n = basis.size()
/// = nu or nv), each O(degree) local work (Boehm insertion has local
/// support) -- O(n) total setup, against the O(nu*nv) columns (each an
/// O(nu*nv)-sized Bezier grid) a directly-formed 2-D Au/Av would need.
template<class T>
gsMatrix<T> buildFoldBarrier1DOperator(const gsBSplineBasis<T> & basis, bool differentiate)
{
    const index_t n = static_cast<index_t>(basis.size());
    gsMatrix<T> e = gsMatrix<T>::Zero(n, 1);
    gsMatrix<T> op;
    for (index_t k = 0; k != n; ++k)
    {
        e(k,0) = T(1);
        const gsMatrix<T> col = differentiate
            ? gsFoldBarrierBezier1D(gsFoldBarrierDeriv1D(basis, e))
            : gsFoldBarrierBezier1D(gsBSpline<T>(basis, e));
        if (k == 0) op.resize(col.rows(), n);
        op.col(k) = col.col(0);
        e(k,0) = T(0);
    }
    return op;
}
} // anonymous namespace

template<class T>
void gsFoldBarrier<T>::initCoefficient()
{
    const gsTensorBSplineBasis<2,T> * sbasis =
        dynamic_cast<const gsTensorBSplineBasis<2,T> *>(&m_domain->domain().basis());
    GISMO_ENSURE(sbasis != nullptr, "gsFoldBarrier: the coefficient barrier "
                 "requires sigma's basis to be a gsTensorBSplineBasis<2,T>");

    m_degU = sbasis->degree(0);
    m_degV = sbasis->degree(1);
    // d_u, d_v require degree >= 1 in the differentiated direction; a
    // degree-0 sigma direction cannot happen for a legal bijective
    // reparametrization anyway, but gsFoldBarrierDeriv1D's derivative
    // formula itself needs p >= 1, so make it a hard failure here instead.
    GISMO_ENSURE(m_degU >= 1 && m_degV >= 1, "gsFoldBarrier: the coefficient "
                 "barrier's degree arithmetic (d/du, d/dv of sigma) requires "
                 "degree >= 1 in both directions, got " << m_degU << "x" << m_degV);
    m_nelU = static_cast<index_t>(sbasis->knots(0).numElements());
    m_nelV = static_cast<index_t>(sbasis->knots(1).numElements());

    const gsBSplineBasis<T> & basisU = sbasis->component(0);
    const gsBSplineBasis<T> & basisV = sbasis->component(1);
    m_nu = static_cast<index_t>(basisU.size());
    m_nv = static_cast<index_t>(basisV.size());

    const index_t nb = static_cast<index_t>(sbasis->size());
    // Sigma's basis is a plain tensor product of its two 1-D factors --
    // gsTensorBSplineBasis guarantees nb == nu*nv by construction, so this
    // is a defensive check (the Kronecker factorisation below is invalid
    // without it), not a live branch.
    GISMO_ENSURE(m_nu * m_nv == nb, "gsFoldBarrier: sigma's basis is not a "
                 "plain tensor product of its two 1-D factors (nu*nv != nb) "
                 "-- the Kronecker factorisation this class relies on does "
                 "not hold");

    // Full Bezier grid sizes of the two extraction maps (d_u sigma has
    // bidegree (degU-1,degV); d_v sigma has bidegree (degU,degV-1)).
    const index_t nA0 = m_nelU * m_degU,       nA1 = m_nelV * (m_degV + 1);
    const index_t nB0 = m_nelU * (m_degU + 1), nB1 = m_nelV * m_degV;

    // The four 1-D Kronecker factors (Au = m_Ev(x)m_Du, Av = m_Fv(x)m_Cu --
    // see the derivation comment above and gsFoldBarrier.h's member doc).
    // Built once here from sigma's two knot vectors alone, O(nu+nv) columns
    // total against the O(nu*nv) columns (each an O(nu*nv)-sized Bezier
    // grid) a directly-formed 2-D Au/Av would need.
    m_Du = buildFoldBarrier1DOperator(basisU, /*differentiate*/true);
    m_Ev = buildFoldBarrier1DOperator(basisV, /*differentiate*/false);
    m_Cu = buildFoldBarrier1DOperator(basisU, /*differentiate*/false);
    m_Fv = buildFoldBarrier1DOperator(basisV, /*differentiate*/true);

    // Bezier layout identity, one per factor: catches a
    // knot-multiplicity/element-count mismatch at setup instead of a
    // silent misindex in computeCoefficients()/addCoefficientGrad().
    GISMO_ENSURE(m_Du.rows() == nA0 && m_Du.cols() == m_nu, "gsFoldBarrier: "
                 "Bezier layout identity failed for m_Du (" << m_Du.rows()
                 << "x" << m_Du.cols() << " != " << nA0 << "x" << m_nu << ")");
    GISMO_ENSURE(m_Ev.rows() == nA1 && m_Ev.cols() == m_nv, "gsFoldBarrier: "
                 "Bezier layout identity failed for m_Ev (" << m_Ev.rows()
                 << "x" << m_Ev.cols() << " != " << nA1 << "x" << m_nv << ")");
    GISMO_ENSURE(m_Cu.rows() == nB0 && m_Cu.cols() == m_nu, "gsFoldBarrier: "
                 "Bezier layout identity failed for m_Cu (" << m_Cu.rows()
                 << "x" << m_Cu.cols() << " != " << nB0 << "x" << m_nu << ")");
    GISMO_ENSURE(m_Fv.rows() == nB1 && m_Fv.cols() == m_nv, "gsFoldBarrier: "
                 "Bezier layout identity failed for m_Fv (" << m_Fv.rows()
                 << "x" << m_Fv.cols() << " != " << nB1 << "x" << m_nv << ")");

    // Two binomial weight tables, one per Leibniz term: term 1 is
    // a*b' = (d_u sigma_0)*(d_v sigma_1), operand bidegrees
    // (degU-1,degV) x (degU,degV-1); term 2 is a'*b = (d_v sigma_0)*(d_u sigma_1),
    // operand bidegrees (degU,degV-1) x (degU-1,degV) -- the SAME pair of
    // bidegrees with the roles swapped, which is why the weight differs
    // between the two terms even though both results share the bidegree
    // (2degU-1,2degV-1). Reusing one table for both would be a silent wrong
    // answer whenever degU != degV or the two directions are refined
    // differently.
    buildFoldBarrierWeightTable(m_degU-1, m_degU,   m_w1u);
    buildFoldBarrierWeightTable(m_degV,   m_degV-1, m_w1v);
    buildFoldBarrierWeightTable(m_degU,   m_degU-1, m_w2u);
    buildFoldBarrierWeightTable(m_degV-1, m_degV,   m_w2v);

    m_nC = (2*m_degU) * m_nelU * (2*m_degV) * m_nelV;

    // Scratch buffers, sized ONCE here and reused by every subsequent
    // addObj()/addGrad()/minDetJacobian() call (see the m_aFull etc. doc in
    // gsFoldBarrier.h): turns what would otherwise be O(nelU*nelV) fresh
    // allocations per call into a one-time allocation at construction.
    m_aFull.resize(nA0 * nA1, 1);  m_bFull.resize(nA0 * nA1, 1);
    m_apFull.resize(nB0 * nB1, 1); m_bpFull.resize(nB0 * nB1, 1);
    m_c.resize(m_nC, 1);
    m_muA.resize(nA0 * nA1, 1);  m_muB.resize(nA0 * nA1, 1);
    m_muAp.resize(nB0 * nB1, 1); m_muBp.resize(nB0 * nB1, 1);
    m_cLoc.resize(2*m_degU, 2*m_degV);
    m_lambdaLoc.resize(2*m_degU, 2*m_degV);
    m_gX.resize(nb, 1); m_gY.resize(nb, 1);
}

template<class T>
void gsFoldBarrier<T>::computeExtraction() const
{
    const gsMatrix<T> & coefs = m_domain->domain().coefs();
    // The reinterpretation below aliases coefs's own storage as an m_nu x
    // m_nv grid, so its validity depends on coefs.rows() == m_nu*m_nv -- a
    // statement about coefs specifically, distinct from the constructor's
    // nu*nv == nb ENSURE (initCoefficient(), gsFoldBarrier.hpp:360-363),
    // which is about sigma's BASIS size, not this per-call coefficient
    // matrix. A mismatch here would make Xmat/Ymat read out of bounds.
    GISMO_ENSURE(coefs.rows() == m_nu * m_nv, "gsFoldBarrier::computeExtraction: "
                 "coefs.rows() = " << coefs.rows() << " != m_nu*m_nv = "
                 << m_nu * m_nv);
    // coefs is column-major, so each column is contiguous and can be
    // reinterpreted in place as sigma's nu x nv coefficient grid, column-major
    // with u fastest -- the SAME flattening computeCoefficients()'s flat index
    // formula assumes for its output (see the derivation comment above
    // initCoefficient()). The .col() Block temporaries die at the end of this
    // statement, but .data() points into coefs, which outlives them.
    const gsAsConstMatrix<T> Xmat(coefs.col(0).data(), m_nu, m_nv);
    const gsAsConstMatrix<T> Ymat(coefs.col(1).data(), m_nu, m_nv);

    const index_t nA0 = m_Du.rows(), nA1 = m_Ev.rows();
    const index_t nB0 = m_Cu.rows(), nB1 = m_Fv.rows();

    // Au*X = vec(m_Du*Xmat*m_Ev^T), Av*X = vec(m_Cu*Xmat*m_Fv^T): the
    // Kronecker factor applied as two small dense products, m_aFull etc.
    // reinterpreted in place as the nA0 x nA1 (resp. nB0 x nB1) result --
    // the Kronecker product Au = m_Ev(x)m_Du is NEVER formed.
    m_aFull.reshapeCol(0, nA0, nA1).noalias()  = m_Du * Xmat * m_Ev.transpose();
    m_bFull.reshapeCol(0, nA0, nA1).noalias()  = m_Du * Ymat * m_Ev.transpose();
    m_apFull.reshapeCol(0, nB0, nB1).noalias() = m_Cu * Xmat * m_Fv.transpose();
    m_bpFull.reshapeCol(0, nB0, nB1).noalias() = m_Cu * Ymat * m_Fv.transpose();
}

template<class T>
void gsFoldBarrier<T>::computeCoefficients() const
{
    computeExtraction();

    const index_t degR0 = 2*m_degU, degR1 = 2*m_degV;
    const index_t nA0 = m_nelU*m_degU, nB0 = m_nelU*(m_degU+1);
    const index_t nC0 = m_nelU*degR0;

    for (index_t ev = 0; ev != m_nelV; ++ev)
    for (index_t eu = 0; eu != m_nelU; ++eu)
    {
        m_cLoc.setZero();
        // term 1: a * b'  (a bidegree (degU-1,degV), b' bidegree (degU,degV-1))
        for (index_t i1 = 0; i1 <= m_degV; ++i1)
        for (index_t i0 = 0;  i0 != m_degU;  ++i0)
        {
            const T aVal = m_aFull((eu*m_degU+i0) + nA0*(ev*(m_degV+1)+i1), 0);
            for (index_t j1 = 0; j1 != m_degV;  ++j1)
            for (index_t j0 = 0; j0 <= m_degU;  ++j0)
            {
                const T w = m_w1u(i0, i0+j0) * m_w1v(i1, i1+j1);
                if (w == T(0)) continue;
                const T bpVal = m_bpFull((eu*(m_degU+1)+j0) + nB0*(ev*m_degV+j1), 0);
                m_cLoc(i0+j0, i1+j1) += w * aVal * bpVal;
            }
        }
        // term 2: a' * b, subtracted  (a' bidegree (degU,degV-1), b bidegree (degU-1,degV))
        for (index_t i1 = 0; i1 != m_degV;  ++i1)
        for (index_t i0 = 0; i0 <= m_degU;  ++i0)
        {
            const T apVal = m_apFull((eu*(m_degU+1)+i0) + nB0*(ev*m_degV+i1), 0);
            for (index_t j1 = 0; j1 <= m_degV; ++j1)
            for (index_t j0 = 0; j0 != m_degU;  ++j0)
            {
                const T w = m_w2u(i0, i0+j0) * m_w2v(i1, i1+j1);
                if (w == T(0)) continue;
                const T bVal = m_bFull((eu*m_degU+j0) + nA0*(ev*(m_degV+1)+j1), 0);
                m_cLoc(i0+j0, i1+j1) -= w * apVal * bVal;
            }
        }

        for (index_t k1 = 0; k1 != degR1; ++k1)
            for (index_t k0 = 0; k0 != degR0; ++k0)
                m_c((eu*degR0+k0) + nC0*(ev*degR1+k1), 0) = m_cLoc(k0,k1);
    }
}

template<class T>
void gsFoldBarrier<T>::addCoefficientGrad(gsAsVector<T> & result) const
{
    computeExtraction();
    m_muA.setZero(); m_muB.setZero(); m_muAp.setZero(); m_muBp.setZero();

    const index_t degR0 = 2*m_degU, degR1 = 2*m_degV;
    const index_t nA0 = m_nelU*m_degU, nB0 = m_nelU*(m_degU+1);

    for (index_t ev = 0; ev != m_nelV; ++ev)
    for (index_t eu = 0; eu != m_nelU; ++eu)
    {
        // ---- forward pass: c^e = a*b' - a'*b (same convolution as
        // computeCoefficients(), redone per element here so the active set
        // -- and hence lambda -- is available for the backward pass below
        // without round-tripping through a global coefficient vector). ----
        m_cLoc.setZero();
        for (index_t i1 = 0; i1 <= m_degV; ++i1)
        for (index_t i0 = 0;  i0 != m_degU;  ++i0)
        {
            const T aVal = m_aFull((eu*m_degU+i0) + nA0*(ev*(m_degV+1)+i1), 0);
            for (index_t j1 = 0; j1 != m_degV;  ++j1)
            for (index_t j0 = 0; j0 <= m_degU;  ++j0)
            {
                const T w = m_w1u(i0, i0+j0) * m_w1v(i1, i1+j1);
                if (w == T(0)) continue;
                const index_t bpIdx = (eu*(m_degU+1)+j0) + nB0*(ev*m_degV+j1);
                m_cLoc(i0+j0, i1+j1) += w * aVal * m_bpFull(bpIdx,0);
            }
        }
        for (index_t i1 = 0; i1 != m_degV;  ++i1)
        for (index_t i0 = 0; i0 <= m_degU;  ++i0)
        {
            const T apVal = m_apFull((eu*(m_degU+1)+i0) + nB0*(ev*m_degV+i1), 0);
            for (index_t j1 = 0; j1 <= m_degV; ++j1)
            for (index_t j0 = 0; j0 != m_degU;  ++j0)
            {
                const T w = m_w2u(i0, i0+j0) * m_w2v(i1, i1+j1);
                if (w == T(0)) continue;
                const index_t bIdx = (eu*m_degU+j0) + nA0*(ev*(m_degV+1)+j1);
                m_cLoc(i0+j0, i1+j1) -= w * apVal * m_bFull(bIdx,0);
            }
        }

        // lambda_k = -2*mu*max(0,eps-c_k) -- zero off the active set, so an
        // element with no coefficient below eps contributes nothing and is
        // skipped entirely.
        bool anyActive = false;
        for (index_t k1 = 0; k1 != degR1; ++k1)
            for (index_t k0 = 0; k0 != degR0; ++k0)
            {
                const T v = std::max((T)0, m_eps - m_cLoc(k0,k1));
                m_lambdaLoc(k0,k1) = T(-2) * m_mu * v;
                if (v > 0) anyActive = true;
            }
        if (!anyActive) continue;

        // ---- backward pass: differentiate the bilinear form directly.
        // Term 1 (a*b') contributes to mu_a (paired with Au, no sign flip)
        // and mu_b' (paired with Av, no sign flip); term 2 (a'*b, entering c
        // with a MINUS) contributes to mu_a' and mu_b with the sign already
        // folded in, matching addObj()'s c = c1 - c2. ----
        for (index_t i1 = 0; i1 <= m_degV; ++i1)
        for (index_t i0 = 0;  i0 != m_degU;  ++i0)
        {
            const index_t aIdx = (eu*m_degU+i0) + nA0*(ev*(m_degV+1)+i1);
            T acc = 0;
            for (index_t j1 = 0; j1 != m_degV;  ++j1)
            for (index_t j0 = 0; j0 <= m_degU;  ++j0)
            {
                const T w = m_w1u(i0, i0+j0) * m_w1v(i1, i1+j1);
                if (w == T(0)) continue;
                const index_t bpIdx = (eu*(m_degU+1)+j0) + nB0*(ev*m_degV+j1);
                const T lam = m_lambdaLoc(i0+j0, i1+j1);
                acc             += w * lam * m_bpFull(bpIdx,0);
                m_muBp(bpIdx,0) += w * lam * m_aFull(aIdx,0);
            }
            m_muA(aIdx,0) += acc;
        }
        for (index_t i1 = 0; i1 != m_degV;  ++i1)
        for (index_t i0 = 0; i0 <= m_degU;  ++i0)
        {
            const index_t apIdx = (eu*(m_degU+1)+i0) + nB0*(ev*m_degV+i1);
            T acc = 0;
            for (index_t j1 = 0; j1 <= m_degV; ++j1)
            for (index_t j0 = 0; j0 != m_degU;  ++j0)
            {
                const T w = m_w2u(i0, i0+j0) * m_w2v(i1, i1+j1);
                if (w == T(0)) continue;
                const index_t bIdx = (eu*m_degU+j0) + nA0*(ev*(m_degV+1)+j1);
                const T lam = m_lambdaLoc(i0+j0, i1+j1);
                acc           -= w * lam * m_bFull(bIdx,0);
                m_muB(bIdx,0) -= w * lam * m_apFull(apIdx,0);
            }
            m_muAp(apIdx,0) += acc;
        }
    }

    // Pull back through the cached Kronecker factors (a=Au*X, a'=Av*X,
    // b=Au*Y, b'=Av*Y are all linear, so the adjoint is just the transpose
    // maps: Au^T*v = vec(m_Du^T*Vmat*m_Ev), Av^T*v = vec(m_Cu^T*Vmat*m_Fv)
    // -- see the derivation comment above initCoefficient() for the
    // Kronecker-transpose identity this uses; m_muA/m_muB/m_muAp/m_muBp
    // reshaped in place, same shapes as m_aFull/m_bFull/m_apFull/m_bpFull)
    // and scatter through the dof mapper. Accumulate (+=), never assign:
    // matchDof (gsSquareDomain::_initMapper) can map two active coefficients
    // to the same free index -- same reasoning as
    // gsSquareDomain::detJacobianDeriv_into's own scatter loop.
    // nA0, nB0 already declared above (same values as m_Du.rows()/m_Cu.rows()
    // respectively -- both equal m_nelU*m_degU / m_nelU*(m_degU+1) by the
    // Bezier layout identity GISMO_ENSURE'd in initCoefficient()).
    const index_t nA1 = m_Ev.rows(), nB1 = m_Fv.rows();
    const gsAsConstMatrix<T> muAmat(m_muA.data(), nA0, nA1);
    const gsAsConstMatrix<T> muBmat(m_muB.data(), nA0, nA1);
    const gsAsConstMatrix<T> muApMat(m_muAp.data(), nB0, nB1);
    const gsAsConstMatrix<T> muBpMat(m_muBp.data(), nB0, nB1);
    m_gX.reshapeCol(0, m_nu, m_nv).noalias() =
        m_Du.transpose() * muAmat * m_Ev + m_Cu.transpose() * muApMat * m_Fv;
    m_gY.reshapeCol(0, m_nu, m_nv).noalias() =
        m_Du.transpose() * muBmat * m_Ev + m_Cu.transpose() * muBpMat * m_Fv;

    const gsDofMapper & mapper = m_domain->mapper();
    const index_t nb = m_gX.rows();
    for (index_t k = 0; k != nb; ++k)
    {
        if (mapper.is_free(k, 0, 0)) result[mapper.index(k, 0, 0)] += m_gX(k,0);
        if (mapper.is_free(k, 0, 1)) result[mapper.index(k, 0, 1)] += m_gY(k,0);
    }
}

template<class T>
void gsFoldBarrier<T>::addSampledGrad(gsAsVector<T> & result) const
{
    m_domain->detJacobianDeriv_into(m_barrierPts, m_dDet, &m_dsig);
    gsSquareDomain<T>::detFromJacobian_into(m_dsig, m_domain->domainDim(), m_det);
    const index_t M = m_barrierPts.cols();
    for (index_t q = 0; q != M; q++)
    {
        const T v = std::max((T)0, m_eps - m_det(0,q));
        if (v > 0)
            for (index_t j = 0; j != result.rows(); j++)
                result[j] -= 2.0 * m_mu * m_barrierWts[q] * v * m_dDet(j,q);
    }
}

// ---------------------------------------------------------------------------
// Public interface
// ---------------------------------------------------------------------------

template<class T>
void gsFoldBarrier<T>::addObj(T & E) const
{
    if (m_mu <= 0) return;

    if (m_mode == gsFoldBarrierMode::Sampled)
    {
        gsBasis<T>::linearCombination_into(m_domain->domain().coefs(), m_bpActives, m_bpDerivs, m_dsig);
        gsSquareDomain<T>::detFromJacobian_into(m_dsig, m_domain->domainDim(), m_det);
        const index_t M = m_barrierPts.cols();
        for (index_t q = 0; q != M; q++)
        {
            const T v = std::max((T)0, m_eps - m_det(0,q));
            E += m_mu * m_barrierWts[q] * v * v;
        }
    }
    else // Coefficient
    {
        computeCoefficients();
        for (index_t i = 0; i != m_c.rows(); ++i)
        {
            const T v = std::max((T)0, m_eps - m_c(i,0));
            E += m_mu * v * v;
        }
    }

    const gsVector<T> u = m_domain->getControls();   // soft BOX penalty on the controls
    for (index_t j = 0; j != u.rows(); j++)
    {
        const T lo = std::max((T)0, -u[j]);
        const T hi = std::max((T)0, u[j] - (T)1);
        E += m_mu * (lo*lo + hi*hi);
    }
}

template<class T>
void gsFoldBarrier<T>::addGrad(gsAsVector<T> & result) const
{
    if (m_mu <= 0) return;

    if (m_mode == gsFoldBarrierMode::Sampled)
        addSampledGrad(result);
    else
        addCoefficientGrad(result);

    const gsVector<T> u = m_domain->getControls();
    for (index_t j = 0; j != u.rows(); j++)
        result[j] += 2.0 * m_mu * (std::max((T)0, u[j]-(T)1) - std::max((T)0, -u[j]));
}

template<class T>
T gsFoldBarrier<T>::minDetJacobian() const
{
    // Defensive default for a default-constructed gsFoldBarrier (m_domain == nullptr).
    if (m_domain == nullptr) return std::numeric_limits<T>::max();

    if (m_mode == gsFoldBarrierMode::Sampled)
    {
        gsBasis<T>::linearCombination_into(m_domain->domain().coefs(), m_bpActives, m_bpDerivs, m_dsig);
        gsSquareDomain<T>::detFromJacobian_into(m_dsig, m_domain->domainDim(), m_det);
        return m_det.minCoeff();
    }
    else
    {
        computeCoefficients();
        return m_c.minCoeff();
    }
}

template<class T>
index_t gsFoldBarrier<T>::nPoints() const
{
    return m_mode == gsFoldBarrierMode::Sampled ? m_barrierPts.cols() : m_nC;
}

} // namespace gismo
