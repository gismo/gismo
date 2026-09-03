/** @file gsFoldBarrier.h

    @brief Fold-avoidance barrier over the controls of a \a gsSquareDomain
    "sigma" map: a differentiable penalty on det(J_sigma) plus a soft box
    penalty on the controls, used by the D-step objectives (\a gsOptFit,
    \a gsOptL2) to keep sigma bijective during optimization.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-)
*/

#pragma once

#include <gsNurbs/gsSquareDomain.h>

namespace gismo
{

/// Which representation of det(J_sigma) the fold barrier penalises.
enum class gsFoldBarrierMode
{
    /// mu * sum_q w_q * max(0, eps - det J_sigma(x_q))^2, with (x_q,w_q) a
    /// Gauss rule on sigma's own knot mesh (order set by quB). The
    /// historical mode.
    Sampled     = 0,
    /// mu * sum_i max(0, eps - c_i)^2 over the per-element Bezier
    /// coefficients c_i of the EXACT det(J_sigma) spline. No quadrature, no
    /// quB.
    Coefficient = 1
};

/**
 * @brief Fold-avoidance barrier plus soft box penalty over the controls of a
 *   \a gsSquareDomain map sigma.
 *
 * Owns the sigma map (a non-owning pointer -- see the hazards below), mu,
 * eps, and whatever per-mode data the chosen \a gsFoldBarrierMode needs to
 * evaluate the barrier and its analytic gradient without repeating a basis
 * evaluation or a spline-arithmetic construction on every optimizer
 * iteration. Mode design and per-call cost comparison: \ref adaptparam_foldbarrier.
 *
 * addObj()/addGrad() read the held domain's CURRENT state directly
 * (coefficients for the fold term, getControls() for the box term) instead
 * of taking a design vector, so a setCtrls()-before-addObj() ordering
 * requirement is structurally impossible rather than merely documented.
 */
template <class T>
class gsFoldBarrier
{
public:

    /// Inactive barrier (mu <= 0): addObj()/addGrad() are no-ops,
    /// minDetJacobian() returns +infinity, nPoints() returns 0.
    gsFoldBarrier();

    /**
     * @param domain sigma; HELD as a member (non-owning pointer) -- addObj()
     *   / addGrad() read its CURRENT coefficients directly, so there is no
     *   setCtrls()-before-addObj() ordering contract for a caller to violate.
     * @param mu     barrier weight; mu <= 0 disables both the fold term and
     *   the soft box penalty.
     * @param eps    det(J_sigma) floor.
     * @param mode   see \a gsFoldBarrierMode.
     * @param quB    Gauss order for \a Sampled mode; IGNORED in \a Coefficient
     *   mode (the coefficient barrier has no quadrature at all).
     * @pre quB >= 0 (GISMO_ENSURE'd here): a negative quB reaches
     *   gsQuadrature::get() unguarded and segfaults, including under
     *   -DNDEBUG release builds. Checked in both modes.
     * @pre domain.domainDim() == 2 (GISMO_ENSURE'd here): the determinant
     *   formula (both modes) and the Coefficient mode's degree arithmetic
     *   are hard-coded for a 2D sigma domain.
     *
     * @warning The basis data cached here goes STALE if sigma's basis is
     *   refined AFTER this call: a caller that refines sigma after building
     *   the D-step objective hits it silently (no crash, a barrier evaluated
     *   on the wrong basis).
     */
    gsFoldBarrier(gsSquareDomain<T> & domain, T mu, T eps,
                  gsFoldBarrierMode mode, index_t quB);

    /// Adds the fold term and the soft box penalty on the controls to \a E.
    /// Reads the held domain's CURRENT coefficients and controls directly --
    /// no design-vector argument, hence no ordering contract with a caller's
    /// setCtrls(u) (see the class doc).
    void addObj(T & E) const;

    /// Adds their analytic gradient to \a result (size == the held domain's
    /// number of free controls).
    void addGrad(gsAsVector<T> & result) const;

    /// Sampled: minimum over the barrier quadrature points of det(J_sigma).
    /// Coefficient: min_i c_i, a certified lower bound on det(J_sigma)
    /// (same partition-of-unity argument as gsSquareDomain::minDetJCoefficient()).
    /// Available regardless of \a active() -- the cache this reads is built
    /// unconditionally by the constructor.
    T minDetJacobian() const;

    /// Number of evaluation points (Sampled: quadrature points; Coefficient:
    /// per-element Bezier coefficients of det(J_sigma)). 0 for a
    /// default-constructed (inactive, no domain) barrier.
    index_t nPoints() const;

    gsFoldBarrierMode mode() const { return m_mode; }

    /// True iff mu > 0, i.e. addObj()/addGrad() actually contribute.
    bool active() const { return m_mu > 0; }

private:

    /// Builds the Gauss rule on sigma's knot mesh and caches sigma's basis
    /// derivatives at the quadrature points (Sampled mode).
    void initSampled(index_t quB);

    /// Builds the two extraction operators and the two binomial weight
    /// tables of the Leibniz product (Coefficient mode). Kronecker/Leibniz
    /// derivation: \ref adaptparam_foldbarrier.
    void initCoefficient();

    /// Coefficient mode: the four mat-vecs a=Au*X, b=Au*Y, a'=Av*X, b'=Av*Y
    /// that turn sigma's CURRENT control coefficients into the Bezier-form
    /// derivative fields the per-element convolution consumes. Writes into
    /// the m_aFull/m_bFull/m_apFull/m_bpFull scratch members instead of
    /// returning by value, so repeated calls reuse those buffers.
    /// @note Allocates O(nb) per call (copies sigma's two coefficient
    ///   columns into local nb-sized vectors X, Y).
    void computeExtraction() const;

    /// Coefficient mode: the per-element Bernstein convolution (Leibniz
    /// product of the extracted derivative fields), assembled into the m_c
    /// scratch member (the flat per-element Bezier coefficient vector of
    /// det(J_sigma)). Calls computeExtraction() itself.
    void computeCoefficients() const;

    /// Sampled mode: the fold-term contribution to addGrad().
    void addSampledGrad(gsAsVector<T> & result) const;

    /// Coefficient mode: the fold-term contribution to addGrad() -- the
    /// adjoint of the per-element convolution, pulled back through the
    /// cached extraction operators.
    void addCoefficientGrad(gsAsVector<T> & result) const;

    gsSquareDomain<T> * m_domain = nullptr;

    T m_mu  = 0; // <=0 disables the fold term and the box penalty
    T m_eps = 0;

    gsFoldBarrierMode m_mode = gsFoldBarrierMode::Sampled;

    // --- Sampled mode ---
    gsMatrix<T>       m_barrierPts;
    gsVector<T>       m_barrierWts;
    gsMatrix<index_t> m_bpActives;
    gsMatrix<T>       m_bpDerivs;

    // Sampled-mode per-call scratch (dsig/det/dDet in addObj/addGrad/
    // minDetJacobian). MUTABLE and reused across calls: addObj/addGrad are
    // always called single-threaded, outside any OpenMP region, by both
    // gsOptFit and gsOptL2 (the parallel element sweeps are the DATA-term
    // sweeps in those classes; the barrier itself is never split across
    // threads), so a shared scratch buffer introduces no race. The various
    // *_into() helpers below call resize() on these internally, which is a
    // no-op once the size stabilises after the first call, turning what
    // would otherwise be a per-call heap allocation into a one-time
    // allocation at construction.
    mutable gsMatrix<T> m_dsig, m_det, m_dDet;

    // --- Coefficient mode ---
    // Sigma's degrees and per-direction element counts, cached for the
    // closed-form element indexing of the Bezier-full extraction operators.
    // Derivation: \ref adaptparam_foldbarrier.
    short_t m_degU = 0, m_degV = 0;
    index_t m_nelU = 0, m_nelV = 0;
    index_t m_nC   = 0; // total number of per-element Bezier coefficients of det(J_sigma)
    index_t m_nu = 0, m_nv = 0; // sizes of sigma's two 1-D basis factors (nu*nv == nb)

    // The two extraction operators Au (sigma-coeffs -> Bezier d(sigma)/du)
    // and Av (-> Bezier d(sigma)/dv) are Kronecker products
    // Au = m_Ev (x) m_Du, Av = m_Fv (x) m_Cu. Only the four 1-D factors are
    // stored, never the nC x nb 2-D operator. Derivation: \ref adaptparam_foldbarrier.
    //   m_Du: (nelU*degU)     x nu -- Bezier extraction of d/du (degU -> degU-1)
    //   m_Ev: (nelV*(degV+1)) x nv -- Bezier extraction, no derivative
    //   m_Cu: (nelU*(degU+1)) x nu -- Bezier extraction, no derivative
    //   m_Fv: (nelV*degV)     x nv -- Bezier extraction of d/dv (degV -> degV-1)
    gsMatrix<T> m_Du, m_Ev, m_Cu, m_Fv;

    // Per-direction binomial weight tables of the Leibniz product's two
    // terms (a*b' and a'*b have DIFFERENT operand bidegrees, hence two
    // tables). Derivation: \ref adaptparam_foldbarrier.
    gsMatrix<T> m_w1u, m_w1v, m_w2u, m_w2v;

    // Coefficient-mode per-call/per-element scratch, sized ONCE in
    // initCoefficient() and reused thereafter. m_aFull/m_bFull/m_apFull/
    // m_bpFull hold computeExtraction()'s output; m_c holds
    // computeCoefficients()'s flat per-element coefficient vector;
    // m_muA/m_muB/m_muAp/m_muBp are addCoefficientGrad()'s adjoint
    // accumulators; m_cLoc/m_lambdaLoc are the per-element local blocks;
    // m_gX/m_gY hold addCoefficientGrad()'s pulled-back gradient before the
    // dof-mapper scatter.
    mutable gsMatrix<T> m_aFull, m_bFull, m_apFull, m_bpFull;
    mutable gsMatrix<T> m_c;
    mutable gsMatrix<T> m_muA, m_muB, m_muAp, m_muBp;
    mutable gsMatrix<T> m_cLoc, m_lambdaLoc;
    mutable gsMatrix<T> m_gX, m_gY;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFoldBarrier.hpp)
#endif
