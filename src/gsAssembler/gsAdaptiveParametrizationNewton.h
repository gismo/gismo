/** @file gsAdaptiveParametrizationNewton.h

    @brief Newton and Picard solvers for the composite spline relocation
    objective via exact analytic Hessian.

    See doc/derivation_hessian.md for the mathematical derivation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-)
*/

#pragma once

#include <gsAssembler/gsAdaptiveParametrization.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsIO/gsOptionList.h>
#include <gsCore/gsComposedGeometry.h>
#include <gsCore/gsMultiBasis.h>
#include <gsAssembler/gsQuadrature.h>
#include <gsMatrix/gsFiberMatrix.h>
#include <gsUtils/gsThreaded.h>

namespace gismo
{

/**
 * @brief Newton and Picard solvers for the composite r-adaptivity energy.
 *
 * Minimises the same energy as gsOptMesh (evalObj / gradObj_into) but via:
 * - **Newton**: exact sparse Hessian K assembled analytically (see [H-4]),
 *   modified-Newton (Levenberg shift when K is indefinite) + Armijo line search.
 * - **Picard**: SPD frozen-coefficient surrogate B (see [H-5]), same loop structure.
 *
 * Both solvers share:
 * - the same residual `R(σ) = gradObj_into(u)` as ground truth,
 * - the same per-iteration logger → CSV,
 * - fold detection via `computeMinJacobian`.
 *
 * The gsExprAssembler-based assembleResidual() independently assembles the EL
 * weak form and is validated against gradObj_into() to machine precision.
 *
 * Reference: doc/derivation_hessian.md  §[H-4], [H-5]
 *
 * @tparam T       Real scalar type
 * @tparam MODE    MonitorMode::ValueBased or MonitorMode::GradientBased
 */
template<class T, enum MonitorMode MODE>
class gsAdaptiveParametrizationNewton
{
public:

    /**
     * @brief Constructor without monitor function.
     *
     * Builds the same internal state as gsOptMesh (composition, geometry, basis)
     * so that gradObj_into can be called for residual evaluation.
     *
     * @note \a integrationBasis must not be coarser than the sigma basis
     * (\c composition.domain().basis()): assembly assumes one active set per
     * integration element (SAME_ELEMENT). Use
     * gsAdaptiveParametrization::makeIntegrationBasis to build a valid
     * union basis (geometry knots + composition knots).
     */
    gsAdaptiveParametrizationNewton(
              gsSquareDomain<T>  & composition,
        const gsGeometry<T>    & geometry,
        const gsBasis<T>       & integrationBasis);

    /**
     * @brief Constructor with monitor function.
     *
     * @param parametric  If true, the monitor f is evaluated in parametric
     *   coordinates (m_parametric); if false, in physical coordinates.
     */
    gsAdaptiveParametrizationNewton(
              gsSquareDomain<T>  & composition,
        const gsGeometry<T>    & geometry,
        const gsFunction<T>    & fun,
        const gsBasis<T>       & integrationBasis,
              bool               parametric);

    // ------------------------------------------------------------------
    // Options
    // ------------------------------------------------------------------

    gsOptionList & options() { return m_options; }

    static gsOptionList defaultOptions();

    // ------------------------------------------------------------------
    // Assembly
    // ------------------------------------------------------------------

    /**
     * @brief Assemble the EL residual as a weak form (flux + reaction).
     *
     * Hand-rolled assembly (OpenMP over elements, mirrors gsExprAssembler's
     * element sweep) of
     *   R_A^ℓ = ∫ Π_ℓ·∇β̂_A + r_ℓ β̂_A dΩ̂           [H-3.2]
     * The result equals gradObj_into() to machine precision (verified in
     * unittests) but is computed *independently* from the flux/reaction
     * decomposition — validating the weak form reused by the Hessian.
     *
     * @param u      Current free-DoF control vector
     * @param R_out  Output residual vector (size = nControls)
     */
    void assembleResidual(const gsVector<T> & u, gsVector<T> & R_out) const;

    /**
     * @brief Assemble the exact analytic sparse Newton tangent K = ∂²E/∂α².
     *
     * Differentiates the discrete gradient of gradObj_into() exactly:
     * material [H-4.1], geometric/curvature [H-4.2] (third derivatives of
     * the geometry; FD-of-deriv2 fallback for NURBS) and monitor [H-4.3]
     * blocks. Implemented for dd==2 (planar and surface); for dd==3 a
     * finite-difference fallback of the analytic gradient is used.
     *
     * Storage: two-phase gsFiberMatrix (cached sparsity pattern + OpenMP
     * atomic scatter), complexity O(n_qp · (dd·k)²) per assembly with
     * k = (deg+1)^dd local basis functions.
     *
     * Symmetry ‖K−Kᵀ‖ < 1e-12·‖K‖ is verified in unittests.
     *
     * @param u      Current free-DoF control vector
     * @param K_out  Output symmetric sparse matrix (nControls × nControls)
     */
    void assembleHessian(const gsVector<T> & u, gsSparseMatrix<T> & K_out) const;

    /**
     * @brief Assemble the Picard surrogate B (SPD frozen-coefficient operator).
     *
     * For q=2 (no-monitor): B is the exact linearised EL operator (SPD).
     * For q=1 (with-monitor): B is the 2ω/g-weighted anisotropic Laplacian
     *   surrogate — SPD, but a fixed-point approximation.
     * See doc/derivation_hessian.md [H-5].
     *
     * @param u      Current free-DoF control vector
     * @param B_out  Output SPD sparse matrix (nControls × nControls)
     */
    void assemblePicard(const gsVector<T> & u, gsSparseMatrix<T> & B_out) const;

    // ------------------------------------------------------------------
    // Solvers
    // ------------------------------------------------------------------

    /**
     * @brief Newton solver with modified-Newton inertia check and Armijo line search.
     *
     * Terminates when ‖R‖ < tol (option "Tolerance") or maxIter (option "MaxIter").
     * Logs energy, ‖R‖, min det J_σ, step α, wall time per iteration.
     * Returns number of iterations taken.
     */
    index_t solveNewton();

    /**
     * @brief Picard solver (frozen-coefficient linear solve per iterate).
     *
     * Same loop structure as solveNewton(), same logger, same line search.
     * Convergence is linear (quadratic only for Newton).
     */
    index_t solvePicard();

    // ------------------------------------------------------------------
    // Queries (delegated to the underlying gsOptMesh)
    // ------------------------------------------------------------------

    T evalObj() const;

    T computeMinJacobian() const;

    /// Iteration log: rows = [energy, ||R||, minDet, stepAlpha, wallTime]
    const gsMatrix<T> & iterLog() const { return m_iterLog; }

    // ------------------------------------------------------------------

private:
    // Shared solve loop: Newton (useHessian=true) or Picard (useHessian=false)
    index_t _solveLoop(bool useHessian);

    // Evaluate residual using gsOptMesh::gradObj_into
    void _residual(const gsVector<T> & u, gsVector<T> & R) const;

    // Fused single-sweep evaluation of the energy E(u) AND the minimum
    // sigma-Jacobian det over all quadrature points. Equivalent to
    // gsOptMesh::evalObj + gsOptMesh::computeMinJacobian, but in one OpenMP
    // element pass: det(J_sigma) is essentially free since J_sigma is already
    // formed for the energy integrand. The energy formula mirrors evalObj
    // exactly (validated to machine precision in unittests); the line search
    // uses this so each backtrack costs one sweep instead of two.
    void _evalObjAndMinJ(const gsVector<T> & u, T & energy, T & minJ) const;

    // FD Hessian of the analytic gradient (validation oracle + dd==3 fallback)
    void _hessianFD(const gsVector<T> & u, gsSparseMatrix<T> & K_out) const;

    // Apply Levenberg shift until K has correct inertia (all positive eigenvalues)
    void _makeDefinite(gsSparseMatrix<T> & K, T & tau) const;

    // Armijo back-tracking line search returning accepted step length. Each
    // trial step costs ONE fused sweep (_evalObjAndMinJ) yielding both the
    // energy and the fold guard. If \a nBacktracks is non-null, the number of
    // trial steps tried (one fused sweep each) is added to it.
    T _armijoStep(const gsVector<T> & u0, const gsVector<T> & R,
                  const gsVector<T> & direction, T f0,
                  index_t * nBacktracks = nullptr) const;

    // Helper: copy free DoFs between gsSquareDomain and gsVector
    void _getControls(gsVector<T> & u) const;
    void _setControls(const gsVector<T> & u) const;

    // ------------------------------------------------------------------
    // Hand-rolled assembly infrastructure (mirrors gsExprAssembler)
    // ------------------------------------------------------------------

    // Build (once) and cache the free-DoF sparsity pattern of K/B in
    // m_fpattern: all (ii,jj) pairs of free DoFs sharing an element, all
    // component pairs (d,m). During numeric assembly each thread only
    // accumulates into its own value buffer via _fpatternPos (no insertion,
    // no cross-thread writes), merged deterministically by _mergePartialK.
    void _buildPattern() const;

    // Flat position of the fixed pattern entry (ii,jj) within a per-thread
    // value buffer of size m_fpattern.nonZeros(), laid out column-by-column
    // (matching m_fpattern's ColMajor fiber order) with rows sorted within
    // each column. Valid only after _buildPattern() has inserted (ii,jj).
    // Returns -1 if (ii,jj) is not in the pattern (GISMO_ASSERT still fires
    // in Debug); callers running inside an OpenMP region must check this
    // instead of letting a missing entry silently overrun the value buffer,
    // since a GISMO_ENSURE cannot be thrown across a parallel region.
    index_t _fpatternPos(index_t ii, index_t jj) const;

    // Merge per-thread value buffers (indexed via _fpatternPos, thread-id
    // order) into m_fpattern column by column, then export to \a out.
    // Shared by assembleHessian and assemblePicard. Pointers rather than
    // values: the buffers are the persistent per-thread scratch (thPartialK),
    // not per-call copies.
    void _mergePartialK(const std::vector<const gsVector<T> *> & partial,
                         gsSparseMatrix<T> & out) const;

    // Third derivatives of the geometry S at (geometry-)parametric points.
    // out layout matches evalAllDers_into()[3]: row = nC3*a + c1 with
    // c1 = #ones among the derivative triple (dd==2: nC3 = 4).
    // Exact via evalAllDers_into(n=3) for spline geometries; central FD of
    // deriv2 (h=1e-5) for bases without third derivatives (e.g. NURBS).
    void _geomDeriv3_into(const gsMatrix<T> & pts, gsMatrix<T> & out) const;

    // Third derivatives of the monitor f at points (parametric or physical
    // according to m_parametric). Always central FD of deriv2 (generic
    // gsFunction does not provide n=3). Same layout as _geomDeriv3_into.
    void _funDeriv3_into(const gsMatrix<T> & pts, gsMatrix<T> & out) const;

private:
    gsSquareDomain<T>         * m_comp;
    const gsGeometry<T>       * m_geom;
    const gsFunction<T>       * m_fun;
    const gsBasis<T>          * m_ib;
    bool                        m_parametric;
    gsMultiBasis<T>             m_mb;

    // gsOptMesh instance: reused for evalObj, gradObj_into, computeMinJacobian
    mutable gsOptMesh<T,MODE>   m_optMesh;

    gsOptionList                m_options;

    /** @brief Per-thread scratch space for the _evalObjAndMinJ() element sweep.
     *   Sibling of gsOptMesh::EvalScratch. Rationale: \ref adaptparam_performance.
     * @warning Sized with omp_get_max_threads() at construction; raising
     *   the thread count after this object has been built is not supported.
     */
    struct EvalScratch
    {
        gsFuncData<T> compData, geomData, funData;
        gsMatrix<T>   Js, Jg, Jc, C, Cinv, Cg, Cg_inv, grad_xi_f;
        gsQuadRule<T> QuRule;
        index_t       QuPatch = -1;
        gsMatrix<T>   uvPoints, Jsigma_flat, Jgeom_flat, monVals, monDerivs;
        gsVector<T>   tmpWeights;
        bool          ready = false;

        void init(index_t dd, index_t td);
    };

    /// Per-thread scratch space for the assembleResidual() element sweep.
    /// See EvalScratch for the rationale.
    struct ResidualScratch
    {
        gsFuncData<T> geomData, funData, sigmaData, sigmaBasisData;
        gsQuadRule<T> QuRule;
        index_t       QuPatch = -1;
        gsMatrix<T>   uvPoints;
        gsVector<T>   tmpWeights;
        gsMatrix<T>   Jsigma_flat, Jgeom_flat, deriv2_geom;
        gsMatrix<T>   monVals, monDerivs, monDeriv2;
        gsMatrix<index_t> actives;
        gsMatrix<T>   basisVals, basisDerivs;
        gsMatrix<T>   Js, Jg, Jc, C, Q;
        gsMatrix<T>   D_d, Md, A_d, QA;
        gsMatrix<T>   adj, AdjJg;
        gsVector<T>   b_d, Qb, Q2b;
        gsMatrix<T>   Pi;
        gsVector<T>   r;
        gsVector<T>   gA, v;
        gsMatrix<T>   Cg, Cg_inv;
        gsMatrix<T>   gradMon, grad_xi_f, grad_x_f;
        gsMatrix<T>   Hess_f;
        gsVector<T>   mon_scalar_d;
        gsMatrix<T>   E_d, HfJgd, Dg_d, JtHJd;
        gsVector<T>   vs;
        gsVector<T>   dm2coef;
        bool          ready = false;

        void init(index_t dd, index_t td);
    };

    /// Per-thread scratch space for the assembleHessian() element sweep.
    /// See EvalScratch for the rationale. \a thPartialK accumulates this
    /// thread's contribution to the Newton tangent; it is re-zeroed at the
    /// start of every assembleHessian() call (it must not carry over between
    /// calls) but its Eigen buffer is kept allocated across calls.
    struct HessianScratch
    {
        gsFuncData<T> geomData, funData, sigmaData, sigmaBasisData;
        gsQuadRule<T> QuRule;
        index_t       QuPatch = -1;
        gsMatrix<T>   uvPoints;
        gsVector<T>   tmpWeights;
        gsMatrix<T>   Jsigma_flat, Jgeom_flat, deriv2_geom, deriv3_geom;
        gsMatrix<T>   monVals, monDerivs, monDeriv2, monDeriv3;
        gsMatrix<index_t> actives;
        gsMatrix<T>   basisVals, basisDerivs;
        gsMatrix<T>   Js, Jg, Jc, C, Q;
        gsMatrix<T>   Cg, Cg_inv;
        std::vector<gsMatrix<T> > D, Md, A, QAQ;
        std::vector<gsVector<T> > b, Qb, Q2b;
        gsVector<T>   trCAC, ad;
        gsMatrix<T>   AdjJg;
        std::vector<std::vector<gsMatrix<T> > > T3, E3;
        std::vector<std::vector<gsVector<T> > > Ddm, cdm, Q2c, edm, amg, adjD;
        gsMatrix<T>   M2, trQHQ, S1, Ss2, cJg, tT, bdMM;
        gsMatrix<T>   adj, adjMd, tmpM, Htmp, QQ;
        gsVector<T>   gA, gB, qA, qB, v, v_m;
        gsMatrix<T>   Hess_f, gradMon, grad_xi_f, grad_x_f;
        gsVector<T>   ms;
        gsMatrix<T>   dms, hEta;
        gsMatrix<T>   localK;
        std::vector<gsMatrix<T> > E1;
        std::vector<gsVector<T> > gfm;
        gsVector<T>   dgf;
        gsMatrix<T>   dHxf;
        gsVector<T>   trCginvEd_vec;
        gsMatrix<T>   trCginvEdm_mat, trCginvEmCginvEd_mat;
        gsVector<T>   dm2coef;
        gsMatrix<T>   d2m2coef;
        gsVector<T>   thPartialK;
        // Set true (per call, per thread) if a scatter target (ii,jj) is not
        // found in m_fpattern; checked and reported after the parallel region
        // (see the comment at the two scatter sites in .hpp).
        bool          patternMiss = false;
        bool          ready = false;

        void init(index_t dd, index_t td);
    };

    /// Per-thread scratch space for the assemblePicard() element sweep.
    /// See HessianScratch for the \a thPartialK re-zeroing rule and the
    /// \a patternMiss guard.
    struct PicardScratch
    {
        gsFuncData<T> geomData, funData, sigmaData, sigmaBasisData;
        gsQuadRule<T> QuRule;
        index_t       QuPatch = -1;
        gsMatrix<T>   uvPoints;
        gsVector<T>   tmpWeights;
        gsMatrix<T>   Jsigma_flat, Jgeom_flat;
        gsMatrix<T>   monVals, monDerivs;
        gsMatrix<index_t> actives;
        gsMatrix<T>   basisVals, basisDerivs;
        gsMatrix<T>   Js, Jg, Jc, C, Q;
        gsMatrix<T>   Cg, Cg_inv, Bc;
        gsMatrix<T>   grad_xi_f;
        gsVector<T>   gA, gB, BcgB;
        gsMatrix<T>   localB;
        gsVector<T>   thPartialK;
        bool          patternMiss = false;
        bool          ready = false;

        void init(index_t dd, index_t td);
    };

    mutable util::gsThreaded<EvalScratch>     m_evalScratch;
    mutable util::gsThreaded<ResidualScratch> m_residualScratch;
    mutable util::gsThreaded<HessianScratch>  m_hessianScratch;
    mutable util::gsThreaded<PicardScratch>   m_picardScratch;

    // Iteration log: columns = iterations, rows = [E, ||R||, minDet, alpha, time]
    mutable gsMatrix<T>         m_iterLog;

    // Cached sparsity pattern for matrix assembly (built once; the support
    // graph of the sigma basis is iteration-independent)
    mutable gsFiberMatrix<T,ColMajor> m_fpattern;
    mutable bool                m_patternReady;

    // Exclusive prefix sum of m_fpattern.nonZerosPerFiber(): column j's
    // entries occupy [m_fpatternColOffset[j], m_fpatternColOffset[j] +
    // nonZeros(column j)) in the flat value buffer used by _fpatternPos().
    // Cached alongside m_fpattern, valid whenever m_patternReady is true.
    mutable gsVector<index_t>   m_fpatternColOffset;

    // -1: unknown (probe on first use), 1: exact evalAllDers(3), 0: FD of deriv2
    mutable short_t             m_geomDeriv3Mode;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAdaptiveParametrizationNewton.hpp)
#endif
