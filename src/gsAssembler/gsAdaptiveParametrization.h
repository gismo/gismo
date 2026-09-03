/** @file gsAdaptiveParametrization.h

    @brief Provides class for adaptive parametrization (r-adaptivity).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-)
*/

#pragma once

#include <gsOptimizer/gsOptProblem.h>
#include <gsOptimizer/gsOptimizer.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsIO/gsOptionList.h>
#include <gsCore/gsComposedGeometry.h>
#include <gsUtils/gsThreaded.h>
#include <gsAssembler/gsQuadrature.h>
#include <gsAssembler/gsFoldBarrier.h>
#include <gsMatrix/gsSparseMatrix.h>
#include <gsHSplines/gsTHBSplineBasis.h>

enum MonitorMode
{
    // We can add the option "Runtime" where the mode is determined at runtime using an option
    ValueBased = 0,
    GradientBased = 1
};

namespace gismo
{

/**
 * @brief Base class for optimization problems over the controls of a
 *   \a gsSquareDomain "sigma" map.
 *
 * Owns the design-variable bookkeeping shared by all sigma optimizers
 * (composition(), fillCurDesign(), setCtrls(), m_domain,
 * m_curDesign/m_numDesignVars). Derivation: \ref adaptparam_overview.
 */
template <class T>
class gsOptSigma : public gsOptProblem<T>
{
    using Base = gsOptProblem<T>;

public:

    /// Default empty constructor (m_domain = nullptr)
    gsOptSigma();

    /// Constructs the design-variable bookkeeping for \a domain: sets
    /// m_numDesignVars and resizes m_curDesign (does NOT fill it -- see
    /// fillCurDesign()).
    explicit gsOptSigma(gsSquareDomain<T> & domain);

    virtual ~gsOptSigma() {}

    /// Returns the sigma domain this optimization problem controls.
    gsSquareDomain<T> & composition();

protected:

    /// Fill m_curDesign with the domain's current controls. gsOptMesh does
    /// NOT call this (it only resizes m_curDesign); gsOptFit/gsOptL2 do.
    void fillCurDesign();

    /// Pushes the design vector \a u back into the sigma domain's controls.
    void setCtrls(const gsAsConstVector<T> & u) const;

    // From gsOptProblem
    using Base::m_curDesign;
    using Base::m_numDesignVars;

    gsSquareDomain<T> * m_domain;
};

/// Central finite-difference validation of an objective's analytic gradient over
/// the free controls of \a domain. Diagnostic only; returns the maximum relative
/// error, skipping components where both the gradient and the difference quotient
/// are below \a floor. Restores \a domain's controls before returning.
template <class T>
T gsCheckSigmaGradient(gsOptProblem<T> & problem, gsSquareDomain<T> & domain,
                       T h = (T)1e-6, T floor = (T)1e-10);

/**
 * @brief Winslow/monitor-based mesh-relocation objective (R step) over the
 *   controls of a \a gsSquareDomain map \f$\sigma\f$ -- see evalObj() for the
 *   energy.
 *
 * Exposes "Smoothing", "Penalty", "quA", "quB" (the element-sweep quadrature)
 * as options. Carries no additive \a gsFoldBarrier (unlike \a gsOptFit /
 * \a gsOptL2, the D step). Derivation: \ref adaptparam_rstep.
 *
 * \pre "Penalty" > 0 (fold-barrier regularisation).
 * \pre "Smoothing" (theta) >= 0.
 * \pre With a ValueBased monitor, 1+theta*f > 0 at every quadrature point.
 */
template <typename T,enum MonitorMode MODE>
class gsOptMesh : public gsOptSigma<T>
{
    using Base = gsOptSigma<T>;

public:

    gsOptMesh();

    gsOptMesh(        gsSquareDomain<T> & composition,
                const gsGeometry<T> & geometry,
                const gsBasis<T>    * integrationBasis);

    gsOptMesh(        gsSquareDomain<T> & composition,
                const gsGeometry<T> & geometry,
                const gsFunction<T> * fun,
                const gsBasis<T>    * integrationBasis,
                const bool            parametric);

    gsOptionList & options();

    /**
     * @brief Evaluates the objective function at the given control point \a u.
     *
     * Without a monitor function, this is the harmonic-map Dirichlet energy
     * \f$ E(\alpha) = \int_{\hat\Omega} \operatorname{tr}(C^{-1}) \, d\hat\Omega
     * = \int_{\hat\Omega} T/g^2 \, d\hat\Omega \f$ of the inverse
     * parameterization (Radó–Kneser–Choquet bijectivity theory), with
     * composed Jacobian \f$J_c=J_g J_\sigma\f$ and \f$C=J_c^\top J_c\f$. With
     * a monitor function, \f$ E(\alpha) = \int_{\hat\Omega} \omega(f)\,T/g \,
     * d\hat\Omega \f$, \f$\omega\f$ depending on the \c MonitorMode.
     * Derivation (both energies, the per-mode weight, and the surface-vs-planar
     * area element): \ref adaptparam_rstep.
     *
     * @note For GradientBased mode, the derivatives of the monitor
     *   function \a m_fun must be exact (e.g. via automatic
     *   differentiation). Finite-difference derivatives introduce
     *   errors that propagate into the gradient computation.
     */
    T evalObj(const gsAsConstVector<T> &u) const;

    /**
     * @brief Analytical gradient of the objective function w.r.t.
     *   the control variables \a u.
     *
     * Differentiates the integrand of evalObj() with respect to each
     * control variable \f$ \alpha_{k,d} \f$ (the \a d-th coordinate
     * of the \a k-th basis function coefficient of \f$\sigma\f$), using
     * \f{align*}{
     *   \frac{\partial\xi_i}{\partial\alpha_{k,d}} &= N_k\,\delta_{id}, \\
     *   \frac{\partial J_\sigma}{\partial\alpha_{k,d}}
     *     &= e_d\,(\nabla_{\hat u} N_k)^\top, \\
     *   \frac{\partial J_g}{\partial\alpha_{k,d}}(a,j)
     *     &= N_k \frac{\partial^2 G_a}{\partial\xi_j\,\partial\xi_d}.
     * \f}
     * Full per-mode derivation: \ref adaptparam_rgrad.
     */
    void gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result) const override;

    /**
     * @brief Finite-difference gradient (central differences) for validation
     *   of gradObj_into(). Per component, the step is \a h scaled by
     *   \c max(1,|u_i|), matching gsCheckSigmaGradient()'s convention.
     * @param h base step size, same default as gsCheckSigmaGradient().
     */
    void gradObj_FD_into( const gsAsConstVector<T> & u, gsAsVector<T> & result,
                          T h = (T)1e-6) const;

    /**
     * @brief Returns the minimum \f$\det(J_\sigma)\f$ over all quadrature points.
     *
     * A negative return value means at least one element is folded.
     * Useful for post-optimization fold detection.
     */
    T computeMinJacobian(const gsAsConstVector<T> & u) const;

protected:

      // From gsOptProblem
      using Base::m_curDesign;
      using Base::m_numDesignVars;
      // From gsOptSigma
      using Base::m_domain;

      const gsGeometry<T>       *   m_geom;
      const gsFunction<T>       *   m_fun;
      const gsBasis<T>          *   m_ib;
      gsMultiBasis<T>               m_mb;
      gsComposedGeometry<T>         m_cgeom;
      gsMultiPatch<T>               m_mp;
      bool                          m_parametric;   ///< true: the monitor is evaluated at the parametric point xi = sigma(u); false: at the physical point x = G(xi).


      // Controls of the composition
      // NOTE: Different from m_curDesign, since m_curDesign updates every time
      gsMatrix<T>                   m_controls;

      gsOptionList                  m_options;

      /** @brief Per-thread scratch space for the evalObj() element sweep;
       *   allocated once and reused for the whole solve. Rationale:
       *   \ref adaptparam_performance.
       * @warning Sized with omp_get_max_threads() at construction; raising
       *   the thread count after a gsOptMesh has been built is not supported.
       */
      struct EvalScratch
      {
          gsFuncData<T> compData, geomData, funData;
          gsMatrix<T>   Js, Jg, Jc, Cg, Cg_inv, grad_xi_f;
          gsQuadRule<T> QuRule;
          index_t       QuPatch = -1;
          gsMatrix<T>   uvPoints;
          gsVector<T>   tmpWeights;
          gsMatrix<T>   monVals, monDerivs_eval;
          gsMatrix<T>   Jgeom_flat, Jsigma_flat;
          bool          ready = false;

          void init(index_t dd, index_t td);
      };

      /// Per-thread scratch space for the gradObj_into() element sweep.
      /// See EvalScratch for the rationale.
      struct GradScratch
      {
          gsFuncData<T> geomData, funData, sigmaData, sigmaBasisData;
          gsMatrix<T>   Js, Jg, Jc, JcT, Cg, Cg_inv;
          gsMatrix<T>   gradMon, grad_xi_f, grad_x_f, Hess_f;
          gsVector<T>   gradNk, v;
          gsMatrix<T>   D_d, DdJs;
          gsVector<T>   b_d, trA;
          gsMatrix<T>   b_all;
          gsVector<T>   mon_scalar_d;
          gsVector<T>   trAdjJcDdJs, trCginvE_d;
          gsMatrix<T>   adjJcT_precomp, adjJc, E_d;
          gsMatrix<T>   HfJgd, Dg_d, JtHJd;
          gsQuadRule<T> QuRule;
          index_t       QuPatch = -1;
          gsMatrix<T>   uvPoints;
          gsVector<T>   tmpWeights;
          gsMatrix<T>   monVals, monDerivs, monDeriv2;
          gsMatrix<T>   Jsigma_flat, Jgeom_flat, deriv2_geom;
          gsMatrix<index_t> actives;
          gsMatrix<T>   basisVals, basisDerivs;
          gsVector<T>   thResult;
          bool          ready = false;

          void init(index_t dd, index_t td, bool parametric, index_t nc);
      };

      mutable util::gsThreaded<EvalScratch> m_evalScratch;
      mutable util::gsThreaded<GradScratch> m_gradScratch;
};

/**
 * @brief Direct point-cloud fitting objective over the controls of a
 *   \a gsSquareDomain map \f$\sigma\f$, for a composed fit \f$G=S\circ\sigma\f$
 *   (the D step). Objective and gradient: \ref adaptparam_dstep.
 *
 * @note This class is written for a planar (domainDim==targetDim==2) sigma
 *   domain; the constructor asserts it.
 */
template <class T>
class gsOptFit : public gsOptSigma<T>
{
    using Base = gsOptSigma<T>;

public:

    /// Constructs the fitting problem. \a mu, \a eps, \a mode, \a quB are
    /// forwarded verbatim to the \a gsFoldBarrier constructor (see there for
    /// their meaning). \a uv are the parametric fitting points (on sigma's
    /// reference square) and \a xyz the corresponding target points in the
    /// physical space of \a S.
    /// @param quB Barrier quadrature order for \a gsFoldBarrierMode::Sampled
    ///   (ignored in \a gsFoldBarrierMode::Coefficient); must be >= 0 (a
    ///   negative value segfaults inside gsQuadrature::get() and is rejected
    ///   with GISMO_ENSURE here).
    gsOptFit(      gsSquareDomain<T> & domain,
             const gsGeometry<T>     & S,
             const gsMatrix<T>       & uv,
             const gsMatrix<T>       & xyz,
                   T mu, T eps, gsFoldBarrierMode mode, index_t quB);

    /// Evaluates E(u) (data term + barrier/box, see class doc).
    T evalObj(const gsAsConstVector<T> & u) const override;

    /// Analytic gradient of evalObj() (data term + barrier/box).
    void gradObj_into(const gsAsConstVector<T> & u, gsAsVector<T> & result) const override;

    /// Number of fold-barrier evaluation points (see gsFoldBarrier::nPoints()).
    index_t nBarrierPts() const { return m_barrier.nPoints(); }

protected:

    using Base::m_domain;

    const gsGeometry<T> * m_S;
    gsMatrix<T> m_uv, m_xyz;

    gsFoldBarrier<T> m_barrier;

    /// Collocation matrix of sigma's basis at m_uv (N x nb), fixed for the
    /// whole optimization -- see the class doc for why this makes the
    /// gradient assembly cheap.
    gsSparseMatrix<T> m_colloc;

    /// Cached basis data of sigma at m_uv (active function ids + basis
    /// values), reused every iteration via linearCombination_into with the
    /// current control coefficients.
    gsMatrix<index_t> m_uvActives;
    gsMatrix<T>       m_uvVals;
};

/**
 * @brief L2-projection-error objective over the controls of a
 *   \a gsSquareDomain map \f$\sigma\f$, with the analysis-space (projection)
 *   coefficients FROZEN (the D step). Objective and gradient derivation:
 *   \ref adaptparam_dstep.
 *
 * @warning Had the analysis/projection basis instead been COMPOSED with
 *   \f$\sigma\f$ too, a change of variables would make this frozen-coefficient
 *   error \f$\sigma\f$-INDEPENDENT and the gradient identically zero -- a
 *   silent no-op. \a gsOptL2 exists to avoid exactly that trap: its
 *   integration basis and \a u_h must both be expressed on \f$\hat\Omega\f$
 *   only.
 */
template <class T>
class gsOptL2 : public gsOptSigma<T>
{
    using Base = gsOptSigma<T>;

public:

    /// Constructs the L2-projection-error problem. \a solution is the frozen
    /// discrete solution \f$u_h\f$, a function of the reference variable
    /// \a v on \f$[0,1]^2\f$; \a fun is the target \a f. \a parametric
    /// follows gsOptMesh's convention: \c true means \a fun is a function of
    /// \f$\xi=\sigma(v)\f$, \c false (default) means it is a function of the
    /// physical point \f$G(v)\f$. \a mu, \a eps, \a mode, \a quB are
    /// forwarded verbatim to the \a gsFoldBarrier constructor; \a quA/\a quB
    /// of the ELEMENT SWEEP (distinct from the barrier's \a quB) are exposed
    /// as the "quA"/"quB" options, mirroring gsOptMesh's constructor.
    /// @param quB Barrier quadrature order for \a gsFoldBarrierMode::Sampled
    ///   (ignored in \a gsFoldBarrierMode::Coefficient); must be >= 0 (a
    ///   negative value segfaults inside gsQuadrature::get() and is rejected
    ///   with GISMO_ENSURE here).
    /// @note \a domain and \a geometry must both have \c domainDim()==2
    ///   (sigma is the unit square); \a geometry may be planar
    ///   (\c targetDim()==2) or a surface (\c targetDim()>2). \a fun's own
    ///   \c domainDim() must match how it is evaluated: \c geometry's
    ///   \c targetDim() when \c parametric==false (\a fun is evaluated at
    ///   the physical point \f$G(v)\f$), or 2 when \c parametric==true
    ///   (\a fun is evaluated at \f$\xi=\sigma(v)\f$). A mismatch here is
    ///   otherwise a silent out-of-bounds read inside \c gradObj_into's
    ///   chain-rule term.
    gsOptL2(      gsSquareDomain<T> & domain,      // sigma
            const gsGeometry<T>     & geometry,    // S, fixed
            const gsFunction<T>     & solution,    // frozen u_h, function of v on [0,1]^2
            const gsFunction<T>     & fun,         // target f
            const gsBasis<T>        * integrationBasis,
                  T mu, T eps, gsFoldBarrierMode mode, index_t quB,
            const bool                parametric = false);

    gsOptionList & options() { return m_options; }

    /// Evaluates E(u) (L2 error term + barrier/box, see class doc). The
    /// element sweep is OpenMP-parallel; per-thread partials are summed in
    /// thread-id order, so the value is reproducible run-to-run for a given
    /// thread count.
    T evalObj(const gsAsConstVector<T> & u) const override;

    /// Analytic gradient of evalObj() (see class doc for the derivation).
    void gradObj_into(const gsAsConstVector<T> & u, gsAsVector<T> & result) const override;

    /// Number of fold-barrier evaluation points (see gsFoldBarrier::nPoints()).
    index_t nBarrierPts() const { return m_barrier.nPoints(); }

protected:

    using Base::m_domain;

    const gsGeometry<T> * m_geom;     // S, fixed
    const gsFunction<T> * m_solution; // frozen u_h(v)
    const gsFunction<T> * m_fun;      // target f
    const gsBasis<T>    * m_ib;
    gsMultiBasis<T>        m_mb;
    bool                   m_parametric;

    gsFoldBarrier<T> m_barrier;

    gsOptionList m_options; // element-sweep quA/quB (distinct from the barrier's)

    /** @brief Per-thread scratch space for the evalObj() element sweep.
     *   Sibling of gsOptMesh::EvalScratch. Rationale: \ref adaptparam_performance.
     * @warning Sized with omp_get_max_threads() at construction; raising
     *   the thread count after a gsOptL2 has been built is not supported.
     */
    struct EvalScratch
    {
        gsFuncData<T> compData, geomData, funData;
        gsMatrix<T>   solVals;
        gsMatrix<T>   Js, JS, Cg;
        gsQuadRule<T> QuRule;
        index_t       QuPatch = -1;
        gsMatrix<T>   uvPoints;
        gsVector<T>   tmpWeights;
        bool          ready = false;

        void init(index_t dd, index_t tdS);
    };

    /// Per-thread scratch space for the gradObj_into() element sweep.
    /// See EvalScratch for the rationale.
    struct GradScratch
    {
        gsFuncData<T> compData, geomData, sigmaBasisData, funData;
        gsMatrix<T>   solVals, dDetJs;
        gsMatrix<T>   Js, JS, Cg, Cg_inv, adjJS, D_d, E_d;
        gsVector<T>   r;
        gsMatrix<T>   dFdxi;
        gsQuadRule<T> QuRule;
        index_t       QuPatch = -1;
        gsMatrix<T>   uvPoints;
        gsVector<T>   tmpWeights;
        gsVector<T>   thResult;
        bool          ready = false;

        void init(index_t dd, index_t tdS, index_t td_sol, index_t nc);
    };

    mutable util::gsThreaded<EvalScratch> m_evalScratch;
    mutable util::gsThreaded<GradScratch> m_gradScratch;
};

/**
 * @brief The gsAdaptiveParametrizationBase class is independent of the monitor mode
 */
template<class T>
class gsAdaptiveParametrizationBase
{

public:

    /// Default empty constructor
    gsAdaptiveParametrizationBase() {};

    /// Default deconstructor
    virtual ~gsAdaptiveParametrizationBase() {};

public:
      virtual gsOptionList & options() = 0;

      virtual void defaultOptions() = 0;

      virtual void solve() = 0;

      /**
       * @brief Returns the minimum \f$\det(J_\sigma)\f$ over all quadrature points.
       *
       * A negative return value indicates at least one folded element.
       */
      virtual T computeMinJacobian() const = 0;

};


/// Tag: the supplied basis is ALREADY the integration mesh (union + degree
/// raise already applied by the caller, e.g. via makeIntegrationBasis).
/// The constructor stores it verbatim instead of building the union itself -
/// passing a pre-built basis to the ordinary constructors would raise the
/// degree a second time.
struct integrationBasisIsFinal_t {};
static constexpr integrationBasisIsFinal_t integrationBasisIsFinal {};

/**
 * @brief The gsAdaptiveParametrization class is a template class that can be used to perform adaptive parametrization.
 *
 * @tparam T the type of the scalar values
 * @tparam MODE the monitor mode (ValueBased or GradientBased)
 */
template<class T, enum MonitorMode MODE=MonitorMode::ValueBased>
class gsAdaptiveParametrization : public gsAdaptiveParametrizationBase<T>
{
protected:

public:

    /**
     * @brief Constructs a gsAdaptiveParametrization object.
     *
     * This constructor takes a \ref gsSquareDomain as composition and integration basis among other inputs, and computes integrals based on the union of the integration basis and the composition
     *
     * @param composition       a \a gsSquareDomain object representing the composition of the parametrization. This object will be exactly integrated by inserting its integration points in the integration basis.
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param integrationBasis  a \a gsBasis object used to define the integration points. The integration points from the composition will be added to this basis.
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
    gsAdaptiveParametrization(        gsSquareDomain<T>     & composition,
                                const gsGeometry<T>         & geometry,
                                const gsBasis<T>            & integrationBasis,
                                      gsOptimizer<T>        & optimizer,
                                const bool                    parametric=true);


    /**
     * @brief Constructs a gsAdaptiveParametrization object.
     *
     * This constructor takes a \ref gsSquareDomain as composition and integration basis among other inputs, and computes integrals based on the union of the integration basis and the composition
     *
     * @param composition       a \a gsSquareDomain object representing the composition of the parametrization. This object will be exactly integrated by inserting its integration points in the integration basis.
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param function          a \a gsFunction object representing the indicator function
     * @param integrationBasis  a \a gsBasis object used to define the integration points. The integration points from the composition will be added to this basis.
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
    gsAdaptiveParametrization(        gsSquareDomain<T>     & composition,
                                const gsGeometry<T>         & geometry,
                                const gsFunction<T>         & function,
                                const gsBasis<T>            & integrationBasis,
                                      gsOptimizer<T>        & optimizer,
                                const bool                    parametric=true);

    /**
     * @brief Constructs a gsAdaptiveParametrization object.
     *
     * This constructor takes a \ref gsSquareDomain as composition among other inputs, and computes integrals based on the union of the integration basis and the composition
     * The integration basis is set to the geometry's basis.
     *
     * @param composition       a \a gsSquareDomain object representing the composition of the parametrization. This object will be exactly integrated by inserting its integration points in the integration basis.
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
    gsAdaptiveParametrization(        gsSquareDomain<T>     & composition,
                                const gsGeometry<T>         & geometry,
                                      gsOptimizer<T>        & optimizer,
                                const bool                    parametric=true);

    /**
     * @brief Constructs a gsAdaptiveParametrization object.
     *
     * This constructor takes a \ref gsSquareDomain as composition among other inputs, and computes integrals based on the union of the integration basis and the composition
     * The integration basis is set to the geometry's basis.
     *
     * @param composition       a \a gsSquareDomain object representing the composition of the parametrization. This object will be exactly integrated by inserting its integration points in the integration basis.
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param function          a \a gsFunction object representing the indicator function
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
    gsAdaptiveParametrization(        gsSquareDomain<T>     & composition,
                                const gsGeometry<T>         & geometry,
                                const gsFunction<T>         & function,
                                      gsOptimizer<T>        & optimizer,
                                const bool                    parametric=true);

/**
     * @brief Constructs a gsAdaptiveParametrization object.
     *
     * This constructor takes a \ref gsSquareDomain as composition and integration basis among other inputs, and computes integrals based on the union of the integration basis and the composition
     *
     * @param composition       a \a gsSquareDomain object representing the composition of the parametrization. This object will be exactly integrated by inserting its integration points in the integration basis.
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param function          a \a gsFunction object representing the indicator function
     * @param integrationBasis  a \a gsBasis object used to define the integration points. The integration points from the composition will be added to this basis.
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
    gsAdaptiveParametrization(        gsSquareDomain<T>     & composition,
                                const gsGeometry<T>         & geometry,
                                const gsFunction<T>         * function,
                                const gsBasis<T>            & integrationBasis,
                                      gsOptimizer<T>        & optimizer,
                                const bool                    parametric=true);

    /**
     * @brief Constructs a gsAdaptiveParametrization object from an
     *   ALREADY-BUILT integration basis. Stores \a integrationBasis verbatim
     *   instead of building the union mesh the other constructors build via
     *   makeIntegrationBasis() -- for a caller that already built the union
     *   itself (e.g. via the THB overload of makeIntegrationBasis()), which
     *   would otherwise have the degree raise applied a second time.
     *   Validation and admissibility rules: \ref adaptparam_supermesh.
     *
     * @param composition       a \a gsSquareDomain object representing the composition of the parametrization.
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param function          a \a gsFunction object representing the indicator function (may be \c nullptr)
     * @param integrationBasis  the FINAL integration basis, stored as-is (no union, no degree raise)
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain
     */
    gsAdaptiveParametrization(        gsSquareDomain<T>     & composition,
                                const gsGeometry<T>         & geometry,
                                const gsFunction<T>         * function,
                                const gsBasis<T>            & integrationBasis,
                                      gsOptimizer<T>        & optimizer,
                                const bool                    parametric,
                                      integrationBasisIsFinal_t);



      gsOptionList & options() override;

      void defaultOptions() override;

      void solve() override;

      T computeMinJacobian() const override;

public:

      template <short_t d>
      static gsTensorBSplineBasis<d,T> makeIntegrationBasis(const gsTensorBSplineBasis<d,T> & basis1,
                                                            const gsTensorBSplineBasis<d,T> & basis2);

      /// Builds the super mesh on the analysis THB's ACTUAL hierarchical
      /// element partition. \a basis2's knot mesh MUST be exactly a dyadic
      /// level of \a basis1's hierarchy (see sigmaLevelInHierarchy()) --
      /// otherwise \c GISMO_ENSURE throws; use the tensor overload above for
      /// a non-nested pair. Construction and preconditions: \ref adaptparam_supermesh.
      template <short_t d>
      static typename gsBasis<T>::uPtr
      makeIntegrationBasis(const gsTHBSplineBasis<d,T>     & basis1,
                           const gsTensorBSplineBasis<d,T> & basis2);

      /// True iff sigma's knot mesh (\a comp) is exactly a refinement LEVEL of
      /// the dyadic hierarchy of \a analysis; then \a sigmaLevel is set to that
      /// level, else \a reason receives a human-readable failure cause.
      /// Preconditions and the pass-through-constructor relationship:
      /// \ref adaptparam_supermesh.
      template <short_t d>
      static bool sigmaLevelInHierarchy(const gsTHBSplineBasis<d,T>      & analysis,
                                        const gsTensorBSplineBasis<d,T>  & comp,
                                        index_t                          & sigmaLevel,
                                        std::string                      & reason);

      /// True iff \a candidate is admissible as an integration basis for the
      /// composition basis \a comp, per direction: same parameter range,
      /// degree at least \a comp's, and every interior knot of \a comp
      /// present up to tolerance. Used by the pass-through
      /// (integrationBasisIsFinal_t) constructor. Coverage note: \ref adaptparam_supermesh.
      template <short_t d>
      static bool tensorBasisAdmissible(const gsTensorBSplineBasis<d,T>  & comp,
                                        const gsTensorBSplineBasis<d,T>  & candidate);

protected:

    gsSquareDomain<T>          & m_comp;
    const gsGeometry<T>       & m_geom;
    const gsFunction<T>       * m_fun;

    gsOptMesh<T,MODE>           m_optProblem;
    gsOptimizer<T>            & m_optimizer;
    typename gsBasis<T>::uPtr   m_integrationBasis;
    gsOptionList                m_options;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAdaptiveParametrization.hpp)
#endif