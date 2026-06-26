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

enum MonitorMode
{
    // We can add the option "Runtime" where the mode is determined at runtime using an option
    ValueBased = 0,
    GradientBased = 1
};

namespace gismo
{

template <typename T,enum MonitorMode MODE>
class gsOptMesh : public gsOptProblem<T>
{
    using Base = gsOptProblem<T>;

private:
//     typedef typename gsExprAssembler<T>::space space;
//     typedef typename gsExprAssembler<T>::solution solution;

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

    gsSquareDomain<T> & composition();

    gsOptionList & options();

    /**
     * @brief Evaluates the objective function at the given control point \a u.
     *
     * The domain chain is \f$ \hat\Omega \xrightarrow{\sigma(\cdot;\alpha)}
     * \tilde\Omega \xrightarrow{G} \Omega \f$, with composed Jacobian
     * \f$ J_c = J_g\,J_\sigma \f$ and metric tensor
     * \f$ C = J_c^\top J_c \f$.
     *
     * **Without monitor** (\c m_fun == nullptr, \f$\omega=1/g\f$, paper Eq.(14)):
     * \f[
     *   E(\alpha) = \int_{\hat\Omega} \operatorname{tr}(C^{-1}) \; d\hat\Omega
     *             = \int_{\hat\Omega} \frac{T}{g^2} \; d\hat\Omega
     * \f]
     * This is the harmonic-map Dirichlet energy of the inverse parameterization
     * \f$G^{-1}:\Omega\to\hat\Omega\f$ (Radó–Kneser–Choquet bijectivity theory).
     *
     * **With monitor** (\c m_fun != nullptr), paper Eq.(18) [surface] / Eq.(17) [planar]:
     * \f[
     *   E(\alpha) = \int_{\hat\Omega} \omega(f)\,\frac{T}{g} \; d\hat\Omega
     * \f]
     * The discrete integrand uses \f$ \omega = \sqrt{m^2} \f$ (the paper's monitor
     * weight, NOT \f$m^2\f$): with the fold barrier \f$\chi\f$ and (penalty\f$\to0\f$)
     * \f$ |\det|^3/\chi^2 \to g_\sigma \f$ this reads
     * \f$ \omega\,\operatorname{tr}(C^{-1})\,|\det|^3/\chi^2 \to \omega\,T/g \f$.
     * The weight \f$\omega\f$ depends on the \c MonitorMode:
     * - **ValueBased** (paper Eq.13): \f$ \omega = 1/\sqrt{1+\theta f} \f$,
     *   i.e. \f$ m^2 = 1/(1+\theta f) \f$ (linear in \f$f\f$, NOT \f$f^2\f$).
     * - **GradientBased** (paper Eq.12, implemented form): \f$ \omega =
     *   1/\sqrt{1+\theta\,\eta^2} \f$ with \f$ \eta^2 = (\nabla_\xi f)^\top
     *   C_g^{-1}\,(\nabla_\xi f) = \|\nabla_x f\|^2 \f$ (squared physical gradient
     *   norm), \f$ C_g = J_g^\top J_g \f$.
     *
     * @note Surface case (td>dd): the integrand carries an extra surface area
     *   element \f$ g_S = \sqrt{\det C_g} \f$ so that
     *   \f$ \omega\,T/(g_\sigma g_S^2)\cdot g_S = \omega\,T/g_c \f$ matches Eq.(18)
     *   exactly (\f$g_c = g_\sigma g_S\f$ is the composite area element). Planar
     *   (td==dd) needs no such factor (\f$\det J_c\f$ already carries the full area).
     *
     * @note The \f$\det^2\f$ (no-monitor) vs \f$\det^3\f$ (with-monitor) exponent
     *   tracks MONITOR PRESENCE, not planar-vs-surface.
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
     * of the \a k-th basis function coefficient of \f$\sigma\f$).
     *
     * The key kinematic derivatives are:
     * \f{align*}{
     *   \frac{\partial\xi_i}{\partial\alpha_{k,d}} &= N_k\,\delta_{id}, \\
     *   \frac{\partial J_\sigma}{\partial\alpha_{k,d}}
     *     &= e_d\,(\nabla_{\hat u} N_k)^\top, \\
     *   \frac{\partial J_g}{\partial\alpha_{k,d}}(a,j)
     *     &= N_k \frac{\partial^2 G_a}{\partial\xi_j\,\partial\xi_d}.
     * \f}
     *
     * See gradObj_into() implementation for the full derivation per mode.
     */
    void gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result) const override;

    /**
     * @brief Finite-difference gradient (central differences, step \f$h=10^{-7}\f$)
     *   for validation of gradObj_into().
     */
    void gradObj_FD_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const;

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

      gsSquareDomain<T>          *   m_comp;
      const gsGeometry<T>       *   m_geom;
      const gsFunction<T>       *   m_fun;
      const gsBasis<T>          *   m_ib;
      gsMultiBasis<T>               m_mb;
      gsComposedGeometry<T>         m_cgeom;
      mutable gsComposedFunction<T> m_cfun; // MUTABLE BECAUSE OF TEMPORARY FIX
      gsMultiPatch<T>               m_mp;
      bool                          m_parametric;   // CAN BE REMOVED IF TEMPORARY FIX IS REMOVED


      // Controls of the composition
      // NOTE: Different from m_curDesign, since m_curDesign updates every time
      gsMatrix<T>                   m_controls;

      gsOptionList                  m_options;
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
     * This constructor takes any composition and integration basis among other inputs, and computes integrals based on the integration rule on the integration basis.
     *
     * @param composition       a \a gsFunction object representing the composition of the parametrization
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param integrationBasis  a \a gsBasis object used to define the integration points
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
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



      gsOptionList & options() override;

      void defaultOptions() override;

      void solve() override;

      T computeMinJacobian() const override;

public:

      template <short_t d>
      static gsTensorBSplineBasis<d,T> makeIntegrationBasis(const gsTensorBSplineBasis<d,T> & basis1,
                                                            const gsTensorBSplineBasis<d,T> & basis2);

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