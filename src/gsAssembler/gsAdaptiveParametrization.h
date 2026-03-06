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

    gsOptMesh(        gsFunction<T> & composition,
                const gsGeometry<T> & geometry,
                const gsBasis<T>    * integrationBasis);

    gsOptMesh(        gsFunction<T> & composition,
                const gsGeometry<T> & geometry,
                const gsFunction<T> * fun,
                const gsBasis<T>    * integrationBasis,
                const bool            parametric);

    gsFunction<T> & composition();

    gsOptionList & options();

    /// Evaluates the objective function at the given point u.
    T evalObj(const gsAsConstVector<T> &u) const;

    /// Gradient evaluation of the objective function at the given point u.
    void gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result) const override;

    /// Finite-difference gradient (central differences) for testing.
    void gradObj_FD_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const;

protected:

      // From gsOptProblem
      using Base::m_curDesign;
      using Base::m_numDesignVars;

      gsFunction<T>             *   m_comp;
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
    gsAdaptiveParametrization(        gsFunction<T>  & composition,
                                const gsGeometry<T>  & geometry,
                                const gsBasis<T>     & integrationBasis,
                                      gsOptimizer<T> & optimizer,
                                const bool             parametric=true);

    /**
     * @brief Constructs a gsAdaptiveParametrization object.
     *
     * This constructor takes any composition and integration basis among other inputs, and computes integrals based on the integration rule on the integration basis.
     *
     * @param composition       a \a gsFunction object representing the composition of the parametrization
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param function          a \a gsFunction object representing the indicator function
     * @param integrationBasis  a \a gsBasis object used to define the integration points
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
    gsAdaptiveParametrization(        gsFunction<T>  & composition,
                                const gsGeometry<T>  & geometry,
                                const gsFunction<T>  & function,
                                const gsBasis<T>     & integrationBasis,
                                      gsOptimizer<T> & optimizer,
                                const bool             parametric=true);

    /**
     * @brief Constructs a gsAdaptiveParametrization object.
     *
     * This constructor takes any composition among other inputs, and computes integrals based on the integration rule on the integration basis.
     *
     * @param composition       a \a gsFunction object representing the composition of the parametrization
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
    gsAdaptiveParametrization(        gsFunction<T>  & composition,
                                const gsGeometry<T>  & geometry,
                                      gsOptimizer<T> & optimizer,
                                const bool             parametric=true);

    /**
     * @brief Constructs a gsAdaptiveParametrization object.
     *
     * This constructor takes any composition among other inputs, and computes integrals based on the integration rule on the integration basis.
     *
     * @param composition       a \a gsFunction object representing the composition of the parametrization
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param function          a \a gsFunction object representing the indicator function
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
    gsAdaptiveParametrization(        gsFunction<T>  & composition,
                                const gsGeometry<T>  & geometry,
                                const gsFunction<T>  & function,
                                      gsOptimizer<T> & optimizer,
                                const bool             parametric=true);

    /**
     * @brief Constructs a gsAdaptiveParametrization object.
     *
     * This constructor takes any composition and integration basis among other inputs, and computes integrals based on the integration rule on the integration basis.
     *
     * @param composition       a \a gsFunction object representing the composition of the parametrization
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param function          a \a gsFunction object representing the indicator function
     * @param integrationBasis  a \a gsBasis object used to define the integration points
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
    gsAdaptiveParametrization(        gsFunction<T>  & composition,
                                const gsGeometry<T>  & geometry,
                                const gsFunction<T>  * function,
                                const gsBasis<T>     & integrationBasis,
                                      gsOptimizer<T> & optimizer,
                                const bool             parametric=true);

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

public:

      template <short_t d>
      static gsTensorBSplineBasis<d,T> makeIntegrationBasis(const gsTensorBSplineBasis<d,T> & basis1,
                                                            const gsTensorBSplineBasis<d,T> & basis2);

protected:

    gsFunction<T>             & m_comp;
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