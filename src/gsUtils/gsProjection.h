/** @file gsProjection.h

    @brief Class that performs a projection

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/


#pragma once

#include <gsAssembler/gsExprAssembler.h>

namespace gismo {

/**
 * \brief Enumeration of projection norms

    This enumeration defines the projection norm identifiers exposed by gsProjection.
    The choice of norm affects the definition of the projection and the resulting
    coefficients.

    Currently implemented in gsProjection:
    - L2: Standard L2 projection, which minimizes the L2 norm of the error
      between the source function and its projection.
    - H1: H1 projection, which minimizes the H1 norm of the error, taking
      into account both the function values and their gradients.
    - H2: H2 projection, which takes into account function values and higher
      derivatives as supported by gsProjection.

    The choice of norm should be based on the specific requirements of the
    problem being solved and on the set of projection norms currently
    implemented by gsProjection.
 */
enum ProjectionNorm
{
    L2,
    H1,
    H2,
};

/** \brief Class that performs a projection

    \tparam T coefficient type
 */

template <enum ProjectionNorm Norm, class T>
struct gsProjection
{

protected:
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    typedef gsExprAssembler<>::element     element;

    /**
     * \brief Projects a source function onto a projection basis using a geometry map.
     *
     * This function computes the coefficients of the projection of a given source function onto a projection basis.
     * The projection is performed using a geometry map, which maps the integration domain to the physical domain.
     *
     * \param integrationBasis The basis used for numerical integration.
     * \param projectionBasis The basis functions used for the projection.
     * \param geometryMap The geometry map that maps the integration domain to the physical domain.
     * \param sourceFunction The source function to be projected.
     * \param coefs The output matrix that stores the computed coefficients of the projection.
     * \param options The options that control the projection process.
     *
     * \return the projection error.
     */
    static T _project(  const gsMultiBasis<T>         & integrationBasis,
                        const gsFunctionSet<T>        & projectionBasis,
                        const gsFunctionSet<T>        & geometryMap,
                        const gsFunctionSet<T>        & sourceFunction,
                              gsMatrix<T>             & coefs,
                        const gsBoundaryConditions<T> & bc,
                        const gsOptionList            & options,
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0);

	/** 
	 * \brief      Obtain the system matrix and right-hand side for the projection of a function onto a basis
	 * 
	 * @param[in]  integrationBasis The basis used for numerical integration.
	 * @param[in]  projectionBasis The basis functions used for the projection.
	 * @param[in]  geometryMap The geometry map that maps the integration domain to the physical domain.
	 * @param[in]  sourceFunction The source function to be projected.
	 * @param      systemMatrix The output system matrix of the projection problem.
	 * @param      rhs The output right-hand side vector of the projection problem.
	 * @param      options The options that control the projection process.
	 */
    static void _system(const gsMultiBasis<T>         & integrationBasis,
                        const gsFunctionSet<T>        & projectionBasis,
                        const gsFunctionSet<T>        & geometryMap,
                        const gsFunctionSet<T>        & sourceFunction,
                              gsSparseMatrix<T>       & systemMatrix,
                              gsMatrix<T>             & rhs,
                        const gsBoundaryConditions<T> & bc,
                        const gsOptionList            & options,
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0);

	/** 
	 * \brief      Obtain the system matrix for the projection of a function onto a basis
	 * 
	 * @param[in]  integrationBasis The basis used for numerical integration.
	 * @param[in]  projectionBasis The basis functions used for the projection.
	 * @param[in]  geometryMap The geometry map that maps the integration domain to the physical domain.
	 * @param      systemMatrix The output system matrix of the projection problem.
	 * @param      targetDim The target dimension of the projection.
	 * @param      options The options that control the projection process.
	 */
    static void _matrix(const gsMultiBasis<T>         & integrationBasis,
                        const gsFunctionSet<T>        & projectionBasis,
                        const gsFunctionSet<T>        & geometryMap,
                              gsSparseMatrix<T>       & systemMatrix,
                              short_t                   targetDim,
                        const gsBoundaryConditions<T> & bc,
                        const gsOptionList            & options,
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0);

	/**
	 * \brief      Obtain the right-hand side for the projection of a function onto a basis
	 * 
	 * @param[in]  integrationBasis The basis used for numerical integration.
	 * @param[in]  projectionBasis The basis functions used for the projection.
	 * @param[in]  geometryMap The geometry map that maps the integration domain to the physical domain.
	 * @param[in]  sourceFunction The source function to be projected.
	 * @param      rhs The output right-hand side vector of the projection problem.
	 * @param      options The options that control the projection process.
	 */
    static void _rhs   (const gsMultiBasis<T>         & integrationBasis,
                        const gsFunctionSet<T>        & projectionBasis,
                        const gsFunctionSet<T>        & geometryMap,
                        const gsFunctionSet<T>        & sourceFunction,
                              gsMatrix<T>             & rhs,
                        const gsBoundaryConditions<T> & bc,
                        const gsOptionList            & options,
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0);

	template<enum ProjectionNorm _Norm>
	static typename std::enable_if<(_Norm == L2), void>::type
	_assembleMatrix(gsExprAssembler<T> & A, space & u, geometryMap & G, T alpha, T beta, T gamma)
	{
		GISMO_UNUSED(beta);
		GISMO_UNUSED(gamma);
		A.assemble(alpha * (u * u.tr()) * meas(G));
	}

	template<enum ProjectionNorm _Norm>
	static typename std::enable_if<(_Norm == H1), void>::type
	_assembleMatrix(gsExprAssembler<T> & A, space & u, geometryMap & G, T alpha, T beta, T gamma)
	{
		GISMO_UNUSED(gamma);
		A.assemble((alpha * u * u.tr() + beta * igrad(u,G) * igrad(u,G).tr()) * meas(G));
	}

	template<enum ProjectionNorm _Norm>
	static typename std::enable_if<(_Norm == H2), void>::type
	_assembleMatrix(gsExprAssembler<T> & A, space & u, geometryMap & G, T alpha, T beta, T gamma)
	{
		A.assemble((alpha * u * u.tr() + beta * igrad(u,G) * igrad(u,G).tr() + gamma * ilapl(u,G) * ilapl(u,G).tr()) * meas(G));
	}

	template<enum ProjectionNorm _Norm, typename _F>
	static typename std::enable_if<(_Norm == L2), void>::type
	_assembleRhs(gsExprAssembler<T> & A, space & u, _F & f, geometryMap & G, T alpha, T beta, T gamma)
	{
		GISMO_UNUSED(beta);
		GISMO_UNUSED(gamma);
		A.assemble(alpha * (u * f) * meas(G));
	}

	template<enum ProjectionNorm _Norm, typename _F>
	static typename std::enable_if<(_Norm == H1), void>::type
	_assembleRhs(gsExprAssembler<T> & A, space & u, _F & f, geometryMap & G, T alpha, T beta, T gamma)
	{
		GISMO_UNUSED(gamma);
		A.assemble((alpha * u * f + beta * igrad(u,G) * igrad(f,G).tr()) * meas(G));
	}

	template<enum ProjectionNorm _Norm, typename _F>
	static typename std::enable_if<(_Norm == H2), void>::type
	_assembleRhs(gsExprAssembler<T> & A, space & u, _F & f, geometryMap & G, T alpha, T beta, T gamma)
	{
		A.assemble((alpha * u * f + beta * igrad(u,G) * igrad(f,G).tr() + gamma * ilapl(u,G) * ilapl(f,G)) * meas(G));
	}

	template<enum ProjectionNorm _Norm, typename _F>
	static typename std::enable_if<(_Norm == L2), T>::type
	_computeError(gsExprEvaluator<T> & ev, solution & s, _F & f, geometryMap & G, T alpha, T beta, T gamma)
	{
		GISMO_UNUSED(beta);
		GISMO_UNUSED(gamma);
		return ev.integral(alpha * ((s-f).sqNorm()) * meas(G));
	}

	template<enum ProjectionNorm _Norm, typename _F>
	static typename std::enable_if<(_Norm == H1), T>::type
	_computeError(gsExprEvaluator<T> & ev, solution & s, _F & f, geometryMap & G, T alpha, T beta, T gamma)
	{
		GISMO_UNUSED(gamma);
		return ev.integral((alpha * (s-f).sqNorm() + beta * (igrad(s,G)-igrad(f,G)).sqNorm()) * meas(G));
	}

	template<enum ProjectionNorm _Norm, typename _F>
	static typename std::enable_if<(_Norm == H2), T>::type
	_computeError(gsExprEvaluator<T> & ev, solution & s, _F & f, geometryMap & G, T alpha, T beta, T gamma)
	{
		return ev.integral((alpha * (s-f).sqNorm() + beta * (igrad(s,G)-igrad(f,G)).sqNorm() + gamma * (ilapl(s,G)-ilapl(f,G)).sqNorm()) * meas(G));
	}

public:

    /**
     * @brief      Project a geometry onto a basis (multi-patch)
     *
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param      coefs            The coefficients of the new geometry on \a projectionBasis
     *
     * @return     The L2 error of the projection
     */
    static T project(   const gsMultiBasis<T>         & projectionBasis,
                        const gsMultiPatch<T>         & geometryMap,
                              gsMatrix<T>             & coefs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        return _project(projectionBasis, projectionBasis, geometryMap, geometryMap, coefs, bc, options, alpha, beta, gamma);
    }
    
    /**
     * @brief      Project a geometry onto a basis (multi-patch)
     *
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param      coefs            The coefficients of the new geometry on \a projectionBasis
     *
     * @return     The L2 error of the projection
     */
    static T project(   const gsMultiBasis<T>         & integrationBasis,
                        const gsFunctionSet<T>        & projectionBasis,
                        const gsMultiPatch<T>         & geometryMap,
                              gsMatrix<T>             & coefs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        return _project(integrationBasis, projectionBasis, geometryMap, geometryMap, coefs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Project a geometry onto a basis (single patch)
     *
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param      coefs            The coefficients of the new geometry on \a projectionBasis
     *
     * @return     The L2 error of the projection
     */             
    static T project(   const gsBasis<T>              & projectionBasis,
                        const gsGeometry<T>           & geometryMap,
                              gsMatrix<T>             & coefs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiPatch<T> geometry(geometryMap);
        return _project(basis, basis, geometry, geometry, coefs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Project a geometry onto a basis (single patch)
     *
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param      coefs            The coefficients of the new geometry on \a projectionBasis
     *
     * @return     The L2 error of the projection
     */        
    static T project(   const gsBasis<T>              & projectionBasis,
                        const gsBasis<T>              & integrationBasis,
                        const gsGeometry<T>           & geometryMap,
                              gsMatrix<T>             & coefs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiBasis<T> intbasis(integrationBasis);
        gsMultiPatch<T> geometry(geometryMap);
        return _project(intbasis, basis, geometry, geometry, coefs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Project a function onto a basis
     *
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      coefs            The coefficients of the new geometry on \a projectionBasis
     *
     * @return     The L2 error of the projection
     */
    static T project(   const gsMultiBasis<T>         & projectionBasis,
                        const gsMultiPatch<T>         & geometryMap,
                        const gsFunctionSet<T>        & sourceFunction,
                              gsMatrix<T>             & coefs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        return _project(projectionBasis, projectionBasis, geometryMap, sourceFunction, coefs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Project a function onto a basis
     *
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      coefs            The coefficients of the new geometry on \a projectionBasis
     *
     * @return     The L2 error of the projection
     */
    static T project(   const gsFunctionSet<T>        & projectionBasis,
                        const gsMultiBasis<T>         & integrationBasis,
                        const gsMultiPatch<T>         & geometryMap,
                        const gsFunctionSet<T>        & sourceFunction,
                              gsMatrix<T>             & coefs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        return _project(integrationBasis, projectionBasis, geometryMap, sourceFunction, coefs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Project a function onto a basis
     *
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      coefs            The coefficients of the new geometry on \a projectionBasis
     *
     * @return     The L2 error of the projection
     */
    static T project(   const gsBasis<T>              & projectionBasis,
                        const gsGeometry<T>           & geometryMap,
                        const gsFunction<T>           & sourceFunction,
                              gsMatrix<T>             & coefs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiPatch<T> geometry(geometryMap);
        return _project(basis, basis, geometry, sourceFunction, coefs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Project a function onto a basis
     *
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      coefs            The coefficients of the new geometry on \a projectionBasis
     *
     * @return     The L2 error of the projection
     */
    static T project(   const gsBasis<T>              & projectionBasis,
                        const gsBasis<T>              & integrationBasis,
                        const gsGeometry<T>           & geometryMap,
                        const gsFunction<T>           & sourceFunction,
                              gsMatrix<T>             & coefs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiBasis<T> intbasis(integrationBasis);
        gsMultiPatch<T> geometry(geometryMap);
        return _project(intbasis, basis, geometry, sourceFunction, coefs, bc, options, alpha, beta, gamma);
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Implementations to obtain the system
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief      Obtain the system matrix and right-hand side for the L2 projection onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param      systemMatrix     The output system matrix
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     */
    static void system( const gsMultiBasis<T>         & projectionBasis,
                        const gsMultiPatch<T>         & geometryMap,
                              gsSparseMatrix<T>       & systemMatrix,
                              gsMatrix<T>             & rhs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        _system(projectionBasis, projectionBasis, geometryMap, geometryMap, systemMatrix, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the system matrix and right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      systemMatrix     The output system matrix
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void system( const gsMultiBasis<T>         & integrationBasis,
                        const gsFunctionSet<T>        & projectionBasis,
                        const gsMultiPatch<T>         & geometryMap,
                              gsSparseMatrix<T>       & systemMatrix,
                              gsMatrix<T>             & rhs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        _system(integrationBasis, projectionBasis, geometryMap, geometryMap, systemMatrix, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the system matrix and right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      systemMatrix     The output system matrix
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void system( const gsBasis<T>              & projectionBasis,
                        const gsGeometry<T>           & geometryMap,
                              gsSparseMatrix<T>       & systemMatrix,
                              gsMatrix<T>             & rhs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiPatch<T> geometry(geometryMap);
        _system(basis, basis, geometry, geometry, systemMatrix, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the system matrix and right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      systemMatrix     The output system matrix
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void system( const gsBasis<T>              & projectionBasis,
                        const gsBasis<T>              & integrationBasis,
                        const gsGeometry<T>           & geometryMap,
                              gsSparseMatrix<T>       & systemMatrix,
                              gsMatrix<T>             & rhs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiBasis<T> intbasis(integrationBasis);
        gsMultiPatch<T> geometry(geometryMap);
        _system(intbasis, basis, geometry, geometry, systemMatrix, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the system matrix and right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      systemMatrix     The output system matrix
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void system( const gsMultiBasis<T>         & projectionBasis,
                        const gsMultiPatch<T>         & geometryMap,
                        const gsFunctionSet<T>        & sourceFunction,
                              gsSparseMatrix<T>       & systemMatrix,
                              gsMatrix<T>             & rhs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        _system(projectionBasis, projectionBasis, geometryMap, sourceFunction, systemMatrix, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the system matrix and right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      systemMatrix     The output system matrix
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void system( const gsFunctionSet<T>        & projectionBasis,
                        const gsMultiBasis<T>         & integrationBasis,
                        const gsMultiPatch<T>         & geometryMap,
                        const gsFunctionSet<T>        & sourceFunction,
                              gsSparseMatrix<T>       & systemMatrix,
                              gsMatrix<T>             & rhs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        _system(integrationBasis, projectionBasis, geometryMap, sourceFunction, systemMatrix, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the system matrix and right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      systemMatrix     The output system matrix
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void system( const gsBasis<T>              & projectionBasis,
                        const gsGeometry<T>           & geometryMap,
                        const gsFunction<T>           & sourceFunction,
                              gsSparseMatrix<T>       & systemMatrix,
                              gsMatrix<T>             & rhs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiPatch<T> geometry(geometryMap);
        _system(basis, basis, geometry, sourceFunction, systemMatrix, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the system matrix and right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      systemMatrix     The output system matrix
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void system( const gsBasis<T>              & projectionBasis,
                        const gsBasis<T>              & integrationBasis,
                        const gsGeometry<T>           & geometryMap,
                        const gsFunction<T>           & sourceFunction,
                              gsSparseMatrix<T>       & systemMatrix,
                              gsMatrix<T>             & rhs,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiBasis<T> intbasis(integrationBasis);
        gsMultiPatch<T> geometry(geometryMap);
        _system(intbasis, basis, geometry, sourceFunction, systemMatrix, rhs, bc, options, alpha, beta, gamma);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Implementations to obtain the matrix
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief      Obtain the system matrix for the L2 projection of a geometry onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param      systemMatrix     The output system matrix
     * @param      options          The options that control the projection process
     *  
     * @return     The L2 error of the projection
     */
    static void matrix( const gsMultiBasis<T>         & projectionBasis,
                        const gsMultiPatch<T>         & geometryMap,
                              gsSparseMatrix<T>       & systemMatrix,
                              short_t                   targetDim = 1,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        _matrix(projectionBasis, projectionBasis, geometryMap, systemMatrix, targetDim, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the system matrix for the L2 projection of a geometry onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param      systemMatrix     The output system matrix
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void matrix( const gsMultiBasis<T>         & integrationBasis,
                        const gsFunctionSet<T>        & projectionBasis,
                        const gsMultiPatch<T>         & geometryMap,
                              gsSparseMatrix<T>       & systemMatrix,
                              short_t                   targetDim = 1,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        _matrix(integrationBasis, projectionBasis, geometryMap, systemMatrix, targetDim, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the system matrix for the L2 projection of a geometry onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param      systemMatrix     The output system matrix
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void matrix( const gsBasis<T>              & projectionBasis,
                        const gsGeometry<T>           & geometryMap,
                              gsSparseMatrix<T>       & systemMatrix,
                              short_t                   targetDim = 1,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiPatch<T> geometry(geometryMap);
        _matrix(basis, basis, geometry, systemMatrix, targetDim, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the system matrix for the L2 projection of a geometry onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param      systemMatrix     The output system matrix
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void matrix( const gsBasis<T>              & projectionBasis,
                        const gsBasis<T>              & integrationBasis,
                        const gsGeometry<T>           & geometryMap,
                              gsSparseMatrix<T>       & systemMatrix,
                              short_t                   targetDim = 1,
                        const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                        const gsOptionList            & options = gsOptionList(),
                        T alpha = 1.0,
                        T beta = 1.0,
                        T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiBasis<T> intbasis(integrationBasis);
        gsMultiPatch<T> geometry(geometryMap);
        _matrix(intbasis, basis, geometry, systemMatrix, targetDim, bc, options, alpha, beta, gamma);
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Implementations to obtain the rhs
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief      Obtain the right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void rhs(const gsMultiBasis<T>         & projectionBasis,
                    const gsMultiPatch<T>         & geometryMap,
                          gsMatrix<T>             & rhs,
                    const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                    const gsOptionList            & options = gsOptionList(),
                    T alpha = 1.0,
                    T beta = 1.0,
                    T gamma = 1.0)
    {
        _rhs(projectionBasis, projectionBasis, geometryMap, geometryMap, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void rhs(const gsMultiBasis<T>         & integrationBasis,
                    const gsFunctionSet<T>        & projectionBasis,
                    const gsMultiPatch<T>         & geometryMap,
                          gsMatrix<T>             & rhs,
                    const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                    const gsOptionList            & options = gsOptionList(),
                    T alpha = 1.0,
                    T beta = 1.0,
                    T gamma = 1.0)
    {
        _rhs(integrationBasis, projectionBasis, geometryMap, geometryMap, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void rhs(const gsBasis<T>              & projectionBasis,
                    const gsGeometry<T>           & geometryMap,
                          gsMatrix<T>             & rhs,
                    const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                    const gsOptionList            & options = gsOptionList(),
                    T alpha = 1.0,
                    T beta = 1.0,
                    T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiPatch<T> geometry(geometryMap);
        _rhs(basis, basis, geometry, geometry, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void rhs(const gsBasis<T>              & projectionBasis,
                    const gsBasis<T>              & integrationBasis,
                    const gsGeometry<T>           & geometryMap,
                          gsMatrix<T>             & rhs,
                    const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                    const gsOptionList            & options = gsOptionList(),
                    T alpha = 1.0,
                    T beta = 1.0,
                    T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiBasis<T> intbasis(integrationBasis);
        gsMultiPatch<T> geometry(geometryMap);
        _rhs(intbasis, basis, geometry, geometry, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void rhs(const gsMultiBasis<T>         & projectionBasis,
                    const gsMultiPatch<T>         & geometryMap,
                    const gsFunctionSet<T>        & sourceFunction,
                          gsMatrix<T>             & rhs,
                    const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                    const gsOptionList            & options = gsOptionList(),
                    T alpha = 1.0,
                    T beta = 1.0,
                    T gamma = 1.0)
    {
        _rhs(projectionBasis, projectionBasis, geometryMap, sourceFunction, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void rhs(const gsFunctionSet<T>        & projectionBasis,
                    const gsMultiBasis<T>         & integrationBasis,
                    const gsMultiPatch<T>         & geometryMap,
                    const gsFunctionSet<T>        & sourceFunction,
                          gsMatrix<T>             & rhs,
                    const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                    const gsOptionList            & options = gsOptionList(),
                    T alpha = 1.0,
                    T beta = 1.0,
                    T gamma = 1.0)
    {
        _rhs(integrationBasis, projectionBasis, geometryMap, sourceFunction, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void rhs(const gsBasis<T>              & projectionBasis,
                    const gsGeometry<T>           & geometryMap,
                    const gsFunction<T>           & sourceFunction,
                          gsMatrix<T>             & rhs,
                    const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                    const gsOptionList            & options = gsOptionList(),
                    T alpha = 1.0,
                    T beta = 1.0,
                    T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiPatch<T> geometry(geometryMap);
        _rhs(basis, basis, geometry, sourceFunction, rhs, bc, options, alpha, beta, gamma);
    }

    /**
     * @brief      Obtain the right-hand side for the L2 projection of a function onto a basis
     * 
     * @param[in]  projectionBasis  The basis to project on
     * @param[in]  integrationBasis The basis used for numerical integration.
     * @param[in]  geometryMap      The geometry
     * @param[in]  sourceFunction   The source function
     * @param      rhs              The output right-hand side vector
     * @param      options          The options that control the projection process
     * 
     * @return     The L2 error of the projection
     */
    static void rhs(const gsBasis<T>              & projectionBasis,
                    const gsBasis<T>              & integrationBasis,
                    const gsGeometry<T>           & geometryMap,
                    const gsFunction<T>           & sourceFunction,
                          gsMatrix<T>             & rhs,
                    const gsBoundaryConditions<T> & bc = gsBoundaryConditions<T>(),
                    const gsOptionList            & options = gsOptionList(),
                    T alpha = 1.0,
                    T beta = 1.0,
                    T gamma = 1.0)
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiBasis<T> intbasis(integrationBasis);
        gsMultiPatch<T> geometry(geometryMap);
        _rhs(intbasis, basis, geometry, sourceFunction, rhs, bc, options, alpha, beta, gamma);
    }

}; //struct

// Type aliases for convenience
template <class T>
using gsL2Projection = gsProjection<ProjectionNorm::L2, T>;

template <class T>
using gsH1Projection = gsProjection<ProjectionNorm::H1, T>;

template <class T>
using gsH2Projection = gsProjection<ProjectionNorm::H2, T>;

#ifdef GISMO_WITH_PYBIND11

    /**
     * @brief Initializes the Python wrapper for the ProjectionNorm enum
     */
    void pybind11_enum_gsProjectionNorm(pybind11::module &m);

    /**
     * @brief Initializes the Python wrapper for the class: gsProjection
     */
    template<ProjectionNorm Norm>
    void pybind11_init_gsProjection(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsProjection.hpp)
#endif
