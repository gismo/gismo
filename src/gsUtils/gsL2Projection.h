/** @file gsL2Projection.h

    @brief Class that performs an L2 projection

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/


#pragma once

#include <gsAssembler/gsExprAssembler.h>

namespace gismo {

/** \brief Class that performs an L2 projection

    \tparam T coefficient type
 */

template <class T>
struct gsL2Projection
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
    static T _project(  const gsMultiBasis<T>  & integrationBasis,
                        const gsFunctionSet<T> & projectionBasis,
                        const gsFunctionSet<T> & geometryMap,
                        const gsFunctionSet<T> & sourceFunction,
                              gsMatrix<T>      & coefs,
                        const gsOptionList     & options);


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
    static T project(   const gsMultiBasis<T> & projectionBasis,
                        const gsMultiPatch<T> & geometryMap,
                              gsMatrix<T>     & coefs,
                        const gsOptionList    & options = gsOptionList())
    {
        return _project(projectionBasis, projectionBasis, geometryMap, geometryMap, coefs, options);   
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
    static T project(   const gsFunctionSet<T>& integrationBasis,
                        const gsMultiBasis<T> & projectionBasis,
                        const gsMultiPatch<T> & geometryMap,
                              gsMatrix<T>     & coefs,
                        const gsOptionList    & options = gsOptionList())
    {
        return _project(integrationBasis, projectionBasis, geometryMap, geometryMap, coefs, options);   
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
    static T project(   const gsBasis<T>      & projectionBasis,
                        const gsGeometry<T>   & geometryMap,
                              gsMatrix<T>     & coefs,
                        const gsOptionList    & options = gsOptionList())
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiPatch<T> geometry(geometryMap);
        return _project(basis, basis, geometry, geometry, coefs, options);   
        
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
    static T project(   const gsBasis<T>      & projectionBasis,
                        const gsBasis<T>      & integrationBasis,
                        const gsGeometry<T>   & geometryMap,
                              gsMatrix<T>     & coefs,
                        const gsOptionList    & options = gsOptionList())
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiBasis<T> intbasis(integrationBasis);
        gsMultiPatch<T> geometry(geometryMap);
        return _project(intbasis, basis, geometry, geometry, coefs, options);
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
    static T project(   const gsMultiBasis<T>  & projectionBasis,
                        const gsMultiPatch<T>  & geometryMap,
                        const gsFunctionSet<T> & sourceFunction,
                              gsMatrix<T>      & coefs,
                        const gsOptionList     & options = gsOptionList())
    {
        return _project(projectionBasis, projectionBasis, geometryMap, sourceFunction, coefs, options);
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
    static T project(   const gsFunctionSet<T> & projectionBasis,
                        const gsMultiBasis<T>  & integrationBasis,
                        const gsMultiPatch<T>  & geometryMap,
                        const gsFunctionSet<T> & sourceFunction,
                              gsMatrix<T>      & coefs,
                        const gsOptionList     & options = gsOptionList())
    {
        return _project(projectionBasis, projectionBasis, geometryMap, sourceFunction, coefs, options);
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
    static T project(   const gsBasis<T>       & projectionBasis,
                        const gsGeometry<T>    & geometryMap,
                        const gsFunction<T>    & sourceFunction,
                              gsMatrix<T>      & coefs,
                        const gsOptionList     & options = gsOptionList())
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiPatch<T> geometry(geometryMap);
        return _project(basis, basis, geometry, sourceFunction, coefs, options);
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
    static T project(   const gsBasis<T>       & projectionBasis,
                        const gsBasis<T>       & integrationBasis,
                        const gsGeometry<T>    & geometryMap,
                        const gsFunction<T>    & sourceFunction,
                              gsMatrix<T>      & coefs,
                        const gsOptionList     & options = gsOptionList())
    {
        gsMultiBasis<T> basis(projectionBasis);
        gsMultiBasis<T> intbasis(integrationBasis);
        gsMultiPatch<T> geometry(geometryMap);
        return _project(intbasis, basis, geometry, sourceFunction, coefs, options);
    }

    /**
     * @brief      Projects a \a source geometry onto \a basis and returns it in
     *             \a result
     *
     * @param[in]  basis     The basis to project on
     * @param[in]  geometry  The geometry
     * @param      result    The coefficients of the new geometry on \a basis
     *
     * @return     The L2 error of the projection
     */
    GISMO_DEPRECATED
    static T projectGeometry(   const gsBasis<T> & basis,
                                const gsGeometry<T> & geometry,
                                gsMatrix<T> & result);

    /**
     * @brief      Projects a \a source geometry onto \a basis and returns it in
     *             \a result
     *
     * @param[in]  basis     The basis to project on
     * @param[in]  geometry  The geometry
     * @param      result    The coefficients of the new geometry on \a basis
     *
     * @return     The L2 error of the projection
     */
    GISMO_DEPRECATED
    static T projectGeometry(   const gsMultiBasis<T> & basis,
                                const gsFunctionSet<T> & geometry,
                                gsMatrix<T> & result);

    /**
     * @brief      Projects a \a source geometry onto \a basis and returns it in
     *             \a result
     *
     * @param[in]  basis     The basis to project on
     * @param[in]  geometry  The geometry
     * @param      result    The new geometry
     *
     * @return     The L2 error of the projection
     */
    GISMO_DEPRECATED
    static T projectGeometry(   const gsMultiBasis<T> & basis,
                                const gsFunctionSet<T> & geometry,
                                gsMultiPatch<T> & result);

    /**
     * @brief      Projects a \a source geometry onto \a basis and returns it in
     *             \a result
     *
     * @param[in]  intbasis  The basis used for quadrature
     * @param[in]  basis     The mapped basis to project on
     * @param[in]  geometry  The geometry
     * @param      result    The coefficients of the new geometry on \a basis
     *
     * @return     The L2 error of the projection
     */
    GISMO_DEPRECATED
    static T projectGeometry(   const gsMultiBasis<T> & intbasis,
                                const gsMappedBasis<2,T> & basis,
                                const gsFunctionSet<T> & geometry,
                                gsMatrix<T> & result);

    /**
     * @brief      Projects a function on a basis
     *
     * @param[in]  basis     The basis to project on
     * @param[in]  source    The source function
     * @param[in]  geometry  The geometry to evaluate the function on
     * @param      result    The coefficients of the function
     *
     * @return     The L2 error of the projection
     */
    GISMO_DEPRECATED
    static T projectFunction(   const gsMultiBasis<T> & basis,
                                const gsFunctionSet<T> & source,
                                const gsMultiPatch<T> & geometry,
                                gsMatrix<T> & result);

    /**
     * @brief      Projects a function on a basis
     *
     * @param[in]  basis     The basis to project on
     * @param[in]  source    The source function
     * @param[in]  geometry  The geometry to evaluate the function on
     * @param      result    The function as a multipatch
     *
     * @return     The L2 error of the projection
     */
    GISMO_DEPRECATED
    static T projectFunction(   const gsMultiBasis<T> & basis,
                                const gsFunctionSet<T> & source,
                                const gsMultiPatch<T> & geometry,
                                gsMultiPatch<T> & result);

    /**
     * @brief      Projects a function on a basis
     *
     * @param[in]  intbasis  The basis used for quadrature
     * @param[in]  basis     The basis to project on
     * @param[in]  source    The source function
     * @param[in]  geometry  The geometry to evaluate the function on
     * @param      result    The function as a multipatch
     *
     * @return     The L2 error of the projection
     */
    GISMO_DEPRECATED
    static T projectFunction(   const gsMultiBasis<T>   & intbasis,
                                const gsMappedBasis<2,T>& basis,
                                const gsFunctionSet<T>  & source,
                                const gsMultiPatch<T>   & geometry,
                                gsMatrix<T> & result);

    /**
     * @brief      Projects a \a source geometry onto \a basis and returns it in
     *             \a result. Fixes the boundaries
     *
     * @param[in]  basis     The basis to project on
     * @param[in]  geometry  The geometry
     * @param      result    The coefficients of the new geometry on \a basis
     *
     * @return     The L2 error of the projection
     */
    GISMO_DEPRECATED
    static T projectGeometryBoundaries( const gsMultiBasis<T> & basis,
                                        const gsMultiPatch<T> & geometry,
                                        gsMultiPatch<T> & result);


    /**
     * @brief      Projects a \a source geometry onto \a basis and returns it in
     *             \a result. Penalizes interfaces and boundaries
     *
     * @param[in]  basis     The basis to project on
     * @param[in]  geometry  The geometry
     * @param      result    The coefficients of the new geometry on @a basis
     * @param[in]  penalty   The penalty factor
     *
     * @return     The L2 error of the projection
     */
    GISMO_DEPRECATED
    static T projectGeometryPenalty(    const gsMultiBasis<T> & basis,
                                        const gsMultiPatch<T> & geometry,
                                        gsMultiPatch<T> & result,
                                        T penalty = 1e3);


}; //struct

} // gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsL2Projection.hpp)
#endif
