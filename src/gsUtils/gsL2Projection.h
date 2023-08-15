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

public:
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
    static T projectGeometryPenalty(    const gsMultiBasis<T> & basis,
                                        const gsMultiPatch<T> & geometry,
                                        gsMultiPatch<T> & result,
                                        T penalty = 1e3);


}; //struct

} // gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsL2Projection.hpp)
#endif
