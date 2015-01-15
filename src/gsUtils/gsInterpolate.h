/** @file gsInterpolate.h

    @brief Provides functions to interpolate data

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris

*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>

namespace gismo {

/** 
    Compute an interpolating function in basis \a g which has the given values \a vals
    at the given interpolation points \a pts in the parameter space.

    \tparam T   coefficient type
    \param  g   The basis in which to interpolate.
    \param pts  The interpolation nodes in parameter space.
    \param vals The values to be interpolated.
    \return     A new geometry object defined over \a g which interpolates the given values.
    
    \ingroup Utils
*/
template<class T>
gsGeometry<T> * gsInterpolate( gsBasis<T> const& g, gsMatrix<T> 
                               const& pts, gsMatrix<T> const& vals );


/** 
    Perform Lagrangian interpolation for the function \a f using the basis' anchors as interpolation points.

    \tparam T   coefficient type
    \param  g   The basis in which to interpolate the function.
    \param  f   The function to be interpolated.
    \return     A new geometry object defined over \a g which interpolates \a f
    
    \ingroup Utils
*/
template<class T>
gsGeometry<T> * gsInterpolate( const gsBasis<T>& g, const gsFunction<T>& f );


/** \brief Interpolate a function \em f at the boundary of a domain in the basis \em basis.
 *
 * The driving motivation for this function is to provide data for the strong
 * incorporation of Dirichlet boundary conditions when solving
 * partial differential equations.
 *
 * For the given basis <b>\em basis</b>, the boundaries/sides of the parameter domain
 * are mapped by <b>\em geo</b> to a physical domain. At this physical domain, an
 * L2-projection of the function <b>\em f</b> into the (push-forward of) \em basis
 * is performed.
 *
 * With input argument <b>\em Sides</b>, you specify which/how many sides of the
 * parameter domain should be used.\n
 * E.g., in 2D, \n
 * if <em>Sides = (1,2,3,4)</em>, the L2-projection will be done
 * on all boundaries, \n
 * if <em>Sides = (1,3)</em>, the L2-projection will be
 * done at the western and southern sides (see also gismo::boundary).
 *
 * In/out-parameter <b>\em vecIdx</b> and <b>\em vecCoeff</b> will be
 * overwritten with the
 * indices of the degrees of freedom (DOF) at the boundary sides and with the
 * corresponding coefficient values, respectively.
 *
 * The flag <b>\em getGlobalData</b> specifies, whether the returned
 * data should be global in the following sense:\n
 * If <em>getGlobalData = false (default)</em>, only those indices and values
 * are returned which correspond to the specified boundary sides.
 * This data can be used to strongly enforce Dirichlet boundary conditions
 * using the function gsPoissonAssembler::boundaryFixDofs.\n
 * If <em>getGlobalData = true</em>, the size of \em vecCoeff equals the
 * number of basis functions in \em basis. All coefficients which do not
 * correspond to the specified boundary sides are set to zero.
 * This can be used to create a gsGeometry that interpolates the
 * function \em f at the specified boundary sides and is zero everywhere else.
 *
 *
 * \param[in] basis The basis in which the function should be interpolated
 * \param[in] f Function to be interpolated (defined on the physical domain).
 * \param[in] geo Geometry of the domain.
 * \param[in] Sides Specifies which sides of the domain should be used.\n
 * Entries should be from 1,2,3,4 in 2D \n and from 1,...,6 in 3D
 * (see also gismo::boundary).
 * \param[in,out] vecIdx Gets \b overwritten by the (global) indices of the
 * DOF at the given sides.
 * \param[in,out] vecCoeff Gets \b overwritten by the coefficients corresponding
 * to the DOFs in \em vecIdx such that they represent the projection of \em f.
 * \param[in] getGlobalData When set to \em true, \em vecIdx and \em vecCoeff
 * will include values for \em all DOF, where the DOFS which are not used for
 * the interpolation are set to zero (useful when you want to create a
 * gsGeometry with \em basis and \em vecCoeff).
 *
 *
 * \ingroup Utils
 */
template <class T>
void gsL2ProjectOnBoundary(const gsBasis<T> & basis,
                           const gsFunction<T> & f,
                           const gsGeometry<T> & geo,
                           const gsVector<int> & Sides,
                           gsVector<unsigned> & vecIdx,
                           gsMatrix<T> & vecCoeff,
                           bool getGlobalData = false);


}// namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsInterpolate.hpp)
#endif
