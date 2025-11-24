/** @file _expr_macros.h

    @brief This file provides a macros for expressions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
               H.M. Verhelst
*/

#pragma once

namespace gismo
{
namespace expr
{

// Shortcuts for common quantities, for instance function
// transformations by the geometry map \a G
#define GISMO_SHORTCUT_VAR_EXPRESSION(name,impl,docstring)  \
/**                                                         \
     * @brief   docstring                                   \
     * @param u The variable                                \
     * @ingroup Expressions                                 \
*/                                                          \
    template<class E> EIGEN_STRONG_INLINE                   \
    auto name(const E & u) -> decltype(impl) { return impl; }
#define GISMO_SHORTCUT_MAP_EXPRESSION(name,impl,docstring)  \
/**                                                         \
     * @brief   docstring                                   \
     * @param G The geometry map                            \
     * @ingroup Expressions                                 \
*/                                                          \
    template<class T> EIGEN_STRONG_INLINE                   \
    auto name(const gsGeometryMap<T> & G)  -> decltype(impl) { return impl; }
#define GISMO_SHORTCUT_PHY_EXPRESSION(name,impl,docstring)  \
/**                                                         \
     * @brief   docstring                                   \
     * @param u The variable                                \
     * @param G The geometry map                            \
     * @ingroup Expressions                                 \
*/                                                          \
    template<class E> EIGEN_STRONG_INLINE                   \
    auto name(const E & u, const gsGeometryMap<typename E::Scalar> & G)  -> decltype(impl) { return impl; }

GISMO_SHORTCUT_VAR_EXPRESSION(  div, gismo::expr::jac(u).trace(),
                                Divergence of \a u with respect to the parametric domain )



GISMO_SHORTCUT_MAP_EXPRESSION(  unv, gismo::expr::nv(G).normalized(),
                                The normalized normal vector on the boundary)

GISMO_SHORTCUT_MAP_EXPRESSION(  usn, gismo::expr::sn(G).normalized(),
                                The normalized surface normal vector)

GISMO_SHORTCUT_VAR_EXPRESSION(  igrad, gismo::expr::grad(u),
                                Gradient of \a u with respect to the parametric domain) // u is presumed to be defined over G

GISMO_SHORTCUT_PHY_EXPRESSION(  igrad, gismo::expr::grad(u)*gismo::expr::jac(G).ginv(),
                                Gradient of \a u with respect to the physical domain ) // transpose() problem ??

GISMO_SHORTCUT_PHY_EXPRESSION(  ijac, gismo::expr::jac(u) * gismo::expr::jac(G).ginv(),
                                Jacobian of \a u with respect to the physical domain)

GISMO_SHORTCUT_VAR_EXPRESSION(  ihess, gismo::expr::hess(u),
                                Hessian of \a u with respect to the parametric domain )

GISMO_SHORTCUT_PHY_EXPRESSION(  idiv, gismo::expr::ijac(u,G).trace(),
                                Divergence of \a u with respect to the physical domain    )

// @note: does this work for non-scalar solutions? (todo)
GISMO_SHORTCUT_PHY_EXPRESSION(ihess,
                              gismo::expr::jac(G).ginv().tr()*( gismo::expr::hess(u) - summ(igrad(u,G),gismo::expr::hess(G)) ) * gismo::expr::jac(G).ginv(),
                              Hessian of \a u with respect to the physical domain)

                              GISMO_SHORTCUT_VAR_EXPRESSION(  ilapl, hess(u).trace(),
                                Hessian of \a u with respect to the parametric domain )

GISMO_SHORTCUT_PHY_EXPRESSION(  ilapl, gismo::expr::ihess(u,G).trace(),
                                Hessian of \a u with respect to the physical domain   )

GISMO_SHORTCUT_VAR_EXPRESSION(  fform, gismo::expr::jac(u).tr()*gismo::expr::jac(u),
                                Fundamental form )

GISMO_SHORTCUT_VAR_EXPRESSION(shapeop, gismo::expr::fform(u).inv() * gismo::expr::fform2nd(u),
                              ?? @note todo )

#undef GISMO_SHORTCUT_PHY_EXPRESSION
#undef GISMO_SHORTCUT_VAR_EXPRESSION
#undef GISMO_SHORTCUT_MAP_EXPRESSION

}// namespace expr
}// namespace gismo