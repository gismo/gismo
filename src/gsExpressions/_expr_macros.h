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
    auto name(const E & u) { return impl; }
#define GISMO_SHORTCUT_MAP_EXPRESSION(name,impl,docstring)  \
/**                                                         \
     * @brief   docstring                                   \
     * @param G The geometry map                            \
     * @ingroup Expressions                                 \
*/                                                          \
    template<class T> EIGEN_STRONG_INLINE                   \
    auto name(const gsGeometryMap<T> & G) { return impl; }
#define GISMO_SHORTCUT_PHY_EXPRESSION(name,impl,docstring)  \
/**                                                         \
     * @brief   docstring                                   \
     * @param u The variable                                \
     * @param G The geometry map                            \
     * @ingroup Expressions                                 \
*/                                                          \
    template<class E> EIGEN_STRONG_INLINE                   \
    auto name(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return impl; }

GISMO_SHORTCUT_VAR_EXPRESSION(  div, jac(u).trace(),
                                Divergence of \a u with respect to the parametric domain )

GISMO_SHORTCUT_PHY_EXPRESSION(  idiv, ijac(u,G).trace(),
                                Divergence of \a u with respect to the physical domain    )

// Normal vector (normalized). Use a fallback implementation when compiling
// with NVCC (Kokkos/NVCC wrapper) because some NVCC host/device parsing
// paths have trouble with the direct `.normalized()` call on expression
// templates. The fallback computes normalization via division by norm().
template<class T> EIGEN_STRONG_INLINE
auto unv(const gsGeometryMap<T> & G)
{
#ifdef __CUDACC__
    auto tmp = onormal_expr<T>(G);
    return tmp / tmp.norm();
#else
    return onormal_expr<T>(G).normalized();
#endif
}

template<class T> EIGEN_STRONG_INLINE
auto usn(const gsGeometryMap<T> & G)
{
#ifdef __CUDACC__
    auto tmp = normal_expr<T>(G);
    return tmp / tmp.norm();
#else
    return normal_expr<T>(G).normalized();
#endif
}

GISMO_SHORTCUT_VAR_EXPRESSION(  igrad, grad(u),
                                Gradient of \a u with respect to the parametric domain) // u is presumed to be defined over G

GISMO_SHORTCUT_PHY_EXPRESSION(  igrad, grad(u)*jac(G).ginv(),
                                Gradient of \a u with respect to the physical domain ) // transpose() problem ??

GISMO_SHORTCUT_PHY_EXPRESSION(  ijac, jac(u) * jac(G).ginv(),
                                Jacobian of \a u with respect to the physical domain)

GISMO_SHORTCUT_VAR_EXPRESSION(  ihess, hess(u),
                                Hessian of \a u with respect to the parametric domain )

// @note: does this work for non-scalar solutions? (todo)
GISMO_SHORTCUT_PHY_EXPRESSION(ihess,
                              jac(G).ginv().tr()*( hess(u) - summ(igrad(u,G),hess(G)) ) * jac(G).ginv(),
                              Hessian of \a u with respect to the physical domain)

                              GISMO_SHORTCUT_VAR_EXPRESSION(  ilapl, hess(u).trace(),
                                Hessian of \a u with respect to the parametric domain )

GISMO_SHORTCUT_PHY_EXPRESSION(  ilapl, ihess(u,G).trace(),
                                Hessian of \a u with respect to the physical domain   )

GISMO_SHORTCUT_VAR_EXPRESSION(  fform, jac(u).tr()*jac(u),
                                Fundamental form )

GISMO_SHORTCUT_VAR_EXPRESSION(shapeop, fform(u).inv() * fform2nd(u),
                              ?? @note todo )

#undef GISMO_SHORTCUT_PHY_EXPRESSION
#undef GISMO_SHORTCUT_VAR_EXPRESSION
#undef GISMO_SHORTCUT_MAP_EXPRESSION

}// namespace expr
}// namespace gismo