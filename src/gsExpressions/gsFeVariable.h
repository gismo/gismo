/** @file gsFeVariable.h

    @brief Defines an expression for a variable

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

/**
 * @brief Expression for a finite element variables or PDE coefficient functionals.
 *        This can be e.g. a diffusion coefficient, or an isogeometric function.
 * @ingroup Expressions
 * @tparam T The expression type
 */
template<class T>
class gsFeVariable  : public symbol_expr< gsFeVariable<T> >
{
    friend class gismo::gsExprHelper<T>;
    typedef symbol_expr< gsFeVariable<T> > Base;
protected:
    explicit gsFeVariable(index_t _d = 1) : Base(_d) { }
public:
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};
};


template <typename T>
struct expr_traits<gsFeVariable<T> >
{
    typedef T Scalar;
    typedef const gsFeVariable<T> Nested_t;
};


}// namespace expr
}// namespace gismo