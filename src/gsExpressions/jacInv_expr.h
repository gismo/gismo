/** @file jacInv_expr.h

    @brief Defines the inverse jacobian expression

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
 * @brief Expression for the inverse of the Jacobian matrix
 * @ingroup Expressions
 * @tparam T The expression type
 */
template<class T>
class jacInv_expr  : public _expr<jacInv_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;
public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

    jacInv_expr(const gsGeometryMap<T> & G) : _G(G)
    {
        // Note: for non-square Jacobian matrices, generalized inverse, i.e.: (J^t J)^{-t} J^t
        //GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
    }

    MatExprType eval(const index_t k) const { return _G.data().jacInvTr.reshapeCol(k,cols(),rows()).transpose(); }

    index_t rows() const { return _G.source().domainDim(); }
    index_t cols() const { return _G.source().targetDim(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_GRAD_TRANSFORM;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<T>::get();}

    // todo mat_expr ?
    // tr() const --> _G.data().fundForm(k)

    void print(std::ostream &os) const { os << "jacInv("; _G.print(os); os <<")"; }
};

}// namespace expr
}// namespace gismo