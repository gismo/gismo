/** @file dJacdc_expr.h

    @brief Defines the dJacdc expression

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
 * @brief Expression for the derivative of the jacobian of a geometry map
 *        with respect to a coordinate c. It returns a matrix with the gradient of u in row d.
 * @ingroup Expressions
 * @tparam E The expression type
 */
template<class E>
class dJacdc_expr : public _expr<dJacdc_expr<E> >
{
    typename E::Nested_t _u;
public:
    enum{ Space = E::Space, ScalarValued = 0, ColBlocks = (1==E::Space?1:0)};

    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;
    index_t _c;

    dJacdc_expr(const E & u, index_t c) : _u(u), _c(c)
    { GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");}

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        index_t dd = _u.source().domainDim();
        index_t n = _u.rows();
        res.setZero(dd, dd*n);

        gsMatrix<Scalar> grad = _u.data().values[1].reshapeCol(k, dd, n);
        for(index_t i = 0; i < n; i++){
            res.row(_c).segment(i*dd,dd) = grad.col(i);
        }
        return res;
    }

    index_t rows() const { return _u.source().domainDim(); }

    index_t cols() const { return _u.source().domainDim()*_u.rows(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_GRAD;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "dJacdc("; _u.print(os); os <<")"; }
};

/**
 * @brief Returns the derivative of the jacobian of a geometry map with respect to a coordinate
 * @param u The expression
 * @param c The coordinate
 */
template<class E> EIGEN_STRONG_INLINE
dJacdc_expr<E> dJacdc(const E & u, index_t c) { return dJacdc_expr<E>(u,c); }

}// namespace expr
}// namespace gismo