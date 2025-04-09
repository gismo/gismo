/** @file rowsum_expr.h

    @brief Defines the rowsum expression

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
 * @brief Expression for the row sum of a matrix
 * @ingroup Expressions
 * @tparam E The type of the expression
 */
template<class E>
class rowsum_expr  : public _expr<rowsum_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, Space = E::Space, ColBlocks = E::ColBlocks};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> tmp;

public:

    rowsum_expr(_expr<E> const& u) : _u(u)
    {
        //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        tmp = _u.eval(k).rowwise().sum();
        return tmp;
    }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return 1; }
    void setFlag() const { _u.setFlag(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "rowsum("; _u.print(os); os<<")"; }
};

}// namespace expr
}// namespace gismo