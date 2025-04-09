/** @file asdiag_expr.h

    @brief Defines the asdiag expression

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
 * @brief   Expression for turning a vector into a diagonal matrix
 * @ingroup Expressions
 * @tparam E The expression type
 */
template<class E>
class asdiag_expr : public _expr<asdiag_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;

public:
    enum{Space = E::Space, ScalarValued= 0, ColBlocks= E::ColBlocks};

    asdiag_expr(_expr<E> const& u) : _u(u) { }

public:

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto m = _u.eval(k);
        const index_t r = m.rows();
        const index_t c = m.cols();
        res.resize(r,r*c);
        for (index_t i = 0; i!=c; ++i)
            res.middleCols(i*r,r) = m.col(i).asDiagonal();
        return res;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.rows() * _u.cols(); }
    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    void print(std::ostream &os) const { os << "diag("; _u.print(os); os <<")";}
};

}// namespace expr
}// namespace gismo