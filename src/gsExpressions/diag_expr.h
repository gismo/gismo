/** @file diag_expr.h

    @brief Defines the diag expression

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
  * @brief Expression for the diagonal of a matrix expression
  * @ingroup Expressions
  * @tparam E The expression type
  */
template<class E>
class diag_expr  : public _expr<diag_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {Space=0, ColBlocks=E::ColBlocks, ScalarValued = 0};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;

public:
    diag_expr(_expr<E> const& u) : _u(u)
    {
        GISMO_ASSERT(0== _u.cols()%_u.rows(), "Expecting square-block expression, got "
        << _u.rows() <<" x "<< _u.cols() );
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        // Assume mat ??
        MatExprType tmp = _u.eval(k);
        const index_t cb = _u.rows();
        const index_t r  = _u.cols() / cb;
        res.resize(r, cb);
        for (index_t i = 0; i!=r; ++i)
            res.row(i) = tmp.middleCols(i*cb,cb).diagonal();
        return res;
    }

    index_t rows() const { return _u.cols() / _u.rows(); }
    index_t cols() const { return _u.rows(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }


    void print(std::ostream &os) const { os << "diag("; _u.print(os); os<<")"; }
};

/**
 * @brief Returns the diagonal of a matrix expression
 * @param u The expression
 * @ingroup Expressions
 * @return A diagonal matrix expression
 */
template <typename E> EIGEN_STRONG_INLINE
diag_expr<E> const diagonal(E const & u)
{ return diag_expr<E>(u); }

}// namespace expr
}// namespace gismo