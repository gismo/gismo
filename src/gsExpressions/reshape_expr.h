/** @file reshape_expr.h

    @brief Defines the reshape expression

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
 * @brief Expression for reshaping an expression
 * @ingroup Expressions
 * @tparam E The type of the expression
 */
template<class E>
class reshape_expr  : public _expr<reshape_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, ColBlocks = E::ColBlocks};
    enum {Space = E::Space};
private:
    typename E::Nested_t _u;
    index_t _n, _m;
    mutable gsMatrix<Scalar> tmp;

public:

    //the reshaping is done column-wise
    reshape_expr(_expr<E> const& u, index_t n, index_t m) : _u(u), _n(n), _m(m)
    {
        //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
    }

    const gsAsConstMatrix<Scalar> eval(const index_t k) const
    {
        // Note: this assertion would fail in the constructor!
        GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension "
                      << _u.rows() << " x "<<_u.cols() << "!=" << _n << " * "<< _m );
        tmp = _u.eval(k);
        return gsAsConstMatrix<Scalar>(tmp.data(),_n,_m);
    }

    index_t rows() const { return _n; }
    index_t cols() const { return _m; }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    void print(std::ostream &os) const { os << "reshape("; _u.print(os); os<<","<<_n<<","<<_m<<")"; }
};


/**
 * @brief Reshape an expression
 * @ingroup Expressions
 * @tparam E The type of the expression
 * @param u The expression to reshape
 * @param n The number of rows in the reshaped expression
 * @param m The number of columns in the reshaped expression
 */
template <class E> EIGEN_STRONG_INLINE
reshape_expr<E> const reshape(E const & u, index_t n, index_t m)
{ return reshape_expr<E>(u, n, m); }

}// namespace expr
}// namespace gismo