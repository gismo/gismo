/** @file replicate_expr.h

    @brief Defines the replicate expression

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
 * @brief Expression for the replicate the evaluation of an expression \f$n\times m\f$ times
 * @ingroup Expressions
 * @tparam E The type of the expression
 */
template<class E>
class replicate_expr  : public _expr<replicate_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, Space = E::Space, ColBlocks= E::ColBlocks};
private:
    typename E::Nested_t _u;
    index_t _n, _m;
    mutable gsMatrix<Scalar> tmp;

public:

    //the replicate is done nxm times
    replicate_expr(_expr<E> const& u, index_t n, index_t m) : _u(u), _n(n), _m(m)
    {
    }

    auto eval(const index_t k) const -> decltype(tmp.replicate(_n,_m))
    {
        tmp = _u.eval(k);
        return tmp.replicate(_n,_m);
    }

    index_t rows() const { return _n*_u.rows(); }
    index_t cols() const { return _m*_u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "replicate("; _u.print(os); os<<","<<_n<<","<<_m<<")"; }
};

/**
 * @brief Creates a replicate expression by replicating the given expression
 *        \f[ \text{replicate}(u) = u \otimes I_{n \times m} \f]
 * @ingroup Expressions
 * @tparam E The type of the expression being replicated
 * @param u The input expression to be replicated.
 */
template <typename E> EIGEN_STRONG_INLINE
replicate_expr<E> const replicate(E const & u, index_t n, index_t m = 1)
{ return replicate_expr<E>(u, n, m); }

}// namespace expr
}// namespace gismo