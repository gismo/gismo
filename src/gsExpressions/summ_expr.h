/** @file summ_expr.h

    @brief Defines the summ expression

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
 * @brief Expression for the summation operation
 *        M [r x r*k] is a list of matrices
 *        Summation is done over k,
 *        [M1 M2 .. Mk]
 *        u [s x k] is a list of vectors
 *        Computed quantity is of size [r x r*s] and contains
 *        [ ... sum_k(Mk * u(s,k) ) ... ]_s
 * @ingroup Expressions
 * @tparam E1 the type of the first expression
 */
template <typename E1, typename E2>
class summ_expr : public _expr<summ_expr<E1,E2> >
{
public:
    typedef typename E1::Scalar Scalar;

    enum {Space = E1::Space, ScalarValued= 0, ColBlocks= E2::ColBlocks};

    summ_expr(E1 const& u, E2 const& M) : _u(u), _M(M) { }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto sl   = _u.eval(k);
        const index_t sr = sl.rows();
        auto ml   = _M.eval(k);
        const index_t mr = ml.rows();
        const index_t mb = _M.cardinality();

        GISMO_ASSERT(_M.cols()==_M.rows(),"Matrix must be square: "<< _M.rows()<<" x "<< _M.cols() << " expr: "<< _M );
        GISMO_ASSERT(mb==_u.cols(),"cardinality must match vector, but card(M)="<<_M.cardinality()<<" and cols(u)="<<_u.cols());

        res.setZero(mr, sr * mr);
        for (index_t i = 0; i!=sr; ++i)
            for (index_t t = 0; t!=mb; ++t) // lc
                res.middleCols(i*mr,mr) += sl(i,t) * ml.middleCols(t*mr,mr);
        return res;
    }

    index_t rows() const { return _M.rows(); }
    index_t cols() const { return _M.rows(); } //_u.rows()

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _M.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    index_t cardinality_impl() const
    { return _u.cardinality(); }

    void print(std::ostream &os) const
    { os << "sum("; _M.print(os); os<<","; _u.print(os); os<<")"; }

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _M;

    mutable gsMatrix<Scalar> res;
};

/**
 * @brief
 * @todo this expression is not clear
 * @ingroup Expressions
 * @param u The first expression
 * @param M The second expression
 */
template <typename E1, typename E2> EIGEN_STRONG_INLINE
summ_expr<E1,E2> const summ(E1 const & u, E2 const& M)
{ return summ_expr<E1,E2>(u, M); }

}// namespace expr
}// namespace gismo