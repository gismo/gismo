/** @file max_expr.h

    @brief Defines the max expression

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
 * @brief Expression for the max of a vector
 * @ingroup Expressions
 * @tparam E The expression type
 */
template<class E>
class max_expr  : public _expr<max_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, Space = E::Space, ColBlocks = 1};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> tmp;
    mutable gsMatrix<Scalar> res;

public:

    max_expr(_expr<E> const& u) : _u(u)
    {
        //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
    }

    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,k); }

    index_t rows() const { return 1; }
    index_t cols() const { return 1; }
    void setFlag() const { _u.setFlag(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "max("; _u.print(os); os<<")"; }
private:
    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        tmp = u.eval(k);

        res.resize(1,u.cardinality());
        if (E::ColBlocks)
            for (index_t c=0; c!=_u.cardinality(); c++)
                res(0,c) = tmp.block(0,c*u.cols(),u.rows(),u.cols()).maxCoeff();
        else
            for (index_t c=0; c!=_u.rows(); c++)
                res(0,c) = tmp.block(c*u.rows(),0,u.rows(),u.cols()).maxCoeff();
        return res;
    }


    template<class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        res = u.eval(k).colwise().maxCoeff();
        return res;
    }
};

}// namespace expr
}// namespace gismo