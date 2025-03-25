/** @file grad_expr.h

    @brief Defines the grad expression

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

/*
  Expression for the gradient of a finite element variable

  Transposed gradient vectors are returned as a matrix
*/
template<class E>
class grad_expr : public _expr<grad_expr<E> >
{
    typename E::Nested_t _u;
public:
    enum {Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    typedef typename E::Scalar Scalar;
    mutable gsMatrix<Scalar> tmp;

    grad_expr(const E & u) : _u(u)
    { GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");}

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        // Assumes: derivatives are in _u.data().values[1]
        // gsExprHelper acounts for compositions/physical expressions
        // so that derivs are directly computed
        tmp = _u.data().values[1].reshapeCol(k, cols(), cardinality_impl()).transpose();
        return tmp;
    }

    index_t rows() const { return 1 /*==u.dim()*/; }

    index_t cols() const { return _u.source().domainDim(); }

    index_t cardinality_impl() const
    { return _u.data().values[1].rows() / cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_GRAD;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; _u.print(os); os <<")"; }
private:

    template<class U> static inline
    typename util::enable_if<util::is_same<U,gsComposition<Scalar> >::value,
                             const gsMatrix<Scalar> &>::type
    eval_impl(const U & u, const index_t k)
    {
        return u.eval(k);
    }
};

/*
  \brief Expression for the gradient of a finite element variable

  Transposed gradient vectors are returned as a matrix.
  This specialization is for a gsFeSolution object
*/

template<class T>
class grad_expr<gsFeSolution<T> > : public _expr<grad_expr<gsFeSolution<T> > >
{
protected:
    const gsFeSolution<T> _u;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    explicit grad_expr(const gsFeSolution<T> & u) : _u(u) { }

    mutable gsMatrix<T> res;
    const gsMatrix<T> & eval(index_t k) const
    {
        GISMO_ASSERT(_u.check(), "Invalid state in gsFeSolution");
        const gsDofMapper & map = _u.mapper();
        auto & act = _u.data().actives.col(1 == _u.data().actives.cols() ? 0:k );
        res.setZero(_u.dim(), _u.parDim());
        for (index_t c = 0; c!= _u.dim(); c++)
        {
            for (index_t i = 0; i!=_u.data().actives.rows(); ++i)
            {
                const index_t ii = map.index(act[i], _u.data().patchId, c);
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                {
                    res.row(c) += _u.coefs().at(ii) *
                        _u.data().values[1].col(k).segment(i*_u.parDim(), _u.parDim()).transpose();
                }
                else
                {
                    res.row(c) +=
                        _u.fixedPart().at( map.global_to_bindex(ii) ) *
                        _u.data().values[1].col(k).segment(i*_u.parDim(), _u.parDim()).transpose();
                }
            }
        }
        return res;
    }

    index_t rows() const {return _u.dim();}
    index_t cols() const {return _u.parDim(); }

    const gsFeSpace<Scalar> & rowVar() const
    {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList);                         // add symbol
        evList.add(_u.space());
        _u.data().flags |= NEED_GRAD|NEED_ACTIVE; // define flags
    }

    void print(std::ostream &os) const { os << "\u2207(s)"; }
};

}// namespace expr
}// namespace gismo