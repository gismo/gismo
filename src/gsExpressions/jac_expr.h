/** @file jac_expr.h

    @brief Defines the jacobian expression

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
               H.M. Verhelst
*/

#include <gsExpressions/jacInv_expr.h>

#pragma once

namespace gismo
{
namespace expr
{

/*
  Expression for the Jacobian matrix of a FE variable
*/
template<class E>
class jac_expr : public _expr<jac_expr<E> >
{
    typename E::Nested_t _u;
public:
    enum {ColBlocks = (1==E::Space?1:0) };
    enum {Space = E::Space, ScalarValued= 0 };

    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;

    jac_expr(const E & _u_)
    : _u(_u_) { }

    MatExprType eval(const index_t k) const
    {
        if (0!=Space)
        {
            // Dim x (numActive*Dim)
            res = _u.data().values[1].col(k).transpose().blockDiag(_u.dim());
        }
        else
        {
            res = _u.data().values[1]
                .reshapeCol(k, _u.parDim(), _u.targetDim()).transpose()
                .blockDiag(_u.dim());
        }
        return res;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    //const gsFeSpace<Scalar> & rowVar() const { return rowVar_impl<E>(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    index_t rows() const { return rows_impl(_u); }
    index_t cols() const { return cols_impl(_u); }

    // index_t rows() const { return _u.dim(); }
    // index_t cols() const { return _u.source().domainDim(); }

    index_t cardinality_impl() const
    {
        return _u.dim() * _u.data().actives.rows();
    }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_DERIV|NEED_ACTIVE;
        //note: cardinality() depends on actives
    }

    void print(std::ostream &os) const { os << "\u2207("; _u.print(os);os <<")"; }

private:

    // The jacobian is different for gsFeVariable, gsFeSolution and gsFeSpace
    // gsFeSolution: Does not work
    // gsFeVariable: dim()=1 and source().targetDim()=d
    // gsFeSpace: dim()=d and source().targetDim()=1
    template<class U> inline
    typename util::enable_if<!(util::is_same<U,gsFeSpace<Scalar> >::value), index_t >::type  // What about solution??
    rows_impl(const U & u)  const
    {
        return u.source().targetDim();
    }

    template<class U> inline
    typename util::enable_if< (util::is_same<U,gsFeSpace<Scalar> >::value), index_t >::type
    rows_impl(const U & u) const
    {
        return u.dim();
    }

    template<class U> inline
    typename util::enable_if<!(util::is_same<U,gsFeSpace<Scalar> >::value), index_t >::type
    cols_impl(const U & u)  const
    {
        return u.source().domainDim();
    }

    template<class U> inline
    typename util::enable_if< (util::is_same<U,gsFeSpace<Scalar> >::value), index_t >::type
    cols_impl(const U & u) const
    {
        return u.source().domainDim();
    }

    // The jacobian is different for gsFeVariable, gsFeSolution and gsFeSpace
    // gsFeSolution: Does not work
    // gsFeVariable: rowVar = NULL
    // gsFeSpace: rowVar = u
    template<class U> inline
    typename util::enable_if<!(util::is_same<U,gsFeSpace<Scalar> >::value), const gsFeSpace<Scalar> &  >::type
    rowVar_impl() const
    {
        return gsNullExpr<Scalar>::get();
    }

    template<class U> inline
    typename util::enable_if<(util::is_same<U,gsFeSpace<Scalar> >::value), const gsFeSpace<Scalar> &  >::type
    rowVar_impl() const
    {
        return _u;
    }
};

/*
  Expression for the Jacobian matrix of a geometry map
*/
template<class T>
class jac_expr<gsGeometryMap<T> > : public _expr<jac_expr<gsGeometryMap<T> > >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    jac_expr(const gsGeometryMap<T> & G) : _G(G) { }
    MatExprType eval(const index_t k) const
    {
        // TarDim x ParDim
        return _G.data().values[1]
            .reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
    }

    index_t rows() const { return _G.source().targetDim(); }

    index_t cols() const { return _G.source().domainDim(); }

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_DERIV;
    }

    meas_expr<T> absDet() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return meas_expr<T>(_G);
    }

    jacInv_expr<T> inv() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return jacInv_expr<T>(_G);
    }

    /// The generalized Jacobian matrix inverse, i.e.: (J^t J)^{-t} J^t
    jacInv_expr<T> ginv() const { return jacInv_expr<T>(_G); }

    void print(std::ostream &os) const { os << "\u2207("; _G.print(os); os <<")"; }
};

/// The Jacobian matrix of a FE variable
template<class E> EIGEN_STRONG_INLINE
jac_expr<E> jac(const symbol_expr<E> & u) { return jac_expr<E>(u); }

/// The Jacobian matrix of a geometry map
template<class T> EIGEN_STRONG_INLINE
jac_expr<gsGeometryMap<T> > jac(const gsGeometryMap<T> & G) {return jac_expr<gsGeometryMap<T> >(G);}

/// Jacobian matrix for a solution expression
template<class T> EIGEN_STRONG_INLINE
grad_expr<gsFeSolution<T> > jac(const gsFeSolution<T> & s) {return grad_expr<gsFeSolution<T> >(s);}

}// namespace expr
}// namespace gismo