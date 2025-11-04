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
#include <gsExpressions/add_expr.h>
#include <gsExpressions/mult_expr.h>
#include <gsExpressions/sub_expr.h>
#include <gsExpressions/divide_expr.h>

#pragma once

namespace gismo
{
namespace expr
{

/**
 * @brief Expression for the Jacobian matrix of a finite element variable
 * @ingroup Expressions
 * @tparam E The expression type
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

/**
 * @brief Expression for the Jacobian matrix of a geometry map
 * @ingroup Expressions
 * @tparam T The expression type
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

/**
 * @brief Jacobian matrix for an addition
 * @ingroup Expressions
 * @tparam E1 The first expression
 * @tparam E2 The second expression
 * @details The Jacobian of the addition is the sum of the Jacobians.
 */
template <typename E1, typename E2>
class jac_expr<add_expr<E1, E2> > : public _expr<jac_expr<add_expr<E1, E2> > >
{
    const typename E1::Nested_t _u;
    const typename E2::Nested_t _v;

public:
    enum {Space = E1::Space, ScalarValued= E1::ScalarValued, ColBlocks= 0};

    typedef typename E1::Scalar Scalar;
    mutable gsMatrix<Scalar> uVals, uGrads, vVals, vGrads, tmp;

    jac_expr(const add_expr<E1, E2> & expr)
    :
    _u(expr.first()),
    _v(expr.second())
    {
        // GISMO_ASSERT(E1::Space == E2::Space,"Error: grad(x+v) requires u and v to have the same space.");
        // GISMO_ASSERT()
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto expr = jac(_u) + jac(_v);
        tmp = expr.eval(k);
        return tmp;
    }

    index_t rows() const { return 1 /*==u.dim()*/; }
    index_t cols() const { return _u.source().domainDim(); }

    index_t cardinality_impl() const
    { return _u.data().values[1].rows() / cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList); // WHY NEEDED??
        _v.parse(evList); // WHY NEEDED??
        jac(_u).parse(evList);
        jac(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; _u.print(os); os <<")"; }
};

/**
 * @brief Jacobian of the subtraction of two expressions
 * @ingroup Expressions
 * @tparam E1 The first expression
 * @tparam E2 The second expression
 * @details The Jacobian of the subtraction is the difference of the Jacobians.
 */
template <typename E1, typename E2>
class jac_expr<sub_expr<E1, E2> > : public _expr<jac_expr<sub_expr<E1, E2> > >
{
    const typename E1::Nested_t _u;
    const typename E2::Nested_t _v;

public:
    enum {Space = E1::Space, ScalarValued= E1::ScalarValued, ColBlocks= 0};

    typedef typename E1::Scalar Scalar;
    mutable gsMatrix<Scalar> uVals, uGrads, vVals, vGrads, tmp;

    jac_expr(const sub_expr<E1, E2> & expr)
    :
    _u(expr.first()),
    _v(expr.second())
    {
        // GISMO_ASSERT(E1::Space == E2::Space,"Error: grad(x+v) requires u and v to have the same space.");
        // GISMO_ASSERT()
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto expr = jac(_u) - jac(_v);
        tmp = expr.eval(k);
        return tmp;
    }

    index_t rows() const { return 1 /*==u.dim()*/; }
    index_t cols() const { return _u.source().domainDim(); }

    index_t cardinality_impl() const
    { return _u.data().values[1].rows() / cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList); // WHY NEEDED??
        _v.parse(evList); // WHY NEEDED??
        jac(_u).parse(evList);
        jac(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; _u.print(os); os <<")"; }
};

/**
 * @brief Gradient of the multiplication of two expressions
 * @ingroup Expressions
 * @tparam E1 The first expression
 * @tparam E2 The second expression
 * @details The Jacobian of the multiplication is the product rule:
 * \f[
 * \nabla (u v) = u \nabla v + v \nabla u
 * \f]
 */
template <typename E1, typename E2>
class jac_expr<mult_expr<E1, E2> > : public _expr<jac_expr<mult_expr<E1, E2> > >
{
    const typename E1::Nested_t _u;
    const typename E2::Nested_t _v;

public:
    enum {Space = E1::Space, ScalarValued= E1::ScalarValued, ColBlocks= 0};

    typedef typename E1::Scalar Scalar;
    mutable gsMatrix<Scalar> uVals, uGrads, vVals, vGrads, tmp;

    jac_expr(const mult_expr<E1, E2> & expr)
    :
    _u(expr.first()),
    _v(expr.second())
    {
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto expr = _v * jac(_u) + _u * jac(_v);
        tmp = expr.eval(k);
        return tmp;
    }

    index_t rows() const { return 1 /*==u.dim()*/; }
    index_t cols() const { return _u.source().domainDim(); }

    index_t cardinality_impl() const
    { return _u.data().values[1].rows() / cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList);
        _v.parse(evList);
        jac(_u).parse(evList);
        jac(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; _u.print(os); os <<")"; }
};

/**
 * @brief Gradient of the division of two expressions
 * @ingroup Expressions
 * @tparam E1 The first expression
 * @tparam E2 The second expression
 * @details The Jacobian of the division is the quotient rule:
 * \f[
 * \nabla (u / v) = \frac{v \nabla u - u \nabla v}{v^2}
 * \f]
 */
template <typename E1, typename E2>
class jac_expr<divide_expr<E1, E2> > : public _expr<jac_expr<divide_expr<E1, E2> > >
{
    const typename E1::Nested_t _u;
    const typename E2::Nested_t _v;

public:
    enum {Space = E1::Space, ScalarValued= E1::ScalarValued, ColBlocks= 0};

    typedef typename E1::Scalar Scalar;
    mutable gsMatrix<Scalar> uVals, uGrads, vVals, vGrads, tmp;

    jac_expr(const divide_expr<E1, E2> & expr)
    :
    _u(expr.first()),
    _v(expr.second())
    {
        GISMO_ASSERT(E2::ScalarValued, "The denominator needs to be scalar valued.");
        GISMO_ASSERT(E2::Space == 0, "The gradient expression is not implemented for spaces in the denominator.");
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto expr = (jac(_u) * _v - _u * jac(_v)) / (_v * _v);
        tmp = expr.eval(k);
        return tmp;
    }

    index_t rows() const { return 1 /*==u.dim()*/; }
    index_t cols() const { return _u.source().domainDim(); }

    index_t cardinality_impl() const
    { return _u.data().values[1].rows() / cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList); // WHY NEEDED??
        _v.parse(evList); // WHY NEEDED??
        jac(_u).parse(evList);
        jac(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; _u.print(os); os <<")"; }
};

/// The Jacobian matrix of a FE variable

/**
 * @brief Jacobian matrix of a finite element variable
 * @ingroup Expressions
 * @tparam E The expression type
 * @param u The finite element variable
 */
template<class E> EIGEN_STRONG_INLINE
jac_expr<E> jac(const symbol_expr<E> & u) { return jac_expr<E>(u); }

/**
 * @brief Jacobian matrix of a geometry map
 * @ingroup Expressions
 * @tparam T The expression type
 * @param G The geometry map
 */
template<class T> EIGEN_STRONG_INLINE
jac_expr<gsGeometryMap<T> > jac(const gsGeometryMap<T> & G) {return jac_expr<gsGeometryMap<T> >(G);}

/**
 * @brief Jacobian matrix of a finite element solution
 * @ingroup Expressions
 * @tparam T The expression type
 * @param s The finite element solution
 */
template<class T> EIGEN_STRONG_INLINE
grad_expr<gsFeSolution<T> > jac(const gsFeSolution<T> & s) {return grad_expr<gsFeSolution<T> >(s);}

}// namespace expr
}// namespace gismo