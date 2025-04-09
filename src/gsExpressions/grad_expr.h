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

#include <gsExpressions/add_expr.h>
#include <gsExpressions/mult_expr.h>
#include <gsExpressions/sub_expr.h>
#include <gsExpressions/divide_expr.h>

namespace gismo
{
namespace expr
{

/**
 * @brief Expression for the gradient of a finite element variable
 * @ingroup Expressions
 * @note  Transposed gradient vectors are returned as a matrix
 * @tparam E The expression type
 */
template<class E>
class grad_expr : public _expr<grad_expr<E> >
{
    typename E::Nested_t _u;
public:
    enum {Space = E::Space, ScalarValued= 0, ColBlocks= 0}; // Order = E::Order+1

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

/**
 * @brief Expression for the gradient of a finite element solution
 * @ingroup Expressions
 * @note  Transposed gradient vectors are returned as a matrix
 * @tparam T The scalar type
 */
template<class T>
class grad_expr<gsFeSolution<T> > : public _expr<grad_expr<gsFeSolution<T> > >
{
protected:
    const gsFeSolution<T> _u;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0}; // ORDER IS DIFFICULT! Maybe make gsScalarSolution and gsVectorSolution?

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

/**
 * @brief Gradient of the addition of two expressions
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 */
template <typename E1, typename E2>
class grad_expr<add_expr<E1, E2> > : public _expr<grad_expr<add_expr<E1, E2> > >
{
    const typename E1::Nested_t _u;
    const typename E2::Nested_t _v;
    typedef add_expr<E1,E2> op_t;

public:
    enum {Space = E1::Space, ScalarValued= E1::ScalarValued, ColBlocks= 0};

    typedef typename E1::Scalar Scalar;
    mutable gsMatrix<Scalar> tmp;

    grad_expr(const op_t & expr)
    :
    _u(expr.first()),
    _v(expr.second())
    {
        // GISMO_ASSERT(E1::Space == E2::Space,"Error: grad(x+v) requires u and v to have the same space.");
        // GISMO_ASSERT()
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto expr = grad(_u) + grad(_v);
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
        grad(_u).parse(evList);
        grad(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; op_t::print(os); os <<")"; }
};

/**
 * @brief Gradient of the subtraction of two expressions
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 */
template <typename E1, typename E2>
class grad_expr<sub_expr<E1, E2> > : public _expr<grad_expr<sub_expr<E1, E2> > >
{
    const typename E1::Nested_t _u;
    const typename E2::Nested_t _v;
    typedef sub_expr<E1, E2> op_t;

public:
    enum {Space = E1::Space, ScalarValued= E1::ScalarValued, ColBlocks= 0};

    typedef typename E1::Scalar Scalar;
    mutable gsMatrix<Scalar> tmp;

    grad_expr(const sub_expr<E1, E2> & expr)
    :
    _u(expr.first()),
    _v(expr.second())
    {
        // GISMO_ASSERT(E1::Space == E2::Space,"Error: grad(x+v) requires u and v to have the same space.");
        // GISMO_ASSERT()
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto expr = grad(_u) - grad(_v);
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
        grad(_u).parse(evList);
        grad(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; op_t::print(os); os <<")"; }
};

/**
 * @brief Gradient of the multiplication of two expressions
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 */
// Implementation for a multiplication
template <typename E1, typename E2>
class grad_expr<mult_expr<E1, E2> > : public _expr<grad_expr<mult_expr<E1, E2> > >
{
    const typename E1::Nested_t _u;
    const typename E2::Nested_t _v;
    typedef mult_expr<E1, E2> op_t;

public:
    enum {Space = E1::Space, ScalarValued= E1::ScalarValued, ColBlocks= 0};

    typedef typename E1::Scalar Scalar;
    mutable gsMatrix<Scalar> tmp;

    grad_expr(const mult_expr<E1, E2> & expr)
    :
    _u(expr.first()),
    _v(expr.second())
    {
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto expr = _v * grad(_u) + _u * grad(_v);
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
        grad(_u).parse(evList);
        grad(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; op_t::print(os); os <<")"; }
};

/**
 * @brief Gradient of the multiplication of two expressions (scalar)
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 */
template <typename E2>
class grad_expr<mult_expr<typename E2::Scalar, E2,false> > : public _expr<grad_expr<mult_expr<typename E2::Scalar, E2,false> > >
{

public:
    enum {Space = E2::Space, ScalarValued= E2::ScalarValued, ColBlocks= 0};

    typedef typename E2::Scalar Scalar;

private:
    Scalar const _c;
    typename E2::Nested_t _v;
    typedef mult_expr<typename E2::Scalar, E2,false> op_t;

public:
    mutable gsMatrix<Scalar> tmp;

    grad_expr(const mult_expr<typename E2::Scalar, E2,false> & expr)
    :
    _c(expr.first()),
    _v(expr.second())
    {
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto expr = _c * grad(_v);
        tmp = expr.eval(k);
        return tmp;
    }

    index_t rows() const { return 1 /*==u.dim()*/; }
    index_t cols() const { return _v.source().domainDim(); }

    index_t cardinality_impl() const
    { return _v.data().values[1].rows() / cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        // _u.parse(evList);
        _v.parse(evList);
        // grad(_u).parse(evList);
        grad(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; op_t::print(os); os <<")"; }
};

/**
 * @brief Gradient of the division of two expressions
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 */
template <typename E1, typename E2>
class grad_expr<divide_expr<E1, E2> > : public _expr<grad_expr<divide_expr<E1, E2> > >
{
    const typename E1::Nested_t _u;
    const typename E2::Nested_t _v;
    typedef divide_expr<E1, E2> op_t;

public:
    enum {Space = E1::Space, ScalarValued= E1::ScalarValued, ColBlocks= 0};

    typedef typename E1::Scalar Scalar;
    mutable gsMatrix<Scalar> tmp;

    grad_expr(const divide_expr<E1, E2> & expr)
    :
    _u(expr.first()),
    _v(expr.second())
    {
        GISMO_ASSERT(E2::ScalarValued, "The denominator needs to be scalar valued.");
        GISMO_ASSERT(E2::Space == 0, "The gradient expression is not implemented for spaces in the denominator.");
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto expr = (grad(_u) * _v - _u * grad(_v)) / (_v * _v);
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
        grad(_u).parse(evList);
        grad(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; op_t::print(os); os <<")"; }
};

/**
 * @brief Gradient of the division of two expressions
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 */
template <typename E1>
class grad_expr<divide_expr<E1, typename E1::Scalar> > : public _expr<grad_expr<divide_expr<E1, typename E1::Scalar> > >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    const typename E1::Nested_t _u;
    Scalar  const   _c;
    typedef divide_expr<E1, typename E1::Scalar> op_t;

public:
    enum {Space = E1::Space, ScalarValued= E1::ScalarValued, ColBlocks= 0};

    mutable gsMatrix<Scalar> tmp;

    grad_expr(const divide_expr<E1, typename E1::Scalar> & expr)
    :
    _u(expr.first()),
    _c(expr.second())
    {
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto expr = grad(_u) / _c;
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
        // _v.parse(evList); // WHY NEEDED??
        grad(_u).parse(evList);
        // grad(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; op_t::print(os); os <<")"; }
};

/**
 * @brief Gradient of the division of two expressions
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 */
template <typename E2>
class grad_expr<divide_expr<typename E2::Scalar,E2> > : public _expr<grad_expr<divide_expr<typename E2::Scalar,E2 > > >
{
public:
    typedef typename E2::Scalar Scalar;

private:
    Scalar  const   _c;
    const typename E2::Nested_t _u;
    typedef divide_expr< typename E2::Scalar,E2 > op_t;

public:
    enum {Space = E2::Space, ScalarValued= E2::ScalarValued, ColBlocks= 0};

    mutable gsMatrix<Scalar> tmp;

    grad_expr(const divide_expr<typename E2::Scalar,E2> & expr)
    :
    _u(expr.first()),
    _c(expr.second())
    {
        GISMO_ASSERT(E2::ScalarValued, "The denominator needs to be scalar valued.");
        GISMO_ASSERT(E2::Space == 0, "The gradient expression is not implemented for spaces in the denominator.");
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        GISMO_NO_IMPLEMENTATION;
        // auto expr = grad(_u) / _c;
        // tmp = expr.eval(k);
        return tmp;
    }

    index_t rows() const { return 1 /*==u.dim()*/; }
    index_t cols() const { return _u.source().domainDim(); }

    index_t cardinality_impl() const
    { return _u.data().values[1].rows() / cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList); // WHY NEEDED??
        // _v.parse(evList); // WHY NEEDED??
        grad(_u).parse(evList);
        // grad(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; op_t::print(os); os <<")"; }
};

/**
 * @brief Gradient of the exponentiation of an expression
 * @ingroup Expressions
 * @tparam E The expression type
 * @param u The expression
 */
template<class E> EIGEN_STRONG_INLINE
grad_expr<E> grad(const E & u) { return grad_expr<E>(u); }

}// namespace expr
}// namespace gismo