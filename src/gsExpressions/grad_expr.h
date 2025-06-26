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
 * @brief Expression for the gradient of a constant value
 * @ingroup Expressions
 * @note  Transposed gradient vectors are returned as a matrix
 * @tparam E The expression type
 */
template<class T>
class grad_expr<_expr<T, true>> : public _expr<grad_expr<_expr<T, true>> >
{
public:
    enum {Space = 0, ScalarValued= 1, ColBlocks= 0}; // Order = E::Order+1

    typedef T Scalar;
    grad_expr(const _expr<T, true> & )
    {}

    static Scalar eval(const index_t k)
    {
        return Scalar(0.0);
    }

    index_t rows() const { return 0; }
    index_t cols() const { return 0; }

    void parse(gsExprHelper<Scalar> &) const { }

    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>::get(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>::get(); }

    void print(std::ostream &os) const { os << "0"; }
};

/**
 * @brief Gradient of an expression
 * @ingroup Expressions
 * @tparam E The expression type
 * @param u The expression
 */
template<class E> EIGEN_STRONG_INLINE
grad_expr<E> grad(const E & u) { return grad_expr<E>(u); }

template<class T> EIGEN_STRONG_INLINE
grad_expr<_expr<T,true>> grad(const _expr<T,true> & u) { return grad_expr<_expr<T,true>>(u); }


/**
 * @brief Gradient of an addition of two expressions
 * @ingroup Expressions
 */
template <typename E1, typename E2> EIGEN_STRONG_INLINE
add_expr<grad_expr<E1>, grad_expr<E2>>
grad(const add_expr<E1,E2> & e)
{ return grad(e.first()) + grad(e.second()); }

template <typename E1> EIGEN_STRONG_INLINE
grad_expr<E1>
grad(const add_expr<E1,_expr<typename E1::Scalar,true>> & e)
{ return grad(e.first()); }

template <typename E2> EIGEN_STRONG_INLINE
grad_expr<E2>
grad(const add_expr<_expr<typename E2::Scalar,true>,E2> & e)
{ return grad(e.second()); }

/**
 * @brief Gradient of a subtraction of two expressions
 * @ingroup Expressions
 */
template <typename E1, typename E2> EIGEN_STRONG_INLINE
sub_expr<grad_expr<E1>, grad_expr<E2>>
grad(const sub_expr<E1,E2> & e)
{ return grad(e.first()) - grad(e.second()); }

template <typename E1> EIGEN_STRONG_INLINE
grad_expr<E1>
grad(const sub_expr<E1,_expr<typename E1::Scalar,true>> & e)
{ return grad(e.first()); }

template <typename E2> EIGEN_STRONG_INLINE
sub_expr<_expr<typename E2::Scalar,true>,grad_expr<E2>>
grad(const sub_expr<_expr<typename E2::Scalar,true>,E2> & e)
{ return -grad(e.second()); }

/**
 * @brief Gradient of a multiplication of two expressions
 * @ingroup Expressions
 */
template <typename E1, typename E2> EIGEN_STRONG_INLINE
add_expr<mult_expr<E2,grad_expr<E1>>,mult_expr<E1,grad_expr<E2>>>
grad(const mult_expr<_expr<E1,false>,E2> & e)
{ return e.second() * grad(e.first()) + e.first() * grad(e.second()); }

template <typename E1> EIGEN_STRONG_INLINE
mult_expr<_expr<typename E1::Scalar,true>,grad_expr<E1>,false>
grad(const mult_expr<_expr<typename E1::Scalar,true>,E1,false> & e)
{ return e.first() * grad(e.second()); }

/**
 * @brief Gradient of an expression
 * @ingroup Expressions
 * @tparam E The expression type
 * @param u The expression
 */
template <typename E1, typename E2> EIGEN_STRONG_INLINE
divide_expr<sub_expr<mult_expr<grad_expr<E1>,E2>,mult_expr<E1,grad_expr<E2>>>,mult_expr<E2,E2>>
grad(const divide_expr<_expr<E1,false>,_expr<E2,false>> & e)
{ return (grad(e.first()) * e.second() - e.first() * grad(e.second())) / (e.second()*e.second()); }

template <typename E1> EIGEN_STRONG_INLINE
mult_expr<_expr<typename E1::Scalar,true>,grad_expr<E1> >
grad(const divide_expr<E1,_expr<typename E1::Scalar,true>> & e)
{ return ((typename E1::Scalar)1/e.second().val().eval()) * grad(e.first()); }

template <typename E2> EIGEN_STRONG_INLINE
mult_expr<_expr<typename E2::Scalar,true>,
          divide_expr<grad_expr<E2>,mult_expr<E2,E2> > >
grad(const divide_expr<_expr<typename E2::Scalar,true>,_expr<E2,false> > & e)
{ return (-e.first().val()) * ( grad(e.second()) / (e.second()*e.second()) ); }

// /**
//  * @brief Gradient of an expression
//  * @ingroup Expressions
//  * @tparam E The expression type
//  * @param u The expression
//  */
// template <typename E> EIGEN_STRONG_INLINE
// auto grad(const div_expr<E> & u) { return lapl(u); }


// /**
//  * @brief Gradient of an expression
//  * @ingroup Expressions
//  * @tparam E The expression type
//  * @param u The expression
//  */
// template <typename E> EIGEN_STRONG_INLINE
// auto div(const div_expr<E> & u) { GISMO_ERROR; }

// /**
//  * @brief Gradient of an expression
//  * @ingroup Expressions
//  * @tparam E The expression type
//  * @param u The expression
//  */
// template <typename E> EIGEN_STRONG_INLINE
// auto jac(const U & u) { return grad(u).tr(); }



}// namespace expr
}// namespace gismo
