/** @file add_expr.h

    @brief Defines the add expression

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
 *  \brief Expression for the addition of two expressions
 *  \ingroup Expressions
 *  \tparam E1 The first expression type
 */
template <typename E1, typename E2>
class add_expr : public _expr<add_expr<E1, E2> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = E1::ScalarValued && E2::ScalarValued,
        ColBlocks = E1::ColBlocks || E2::ColBlocks };
    enum {Space = E1::Space}; // == E2::Space

    typedef typename E1::Scalar Scalar;

    add_expr(_expr<E1> const& u, _expr<E2> const& v)
    : _u(u), _v(v)
    {
        GISMO_ENSURE((int)E1::Space == (int)E2::Space &&
                     _u.rowVar()==_v.rowVar() && _u.colVar()==_v.colVar(),
                     "Error: adding apples and oranges (use comma instead),"
                     " namely:\n" << _u <<"\n"<<_v<<
                     " \nvars:\n" << _u.rowVar().id()<<"!="<<_v.rowVar().id() <<", "<< _u.colVar().id()<<"!="<<_v.colVar().id()<<
                     " \nspaces:\n" << (int)E1::Space<< "!="<< (int)E2::Space
            );
    }

    mutable Temporary_t res;
    const Temporary_t & eval(const index_t k) const
    {
        GISMO_ASSERT(_u.rows() == _v.rows(),
                     "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in + operation:\n"
                     << _u <<" plus \n" << _v );
        GISMO_ASSERT(_u.cols() == _v.cols(),
                     "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in + operation:\n"
                     << _u <<" plus \n" << _v );
        res = _u.eval(k) + _v.eval(k);
        return res;
    }

    typename E1::Nested_t const & first() const { return _u; }
    typename E2::Nested_t const & second() const { return _v; }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }


    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<" + ";_v.print(os);os << ")"; }
};

/**
 *  @brief Returns the sum of two expressions
 * @param u The first expression
 * @param v The second expression
 * @ingroup Expressions
 */
template <typename E1, typename E2> EIGEN_STRONG_INLINE
add_expr<E1,E2> const operator+(_expr<E1> const& u, _expr<E2> const& v)
{ return add_expr<E1, E2>(u, v); }

/**
 *  @brief Returns the sum of an expression and a scalar
 * @param u The expression
 * @param v The scalar
 * @ingroup Expressions
 */
template <typename E> EIGEN_STRONG_INLINE
add_expr< E, _expr<typename E::Scalar, true> >
operator+(_expr<E> const& u, const typename E::Scalar v)
{ return add_expr<E,_expr<typename E::Scalar>>(u, _expr<typename E::Scalar,true>(v)); }

/**
 *  @brief Returns the sum of a scalar and an expression
 * @param v The scalar
 * @param u The expression
 * @ingroup Expressions
 */
template <typename E> EIGEN_STRONG_INLINE
add_expr< E, _expr<typename E::Scalar, true> >
operator+(const typename E::Scalar v, _expr<E> const& u)
{ return add_expr<E,_expr<typename E::Scalar>>(u, _expr<typename E::Scalar,true>(v)); }

}// namespace expr
}// namespace gismo