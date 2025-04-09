/** @file divide_expr.h

    @brief Defines the divide expression

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
  Expression for scalar division operation (first version)
*/

/**
 * @brief Expression for the division of two expressions (first version)
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 */
template <typename E1, typename E2>
class divide_expr : public _expr<divide_expr<E1,E2> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    typedef typename E1::Scalar Scalar;

public:
    enum {ScalarValued = E1::ScalarValued, ColBlocks= E2::ColBlocks};
    enum {Space = E1::Space}; // The denominator E2 has to be scalar.

    divide_expr(_expr<E1> const& u, _expr<E2> const& v)
    : _u(u), _v(v)
    {
        GISMO_STATIC_ASSERT(E2::ScalarValued, "The denominator needs to be scalar valued.");
    }

    AutoReturn_t eval(const index_t k) const
    { return ( _u.eval(k) / _v.eval(k) ); }

    typename E1::Nested_t const & first() const { return _u; }
    typename E2::Nested_t const & second() const { return _v; }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }


    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<" / ";_v.print(os);os << ")"; }
};

/**
 * @brief Expression for the division of an expression by a constant (second version)
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 *
 * This class is a partial specialization for the case where the second
 * expression is a constant value.
 */
template <typename E1>
class divide_expr<E1,typename E1::Scalar>
    : public _expr<divide_expr<E1,typename E1::Scalar> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    Scalar  const   _c;

public:
    enum {Space= E1::Space, ScalarValued = E1::ScalarValued, ColBlocks= E1::ColBlocks};

    divide_expr(_expr<E1> const& u, Scalar const  c)
    : _u(u), _c(c) { }

    AutoReturn_t eval(const index_t k) const
    { return ( _u.eval(k) / _c ); }


    typename E1::Nested_t const & first() const { return _u; }
                  Scalar const & second() const { return _c; }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }


    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<"/"<< _c << ")"; }
};

/**
 * @brief Expression for the division of a constant by an expression (third version)
 * @ingroup Expressions
 * @tparam E2 The second expression type
 *
 * This class is a partial specialization for the case where the first
 * expression is a constant value.
 */
template <typename E2>
class divide_expr<typename E2::Scalar,E2>
    : public _expr<divide_expr<typename E2::Scalar,E2> >
{
public:
    typedef typename E2::Scalar Scalar;

private:
    Scalar  const   _c;
    typename E2::Nested_t _u;
public:
    enum {Space= 0, ScalarValued = 1, ColBlocks= 0};

    divide_expr(Scalar const c, _expr<E2> const& u)
    : _c(c), _u(u)
    { GISMO_STATIC_ASSERT(E2::ScalarValued, "The denominator needs to be scalar valued."); }

    Scalar eval(const index_t k) const
    { return ( _c / _u.val().eval(k) ); }

                    Scalar const & first() const { return _c; }
    typename E2::Nested_t const & second() const { return _u; }

    index_t rows() const { return 0; }
    index_t cols() const { return 0; }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }


    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    void print(std::ostream &os) const
    { os << "("<< _c <<"/";_u.print(os);os << ")";}
};

/**
 * @brief Returns the division of two expressions
 * @param u The first expression
 * @param v The second expression
 * @ingroup Expressions
 */
template <typename E1, typename E2> EIGEN_STRONG_INLINE
divide_expr<E1,E2> const operator/(_expr<E1> const& u, _expr<E2> const& v)
{ return divide_expr<E1,E2>(u, v); }

/**
 * @brief Returns the division of an expression by a constant
 * @param u The expression
 * @param v The constant
 * @ingroup Expressions
 */
template <typename E> EIGEN_STRONG_INLINE
divide_expr<E,typename E::Scalar> const
operator/(_expr<E> const& u, const typename E::Scalar v)
{ return divide_expr<E,typename E::Scalar>(u, v); }

/**
 * @brief Returns the division of a constant by an expression
 * @param u The constant
 * @param v The expression
 * @ingroup Expressions
 */
template <typename E> EIGEN_STRONG_INLINE
divide_expr<typename E::Scalar,E> const
operator/(const typename E::Scalar u, _expr<E> const& v)
{ return divide_expr<typename E::Scalar,E>(u, v); }

}// namespace expr
}// namespace gismo