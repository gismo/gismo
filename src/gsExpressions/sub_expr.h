/** @file sub_expr.h

    @brief Defines the sub expression

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
  Expression for subtraction operation
*/
template <typename E1, typename E2>
class sub_expr : public _expr<sub_expr<E1, E2> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = E1::ScalarValued && E2::ScalarValued,
        ColBlocks = E1::ColBlocks || E2::ColBlocks };
    enum {Space = E1::Space}; // == E2::Space

    typedef typename E1::Scalar Scalar;

    sub_expr(_expr<E1> const& u, _expr<E2> const& v)
    : _u(u), _v(v)
    {
        GISMO_ENSURE((int)E1::Space == (int)E2::Space &&
                     _u.rowVar()==_v.rowVar() && _u.colVar()==_v.colVar(),
                     "Error: substracting apples from oranges (use comma instead),"
                     " namely:\n" << _u <<"\n"<<_v);
    }

    mutable Temporary_t res;
    const Temporary_t & eval(const index_t k) const
    {
        // GISMO_ASSERT(_u.rowVar().id()==_v.rowVar().id() && _u.rowVar().isAcross()==_v.rowVar().isAcross(),
        //     "The row spaces are not split compatibly.");
        // GISMO_ASSERT(_u.colVar().id()==_v.colVar().id() && _u.colVar().isAcross()==_v.colVar().isAcross(),
        //     "The col spaces are not split compatibly.");
        GISMO_ASSERT(_u.rows() == _v.rows(),
                     "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in - operation:\n" << _u <<" minus \n" << _v );
        GISMO_ASSERT(_u.cols() == _v.cols(),
                     "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in - operation:\n" << _u <<" minus \n" << _v );
        GISMO_ASSERT(_u.cardinality() == _u.cardinality(),
                     "Cardinality "<< _u.cardinality()<<" != "<< _v.cardinality());
        //return (_u.eval(k) - _v.eval(k) ).eval();
        //return (_u.eval(k) - _v.eval(k) ); // any temporary matrices eval(.) will leak mem.
        res = _u.eval(k) - _v.eval(k);
        return res;
    }

    typename E1::Nested_t const & first() const { return _u; }
    typename E2::Nested_t const & second() const { return _v; }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    index_t cardinality_impl() const
    {
        GISMO_ASSERT(_u.cardinality() == _u.cardinality(),
                     "Cardinality "<< _u.cardinality()<<" != "<< _v.cardinality());
        return _u.cardinality();
    }

    void print(std::ostream &os) const
    { os << "("; _u.print(os); os<<" - ";_v.print(os); os << ")";}
};

/// Subtraction operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
sub_expr<E1,E2> const operator-(_expr<E1> const& u, _expr<E2> const& v)
{ return sub_expr<E1, E2>(u, v); }

template <typename E> EIGEN_STRONG_INLINE
sub_expr<_expr<typename E::Scalar>,E> const
operator-(typename E::Scalar const& s, _expr<E> const& v)
{ return sub_expr<_expr<typename E::Scalar>, E>(_expr<typename E::Scalar>(s), v); }

template <typename E> EIGEN_STRONG_INLINE
sub_expr<E, _expr<typename E::Scalar>> const
operator-(_expr<E> const& u, typename E::Scalar const& s)
{ return sub_expr<E, _expr<typename E::Scalar>>(u, _expr<typename E::Scalar>(s)); }

}// namespace expr
}// namespace gismo