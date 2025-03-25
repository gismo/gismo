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

/*
  Expression for addition operation
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

}// namespace expr
}// namespace gismo