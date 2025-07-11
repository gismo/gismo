/** @file gsFeSpace.h

    @brief Defines a space expression

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
    @brief Null expression is a compatibility expression invalid at runtime
    @ingroup Expressions
    @tparam T The scalar type
*/
template<class T>
class gsNullExpr : public _expr<gsNullExpr<T> >
{
public:

    operator const gsFeSpace<T> & () const
    {
        static gsFeSpace<T> vv(-1);
        return vv;
    }

    typedef T Scalar;
    gsMatrix<T> eval(const index_t) const { GISMO_ERROR("gsNullExpr"); }
    inline index_t rows() const { GISMO_ERROR("gsNullExpr"); }
    inline index_t cols() const { GISMO_ERROR("gsNullExpr"); }
    void parse(gsExprHelper<T> &) const { }

    const gsFeSpace<T> & rowVar() const { GISMO_ERROR("gsNullExpr"); }
    const gsFeSpace<T> & colVar() const { GISMO_ERROR("gsNullExpr"); }

    void print(std::ostream &os) const { os << "NullExpr"; }

    static const gsNullExpr & get()
    {
        static gsNullExpr o;
        return o;
    }
//private:
    gsNullExpr() {}
};

}// namespace expr
}// namespace gismo