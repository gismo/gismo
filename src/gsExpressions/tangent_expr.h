/** @file tangent_expr.h

    @brief Defines the tangent expression

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
 * @brief Expression for the tangent vector of a geometry map
 *        This expression is valid only at the boundaries of a geometric patch
 * @ingroup Expressions
 * @tparam T The type of the expression
 */
template<class T>
class tangent_expr : public _expr<tangent_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    tangent_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsVector<Scalar> res;
    const gsVector<Scalar> & eval(const index_t k) const
    {
        if (_G.targetDim()==2)
        {
            res = _G.data().outNormals.col(k);//2x1
            std::swap( res(0,0), res(1,0) );
            res(0,0) *= -1;
            return res;
        }
        else if (_G.targetDim()==3)
        {
            res.resize(3);
            res.col3d(0) = _G.data().normals.col3d(k)
                .cross( _G.data().outNormals.col3d(k) );
            return res;
        }
        else
            GISMO_ERROR("Function not implemented for dimension"<<_G.targetDim());

    }

    index_t rows() const { return _G.source().targetDim(); }
    index_t cols() const { return 1; }

    static const gsFeSpace<Scalar> & rowVar() {return gsNullExpr<Scalar>::get();}
    static const gsFeSpace<Scalar> & colVar() {return gsNullExpr<Scalar>::get();}

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL;
        _G.data().flags |= NEED_OUTER_NORMAL;
    }

    void print(std::ostream &os) const { os << "tv("; _G.print(os); os <<")"; }
};

/**
 * @brief Expression for the tangent vector of a geometry map
 * @ingroup Expressions
 * @tparam T The type of the expression
 */
template<class T> EIGEN_STRONG_INLINE
tangent_expr<T> tv(const gsGeometryMap<T> & u) { return tangent_expr<T>(u); }

}// namespace expr
}// namespace gismo