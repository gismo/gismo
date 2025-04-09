/** @file fform2nd_expr.h

    @brief Defines the second fundamental form expression

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
  Expression for the (precomputed) second fundamental form of a surface
*/

/**
 * @brief Expression for the second fundamental form of a geometry map,
 *        which is pre-computed by the \ref gsFuncData
 * @ingroup Expressions
 * @tparam T The scalar type
 */
template<class T>
class fform2nd_expr  : public _expr<fform2nd_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;
public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

    fform2nd_expr(const gsGeometryMap<T> & G) : _G(G) { }

    const gsAsConstMatrix<Scalar> eval(const index_t k) const
    {
        return gsAsConstMatrix<Scalar>(_G.data().fundForms.col(k).data(),rows(),cols());
    }

    index_t rows() const { return _G.source().domainDim() ; }
    index_t cols() const { return _G.source().domainDim() ; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_2ND_FFORM;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<T>::get();}

    void print(std::ostream &os) const { os << "fform2nd("; _G.print(os); os <<")"; }
};

/**
 * @brief Returns the second fundamental form of a geometry map
 * @ingroup Expressions
 * @param G The geometry map
 */
template<class T> EIGEN_STRONG_INLINE fform2nd_expr<T> fform2nd(const gsGeometryMap<T> & G)
{ return fform2nd_expr<T>(G); }

}// namespace expr
}// namespace gismo