/** @file dJacG_expr.h

    @brief Defines the dJacG expression

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
 * @brief Expression for the partial derivative of the Jacobian of a geometry map
 * @ingroup Expressions
 * @tparam T The expression type
 */
template<class T>
class dJacG_expr : public _expr<dJacG_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

    mutable gsMatrix<T> res;
public:
    typedef T Scalar;

    dJacG_expr(const gsGeometryMap<T> & G) : _G(G) { }

    MatExprType eval(const index_t k) const
    {
        const index_t sz = _G.data().values[0].rows();
        const index_t s = _G.data().derivSize(); //dim.first*(_G.data().dim.first+1)/2;
        (void)s;
        res.resize(_G.data().dim.second, sz*_G.data().dim.first);
        res.setOnes();//_G.data().values[2].segment(i*k,k); // todo
        return res;
    }

    index_t rows() const { return _G.source().targetDim(); }
    index_t cols() const { return _G.source().domainDim(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_2ND_DER;
    }
};

/**
 * @brief Expression for the partial derivative of the Jacobian of a geometry map
 * @ingroup Expressions
 * @tparam T The expression type
 * @param G The geometry map
 */
template<class T> EIGEN_STRONG_INLINE
dJacG_expr<T> dJac(const gsGeometryMap<T> & G) { return dJacG_expr<T>(G); }

}// namespace expr
}// namespace gismo