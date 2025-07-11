/** @file idMat_expr.h

    @brief Defines the idMat expression

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
 * @brief Expression for the identity matrix
 * @ingroup Expressions
 * @tparam E The expression type
 */
class idMat_expr : public _expr<idMat_expr >
{
public:
    typedef real_t Scalar;
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};
private:
    index_t _dim;

public:
    idMat_expr(const index_t dim) : _dim(dim) { }

public:

    gsMatrix<Scalar>::IdentityReturnType eval(const index_t) const
    {
        return gsMatrix<Scalar>::Identity(_dim,_dim);
    }

    index_t rows() const { return _dim; }
    index_t cols() const { return  _dim; }
    void parse(gsExprHelper<Scalar> & ) const {  }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "id("<<_dim <<")";}
};

/**
 * @brief Function to create an identity matrix expression
 * @ingroup Expressions
 * @param dim The dimension of the identity matrix
 */
EIGEN_STRONG_INLINE idMat_expr id(const index_t dim) { return idMat_expr(dim); }

}// namespace expr
}// namespace gismo