/** @file constMat_expr.h

    @brief Defines the constMat expression

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
 * @brief Expression for a constant matrix
 * @ingroup Expressions
 * @tparam E The expression type
 */
class constMat_expr : public _expr<constMat_expr >
{
public:
    typedef real_t Scalar;
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};
private:
    gsMatrix<Scalar> _mat;

public:
    constMat_expr(const gsMatrix<Scalar> mat) : _mat(mat) { }

public:

    gsMatrix<Scalar> eval(const index_t) const
    {
        return _mat;
    }

    index_t rows() const { return _mat.rows(); }
    index_t cols() const { return  _mat.cols(); }
    void parse(gsExprHelper<Scalar> & ) const {  }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "constMat";}
};

/**
 * @brief Expression for a constant matrix of ones
 * @ingroup Expressions
 * @param dim The dimension of the matrix
 * @return A constant matrix of ones
 */
EIGEN_STRONG_INLINE constMat_expr ones(const index_t dim)
{
    gsMatrix<real_t> ones(dim, dim);
    ones.fill(1);
    return constMat_expr(ones);
}

/**
 * @brief Expression that turns a matrix into an expression
 * @ingroup Expressions
 * @param mat The matrix
 * @return A constant matrix of zeros
 */
EIGEN_STRONG_INLINE constMat_expr mat(const gsMatrix<real_t> mat) { return constMat_expr(mat); }


}// namespace expr
}// namespace gismo