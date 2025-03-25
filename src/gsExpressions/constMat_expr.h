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
   Wrapper expression for constant matrices
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

}// namespace expr
}// namespace gismo