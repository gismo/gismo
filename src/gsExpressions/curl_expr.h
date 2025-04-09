/** @file curl_expr.h

    @brief Defines the curl expression

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
 * @brief Expression for the curl of a vector field
 * @ingroup Expressions
 * @tparam T The scalar type
 * @todo finish this expression
 */
template<class T>
class curl_expr : public _expr<curl_expr<T> >
{
public:
    typedef T Scalar;
private:
    typename gsFeVariable<T>::Nested_t _u;
    mutable gsMatrix<Scalar> res;
public:
    enum{ Space = 1, ScalarValued= 0, ColBlocks= 0};

    curl_expr(const gsFeVariable<T> & u) : _u(u)
    { GISMO_ASSERT(3==u.dim(),"curl(.) requires 3D variable."); }

    const gsMatrix<T> & eval(const index_t k) const
    {
        res.setZero( rows(), _u.dim());
        const index_t na = _u.data().values[0].rows();
        gsAsConstMatrix<T, Dynamic, Dynamic> pd =
            _u.data().values[1].reshapeCol(k, cols(), na);

        res.col(0).segment(na  ,na) = -pd.row(2);
        res.col(0).segment(2*na,na) =  pd.row(1);
        res.col(1).segment(0   ,na) =  pd.row(2);
        res.col(1).segment(2*na,na) = -pd.row(0);
        res.col(2).segment(0   ,na) = -pd.row(1);
        res.col(2).segment(na  ,na) =  pd.row(0);
        return res;
    }

    index_t rows() const { return _u.dim() * _u.data().values[0].rows(); }
    index_t cols() const { return _u.data().dim.first; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_GRAD;
    }

    const gsFeSpace<T> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<T> & colVar() const {return gsNullExpr<T>::get();}

    void print(std::ostream &os) const { os << "curl("; _u.print(os); os <<")"; }
};

/// The curl of a finite element variable
template<class T> EIGEN_STRONG_INLINE
curl_expr<T> curl(const gsFeVariable<T> & u) { return curl_expr<T>(u); }

}// namespace expr
}// namespace gismo