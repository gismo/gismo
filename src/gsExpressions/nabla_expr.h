/** @file nabla_expr.h

    @brief Defines the nabla expression

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
 * @brief Expression for the nabla (\f$\nabla\f$) of a finite element variable
 * @ingroup Expressions
 * @tparam T The type of the expression
 */
template<class T>
class nabla_expr : public _expr<nabla_expr<T> >
{
    typename gsFeVariable<T>::Nested_t u;

public:
    typedef T Scalar;
    enum{Space = 1};

    /* // todo
       nabla_expr(const gsGeometryMap<T> & G)
       : m_data(G.data()) { }
    */

    nabla_expr(const gsFeVariable<T> & _u) : u(_u)
    {
        //GISMO_ASSERT(u.parDim()==u.dim(),"nabla(.) requires tarDim==parDim:"
        //             << u.parDim() <<"!="<< u.dim() <<"\n" );
    }

    mutable gsMatrix<Scalar> res;

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t d = u.cols();
        const index_t n = rows() / d;
        res.setZero(rows(), d);

        for (index_t i = 0; i!=d; ++i)
            res.col(i).segment(i*n,n) = u.data().values[1].reshapeCol(k, d, n).row(i);
        return res;
    }

    index_t rows() const { return u.rows(); }
    index_t cols() const { return u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(u);
        u.data().flags |= NEED_GRAD;
    }

    const gsFeSpace<T> & rowVar() const { return u.rowVar(); }
    const gsFeSpace<T> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "nabla("; u.print(os); os <<")"; }
};

/**
 * @brief nabla operator for finite element variables
 * @ingroup Expressions
 * @param u The finite element variable
 */
template<class T> EIGEN_STRONG_INLINE
nabla_expr<T> nabla(const gsFeVariable<T> & u) { return nabla_expr<T>(u); }

}// namespace expr
}// namespace gismo