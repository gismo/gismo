/** @file collapse_expr.h

    @brief Defines the collapse expression

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
 * @brief
 * @todo document
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 */
template <typename E1, typename E2>
class collapse_expr : public _expr<collapse_expr<E1, E2> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = 0, ColBlocks = 0};
    enum { Space = (int)E1::Space + (int)E2::Space };

    typedef typename E1::Scalar Scalar;

    mutable gsMatrix<Scalar> res;

    collapse_expr(_expr<E1> const& u,
                    _expr<E2> const& v)
    : _u(u), _v(v) { }

    //EIGEN_STRONG_INLINE MatExprType
    const gsMatrix<Scalar> &
    eval(const index_t k) const
    {
        const index_t nb = rows();
        const auto tmpA = _u.eval(k);
        const auto tmpB = _v.eval(k);

        if (E1::ColBlocks)
        {
            const index_t ur = _v.rows();
            res.resize(nb, ur);
            for (index_t i = 0; i!=nb; ++i)
            {
                res.row(i).transpose().noalias() = tmpA.middleCols(i*ur,ur) * tmpB;
            }
        }
        else if (E2::ColBlocks)
        {
            const index_t ur = _u.cols();
            res.resize(nb, ur);
            for (index_t i = 0; i!=nb; ++i)
            {
                res.row(i).noalias() = tmpA * tmpB.middleCols(i*ur,ur);
            }
        }

        return res;
    }

    index_t rows() const { return E1::ColBlocks ? _u.cols() / _v.rows() : _v.cols() / _u.cols() ; }
    index_t cols() const { return E1::ColBlocks ? _v.rows()  : _u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const
    { return E1::ColBlocks ? _u.rowVar() : _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {
        GISMO_ERROR("none");
    }

    void print(std::ostream &os) const { _u.print(os); os<<"~"; _v.print(os); }
};

// Multi-matrix collapsed by a vector
template <typename E1, typename E2> //EIGEN_STRONG_INLINE
//collapse_expr<E1,E2> const  operator&(<E1> const& u, _expr<E2> const& v)
collapse_expr<E1,E2> collapse( _expr<E1> const& u, _expr<E2> const& v)
{ return collapse_expr<E1, E2>(u, v); }

}// namespace expr
}// namespace gismo