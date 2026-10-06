/** @file matrix_by_space_tr_expr.h

    @brief Defines the matrix_by_space expression

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
 * @brief Expression for computing the outer products of a matrix by a space of dimension > 1
 * @ingroup Expressions
 * @tparam E1 The first expression type
 * @tparam E2 The second expression type
 */
template <typename E1, typename E2>
class matrix_by_space_tr_expr  : public _expr<matrix_by_space_tr_expr<E1,E2> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0, ColBlocks = 1};
    enum {Space = E2::Space};

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    mutable gsMatrix<Scalar> res;

public:
    matrix_by_space_tr_expr(E1 const& u, E2 const& v) : _u(u), _v(v) { }


    // choose if ColBlocks
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t r  = _u.rows();
        const index_t N  = _v.cols() / (r*r);

        const auto uEv  = _u.eval(k);
        const auto vEv  = _v.eval(k);

        res.resize(r, N*r*r);
        for (index_t s = 0; s!=r; ++s)
            for (index_t i = 0; i!=N; ++i)
            {
                res.middleCols((s*N + i)*r,r).noalias() =
                    uEv.transpose()*vEv.middleCols((s*N + i)*r,r).transpose();
            }
        //meaning: [Jg Jg Jg] * Jb ..
        return res;
    }

    index_t rows() const { return _u.cols(); }
    index_t cols() const { return _v.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const { os << "matrix_by_space_tr("; _u.print(os); os<<")"; }
};

/**
 * @brief Matrix by space expression
 * @ingroup Expressions
 * @tparam E1 The first expression
 * @tparam E2 The second expression
 * @param u The first expression
 * @param v The second expression
 * @note Matrix by space TODO: find better name and/or description? And is this the best place?
 */
template <typename E1, typename E2> EIGEN_STRONG_INLINE
matrix_by_space_tr_expr<E1,E2> const matrix_by_space_tr(E1 const & u, E2 const& v)
{ return matrix_by_space_tr_expr<E1,E2>(u, v); }


}// namespace expr
}// namespace gismo