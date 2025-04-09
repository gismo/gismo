/** @file colBlocks_expr.h

    @brief Defines the colBlocks expression

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
 * @brief Expression to make an expression colblocks
 * @ingroup Expressions
 * @tparam E The expression type
 */
template<class E>
class colBlocks_expr : public _expr<colBlocks_expr<E> >
{
    typename E::Nested_t _u;

public:

    typedef typename E::Scalar Scalar;

    colBlocks_expr(_expr<E> const& u)
    : _u(u) { }

public:
    enum {ColBlocks = 1, ScalarValued=E::ScalarValued};
    enum {Space = E::Space};

    mutable gsMatrix<Scalar> ev, res;

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        //return _u.eval(k).transpose();
        // /*
        ev = _u.eval(k);
        if (E::ColBlocks)
        {
            return ev;
        }
        else
        {
            res = _u.eval(k);

            GISMO_ASSERT(res.rows() % _u.rows() == 0 && res.cols() % _u.cols() == 0,"Result is not a multiple of the space dimensions?");

            index_t cardinality;
            if ( (cardinality = res.rows() / _u.rows()) >= 1 && res.cols() / _u.cols() == 1 ) // stored in rows
            {
                res.resize(_u.rows(), cardinality * _u.cols());
                for (index_t r = 0; r!=cardinality; r++)
                    res.block(0 , r * _u.cols(), _u.rows(), _u.cols()) = ev.block( r * _u.rows(), 0, _u.rows(), _u.cols() );
            }
            else if ( (cardinality = res.rows() / _u.rows()) == 1 && res.cols() / _u.cols() >= 1 ) // stored in cols ----->>>> This is already colBlocks???
            {
                res.resize(_u.rows(), cardinality * _u.cols());
                for (index_t r = 0; r!=cardinality; r++)
                    res.block(0 , r * _u.cols(), _u.rows(), _u.cols()) = ev.block( 0, r * _u.cols(), _u.rows(), _u.cols() );
            }
        }
        return res;
    }

    index_t rows() const { return _u.rows(); }

    index_t cols() const { return _u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    index_t cardinality_impl() const
    {
        res = _u.eval(0);
        index_t cardinality;
        if ( res.rows() / _u.rows() >= 1 && res.cols() / _u.cols() == 1 ) // stored in rows
            cardinality = res.rows() / _u.rows();
        else if ( res.rows() / _u.rows() == 1 && res.cols() / _u.cols() >= 1 )
            cardinality = res.cols() / _u.cols();
        else
            GISMO_ERROR("Cardinality for colBlocks_expr cannot be determined.");

        return cardinality;
    }

    void print(std::ostream &os) const { os<<"{"; _u.print(os); os <<"}"; }
};

}// namespace expr
}// namespace gismo