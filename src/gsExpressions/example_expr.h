/** @file example_expr.h

    @brief This file provides a template for expressions

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
 */
template<class E>
class example_expr : public _expr<example_expr<E> >
{
    typename E::Nested_t _u;

public:

    typedef typename E::Scalar Scalar;

    example_expr(_expr<E> const& u)
    : _u(u) { }

public:

    enum
    {
        ColBlocks = 1,
        ScalarValued=E::ScalarValued,
        Space = E::Space
    };

    AutoReturn_t & eval(const index_t k) const
    {
        return 0.0;
    }

    /**
     * @brief
     */
    index_t rows() const { return _u.rows(); }

    /**
     * @brief
     */
    index_t cols() const { return _u.cols(); }

    /**
     * @brief
     */
    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    /**
     * @brief
     */
    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    /**
     * @brief
     */
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    /**
     * @brief
     */
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

    /**
     * @brief
     */
    void print(std::ostream &os) const { os<<"{"; _u.print(os); os <<"}"; }
};

}// namespace expr
}// namespace gismo