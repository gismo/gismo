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
 * @brief An example expression, for illustration purposes
 * @ingroup Expressions
 * @tparam E The expression type
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

    /// @brief Expression properties
    enum
    {
        ColBlocks = 1, /// @brief The expression is stored block-wise in columns
        ScalarValued=E::ScalarValued, /// @brief The expression is scalar valued
        Space = E::Space /// @brief The expression is a space (1: trail, 2: test, 3: combination)
    };

    /**
     * @brief Evaluates the expression at pre-computed point \a k
     */
    AutoReturn_t & eval(const index_t k) const
    {
        return 0.0;
    }

    /**
     * @brief Returns the number of rows of the expression
     */
    index_t rows() const { return _u.rows(); }

    /**
     * @brief Returns the number of columns of the expression
     */
    index_t cols() const { return _u.cols(); }

    /**
     * @brief Sets the flag for the expression
     */
    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    /**
     * @brief Returns the row variable of the expression
     * @todo  Elaborate more
     */
    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    /**
     * @brief Returns the column variable of the expression
     * @todo  Elaborate more
     */
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    /**
     * @brief Returns the cardinality of the expression, i.e. the number of basis functions at the point
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