/** @file flat_expr.h

    @brief Defines the flat expression

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

*/

/**
 * @brief Expression that transforms a matrix expression into a vector expression (Voigt-like)
 *
 *          Transforms a matrix expression into a vector expression by computing the vector
 *          [ a b c+d]^T
 *          for each matrix block
 *          [ a d ]
 *          [ c b ]
 *
 * @ingroup Expressions
 * @tparam E The expression type
 * @todo   Rename this to voigt_expr
 */
template<class E>
class flat_expr  : public _expr<flat_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, Space = E::Space, ColBlocks= 0}; // to do: ColBlocks
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> tmp;

public:

    flat_expr(_expr<E> const& u) : _u(u)
    {
        //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        tmp = _u.eval(k);
        const index_t numActives = _u.cardinality();

        for (index_t i = 0; i<numActives; ++i)
        {
            tmp(0,2*i+1) += tmp(1,2*i);
            std::swap(tmp(1,2*i), tmp(1,2*i+1));
        }

        tmp.resize(4,numActives);
        tmp.conservativeResize(3,numActives);

        if ( 1==Space )
            tmp.transposeInPlace();
        else if (2!=Space) // if not colSpan and not rowSpan
            tmp.transposeInPlace();

        return tmp;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "flat("; _u.print(os); os<<")"; }
};

/**
 * @brief Returns the flat expression of a matrix expression
 * @param u The expression
 * @ingroup Expressions
 */
template <typename E> EIGEN_STRONG_INLINE
flat_expr<E> const flat(E const & u)
{ return flat_expr<E>(u); }

}// namespace expr
}// namespace gismo