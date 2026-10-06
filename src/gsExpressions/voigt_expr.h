/** @file voigt_expr.h

    @brief Defines the voigt expression

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
 * @brief Expression for the Voigt notation of a matrix expression
 *        Evaluates
 *       \f[
 *          \text{voigt}(u) = [ u_{11} u_{22} u_{12}+u_{21} ]^T
 *       \f]
 * @todo finish
 * @ingroup Expressions
 * @tparam E The type of the expression
 */
template<class E>
class voigt_expr  : public _expr<voigt_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, Space = E::Space, ColBlocks= 0}; // to do: ColBlocks
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> tmp;

public:

    voigt_expr(_expr<E> const& u) : _u(u)
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

}// namespace expr
}// namespace gismo