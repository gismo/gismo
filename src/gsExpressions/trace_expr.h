/** @file trace_expr.h

    @brief Defines the trace expression

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
 * @brief Expression for the trace of a (matrix) expression
 * @ingroup Expressions
 * @tparam E The type of the expression
 */
template<class E>
class trace_expr  : public _expr<trace_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, Space = E::Space, ColBlocks= 0};

private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;

public:
    trace_expr(_expr<E> const& u) : _u(u)
    {
        // gcc 4.8.4: invalid read due to _u.rows() using gsFuncData
        //GISMO_ASSERT(0== _u.cols()%_u.rows(), "Expecting square-block expression, got " << _u.rows() <<" x "<< _u.cols() );
    }

    // choose if ColBlocks
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto tmp = _u.eval(k);
        const index_t cb = _u.rows();
        const index_t r  = _u.cardinality();
        if (Space==1)
            res.resize(r, 1);
        else
            res.resize(1, r);

        for (index_t i = 0; i!=r; ++i)
            res.at(i) = tmp.middleCols(i*cb,cb).trace();
        return res;
    }

    // choose if !ColBlocks
    //todo: Scalar eval(const index_t k) const

    index_t rows() const { return _u.cols() / _u.rows(); } //_u.cardinality()?
    index_t cols() const { return 1; }

    index_t cardinality_impl() const { return _u.cardinality(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    void print(std::ostream &os) const { os << "trace("; _u.print(os); os<<")"; }
};

}// namespace expr
}// namespace gismo