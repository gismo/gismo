/** @file adjugate_expr.h

    @brief Defines the adjugate expression

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
 * @brief Expression for the adjugate of a matrix
 * @ingroup Expressions
 * @tparam E The expression type
 */
template<class E>
class adjugate_expr  : public _expr<adjugate_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, ColBlocks = E::ColBlocks};
    enum {Space = E::Space};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;

public:
    adjugate_expr(_expr<E> const& u) : _u(u)
    {
        // gcc 4.8.4: invalid read due to _u.rows() using gsFuncData
        //GISMO_ASSERT(0== _u.cols()%_u.rows(), "Expecting square-block expression, got " << _u.rows() <<" x "<< _u.cols() );
    }

    // choose if ColBlocks
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto tmp = _u.eval(k);
        const index_t cb = _u.rows();
        const index_t r  = _u.cols() / cb;
        res.resize(_u.rows(),_u.cols());
        for (index_t i = 0; i!=r; ++i){
            res.middleCols(i*cb,cb) = tmp.middleCols(i*cb,cb).adjugate();
        }
        return res;
    }

    // choose if !ColBlocks
    //todo: Scalar eval(const index_t k) const

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    void print(std::ostream &os) const { os << "adj("; _u.print(os); os<<")"; }
};

}// namespace expr
}// namespace gismo