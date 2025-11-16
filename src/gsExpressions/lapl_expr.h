/** @file lapl_expr.h

    @brief Defines the lapl expression

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
 * @brief Expression for the Laplacian of a finite element variable
 * @ingroup Expressions
 * @tparam E The expression type
 */
template<class E>
class lapl_expr : public _expr<lapl_expr<E> >
{
    typename E::Nested_t _u;

public:
    typedef typename E::Scalar Scalar;
    enum {Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    lapl_expr(const E & u) : _u(u) { }

    auto eval(const index_t k) const -> decltype(_u.data().laplacians.col(k))
    {
        // numActive x 1
        return _u.data().laplacians.col(k);
        //todo: replace by
        // NEED_DERIV2
        // ..nabla2.sum()
    }

    index_t rows() const { return _u.data().laplacians.rows(); }
    index_t cols() const { return 1; }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_LAPLACIAN;
    }

    static const gsFeSpace<Scalar> & rowVar() {return E::rowVar();}
    static const gsFeSpace<Scalar> & colVar() {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2206("; _u.print(os); os <<")"; } //or \u0394
};

/**
 * @brief Expression for the Laplacian of a finite element solution
 * @ingroup Expressions
 * @tparam T The expression type
 */
template<class T>
class lapl_expr<gsFeSolution<T> > : public _expr<lapl_expr<gsFeSolution<T> > >
{
protected:
    const gsFeSolution<T> _u;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    lapl_expr(const gsFeSolution<T> & u) : _u(u) { }

    mutable gsMatrix<T> res;
    const gsMatrix<T> & eval(const index_t k) const
    {
        GISMO_ASSERT(_u.check(), "Invalid state in gsFeSolution");
        const gsDofMapper & map = _u.mapper();
        res.setZero(_u.dim(), 1); //  scalar, but per component

        index_t numActs = _u.data().values[0].rows();
        index_t numDers = _u.parDim() * (_u.parDim() + 1) / 2;
        gsMatrix<T> deriv2;

        auto & act = _u.data().actives.col(1 == _u.data().actives.cols() ? 0:k );
        for (index_t c = 0; c!= _u.dim(); c++)
            for (index_t i = 0; i!=numActs; ++i)
            {
                const index_t ii = map.index(act[i], _u.data().patchId, c);
                deriv2 = _u.data().values[2].block(i*numDers,k,_u.parDim(),1); // this only takes d11, d22, d33 part. For all the derivatives [d11, d22, d33, d12, d13, d23]: col.block(i*numDers,k,numDers,1)
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                    res.at(c) += _u.coefs().at(ii) * deriv2.sum();
                else
                    res.at(c) +=_u.fixedPart().at( map.global_to_bindex(ii) ) * deriv2.sum();
            }
        return res;
    }

    index_t rows() const { return _u.dim(); }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u.space());
        _u.data().flags |= NEED_ACTIVE | NEED_DERIV2;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<T>::get();}

    void print(std::ostream &os) const { os << "\u2206(s)"; }
};

/**
 * @brief Returns the Laplacian of an addition
 * @ingroup Expressions
 */


/**
 * @brief Returns the Laplacian of an expression
 * @ingroup Expressions
 * @param u The expression
 */
template<class E> EIGEN_STRONG_INLINE
lapl_expr<E> lapl(const symbol_expr<E> & u) { return lapl_expr<E>(u); }

/**
 * @brief Returns the Laplacian of a geometry map
 * @ingroup Expressions
 * @param G The geometry map
 */
template<class T> EIGEN_STRONG_INLINE
lapl_expr<gsFeSolution<T> > lapl(const gsFeSolution<T> & u)
{ return lapl_expr<gsFeSolution<T> >(u); }

}// namespace expr
}// namespace gismo