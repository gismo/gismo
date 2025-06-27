/** @file hess_expr.h

    @brief Defines the hessian expression

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

// Adaptor to compute Hessian
template <typename Derived>
void secDerToHessian(const gsEigen::DenseBase<Derived> &  secDers,
                     const index_t dim,
                     gsMatrix<typename Derived::Scalar> & hessian)
{
    const index_t sz = dim*(dim+1)/2;
    const gsAsConstMatrix<typename Derived::Scalar>
        ders(secDers.derived().data(), sz, secDers.size() / sz );
    hessian.resize(dim*dim, ders.cols() );

    switch ( dim )
    {
    case 1:
        hessian = secDers.transpose(); //==ders
        break;
    case 2:
        hessian.row(0)=ders.row(0);//0,0
        hessian.row(1)=//1,0
            hessian.row(2)=ders.row(2);//0,1
        hessian.row(3)=ders.row(1);//1,1
        break;
    case 3:
        hessian.row(0)=ders.row(0);//0,0
        hessian.row(3)=//0,1
            hessian.row(1)=ders.row(3);//1,0
        hessian.row(6)=//0,2
            hessian.row(2)=ders.row(4);//2,0
        hessian.row(4)=ders.row(1);//1,1
        hessian.row(7)=//1,2
            hessian.row(5)=ders.row(5);//2,1
        hessian.row(8)=ders.row(2);//2,2
        break;
    default:
        break;
    }
}

/**
 * @brief Expression for the Hessian of a function
 * @ingroup Expressions
 * @tparam E The expression type
 */
template<class E>
class hess_expr : public _expr<hess_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;
public:
    enum {ScalarValued = 0, ColBlocks = (1==E::Space?1:0) };
    enum {Space = E::Space };

public:
    hess_expr(const E & u) : _u(u)
    {
        //gsInfo << "\n-expression is space ? "<<E::Space <<"\n"; _u.print(gsInfo);
        //GISMO_ASSERT(1==_u.dim(),"hess(.) requires 1D variable");
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const gsFuncData<Scalar> & dd = _u.data();
        const index_t sz = cardinality_impl();
        res.resize(dd.dim.first, sz*dd.dim.first);
        secDerToHessian(dd.values[2].col(k), dd.dim.first, res);
        res.resize(dd.dim.first, res.cols()*dd.dim.first);
        // Note: auto returns by value here,
        // in C++11 we may add in -> decltype(res) &
        return res;
    }

    index_t rows() const { return _u.data().dim.first; }
    index_t cols() const
    {   return rows();
        //return 2*_u.data().values[2].rows() / (1+_u.data().dim.first);
    }

    index_t cardinality_impl() const
    {
        return 2*_u.data().values[2].rows()/
            (_u.data().dim.first*(1+_u.data().dim.first));
        //gsDebugVar(_u.data().values.front().rows());//empty!
    }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_2ND_DER;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    void print(std::ostream &os) const
    //    { os << "hess("; _u.print(os);os <<")"; }
    { os << "\u210D(U)"; }
};

/**
 * @brief Expression for the Hessian of a FE solution
 * @ingroup Expressions
 * @tparam T The expression type
 */
template<class T>
class hess_expr<gsFeSolution<T> > : public _expr<hess_expr<gsFeSolution<T> > >
{
protected:
    const gsFeSolution<T> _u;

public:
    typedef T Scalar;
    enum{Space = 0, ScalarValued = 0, ColBlocks = 0 };

    hess_expr(const gsFeSolution<T> & u) : _u(u) { }

    mutable gsMatrix<T> deriv2, res;
    const gsMatrix<T> & eval(const index_t k) const
    {
        GISMO_ASSERT(_u.check(), "Invalid state in gsFeSolution");
        const gsDofMapper & map = _u.mapper();
        const index_t numActs = _u.data().values[0].rows();
        const index_t pdim = _u.parDim();
        index_t numDers = pdim*(pdim+1)/2;
        auto & act = _u.data().actives.col(1 == _u.data().actives.cols() ? 0:k );

        // In the scalar case, the hessian is returned as a pdim x pdim matrix
        if (1==_u.dim())
        {
            res.setZero(numDers,1);
            for (index_t i = 0; i!=numActs; ++i)
            {
                const index_t ii = map.index(act[i], _u.data().patchId, 0);
                deriv2 = _u.data().values[2].block(i*numDers,k,numDers,1);
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                    res += _u.coefs().at(ii) * deriv2;
                else
                    res +=_u.fixedPart().at( map.global_to_bindex(ii) ) * deriv2;
            }
            secDerToHessian(res, pdim, deriv2);
            res.swap(deriv2);
            res.resize(pdim,pdim);
        }
        // In the vector case, the hessian is returned as a matrix where each row corresponds to the component of the solution and contains the derivatives in the columns
        else
        {
            res.setZero(rows(), numDers);
            for (index_t c = 0; c != _u.dim(); c++)
                for (index_t i = 0; i != numActs; ++i)
                {
                    const index_t ii = map.index(act[i], _u.data().patchId, c);
                    deriv2 = _u.space().data().values[2].block(i * numDers, k, numDers,
                                                                1).transpose(); // start row, start col, rows, cols
                    if (map.is_free_index(ii)) // DoF value is in the solVector
                        res.row(c) += _u.coefs().at(ii) * deriv2;
                    else
                        res.row(c) += _u.fixedPart().at(map.global_to_bindex(ii)) * deriv2;
                }
        }
        return res;
    }

    index_t rows() const
    {
        if (1==_u.dim())
            return _u.parDim();
        else
            return _u.dim(); //  number of components
    }
    index_t cols() const
    {
        if (1==_u.dim())
            return _u.parDim();
        // second derivatives in the columns; i.e. [d11, d22, d33, d12, d13, d23]
        else
            return _u.parDim() * (_u.parDim() + 1) / 2;
    }

    const gsFeSpace<Scalar> & rowVar() const { return gsNullExpr<Scalar>::get(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList);                         // add symbol
        evList.add(_u.space());
        _u.data().flags |= NEED_ACTIVE | NEED_VALUE | NEED_DERIV2;
    }

    void print(std::ostream &os) const { os << "\u210D(s)"; }
};

/**
 * @brief Returns the Hessian of an expression
 * @ingroup Expressions
 * @param u The expression
 */
template<class E> EIGEN_STRONG_INLINE
hess_expr<E> hess(const symbol_expr<E> & u) { return hess_expr<E>(u); }

/**
 * @brief Returns the Hessian of a geometry map
 * @ingroup Expressions
 * @param u The geometry map
 */
template<class T> EIGEN_STRONG_INLINE
hess_expr<gsGeometryMap<T> > hess(const gsGeometryMap<T> & u) { return hess_expr<gsGeometryMap<T> >(u); }

/**
 * @brief Returns the Hessian of a \ref gsFeSolution
 * @ingroup Expressions
 * @param u The \ref gsFeSolution
 */
template<class T> EIGEN_STRONG_INLINE
hess_expr<gsFeSolution<T> > hess(const gsFeSolution<T> & u) { return hess_expr<gsFeSolution<T> >(u); }

}// namespace expr
}// namespace gismo