/** @file poisson2_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <typeinfo>
#include <gismo.h>
#include <gsAssembler/gsAdaptiveRefUtils.h>

#  define MatExprType  auto

namespace gismo{
namespace expr{

class unitVec_expr : public _expr<unitVec_expr >
{
public:
    typedef real_t Scalar;
private:
    index_t _dim;
    index_t _index;

public:
    unitVec_expr(const index_t index, const index_t dim) : _index(index), _dim(dim) { }

public:
    enum{ Space = 0 };

    gsMatrix<Scalar> eval(const index_t) const
    {
        gsMatrix<Scalar> vec = gsMatrix<Scalar>::Zero(_dim,1);
        // vec.setZero();
        vec(_index,0) = 1;
        return vec;
    }

    index_t rows() const { return _dim; }
    index_t cols() const { return  1; }
    void setFlag() const { }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & ) const {  }

    static constexpr bool rowSpan() {return false;}
    static bool colSpan() {return false;}

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "uv("<<_dim <<")";}
};

// Comments for var1:
// - TODO: dimensionm indep. later on
template<class E>
class var1_expr : public _expr<var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:

    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _G;

public:
    enum{ Space = E::Space };

    var1_expr(const E & u, const gsGeometryMap<Scalar> & G) : _u(u), _G(G) { }

    mutable gsMatrix<Scalar> res;

    mutable gsMatrix<Scalar> bGrads, cJac;
    mutable gsVector<Scalar,3> m_v, normal;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,k); }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD | NEED_ACTIVE;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    static constexpr bool rowSpan() {return E::rowSpan(); }
    static constexpr bool colSpan() {return false;}

    void print(std::ostream &os) const { os << "var1("; _u.print(os); os <<")"; }

private:
    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        const index_t A = _u.cardinality()/_u.targetDim();
        res.resize(_u.cardinality(), cols()); // rows()*

        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();
        bGrads = _u.data().values[1].col(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        // gsDebugVar(_G.data().values[0].col(k).transpose());

        for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all actives (only one component)
            {
                // Jac(u) ~ Jac(G) with alternating signs ?..
                m_v.noalias() = (vecFun(d, bGrads.at(2*j  ) ).cross( cJac.col(1).template head<3>() )
                              - vecFun(d, bGrads.at(2*j+1) ).cross( cJac.col(0).template head<3>() )) / measure;

                // ---------------  First variation of the normal
                // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
                res.row(s+j).noalias() = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
            }
        }
        return res;
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        res.resize(rows(), cols()); // rows()*
        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();
        grad_expr<U> vGrad = grad_expr(_u);

        bGrads = vGrad.eval(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        // gsDebugVar(_G.data().values[0].col(k).transpose());

        m_v.noalias() = ( ( bGrads.col(0).template head<3>() ).cross( cJac.col(1).template head<3>() )
                      -   ( bGrads.col(1).template head<3>() ).cross( cJac.col(0).template head<3>() ) ) / measure;

        // ---------------  First variation of the normal
        // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
        res = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
        return res;
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");
        solGrad_expr<Scalar> sGrad =  solGrad_expr(_u);
        res.resize(rows(), cols()); // rows()*

        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();
        bGrads = sGrad.eval(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        // gsDebugVar(_G.data().values[0].col(k).transpose());

        m_v.noalias() = ( ( bGrads.col(0).template head<3>() ).cross( cJac.col(1).template head<3>() )
                      -   ( bGrads.col(1).template head<3>() ).cross( cJac.col(0).template head<3>() ) ) / measure;

        // ---------------  First variation of the normal
        // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
        res = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
        return res;
    }
};

template<class E1, class E2>
class var1dif_expr : public _expr<var1dif_expr<E1, E2> >
{
public:
    typedef typename E1::Scalar Scalar;

private:

    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _G;

public:
    enum{ Space = E1::Space };

    var1dif_expr(const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & G) : _u(u), _v(v), _G(G) { }

    mutable gsMatrix<Scalar> res;

    mutable gsMatrix<Scalar> bGrads, cJac;
    mutable gsVector<Scalar,3> m_v, normal;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,_v,k); }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD | NEED_ACTIVE;
        _v.data().flags |= NEED_GRAD | NEED_ACTIVE;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        _v.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    static constexpr bool rowSpan() {return E1::rowSpan(); }
    static constexpr bool colSpan() {return false;}

    void print(std::ostream &os) const { os << "var1("; _u.print(os); os <<")"; }

private:
    template<class U, class V> inline
    typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value && util::is_same<V,gsFeVariable<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const V & v, const index_t k)  const
    {
        res.resize(rows(), cols()); // rows()*
        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();
        grad_expr<U> uGrad = grad_expr(_u);
        grad_expr<V> vGrad = grad_expr(_v);

        bGrads = uGrad.eval(k) - vGrad.eval(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        // gsDebugVar(_G.data().values[0].col(k).transpose());

        m_v.noalias() = ( ( bGrads.col(0).template head<3>() ).cross( cJac.col(1).template head<3>() )
                      -   ( bGrads.col(1).template head<3>() ).cross( cJac.col(0).template head<3>() ) ) / measure;

        // ---------------  First variation of the normal
        // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
        res = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
        return res;
    }

    template<class U, class V> inline
     typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value && util::is_same<V,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const V & v, const index_t k)  const
    {
        GISMO_ASSERT(1==_v.data().actives.cols(), "Single actives expected");
        grad_expr<U> uGrad = grad_expr(_u);
        solGrad_expr<Scalar> vGrad =  solGrad_expr(_v);
        res.resize(rows(), cols()); // rows()*

        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();
        bGrads = uGrad.eval(k) - vGrad.eval(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        // gsDebugVar(_G.data().values[0].col(k).transpose());

        m_v.noalias() = ( ( bGrads.col(0).template head<3>() ).cross( cJac.col(1).template head<3>() )
                      -   ( bGrads.col(1).template head<3>() ).cross( cJac.col(0).template head<3>() ) ) / measure;

        // ---------------  First variation of the normal
        // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
        res = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
        return res;
    }
};

// Comments for var2:
// - TODO: dimensionm indep. later on
// - TODO: how to structure this matrix
template<class E1, class E2, class E3>
class var2_expr : public _expr<var2_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:

    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _G;
    typename E3::Nested_t _Ef;

public:
    enum{ Space = E1::Space };

    var2_expr( const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & G, _expr<E3> const& Ef) : _u(u),_v(v), _G(G), _Ef(Ef) { }

    mutable gsMatrix<Scalar> res;

    mutable gsMatrix<Scalar> uGrads, vGrads, cJac, cDer2, evEf, result;
    mutable gsVector<Scalar> m_u, m_v, normal, m_uv, m_u_der, n_der, n_der2, tmp; // memomry leaks when gsVector<T,3>, i.e. fixed dimension
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        res.resize(_u.cardinality(), _u.cardinality());

        normal = _G.data().normal(k);
        normal.normalize();
        uGrads = _u.data().values[1].col(k);
        vGrads = _v.data().values[1].col(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        cDer2 = _G.data().values[2].reshapeCol(k, _G.data().dim.second, _G.data().dim.second);

        const index_t cardU = _u.data().values[0].rows(); // number of actives per component of u
        const index_t cardV = _v.data().values[0].rows(); // number of actives per component of v
        const Scalar measure =  _G.data().measures.at(k);

        evEf = _Ef.eval(k);

        for (index_t j = 0; j!= cardU; ++j) // for all basis functions u (1)
        {
            for (index_t i = 0; i!= cardV; ++i) // for all basis functions v (1)
            {
                for (index_t d = 0; d!= _u.dim(); ++d) // for all basis functions u (2)
                {
                    m_u.noalias() = ( vecFun(d, uGrads.at(2*j  ) ).cross( cJac.col(1).template head<3>() )
                                     -vecFun(d, uGrads.at(2*j+1) ).cross( cJac.col(0).template head<3>() ))
                                    / measure;

                    const short_t s = d*cardU;

                    for (index_t c = 0; c!= _v.dim(); ++c) // for all basis functions v (2)
                    {
                        const short_t r = c*cardV;
                        m_v.noalias() = ( vecFun(c, vGrads.at(2*i  ) ).cross( cJac.col(1).template head<3>() )
                                         -vecFun(c, vGrads.at(2*i+1) ).cross( cJac.col(0).template head<3>() ))
                                        / measure;

                        // n_der.noalias() = (m_v - ( normal.dot(m_v) ) * normal);
                        n_der.noalias() = (m_v - ( normal*m_v.transpose() ) * normal); // outer-product version

                        m_uv.noalias() = ( vecFun(d, uGrads.at(2*j  ) ).cross( vecFun(c, vGrads.at(2*i+1) ) )
                                          +vecFun(c, vGrads.at(2*i  ) ).cross( vecFun(d, uGrads.at(2*j+1) ) ))
                                          / measure; //check

                        m_u_der.noalias() = (m_uv - ( normal.dot(m_v) ) * m_u);
                        // m_u_der.noalias() = (m_uv - ( normal*m_v.transpose() ) * m_u); // outer-product version TODO

                        // ---------------  Second variation of the normal
                        tmp = m_u_der - (m_u.dot(n_der) + normal.dot(m_u_der) ) * normal - (normal.dot(m_u) ) * n_der;
                        // tmp = m_u_der - (m_u.dot(n_der) + normal.dot(m_u_der) ) * normal - (normal.dot(m_u) ) * n_der;

                        // Evaluate the product
                        tmp = cDer2 * tmp; // E_f_der2, last component
                        tmp.row(2) *= 2.0;
                        result = evEf * tmp;

                        res(s + j, r + i ) = result(0,0);
                    }
                }
            }
        }
        return res;
    }

    index_t rows() const
    {
        return 1; // because the resulting matrix has scalar entries for every combination of active basis functions
    }

    index_t cols() const
    {
        return 1; // because the resulting matrix has scalar entries for every combination of active basis functions
    }

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_2ND_DER | NEED_MEASURE;
        _Ef.setFlag();
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        _G.data().flags |=  NEED_NORMAL | NEED_DERIV | NEED_2ND_DER | NEED_MEASURE;
        _Ef.setFlag();
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    static constexpr bool rowSpan() {return true; }
    static constexpr bool colSpan() {return true; }

    void print(std::ostream &os) const { os << "var2("; _u.print(os); os <<")"; }
};


template<class E1, class E2>
class deriv2dot_expr : public _expr<deriv2dot_expr<E1, E2> >
{
public:
    typedef typename E1::Scalar Scalar;

private:

    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum{ Space = E1::Space };

    deriv2dot_expr(const E1 & u, const E2 & v) : _u(u), _v(v) { }

    mutable gsMatrix<Scalar> res,tmp, vEv;

    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,k); }

    index_t rows() const
    {
        return 1; //since the product with another vector is computed
    }

    index_t cols() const
    {
        return cols_impl(_u);
    }

    void setFlag() const
    {
        _u.data().flags |= NEED_DERIV2;
        _v.setFlag();
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_DERIV2;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    static constexpr bool rowSpan() {return E1::rowSpan(); }
    static constexpr bool colSpan() {return E2::rowSpan(); }

    void print(std::ostream &os) const { os << "deriv2("; _u.print(os); _v.print(os); os <<")"; }

private:
    template<class U> inline
    typename util::enable_if< util::is_same<U,gsGeometryMap<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        /*
            Here, we multiply the hessian of the geometry map by a vector, which possibly has multiple actives.
            The hessian of the geometry map c has the form: hess(c)
            [d11 c1, d11 c2, d11 c3]
            [d22 c1, d22 c2, d22 c3]
            [d12 c1, d12 c2, d12 c3]
            And we want to compute [d11 c .v; d22 c .v;  d12 c .v] ( . denotes a dot product and c and v are both vectors)
            So we simply evaluate for every active basis function v_k the product hess(c).v_k
        */


        // evaluate the geometry map of U
        tmp =_u.data().values[2].reshapeCol(k, cols(), _u.data().dim.second );
        vEv = _v.eval(k);
        res = vEv * tmp.transpose();
        return res;
    }

    template<class U> inline
    typename util::enable_if<util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k) const
    {
        /*
            We assume that the basis has the form v*e_i where e_i is the unit vector with 1 on index i and 0 elsewhere
            This implies that hess(v) = [hess(v_1), hess(v_2), hess(v_3)] only has nonzero entries in column i. Hence,
            hess(v) . normal = hess(v_i) * n_i (vector-scalar multiplication. The result is then of the form
            [hess(v_1)*n_1 .., hess(v_2)*n_2 .., hess(v_3)*n_3 ..]. Here, the dots .. represent the active basis functions.
        */
        const index_t numAct = u.data().values[0].rows();   // number of actives of a basis function
        const index_t cardinality = u.cardinality();        // total number of actives (=3*numAct)
        res.resize(rows()*cardinality, cols() );
        tmp.transpose() =_u.data().values[2].reshapeCol(k, cols(), numAct );
        vEv = _v.eval(k);

        for (index_t i = 0; i!=_u.dim(); i++)
            res.block(i*numAct, 0, numAct, cols() ) = tmp * vEv.at(i);

        return res;
    }

    template<class U> inline
    typename util::enable_if<util::is_same<U,gsFeVariable<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k) const
    {
        /*
            Here, we multiply the hessian of the geometry map by a vector, which possibly has multiple actives.
            The laplacian of the variable has the form: hess(v)
            [d11 c1, d22 c1, d12 c1]
            [d11 c2, d22 c2, d12 c2]
            [d11 c3, d22 c3, d12 c3]
            And we want to compute [d11 c .v; d22 c .v;  d12 c .v] ( . denotes a dot product and c and v are both vectors)
            So we simply evaluate for every active basis function v_k the product hess(c).v_k
        */

        gsMatrix<> tmp2;
        tmp =  u.data().values[2].col(k);
        index_t nDers = _u.source().domainDim() * (_u.source().domainDim() + 1) / 2;
        index_t dim = _u.source().targetDim();
        tmp2.resize(nDers,dim);
        for (index_t comp = 0; comp != u.source().targetDim(); comp++)
            tmp2.col(comp) = tmp.block(comp*nDers,0,nDers,1); //star,length

        vEv = _v.eval(k);
        res = vEv * tmp2.transpose();
        return res;
    }

    template<class U> inline
    typename util::enable_if<util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k) const
    {
        /*
            Here, we multiply the hessian of the geometry map by a vector, which possibly has multiple actives.
            The hessian of the geometry map c has the form: hess(c)
            [d11 c1, d22 c1, d12 c1]
            [d11 c2, d22 c2, d12 c2]
            [d11 c3, d22 c3, d12 c3]
            And we want to compute [d11 c .v; d22 c .v;  d12 c .v] ( . denotes a dot product and c and v are both vectors)
            So we simply evaluate for every active basis function v_k the product hess(c).v_k
        */

        solHess_expr<Scalar> sHess = solHess_expr(u);
        tmp = sHess.eval(k).transpose();
        vEv = _v.eval(k);
        res = vEv * tmp.transpose();
        return res;
    }


    template<class U> inline
    typename util::enable_if< util::is_same<U,gsGeometryMap<Scalar> >::value, index_t >::type
    cols_impl(const U & u)  const
    {
        return _u.data().dim.second;
    }

    template<class U> inline
    typename util::enable_if<util::is_same<U,gsFeSpace<Scalar> >::value || util::is_same<U,gsFeSolution<Scalar> >::value, index_t >::type
    cols_impl(const U & u) const
    {
        return _u.dim();
    }

    template<class U> inline
    typename util::enable_if<util::is_same<U,gsFeVariable<Scalar> >::value, index_t >::type
    cols_impl(const U & u) const
    {
        return _u.source().targetDim();
    }

};


/*
    The deriv2_expr computes the hessian of a basis.
    It assumes that the vector of basis functions is of the form v = u*e_i where u
    is the scalar basis function u: [0,1]^3 -> R^1 and e_i is the unit vector with a 1 on index i and a 0 elsewhere.
    Let us define the following blocks
    hess1(u) =              hess2(u) =              hess3(u) =
    [d11 u , 0 , 0 ]    |   [0 , d11 u , 0 ]     |  [0 , 0 , d11 u ]
    [d22 u , 0 , 0 ]    |   [0 , d22 u , 0 ]     |  [0 , 0 , d22 u ]
    [d12 u , 0 , 0 ]    |   [0 , d12 u , 0 ]     |  [0 , 0 , d12 u ]

    Then the deriv2(u) is defined as follows (for k number of actives)
    [hess1(u)_1]
    ...
    [hess1(u)_k]
    [hess2(u)_1]
    ...
    [hess2(u)_k]
    [hess3(u)_1]
    ...
    [hess3(u)_k]
**/
template<class E>
class deriv2_expr : public _expr<deriv2_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:

    typename E::Nested_t _u;

public:
    enum {ColBlocks = E::rowSpan() };
    enum{ Space = E::Space };

    deriv2_expr(const E & u) : _u(u) { }

    mutable gsMatrix<Scalar> res, tmp;

    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,k); }

    index_t rows() const //(components)
    {
        // return _u.data().values[2].rows() / _u.data().values[0].rows(); // numHessian dimensions
        // return _u.source().targetDim(); // no. dimensions should be 3
        return rows_impl(_u); // _u.dim() for space or targetDim() for geometry
    }

    index_t cols() const // number of function components (targetiDim)
    {
        return 3;
    }

    void setFlag() const
    {
        _u.data().flags |= NEED_DERIV2|NEED_ACTIVE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_DERIV2;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    static constexpr bool rowSpan() {return E::rowSpan();}
    static constexpr bool colSpan() {return E::colSpan();}

    void print(std::ostream &os) const { os << "deriv2("; _u.print(os); os <<")"; }

    private:
        template<class U> inline
        typename util::enable_if< util::is_same<U,gsGeometryMap<Scalar> >::value, const gsMatrix<Scalar> & >::type
        eval_impl(const U & u, const index_t k)  const
        {
            /*
                Here, we compute the hessian of the geometry map.
                The hessian of the geometry map c has the form: hess(c)
                [d11 c1, d11 c2, d11 c3]
                [d22 c1, d22 c2, d22 c3]
                [d12 c1, d12 c2, d12 c3]

                The geometry map has components c=[c1,c2,c3]
            */
            // evaluate the geometry map of U
            res = _u.data().values[2].reshapeCol(k, cols(), _u.data().dim.second );
            return res;
        }

        template<class U> inline
        typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value, const gsMatrix<Scalar> & >::type
        eval_impl(const U & u, const index_t k)  const
        {
            /*
                Here, we compute the hessian of the geometry map.
                The hessian of the geometry map c has the form: hess(c)
                [d11 c1, d11 c2, d11 c3]
                [d22 c1, d22 c2, d22 c3]
                [d12 c1, d12 c2, d12 c3]

                The geometry map has components c=[c1,c2,c3]
            */
            // evaluate the geometry map of U
            tmp =  _u.data().values[2];
            res.resize(rows(),cols());
            for (index_t comp = 0; comp != _u.source().targetDim(); comp++)
                res.col(comp) = tmp.block(comp*rows(),0,rows(),1); //star,length
            return res;
        }

        template<class U> inline
        typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
        eval_impl(const U & u, const index_t k)  const
        {
            /*
                Here, we compute the hessian of the geometry map.
                The hessian of the geometry map c has the form: hess(c)
                [d11 c1, d11 c2, d11 c3]
                [d22 c1, d22 c2, d22 c3]
                [d12 c1, d12 c2, d12 c3]

                The geometry map has components c=[c1,c2,c3]
            */
            // evaluate the geometry map of U
            solHess_expr<Scalar> sHess = solHess_expr(_u);
            res = sHess.eval(k).transpose();
            return res;
        }

        template<class U> inline
        typename util::enable_if<util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
        eval_impl(const U & u, const index_t k) const
        {
            /*
                Here, we compute the hessian of the basis with n actives.
                The hessian of the basis u has the form: hess(u)
                    active 1                active 2                        active n = cardinality
                [d11 u1, d11 u2, d11 u3] [d11 u1, d11 u2, d11 u3] ... [d11 u1, d11 u2, d11 u3]
                [d22 u1, d22 u2, d22 u3] [d22 u1, d22 u2, d22 u3] ... [d22 u1, d22 u2, d22 u3]
                [d12 u1, d12 u2, d12 u3] [d12 u1, d12 u2, d12 u3] ... [d12 u1, d12 u2, d12 u3]

                Here, the basis function has components u = [u1,u2,u3]. Since they are evaluated for scalars
                we use blockDiag to make copies for all components ui

                    active 1     active 2     active k = cardinality/dim   active 1           active 2k       active 1           active 2k
                [d11 u, 0, 0] [d11 u, 0, 0] ... [d11 u, 0, 0]            [0, d11 u, 0]  ... [0, d11 u, 0]  [0, d11 u, 0]  ... [0, d11 u, 0]
                [d22 u, 0, 0] [d22 u, 0, 0] ... [d22 u, 0, 0]            [0, d22 u, 0]  ... [0, d22 u, 0]  [0, d22 u, 0]  ... [0, d22 u, 0]
                [d12 u, 0, 0] [d12 u, 0, 0] ... [d12 u, 0, 0]            [0, d12 u, 0]  ... [0, d12 u, 0]  [0, d12 u, 0]  ... [0, d12 u, 0]

            */
            const index_t numAct = u.data().values[0].rows();   // number of actives of a basis function
            const index_t cardinality = u.cardinality();        // total number of actives (=3*numAct)

            res.resize(rows(), _u.dim() *_u.cardinality()); // (3 x 3*cardinality)
            res.setZero();

            tmp = _u.data().values[2].reshapeCol(k, cols(), numAct );
            for (index_t d = 0; d != cols(); ++d)
            {
                const index_t s = d*(cardinality + 1);
                for (index_t i = 0; i != numAct; ++i)
                    res.col(s+i*_u.cols()) = tmp.col(i);
            }

            // res = _u.data().values[2].col(k).transpose().blockDiag(_u.targetDim()); // analoguous to jacobian..
            // res = res.transpose();
            // gsDebugVar(res);
            return res;
        }

        template<class U> inline
        typename util::enable_if< util::is_same<U,gsFeVariable<Scalar> >::value || util::is_same<U,gsGeometryMap<Scalar> >::value || util::is_same<U,gsFeSpace<Scalar> >::value, index_t >::type
        rows_impl(const U & u)  const
        {
            return _u.source().domainDim() * ( _u.source().domainDim() + 1 ) / 2;
        }

        template<class U> inline
        typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, index_t >::type
        rows_impl(const U & u) const
        {
            return _u.dim();
        }

};



// template<class E1, class E2>
// class hessdot_expr : public _expr<hessdot_expr<E1,E2> >
// {
//     typename E1::Nested_t _u;
//     typename E2::Nested_t _v;

// public:
//     enum{ Space = E1::Space };

//     typedef typename E1::Scalar Scalar;

//     hessdot_expr(const E1 & u, const E2 & v) : _u(u), _v(v) {}

//     mutable gsMatrix<Scalar> res, hess, tmp;
//     mutable gsMatrix<Scalar> normalMat;

//     MatExprType eval(const index_t k) const
//     {
//         const gsFuncData<Scalar> & udata = _u.data(); // udata.values[2].col(k)
//         const index_t numAct = udata.values[0].rows();
//         const gsAsConstMatrix<Scalar> ders(udata.values[2].col(k).data(), 3, numAct );

//         tmp = _v.eval(k);

//         res.resize(rows(), cols() );


//             for (index_t i = 0; i!=tmp.rows(); ++i)
//             {
//                 res.block(i*numAct, 0, numAct, 3).noalias() = ders.transpose() * tmp.at(i);
//             }

//         return res;
//     }

//     index_t rows() const
//     {
//         return _u.dim() * _u.data().values[0].rows();
//     }

//     index_t cols() const
//     {
//         return // ( E2::rowSpan() ? rows() : 1 ) *
//             _u.data().values[2].rows() / _u.data().values[0].rows();//=3
//     }

//     void setFlag() const
//     {
//         _u.data().flags |= NEED_2ND_DER;
//         _v.setFlag();
//     }

//     void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
//     {
//         //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
//         evList.push_sorted_unique(&_u.source());
//         _u.data().flags |= NEED_2ND_DER;
//     }

//     const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
//     const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

//     static constexpr bool rowSpan() {return E1::rowSpan(); }
//     static constexpr bool colSpan() {return E2::rowSpan(); }

//     void print(std::ostream &os) const { os << "hessdot("; _u.print(os); os <<")"; }
// };


/**
   TO ADD
 */
template<class E1, class E2, class E3>
class flatdot_expr  : public _expr<flatdot_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0, Space = E1::Space};
private:
    typename E1::Nested_t _A;
    typename E2::Nested_t _B;
    typename E3::Nested_t _C;
    mutable gsMatrix<Scalar> eA, eB, eC, tmp, res;

public:

    flatdot_expr(_expr<E1> const& A, _expr<E2> const& B, _expr<E3> const& C) : _A(A),_B(B),_C(C)
    {
        //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t Ac = _A.cols();
        const index_t An = _A.cardinality();
        const index_t Bc = _B.cols();
        const index_t Bn = _B.cardinality();

        eA = _A.eval(k);
        eB = _B.eval(k);
        eC = _C.eval(k);

        // gsDebugVar(eA);
        // gsDebugVar(eB);
        // gsDebugVar(eC);

        // gsDebugVar(eA.rows());
        // gsDebugVar(eA.cols());
        // gsDebugVar(_A.rows());
        // gsDebugVar(_A.cols());
        // gsDebugVar(eB.rows());
        // gsDebugVar(eB.cols());
        // gsDebugVar(_B.rows());
        // gsDebugVar(_B.cols());

        GISMO_ASSERT(Bc==_A.rows(), "Dimensions: "<<Bc<<","<< _A.rows()<< "do not match");
        GISMO_ASSERT(_A.rowSpan(), "First entry should be rowSpan");
        GISMO_ASSERT(_B.colSpan(), "Second entry should be colSpan.");

        res.resize(An, Bn);
        for (index_t i = 0; i!=An; ++i)
            for (index_t j = 0; j!=Bn; ++j)
            {
                tmp.noalias() = eB.middleCols(i*Bc,Bc) * eA.middleCols(j*Ac,Ac);
                tmp(0,0) *= eC.at(0);
                tmp(0,1) *= eC.at(2);
                tmp(1,0) *= eC.at(2);
                tmp(1,1) *= eC.at(1);
                res(i,j) = tmp.sum();
            }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 1; }
    void setFlag() const { _A.setFlag();_B.setFlag();_C.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _A.parse(evList);_B.parse(evList);_C.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _A.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _B.colVar(); }
    index_t cardinality_impl() const { return _A.cardinality_impl(); }

    static constexpr bool rowSpan() {return E1::rowSpan();}
    static constexpr bool colSpan() {return E2::colSpan();}

    void print(std::ostream &os) const { os << "flatdot("; _A.print(os);_B.print(os);_C.print(os); os<<")"; }
};

/**
   To Do:
   *    Improve by inputting u instead of deriv2(u)
 */
template<class E1, class E2, class E3>
class flatdot2_expr  : public _expr<flatdot2_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0, Space = E1::Space};
    enum {ColBlocks = 0 };

private:
    typename E1::Nested_t _A;
    typename E2::Nested_t _B;
    typename E3::Nested_t _C;
    mutable gsMatrix<Scalar> eA, eB, eC, res, tmp;

public:

    flatdot2_expr(_expr<E1> const& A, _expr<E2> const& B, _expr<E3> const& C) : _A(A),_B(B),_C(C)
    {
        //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t Ac = _A.cols();
        const index_t An = _A.cardinality();
        const index_t Bn = _B.cardinality();

        eA = _A.eval(k);
        eB = _B.eval(k);
        eC = _C.eval(k);

        // gsDebugVar(eA);
        // gsDebugVar(eB);
        // gsDebugVar(eC);

        // gsDebugVar(eA.rows());
        // gsDebugVar(eA.cols());
        // gsDebugVar(_A.rows());
        // gsDebugVar(_A.cols());
        // gsDebugVar(eB.rows());
        // gsDebugVar(eB.cols());
        // gsDebugVar(_B.rows());
        // gsDebugVar(_B.cols());

        GISMO_ASSERT(_B.rows()==_A.cols(), "Dimensions: "<<_B.rows()<<","<< _A.cols()<< "do not match");
        GISMO_ASSERT(_A.rowSpan(), "First entry should be rowSpan");
        GISMO_ASSERT(_B.colSpan(), "Second entry should be colSpan.");
        GISMO_ASSERT(_C.cols()==_B.rows(), "Dimensions: "<<_C.rows()<<","<< _B.rows()<< "do not match");

        // NOTE: product moved to the loop since that is more consistent with the formulations
        // for (index_t i = 0; i!=An; ++i)
        //     for (index_t j = 0; j!=Ac; ++j)
        //         eA.middleCols(i*Ac,Ac).row(j) *= eC(j);

        res.resize(An, Bn);
        for (index_t i = 0; i!=An; ++i)
            for (index_t j = 0; j!=Bn; ++j)
            {
                tmp = eA.middleCols(i*Ac,Ac) * eB.col(j);   // E_f_der2
                tmp.row(2) *= 2.0;                          // multiply the third row of E_f_der2 by 2 for voight notation
                res(i,j) = eC.row(0) * tmp.col(0);          // E_f^T * mm * E_f_der2
            }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 1; }
    void setFlag() const { _A.setFlag();_B.setFlag();_C.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _A.parse(evList);_B.parse(evList);_C.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _A.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _B.colVar(); }
    index_t cardinality_impl() const { return _A.cardinality_impl(); }

    static constexpr bool rowSpan() {return E1::rowSpan();}
    static constexpr bool colSpan() {return E2::colSpan();}

    void print(std::ostream &os) const { os << "flatdot2("; _A.print(os);_B.print(os);_C.print(os); os<<")"; }
};

/*
   Expression for the transformation matrix between local cartesian and covariant bases, based on a geometry map
 */
template<class T> class cartcovinv_expr ;

template<class T>
class cartcov_expr : public _expr<cartcov_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    cartcov_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsMatrix<Scalar,3,3> covBasis, conBasis, covMetric, conMetric, cartBasis, result;
    mutable gsVector<Scalar,3> normal, tmp;
    mutable gsVector<Scalar,3> e1, e2, a1, a2;

    MatExprType eval(const index_t k) const
    {
        // Compute covariant bases in deformed and undeformed configuration
        normal = _G.data().normals.col(k);
        normal.normalize();
        covBasis.leftCols(2) = _G.data().jacobian(k);
        covBasis.col(2)      = normal;
        covMetric = covBasis.transpose() * covBasis;

        conMetric = covMetric.inverse();

        // conBasis.col(0) = conMetric(0,0)*covBasis.col(0)+conMetric(0,1)*covBasis.col(1)+conMetric(0,2)*covBasis.col(2);
        conBasis.col(1) = conMetric(1,0)*covBasis.col(0)+conMetric(1,1)*covBasis.col(1)+conMetric(1,2)*covBasis.col(2);
        // conBasis.col(2) = conMetric(2,0)*covBasis.col(0)+conMetric(2,1)*covBasis.col(1)+conMetric(2,2)*covBasis.col(2);

        e1 = covBasis.col(0); e1.normalize();
        e2 = conBasis.col(1); e2.normalize();
        // e3 = normal;

        a1 = covBasis.col(0);
        a2 = covBasis.col(1);

        result(0,0) = (e1.dot(a1))*(a1.dot(e1));
        result(0,1) = (e1.dot(a2))*(a2.dot(e2));
        result(0,2) = 2*(e1.dot(a1))*(a2.dot(e1));
        // Row 1
        result(1,0) = (e2.dot(a1))*(a1.dot(e2));
        result(1,1) = (e2.dot(a2))*(a2.dot(e2));
        result(1,2) = 2*(e2.dot(a1))*(a2.dot(e2));
        // Row 2
        result(2,0) = (e1.dot(a1))*(a1.dot(e2));
        result(2,1) = (e1.dot(a2))*(a2.dot(e2));
        result(2,2) = (e1.dot(a1))*(a2.dot(e2)) + (e1.dot(a2))*(a1.dot(e2));

        // return result.inverse(); // !!!!
        return result;
    }

    cartcovinv_expr<T> inv() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return cartcovinv_expr<T>(_G);
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    static constexpr bool rowSpan() {return false; }
    static bool colSpan() {return false;}

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void setFlag() const { _G.data().flags |= NEED_NORMAL|NEED_DERIV; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_NORMAL|NEED_DERIV;
    }

    void print(std::ostream &os) const { os << "cartcov("; _G.print(os); os <<")"; }
};

template<class T>
class cartcovinv_expr : public _expr<cartcovinv_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    cartcovinv_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsMatrix<T> temp;

    MatExprType eval(const index_t k) const
    {
        return  (cartcov_expr(_G).eval(k)).reshape(3,3).inverse();
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    static constexpr bool rowSpan() {return false; }
    static bool colSpan() {return false;}

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void setFlag() const { _G.data().flags |= NEED_NORMAL|NEED_DERIV; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_NORMAL|NEED_DERIV;
    }

    void print(std::ostream &os) const { os << "cartcovinv("; _G.print(os); os <<")"; }
};


/*
   Expression for the transformation matrix between local cartesian and contravariant bases, based on a geometry map
 */
template<class T> class cartconinv_expr ;

template<class T>
class cartcon_expr : public _expr<cartcon_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    cartcon_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsMatrix<Scalar,3,3> covBasis, conBasis, covMetric, conMetric, cartBasis, result;
    mutable gsVector<Scalar,3> normal, tmp;
    mutable gsVector<Scalar,3> e1, e2, ac1, ac2;

    MatExprType eval(const index_t k) const
    {
        // Compute covariant bases in deformed and undeformed configuration
        normal = _G.data().normals.col(k);
        normal.normalize();
        covBasis.leftCols(2) = _G.data().jacobian(k);
        covBasis.col(2)      = normal;
        covMetric = covBasis.transpose() * covBasis;

        conMetric = covMetric.inverse();

        conBasis.col(0) = conMetric(0,0)*covBasis.col(0)+conMetric(0,1)*covBasis.col(1)+conMetric(0,2)*covBasis.col(2);
        conBasis.col(1) = conMetric(1,0)*covBasis.col(0)+conMetric(1,1)*covBasis.col(1)+conMetric(1,2)*covBasis.col(2);
        // conBasis.col(2) = conMetric(2,0)*covBasis.col(0)+conMetric(2,1)*covBasis.col(1)+conMetric(2,2)*covBasis.col(2);

        e1 = covBasis.col(0); e1.normalize();
        e2 = conBasis.col(1); e2.normalize();
        // e3 = normal;

        ac1 = conBasis.col(0);
        ac2 = conBasis.col(1);

        result(0,0) = (e1.dot(ac1))*(ac1.dot(e1));
        result(0,1) = (e1.dot(ac2))*(ac2.dot(e2));
        result(0,2) = 2*(e1.dot(ac1))*(ac2.dot(e1));
        // Row 1
        result(1,0) = (e2.dot(ac1))*(ac1.dot(e2));
        result(1,1) = (e2.dot(ac2))*(ac2.dot(e2));
        result(1,2) = 2*(e2.dot(ac1))*(ac2.dot(e2));
        // Row 2
        result(2,0) = (e1.dot(ac1))*(ac1.dot(e2));
        result(2,1) = (e1.dot(ac2))*(ac2.dot(e2));
        result(2,2) = (e1.dot(ac1))*(ac2.dot(e2)) + (e1.dot(ac2))*(ac1.dot(e2));

        return result;
    }

    cartconinv_expr<T> inv() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return cartconinv_expr<T>(_G);
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    static constexpr bool rowSpan() {return false; }
    static bool colSpan() {return false;}

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void setFlag() const { _G.data().flags |= NEED_NORMAL|NEED_DERIV; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_NORMAL|NEED_DERIV;
    }

    void print(std::ostream &os) const { os << "cartcon("; _G.print(os); os <<")"; }
};

template<class T>
class cartconinv_expr : public _expr<cartconinv_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    cartconinv_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsMatrix<T> temp;

    MatExprType eval(const index_t k) const
    {
        return  (cartcon_expr(_G).eval(k)).reshape(3,3).inverse();
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    static constexpr bool rowSpan() {return false; }
    static bool colSpan() {return false;}

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void setFlag() const { _G.data().flags |= NEED_NORMAL|NEED_DERIV; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_NORMAL|NEED_DERIV;
    }

    void print(std::ostream &os) const { os << "cartconinv("; _G.print(os); os <<")"; }
};

EIGEN_STRONG_INLINE
unitVec_expr uv(const index_t index, const index_t dim) { return unitVec_expr(index,dim); }

template<class E> EIGEN_STRONG_INLINE
var1_expr<E> var1(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return var1_expr<E>(u, G); }
template<class E1, class E2> EIGEN_STRONG_INLINE
var1dif_expr<E1,E2> var1dif(const E1 & u,const E2 & v, const gsGeometryMap<typename E1::Scalar> & G) { return var1dif_expr<E1,E2>(u,v, G); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
var2_expr<E1,E2,E3> var2(const E1 & u, const E2 & v, const gsGeometryMap<typename E1::Scalar> & G, const E3 & Ef)
{ return var2_expr<E1,E2,E3>(u,v, G, Ef); }

// template<class E1, class E2> EIGEN_STRONG_INLINE
// hessdot_expr<E1,E2> hessdot(const E1 & u, const E2 & v) { return hessdot_expr<E1,E2>(u, v); }

template<class E> EIGEN_STRONG_INLINE
deriv2_expr<E> deriv2(const E & u) { return deriv2_expr<E>(u); }

template<class E1, class E2> EIGEN_STRONG_INLINE
deriv2dot_expr<E1, E2> deriv2(const E1 & u, const E2 & v) { return deriv2dot_expr<E1, E2>(u,v); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
flatdot_expr<E1,E2,E3> flatdot(const E1 & u, const E2 & v, const E3 & w)
{ return flatdot_expr<E1,E2,E3>(u, v, w); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
flatdot2_expr<E1,E2,E3> flatdot2(const E1 & u, const E2 & v, const E3 & w)
{ return flatdot2_expr<E1,E2,E3>(u, v, w); }

template<class E> EIGEN_STRONG_INLINE
cartcov_expr<E> cartcov(const gsGeometryMap<E> & G) { return cartcov_expr<E>(G); }

template<class E> EIGEN_STRONG_INLINE
cartcon_expr<E> cartcon(const gsGeometryMap<E> & G) { return cartcon_expr<E>(G); }

}
}


using namespace gismo;

template <class T>
void evaluateFunction(gsExprEvaluator<T> ev, auto expression, gsVector<T> pt);
template <class T>
void evaluateFunction(gsExprEvaluator<T> ev, auto expression, gsMatrix<T> pt);
template <class T>
void constructSolution(gsExprAssembler<T> assembler, gsMultiPatch<T> mp_def);

// To Do:
// * struct for material model
// Input is parametric coordinates of the surface \a mp
template <class T>
class gsMaterialMatrix : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
  //
protected:
    const gsFunctionSet<T> * _mp;
    const gsFunction<T> * _YoungsModulus;
    const gsFunction<T> * _PoissonRatio;
    mutable gsMapData<T> _tmp;
    mutable gsMatrix<real_t,3,3> F0;
    mutable gsMatrix<T> Emat,Nmat;
    mutable real_t lambda, mu, E, nu, C_constant;

public:
    /// Shared pointer for gsMaterialMatrix
    typedef memory::shared_ptr< gsMaterialMatrix > Ptr;

    /// Unique pointer for gsMaterialMatrix
    typedef memory::unique_ptr< gsMaterialMatrix > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsMaterialMatrix(const gsFunctionSet<T> & mp, const gsFunction<T> & YoungsModulus,
                   const gsFunction<T> & PoissonRatio) :
    _mp(&mp), _YoungsModulus(&YoungsModulus), _PoissonRatio(&PoissonRatio), _mm_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
    }

    ~gsMaterialMatrix() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrix)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 9;}

    mutable gsMaterialMatrix<T> * _mm_piece; // todo: improve the way pieces are accessed

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _mm_piece;
        _mm_piece = new gsMaterialMatrix(_mp->piece(k), *_YoungsModulus, *_PoissonRatio);
        return *_mm_piece;
    }

    //class .. matMatrix_z
    // should contain eval_into(thickness variable)

    // Input is parametric coordinates of the surface \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // NOTE 1: if the input \a u is considered to be in physical coordinates
        // then we first need to invert the points to parameter space
        // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
        // otherwise we just use the input paramteric points
        _tmp.points = u;

        static_cast<const gsFunction<T>&>(_mp->piece(0)).computeMap(_tmp); // the piece(0) here implies that if you call class.eval_into, it will be evaluated on piece(0). Hence, call class.piece(k).eval_into()

        // NOTE 2: in the case that parametric value is needed it suffices
        // to evaluate Youngs modulus and Poisson's ratio at
        // \a u instead of _tmp.values[0].
        _YoungsModulus->eval_into(_tmp.values[0], Emat);
        _PoissonRatio->eval_into(_tmp.values[0], Nmat);

        result.resize( targetDim() , u.cols() );
        for( index_t i=0; i< u.cols(); ++i )
        {
            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(i,3,3);

            F0.leftCols(2) = _tmp.jacobian(i);
            F0.col(2)      = _tmp.normal(i).normalized();
            F0 = F0.inverse();
            F0 = F0 * F0.transpose(); //3x3

            // Evaluate material properties on the quadrature point
            E = Emat(0,i);
            nu = Nmat(0,i);
            lambda = E * nu / ( (1. + nu)*(1.-2.*nu)) ;
            mu     = E / (2.*(1. + nu)) ;

            C_constant = 2*lambda*mu/(lambda+2*mu);

            C(0,0) = C_constant*F0(0,0)*F0(0,0) + 1*mu*(2*F0(0,0)*F0(0,0));
            C(1,1) = C_constant*F0(1,1)*F0(1,1) + 1*mu*(2*F0(1,1)*F0(1,1));
            C(2,2) = C_constant*F0(0,1)*F0(0,1) + 1*mu*(F0(0,0)*F0(1,1) + F0(0,1)*F0(0,1));
            C(1,0) =
            C(0,1) = C_constant*F0(0,0)*F0(1,1) + 1*mu*(2*F0(0,1)*F0(0,1));
            C(2,0) =
            C(0,2) = C_constant*F0(0,0)*F0(0,1) + 1*mu*(2*F0(0,0)*F0(0,1));
            C(2,1) = C(1,2) = C_constant*F0(0,1)*F0(1,1) + 1*mu*(2*F0(0,1)*F0(1,1));

            //gsDebugVar(C);
        }
    }

    // piece(k) --> for patch k

}; //! [Include namespace]

template <class T>
class gsMaterialMatrixD : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
  //
protected:
    const gsFunctionSet<T> * _mp;
    const gsFunction<T> * _YoungsModulus;
    const gsFunction<T> * _PoissonRatio;
    mutable gsMapData<T> _tmp;
    mutable gsMatrix<T> Emat,Nmat;
    mutable real_t E, nu, D_constant;

public:
    /// Shared pointer for gsMaterialMatrixD
    typedef memory::shared_ptr< gsMaterialMatrixD > Ptr;

    /// Unique pointer for gsMaterialMatrixD
    typedef memory::unique_ptr< gsMaterialMatrixD > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsMaterialMatrixD(const gsFunctionSet<T> & mp, const gsFunction<T> & YoungsModulus,
                   const gsFunction<T> & PoissonRatio) :
    _mp(&mp), _YoungsModulus(&YoungsModulus), _PoissonRatio(&PoissonRatio), _mm_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
    }

    ~gsMaterialMatrixD() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrixD)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 9;}

    mutable gsMaterialMatrixD<T> * _mm_piece; // todo: improve the way pieces are accessed

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _mm_piece;
        _mm_piece = new gsMaterialMatrixD(_mp->piece(k), *_YoungsModulus, *_PoissonRatio);
        return *_mm_piece;
    }

    //class .. matMatrix_z
    // should contain eval_into(thickness variable)

    // Input is parametric coordinates of the surface \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // NOTE 1: if the input \a u is considered to be in physical coordinates
        // then we first need to invert the points to parameter space
        // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
        // otherwise we just use the input paramteric points
        _tmp.points = u;

        static_cast<const gsFunction<T>*>(_mp)->computeMap(_tmp);

        // NOTE 2: in the case that parametric value is needed it suffices
        // to evaluate Youngs modulus and Poisson's ratio at
        // \a u instead of _tmp.values[0].
        _YoungsModulus->eval_into(_tmp.values[0], Emat);
        _PoissonRatio->eval_into(_tmp.values[0], Nmat);

        result.resize( targetDim() , u.cols() );
        for( index_t i=0; i< u.cols(); ++i )
        {
            gsAsMatrix<T, Dynamic, Dynamic> D = result.reshapeCol(i,3,3);

            // Evaluate material properties on the quadrature point
            E = Emat(0,i);
            nu = Nmat(0,i);

            D_constant = E/(1-nu*nu);
            D(0,0) = D(1,1) = 1;
            D(2,2) = (1-nu)/2;
            D(0,1) = D(1,0) = nu;
            D(2,0) = D(0,2) = D(2,1) = D(1,2) = 0;
            D *= D_constant;

            //gsDebugVar(C);
        }
    }

    // piece(k) --> for patch k

}; //! [Include namespace]

/*
    Todo:
        * Improve for mu, E, phi as gsFunction instead of reals!!
*/
template <class T, int mat>
class gsMaterialMatrixComp : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
  // NOTE: This material matrix is in local Cartesian coordinates and should be transformed!
  // NOTE: To make this efficient, we can output all matrices stacked and use expressions to pick the 1st, second and third
protected:
    // const gsFunctionSet<T> * _mp;
    const std::vector<std::pair<T,T>> _YoungsModuli;
    const std::vector<T> _ShearModuli;
    const std::vector<std::pair<T,T>> _PoissonRatios;
    const std::vector<T> _thickness;
    const std::vector<T> _phi;
    mutable gsMapData<T> _tmp;
    mutable real_t E1, E2, G12, nu12, nu21, t, t_tot, t_temp, z, z_mid, phi;
    mutable gsMatrix<T> Tmat, Dmat;

public:
    /// Shared pointer for gsMaterialMatrixComp
    typedef memory::shared_ptr< gsMaterialMatrixComp > Ptr;

    /// Unique pointer for gsMaterialMatrixComp
    typedef memory::unique_ptr< gsMaterialMatrixComp > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsMaterialMatrixComp(  //const gsFunctionSet<T> & mp,
                        const std::vector<std::pair<T,T>> & YoungsModuli,
                        const std::vector<T> & ShearModuli,
                        const std::vector<std::pair<T,T>> & PoissonRatios,
                        const std::vector<T> thickness,
                        const std::vector<T> phi) :
    // _mp(&mp),
    _YoungsModuli(YoungsModuli),
    _ShearModuli(ShearModuli),
    _PoissonRatios(PoissonRatios),
    _thickness(thickness),
    _phi(phi)//,
    // _mm_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
    }

    // ~gsMaterialMatrixComp() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrixComp)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 9;}

    // mutable gsMaterialMatrixComp<T> * _mm_piece; // todo: improve the way pieces are accessed

    // const gsFunction<T> & piece(const index_t k) const
    // {
    //     delete _mm_piece;
    //     _mm_piece = new gsMaterialMatrixComp(_mp->piece(k), _YoungsModuli, _ShearModuli, _PoissonRatios, _thickness, _phi);
    //     return *_mm_piece;
    // }

    // Input is parametric coordinates of the surface \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // static_cast<const gsFunction<T>*>(_mp)->computeMap(_tmp);
        GISMO_ASSERT(_YoungsModuli.size()==_PoissonRatios.size(),"Size of vectors of Youngs Moduli and Poisson Ratios is not equal: " << _YoungsModuli.size()<<" & "<<_PoissonRatios.size());
        GISMO_ASSERT(_YoungsModuli.size()==_ShearModuli.size(),"Size of vectors of Youngs Moduli and Shear Moduli is not equal: " << _YoungsModuli.size()<<" & "<<_ShearModuli.size());
        GISMO_ASSERT(_thickness.size()==_phi.size(),"Size of vectors of thickness and angles is not equal: " << _thickness.size()<<" & "<<_phi.size());
        GISMO_ASSERT(_YoungsModuli.size()==_thickness.size(),"Size of vectors of material properties and laminate properties is not equal: " << _YoungsModuli.size()<<" & "<<_thickness.size());
        GISMO_ASSERT(_YoungsModuli.size()!=0,"No laminates defined");

        // Compute total thickness (sum of entries)
        t_tot = std::accumulate(_thickness.begin(), _thickness.end(), 0.0);

        // compute mid-plane height of total plate
        z_mid = t_tot / 2.0;

        // now we use t_temp to add the thickness of all plies iteratively
        t_temp = 0.0;

        // Initialize material matrix and result
        Dmat.resize(3,3);
        result.resize( targetDim(), 1 );

        // Initialize transformation matrix
        Tmat.resize(3,3);


        for (size_t i = 0; i != _phi.size(); ++i) // loop over laminates
        {
            // Compute all quantities
            E1 = _YoungsModuli[i].first;
            E2 = _YoungsModuli[i].second;
            G12 = _ShearModuli[i];
            nu12 = _PoissonRatios[i].first;
            nu21 = _PoissonRatios[i].second;
            t = _thickness[i];
            phi = _phi[i];

            GISMO_ASSERT(nu21*E1 == nu12*E2, "No symmetry in material properties for ply "<<i<<". nu12*E2!=nu21*E1:\n"<<
                    "\tnu12 = "<<nu12<<"\t E2 = "<<E2<<"\t nu12*E2 = "<<nu12*E2<<"\n"
                  <<"\tnu21 = "<<nu21<<"\t E1 = "<<E1<<"\t nu21*E1 = "<<nu21*E1);

            // Fill material matrix
            Dmat(0,0) = E1 / (1-nu12*nu21);
            Dmat(1,1) = E2 / (1-nu12*nu21);;
            Dmat(2,2) = G12;
            Dmat(0,1) = nu21*E1 / (1-nu12*nu21);
            Dmat(1,0) = nu12*E2 / (1-nu12*nu21);
            Dmat(2,0) = Dmat(0,2) = Dmat(2,1) = Dmat(1,2) = 0.0;

            // Make transformation matrix
            Tmat(0,0) = Tmat(1,1) = math::pow(math::cos(phi),2);
            Tmat(0,1) = Tmat(1,0) = math::pow(math::sin(phi),2);
            Tmat(2,0) = Tmat(0,2) = Tmat(2,1) = Tmat(1,2) = math::sin(phi) * math::cos(phi);
            Tmat(2,0) *= -2.0;
            Tmat(2,1) *= 2.0;
            Tmat(1,2) *= -1.0;
            Tmat(2,2) = math::pow(math::cos(phi),2) - math::pow(math::sin(phi),2);

            // Compute laminate stiffness matrix
            Dmat = Tmat.transpose() * Dmat * Tmat;

            z = math::abs(z_mid - (t/2.0 + t_temp) ); // distance from mid-plane of plate

            // Make matrices A, B and C
            // [NOTE: HOW TO DO THIS NICELY??]
            // result.reshape(3,3) += Dmat * t; // A
            // result.reshape(3,3) += Dmat * t*z; // B
            // result.reshape(3,3) += Dmat * ( t*z*z + t*t*t/12.0 ); // D
            switch (mat)
            {
                case 0:
                    result.reshape(3,3) += Dmat * t; // A
                    break;
                case 1:
                    result.reshape(3,3) += Dmat * t; // A
                    break;
                case 2:
                    result.reshape(3,3) += Dmat * t; // A
                    break;
            }

            t_temp += t;
        }

        GISMO_ASSERT(t_tot==t_temp,"Total thickness after loop is wrong. t_temp = "<<t_temp<<" and sum(thickness) = "<<t_tot);

        // Replicate for all points since the quantities are equal over the whole domain
        result.replicate(1, u.cols());
    }

    // piece(k) --> for patch k

};


template<typename T>
class gsElementErrorPlotter : public gsFunction<T>
{
public:
    gsElementErrorPlotter(const gsBasis<T>& mp, const std::vector<T>& errors ) : m_mp(mp),m_errors(errors)
    {

    }

    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& res) const
    {
        // Initialize domain element iterator -- using unknown 0
        res.setZero(1,u.cols());
        for(index_t i=0; i<u.cols();++i)
        {
            int iter =0;
            // Start iteration over elements

            typename gsBasis<T>::domainIter domIt = m_mp.makeDomainIterator();
            for (; domIt->good(); domIt->next() )
            {
                 bool flag = true;
                const gsVector<T>& low = domIt->lowerCorner();
                const gsVector<T>& upp = domIt->upperCorner();


                for(int d=0; d<domainDim();++d )
                {
                    if(low(d)> u(d,i) || u(d,i) > upp(d))
                    {
                        flag = false;
                        break;
                    }
                }
                if(flag)
                {
                     res(0,i) = m_errors.at(iter);
                     break;
                }
                iter++;
            }
        }
    }

    short_t domainDim() const { return m_mp.dim();}

private:
    const gsBasis<T>& m_mp;
    const std::vector<T>& m_errors;
};

int main(int argc, char *argv[])
{
    // Number of adaptive refinement loops
    index_t RefineLoopMax;
    // Flag for refinemet criterion
    // (see doxygen documentation of the free function
    // gsMarkElementsForRef explanation)
    index_t refCriterion;
    // Parameter for computing adaptive refinement threshold
    // (see doxygen documentation of the free function
    // gsMarkElementsForRef explanation)
    real_t refParameter;  // ...specified below with the examples


    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t goal = 1;
    bool nonlinear = false;
    std::string fn;

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t thickness = 1.0;

    refCriterion = GARU;
    refParameter = 0.85;
    RefineLoopMax = 1;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt("R", "refine", "Maximum number of adaptive refinement steps to perform",
        RefineLoopMax);
    cmd.addReal( "T", "thickness", "thickness",  thickness );
    cmd.addInt( "g", "goal", "Goal function to use", goal );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    gsMultiPatch<> mp_ex;
    gsReadFile<>("deformed_plate_T=" + std::to_string(thickness) + ".xml",mp_ex);
    gsMultiBasis<> basisR(mp_ex);

    // Unit square
    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.addAutoBoundaries();
    mp.embed(3);
    E_modulus = 1.0;
    // thickness = 1.0;

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    gsMultiBasis<> dbasis(mp);
    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";


    // Cast all patches of the mp object to THB splines
    gsTHBSpline<2,real_t> thb;
    for (index_t k=0; k!=mp.nPatches(); ++k)
    {
        gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp.patch(k));
        thb = gsTHBSpline<2,real_t>(*geo);
        mp.patch(k) = thb;
    }

    mp_def = mp;

    gsMultiBasis<> basisL(mp);
    gsMultiBasis<> basisH = basisL;
    basisH.degreeElevate();

    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    gsVector<> neu(3);
    neu << 0, 0, 0;

    gsFunctionExpr<> displ("1",3);
    real_t load = 1.0;
    real_t D = E_modulus * math::pow(thickness,3) / ( 12 * ( 1- math::pow(PoissonRatio,2) ) );

    gsConstantFunction<> neuData(neu,3);

    for (index_t i=0; i!=3; ++i)
    {
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i ); // unknown 0 - x
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i ); // unknown 1 - y
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
    }
    tmp << 0,0,-load;
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> exL;
    gsExprAssembler<> exH;
    gsExprAssembler<> exRef;

    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    exL.setIntegrationElements(basisH);
    exH.setIntegrationElements(basisH);
    exRef.setIntegrationElements(basisR);
    gsExprEvaluator<> evL(exL);
    gsExprEvaluator<> evH(exH);
    gsExprEvaluator<> evRef(exRef);

    // Set the geometry map
    geometryMap mapL = exL.getMap(mp); // the last map counts
    geometryMap defL = exL.getMap(mp_def);
    geometryMap mapH = exH.getMap(mp); // the last map counts
    geometryMap defH = exH.getMap(mp_def);
    geometryMap mapRef = exRef.getMap(mp);
    geometryMap defRef = exRef.getMap(mp_ex);

    // Set the discretization spaces
    space uL = exL.getSpace(basisL, 3); //primal space on L
    space zH = exH.getSpace(basisH, 3); // dual space on H

    uL.setup(bc, dirichlet::interpolation, 0);
    zH.setup(bc, dirichlet::interpolation, 0);

    // Solution vector and solution variable
    gsMatrix<> random;
    solution uL_sol = exL.getSolution(uL,random);
    solution zL_sol = exL.getSolution(uL,random);
    solution zH_sol = exH.getSolution(zH,random);

    // zH2
    gsMultiPatch<> zH2_mp(mp);//just initialize for not being empty
    variable zH2 = exL.getCoeff(zH2_mp);
    // zL2
    gsMultiPatch<> zL2_mp(mp);//just initialize for not being empty
    variable zL2 = exH.getCoeff(zL2_mp);
    // uL2
    gsMultiPatch<> uL2_mp(mp);//just initialize for not being empty
    variable uL2 = exH.getCoeff(uL2_mp);

    gsMultiPatch<> uH2_mp(mp);//just initialize for not being empty
    variable uH2 = exL.getCoeff(uH2_mp);

    // gsFunctionExpr<> materialMat("1","0","0","0","1","0","0","0","1",3);
    // variable mm = A.getCoeff(materialMat, G);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsMaterialMatrix materialMat(mp, E, nu);
    variable mmL = exL.getCoeff(materialMat); // evaluates in the parametric domain, but the class transforms E and nu to physical
    variable mmH = exH.getCoeff(materialMat); // evaluates in the parametric domain, but the class transforms E and nu to physical
    variable mmRef = exRef.getCoeff(materialMat);

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",3);
    variable m2L = exL.getCoeff(mult2t, mapL); // evaluates in the physical domain
    variable m2H = exH.getCoeff(mult2t, mapH); // evaluates in the physical domain
    variable m2Ref = exRef.getCoeff(mult2t, mapRef);

    gsFunctionExpr<> t(std::to_string(thickness), 3);
    variable ttL = exL.getCoeff(t, mapL); // evaluates in the physical domain
    variable ttH = exH.getCoeff(t, mapH); // evaluates in the physical domain
    // TEMPORARILY!!
    // real_t tt = thickness;

    // gsFunctionExpr<> force("0","0","1", 3);
    gsConstantFunction<> force(tmp,3);
    variable ffL = exL.getCoeff(force,mapL); // evaluates in the physical domain
    variable ffH = exH.getCoeff(force,mapH); // evaluates in the physical domain

    gsSparseSolver<>::CGDiagonal solver;

    //! [Problem setup]

    //! [Solver loop]

    // Set Dirichlet values
    //A.options().setInt("DirichletValues", dirichlet::homogeneous);

    // Initialize the system
    exL.initSystem(false);
    exH.initSystem(false);

    gsInfo  <<"Lower order basis:\n"
            <<"\t Order: "<<basisL.maxCwiseDegree()<<"\n"
            <<"\t Number of elements: "<<basisL.totalElements()<<"\n"
            <<"\t Number of DoFs: "<<exL.numDofs()<<"\n";
    gsInfo  <<"Higher order basis:\n"
            <<"\t Order: "<<basisH.maxCwiseDegree()<<"\n"
            <<"\t Number of elements: "<<basisH.totalElements()<<"\n"
            <<"\t Number of DoFs: "<<exH.numDofs()<<"\n";

    /*
        We provide the following functions:                         checked with previous assembler
            E_m         membrane strain tensor.                             V
            E_m_der     first variation of E_m.                             V
            E_m_der2    second variation of E_m MULTIPLIED BY S_m.          V
            E_f         flexural strain tensor.                             V
            E_f_der     second variation of E_f.                            V
            E_f_der2    second variation of E_f MULTIPLIED BY S_f.          X NOTE: var1 in E_f_der2 in previous assembler is computed with both G and defG

        Where:
            G       the undeformed geometry,
            defG    the deformed geometry,
            mm      the material matrix,
            m2      an auxillary matrix to multiply the last row of a tensor with 2
    **/

    // --------------------------------------------------Lower-order basis----------------------------------------------------------------- //
    // Membrane components
    auto E_mL = 0.5 * ( flat(jac(defL).tr()*jac(defL)) - flat(jac(mapL).tr()* jac(mapL)) ) ; //[checked]
    auto S_mL = E_mL * reshape(mmL,3,3);
    auto N_L       = ttL.val() * S_mL;

    auto E_m_derL = flat( jac(defL).tr() * jac(uL) ) ; //[checked]
    auto S_m_derL = E_m_derL * reshape(mmL,3,3);
    auto N_derL   = ttL.val() * S_m_derL;

    auto E_m_der2_L = flatdot( jac(uL),jac(uL).tr(), N_L ); //[checked]

    // Flexural components
    auto E_fL = ( deriv2(mapL,sn(mapL).normalized().tr()) - deriv2(defL,sn(defL).normalized().tr()) ) * reshape(m2L,3,3) ; //[checked]
    auto S_fL = E_fL * reshape(mmL,3,3);
    auto M_L   = ttL.val() * ttL.val() * ttL.val() / 12.0 * S_fL;

    auto E_f_derL = ( deriv2(uL,sn(defL).normalized().tr() ) + deriv2(defL,var1(uL,defL) ) ) * reshape(m2L,3,3); //[checked]
    auto S_f_derL = E_f_derL * reshape(mmL,3,3);
    auto M_derL   = ttL.val() * ttL.val() * ttL.val() / 12.0 * S_f_derL;

    auto E_f_der2_L = flatdot2( deriv2(uL), var1(uL,defL).tr(), M_L  ).symmetrize() + var2(uL,uL,defL, M_L );

    auto F_L        = ffL;

    // --------------------------------------------------Higher-order basis---------------------------------------------------------------- //
    // Membrane components
    auto E_mH = 0.5 * ( flat(jac(defH).tr()*jac(defH)) - flat(jac(mapH).tr()* jac(mapH)) ) ; //[checked]
    auto S_mH = E_mH * reshape(mmH,3,3);
    auto N_H       = ttH.val() * S_mH;

    auto E_m_derH = flat( jac(defH).tr() * jac(zH) ) ; //[checked]
    auto S_m_derH = E_m_derH * reshape(mmH,3,3);
    auto N_derH   = ttH.val() * S_m_derH;

    auto E_m_der2_H = flatdot( jac(zH),jac(zH).tr(), N_H ); //[checked]

    // Flexural components
    auto E_fH = ( deriv2(mapH,sn(mapH).normalized().tr()) - deriv2(defH,sn(defH).normalized().tr()) ) * reshape(m2H,3,3) ; //[checked]
    auto S_fH = E_fH * reshape(mmH,3,3);
    auto M_H   = ttH.val() * ttH.val() * ttH.val() / 12.0 * S_fH;

    auto E_f_derH = ( deriv2(zH,sn(defH).normalized().tr() ) + deriv2(defH,var1(zH,defH) ) ) * reshape(m2H,3,3); //[checked]
    auto S_f_derH = E_f_derH * reshape(mmH,3,3);
    auto M_derH   = ttH.val() * ttH.val() * ttH.val() / 12.0 * S_f_derH;

    auto E_f_der2_H = flatdot2( deriv2(zH), var1(zH,defH).tr(), M_H  ).symmetrize() + var2(zH,zH,defH, M_H );

    auto F_H        = ffH;

    // ! [Solve linear problem]

    // [Pre-work]


    gsMatrix<> solVectorL;
    uL_sol.setSolutionVector(solVectorL);

    gsMatrix<> solVectorDualL;
    zL_sol.setSolutionVector(solVectorDualL);

    gsMatrix<> solVectorDualH;
    std::vector<bool> elMarked;
    // So, ready to start the adaptive refinement loop:
    for( int RefineLoop = 1; RefineLoop <= RefineLoopMax ; RefineLoop++ )
    {
        gsInfo  << "\n ====== Loop " << RefineLoop << " of " << RefineLoopMax << " ======" << "\n" << "\n"
                << "Number of elements: "<<basisL.totalElements()<<"\n";

        // --------------------------------------------------------------------------------------------------------------------------------- //
        // ---------------------------------------------------------Solving primal and dual problems---------------------------------------- //
        // --------------------------------------------------------------------------------------------------------------------------------- //

        mp_def = mp;

        // Redefine basis
        basisL = gsMultiBasis<>(mp);
        basisH = basisL;
        basisH.degreeElevate(1);

        // Assemble matrix and rhs
        exL.initSystem(false);
        gsInfo << "Assembling primal, size ="<<exL.matrix().rows()<<","<<exL.matrix().cols()<<"... "<< std::flush;

        exL.assemble(
                // (N_der * cartcon(G) * (E_m_der * cartcon(G)).tr() + M_der * cartcon(G) * (E_f_der * cartcon(G)).tr()) * meas(G)
                (N_derL * (E_m_derL).tr() + M_derL * (E_f_derL).tr()) * meas(mapL)
                ,
                uL * F_L  * meas(mapL)
            );
        gsInfo << "done." << "\n";

        // Solve system
        gsInfo << "Solving primal, size ="<<exL.matrix().rows()<<","<<exL.matrix().cols()<<"... "<< std::flush;
        solver.compute(exL.matrix());
        solVectorL = solver.solve(exL.rhs());

        gsInfo << "done." << " --> ";

        gsInfo << "Solving primal (nonlinear), size ="<<exL.matrix().rows()<<","<<exL.matrix().cols()<<"... "<< "\n";

        // Deform mps
        gsMatrix<> cc;
        for ( size_t k =0; k!=mp.nPatches(); ++k) // Deform the geometry
        {
            // extract deformed geometry
            uL_sol.extract(cc, k);
            mp_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
        }

        index_t itMax = 100;
        real_t tol = 1e-14;
        real_t residual = exL.rhs().norm();
        real_t residual0 = residual;
        real_t residualOld = residual;
        for (index_t it = 0; it != itMax; ++it)
        {
            exL.initSystem(false);
            // assemble system
            exL.assemble(
                ( N_derL * E_m_derL.tr() + E_m_der2_L + M_derL * E_f_derL.tr() - E_f_der2_L ) * meas(mapL)
                , uL * F_L * meas(mapL) - ( ( N_L * E_m_derL.tr() - M_L * E_f_derL.tr() ) * meas(mapL) ).tr()
                );
            solver.compute( exL.matrix() );
            gsMatrix<> updateVectorL = solver.solve(exL.rhs()); // this is the UPDATE
            solVectorL += updateVectorL;
            residual = exL.rhs().norm();

            gsInfo<<"Iteration: "<< it
                   <<", residue: "<< residual
                   <<", update norm: "<<updateVectorL.norm()
                   <<", log(Ri/R0): "<< math::log10(residualOld/residual0)
                   <<", log(Ri+1/R0): "<< math::log10(residual/residual0)
                   <<"\n";

            residualOld = residual;


            // update deformed patch
            uL_sol.setSolutionVector(updateVectorL);
            for ( size_t k =0; k!=mp_def.nPatches(); ++k) // Deform the geometry
            {
                // // extract deformed geometry
                uL_sol.extract(cc, k);
                mp_def.patch(k).coefs() += cc;  // defG points to mpL_def, therefore updated
            }


            if (updateVectorL.norm() < tol)
                break;
        }

        uL_sol.setSolutionVector(solVectorL);
        uL_sol.extract(uL2_mp);
        uL_sol.extract(uH2_mp);
        gsInfo << "done." << " --> ";

        // Assemble matrix and rhs
        exL.initSystem(false);
        exL.assemble( ( N_derL * E_m_derL.tr() + E_m_der2_L + M_derL * E_f_derL.tr() - E_f_der2_L ) * meas(mapL) );
        exL.initVector(1,false);
        gsInfo << "Assembling dual (low), size = "<<exL.matrix().rows()<<","<<exL.matrix().cols()<<"... "<< std::flush;
        // NOTE, we assume that the matrix in space uL is equal to that in space zL, hence it is not re-assembled!
        if (goal == 1)
            exL.assemble( uL * gismo::expr::uv(2,3) * meas(mapL) );
        else if (goal == 2)
            exL.assemble( 2 * uL * uL_sol * meas(mapL) );
        else if (goal == 3)
            exL.assemble( 2 * E_m_derL * E_mL.tr() * meas(mapL) );
        else if (goal == 4)
            exL.assemble( S_m_derL * gismo::expr::uv(0,3) * meas(mapL) );
        gsInfo << "done." << "\n";


        // Solve system
        gsInfo << "Solving dual (low), size = "<<exL.matrix().rows()<<","<<exL.matrix().cols()<<"... "<< std::flush;
        gsSparseMatrix<> matrixL = exL.matrix().transpose();
        gsDebug<<"matrix is transposed\n";
        solver.compute(matrixL); // not needed
        solVectorDualL = solver.solve(exL.rhs());
        // gsDebugVar(exL.rhs().transpose());
        // gsDebugVar(exL.matrix().toDense());

        zL_sol.setSolutionVector(solVectorDualL);
        zL_sol.extract(zL2_mp);

        gsInfo << "done." << " --> ";

        // Assemble matrix and rhs
        exH.initSystem(false);
        gsInfo << "Assembling dual (high), size = "<<exH.matrix().rows()<<","<<exH.matrix().cols()<<"... "<< std::flush;
        exH.assemble( ( N_derH * E_m_derH.tr() + E_m_der2_H + M_derH * E_f_derH.tr() - E_f_der2_H ) * meas(mapH) );

        if (goal == 1)
            exH.assemble( zH * gismo::expr::uv(2,3) * meas(mapH) );
        else if (goal == 2)
            exH.assemble( 2 * zH * uL2 * meas(mapH) );
        else if (goal == 3)
            exH.assemble( 2*E_m_derH * E_mH.tr() * meas(mapH) );
        else if (goal == 4)
            exH.assemble( S_m_derH * gismo::expr::uv(0,3) * meas(mapH) );

        gsInfo << "done." << "\n";

        // Solve system
        gsInfo << "Solving dual (high), size = "<<exH.matrix().rows()<<","<<exH.matrix().cols()<<"... "<< std::flush;
        gsSparseMatrix<> matrixH = exH.matrix().transpose();
        gsDebug<<"matrix is transposed\n";
        solver.compute(matrixH);
        solVectorDualH = solver.solve(exH.rhs());
        // gsDebugVar(exH.matrix().toDense());

        zH_sol.setSolutionVector(solVectorDualH);
        zH_sol.extract(zH2_mp);
        gsInfo << "done." << " --> \n";

        // --------------------------------------------------------------------------------------------------------------------------------- //
        // ---------------------------------------------------------Computing DWR error estimate-------------------------------------------- //
        // --------------------------------------------------------------------------------------------------------------------------------- //

        auto E_m = 0.5 * ( flat(jac(defL).tr()*jac(defL)) - flat(jac(mapL).tr()* jac(mapL)) ) ; //[checked]
        auto S_m = E_m * reshape(mmL,3,3);
        auto N   = ttL.val() * S_m;

        auto E_f = ( deriv2(mapL,sn(mapL).normalized().tr()) - deriv2(defL,sn(defL).normalized().tr()) ) * reshape(m2L,3,3) ; //[checked]
        auto S_f = E_f * reshape(mmL,3,3);
        auto M   = ttL.val() * ttL.val() * ttL.val() / 12.0 * S_f;

        auto E_m_der = flat( jac(defL).tr() * (fjac(zH2) - jac(zL_sol)) ) ; //[checked]
        auto S_m_der = E_m_der * reshape(mmL,3,3);
        auto N_der   = ttL.val() * S_m_der;

        auto E_f_der = ( deriv2(zH2,sn(defL).normalized().tr() ) - deriv2(zL_sol,sn(defL).normalized().tr() ) + deriv2(defL,var1(zH2,defL) ) - deriv2(defL,var1(zL_sol,defL) ) ) * reshape(m2L,3,3); //[checked]
        auto S_f_der = E_f_der * reshape(mmL,3,3);
        auto M_der   = ttL.val() * ttL.val() * ttL.val() / 12.0 * S_f_der;

        auto Fint_m = ( N * E_m_der.tr() ) * meas(mapL);
        auto Fint_f = ( M * E_f_der.tr() ) * meas(mapL);
        auto Fint = ( N * E_m_der.tr() - M * E_f_der.tr() ) * meas(mapL);
        auto Fext = ( zH2.tr() * ffL - zL_sol.tr() * ffL ) * meas(mapL);

        auto E_mG = 0.5 * ( flat(jac(defRef).tr()*jac(defRef)) - flat(jac(mapRef).tr()* jac(mapRef)) ) ; //[checked]
        auto E_fG = ( deriv2(mapRef,sn(mapRef).normalized().tr()) - deriv2(defRef,sn(defRef).normalized().tr()) ) * reshape(m2Ref,3,3) ; //[checked]
        auto S_mG = E_mG * reshape(mmRef,3,3);
        auto S_fG = E_fG * reshape(mmRef,3,3);

        real_t fint_m = evL.integral(( N * E_m_der.tr() ) * meas(mapL) );
        real_t fint_f = evL.integral(( M * E_f_der.tr() ) * meas(mapL) );
        real_t fint = fint_m - fint_f;
        real_t fext = evL.integral( Fext  );

        gsInfo<<"Fint_m = "<<fint_m<<"\n";
        gsInfo<<"Fint_f = "<<fint_f<<"\n";

        gsInfo<<"Fint = "<<fint<<"\n";
        gsInfo<<"Fext = "<<fext<<"\n";
        real_t Res = fext-fint;
        gsInfo<<"R = "<<Res<<"\n";

        real_t approx = Res;
        real_t exact = 0;
        if (goal==1)
            exact = evRef.integral((defRef-mapRef).tr() * gismo::expr::uv(2,3)*meas(mapRef)) - evL.integral((defL-mapL).tr() * gismo::expr::uv(2,3)*meas(mapL));
        else if (goal==2)
            exact = evRef.integral( (defRef - mapRef).tr() * (defRef - mapRef) * meas(mapRef)) - evL.integral((defL - mapL).tr() * (defL - mapL) * meas(mapL));
        else if (goal==3)
            exact = evRef.integral(E_mG * E_mG.tr() *meas(mapRef) ) - evL.integral(E_m * E_m.tr() *meas(mapL) );
        else if (goal==4)
            exact = evRef.integral(S_mG * gismo::expr::uv(0,3)*meas(mapRef) ) - evL.integral(S_m * gismo::expr::uv(0,3)*meas(mapL) );


        gsInfo<<"Exact = "<<exact<<"\n";
        gsInfo<<"Efficiency = "<<approx/exact<<"\n";
    }

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsMatrix<> cc;
        uL_sol.setSolutionVector(solVectorL);
        mp_def = mp;

        for ( size_t k =0; k!=mp.nPatches(); ++k) // Deform the geometry
        {
            // extract deformed geometry
            uL_sol.extract(cc, k);
            mp_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
        }
        gsWriteParaview<>( mp_def, "mp_def", 1000, true);

        gsMultiPatch<> deformation = mp_def;
        for (index_t k = 0; k != mp_def.nPatches(); ++k)
            deformation.patch(k).coefs() -= mp.patch(k).coefs();

        gsInfo<<"Plotting in Paraview...\n";
        evL.options().setSwitch("plot.elements", true);
        evH.options().setSwitch("plot.elements", true);

        evL.writeParaview( defL-mapL   , mapL, "solution_primalL");
        evL.writeParaview( zL_sol   , mapL, "solution_dualL");
        evH.writeParaview( zH_sol   , mapH, "solution_dualH");

        evL.writeParaview( zH2-zL_sol   , mapL, "solution_dual");

        gsWriteParaview<>( basisL.basis(0), "basisL", 1000);
        gsWriteParaview<>( basisH.basis(0), "basisH", 1000);


        // ev.options().setSwitch("plot.elements", true);
        // ev.writeParaview( S_f2, G, "stress");
        // evaluateFunction(ev, S_f2[0], pt); // evaluates an expression on a point


        // gsFileManager::open("solution.pvd");
    }

    return EXIT_SUCCESS;

}// end main

template <class T>
void evaluateFunction(gsExprEvaluator<T> ev, auto expression, gsVector<T> pt)
{
    gsMatrix<T> evresult = ev.eval( expression,pt );
    gsInfo << "Eval on point ("<<pt.at(0)<<" , "<<pt.at(1)<<") :\n"<< evresult;
    gsInfo << "\nEnd ("<< evresult.rows()<< " x "<<evresult.cols()<<")\n";
};

template <class T>
void evaluateFunction(gsExprEvaluator<T> ev, auto expression, gsMatrix<T> pt)
{
    gsVector<T> tmp;
    gsMatrix<T> evresult;
    for (index_t k=0; k!= pt.cols(); ++k)
    {
        tmp = pt.col(k);
        evresult = ev.eval( expression,tmp );
        gsInfo<<evresult<<"\n";
    }
};

/*
template<class T>
void gsShellAssembler<T>::applyLoads()
{
    gsMatrix<T>        bVals;
    gsMatrix<unsigned> acts,globalActs;

    for (size_t i = 0; i< m_pLoads.numLoads(); ++i )
    {
        if ( m_pLoads[i].parametric )
        {
            m_bases.front().basis(m_pLoads[i].patch).active_into( m_pLoads[i].point, acts );
            m_bases.front().basis(m_pLoads[i].patch).eval_into  ( m_pLoads[i].point, bVals);
        }
        else
        {
            gsMatrix<> forcePoint;
            m_patches.patch(m_pLoads[i].patch).invertPoints(m_pLoads[i].point,forcePoint);
            u.source().piece(m_pLoads[i].patch).active_into( forcePoint, acts );
            u.source().piece(m_pLoads[i].patch).active_into( forcePoint, bVals);
        }

        // translate patch-local indices to global dof indices
        for (size_t j = 0; j< 3; ++j)
        {
            if (m_pLoads[i].value[j] != 0.0)
            {
                u.dofMappers[j].localToGlobal(acts, m_pLoads[i].patch, globalActs);

                for (index_t k=0; k < globalActs.rows(); ++k)
                {
                    if (int(globalActs(k,0)) < m_dofs)
                        m_rhs(globalActs(k,0), 0) += bVals(k,0) * m_pLoads[i].value[j];
                }
            }
        }
    }
}
*/