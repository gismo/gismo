/** @file kirchhoff-Love_example.cpp

    @brief Solver for kirchhoff-Love shells

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst & A. Mantzaflaris
*/

//! [Include namespace]

#include <gismo.h>

namespace gismo{
namespace expr{

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
    index_t cols() const { return _u.dim(); }

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

    enum{rowSpan = E::rowSpan, colSpan = 0};

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

        for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all actives
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
    typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");
        solGrad_expr<Scalar> sGrad =  solGrad_expr<Scalar>(_u);
        res.resize(rows(), cols()); // rows()*

        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();
        bGrads = sGrad.eval(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        m_v.noalias() = ( ( bGrads.col(0).template head<3>() ).cross( cJac.col(1).template head<3>() )
                      -   ( bGrads.col(1).template head<3>() ).cross( cJac.col(0).template head<3>() ) ) / measure;

        // ---------------  First variation of the normal
        // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
        res = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
        return res;
    }
};

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

    enum{rowSpan = 1, colSpan = 1};

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

    enum{rowSpan = E1::rowSpan, colSpan = E2::rowSpan};

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

        solHess_expr<Scalar> sHess = solHess_expr<Scalar>(_u);
        tmp = sHess.eval(k);
        vEv = _v.eval(k);
        res = vEv * tmp;
        return res;
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsGeometryMap<Scalar> >::value, index_t >::type
    cols_impl(const U & u)  const
    {
        return _u.data().dim.second;
    }

    template<class U> inline
    typename util::enable_if<util::is_same<U,gsFeSpace<Scalar> >::value, index_t >::type
    cols_impl(const U & u) const
    {
        return _u.dim();
    }

    template<class U> inline
    typename util::enable_if<util::is_same<U,gsFeSolution<Scalar> >::value, index_t >::type
    cols_impl(const U & u) const
    {
        return _u.dim();
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
    enum {ColBlocks = E::rowSpan };
    enum{ Space = E::Space };

    deriv2_expr(const E & u) : _u(u) { }

    mutable gsMatrix<Scalar> res, tmp;

    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,k); }

    index_t rows() const //(components)
    {
        return 3; // _u.dim() for space or targetDim() for geometry
    }

    index_t cols() const
    {
        return _u.source().domainDim() * ( _u.source().domainDim() + 1 ) / 2;
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

    enum{rowSpan = E::rowSpan, colSpan = E::colSpan};

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

            return res;
        }

};


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

        GISMO_ASSERT(Bc==_A.rows(), "Dimensions: "<<Bc<<","<< _A.rows()<< "do not match");
        GISMO_STATIC_ASSERT(E1::rowSpan, "First entry should be rowSpan"  );
        GISMO_STATIC_ASSERT(E2::colSpan, "Second entry should be colSpan.");

        res.resize(An, Bn);

        // to do: evaluate the jacobians internally for the basis functions and loop over dimensions as well, since everything is same for all dimensions.
        for (index_t i = 0; i!=An; ++i) // for all actives u
            for (index_t j = 0; j!=Bn; ++j) // for all actives v
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

    enum{rowSpan = E1::rowSpan, colSpan = E2::colSpan};

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
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t Ac = _A.cols();
        const index_t An = _A.cardinality();
        const index_t Bn = _B.cardinality();

        eA = _A.eval(k);
        eB = _B.eval(k);
        eC = _C.eval(k);

        GISMO_ASSERT(_B.rows()==_A.cols(), "Dimensions: "<<_B.rows()<<","<< _A.cols()<< "do not match");
        GISMO_STATIC_ASSERT(E1::rowSpan, "First entry should be rowSpan");
        GISMO_STATIC_ASSERT(E2::colSpan, "Second entry should be colSpan.");
        GISMO_ASSERT(_C.cols()==_B.rows(), "Dimensions: "<<_C.rows()<<","<< _B.rows()<< "do not match");

        res.resize(An, Bn);
        for (index_t i = 0; i!=An; ++i)
            for (index_t j = 0; j!=Bn; ++j)
            {
                tmp = eA.middleCols(i*Ac,Ac) * eB.col(j);   // E_f_der2
                tmp.row(2) *= 2.0;                          // multiply the third row of E_f_der2 by 2 for voight notation
                res(i,j) = (eC.row(0) * tmp.col(0)).value();          // E_f^T * mm * E_f_der2
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

    enum{rowSpan = E1::rowSpan, colSpan = E2::colSpan};

    void print(std::ostream &os) const { os << "flatdot2("; _A.print(os);_B.print(os);_C.print(os); os<<")"; }
};

/*
   Expression for the transformation matrix FROM local covariant TO local cartesian bases, based on a geometry map
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

    gsMatrix<Scalar,3,3> eval(const index_t k) const
    {
        // Compute covariant bases in deformed and undeformed configuration
        normal = _G.data().normals.col(k);
        normal.normalize();
        covBasis.leftCols(2) = _G.data().jacobian(k);
        covBasis.col(2)      = normal;
        covMetric = covBasis.transpose() * covBasis;

        conMetric = covMetric.inverse();

        conBasis.col(1) = conMetric(1,0)*covBasis.col(0)+conMetric(1,1)*covBasis.col(1)+conMetric(1,2)*covBasis.col(2);

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

    enum{rowSpan = 0, colSpan = 0};

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

    gsMatrix<T> eval(const index_t k) const
    {
        temp = (cartcov_expr<gsGeometryMap<T> >(_G).eval(k)).reshape(3,3).inverse();
        return temp;
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    enum{rowSpan = 0, colSpan = 0};

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
   Expression for the transformation matrix FROM local contravariant TO local cartesian bases, based on a geometry map
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

    gsMatrix<Scalar,3,3> eval(const index_t k) const
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

    enum{rowSpan = 0, colSpan = 0};

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

    gsMatrix<T> eval(const index_t k) const
    {
        temp = (cartcon_expr<gsGeometryMap<T> >(_G).eval(k)).reshape(3,3).inverse();
        return temp;
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 3; }

    enum{rowSpan = 0, colSpan = 0};

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

template<class E> EIGEN_STRONG_INLINE
var1_expr<E> var1(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return var1_expr<E>(u, G); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
var2_expr<E1,E2,E3> var2(const E1 & u, const E2 & v, const gsGeometryMap<typename E1::Scalar> & G, const E3 & Ef)
{ return var2_expr<E1,E2,E3>(u,v, G, Ef); }

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
//! [Include namespace]
using namespace gismo;
//! [Include namespace]

// Input is parametric coordinates of the surface \a mp
template <class T>
class gsMaterialMatrix : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
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
        }
    }

    // piece(k) --> for patch k

};

//! [Forward declarations]
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B);
template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B);
//! [Forward declarations]


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool weak = false;
    bool nonlinear = false;
    std::string fn;

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t thickness = 1.0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("weak", "Weak BCs", weak);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [set test case data]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;

    if (testCase == 2  || testCase == 3)
    {
        thickness = 0.25;
        E_modulus = 4.32E8;
        fn = "surfaces/scordelis_lo_roof.xml";
        gsReadFile<>(fn, mp);
        PoissonRatio = 0.0;
    }
    else if (testCase == 6)
    {
        thickness = 0.1;
        real_t mu = 4.225e5;
        PoissonRatio = 0.3;
        E_modulus = (2+PoissonRatio)*mu;
        fn = "surfaces/quarter_sphere.xml";
        gsReadFile<>(fn, mp);
    }
    else
    {
        // Unit square
        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.addAutoBoundaries();
        mp.embed(3);
        E_modulus = 1e0;
        thickness = 1e0;
        PoissonRatio = 0.0;
    }
    //! [set test case data]

    //! [Refinement]
    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();
    mp_def = mp;
    gsWriteParaview(mp,"mp",1000,true);

    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";
    //! [Refinement]

    //! [pre-define boundary conditions]
    gsBoundaryConditions<> bc;
    bc.setGeoMap(mp);
    gsVector<> tmp(3);
    tmp << 0, 0, 0;

    gsVector<> neu(3);
    neu << 0, 0, 0;

    gsConstantFunction<> displx(0.09317696,3);

    gsConstantFunction<> neuData(neu,3);
    gsConstantFunction<> weakBC(neu,3);

    real_t pressure = 0.0;
    //! [pre-define boundary conditions]

    //! [Boundary condition case 1]
    if (testCase == 1)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i ); // unknown 0 - x
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }
        tmp << 0,0,-1;
    }
    //! [Boundary condition case 1]
    //! [Boundary condition case 2]
    else if (testCase == 2)
    {
        // Diaphragm conditions
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        bc.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Surface forces
        tmp << 0, 0, -90;
    }
    //! [Boundary condition case 2]
    else if (testCase == 3)
    {
        neu << 0, 0, -90;
        neuData.setValue(neu,3);
        // Diaphragm conditions
        bc.addCondition(boundary::west, condition_type::neumann, &neuData );

        // ORIGINAL
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 1 ); // unknown 1 - y
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // NOT ORIGINAL
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 ); // unknown 1 - x

        // Surface forces
        tmp << 0, 0, 0;
    }
    else if (testCase == 4)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i ); // unknown 0 - x
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
            if (!weak)
                bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 0 - x

        }
        if (weak)
            bc.addCondition(boundary::west, condition_type::weak_dirichlet, &weakBC ); // unknown 2 - z

        tmp << 0, 0, -1;
    }

    else if (testCase == 6)
    {
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        // Symmetry in x-direction:
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, 0 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 1 );
        bc.addCondition(boundary::east, condition_type::clamped, 0, 0, false, 2 );

        // Symmetry in y-direction:
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 0 );
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, 1 );
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0, false, 2 );

        // Pressure
        pressure = 5e2;
    }

    else if (testCase == 10)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }
        bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z
        bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, 2 ); // unknown 2 - z

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,1);
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ,false,2);

        if (!weak)
            bc.addCondition(boundary::east, condition_type::dirichlet, &displx, 0 ,false,0);
        else
        {
            neu << -0.09317696, 0, 0;
            weakBC.setValue(neu,3);
            bc.addCondition(boundary::east, condition_type::weak_dirichlet, &weakBC);
        }

        // Surface forces
        tmp.setZero();
    }
    else if (testCase == 11)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i ); // unknown 1 - y
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }

        bc.addCondition(boundary::east, condition_type::clamped, 0, 0 ,false,2);
        bc.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);

        tmp<<0,0,-1;
    }
    else if (testCase == 12)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
        }

        bc.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);

        neu << 0, 0, -0.1;
        neuData.setValue(neu,3);

        bc.addCondition(boundary::east, condition_type::neumann, &neuData );
    }
    else if (testCase == 13)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, i ); // unknown 0 - x
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0, false, i ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, 0, false, i ); // unknown 2 - z
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0, false, i ); // unknown 0 - x
        }

        if (!weak)
        {
            bc.addCondition(boundary::north, condition_type::clamped, 0, 0 ,false,2);
            bc.addCondition(boundary::east, condition_type::clamped, 0, 0 ,false,2);
            bc.addCondition(boundary::south, condition_type::clamped, 0, 0 ,false,2);
            bc.addCondition(boundary::west, condition_type::clamped, 0, 0 ,false,2);
        }
        else
        {
            bc.addCondition(boundary::north, condition_type::weak_clamped, &weakBC );
            bc.addCondition(boundary::east, condition_type::weak_clamped, &weakBC );
            bc.addCondition(boundary::south, condition_type::weak_clamped, &weakBC );
            bc.addCondition(boundary::west, condition_type::weak_clamped, &weakBC ); // unknown 2 - z
        }


        tmp << 0,0,-1;
    }

    //! [Assembler setup]
    gsExprAssembler<> A;

    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp); // the last map counts
    geometryMap defG = A.getMap(mp_def);


    // Set the discretization space
    space u = A.getSpace(dbasis, 3);
    u.setup(bc, dirichlet::interpolation, 0);

    // Solution vector and solution variable
    gsMatrix<> random;
    solution u_sol = A.getSolution(u,random);

    // gsFunctionExpr<> materialMat("1","0","0","0","1","0","0","0","1",3);
    gsFunctionExpr<> E(util::to_string(E_modulus),3);
    gsFunctionExpr<> nu(util::to_string(PoissonRatio),3);
    gsMaterialMatrix<real_t> materialMat(mp, E, nu);
    variable mm = A.getCoeff(materialMat); // evaluates in the parametric domain, but the class transforms E and nu to physical

    gsFunctionExpr<> t(util::to_string(thickness), 3);
    variable tt = A.getCoeff(t, G); // evaluates in the physical domain

    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",3);
    variable m2 = A.getCoeff(mult2t, G); // evaluates in the physical domain

    gsConstantFunction<> force(tmp,3);
    variable ff = A.getCoeff(force,G); // evaluates in the physical domain

    gsSparseSolver<>::CGDiagonal solver;

    //! [Assembler setup]

    //! [System assembly]

    // Initialize the system
    A.initSystem(false);

    gsInfo<<"Number of degrees of freedom: "<< A.numDofs() <<"\n"<<std::flush;

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


    // Membrane components
    auto E_m = 0.5 * ( flat(jac(defG).tr()*jac(defG)) - flat(jac(G).tr()* jac(G)) ) ; //[checked]
    auto S_m = E_m * reshape(mm,3,3);
    auto N   = tt.val() * S_m;

    auto E_m_der = flat( jac(defG).tr() * jac(u) ) ; //[checked]
    auto S_m_der = E_m_der * reshape(mm,3,3);
    auto N_der   = tt.val() * S_m_der;

    auto E_m_der2 = flatdot( jac(u),jac(u).tr(), N ); //[checked]

    // Flexural components
    auto E_f = ( deriv2(G,sn(G).normalized().tr()) - deriv2(defG,sn(defG).normalized().tr()) ) * reshape(m2,3,3) ; //[checked]
    auto S_f = E_f * reshape(mm,3,3);
    auto M   = tt.val() * tt.val() * tt.val() / 12.0 * S_f;

    auto E_f_der = -( deriv2(u,sn(defG).normalized().tr() ) + deriv2(defG,var1(u,defG) ) ) * reshape(m2,3,3); //[checked]
    auto S_f_der = E_f_der * reshape(mm,3,3);
    auto M_der   = tt.val() * tt.val() * tt.val() / 12.0 * S_f_der;

    auto E_f_der2 = - (flatdot2( deriv2(u), var1(u,defG).tr(), M  ).symmetrize() + var2(u,u,defG, M ));

    auto F        = ff;

    auto That   = cartcon(G);
    auto Ttilde = cartcov(G);
    auto E_m_plot = 0.5 * ( flat(jac(defG).tr()*jac(defG)) - flat(jac(G).tr()* jac(G)) ) * That; //[checked]
    auto S_m_plot = E_m_plot * reshape(mm,3,3) * Ttilde; //[checked]

    // // For Neumann (same for Dirichlet/Nitsche) conditions
    variable g_N = A.getBdrFunction();

    real_t alpha_d = 1e3;
    A.assembleLhsRhsBc
    (
        alpha_d * u * u.tr()
        ,
        alpha_d * (u * (defG - G) - u * (g_N) )
        ,
        bc.container("Weak Dirichlet")
    );

    // for weak clamped
    real_t alpha_r = 1e3;
    A.assembleLhsRhsBc
    (
        alpha_r * ( sn(defG).tr()*nv(G) - sn(G).tr()*nv(G) ).val() * ( var2(u,u,defG,nv(G).tr()) )
        +
        alpha_r * ( ( var1(u,defG) * nv(G) ) * ( var1(u,defG) * nv(G) ).tr() )
        ,
        alpha_r * ( sn(defG).tr()*sn(G) - sn(G).tr()*sn(G) ).val() * ( var1(u,defG) * sn(G) )
        ,
        bc.container("Weak Clamped")
    );

    // For Neumann conditions
    A.assembleRhsBc(u * g_N * tv(G).norm(), bc.neumannSides() );

    A.assemble(
        (N_der * (E_m_der).tr() + M_der * (E_f_der).tr()) * meas(G)
        ,
        u * F  * meas(G) + pressure * u * sn(defG).normalized() * meas(G)
        );

    //! [System assembly]

    //! [Linear solve]
    // solve system
    solver.compute( A.matrix() );
    gsMatrix<> solVector = solver.solve(A.rhs());

    // update deformed patch
    gsMatrix<> cc;

    u_sol.setSolutionVector(solVector);
    for ( size_t k =0; k!=mp_def.nPatches(); ++k) // Deform the geometry
    {
        // // extract deformed geometry
        u_sol.extract(cc, k);
        mp_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
    }

    gsMatrix<> result;
    u_sol.extractFull(result);

    gsVector<> pt(2);
    pt.setConstant(0.25);
    ev.options().setInt("plot.npts", 5000);
    ev.writeParaview(S_m_plot,defG,"stress");
    gsInfo<<"Stresses\n"<<ev.eval(E_m_plot,pt)<<"\n";
    //! [Linear solve]

    //! [Nonlinear solve]
    real_t residual = A.rhs().norm();
    real_t residual0 = residual;
    real_t residualOld = residual;
    gsMatrix<> updateVector = solVector;
    if (nonlinear)
    {
        index_t itMax = 25;
        real_t tol = 1e-8;
        for (index_t it = 0; it != itMax; ++it)
        {

            A.initSystem(false);
            // assemble system
            A.assemble(
                (
                  N_der * E_m_der.tr()
                    +
                  E_m_der2
                    +
                  M_der * E_f_der.tr()
                    +
                  E_f_der2
                  ) * meas(G)
                    -
                  pressure * u * var1(u,defG) .tr() * meas(G)
                , u * F * meas(G) + pressure * u * sn(defG).normalized() * meas(G) - ( ( N * E_m_der.tr() + M * E_f_der.tr() ) * meas(G) ).tr()
                );

            A.assembleRhsBc(u * g_N * tv(G).norm(), bc.neumannSides() );

            A.assembleLhsRhsBc
            (
                alpha_d * u * u.tr()
                ,
                -alpha_d * (u * (defG - G) - u * (g_N) )
                ,
                bc.container("Weak Dirichlet")
            );
            A.assembleLhsRhsBc
            (
                alpha_r * ( sn(defG).tr()*nv(G) - sn(G).tr()*nv(G) ).val() * ( var2(u,u,defG,nv(G).tr()) )
                +
                alpha_r * ( ( var1(u,defG) * nv(G) ) * ( var1(u,defG) * nv(G) ).tr() )
                ,
                -alpha_r * ( sn(defG).tr()*nv(G) - sn(G).tr()*nv(G) ).val() * ( var1(u,defG) * nv(G) )
                ,
                bc.container("Weak Clamped")
            );


            // solve system
            solver.compute( A.matrix() );
            updateVector = solver.solve(A.rhs()); // this is the UPDATE


            solVector += updateVector;
            residual = A.rhs().norm();

            gsInfo<<"Iteration: "<< it
                   <<", residue: "<< residual
                   <<", update norm: "<<updateVector.norm()
                   <<", log(Ri/R0): "<< math::log10(residualOld/residual0)
                   <<", log(Ri+1/R0): "<< math::log10(residual/residual0)
                   <<"\n";

            residualOld = residual;

            // update deformed patch
            u_sol.setSolutionVector(updateVector);
            for ( size_t k =0; k!=mp_def.nPatches(); ++k) // Deform the geometry
            {
                // // extract deformed geometry
                u_sol.extract(cc, k);
                mp_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
            }

            if (residual < tol)
                break;
        }
    }
    //! [Nonlinear solve]

    //! [Construct solution]
    u_sol.setSolutionVector(solVector);
    mp_def = mp;
    for ( size_t k =0; k!=mp.nPatches(); ++k) // Deform the geometry
    {
        u_sol.extract(cc, k);
        mp_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
    }

    gsMultiPatch<> deformation = mp_def;
    for (size_t k = 0; k != mp_def.nPatches(); ++k)
        deformation.patch(k).coefs() -= mp.patch(k).coefs();

    gsInfo <<"Maximum deformation coef: "
           << deformation.patch(0).coefs().colwise().maxCoeff() <<".\n";
    gsInfo <<"Minimum deformation coef: "
           << deformation.patch(0).coefs().colwise().minCoeff() <<".\n";
    gsInfo <<"Area (undeformed) = "<<ev.integral(meas(G))<<"\tArea (undeformed) = "<<ev.integral(meas(defG))<<"\n";
    //! [Construct solution]

    //! [Evaluate solution]
    gsMatrix<> coords(2,1);
    if (testCase==1)
      coords<<0,0;
    else if (testCase==2)
      coords<<0.5,0;
    else if (testCase==3)
      coords<<0,0;
    else if (testCase==4)
      coords<<1,1;
    else
      coords<<0,0;

    gsMatrix<> res(3,1);
    mp.patch(0).eval_into(coords,res);
    real_t x=res.at(0);
    real_t y=res.at(1);
    real_t z=res.at(2);

    deformation.patch(0).eval_into(coords,res);
    real_t ux=res.at(0);
    real_t uy=res.at(1);
    real_t uz=res.at(2);

    gsInfo<<"Deformation on point ("<<x<<","<<y<<","<<z<<"):\n";
    gsInfo<<std::setprecision(20)<<"("<<ux<<","<<uy<<","<<uz<<")"<<"\n";
    //! [Evaluate solution]

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsField<> solField(mp, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "solution", 1000, true);
        // gsFileManager::open("solution.pvd");
    }
    //! [Export visualization in ParaView]
    return EXIT_SUCCESS;

}// end main

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int p, T L, T B)
{
  int q = p;
  int m = n;
  gsMultiPatch<T> mp = RectangularDomain(n, m, p, q, L, B);
  return mp;
}

template <class T>
gsMultiPatch<T> RectangularDomain(int n, int m, int p, int q, T L, T B)
{
  // -------------------------------------------------------------------------
  // --------------------------Make beam geometry-----------------------------
  // -------------------------------------------------------------------------
  int dim = 3; //physical dimension
  gsKnotVector<> kv0;
  kv0.initUniform(0,1,0,p+1,1);
  gsKnotVector<> kv1;
  kv1.initUniform(0,1,0,q+1,1);

  for(index_t i = 0; i< n; ++i)
      kv0.uniformRefine();
  for(index_t i = 0; i< m; ++i)
      kv1.uniformRefine();

  // Make basis
  gsTensorBSplineBasis<2,T> basis(kv0,kv1);

  // Initiate coefficient matrix
  gsMatrix<> coefs(basis.size(),dim);
  // Number of control points needed per component
  size_t len0 = basis.component(0).size();
  size_t len1 = basis.component(1).size();
  gsVector<> coefvec0(len0);
  // Uniformly distribute control points per component
  coefvec0.setLinSpaced(len0,0.0,L);
  gsVector<> coefvec1(basis.component(1).size());
  coefvec1.setLinSpaced(len1,0.0,B);

  // Z coordinate is zero
  coefs.col(2).setZero();

  // Define a matrix with ones
  gsVector<> temp(len0);
  temp.setOnes();
  for (index_t k = 0; k < len1; k++)
  {
    // First column contains x-coordinates (length)
    coefs.col(0).segment(k*len0,len0) = coefvec0;
    // Second column contains y-coordinates (width)
    coefs.col(1).segment(k*len0,len0) = temp*coefvec1.at(k);
  }
  // Create gsGeometry-derived object for the patch
  gsTensorBSpline<2,real_t> shape(basis,coefs);

  gsMultiPatch<T> mp;
  mp.addPatch(shape);
  mp.addAutoBoundaries();

  return mp;
}
