/** @file poisson2_example.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>

#  define MatExprType  auto

namespace gismo{
namespace expr{

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

    MatExprType eval(const index_t k) const
    {
        const index_t A = _u.cardinality()/_u.targetDim();
        res.resize(A*_u.targetDim(), cols()); // rows()*

        normal = _G.data().normal(k);// not normalized to unit length
        normal.normalize();
        bGrads = _u.data().values[1].col(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        // gsDebugVar(_G.data().values[0].col(k).transpose());

        for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all actives
            {
                // Jac(u) ~ Jac(G) with alternating signs ?..
                m_v.noalias() = (vecFun(d, bGrads.at(2*j  ) ).cross( cJac.col(1).template head<3>() )
                              - vecFun(d, bGrads.at(2*j+1) ).cross( cJac.col(0).template head<3>() )) / measure;

                // ---------------  First variation of the normal
                res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();

            }
        }
        return res;
    }

    index_t rows() const
    {
        return 1; //cols() * _u.data().values[1].rows() / _u.source().domainDim();
    }

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
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    static constexpr bool rowSpan() {return E::rowSpan(); }
    static constexpr bool colSpan() {return false;}

    void print(std::ostream &os) const { os << "var1("; _u.print(os); os <<")"; }
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

                        n_der.noalias() = (m_v - ( normal.dot(m_v) ) * normal);

                        m_uv.noalias() = ( vecFun(d, uGrads.at(2*j  ) ).cross( vecFun(c, vGrads.at(2*i+1) ) )
                                          +vecFun(c, vGrads.at(2*i  ) ).cross( vecFun(d, uGrads.at(2*j+1) ) ))
                                          / measure; //check

                        m_u_der.noalias() = (m_uv - ( normal.dot(m_v) ) * m_u);

                        // ---------------  Second variation of the normal
                        tmp = m_u_der - (m_u.dot(n_der) + normal.dot(m_u_der) ) * normal - (normal.dot(m_u) ) * n_der;

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
        return _u.source().domainDim() * ( _u.source().domainDim() + 1 ) / 2;
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
            [hess(v_1)*n_1 .., hess(v_1)*n_1 .., hess(v_1)*n_1 ..]. Here, the dots .. represent the active basis functions.
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


template<class E> EIGEN_STRONG_INLINE
var1_expr<E> var1(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return var1_expr<E>(u, G); }

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

}
}


using namespace gismo;

template <class T>
void evaluateFunction(gsExprEvaluator<T> ev, auto expression, gsVector<T> pt);
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

            C_constant = 4*lambda*mu/(lambda+2*mu);

            C(0,0) = C_constant*F0(0,0)*F0(0,0) + 2*mu*(2*F0(0,0)*F0(0,0));
            C(1,1) = C_constant*F0(1,1)*F0(1,1) + 2*mu*(2*F0(1,1)*F0(1,1));
            C(2,2) = C_constant*F0(0,1)*F0(0,1) + 2*mu*(F0(0,0)*F0(1,1) + F0(0,1)*F0(0,1));
            C(1,0) =
            C(0,1) = C_constant*F0(0,0)*F0(1,1) + 2*mu*(2*F0(0,1)*F0(0,1));
            C(2,0) =
            C(0,2) = C_constant*F0(0,0)*F0(0,1) + 2*mu*(2*F0(0,0)*F0(0,1));
            C(2,1) = C(1,2) = C_constant*F0(0,1)*F0(1,1) + 2*mu*(2*F0(0,1)*F0(1,1));

            //gsDebugVar(C);
        }
    }

    // piece(k) --> for patch k

};


// To Do:
// * struct for material model
// Input is parametric coordinates of the surface \a mp
template <class T>
class gsMaterialMatrixIncompressible : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
  //
protected:
    const gsFunctionSet<T> * _mp;
    const gsFunction<T> * _par1;
    // const gsFunction<T> * _par2;
    const gsFunctionSet<T> * _mp_def;
    mutable gsMapData<T> _tmp;
    mutable gsMapData<T> _tmp_def;
    mutable gsMatrix<real_t,3,3> F0, jacGdef, jacGori;
    mutable gsMatrix<T> par1mat,par2mat;
    mutable real_t mu, J0;

public:
    /// Shared pointer for gsMaterialMatrixIncompressible
    typedef memory::shared_ptr< gsMaterialMatrixIncompressible > Ptr;

    /// Unique pointer for gsMaterialMatrixIncompressible
    typedef memory::unique_ptr< gsMaterialMatrixIncompressible > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsMaterialMatrixIncompressible( const gsFunctionSet<T> & mp,
                                    const gsFunction<T> & par1,
                                    // const gsFunction<T> & par2,
                                    const gsFunctionSet<T> & mp_def) : // deformed multipatch
    _mp(&mp),
    _par1(&par1),
    // _par2(&par2),
    _mp_def(&mp_def),
    _mm_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
        _tmp_def.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
    }

    ~gsMaterialMatrixIncompressible() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrixIncompressible)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 9;}

    mutable gsMaterialMatrixIncompressible<T> * _mm_piece; // todo: improve the way pieces are accessed

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _mm_piece;
        // _mm_piece = new gsMaterialMatrixIncompressible(_mp->piece(k), *_par1, *_par2, _mp_def->piece(k) );
        _mm_piece = new gsMaterialMatrixIncompressible(_mp->piece(k), *_par1, _mp_def->piece(k) );
        return *_mm_piece;
    }

    // Input is parametric coordinates of the surface \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // NOTE 1: if the input \a u is considered to be in physical coordinates
        // then we first need to invert the points to parameter space
        // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
        // otherwise we just use the input paramteric points
        _tmp.points = u;
        _tmp_def.points = u;

        real_t z = 0; // parametric coordinate in thickness direction

        static_cast<const gsFunction<T>*>(_mp)->computeMap(_tmp);
        static_cast<const gsFunction<T>*>(_mp_def)->computeMap(_tmp_def);

        // NOTE 2: in the case that parametric value is needed it suffices
        // to evaluate Youngs modulus and Poisson's ratio at
        // \a u instead of _tmp.values[0].
        _par1->eval_into(_tmp.values[0], par1mat);
        // _par2->eval_into(_tmp.values[0], par2mat);

        result.resize( targetDim() , u.cols() );
        for( index_t i=0; i< u.cols(); ++i )
        {
            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(i,3,3);

            // Compute metric tensor gab deformed: jacGdef
            jacGdef.leftCols(2) = _tmp_def.jacobian(i);                 // The in-plane surface tangents and...
            jacGdef.col(2)      = _tmp_def.normal(i).normalized();      // The normal vector form the covariant basis, which is stored in the cols of jacGdef
            jacGdef = jacGdef.inverse();                                // Now we have the inverse of this matrix
            jacGdef = jacGdef * jacGdef.transpose(); //3x3              // And this is the contravariant metric tensor

            // Compute metric tensor gab undeformed (original): jacGori
            jacGori.leftCols(2) = _tmp.jacobian(i);                     // The in-plane surface tangents and...
            jacGori.col(2)      = _tmp.normal(i).normalized();          // The normal vector form the covariant basis, which is stored in the cols of jacGori
            jacGori = jacGori.inverse();                                // Now we have the inverse of this matrix
            jacGori = jacGori * jacGori.transpose(); //3x3              // And this is the contravariant metric tensor

            // Evaluate material properties on the quadrature point
            mu = par1mat(0,i);
            J0 = math::sqrt( jacGdef.determinant() / jacGori.determinant() );
            J0 = math::pow( J0, -2 );

            /*
                C =     C1111,  C1122,  C1112
                        symm,   C2222,  C2212
                        symm,   symm,   C1212
            */

            C(0,0) = 2*F0(0,0)*F0(0,0) + F0(0,0)*F0(0,0) + F0(0,0)*F0(0,0);                 // C1111
            C(1,0) =
            C(0,1) = 2*F0(0,0)*F0(0,0) + F0(0,0)*F0(0,0) + F0(0,0)*F0(0,0);                 // C1122
            C(2,0) =
            C(0,2) = 2*F0(0,0)*F0(0,1) + F0(0,0)*F0(0,1) + F0(0,1)*F0(0,0);                 // C1112
            C(1,1) = 2*F0(1,1)*F0(1,1) + F0(1,1)*F0(1,1) + F0(1,1)*F0(1,1);                 // C2222
            C(2,1) =
            C(1,2) = 2*F0(1,1)*F0(0,1) + F0(1,0)*F0(1,1) + F0(1,1)*F0(1,0);                 // C2212
            C(2,2) = 2*F0(0,1)*F0(0,1) + F0(0,0)*F0(1,1) + F0(0,1)*F0(1,0);                 // C1212

            C *= mu * J0;

            //gsDebugVar(C);



            // for completeness, the computation of variable S^\alpha\beta is given here, commented. It should be implemented later in a separate class or function
             // TRANSFORM JACOBIANS!
            // Sab(0,0) = mu * (jacGori(0,0) - math::pow(J0,-2.) * jacGdef(0,0) );
            // Sab(0,1) =
            // Sab(1,0) = mu * (jacGori(1,0) - math::pow(J0,-2.) * jacGdef(1,0) ); // CHECK SYMMETRIES
            // Sab(1,1) = mu * (jacGori(1,1) - math::pow(J0,-2.) * jacGdef(1,1) );

        }
    }

};

// To Do:
// * struct for material model
// Input is parametric coordinates of the surface \a mp
template <class T>
class gsMaterialMatrixCompressible : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
  //
protected:
    const gsFunctionSet<T> * _mp;
    const gsFunction<T> * _par1;
    const gsFunction<T> * _par2;
    const gsFunctionSet<T> * _mp_def;
    mutable gsMapData<T> _tmp;
    mutable gsMapData<T> _tmp_def;
    mutable gsMatrix<real_t,3,3> F0, jacGdef, jacGori;
    mutable gsMatrix<T> par1mat,par2mat;
    mutable real_t mu, K, J0, J;

public:
    /// Shared pointer for gsMaterialMatrixCompressible
    typedef memory::shared_ptr< gsMaterialMatrixCompressible > Ptr;

    /// Unique pointer for gsMaterialMatrixCompressible
    typedef memory::unique_ptr< gsMaterialMatrixCompressible > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsMaterialMatrixCompressible( const gsFunctionSet<T> & mp,
                                    const gsFunction<T> & par1,
                                    const gsFunction<T> & par2,
                                    const gsFunctionSet<T> & mp_def) : // deformed multipatch
    _mp(&mp),
    _par1(&par1),
    _par2(&par2),
    _mp_def(&mp_def),
    _mm_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
        _tmp_def.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
    }

    ~gsMaterialMatrixCompressible() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrixCompressible)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 9;}

    mutable gsMaterialMatrixCompressible<T> * _mm_piece; // todo: improve the way pieces are accessed

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _mm_piece;
        _mm_piece = new gsMaterialMatrixCompressible(_mp->piece(k), *_par1, *_par2, _mp_def->piece(k) );
        // _mm_piece = new gsMaterialMatrixCompressible(_mp->piece(k), *_par1, _mp_def.piece(k) );
        return *_mm_piece;
    }



//class .. matMatrix_z
// should contain eval_into(thickness variable)


    // void computeThickness()
    // {

    // }

    // Input is parametric coordinates of the surface \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // NOTE 1: if the input \a u is considered to be in physical coordinates
        // then we first need to invert the points to parameter space
        // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
        // otherwise we just use the input paramteric points
        _tmp.points = u;
        _tmp_def.points = u;

        static_cast<const gsFunction<T>*>(_mp)->computeMap(_tmp);
        static_cast<const gsFunction<T>*>(_mp_def)->computeMap(_tmp_def);

        // NOTE 2: in the case that parametric value is needed it suffices
        // to evaluate Youngs modulus and Poisson's ratio at
        // \a u instead of _tmp.values[0].
        _par1->eval_into(_tmp.values[0], par1mat);
        _par2->eval_into(_tmp.values[0], par2mat);

        result.resize( targetDim() , u.cols() );
        for( index_t i=0; i< u.cols(); ++i )
        {
            // Material parameters
            mu = par1mat(0,i);
            K = par1mat(0,i);

            // Define objects
            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(i,3,3);
            gsMatrix<T,3,3> c, cinv;
            T S33, C3333, dc33, traceC;

            // Compute metric tensor gab deformed: jacGdef
            jacGdef.leftCols(2) = _tmp_def.jacobian(i);                 // The in-plane surface tangents and...
            jacGdef.col(2)      = _tmp_def.normal(i).normalized();      // The normal vector form the covariant basis, which is stored in the cols of jacGdef
            jacGdef = jacGdef * jacGdef.transpose(); //3x3              // And this is the covariant metric tensor

            // Compute metric tensor gab undeformed (original): jacGori
            jacGori.leftCols(2) = _tmp.jacobian(i);                     // The in-plane surface tangents and...
            jacGori.col(2)      = _tmp.normal(i).normalized();          // The normal vector form the covariant basis, which is stored in the cols of jacGori
            jacGori = jacGori * jacGori.transpose(); //3x3              // And this is the covariant metric tensor

            // Initialize c
            c.setZero();
            c.block(0,0,2,2) = jacGdef.block(0,0,2,2);
            c(2,2) = 1.0; // c33
            cinv = c.inverse();
            // note: can also just do c = jacGdef because the normal has length one and hence c(2,2) is 1. CHECK!

            J0 = math::sqrt( jacGdef.determinant() / jacGori.determinant() );
            J = J0 * math::sqrt( c(2,2) );

            index_t imax = 20;
            T tol = 1e-6;
            S33 = 0.0;
            C3333 = 1.0;

            // Define lambda function for C
            std::function<T (index_t i, index_t j, index_t k, index_t l)> Cijkl;
            Cijkl = [=](index_t i, index_t j, index_t k, index_t l)
            {
                T res = 1.0 / 9.0 * mu * math::pow( J , -2.0/3.0 ) * ( traceC * ( 2*cinv(i,j)*cinv(k,l) + 3*cinv(i,k)*cinv(j,l) + 3*cinv(i,l)*cinv(j,k) )
                                - 6*jacGori(i,j)*cinv(k,l) + cinv(i,j)*jacGori(k,l) ) + K * ( J*J*cinv(i,j)*cinv(k,l) - 0.5*(J*J-1)*( cinv(i,k)*cinv(j,l) + cinv(i,l)*cinv(j,k) ) );
                return res;
            };

            for (index_t i = 0; i < imax; i++)
            {
                dc33 = -2. * S33 / C3333;
                c(2,2) += dc33;
                cinv(2,2) = 1.0/c(2,2);

                traceC = c.trace();
                J = J0 * math::sqrt( c(2,2) );

                S33     = mu * math::pow( J , -2.0/3.0 ) * ( jacGori(2,2) - 1.0/3.0 * traceC * cinv(2,2) ) + 0.5 * K * ( J*J - 1 ) * cinv(2,2);
                C3333   = Cijkl(2,2,2,2);

                if (S33 < tol)
                {
                    gsInfo<<"Converged, S33 = "<<S33<<" and tolerance = "<<tol<<"\n";

/*
                    S.at(0) = mu * math::pow( J , -2.0/3.0 ) * ( jacGori(0,0) - 1.0/3.0 * traceC * cinv(0,0) ) + 0.5 * K * ( J*J - 1 ) * cinv(0,0); // S11
                    S.at(1) = mu * math::pow( J , -2.0/3.0 ) * ( jacGori(1,1) - 1.0/3.0 * traceC * cinv(1,1) ) + 0.5 * K * ( J*J - 1 ) * cinv(1,1); // S22
                    S.at(2) = mu * math::pow( J , -2.0/3.0 ) * ( jacGori(0,1) - 1.0/3.0 * traceC * cinv(0,1) ) + 0.5 * K * ( J*J - 1 ) * cinv(0,1); // S12
                    // S(i,j) = mu * math::pow( J , -2.0/3.0 ) * ( jacGori(i,j) - 1.0/3.0 * traceC * cinv(i,j) ) + 0.5 * K * ( J*J - 1 ) * cinv(i,j);
*/
                    /*
                        C =     C1111,  C1122,  C1112
                                symm,   C2222,  C2212
                                symm,   symm,   C1212
                        Here, Cabcd = Cijkl - Cab33*C33cd / C3333;
                        a,b,c,d = 1,2; i,j,k,l = 1...3;
                    */
                    C(0,0) = Cijkl(0,0,0,0) - ( Cijkl(0,0,2,2) * Cijkl(2,2,0,0) ) / (Cijkl(2,2,2,2)); // C1111
                    C(0,1) =
                    C(1,0) = Cijkl(0,0,1,1) - ( Cijkl(0,0,2,2) * Cijkl(2,2,1,1) ) / (Cijkl(2,2,2,2)); // C1122
                    C(0,2) =
                    C(2,0) = Cijkl(0,0,0,1) - ( Cijkl(0,0,2,2) * Cijkl(2,2,0,1) ) / (Cijkl(2,2,2,2)); // C1112
                    C(1,1) = Cijkl(1,1,1,1) - ( Cijkl(1,1,2,2) * Cijkl(2,2,1,1) ) / (Cijkl(2,2,2,2)); // C2222
                    C(1,2) =
                    C(2,1) = Cijkl(1,1,0,1) - ( Cijkl(1,1,2,2) * Cijkl(2,2,0,1) ) / (Cijkl(2,2,2,2)); // C2212
                    C(2,2) = Cijkl(0,1,0,1) - ( Cijkl(0,1,2,2) * Cijkl(2,2,0,1) ) / (Cijkl(2,2,2,2)); // C1212
                }
                else if (i == imax - 1)
                {
                    gsInfo<<"Error: Method did not converge";
                    std::terminate();
                }
            }
        }
    }
};



//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    index_t testCase = 1;
    bool nonlinear = false;
    std::string fn("pde/poisson2d_bvp.xml");

    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    real_t thickness = 1.0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "t", "testCase", "Test case to run: 1 = unit square; 2 = Scordelis Lo Roof",  testCase );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("nl", "Solve nonlinear problem", nonlinear);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]
    gsMultiPatch<> mp;
    gsMultiPatch<> mp_def;
    if (testCase==1)
    {
        // Unit square
        mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
        mp.addAutoBoundaries();
        mp.embed(3);
        E_modulus = 1.0;
        thickness = 0.5;

    }
    else if (testCase == 2)
    {
        thickness = 0.125;
        E_modulus = 4.32E8;
        fn = "../extensions/unsupported/filedata/scordelis_lo_roof.xml";
        gsReadFile<>(fn, mp);
    }
    else if (testCase==9)
    {
        gsFileData<> fd(fn);
        gsInfo << "Loaded file "<< fd.lastPath() <<"\n";
        // Annulus
        fd.getId(0, mp); // id=0: Multipatch domain
        mp.embed(3);
        E_modulus = 1.0;
        thickness = 0.5;
    }
    //! [Read input file]

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    mp_def = mp;

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
    gsInfo << dbasis.basis(0)<<"\n";

    gsBoundaryConditions<> bc;
    gsVector<> tmp(3);
    tmp << 0, 0, 0;
    if (testCase == 1)
    {
        for (index_t i=0; i!=3; ++i)
        {
            bc.addCondition(boundary::north, condition_type::dirichlet, 0, i ); // unknown 0 - x
            bc.addCondition(boundary::east, condition_type::dirichlet, 0, i ); // unknown 1 - y
            bc.addCondition(boundary::south, condition_type::dirichlet, 0, i ); // unknown 2 - z
            bc.addCondition(boundary::west, condition_type::dirichlet, 0, i ); // unknown 2 - z
        }
        tmp << 0,0,-1;
    }
    else if (testCase == 2)
    {
        // Diaphragm conditions
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 1 ); // unknown 1 - y
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 2 ); // unknown 2 - z

        // ORIGINAL
        // bc.addCornerValue(boundary::southwest, 0.0, 0, 0); // (corner,value, patch, unknown)

        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 1 ); // unknown 1 - y
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 2 ); // unknown 2 - z

        // NOT ORIGINAL
        bc.addCondition(boundary::west, condition_type::dirichlet, 0, 0 ); // unknown 1 - y
        bc.addCondition(boundary::east, condition_type::dirichlet, 0, 0 ); // unknown 1 - y

        // Surface forces
        tmp << 0, 0, -90;
    }
    else if (testCase == 9)
    {
        for (index_t i=0; i!=3; ++i)
        {
            // patch 0
            bc.addCondition(0,boundary::north, condition_type::dirichlet, 0, i ); // unknown 0 - x
            // bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, i ); // unknown 1 - y
            bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, i ); // unknown 2 - z
            bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, i ); // unknown 2 - z
            // patch 1
            bc.addCondition(1,boundary::north, condition_type::dirichlet, 0, i ); // unknown 0 - x
            bc.addCondition(1,boundary::east, condition_type::dirichlet, 0, i ); // unknown 1 - y
            bc.addCondition(1,boundary::south, condition_type::dirichlet, 0, i ); // unknown 2 - z
            // bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, i ); // unknown 2 - z
        }
        tmp << 0,0,-1;
    }
    //! [Refinement]

    //! [Problem setup]
    gsExprAssembler<> A(1,1);

    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
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
    u.setInterfaceCont(0); // todo: 1 (smooth basis)
    u.addBc( bc.get("Dirichlet") ); // (!) must be called only once

    // Solution vector and solution variable
    gsMatrix<> solVector, vec3(3,1); vec3.setZero();
    solution u_sol = A.getSolution(u, solVector);

    // gsFunctionExpr<> materialMat("1","0","0","0","1","0","0","0","1",3);
    // variable mm = A.getCoeff(materialMat, G);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsMaterialMatrix materialMat(mp, E, nu);
    variable mm = A.getCoeff(materialMat);


    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",3);
    variable m2 = A.getCoeff(mult2t, G);

    gsFunctionExpr<> t(std::to_string(thickness), 3);
    variable tt = A.getCoeff(t,G);
    // TEMPORARILY!!
    // real_t tt = thickness;

    // gsFunctionExpr<> force("0","0","1", 3);
    gsConstantFunction<> force(tmp,3);
    variable ff = A.getCoeff(force, G);

    gsSparseSolver<>::CGDiagonal solver;

    //! [Problem setup]

    //! [Solver loop]

    // Set Dirichlet values
    //A.options().setInt("DirichletValues", dirichlet::homogeneous);

    // Initialize the system
    A.initSystem();

    gsInfo<< A.numDofs() <<"\n"<<std::flush;

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
    auto E_m = 0.5 * ( flat(jac(defG).tr()*jac(defG)) - flat(jac(G).tr()* jac(G)) ) ; //[checked]
    auto E_m_der = flat( jac(defG).tr() * jac(u) ) ; //[checked]
    auto E_m_der2 = flatdot( jac(u),jac(u).tr(), E_m * reshape(mm,3,3) ); //[checked]

    auto E_f = ( deriv2(G,sn(G).normalized().tr()) - deriv2(defG,sn(defG).normalized().tr()) ) * reshape(m2,3,3) ; //[checked]
    auto E_f_der = ( deriv2(u,sn(defG).normalized().tr() ) + deriv2(defG,var1(u,defG) ) ) * reshape(m2,3,3); //[checked]
    auto E_f_der2 = flatdot2( deriv2(u), var1(u,defG).tr(), E_f * reshape(mm,3,3)  ).symmetrize() + var2(u,u,defG,E_f * reshape(mm,3,3) );
    // NOTE: var1(u,G) in E_F_der2 should be var1(u,defG)


    gsVector<> pt(2); pt.setConstant(0.25);
    // evaluateFunction(ev, u * ff * meas(G), pt); // evaluates an expression on a point

    // ! [Solve linear problem]

    // assemble system
    A.assemble(
        (
        (tt.val()) * (E_m_der * reshape(mm,3,3) * E_m_der.tr())
        +
        (tt.val() * tt.val() * tt.val())/3.0 * (E_f_der * reshape(mm,3,3) * E_f_der.tr())
        ) * meas(G)
        ,u * ff * meas(G)
        );

    // solve system
    solver.compute( A.matrix() );
    solVector = solver.solve(A.rhs());

    // update deformed patch
    gsMatrix<> cc;
    for ( size_t k =0; k!=mp_def.nPatches(); ++k) // Deform the geometry
    {
        // extract deformed geometry
        u_sol.extract(cc, k);
        mp_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
    }

    /*Something with Dirichlet homogenization*/

    // ! [Solve linear problem]

    // ! [Solve nonlinear problem]

    real_t residual = A.rhs().norm();
    if (nonlinear)
    {
        index_t itMax = 10;
        real_t tol = 1e-8;
        for (index_t it = 0; it != itMax; ++it)
        {
            A.initSystem();
            // assemble system
            A.assemble(
                (
                (tt.val()) * (E_m_der * reshape(mm,3,3) * E_m_der.tr() + E_m_der2)
                +
                (tt.val() * tt.val() * tt.val())/3.0 * (E_f_der * reshape(mm,3,3) * E_f_der.tr() -  E_f_der2)
                ) * meas(G)
                , u * ff * meas(G)
                -
                (
                 (
                    (tt.val()) *(E_m * reshape(mm,3,3) * E_m_der.tr()) -
                    (tt.val() * tt.val() * tt.val())/3.0 * (E_f * reshape(mm,3,3) * E_f_der.tr())
                 ) * meas(G)
                ).tr()
                );

            // A.assemble(tt.val() * tt.val() * tt.val() / 3.0 * E_f_der2);
            // solve system
            solver.compute( A.matrix() );
            solVector = solver.solve(A.rhs()); // this is the UPDATE
            residual = A.rhs().norm();

            gsInfo<<"Iteration: "<< it
                   <<", residue: "<< residual
                   <<", update norm: "<<solVector.norm()
                   <<"\n";

            // update deformed patch
            gsMatrix<> cc;
            for ( size_t k =0; k!=mp_def.nPatches(); ++k) // Deform the geometry
            {
                // extract deformed geometry
                u_sol.extract(cc, k);
                mp_def.patch(k).coefs() += cc;  // defG points to mp_def, therefore updated
            }


            if (residual < tol)
                break;
        }
    }


    // ! [Solve nonlinear problem]

    // For Neumann (same for Dirichlet/Nitche) conditions
    // variable g_N = A.getBdrFunction();
    // A.assembleRhsBc(u * g_N.val() * nv(G).norm(), bc.neumannSides() );

    // Penalize the matrix? (we need values for the DoFs to be enforced..
    // function/call:  penalize_matrix(DoF_indices, DoF_values)
    // otherwise: should we tag the DoFs inside "u" ?

    // gsInfo<<"RHS rows = "<<A.rhs().rows()<<"\n";
    // gsInfo<<"RHS cols = "<<A.rhs().cols()<<"\n";
    // gsInfo<<"MAT rows = "<<A.matrix().rows()<<"\n";
    // gsInfo<<"MAT cols = "<<A.matrix().cols()<<"\n";

    // gsInfo<< A.rhs().transpose() <<"\n";
    // gsInfo<< A.matrix().toDense()<<"\n";

    // ADD BOUNDARY CONDITIONS! (clamped will be tricky..............)


    //! [Export visualization in ParaView]
    if (plot)
    {
        gsMultiPatch<> deformation = mp_def;
        for (index_t k = 0; k != mp_def.nPatches(); ++k)
            deformation.patch(0).coefs() -= mp.patch(0).coefs();

        gsField<> solField(mp, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "solution", 1000, true);

        // ev.options().setSwitch("plot.elements", true);
        // ev.writeParaview( u_sol   , G, "solution");

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

/*
    to do:
    =  make function for construction of the solution given the space and the mp
*/



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