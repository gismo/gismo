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
                res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
            }
        }

        // gsDebugVar(res.rows());
        // gsDebugVar(res.cols());
        // gsDebugVar(res);
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

    static constexpr bool rowSpan() {return E::rowSpan(); }
    static constexpr bool colSpan() {return false;}

    void print(std::ostream &os) const { os << "var1("; _u.print(os); os <<")"; }
};

// Comments for var2:
// - TODO: dimensionm indep. later on
// - TODO: how to structure this matrix
template<class E1, class E2>
class var2_expr : public _expr<var2_expr<E1,E2> >
{
public:
    typedef typename E1::Scalar Scalar;

private:

    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _G;

public:
    enum{ Space = E1::Space };

    var2_expr( const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & G) : _u(u),_v(v), _G(G) { }

    mutable gsMatrix<Scalar> res;

    mutable gsMatrix<Scalar> uGrads, vGrads, cJac, cDer2;
    mutable gsVector<Scalar,3> m_u, m_v, normal, m_uv, m_u_der, n_der, n_der2, tmp;
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
        res.resize(rows(), cols());

        normal = _G.data().normals.col(k);
        uGrads = _u.data().values[1].col(k);
        vGrads = _v.data().values[1].col(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        cDer2 = _G.data().values[2].reshapeCol(k, _G.data().dim.second, _G.data().dim.second);

        const index_t cardU = _u.data().values[0].rows(); // cardinality of u
        const index_t cardV = _v.data().values[0].rows(); // cardinality of v
        gsDebugVar(cardU);
        const index_t numHess = cDer2.rows();
        const Scalar measure =  _G.data().measures.at(k);


        for (index_t j = 0; j!= cardU; ++j) // for all basis functions u (1)
        {
            for (index_t i = 0; i!= cardV; ++i) // for all basis functions v (1)
            {
                for (index_t d = 0; d!= _u.dim(); ++d) // for all basis functions u (2)
                {
                    m_u.noalias() = vecFun(d, uGrads.at(2*j  ) ).cross( cJac.col(1).template head<3>() )
                                  - vecFun(d, uGrads.at(2*j+1) ).cross( cJac.col(0).template head<3>() ) / measure;
                    const short_t s = d*cardU;
                    for (index_t c = 0; c!= _v.dim(); ++c) // for all basis functions v (2)
                    {
                        const short_t r = c*cardV;
                        m_v.noalias() = vecFun(c, vGrads.at(2*i  ) ).cross( cJac.col(1).template head<3>() )
                                      - vecFun(c, vGrads.at(2*i+1) ).cross( cJac.col(0).template head<3>() ) / measure;

                        n_der = (m_v - ( normal.dot(m_v) ) * normal);

                        m_uv.noalias() = vecFun(d, uGrads.at(2*j  ) ).cross( vecFun(c, vGrads.at(2*i+1) ) )
                                       - vecFun(d, uGrads.at(2*j+1) ).cross( vecFun(c, vGrads.at(2*i  ) ) ) / measure;

                        m_u_der.noalias() = (m_uv - ( normal.dot(m_v) ) * m_u);

                        // ---------------  Second variation of the normal
                        tmp = m_u_der - (m_u.dot(n_der) + normal.dot(m_u_der) ) * normal - (normal.dot(m_u) ) * n_der;

                        // NOTE: This all happens in Colspace. Is that right?
                        // NOTE: We multiply with the column of c. This means that we multiply with [d11 c1, d22 c1, d12 c1]. Should't it be [d11 c1, d11 c2, d11 c3]. We simplyneed to transppse c.
                        for (index_t l=0; l != numHess; l++) // per hessian entry of c
                        {
                            res(s + j, r + i + l*_u.dim()*cardU ) = tmp.dot(cDer2.col(l));
                        }
                    }
                }
            }
        }

        return res;
    }

    index_t rows() const
    {
        return _u.dim() * _u.data().values[0].rows();
    }

    index_t cols() const
    {
        return  _u.source().domainDim() * (_u.source().domainDim()+1) / 2 * _v.dim() * _v.data().values[0].rows(); // hessian of c has dimension d(d+1)/2 and the columnspace of v has dimension v.dim()*numActive
    }

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_2ND_DER | NEED_MEASURE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_OUTER_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    const gsFeSpace<Scalar> & rowVar() const { true; }
    const gsFeSpace<Scalar> & colVar() const { true; }

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

    const gsMatrix<Scalar> & eval(const index_t k) const { return eval_impl(_u,k); }

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
    static constexpr bool colSpan() {return E2::rowSpan();}

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
        tmp =_u.data().values[2].reshapeCol(k, rows(), _u.data().dim.second).transpose();
        vEv = _v.eval(k);

        // gsDebugVar(tmp.rows());
        // gsDebugVar(tmp.cols());
        // gsDebugVar(vEv.rows());
        // gsDebugVar(vEv.cols());

        res = vEv * tmp;

        //..
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
        const index_t numAct = u.data().values[0].rows();
        res.resize(rows()*numAct, cols() );
        tmp.transpose() =_u.data().values[2].reshapeCol(k, rows() , numAct );
        vEv = _v.eval(k);

        for (index_t i = 0; i!=_u.dim(); i++)
            res.block(i*numAct, 0, numAct, cols() ) = tmp * vEv.at(i);

        return res;
    }

};

template<class E>
class deriv2_expr : public _expr<deriv2_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:

    typename E::Nested_t _u;

public:
    enum{ Space = E::Space };

    deriv2_expr(const E & u) : _u(u) { }

    mutable gsMatrix<Scalar> res;

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        //gsDebugVar(_u.targetDim());
        //const index_t c = _u.cardinality() / _u.targetDim();
        res = _u.data().values[2].blockDiag(_u.targetDim());
        return res;
    }

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
    const gsFeSpace<Scalar> & colVar() const { return _u.rowVar(); }

    static constexpr bool rowSpan() {return E::rowSpan(); }
    static constexpr bool colSpan() {return E::rowSpan();}

    void print(std::ostream &os) const { os << "deriv2("; _u.print(os); os <<")"; }
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
   Takes
   A: rowspace
   B: rowspace
   C: vector

   and computes

   A* : B

   where A* = A "hadamard" C
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
    mutable gsMatrix<Scalar> tmpA, tmpB, tmpC, res;

public:

    flatdot_expr(_expr<E1> const& A, _expr<E2> const& B, _expr<E3> const& C) : _A(A),_B(B),_C(C)
    {
        //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t Ac = _A.cols();
        //const index_t Ar = _A.rows();
        const index_t An = _A.cardinality();
        const index_t Bc = _B.cols();
        //const index_t Br = _B.rows();
        const index_t Bn = _B.cardinality();

        MatExprType eA = _A.eval(k);
        const MatExprType eB = _B.eval(k);

        // multiply by C
        const MatExprType eC = _C.eval(k);

        // gsDebugVar(Ac);
        // gsDebugVar(Ar);
        // gsDebugVar(An);
        // gsDebugVar(Bc);
        // gsDebugVar(Br);
        // gsDebugVar(Bn);
        // gsDebugVar(eC);

        for (index_t i = 0; i!=An; ++i)
        {
            // gsDebugVar(eA(0,2*i  ));
            // gsDebugVar(eA(0,2*i+1));
            // gsDebugVar(eC);

            eA(0,2*i  ) *= eC(0);
            eA(1,2*i+1) *= eC(1);
            eA(1,2*i)   *= eC(2);
            eA(0,2*i+1) *= eC(2);
        }

        res.resize(An, Bn);
        for (index_t i = 0; i!=An; ++i)
            for (index_t j = 0; j!=Bn; ++j)
            {
                res(i,j) = ( eA.middleCols(i*Ac,Ac) * eB.middleCols(j*Bc,Bc) ).sum();
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
    static bool colSpan() {return E2::colSpan();}

    void print(std::ostream &os) const { os << "flatdot("; _A.print(os);_B.print(os);_C.print(os); os<<")"; }
};


template<class E> EIGEN_STRONG_INLINE
var1_expr<E> var1(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return var1_expr<E>(u, G); }

template<class E1, class E2> EIGEN_STRONG_INLINE
var2_expr<E1,E2> var2(const E1 & u, const E2 & v, const gsGeometryMap<typename E1::Scalar> & G)
{ return var2_expr<E1,E2>(u,v, G); }

// template<class E1, class E2> EIGEN_STRONG_INLINE
// hessdot_expr<E1,E2> hessdot(const E1 & u, const E2 & v) { return hessdot_expr<E1,E2>(u, v); }

template<class E> EIGEN_STRONG_INLINE
deriv2_expr<E> deriv2(const E & u) { return deriv2_expr<E>(u); }

template<class E1, class E2> EIGEN_STRONG_INLINE
deriv2dot_expr<E1, E2> deriv2(const E1 & u, const E2 & v) { return deriv2dot_expr<E1, E2>(u,v); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
flatdot_expr<E1,E2,E3> flatdot(const E1 & u, const E2 & v, const E3 & w)
{ return flatdot_expr<E1,E2,E3>(u, v, w); }

}
}


using namespace gismo;

//! [Include namespace]



int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 0;
    bool last = false;
    std::string fn("pde/poisson2d_bvp.xml");

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addString( "f", "file", "Input XML file", fn );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> mp;

    // Annulus
    fd.getId(0, mp); // id=0: Multipatch domain

    // Unit square
    // mp.addPatch( gsNurbsCreator<>::BSplineSquare(numElevate+1) ); // degree



    //! [Read input file]

    //! [Refinement]
    gsMultiBasis<> dbasis(mp);

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( dbasis.maxCwiseDegree() + numElevate);

    // h-refine each basis
    if (last)
    {
        for (int r =0; r < numRefine-1; ++r)
            dbasis.uniformRefine();
        numRefine = 0;
    }

    mp.embed(3);
    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< dbasis.minCwiseDegree() <<"\n";
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
    geometryMap G = A.getMap(mp);
    geometryMap defG = G;           // defG ??????

    // Set the discretization space
    space u = A.getSpace(dbasis, 3);
    u.setInterfaceCont(0);

    // Solution vector and solution variable
    gsMatrix<> solVector, vec3(3,1); vec3.setZero();
    solution u_sol = A.getSolution(u, solVector);

    gsFunctionExpr<> materialMat("1","0","0","0","1","0","0","0","1",3);
    variable mm = A.getCoeff(materialMat, G);


    gsFunctionExpr<> mult2t("1","0","0","0","1","0","0","0","2",3);
    variable m2 = A.getCoeff(mult2t, G);

    gsFunctionExpr<> thickness("1", 2);
    variable tt = A.getCoeff(thickness);

    gsFunctionExpr<> force("0","0""1", 3);
    variable ff = A.getCoeff(force, G);

    gsSparseSolver<>::CGDiagonal solver;

    //! [Problem setup]

    //! [Solver loop]
    gsVector<> l2err(numRefine+1), h1err(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=got_error)\n"
        "\nDoFs: ";
    for (int r=0; r<=numRefine; ++r)
    {
        dbasis.uniformRefine();

        // Initialize the system
        A.initSystem();

        gsInfo<< A.numDofs() <<"\n"<<std::flush;

        // Compute the system matrix and right-hand side

        auto E_m = 0.5 * ( flat(jac(defG).tr()*jac(defG)) - flat(jac(G).tr()* jac(G)) ); //[checked]
        auto E_m2 = 0.5 * ( jac(defG).tr()*jac(defG) - jac(G).tr()* jac(G) );

        auto E_m_der = flat( jac(u).tr() * jac(defG) ) ; //[checked]
        // auto E_m_der2 = flat( jac(u).tr() * jac(u) );
        auto E_m_der2 = flatdot( jac(u).tr(),jac(u), E_m.tr() * reshape(mm,3,3) );

        auto E_f = ( deriv2(defG,sn(defG).normalized().tr()) - deriv2(G,sn(G).normalized().tr()) ) * reshape(m2,3,3) ;//[checked]
        auto E_f_der = ( deriv2(u,sn(defG).normalized().tr() ) + deriv2(defG,var1(u,G) ) ) * reshape(m2,3,3);
        //auto E_f_der2 = reshape(m2,3,3) * ( 2 * deriv2(u) * var1(u,defG).tr() + var2(u,u,defG) );

        auto E_f_der2 = flatdot(
            deriv2(u).tr(),
            var1(u,defG).tr(),
            E_f.tr() * reshape(mm,3,3) );
        //    + flatdot( reshape(m2,3,3) * deriv2(u).tr() , var1(u,defG).tr(), E_f.tr() * reshape(mm,3,3) ).tr()

            // + var2*Ef
            ;

        //A.assemble( tt * ( E_m.tr() * reshape(mm,3,3) * E_m_der2 + E_m_der.tr() * reshape(mm,3,3) * E_m_der ) +
        //             (tt*tt*tt/3.0) *( E_f * mm * E_f_der2 + E_f_der * mm * E_f_der ), // Matrix

        //             tt * ( E_m * mm * E_m_der ) + (tt*tt*tt/3.0) *( E_f * mm * E_f_der )  // Rhs
        //             - u * ff.tr()
        //     );

        gsVector<> pt(2); pt.setConstant(0.5);
        gsMatrix<> evresult = ev.eval(
            // E_m              // works; output 3 x 1
            // E_m_der          // works; output 3 x 27
            // E_m_der2         // works; output 27 x 27
//             E_f              // works: output 3 x 1
//             E_f_der          // works output  27 x 3
             E_f_der2         // does not work;

//            deriv2(u,var1(u,G))
//            deriv2(u)

//
//            deriv2(G,var1(u,G))
//            deriv2(G, var1(u,G) )
            // var2(u,u,defG)

            //var1(u,G) * deriv2(defG) //.tr()
            //deriv2(u) * sn(G)

            // deriv2(u)
            // var1(u,G)
            //deriv2(u) //.tr() //* sn(defG).normalized()
            //deriv2(G) //.tr() //* sn(defG).normalized()
                                       ,  pt
            );
        gsInfo << "Eval:\n"<< evresult;
        gsInfo << "\nEnd ("<< evresult.rows()<< " x "<<evresult.cols()<<")\n";

        return 0;

        A.assemble(
            //tt * (

            E_m_der.tr()
              * reshape(mm,3,3)
                * E_m_der
                    // ), - u * ff
            );

        //A.assemble( var1(u,G) * var1(u,G).tr() );

        //A.assemble( var1(u,G)[0] );

        //A.assemble( deriv2(u) * var1(u,G).tr() ); // see comment 1

        //A.assemble( var2(u,u,G)[0] );

        // auto snormal = sn(G);

        // A.assemble( hessdot(u, vff )[0] ); // .normalized()

        //A.assemble( hessdot(u, var1(u,G)) );


        gsInfo<< A.rhs().transpose() <<"\n";
        gsInfo<< A.matrix().toDense().diagonal().transpose() <<"\n";

    } //for loop


    return EXIT_SUCCESS;

}// end main
