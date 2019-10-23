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
        res.resize(rows(), cols());
        const index_t A = rows()/cols(); // note: rows/cols is the number of actives

        normal = _G.data().normal(k);// not normalized to unit length
        bGrads = _u.data().values[1].col(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        for (index_t d = 0; d!= cols(); ++d) // for all basis functions (1)
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all basis functions (2)
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
        return cols() * _u.data().values[1].rows() / _u.source().domainDim();
    }

    index_t cols() const { return _u.dim(); }

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_OUTER_NORMAL | NEED_DERIV | NEED_MEASURE;
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

    var2_expr(const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & G) : _u(u),_v(v), _G(G) { }

    mutable gsMatrix<Scalar> res;

    mutable gsMatrix<Scalar> uGrads, vGrads, cJac, cDer2;
    mutable gsVector<Scalar,3> m_v, m_w, normal, m_vw, m_v_der, n_der, n_der2, tmp;
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

        const index_t numAct = _u.data().values[0].rows();
        const index_t numHess = cDer2.rows();
        const Scalar measure =  _G.data().measures.at(k);

        for (index_t d = 0; d!= _u.dim(); ++d) // for all basis functions u (1)
        {
            const short_t s = d*numAct;
            for (index_t j = 0; j!= numAct; ++j) // for all basis functions u (2)
            {
                m_v.noalias() = vecFun(d, uGrads.at(2*j  ) ).cross( cJac.col(1).template head<3>() )
                              - vecFun(d, uGrads.at(2*j+1) ).cross( cJac.col(0).template head<3>() ) / measure;

                for (index_t c = 0; c!= _v.dim(); ++c) // for all basis functions v (1)
                {
                    const short_t r = c*numAct;
                    for (index_t i = 0; i!= numAct; ++i) // for all basis functions v (2)
                    {

                        m_w.noalias() = vecFun(c, vGrads.at(2*i  ) ).cross( cJac.col(1).template head<3>() )
                                      - vecFun(c, vGrads.at(2*i+1) ).cross( cJac.col(0).template head<3>() ) / measure;

                        n_der = (m_w - ( normal.dot(m_w) ) * normal);

                        m_vw.noalias() = vecFun(d, uGrads.at(2*j  ) ).cross( vecFun(c, vGrads.at(2*i+1) ) )
                                       - vecFun(d, uGrads.at(2*j+1) ).cross( vecFun(c, vGrads.at(2*i  ) ) ) / measure;

                        m_v_der.noalias() = (m_vw - ( normal.dot(m_w) ) * m_v);

                        // ---------------  Second variation of the normal
                        tmp = m_v_der - (m_v.dot(n_der) + normal.dot(m_v_der) ) * normal - (normal.dot(m_v) ) * n_der;

                        for (index_t l=0; l != numHess; l++) // per hessian entry of c
                        {
                            res(s + j, r + i + l*_u.dim()*numAct ) = tmp.dot(cDer2.col(l));
                        }
                    }
                }
            }
        }

//        gsDebugVar(res);
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

    mutable gsMatrix<Scalar> res,tmp;

    const gsMatrix<Scalar> & eval(const index_t k) const { return eval_impl(_u,k); }

    index_t rows() const
    {
        return _u.targetDim() * _u.data().values[0].rows(); // no. dimensions*numAct
    }

    index_t cols() const
    {
        return _u.data().values[2].rows() / _u.data().values[0].rows(); // numHessian dimensions
    }

    void setFlag() const
    {
        _u.data().flags |= NEED_DERIV2;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_DERIV2;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    static constexpr bool rowSpan() {return E::rowSpan(); }
    static constexpr bool colSpan() {return E::colSpan();}

    void print(std::ostream &os) const { os << "deriv2("; _u.print(os); os <<")"; }

private:
    template<class U> inline
    typename util::enable_if< util::is_same<U,gsGeometryMap<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        tmp =_u.data().values[2].reshapeCol(k, _u.data().dim.second, _u.data().dim.second);
        for (index_t i = 0; i!=3; i++)
            res.block(i*3,0,3,3) = tmp;
        return res;
    }

    template<class U> inline
    typename util::enable_if<util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k) const
    {
        res.resize(rows(), cols());
        tmp.resize(rows()/_u.dim(), cols());
        const index_t numAct = _u.data().values[0].rows();
        const index_t numHess = _u.data().values[2].rows() / numAct;
//        gsDebugVar(numHess);
        tmp = gsAsConstMatrix<Scalar>(_u.data().values[2].col(k).data(), 3, numAct ).transpose();

//        gsDebugVar(cols());
//        gsDebugVar(rows());

        for (index_t i = 0; i!=_u.dim(); i++)
            res.block(i*numAct,0,numAct,numHess) = tmp;

        gsDebugVar(res); // 1: WHY DOES IN THE PRODUCT THE NUMBER OF ROWS DOUBLE?
        return res;
    }

};



template<class E1, class E2>
class hessdot_expr : public _expr<hessdot_expr<E1,E2> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum{ Space = E1::Space };

    typedef typename E1::Scalar Scalar;

    hessdot_expr(const E1 & u, const E2 & v) : _u(u), _v(v) {}

    mutable gsMatrix<Scalar> res, hess, tmp;
    mutable gsMatrix<Scalar> normalMat;

    MatExprType eval(const index_t k) const
    {
        const gsFuncData<Scalar> & udata = _u.data(); // udata.values[2].col(k)
        const index_t numAct = udata.values[0].rows();
        const gsAsConstMatrix<Scalar> ders(udata.values[2].col(k).data(), 3, numAct );

        tmp = _v.eval(k);

        res.resize(rows(), cols() );


            for (index_t i = 0; i!=tmp.rows(); ++i)
            {
                res.block(i*numAct, 0, numAct, 3).noalias() = ders.transpose() * tmp.at(i);
            }

        return res;
    }

    index_t rows() const
    {
        return _u.dim() * _u.data().values[0].rows();
    }

    index_t cols() const
    {
        return // ( E2::rowSpan() ? rows() : 1 ) *
            _u.data().values[2].rows() / _u.data().values[0].rows();//=3
    }

    void setFlag() const
    {
        _u.data().flags |= NEED_2ND_DER;
        _v.setFlag();
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_2ND_DER;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    static constexpr bool rowSpan() {return E1::rowSpan(); }
    static constexpr bool colSpan() {return E2::rowSpan(); }

    void print(std::ostream &os) const { os << "hessdot("; _u.print(os); os <<")"; }
};

template<class E> EIGEN_STRONG_INLINE
var1_expr<E> var1(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return var1_expr<E>(u, G); }

template<class E1, class E2> EIGEN_STRONG_INLINE
var2_expr<E1,E2> var2(const E1 & u, const E2 & v, const gsGeometryMap<typename E1::Scalar> & G)
{ return var2_expr<E1,E2>(u,v, G); }

template<class E1, class E2> EIGEN_STRONG_INLINE
hessdot_expr<E1,E2> hessdot(const E1 & u, const E2 & v) { return hessdot_expr<E1,E2>(u, v); }

template<class E> EIGEN_STRONG_INLINE
deriv2_expr<E> deriv2(const E & u) { return deriv2_expr<E>(u); }

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

    gsFunctionExpr<> materialMat("1","1","1","1","1","1","1","1","1",3);
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

        auto E_m = 0.5 * ( flat(jac(defG).tr()*jac(defG)) - flat(jac(G).tr()* jac(G)) );
        auto E_m_der = flat( jac(u).tr() * jac(defG) ) ; 
        auto E_m_der2 = flat( jac(u) * jac(u).tr() );

        auto E_f = reshape(m2,2,2) * ( deriv2(defG) * nv(defG).normalized() - deriv2(G) * nv(G).normalized() ) ;
        auto E_f_der = reshape(m2,2,2) * ( deriv2(u) * nv(defG).normalized() + deriv2(defG) * var1(u,G) );
        auto E_f_der2 = reshape(m2,2,2) * ( 2 * deriv2(u) * var1(u,defG).tr() + var2(u,u,defG) );

        // A.assemble( tt * ( E_m * mm * E_m_der2 + E_m_der * mm * E_m_der ) +
        //             (tt*tt*tt/3.0) *( E_f * mm * E_f_der2 + E_f_der * mm * E_f_der ), // Matrix

        //             tt * ( E_m * mm * E_m_der ) + (tt*tt*tt/3.0) *( E_f * mm * E_f_der )  // Rhs
        //             - u * ff.tr()
        //     );

        A.assemble( tt * ( E_m_der.tr() * reshape(mm,3,3) * E_m_der ), - u * ff );

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
