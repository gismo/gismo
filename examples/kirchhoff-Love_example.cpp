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

template<class E>
class mygrad_expr : public _expr<mygrad_expr<E> >
{
    typename E::Nested_t _u;
public:
    enum{ Space = E::Space };

    typedef typename E::Scalar Scalar;

    mygrad_expr(const E & u) : _u(u)
    { GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");}

    MatExprType eval(const index_t k) const
    {
        // numActive x dim
        return _u.data().values[1].reshapeCol(k, cols(), rows()).transpose();
    }

    index_t rows() const
    {
        //return _u.data().values[0].rows();
        return _u.data().values[1].rows() / cols();
    }
    //index_t rows() const { return _u.data().actives.size(); }
    //index_t rows() const { return _u.rows(); }

    //index_t rows() const { return _u.source().targetDim() is wrong }
    index_t cols() const { return _u.source().domainDim(); }

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    static constexpr bool rowSpan() {return E::rowSpan(); }
    static constexpr bool colSpan() {return false;}

    void print(std::ostream &os) const { os << "grad("; _u.print(os); os <<")"; }
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

        normal = _G.data().outNormals.col(k);// not normalized to unit length
        bGrads = _u.data().values[1].col(k);
        cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        const Scalar measure =  _G.data().measures.at(k);

        for (index_t d = 0; d!= cols(); ++d) // for all basis functions (1)
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all basis functions (2)
            {
                // Jac(u) ~ Jac(G) with alternating signs ?..
                m_v.noalias() = vecFun(d, bGrads.at(2*j  ) ).cross( cJac.col(1).template head<3>() )
                              - vecFun(d, bGrads.at(2*j+1) ).cross( cJac.col(0).template head<3>() );

                // ---------------  First variation of the normal
                res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose() / measure;
            }
        }

        gsDebugVar(res);
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
        _G.data().flags |= NEED_OUTER_NORMAL | NEED_DERIV | NEED_MEASURE;
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

// template<class E>
// class var2_expr : public _expr<var2_expr<E> >
// {
//     typename E::Nested_t _u;
//     typename E::Nested_t _v;
//     typename gsGeometryMap<T>::Nested_t _G;

// public:
//     enum{ Space = E::Space };

//     typedef typename E::Scalar Scalar;

//     var2_expr(const E & u) : _u(u)
//     { GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");}

//     MatExprType eval(const index_t k) const
//     {
//         // numActive x dim
//         return _u.data().values[1].reshapeCol(k, cols(), rows()).transpose();
//     }

//     index_t rows() const
//     {
//         //return _u.data().values[0].rows();
//         return _u.data().values[1].rows() / cols();
//     }
//     //index_t rows() const { return _u.data().actives.size(); }
//     //index_t rows() const { return _u.rows(); }

//     //index_t rows() const { return _u.source().targetDim() is wrong }
//     index_t cols() const { return _u.source().domainDim(); }

//     void setFlag() const
//     {
//         _u.data().flags |= NEED_GRAD;
//         _G.data().flags |= NEED_OUTER_NORMAL;
//         if (_u.composed() )
//             _u.mapData().flags |= NEED_VALUE;
//     }

//     void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
//     {
//         //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
//         evList.push_sorted_unique(&_u.source());
//         _u.data().flags |= NEED_GRAD;
//         if (_u.composed() )
//             _u.mapData().flags |= NEED_VALUE;
//     }

//     const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
//     const gsFeSpace<Scalar> & colVar() const
//     {return gsNullExpr<Scalar>::get();}

//     static constexpr bool rowSpan() {return E::rowSpan(); }
//     static constexpr bool colSpan() {return false;}

//     void print(std::ostream &os) const { os << "grad("; _u.print(os); os <<")"; }
// };

// template<class E>
// class hessdot_expr : public _expr<hessdot_expr<E> >
// {
//     typename E::Nested_t _u;
//     typename gsGeometryMap<T>::Nested_t _G;

// public:
//     enum{ Space = E::Space };

//     typedef typename E::Scalar Scalar;

//     hessdot_expr(const E & u) : _u(u)
//     { GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");}

// !    MatExprType eval(const index_t k) const
//     {
//         // numActive x dim
//         return _u.data().values[1].reshapeCol(k, cols(), rows()).transpose();
//     }

// !    index_t rows() const
//     {
//         //return _u.data().values[0].rows();
//         return _u.data().values[1].rows() / cols();
//     }
//     //index_t rows() const { return _u.data().actives.size(); }
//     //index_t rows() const { return _u.rows(); }

//     //index_t rows() const { return _u.source().targetDim() is wrong }
// !    index_t cols() const { return _u.source().domainDim(); }

//     void setFlag() const
//     {
//         _u.data().flags |= NEED_2ND_DER;
//         _G.data().flags |= NEED_OUTER_NORMAL;
//     }

//     void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
//     {
//         //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
//         evList.push_sorted_unique(&_u.source());
//         _u.data().flags |= NEED_GRAD;
//         if (_u.composed() )
//             _u.mapData().flags |= NEED_VALUE;
//     }

// !    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
// !    const gsFeSpace<Scalar> & colVar() const
//     {return gsNullExpr<Scalar>::get();}

//     static constexpr bool rowSpan() {return E::rowSpan(); }
//     static constexpr bool colSpan() {return false;}

//     void print(std::ostream &os) const { os << "grad("; _u.print(os); os <<")"; }
// };



template<class E> EIGEN_STRONG_INLINE
mygrad_expr<E> mygrad(const E & u) { return mygrad_expr<E>(u); }

template<class E> EIGEN_STRONG_INLINE
var1_expr<E> var1(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return var1_expr<E>(u, G); }

// template<class E> EIGEN_STRONG_INLINE
// mygrad_expr<E> var2(const E & u, const E & v, const gsGeometryMap<T> & G) { return var2_expr<E>(u, v, G); }

// template<class E> EIGEN_STRONG_INLINE
// mygrad_expr<E> hessdot(const E & u, const gsGeometryMap<T> & G) { return hessdot_expr<E>(u, G); }

}}

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
    //fd.getId(0, mp); // id=0: Multipatch domain

    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree

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

    // Set the discretization space
    space u = A.getSpace(dbasis, 3);
    u.setInterfaceCont(0);

    // Solution vector and solution variable
    gsMatrix<> solVector;
    solution u_sol = A.getSolution(u, solVector);

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
//        A.assemble( mygrad(u)*jac(G).ginv() * (mygrad(u)*jac(G).ginv()).tr() * meas(G), u * ff * meas(G) );

        //A.assemble( var1(u,G) * var1(u,G).tr() );

        A.assemble( var1(u,G)[0] );

        gsInfo<< A.rhs().transpose() <<"\n";
        //gsInfo<< A.matrix().toDense().diagonal().transpose() <<"\n";

    } //for loop


    return EXIT_SUCCESS;

}// end main
