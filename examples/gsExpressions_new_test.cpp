/** @file gsExpressions_test.cpp

    @brief Testing integral computation using the expression evaluator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

namespace gismo
{
namespace expr
{


// Implementation for an addition
template <typename E1, typename E2>
class grad_expr<add_expr<E1, E2> > : public _expr<grad_expr<add_expr<E1, E2> > >
{
    const typename E1::Nested_t _u;
    const typename E2::Nested_t _v;

public:
    enum {Space = E1::Space, ScalarValued= E1::ScalarValued, ColBlocks= 0};

    typedef typename E1::Scalar Scalar;
    mutable gsMatrix<Scalar> uVals, uGrads, vVals, vGrads, tmp;

    grad_expr(const add_expr<E1, E2> & expr)
    :
    _u(expr.left()),
    _v(expr.right())
    {
        // GISMO_ASSERT(E1::Space == E2::Space,"Error: grad(x+v) requires u and v to have the same space.");
        // GISMO_ASSERT()
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto expr = grad(_u) + grad(_v);
        tmp = expr.eval(k);
        return tmp;
    }

    index_t rows() const { return 1 /*==u.dim()*/; }
    index_t cols() const { return _u.source().domainDim(); }

    index_t cardinality_impl() const
    { return _u.data().values[1].rows() / cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList); // WHY NEEDED??
        _v.parse(evList); // WHY NEEDED??
        grad(_u).parse(evList);
        grad(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; _u.print(os); os <<")"; }
};


// Implementation for an multiplication
template <typename E1, typename E2>
class grad_expr<mult_expr<E1, E2> > : public _expr<grad_expr<mult_expr<E1, E2> > >
{
    const typename E1::Nested_t _u;
    const typename E2::Nested_t _v;

public:
    enum {Space = E1::Space, ScalarValued= E1::ScalarValued, ColBlocks= 0};

    typedef typename E1::Scalar Scalar;
    mutable gsMatrix<Scalar> uVals, uGrads, vVals, vGrads, tmp;

    grad_expr(const mult_expr<E1, E2> & expr)
    :
    _u(expr.left()),
    _v(expr.right())
    {
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto expr = _v * grad(_u) + _u * grad(_v);
        tmp = expr.eval(k);
        return tmp;
    }

    index_t rows() const { return 1 /*==u.dim()*/; }
    index_t cols() const { return _u.source().domainDim(); }

    index_t cardinality_impl() const
    { return _u.data().values[1].rows() / cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList);
        _v.parse(evList);
        grad(_u).parse(evList);
        grad(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; _u.print(os); os <<")"; }
};

// Implementation for a division
template <typename E1, typename E2>
class grad_expr<divide_expr<E1, E2> > : public _expr<grad_expr<divide_expr<E1, E2> > >
{
    const typename E1::Nested_t _u;
    const typename E2::Nested_t _v;

public:
    enum {Space = E1::Space, ScalarValued= E1::ScalarValued, ColBlocks= 0};

    typedef typename E1::Scalar Scalar;
    mutable gsMatrix<Scalar> uVals, uGrads, vVals, vGrads, tmp;

    grad_expr(const divide_expr<E1, E2> & expr)
    :
    _u(expr.left()),
    _v(expr.right())
    {
        // GISMO_ASSERT()
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto expr = (grad(_u) * _v - _u * grad(_v)) / (_v * _v);
        tmp = expr.eval(k);
        return tmp;
    }

    index_t rows() const { return 1 /*==u.dim()*/; }
    index_t cols() const { return _u.source().domainDim(); }

    index_t cardinality_impl() const
    { return _u.data().values[1].rows() / cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList); // WHY NEEDED??
        _v.parse(evList); // WHY NEEDED??
        grad(_u).parse(evList);
        grad(_v).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; _u.print(os); os <<")"; }
};

// /*
//   Expression for the gradient of a finite element variable

//   Transposed gradient vectors are returned as a matrix
// */
// template<class E>
// class div_expr : public _expr<div_expr<E> >
// {
//     typename E::Nested_t _u;
// public:
//     enum {Space = E::Space, ScalarValued= 0, ColBlocks= 0};

//     typedef typename E::Scalar Scalar;
//     mutable gsMatrix<Scalar> tmp;

//     div_expr(const E & u)
//     :
//     _u(u)
//     {
//         GISMO_ASSERT(_u.source().domainDim() == _u.source().targetDim(),"Error: div(.) requires a vector field.");
//     }

//     const gsMatrix<Scalar> & eval(const index_t k) const
//     {
//         if (0!=Space)
//         {
//             // Dim x (numActive*Dim)
//             res = _u.data().values[1].col(k).transpose().blockDiag(_u.dim());
//         }
//         else
//         {
//             res = _u.data().values[1]
//                 .reshapeCol(k, _u.parDim(), _u.targetDim()).transpose()
//                 .blockDiag(_u.dim());
//         }
//         return res;
//     }

//     index_t rows() const { return 1 /*==u.dim()*/; }

//     index_t cols() const { return _u.source().domainDim(); }

//     index_t cardinality_impl() const
//     { return _u.data().values[1].rows() / cols(); }

//     void parse(gsExprHelper<Scalar> & evList) const
//     {
//         evList.add(_u);
//         _u.data().flags |= NEED_GRAD;
//     }

//     const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
//     const gsFeSpace<Scalar> & colVar() const
//     {return gsNullExpr<Scalar>::get();}

//     void print(std::ostream &os) const { os << "\u2207\u00B7("; _u.print(os); os <<")"; }
// private:

//     template<class U> static inline
//     typename util::enable_if<util::is_same<U,gsComposition<Scalar> >::value,
//                              const gsMatrix<Scalar> &>::type
//     eval_impl(const U & u, const index_t k)
//     {
//         return u.eval(k);
//     }
// };


}
}

using namespace gismo;

int main(int argc, char *argv[])
{



# define M_R  1.0

    bool verbose = false;
    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addSwitch("verbose", "Show result and exact", verbose);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo  <<"To do:\n"
            <<"- add deriv2 expressions to gsExpressions.h and make test case\n"
            <<"\n";

    gsMultiPatch<> mp;
    mp.addPatch(gsNurbsCreator<>::NurbsQuarterAnnulus(M_R,2*M_R));
    mp.degreeElevate();
    gsMultiBasis<> basis(mp);

    gsVector<> physpoint, point(2);
    point.setConstant(0.5);
    physpoint = mp.patch(0).eval(point);
    //b.basis(0).component(0).uniformRefine();

    // Set the expression assembler
    gsExprAssembler<> A(1,1);
    // Set the parameter mesh as the integration mesh
    A.setIntegrationElements(basis);

    gsExprEvaluator<> ev(A);
    gsMatrix<> coefs;

    gsFunctionExpr<> a_("x*y",2);
    gsFunctionExpr<> b_("x/y",2);
    // Define integrant variables
    auto G = A.getMap(mp);
    auto u = A.getSpace(basis);
    auto s = A.getSolution(u,coefs);
    auto a = A.getCoeff(a_, G);
    auto b = A.getCoeff(b_, G);
    gsMatrix<> result, exact, tmp;


    u.setup();
    coefs.setZero(u.mapper().freeSize(),1);



    gsInfo<<"-------------------------------------------------------------------------"<<"\n";
    gsInfo<<"------------------------------Simple stuff ------------------------"<<"\n";
    gsInfo<<"-------------------------------------------------------------------------"<<"\n";

    gsDebugVar(ev.eval(a*b,point));
    // gsDebugVar(ev.eval(a/b,point));
    gsDebugVar(ev.eval(a+b,point));
    gsDebugVar(ev.eval(a-b,point));

    gsDebugVar(ev.eval(1.0*a,point));
    gsDebugVar(ev.eval(a*1.0,point));
    // gsDebugVar(ev.eval(1.0+a,point)); //fails
    // gsDebugVar(ev.eval(a+1.0,point)); //fails
    // gsDebugVar(ev.eval(1.0-a,point)); //fails
    // gsDebugVar(ev.eval(a-1.0,point)); //fails
    // gsDebugVar(ev.eval(1.0/a,point)); //fails
    gsDebugVar(ev.eval(a/1.0,point));

    gsDebugVar(ev.eval(1.0*s,point));
    gsDebugVar(ev.eval(s*1.0,point));
    // gsDebugVar(ev.eval(1.0+s,point)); //fails
    // gsDebugVar(ev.eval(s+1.0,point)); //fails
    // gsDebugVar(ev.eval(1.0-s,point)); //fails
    // gsDebugVar(ev.eval(s-1.0,point)); //fails
    // gsDebugVar(ev.eval(1.0/s,point)); //fails
    gsDebugVar(ev.eval(s/1.0,point));

    // gsDebugVar(ev.eval(a*u,point)); //fails (runtime)
    gsDebugVar(ev.eval(u*a,point));
    // gsDebugVar(ev.eval(a+u,point)); //fails (apples and oranges. Expected behavior: a + u per basis function)
    // gsDebugVar(ev.eval(u+a,point)); //fails (apples and oranges. Expected behavior: u + a per basis function)
    // gsDebugVar(ev.eval(a-u,point)); //fails (apples and oranges. Expected behavior: a - u per basis function)
    // gsDebugVar(ev.eval(u-a,point)); //fails (apples and oranges. Expected behavior: u - a per basis function)
    // gsDebugVar(ev.eval(a/u,point)); //fails
    // gsDebugVar(ev.eval(u/a,point));

    gsDebugVar(ev.eval(1.0*u,point));
    gsDebugVar(ev.eval(u*1.0,point));
    // gsDebugVar(ev.eval(1.0+u,point)); //fails
    // gsDebugVar(ev.eval(u+1.0,point)); //fails
    // gsDebugVar(ev.eval(1.0-u,point)); //fails
    // gsDebugVar(ev.eval(u-1.0,point)); //fails
    // gsDebugVar(ev.eval(1.0/u,point)); //fails
    gsDebugVar(ev.eval(u/1.0,point));

    gsInfo<<"-------------------------------------------------------------------------"<<"\n";
    gsInfo<<"------------------------------Scalar Function operations------------------------"<<"\n";
    gsInfo<<"-------------------------------------------------------------------------"<<"\n";

    auto apb = a + b;
    auto amb = a - b;
    auto atb = a * b;
    auto ap1 = a + 1.0;

    gsDebug<<ev.eval( grad(a), point )<<"\n";
    // gsDebug<<ev.eval( grad(1.0*a), point )<<"\n"; //fails
    gsDebug<<ev.eval( grad(a+b), point )<<"\n";
    gsDebug<<ev.eval( grad(a*b), point )<<"\n";
    // gsDebug<<ev.eval( grad(a/b), point )<<"\n"; //fails

    // gsDebug<<ev.eval( grad(a+u), point )<<"\n"; //fails
    gsDebug<<ev.eval( grad(a*u), point )<<"\n";

    //auto c = ev.getVariable(c_, G);

    /*
      Computes the value of a variable
      Assessment of:
      - gsFeVariable
      = avg(gsFeVariable)
    */
    gsInfo<< "* Value:\t\t";
    result = ev.eval( a, point );
    exact  = a_.eval(physpoint);
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    /*
    gsInfo<< "* Average:\t\t";
    result.resize(1,1); result.at(0) = ev.integralInterface( avg(a) );
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    */

    /*
      Computes the gradient of a variable
      Assessment of:
      - grad_expr(gsFeVariable)
    */
    gsInfo<< "* Gradient:\t\t";
    result = ev.eval( grad(a), point );
    exact.resize(2,2);
    exact<<1,0,0,1;
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";

    /*
      Computes the Jacobian of a variable
      Assessment of:
      - fjac_expr(gsFeVariable)
    */
    gsInfo<< "* Function Jacobian:\t";
    result = ev.eval( jac(a), point );
    //exact.resize(3,2); //done above
    //exact<<2*physpoint(0,0),0,0,2*physpoint(1,0),physpoint(1,0),physpoint(0,0);
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n";
    gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    // /*
    //   Computes the Jacobian of a (1D) variable
    //   Assessment of:
    //   - jac_expr(gsFeVariable)
    // */
    // gsInfo<< "* 1-D Jacobian:\t\t";
    // result = ev.eval( jac(b), point );
    // exact.resize(1,2);
    // exact<<2*physpoint(0,0),2*physpoint(1,0);
    // gsInfo<<( (result-exact).norm() < 1e-10 ? "passed" : "failed" )<<"\n";
    // if (verbose)
    //     gsInfo  <<"Result:\n"<<result<<"\n"
    //             <<"Exact:\n"<<exact<<"\n";

    /*
      Computes the Laplacian of a variable
      Assessment of:
      - lapl_expr(gsFeVariable)
    */
    gsInfo<< "* Laplacian:\t\t";
    result = ev.eval( lapl(b), point );
    exact.resize(1,1);
    exact<<4;
    if (verbose)
        gsInfo  <<"Result:\n"<<result<<"\n"
                <<"Exact:\n"<<exact<<"\n"
                <<"Diff:\n"<< (result-exact) <<"\n";
    gsInfo<<( (result-exact).norm() < 1e-5 ? "passed" : "failed" );//<<"\n";
    gsInfo<<"\t\tnote: with lower tolerance"<<"\n";

    // /*
    //   Computes the Hessian of a variable
    //   Assessment of:
    //   - hess_expr(gsFeVariable)
    // */
    // gsInfo<< "* Hessian:\t\t";
    // result = ev.eval( hess(b), point );
    // exact.resize(2,2);
    // exact<<2,0,0,2;
    // if (verbose)
    //     gsInfo  <<"Result:\n"<<result<<"\n"
    //             <<"Exact:\n"<<exact<<"\n"
    //             <<"Diff:\n"<< (result-exact) <<"\n";
    // gsInfo<<( (result-exact).norm() < 1e-5 ? "passed" : "failed" );//<<"\n";
    // gsInfo<<"\t\tnote: with lower tolerance"<<"\n";

    // DOES NOT COMPILE
    /*
      Computes the curl of a variable
      Assessment of:
      - curl_expr(gsFeVariable)
    */
    /*
      gsInfo<< "* curl:\n";
      result = ev.eval( curl(a), point );
      gsInfo<< "  Result:\n"<< result <<"\n";
      exact.resize(3,1);
      exact<<0,0,1;
      gsInfo<< "  Exact:\n"<< exact <<"\n";
    */



    return EXIT_SUCCESS;
}
