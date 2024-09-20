/** @file monitor_poisson_composed_r-adaptivity.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsHLBFGS/gsHLBFGS.h>
#include <gsModeling/gsBarrierCore.h>
#include <gsOptimizer/gsGradientDescent.h>
#include <gsOptim/gsOptim.h>

//! [Include namespace]

namespace gismo{
namespace expr{

template<class E>
class monitor_expr : public _expr<monitor_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;
    mutable gsMatrix<Scalar> grad, jac, jacInv;
    mutable gsVector<Scalar> ones;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _G;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    monitor_expr(const E & u, const gsGeometryMap<Scalar> & G) : _u(u), _G(G)
    {
    }

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {return eval_impl(_u,k); }

    index_t rows() const { return _u.source().domainDim(); }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        grad_expr<gsFeVariable<Scalar>>(_u).parse(evList);
        jacInv_expr<Scalar>(_G).parse(evList);
        _u.parse(evList);
        _G.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "M("; _u.print(os); os <<")"; }

private:
    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        // const index_t A = _u.cardinality()/_u.dim(); // _u.data().actives.rows()
        // res.resize(_u.cardinality(), cols()); // rows()*

        // normal = _G.data().normal(k);// not normalized to unit length
        // normal.normalize();
        // bGrads = _u.data().values[1].col(k);
        // cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        // const Scalar measure =  (cJac.col3d(0).cross( cJac.col3d(1) )).norm();

        // for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        // {
        //     const short_t s = d*A;
        //     for (index_t j = 0; j!= A; ++j) // for all actives
        //     {
        //         // Jac(u) ~ Jac(G) with alternating signs ?..
        //         m_v.noalias() = (vecFun(d, bGrads.at(2*j  ) ).cross( cJac.col(1).template head<3>() )
        //                       - vecFun(d, bGrads.at(2*j+1) ).cross( cJac.col(0).template head<3>() )) / measure;

        //         // ---------------  First variation of the normal
        //         // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
        //         res.row(s+j).noalias() = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
        //     }
        // }
        // return res;
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        // GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");
        // grad_expr<gsFeSolution<Scalar>> sGrad =  grad_expr<gsFeSolution<Scalar>>(_u);
        // grad = sGrad.eval(k);
        // ones.resize( _u.source().domainDim());
        // ones.setOnes();
        // //ones[0] = 4;
        // res = 1.0 / ( math::sqrt(1.0 + grad.dot(grad) ) ) * ones;
        // //res = ones;
        // return res;
    }

    template<class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        grad_expr<gsFeVariable<Scalar>> gradExpr = grad_expr<gsFeVariable<Scalar>>(_u);
        grad = gradExpr.eval(k);

        jacInv_expr<Scalar> jacInvExpr = jacInv_expr<Scalar>(_G);
        jacInv = jacInvExpr.eval(k);

        grad = grad * jacInv;
        ones.resize( _u.source().domainDim());
        ones.setOnes();

        res.resize( _u.source().domainDim(), 1);

        res = 1.0 / ( math::sqrt(1.0 + grad.squaredNorm() ) ) * ones;
        return res;

        // grad = _u.data().values[1].col(k);
        // ones.resize( _u.source().domainDim());
        // ones.setOnes();

        // res.resize( _u.source().domainDim(), 1);

        // res = 1.0 / ( math::sqrt(1.0 + grad.squaredNorm() ) ) * ones;
        // return res;
    }
};

template<class E> EIGEN_STRONG_INLINE
monitor_expr<E> monitor(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return monitor_expr<E>(u,G); }


}
}

using namespace gismo;

// template<typename T = real_t>
// class gsExprAsFunction : public gsFunction<T>
// {
// public:
//     template<class E>
//     gsExprAsFunction( const expr::_expr<E> & expr)
//     : m_expr(expr)
//     {}

//     short_t domainDim() const { return m_expr.parDim(); }
//     short_t targetDim() const { return m_expr.rows();   }

//     void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
//     {
//         result.resize(u,this->targetDim());
//         gsExprEvaluator<T> ev;
//         for (index_t k = 0; k!=u.cols(); k++)
//             result.col(k) = ev.eval(m_expr,u.col(k));
//     }

// protected:
//     const expr::_expr<E> & m_expr;
// };




template<typename T = real_t>
class gsOptMesh : public gsOptProblem<T> 
{
    using Base = gsOptProblem<T>;

private:
    typedef typename gsExprAssembler<T>::geometryMap geometryMap;
    typedef typename gsExprAssembler<T>::space space;
    typedef typename gsExprAssembler<T>::solution solution;

public:
    gsOptMesh(  gsFunction<T> & composition,
                const gsGeometry<T> & geometry,
                const gsFunctionSet<T> & fun,
                T eps = 1e-4)
    :
    m_comp(&composition),
    m_geom(geometry),
    m_fun(fun),
    m_eps(eps)
    {
        m_numDesignVars = m_comp->nControls();
        m_curDesign.resize(m_numDesignVars,1);
    }

    /// Evaluates the objective function at the given point u.
    T evalObj(const gsAsConstVector<T> &u) const override
    {
        for (index_t k=0; k!=u.rows(); k++)
            m_comp->control(k) = u[k];

        gsComposedGeometry<T> cgeom(*m_comp,m_geom);

        gsMultiPatch<> mp;
        mp.addPatch(cgeom);
        gsMultiBasis<> mb(m_geom.basis());

        m_evaluator.setIntegrationElements(mb);
        geometryMap G = m_evaluator.getMap(mp);

        gsConstantFunction<T> epsilonFunction(m_eps, m_geom.domainDim());
        auto eps = m_evaluator.getVariable(epsilonFunction);

        auto chi = 0.5 * (jac(G).det() + pow(eps.val() + pow(jac(G).det(), 2.0), 0.5));
        auto invJacMat = jac(G).adj()/chi;
        auto eta = m_evaluator.getVariable(m_fun);
        return m_evaluator.integral( (monitor(eta,G).asDiag()*invJacMat).sqNorm()*meas(G));
    }

    // /// Computes the gradient of the objective function at the given point u
    // // and stores it in result.
    // void gradObj_into(const gsAsConstVector<T> &u,
    //                 gsAsVector<T> &result) const override;

protected:

    // From gsOptProblem
    using Base::m_curDesign;
    using Base::m_numDesignVars;


    gsFunction<T> * m_comp;
    const gsGeometry<T> & m_geom;
    const gsFunctionSet<T> & m_fun;
    T m_eps;

    mutable gsExprEvaluator<T> m_evaluator;
};


int main(int arg, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 2;
    index_t numRefineD = 2;
    index_t numElevate = 0;
    index_t numElevateD= 0;
    index_t maxIt = 100;
    real_t tol_g = 5e-5;
    real_t eps = 1e-4;
    bool slide = true;
    index_t expression = 0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "elevAnalysis","Number of degree elevation steps to perform for the analysis", numElevate );
    cmd.addInt( "E", "elevDomain","Number of degree elevation steps to perform for the domain", numElevateD );
    cmd.addInt( "r", "refAnalysis", "Number of Uniform h-refinement loops for the analysis",  numRefine );
    cmd.addInt( "R", "refDomain", "Number of Uniform h-refinement loops for the domain",  numRefineD );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addReal("g", "tolG", "relative tol", tol_g);
    cmd.addInt( "i", "maxIt", "max num iterations",  maxIt );
    cmd.addReal("", "eps", "eps",  eps );
    cmd.addSwitch("noslide", "Do not slide the boundaries",  slide );
    cmd.addInt( "f", "expr", "Problem to be solved: 0 default peak in the center of the domain; 1 peak in the bottom left corner of the domain.",  expression );

    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    GISMO_ASSERT(numElevateD<=numElevate,"The degree of the domain should be lower or equal to the degree of the analysis");
    GISMO_ASSERT(numRefineD<=numRefine,"The number of refinements for the domain should be lower or equal to the number of refinements of the analysis");

    std::string dirname = "r-adaptivity_composed_poisson_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine)+"_E"+util::to_string(numElevateD)+"_R"+util::to_string(numRefineD);
    gsFileManager::mkdir(dirname);
    dirname += gsFileManager::getNativePathSeparator();

    // Basis for the analysis
    gsKnotVector<> kv({0,0,1,1},1);
    gsTensorBSplineBasis<2> tbasis(kv,kv);
        
    tbasis.degreeElevate(numElevate);
    for (index_t i = 0; i < numRefine; i++)
        tbasis.uniformRefine();

    gsInfo<<"Analysis basis:\n"<<tbasis<<"\n";

    gsMatrix<> coefs = tbasis.anchors().transpose();
    // coefs.row(coefs.rows()-1) *= 1.2;
    gsTensorBSpline<2> tbspline(tbasis,coefs);

    // gsSquareDomain<2,real_t> domain(numElevate,numRefine);

    // Basis for the square domain
    gsKnotVector<> kv2({0,0,1,1},1);
    gsTensorBSplineBasis<2> dbasis(kv2,kv2);
    dbasis.degreeElevate(numElevateD);
    for (index_t i = 0; i < numRefineD; i++)
        dbasis.uniformRefine();
    
    gsInfo<<"Mapper basis:\n"<<dbasis<<"\n";

    gsSquareDomain<2,real_t> domain(dbasis);
    domain.options().addSwitch("Slide","",slide);
    domain.applyOptions();

    gsComposedGeometry<real_t> cspline(domain,tbspline);


    gsMultiPatch<> mp;
    mp.addPatch(cspline);
    // mp.embed(3);
    gsMultiBasis<> mb(mp);

/* 
    SOLVE POISSON
 */
    
    // Source function:
    
    index_t dimexpr = mp.geoDim();
    std::string fstring = "";
    std::string msstring = "";

    switch (expression) {
    case 1:
        // bottom left corner
        fstring = "((tanh(20*(x^2 + y^2)^(1/2) - 5)^2 - 1)*(20*x^2 + 20*y^2)*(40*tanh(20*(x^2 + y^2)^(1/2) - 5)*(x^2 + y^2)^(1/2) - 1))/(x^2 + y^2)^(3/2)";
        msstring = "tanh((0.25-sqrt(x^2+y^2))/0.05)+1";
        break;
    default:
        // center of the domain
        fstring = "((tanh(20*((x-0.5)^2 + (y-0.5)^2)^(1/2) - 5)^2 - 1)*(20*(x-0.5)^2 + 20*(y-0.5)^2)*(40*tanh(20*((x-0.5)^2 + (y-0.5)^2)^(1/2) - 5)*((x-0.5)^2 + (y-0.5)^2)^(1/2) - 1))/((x-0.5)^2 + (y-0.5)^2)^(3/2)";
        msstring = "tanh((0.25-sqrt((x-0.5)^2+(y-0.5)^2))/0.05)+1";
        break;
    }

    gsFunctionExpr<> f(fstring, dimexpr);
    gsFunctionExpr<> ms(msstring,dimexpr);

    gsBoundaryConditions<> bc;
    bc.addCondition(boundary::side::west ,condition_type::dirichlet,&ms);
    bc.addCondition(boundary::side::east ,condition_type::dirichlet,&ms);
    bc.addCondition(boundary::side::south,condition_type::dirichlet,&ms);
    bc.addCondition(boundary::side::north,condition_type::dirichlet,&ms);
    bc.setGeoMap(mp);

    //! [Problem setup]
    gsExprAssembler<> A(1,1);

    typedef typename gsExprAssembler<>::geometryMap geometryMap;
    typedef typename gsExprAssembler<>::variable    variable;
    typedef typename gsExprAssembler<>::space       space;
    typedef typename gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    A.setIntegrationElements(mb);
    gsExprEvaluator<> ev(A);
    gsMatrix<> solVector;

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space u = A.getSpace(mb);

    // Set the source term
    auto ff = A.getCoeff(f, G);

    // Recover manufactured solution
    auto u_ex = ev.getVariable(ms, G);

    // Solution vector and solution variable
    solution u_sol = A.getSolution(u, solVector);

    //! [Problem setup]

    typename gsSparseSolver<>::CGDiagonal solver;

    u.setup(bc, dirichlet::l2Projection, 0);

    // Initialize the system
    A.initSystem();
    gsInfo<<"Number of analysis degrees of freedom: "<<A.numDofs()<<"\n";

    // Compute the system matrix and right-hand side
    A.assemble(
        igrad(u, G) * igrad(u, G).tr() * meas(G) //matrix
        ,
        u * ff * meas(G) //rhs vector
        );

    solver.compute( A.matrix() );
    solVector = solver.solve(A.rhs());

    // Convert solution to gsMultiPatch
    gsMultiPatch<> fun;
    u_sol.extract(fun);
    ev.writeParaview(u_sol,G,dirname+"solution");
    ev.writeParaview(ijac(u_sol,G),G,dirname+"solution_gradient");
    ev.writeParaview(ijac(u_sol,G).sqNorm(),G,dirname+"solution_gradient_sqNorm");


/* 
    PERFORM R-ADAPTIVITY
 */

    gsInfo<<"Number of optimizer degrees of freedom: "<<domain.nControls()<<"\n";

    gsOptMesh<> optMesh(domain,tbspline,fun,eps);
    gsVector<> controls(domain.nControls());
     for (size_t k=0; k!=domain.nControls(); k++)
        controls[k] = domain.control(k);


    //gsHLBFGS<real_t> optimizer(&optMesh);
    // gsGradientDescent<real_t> optimizer(&optMesh);
    //optimizer.options().setInt("MaxIterations",maxIt);
    //optimizer.options().setInt("Verbose",2);
    //optimizer.options().setReal("tolRelG",tol_g);
    // // gsDebugVar(optimizer.currentDesign().transpose());
    //optimizer.solve(controls);
    //gsVector<> optSol = optimizer.currentDesign();

    gsVector<> optSol;

    gsOptim<real_t>::LBFGS optim(&optMesh);
    optim.options().setInt("MaxIterations",maxIt);
    optim.options().setInt("Verbose",1);
    optim.options().setReal("GradErrTol",tol_g);
    optim.solve(controls);
    optSol = optim.currentDesign();

    // gsGradientDescent<real_t> hlbfgs(&optMesh);
    // hlbfgs.options().setInt("MaxIterations",maxIt);
    // hlbfgs.options().setInt("Verbose",1);
    // // hlbfgs.options().setReal("tolRelG",tol_g);
    // hlbfgs.solve(optSol);
    // optSol = hlbfgs.currentDesign();
    



    
    
    for (size_t k=0; k!=optSol.rows(); k++)
        domain.control(k) = optSol[k];



    // gsWriteParaview(cspline,"cspline",1000,true);
    // gsWriteParaview(cspline.basis(),"cbasis",1000);
    gsWriteParaview(domain.domain(),dirname+"domain",1000,true,true);

    gsInfo<<"Area = "<<ev.integral(meas(G))<<"\n";
    ev.writeParaview(jac(G).det(),G,dirname+"jacobian_determinant");

    return 0;
}// end main
