/** @file monitor_template_composed_r-adaptivity.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsCore/gsComposedFunction.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsHLBFGS/gsHLBFGS.h>
#include <gsModeling/gsBarrierCore.h>
#include <gsOptimizer/gsGradientDescent.h>
#include <gsOptim/gsOptim.h>

//! [Include namespace]

enum MonitorMode
{
    // We can add the option "Runtime" where the mode is determined at runtime using an option
    ValueBased = 0,
    GradientBased = 1
};

namespace gismo{
namespace expr{

template<enum MonitorMode MODE, class E>
class monitor_expr : public _expr<monitor_expr<MODE,E> >
{
public:
    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;
    mutable gsMatrix<Scalar> grad, jac, jacInv;
    mutable gsMatrix<Scalar> ones;
    mutable Scalar m_theta;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _G;
    const short_t DIM;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    monitor_expr(const E & u, const gsGeometryMap<Scalar> & G, const Scalar theta)
    :
    _u(u),
    _G(G),
    DIM(_u.source().domainDim()),
    m_theta(theta)
    {
    }

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {return eval_impl<MODE>(_u,k); }

    index_t rows() const { return DIM; }

    index_t cols() const { return DIM; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        parse_impl<MODE>(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "M("; _u.print(os); os <<")"; }

private:

    template<enum MonitorMode _MODE> inline
    typename util::enable_if<_MODE == ValueBased, void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        _G.parse(evList);

        evList.add(_u);
        _u.data().flags |= NEED_VALUE; //for value-based
    }

    template<enum MonitorMode _MODE> inline
    typename util::enable_if<_MODE == GradientBased, void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        jacInv_expr<Scalar>(_G).parse(evList);
        _G.parse(evList);

        evList.add(_u);
        _u.data().flags |= NEED_JACOBIAN;
    }

    template<enum MonitorMode _MODE> inline
    typename util::enable_if< _MODE != ValueBased && _MODE != GradientBased , void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

    template<enum MonitorMode _MODE, class U> inline
    typename util::enable_if< util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        GISMO_NO_IMPLEMENTATION;
    }

    template<enum MonitorMode _MODE, class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSpace<Scalar> >::value && _MODE == ValueBased, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        ones.resize(DIM,DIM);
        ones.setIdentity();

        // Value-based:
        Scalar uval = _u.data().values[0].col(k).value();
        res = 1.0 / ( math::sqrt(1.0 +m_theta*math::pow(uval,2) ) ) * ones;
        return res;
    }

    template<enum MonitorMode _MODE, class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSpace<Scalar> >::value && _MODE == GradientBased, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        ones.resize(DIM,DIM);
        ones.setIdentity();

        // Gradient-based
        jac = _u.data().values[1].col(k).transpose();
        jacInv_expr<Scalar> jacInvExpr = jacInv_expr<Scalar>(_G);
        jacInv = jacInvExpr.eval(k);
        grad = jac * jacInv;

        res = 1.0 / ( math::sqrt(1.0 +m_theta*grad.squaredNorm() ) ) * ones;
        return res;
    }

    template<enum MonitorMode _MODE, class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSpace<Scalar> >::value && _MODE != ValueBased && _MODE != GradientBased, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        GISMO_NO_IMPLEMENTATION;
    }

};

template<enum MonitorMode MODE, class E>
EIGEN_STRONG_INLINE
monitor_expr<MODE,E> monitor(const E & u,
                             const gsGeometryMap<typename E::Scalar> & G,
                             const typename E::Scalar theta = 0.1)
                             {
                                return monitor_expr<MODE,E>(u,G,theta);
                             }


}
}

using namespace gismo;


template<typename T, enum MonitorMode MODE>
class gsOptMesh : public gsOptProblem<T>
{
    using Base = gsOptProblem<T>;

private:
    typedef typename gsExprAssembler<T>::geometryMap geometryMap;
    typedef typename gsExprAssembler<T>::space space;
    typedef typename gsExprAssembler<T>::solution solution;

public:
    gsOptMesh(        gsFunction<T> & composition,
                const gsGeometry<T> & geometry,
                const gsFunction<T> & fun,
                const bool            parametric)
    :
    m_comp(composition)
    {


        // Assert dimensions
        // GISMO_ASSERT();

        m_cgeom = gsComposedGeometry<T>(m_comp,geometry);
        if (parametric)
            m_cfun = gsComposedFunction<T>(m_comp,fun);
        else
            m_cfun = gsComposedFunction<T>(m_cgeom,fun);

        m_numDesignVars = m_comp.nControls();
        m_curDesign.resize(m_numDesignVars,1);
        m_controls.resize(m_numDesignVars,1);
        for (index_t k=0; k!=m_numDesignVars; k++)
            m_controls(k,0) = m_comp.control(k);

        m_options.addInt("nRefine","Number of refinement steps for the integration basis",2);
        m_options.addInt("nElevate","Number of elevation steps for the integration basis",1);

        m_options.addReal("Smoothing","Smoothing parameter for the monitor function",0.1);
        m_options.addReal("Penalty","Penalty parameter for the monitor function",1e-2);


    }

    gsFunction<T> & composition() { return m_comp; }

    gsOptionList & options() { return m_options; }

    /// Evaluates the objective function at the given point u.
    T evalObj(const gsAsConstVector<T> &u) const override
    {
        for (index_t k=0; k!=m_numDesignVars; k++)
            m_comp.control(k) = u[k];

        gsMultiPatch<> mpG;
        mpG.addPatch(m_cgeom);

        gsKnotVector<> kv({0,0,1,1},1);
        gsTensorBSplineBasis<2,T> basis(kv,kv);
        basis.degreeElevate(m_options.getInt("nElevate"));
        for (index_t i = 0; i < m_options.getInt("nRefine"); i++)
            basis.uniformRefine();

        gsMultiBasis<> mb(basis);

        m_evaluator.setIntegrationElements(mb);
        geometryMap G = m_evaluator.getMap(mpG);

        gsConstantFunction<T> epsilonFunction(m_options.getReal("Penalty"), m_cgeom.domainDim());
        auto eps = m_evaluator.getVariable(epsilonFunction);

        // auto chi = 0.5 * (jac(G).det() + pow(eps.val() + pow(jac(G).det(), 2.0), 0.5));
        // auto invJacMat = jac(G).adj()/chi;
        // auto eta = m_evaluator.getVariable(m_fun);
        // return m_evaluator.integral( (monitor(eta,G).asDiag()*invJacMat).sqNorm()*meas(G));

        // gsComposedFunction<T> fun(*m_comp,m_fun.function(0));
        // gsComposedFunction<T> fun(mpG.patch(0),m_fun.function(0));


        // gsVector<unsigned> np(m_fun.domainDim());
        // np.setConstant(m_options.getInt("nSamplingPoints"));
        // gsMatrix<T> grid = gsPointGrid<T>(m_cgeom.support().col(0),m_cgeom.support().col(1),np);
        if (m_cgeom.domainDim()==m_cgeom.targetDim())
        {
            auto jacG = jac(G).det();
            auto chi = 0.5 * (jacG + pow(pow(eps.val(),2.0) + pow(jacG, 2.0), 0.5));
            auto invJacMat = jac(G).adj()/chi; // inverse of jacobian matrix with 'determinant' replaced

            // auto eta = m_evaluator.getVariable(fun);
            auto eta = m_evaluator.getVariable(m_cfun);
            // auto eta = m_evaluator.getVariable(m_fun,G);

            auto M   = gismo::expr::monitor<MODE>(eta,G,m_options.getReal("Smoothing"));
            return m_evaluator.integral( (M*invJacMat).sqNorm()*meas(G));
        }
        else
        {
            // auto fform = jac(G).tr()*jac(G);
            // auto jacG = sqrt(fform.det()); //jacobian determinant for a surface, i.e. the measure
            // // auto chi = 0.5 * (jacG + pow(pow(eps.val(),2.0) + pow(jacG, 2.0), 0.5));

            // // Compute the chi part
            // auto chiPPart = eps * ((jacG - eps.val()).exp());

            // // Ternary operation to compute chi and chip
            // auto chi = ternary(eps.val() - jacG, chiPPart.val(), jacG.val());

            // auto invJacMat = fform.adj()/chi; // inverse of jacobian matrix with 'determinant' replaced
            // auto eta = m_evaluator.getVariable(fun);
            // return m_evaluator.integral( (monitor(eta,G)*invJacMat).sqNorm()*meas(G));

            GISMO_ERROR("The dimension of target domain should be 2 or 3.");
            return 0;
        }
    }

    // /// Computes the gradient of the objective function at the given point u
    // // and stores it in result.
    // void gradObj_into(const gsAsConstVector<T> &u,
    //                 gsAsVector<T> &result) const override;

protected:

    // From gsOptProblem
    using Base::m_curDesign;
    using Base::m_numDesignVars;

    // Controls of the composition
    // NOTE: Different from m_curDesign, since m_curDesign updates every time
    gsMatrix<T> m_controls;

    gsFunction<T>         & m_comp;
    gsComposedGeometry<T>   m_cgeom;
    gsComposedFunction<T>   m_cfun;

    gsOptionList m_options;

    mutable gsExprEvaluator<T> m_evaluator;
};

/*
    NOTES:
    * Give the integration basis?? The integration elements should be the finest of the elements of the composition and the basis/geometry
    * There should be a rule of thumb for the number of integration points. Given a composition of degree p and a geometry of degree q, the number of integration points should be at least p*q+1?? or p+q+1??
 */
template<class T, enum MonitorMode MODE=MonitorMode::ValueBased>
class gsComposedRelocation
{
protected:

public:

    /**
     * @brief      Constructs a new instance.
     * @param[in]  composition  The composition to apply to the geometry
     * @param[in]  geometry     The geometry to be relocated
     * @param[in]  fun          The function used for the monitor function
     */
    gsComposedRelocation(      gsFunction<T>  & composition,
                         const gsGeometry<T>  & geometry,
                         const gsFunction<T>  & function,
                               gsOptimizer<T> & optimizer,
                         const bool               parametric=true)
    :
    m_optProblem(composition,geometry,function,parametric),
    m_optimizer(optimizer)
    {
        this->defaultOptions();
    }

    gsOptionList & options() { return m_options; }

    void defaultOptions()
    {
        m_options.addReal("Penalty","Penalization coefficient for Jacobian determinant [default=0.01]",1e-2);
        m_options.addReal("Smoothing","Smoothing parameter in the monitor function [default=0.1]",0.1);
        m_options.addInt("Mode","0: Relocate based on f [default]; 1: Relocate based on grad(f)",0);

        // Integration basis options [HUGO: replace with numGauss]
        m_options.addInt("nRefine","Number of refinement steps for the integration basis",2);
        m_options.addInt("nElevate","Number of elevation steps for the integration basis",1);
    }

    void solve()
    {
        // Set the optimization problem
        m_optProblem.options().setInt("nRefine",m_options.getInt("nRefine"));
        m_optProblem.options().setInt("nElevate",m_options.getInt("nElevate"));
        m_optProblem.options().setReal("Smoothing",m_options.getReal("Smoothing"));
        m_optProblem.options().setReal("Penalty",m_options.getReal("Penalty"));

        m_optimizer.setProblem(&m_optProblem);
        // Solve the optimization problem
        gsVector<T> controls(m_optProblem.composition().nControls());
        for (index_t k=0; k!=m_optProblem.composition().nControls(); k++)
            controls[k] = m_optProblem.composition().control(k);

        gsDebugVar(controls.transpose());
        m_optimizer.solve(controls);
        controls = m_optimizer.currentDesign();
        gsDebugVar(controls.transpose());
        // Write the optimized solution back to the composition
        for (index_t k=0; k!=controls.rows(); k++)
            m_optProblem.composition().control(k) = controls[k];
    }

protected:

    gsOptMesh<T,MODE>       m_optProblem;
    gsOptimizer<T>        & m_optimizer;

    gsOptionList m_options;
};


int main(int arg, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 2;
    index_t numElevate = 1;
    index_t numRefineI = 2;
    index_t numElevateI = 1;
    index_t maxIt = 100;
    real_t tol_g = 5e-5;
    real_t eps = 1e-2;
    bool slide = true;
    index_t testCase = 0;
    index_t opt = 2;
    index_t nSamplingPoints = 50;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "elevAnalysis","Number of degree elevation steps to perform for the analysis", numElevate );
    cmd.addInt( "r", "refAnalysis", "Number of Uniform h-refinement loops for the analysis",  numRefine );
    cmd.addInt( "E", "elevIntegral","Number of degree elevation steps to perform for the integration", numElevateI );
    cmd.addInt( "R", "refIntegral", "Number of Uniform h-refinement loops for the integration",  numRefineI );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addReal("g", "tolG", "relative tol", tol_g);
    cmd.addInt( "i", "maxIt", "max num iterations",  maxIt );
    cmd.addReal("", "eps", "eps",  eps );
    cmd.addSwitch("noslide", "Do not slide the boundaries",  slide );
    cmd.addInt( "t", "testCase", "Function to be used: 0: cosine waves, 1: spiral.",  testCase );
    cmd.addInt( "o", "opt", "Optimizer: 0: gsGradientDescent, 1: gsHLBFGS, 2: gsOptim::LBFGS.",  opt );
    cmd.addInt( "S", "nSamplingPoints", "Number of sampling points in each parametric direction",  nSamplingPoints );

    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // std::string dirname = "template_r-adaptivity_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine)+"_S"+util::to_string(nSamplingPoints)+"_function"+util::to_string(testCase)+"_opt"+util::to_string(opt);
    std::string dirname = "template_r-adaptivity_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine)+"_E"+util::to_string(numElevateI)+"_R"+util::to_string(numRefineI)+"_function"+util::to_string(testCase)+"_opt"+util::to_string(opt);
    gsFileManager::mkdir(dirname);
    dirname += gsFileManager::getNativePathSeparator();

    gsMultiPatch<> geom;
    // geom.addPatch(*gsNurbsCreator<>::BSplineSquare());
    geom.addPatch(*gsNurbsCreator<>::BSplineRectangle(0,0,2,1));
    geom.patch(0).uniformRefine(1,1,0);

    // Basis for the square domain
    gsKnotVector<> kv({0,0,1,1},1);
    gsTensorBSplineBasis<2> tbbasis(kv,kv);
    tbbasis.degreeElevate(numElevate);
    for (index_t i = 0; i < numRefine; i++)
        tbbasis.uniformRefine();

    tbbasis.uniformRefine(1,1,0);


    gsInfo<<"Mapper basis:\n"<<tbbasis<<"\n";

    gsSquareDomain<2,real_t> domain(tbbasis);
    domain.options().addSwitch("Slide","",slide);
    domain.applyOptions();
    // domain.perturb(1e-1);

    gsComposedGeometry<real_t> cspline(domain,geom.patch(0));
    gsFunction<> * fun;
    if      (testCase==0)
    {
        fun = new gsFunctionExpr<>("1 + 5 * exp( -50 * abs( (x-0.5)^2 + (y-0.5)^2 - 0.09 ) ) ",2);
    }
    else if (testCase==1)
    {
        std::string R     = "sqrt( (x-0.5)^2 + (y-0.5)^2 )";
        fun = new gsFunctionExpr<>("1 /(2 + cos( 8 * pi * "+R+"))",2);
    }
    else if (testCase==2)
    {
        std::string R     = "sqrt( (x-0.7)^2 + (y-0.5)^2 )";
        std::string Theta = "atan2((y-0.5),(x-0.7))";
        fun = new gsFunctionExpr<>("1 + 9/(1 + ( 10 * "+R+" * cos(" + Theta +" - 20 * "+R+"^2 ) )^2)",2);
    }
    else
    {
        GISMO_ERROR("Unknown test case");
    }

    gsComposedFunction<real_t> cfun({&domain,fun});
    // gsComposedFunction<real_t> cfun({&cspline,fun});
    gsComposedGeometry<real_t> cgeom(domain,geom.patch(0));

/*
    PERFORM R-ADAPTIVITY
 */

    gsInfo<<"Number of optimizer degrees of freedom: "<<domain.nControls()<<"\n";

    gsOptimizer<real_t> * optimizer;
    if      (opt==0) // gsGradientDescent
    {
        optimizer = new gsGradientDescent<real_t>;
        optimizer->options().setInt("MaxIterations",maxIt);
        optimizer->options().setInt("Verbose",2);
        optimizer->options().setReal("MinGradientLength",tol_g);

    }
    else if (opt==1) // gsHLBFGS
    {
        optimizer = new gsHLBFGS<real_t>;
        optimizer->options().setInt("MaxIterations",maxIt);
        optimizer->options().setInt("Verbose",2);
        optimizer->options().setReal("tolRelG",tol_g);
    }
    else if (opt==2) //gsOptim::LBFGS
    {
        optimizer = new gsOptim<real_t>::LBFGS;
        optimizer->options().setInt("MaxIterations",maxIt);
        optimizer->options().setInt("Verbose",1);
        optimizer->options().setReal("GradErrTol",tol_g);
    }
    else
    {
        GISMO_ERROR("Unknown optimizer");
    }

    for (index_t k=0; k!=domain.nControls(); k++)
        gsInfo<<domain.control(k)<<" ";
    gsInfo<<std::endl;
    // gsComposedRelocation<real_t,MonitorMode::ValueBased> relocator(domain,geom.patch(0),cfun,*optimizer);
    gsComposedRelocation<real_t,MonitorMode::ValueBased> relocator(domain,geom.patch(0),*fun,*optimizer);
    relocator.options().setInt("nRefine",numRefineI);
    relocator.options().setInt("nElevate",numElevateI);
    relocator.solve();
    for (index_t k=0; k!=domain.nControls(); k++)
        gsInfo<<domain.control(k)<<" ";
    gsInfo<<std::endl;

    gsMultiPatch<> mp;
    mp.addPatch(cspline);
    // mp.embed(3);
    gsMultiBasis<> mb(mp);

    gsExprEvaluator<> ev;
    ev.setIntegrationElements(mb);
    auto G = ev.getMap(mp);
    auto f = ev.getVariable(*fun);

    gsField<> cfun_field(mp,cfun);

    gsWriteParaview(mp,*fun,dirname+"fun");
    gsWriteParaview(mp,cfun,dirname+"cfun");
    gsWriteParaview(cfun_field,dirname+"cfun_field");

    // ev.writeParaview(f,G,dirname+"fun");

    // gsWriteParaview(cspline,"cspline",1000,true);
    // gsWriteParaview(cspline.basis(),"cbasis",1000);
    gsWriteParaview(domain.domain(),dirname+"domain",1000,true,true);

    // gsInfo<<"Area = "<<ev.integral(meas(G))<<"\n";
    ev.writeParaview(jac(G).det(),G,dirname+"jacobian_determinant");

    delete fun;
    delete optimizer;
    return 0;
}// end main
