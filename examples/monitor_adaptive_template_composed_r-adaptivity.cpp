/** @file monitor_adaptive_template_composed_r-adaptivity.cpp

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

namespace gismo{
namespace expr{

template<class E>
class monitor_expr : public _expr<monitor_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;
    mutable gsMatrix<Scalar> grad, jac, jacInv;
    mutable gsMatrix<Scalar> ones;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _G;
    const short_t DIM;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    monitor_expr(const E & u, const gsGeometryMap<Scalar> & G) 
    : 
    _u(u), 
    _G(G),
    DIM(_u.source().domainDim())
    {
    }

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {return eval_impl(_u,k); }

    index_t rows() const { return DIM; }

    index_t cols() const { return DIM; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        // jac_expr<gsFeVariable<Scalar>>(_u).parse(evList);
        jacInv_expr<Scalar>(_G).parse(evList);
        // _u.parse(evList);
        _G.parse(evList);

        // OLD
        evList.add(_u);
        _u.data().flags |= NEED_JACOBIAN;
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
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
    }

    template<class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        real_t eps = 0.1; // determines the degree of smoothing
        

        jac = _u.data().values[1].col(k).transpose();

        jacInv_expr<Scalar> jacInvExpr = jacInv_expr<Scalar>(_G);
        jacInv = jacInvExpr.eval(k);

        grad = jac * jacInv;
        ones.resize(DIM,DIM);
        ones.setIdentity();

        res = 1.0 / ( math::sqrt(1.0 +eps*grad.squaredNorm() ) ) * ones;
        return res;

    }
};

template<class E> EIGEN_STRONG_INLINE
monitor_expr<E> monitor(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return monitor_expr<E>(u,G); }

}
}

using namespace gismo;



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

        //m_options.addInt("nSamplingPoints","Number of sampling points in each parametric direction",50);
    
        m_options.addInt("nRefine","Number of refinement steps for the integration basis",2);
        m_options.addInt("nElevate","Number of elevation steps for the integration basis",1);
    }

    gsOptionList & options() { return m_options; }

    /// Evaluates the objective function at the given point u.
    T evalObj(const gsAsConstVector<T> &u) const override
    {
        for (index_t k=0; k!=u.rows(); k++)
            m_comp->control(k) = u[k];

        gsComposedGeometry<T> cgeom(*m_comp,m_geom);

        gsMultiPatch<> mpG;
        mpG.addPatch(cgeom);

        gsKnotVector<> kv({0,0,1,1},1);
        gsTensorBSplineBasis<2,T> basis(kv,kv);
        basis.degreeElevate(m_options.getInt("nElevate"));

        for (index_t i = 0; i < m_options.getInt("nRefine"); i++)
            basis.uniformRefine();

        gsMultiBasis<> mb(basis);

        m_evaluator.setIntegrationElements(mb);
        geometryMap G = m_evaluator.getMap(mpG);

        gsConstantFunction<T> epsilonFunction(m_eps, m_geom.domainDim());
        auto eps = m_evaluator.getVariable(epsilonFunction);

        // auto chi = 0.5 * (jac(G).det() + pow(eps.val() + pow(jac(G).det(), 2.0), 0.5));
        // auto invJacMat = jac(G).adj()/chi;
        // auto eta = m_evaluator.getVariable(m_fun);
        // return m_evaluator.integral( (monitor(eta,G).asDiag()*invJacMat).sqNorm()*meas(G));

        gsComposedFunction<T> fun(*m_comp,m_fun.function(0));

        // gsVector<unsigned> np(m_fun.domainDim());
        // np.setConstant(m_options.getInt("nSamplingPoints"));
        // gsMatrix<T> grid = gsPointGrid<T>(m_geom.support().col(0),m_geom.support().col(1),np);
        if (m_geom.domainDim()==m_geom.targetDim())
        {
            auto jacG = jac(G).det();
            auto chi = 0.5 * (jacG + pow(eps.val() + pow(jacG, 2.0), 0.5));
            auto invJacMat = jac(G).adj()/chi; // inverse of jacobian matrix with 'determinant' replaced
            auto eta = m_evaluator.getVariable(fun);    
            // return m_evaluator.eval((monitor(eta,G)*invJacMat).sqNorm(),grid).sum();
            return m_evaluator.integral( (monitor(eta,G)*invJacMat).sqNorm()*meas(G));
        }
        else
        {
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


    gsFunction<T> * m_comp;
    const gsGeometry<T> & m_geom;
    const gsFunctionSet<T> & m_fun;
    T m_eps;

    gsOptionList m_options;

    mutable gsExprEvaluator<T> m_evaluator;
};


int main(int arg, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numURef  = 3;
    index_t numRefine  = 2;
    index_t numElevate = 1;
    index_t numRefineI = 5;
    index_t numElevateI = 3;
    index_t maxIt = 100;
    real_t tol_g = 5e-5;
    real_t eps = 1e-4;
    bool slide = true;
    index_t testCase = 0;
    index_t opt = 2;
    index_t numAdaptiveLoops = 2;
    index_t rule = 1;
    real_t perc_ref = 0.1;
    

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
    cmd.addInt( "j", "ell", "Maximum number of iterations for the adaptive loop.",  numAdaptiveLoops );
    cmd.addInt( "w", "rule", "Refinement rule: 1 - GARU; 2 - PUKA; 3 - BULK.",  rule );
    cmd.addReal("p", "refp", "Refine percentage",  perc_ref );
    

    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // std::string dirname = "template_r-adaptivity_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine)+"_S"+util::to_string(nSamplingPoints)+"_function"+util::to_string(testCase)+"_opt"+util::to_string(opt);
    std::string dirname = "template_r-adaptivity_e"+util::to_string(numElevate)+"_r"+util::to_string(numRefine)+"_E"+util::to_string(numElevateI)+"_R"+util::to_string(numRefineI)+"_function"+util::to_string(testCase)+"_opt"+util::to_string(opt);
    gsFileManager::mkdir(dirname);
    dirname += gsFileManager::getNativePathSeparator();


    // Geometry (square domain, represented by thb.)
    gsTensorBSpline<2,real_t> nurbs = *gsNurbsCreator<>::BSplineSquare(1,0,0);
    nurbs.uniformRefine( (1<<numURef)-1 );
    gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&nurbs);
    gsTHBSpline<2,real_t> thb = gsTHBSpline<2,real_t>(*geo);
    
    gsMultiPatch<> geom;
    geom.addPatch(thb);
    gsWriteParaview(thb, dirname+"init_geom", 10000, true, true);



    // Basis for the sigma mapper
    gsKnotVector<> kv({0,0,1,1},1);
    gsTensorBSplineBasis<2> tbbasis(kv,kv);
    tbbasis.degreeElevate(numElevate);
    for (index_t i = 0; i < numRefine; i++)
        tbbasis.uniformRefine();
    
    gsInfo<<"Mapper basis:\n"<<tbbasis<<"\n";

    // // sigma mapper
    // gsSquareDomain<2,real_t> domain(tbbasis);
    // domain.options().addSwitch("Slide","",slide);
    // domain.applyOptions();


    // error field to test
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


    /* 
        PERFORM H-ADAPTIVITY
    */
    //! [beginRefLoop]
    for( int refLoop = 0; refLoop <= numAdaptiveLoops; refLoop++)
    {
        gsInfo << "#loop = "<< refLoop + 1<< "\n";
    /* 
        PERFORM R-ADAPTIVITY
    */

    // sigma mapper
    gsSquareDomain<2,real_t> domain(tbbasis);
    domain.options().addSwitch("Slide","",slide);
    domain.applyOptions();


    gsInfo<<"Number of optimizer degrees of freedom: "<<domain.nControls()<<"\n";

    gsWriteParaview(geom,dirname+"geom_in"+std::to_string(refLoop),1,true);
    
    gsMultiPatch<> mp;
    gsComposedGeometry<real_t> cspline(domain,geom.patch(0));
    mp.addPatch(cspline);
    gsWriteParaview(geom,dirname+"cspline_in"+std::to_string(refLoop),1,true);


    gsOptMesh<> optMesh(domain,geom.patch(0),*fun,eps);
    gsVector<> controls(domain.nControls());
    
    optMesh.options().setInt("nRefine",numRefineI);
    optMesh.options().setInt("nElevate",numElevateI);
    for (size_t k=0; k!=domain.nControls(); k++)
        controls[k] = domain.control(k);


    gsOptimizer<real_t> * optimizer;
    if      (opt==0) // gsGradientDescent
    {
        optimizer = new gsGradientDescent<real_t>(&optMesh);
        optimizer->options().setInt("MaxIterations",maxIt);
        optimizer->options().setInt("Verbose",2);
        optimizer->options().setReal("MinGradientLength",tol_g);

    }
    else if (opt==1) // gsHLBFGS
    {
        optimizer = new gsHLBFGS<real_t>(&optMesh);
        optimizer->options().setInt("MaxIterations",maxIt);
        optimizer->options().setInt("Verbose",2);
        optimizer->options().setReal("tolRelG",tol_g);
    }
    else if (opt==2) //gsOptim::LBFGS
    {
        optimizer = new gsOptim<real_t>::LBFGS(&optMesh);
        optimizer->options().setInt("MaxIterations",maxIt);
        optimizer->options().setInt("Verbose",1);
        optimizer->options().setReal("GradErrTol",tol_g);
    }
    else
    {
        GISMO_ERROR("Unknown optimizer");
    }

    optimizer->solve(controls);
    gsVector<> optSol = optimizer->currentDesign();
    
    for (size_t k=0; k!=optSol.rows(); k++)
        domain.control(k) = optSol[k];

    // gsInfo << "Domain after optimization:\n" << domain.domain() << "\n";
    // gsInfo << "Domain #coefficients after optimization:\n" << domain.domain().coefs().rows() << "\n";
    // gsInfo << "Domain #controls after optimization:\n" << domain.nControls() << "\n";

    
    // mp.embed(3);

    gsWriteParaview(cspline,dirname+"cspline_r"+std::to_string(refLoop),1000,true,true); // plots in xi, eta because it uses the support of the basis of cspline (thus the support of the cbasis, thus the support of the original basis)
    //gsWriteParaview(cspline.basis(),dirname+"cbasis_init"+std::to_string(refLoop),1000);

    gsMultiBasis<> mb(geom);

    gsExprEvaluator<> ev;
    ev.setIntegrationElements(mb);
    auto G = ev.getMap(mp); // geometry map
    auto f  = ev.getVariable(*fun,G);

    ev.writeParaview(f,G,dirname+"fun");
   
    gsWriteParaview(domain.domain(),dirname+"domain"+std::to_string(refLoop),1000,true,true);
    ev.writeParaview(jac(G).det(),G,dirname+"jacobian_determinant"+std::to_string(refLoop));
    

    // Get the element-wise errors.
    //ev.integralElWise( f*meas(G) );
    ev.maxElWise(f);

    
    const std::vector<real_t> & elErrs = ev.elementwise();
    gsInfo << elErrs.size() << " element-wise errors\n";

    // to check
    gsElementErrorPlotter<real_t> err_eh(thb.basis(),elErrs);
    // gsComposedFunction<real_t> cerr_eh(domain,err_eh);
    // gsWriteParaview<>(err_eh, cspline.support(), dirname+"error_elem_ref", 10000); 
    gsWriteParaview<>(cspline,err_eh, dirname+"error_elem_ref"+std::to_string(refLoop), 1000); 

    gsAdaptiveMeshing<real_t> mesher(geom);
    mesher.options().setInt("RefineRule",rule);
    mesher.options().setInt("CoarsenRule",rule);
    mesher.options().setSwitch("Admissible",true);
    mesher.options().setReal("RefineParam",perc_ref);
    
    
    mesher.getOptions();
    gsHBoxContainer<2,real_t> refine;
    mesher.markRef_into(elErrs,refine);

    // gsInfo<<"Cells marked for refinement:\n";
    // gsInfo<<refine<<"\n";
    gsWriteParaview(refine,dirname+"marked4ref_"+std::to_string(refLoop));


    mesher.refine(refine);
    gsWriteParaview(geom,dirname+"geom_ref"+std::to_string(refLoop),1,true);

    gsComposedGeometry<real_t> refcspline(domain,geom.patch(0));
    gsWriteParaview(refcspline,dirname+"cspline_ref"+std::to_string(refLoop),1000,true,true);

    delete optimizer;

    }//! [endRefLoop]

    delete fun;
    return 0;
}// end main
