/** @file gsOptim.h

    @brief Provides declaration of an optimization problem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsCore/gsLinearAlgebra.h>
#include <gsOptimizer/gsOptProblem.h>
#include <gsOptimizer/gsOptimizer.h>

#   define Eigen gsEigen
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include <optim/header_only_version/optim.hpp>
#   undef Eigen

using namespace optim;

#pragma once

namespace gismo
{

// Forward declarations
template <typename T> class gsOptimBFGS;
template <typename T> class gsOptimLBFGS;
template <typename T> class gsOptimCG;
template <typename T> class gsOptimGD;
template <typename T> class gsOptimNewton;
template <typename T> class gsOptimNM;
template <typename T> class gsOptimDE;
template <typename T> class gsOptimDEPRMM;
template <typename T> class gsOptimPSO;
template <typename T> class gsOptimPSODV;
template <typename T> class gsOptimSUMT;

/*
 * TODO:
 * - What to do with root finding
 */

template <typename T>
struct gsOptimWrapper
{
    gsOptimWrapper(gsOptProblem<T> & problem)
    :
    m_op(problem)
    { }

    T operator() (const optim::ColVec_t& vals_in, optim::ColVec_t* grad_out, void* opt_data)
    {
        gsAsConstVector<T> in(vals_in.data(),vals_in.rows());
        if (grad_out)
        {
            grad_out->resize(vals_in.size());
            gsAsVector<T> out(grad_out->data(),grad_out->size());
            m_op.gradObj_into(in, out);
        }
        return m_op.evalObj(in);
    }

private:
    gsOptProblem<T> & m_op;
};

template <typename T>
struct gsOptimWrapperConstraint
{
    gsOptimWrapperConstraint(gsOptProblem<T> & problem)
    :
    m_op(problem)
    { }

    optim::ColVec_t operator() (const optim::ColVec_t& vals_in, optim::Mat_t* jacob_out, void* constr_data)
    {
        gsAsConstVector<T> in(vals_in.data(),vals_in.size());
        if (jacob_out)
        {
            std::vector<T> jacTmp(vals_in.size()*m_op.numConstraints());
            gsAsVector<T> jacOut(jacTmp);
            m_op.jacobCon_into(in, jacOut);

            /// WRITE TO MATRIX
            // The matrix is organized as follows (see gsOptProblem.h, computeJacStructure):
            // - row i, constraint i
            // - col j, design variable j
            jacob_out->resize(m_op.numConstraints(),vals_in.size());
            gsAsMatrix<T> out(jacob_out->data(),jacob_out->rows(),jacob_out->cols());
            out.setZero();
            std::vector<index_t> conJacRows = m_op.conJacRows();
            std::vector<index_t> conJacCols = m_op.conJacCols();
            for (size_t k=0; k!=conJacRows.size(); k++)
                out(conJacRows[k],conJacCols[k]) = jacTmp[k];
        }
        gsVector<T> conTmp(m_op.numConstraints());
        gsAsVector<T> conOut(conTmp.data(),conTmp.size());
        m_op.evalCon_into(in,conOut);
        return conTmp;
    }

private:
    gsOptProblem<T> & m_op;
};


template <typename T>
class gsOptim : public gsOptimizer<T>
{
    using Base = gsOptimizer<T>;
public:

    typedef memory::unique_ptr<gsOptim> uPtr;

    typedef gsOptimBFGS<T>              BFGS;
    typedef gsOptimLBFGS<T>             LBFGS;
    typedef gsOptimCG<T>                CG;
    typedef gsOptimGD<T>                GD;
    // typedef gsOptimNewton<T>            Newton;
    typedef gsOptimNM<T>                NM;
    typedef gsOptimDE<T>                DE;
    typedef gsOptimDEPRMM<T>            DEPRMM;
    typedef gsOptimPSO<T>               PSO;
    typedef gsOptimPSODV<T>             PSODV;
    typedef gsOptimSUMT<T>              SUMT;

public:

    gsOptim() {};

    gsOptim(gsOptProblem<T> * problem)
    :
    Base(problem),
    m_success(false)
    {
        this->defaultOptions();
    }

    static uPtr get(const std::string & slv, gsOptProblem<T> * problem)
    {
        if (slv=="BFGS") return uPtr(new BFGS(problem));
        if (slv=="LBFGS") return uPtr(new LBFGS(problem));
        if (slv=="CG") return uPtr(new CG(problem));
        if (slv=="GD") return uPtr(new GD(problem));
        // if (slv=="Newton") return uPtr(new Newton(problem));
        if (slv=="NM") return uPtr(new NM(problem));
        if (slv=="DE") return uPtr(new DE(problem));
        if (slv=="DEPRMM") return uPtr(new DEPRMM(problem));
        if (slv=="PSO") return uPtr(new PSO(problem));
        if (slv=="PSODV") return uPtr(new PSODV(problem));
        if (slv=="SUMT") return uPtr(new SUMT(problem));
        GISMO_ERROR("Solver \'"<< slv << "\' not known to G+Smo");
    }

    virtual void defaultOptions()
    {
        Base::defaultOptions();

        m_options.addInt("ConfFailureSwitch","Switch for convergence failure",0);
        m_options.addReal("GradErrTol","Gradient error tolerance (default: 1e-8)",1e-8);
        m_options.addReal("RelSolChangeTol","Relative solution change tolerance (default: 1e-14)",1e-14);
        m_options.addReal("RelObjFnChangeTol","Relative objective function change tolerance (default: 1e-8)",1e-8);


        // m_options.addReal("MinGradientLength","Minimal gradient length",1e-9);
        // m_options.addReal("MinStepLength","Minimal step length",1e-9);
        // m_options.addInt("LBFGSUpdates","Number of LBFGS updates (typically 3-20, put to 0 for gradient descent)",20);
    }

    virtual void getOptions()
    {
        Base::getOptions();

        // See: https://optimlib.readthedocs.io/en/latest/settings.html?highlight=iterations#main-settings

        // RNG seeding
        // size_t rng_seed_value = std::random_device{}();

        // print and convergence options
        m_optimSettings.print_level = m_verbose;
        m_optimSettings.conv_failure_switch = m_options.getInt("ConfFailureSwitch");

        // error tolerance and maxiumum iterations
        m_optimSettings.iter_max = m_maxIterations;

        m_optimSettings.grad_err_tol = m_options.getReal("GradErrTol");
        m_optimSettings.rel_sol_change_tol  = m_options.getReal("RelSolChangeTol");
        m_optimSettings.rel_objfn_change_tol = m_options.getReal("RelObjFnChangeTol");

        // bounds
        // See solve function
    }

    virtual void solve(const gsMatrix<T> & initialGuess)
    {
        // Get the bounds
        this->getOptions();
        this->setConstraints();
        GISMO_ASSERT(initialGuess.cols()==1,"The initial guess should have vector format");
        gsVector<T> x = initialGuess.col(0);
        m_success = callOptim(x, *m_op, m_optimSettings);
        m_curDesign = x;
        m_numIterations = m_optimSettings.opt_iter;
        m_finalObjective = m_optimSettings.opt_fn_value;
    }

    bool success() { return m_success; }

protected:

    virtual bool callOptim(gsVector<T> & initialGuess, gsOptProblem<T> & op, optim::algo_settings_t & settings) = 0;

    void setConstraints()
    {
        if(m_op->desLowerBounds().size()!=0)
        {
            m_optimSettings.lower_bounds = m_op->desLowerBounds();
            m_optimSettings.upper_bounds = m_op->desUpperBounds();
            m_optimSettings.vals_bound = true;
        }
    }

// Members taken from gsOptimizer
protected:
    using Base::m_op;
    using Base::m_numIterations;
    using Base::m_finalObjective;
    using Base::m_curDesign;
    using Base::m_options;
    using Base::m_verbose;
    using Base::m_maxIterations;

    using Base::defaultOptions;
    using Base::getOptions;

protected:
    // Optim options
    optim::algo_settings_t m_optimSettings;
    bool m_success;
};

template <typename T>
class gsOptimBFGS : public gsOptim<T>
{
public:
    using Base = gsOptim<T>;

public:

    gsOptimBFGS(gsOptProblem<T> * problem) : Base(problem)
    {
        this->defaultOptions();
    }

    bool callOptim( gsVector<T> & x,
                    gsOptProblem<T> & op,
                    optim::algo_settings_t & optimSettings)
    { return optim::bfgs(x, gsOptimWrapper<T>(op), nullptr,optimSettings); }

    void defaultOptions() override
    {
        Base::defaultOptions();
        m_options.addReal("WolfeCons1","Line search tuning parameter that controls the tolerance on the Armijo sufficient decrease condition.",1e-3);
        m_options.addReal("WolfeCons2","Line search tuning parameter that controls the tolerance on the curvature condition.",0.9);
    }

    void getOptions() override
    {
        Base::getOptions();
        m_optimSettings.bfgs_settings.wolfe_cons_1 = m_options.getReal("WolfeCons1");
        m_optimSettings.bfgs_settings.wolfe_cons_2 = m_options.getReal("WolfeCons2");
    }

protected:

    using Base::m_options;
    using Base::m_optimSettings;
};

template <typename T>
class gsOptimLBFGS : public gsOptim<T>
{
public:
    using Base = gsOptim<T>;

public:

    gsOptimLBFGS(gsOptProblem<T> * problem) : Base(problem)
    {
        this->defaultOptions();
    }

    bool callOptim( gsVector<T> & x,
                    gsOptProblem<T> & op,
                    optim::algo_settings_t & optimSettings)
    { return optim::lbfgs(x, gsOptimWrapper<T>(op), nullptr,optimSettings); }

    void defaultOptions() override
    {
        Base::defaultOptions();
        m_options.addInt("ParM","The number of past gradient vectors to use when forming the approximate Hessian matrix.",10);
        m_options.addReal("WolfeCons1","Line search tuning parameter that controls the tolerance on the Armijo sufficient decrease condition.",1e-3);
        m_options.addReal("WolfeCons2","Line search tuning parameter that controls the tolerance on the curvature condition.",0.9);
    }

    void getOptions() override
    {
        Base::getOptions();
        m_optimSettings.lbfgs_settings.par_M = m_options.getInt("ParM");
        m_optimSettings.lbfgs_settings.wolfe_cons_1 = m_options.getReal("WolfeCons1");
        m_optimSettings.lbfgs_settings.wolfe_cons_2 = m_options.getReal("WolfeCons2");
    }

protected:

    using Base::m_options;
    using Base::m_optimSettings;
};

template <typename T>
class gsOptimCG : public gsOptim<T>
{
public:
    using Base = gsOptim<T>;

public:

    gsOptimCG(gsOptProblem<T> * problem) : Base(problem)
    {
        this->defaultOptions();
    }

    bool callOptim( gsVector<T> & x,
                    gsOptProblem<T> & op,
                    optim::algo_settings_t & optimSettings)
    { return optim::cg(x, gsOptimWrapper<T>(op), nullptr,optimSettings); }

    void defaultOptions() override
    {
        Base::defaultOptions();
        m_options.addInt("UpdateMethod","Update method: 1) Fletcher–Reeves, 2) Polak-Ribiere, 3) FR-PR Hybrid, 4) Hestenes-Stiefel, 5) Dai-Yuan, 6) Hager-Zhang.",2);
        m_options.addReal("RestartThreshold","parameter ν from step 2 in the algorithm description.",0.1);
        m_options.addSwitch("UseRelSolChangeCrit","whether to enable the `rel_sol_change_tol` stopping criterion.",false);
        m_options.addReal("WolfeCons1","Line search tuning parameter that controls the tolerance on the Armijo sufficient decrease condition.",1e-3);
        m_options.addReal("WolfeCons2","Line search tuning parameter that controls the tolerance on the curvature condition.",0.9);
    }

    void getOptions() override
    {
        Base::getOptions();
        m_optimSettings.cg_settings.method = m_options.getInt("UpdateMethod");
        m_optimSettings.cg_settings.restart_threshold = m_options.getReal("RestartThreshold");
        m_optimSettings.cg_settings.use_rel_sol_change_crit = m_options.getSwitch("UseRelSolChangeCrit");
        m_optimSettings.cg_settings.wolfe_cons_1 = m_options.getReal("WolfeCons1");
        m_optimSettings.cg_settings.wolfe_cons_2 = m_options.getReal("WolfeCons2");
    }

protected:

    using Base::m_options;
    using Base::m_optimSettings;
};

template <typename T>
class gsOptimGD : public gsOptim<T>
{
public:
    using Base = gsOptim<T>;

public:

    gsOptimGD(gsOptProblem<T> * problem) : Base(problem)
    {
        this->defaultOptions();
    }

    bool callOptim( gsVector<T> & x,
                    gsOptProblem<T> & op,
                    optim::algo_settings_t & optimSettings)
    { return optim::gd(x, gsOptimWrapper<T>(op), nullptr,optimSettings); }

    void defaultOptions() override
    {
        Base::defaultOptions();

        m_options.addInt("Method","Method to use: 0) Vanilla GD, 1) GD with momentum, 2) Nesterov accelerated gradient descent, 3) AdaGrad. 4) RMSProp, 5) AdaDelta. 6) Adam (adaptive moment estimation) and AdaMax, 7) Nadam (adaptive moment estimation) and NadaMax",0);
        m_options.addReal("StepSize","The learning rate.",0.1);
        m_options.addSwitch("StepDecay","Whether to use step decay.",false);
        m_options.addInt("StepDecayPeriods","Number of periods after which to decay the step size.",10);
        m_options.addReal("StepDecayValue","The value by which to decay the step size.",0.5);
        m_options.addReal("Momentum","The momentum parameter.",0.9);
        m_options.addReal("AdaNormTerm","The normalization term for AdaGrad.",1.0e-08);
        m_options.addReal("AdaRho","The rho parameter for AdaGrad.",0.9);
        m_options.addSwitch("AdaMax","Whether to use AdaMax.",false);
        m_options.addReal("AdamBeta1","The beta1 parameter for Adam.",0.9);
        m_options.addReal("AdamBeta2","The beta2 parameter for Adam.",0.999);
        m_options.addSwitch("ClipGrad","Whether to clip the gradient.",false);
        m_options.addSwitch("ClipMaxNorm","Whether to clip the gradient norm.",false);
        m_options.addSwitch("ClipMinNorm","Whether to clip the gradient norm.",false);
        m_options.addInt("ClipNormType","The type of norm to use for clipping.",2);
        m_options.addReal("ClipNormBound","The bound for the norm used for clipping.",5.0);
    }

    void getOptions() override
    {
        Base::getOptions();

        m_optimSettings.gd_settings.method = m_options.getInt("Method");

        // step size, or 'the learning rate'
        m_optimSettings.gd_settings.par_step_size = m_options.getReal("StepSize");

        // decay
        m_optimSettings.gd_settings.step_decay = m_options.getSwitch("StepDecay");
        m_optimSettings.gd_settings.step_decay_periods = m_options.getInt("StepDecayPeriods");
        m_optimSettings.gd_settings.step_decay_val = m_options.getReal("StepDecayValue");

        // momentum parameter
        m_optimSettings.gd_settings.par_momentum = m_options.getReal("Momentum");

        // Ada parameters
        m_optimSettings.gd_settings.par_ada_norm_term = m_options.getReal("AdaNormTerm");
        m_optimSettings.gd_settings.par_ada_rho = m_options.getReal("AdaRho");
        m_optimSettings.gd_settings.ada_max = m_options.getSwitch("AdaMax");

        // Adam parameters
        m_optimSettings.gd_settings.par_adam_beta_1 = m_options.getReal("AdamBeta1");
        m_optimSettings.gd_settings.par_adam_beta_2 = m_options.getReal("AdamBeta2");

        // gradient clipping settings
        m_optimSettings.gd_settings.clip_grad = m_options.getSwitch("ClipGrad");
        m_optimSettings.gd_settings.clip_max_norm = m_options.getSwitch("ClipMaxNorm");
        m_optimSettings.gd_settings.clip_min_norm = m_options.getSwitch("ClipMinNorm");
        m_optimSettings.gd_settings.clip_norm_type = m_options.getInt("ClipNormType");
        m_optimSettings.gd_settings.clip_norm_bound = m_options.getReal("ClipNormBound");
    }

protected:

    using Base::m_options;
    using Base::m_optimSettings;
};

// NOTE: OMMITTED SINCE IT NEEDS HESSIAN!!
//
// template <typename T>
// class gsOptimNewton : public gsOptim<T>
// {
// public:
//     using Base = gsOptim<T>;

// public:

//     gsOptimNewton(gsOptProblem<T> * problem) : Base(problem)
// {
//     this->defaultOptions();
// }

//     bool callOptim( gsVector<T> & x,
//                     gsOptProblem<T> & op,
//                     optim::algo_settings_t & optimSettings)
//     { return optim::newton(x, gsOptimWrapper<T>(op), nullptr,optimSettings); }

// protected:

//     using Base::m_options;
//     using Base::m_optimSettings;
// };

template <typename T>
class gsOptimNM : public gsOptim<T>
{
public:
    using Base = gsOptim<T>;

public:

    gsOptimNM(gsOptProblem<T> * problem) : Base(problem)
    {
        this->defaultOptions();
    }

    bool callOptim( gsVector<T> & x,
                    gsOptProblem<T> & op,
                    optim::algo_settings_t & optimSettings)
    { return optim::nm(x, gsOptimWrapper<T>(op), nullptr,optimSettings); }

    void defaultOptions() override
    {
        Base::defaultOptions();
        m_options.addSwitch("AdaptivePars","scale the contraction, expansion, and shrinkage parameters using the dimension of the optimization problem.",true);
        m_options.addReal("ParAlpha","Reflection parameter.",1.0);
        m_options.addReal("ParBeta","Contraction parameter.",0.5);
        m_options.addReal("ParGamma","Expansion parameter.",2.0);
        m_options.addReal("ParDelta","Shrinkage parameter.",0.5);
        m_options.addSwitch("CustomInitialSimplex","whether to use user-defined values for the initial simplex matrix.",false);
    }

    /// Initial simplex points
    /// From manual:
    ///user-defined values for the initial simplex (optional). Dimensions: (n+1) x n
    void setSimplexPoints(const gsMatrix<T> & points)
    {
        m_optimSettings.nm_settings.initial_simplex_points = points;
    }

    void getOptions() override
    {
        Base::getOptions();
        m_optimSettings.nm_settings.adaptive_pars = m_options.getSwitch("AdaptivePars");

        m_optimSettings.nm_settings.par_alpha = m_options.getReal("ParAlpha");
        m_optimSettings.nm_settings.par_beta  = m_options.getReal("ParBeta");
        m_optimSettings.nm_settings.par_gamma = m_options.getReal("ParGamma");
        m_optimSettings.nm_settings.par_delta = m_options.getReal("ParDelta");

        m_optimSettings.nm_settings.custom_initial_simplex = m_options.getSwitch("CustomInitialSimplex");
    }

protected:

    using Base::m_options;
    using Base::m_optimSettings;
};

template <typename T>
class gsOptimDE : public gsOptim<T>
{
public:
    using Base = gsOptim<T>;

public:

    gsOptimDE(gsOptProblem<T> * problem) : Base(problem)
    {
        this->defaultOptions();
    }

    bool callOptim( gsVector<T> & x,
                    gsOptProblem<T> & op,
                    optim::algo_settings_t & optimSettings)
    { return optim::de(x, gsOptimWrapper<T>(op), nullptr,optimSettings); }

    void defaultOptions() override
    {
        Base::defaultOptions();
        m_options.addInt("nPop","size of population for each generation.",200);
        m_options.addInt("nPopBest","number of best individuals to select.",6);
        m_options.addInt("nGen","number of generations.",1000);
        m_options.addInt("MutationMethod","the mutation strategy, as described in step one of the algorithm description: 1) rand, 2) best.",1);
        m_options.addInt("CheckFreq","how many generations to skip when evaluating whether the best candidate value has improved between generations (i.e., to check for potential convergence).",std::numeric_limits<index_t>::max()-1);
        m_options.addReal("ParF","the mutation parameter.",0.8);
        m_options.addReal("ParCR","the crossover parameter.",0.9);
        m_options.addInt("PMax","Maximum number of perturbations.",4);
        m_options.addInt("MaxFnEval","Maximum number of function evaluations.",100000);
        m_options.addReal("ParFL","Lower bound for the scaling factor.",0.1);
        m_options.addReal("ParFU","Upper bound for the scaling factor.",1.0);
        m_options.addReal("ParTauF","Scaling factor for the scaling factor.",0.1);
        m_options.addReal("ParTauCR","Scaling factor for the crossover probability.",0.1);
        m_options.addReal("ParDEps","Small value for numerical stability.",0.5);
        m_options.addSwitch("ReturnPopulationMat","Whether to return the population matrix.",false);
    }

    /// Set the Upper and lower bounds of the uniform distributions used to generate the initial population
    void setBounds(const gsMatrix<T,Dynamic,2> & bounds)
    {
        m_optimSettings.de_settings.initial_lb = bounds.col(0);
        m_optimSettings.de_settings.initial_ub = bounds.col(1);
    }

    gsMatrix<T> getPopulationMatrix()
    {
        GISMO_ASSERT(m_optimSettings.de_settings.return_population_mat,"The option ReturnPopulationMat was not set to true.");
        return m_optimSettings.de_settings.population_mat;
    }

    void getOptions() override
    {
        Base::getOptions();
        m_optimSettings.de_settings.n_pop = m_options.getInt("nPop");
        m_optimSettings.de_settings.n_pop_best = m_options.getInt("nPopBest");
        m_optimSettings.de_settings.n_gen = m_options.getInt("nGen");

#   ifdef _OPENMP
        m_optimSettings.de_settings.omp_n_threads = omp_get_num_threads(); // numbers of threads to use
#   else
        m_optimSettings.de_settings.omp_n_threads = -1; // numbers of threads to use
#   endif

        m_optimSettings.de_settings.mutation_method = m_options.getInt("MutationMethod");

        m_optimSettings.de_settings.check_freq = m_options.getInt("CheckFreq");

        m_optimSettings.de_settings.par_F = m_options.getReal("ParF");
        m_optimSettings.de_settings.par_CR = m_options.getReal("ParCR");

        // DE-PRMM specific

        m_optimSettings.de_settings.pmax = m_options.getInt("PMax");
        m_optimSettings.de_settings.max_fn_eval = m_options.getInt("MaxFnEval");

        m_optimSettings.de_settings.par_F_l = m_options.getReal("ParFL");
        m_optimSettings.de_settings.par_F_u = m_options.getReal("ParFU");

        m_optimSettings.de_settings.par_tau_F  = m_options.getReal("ParTauF");
        m_optimSettings.de_settings.par_tau_CR = m_options.getReal("ParTauCR");

        m_optimSettings.de_settings.par_d_eps = m_options.getReal("ParDEps");
    }

protected:

    using Base::m_options;
    using Base::m_optimSettings;
};

template <typename T>
class gsOptimDEPRMM : public gsOptimDE<T>
{
    using Base = gsOptimDE<T>;
public:
    gsOptimDEPRMM(gsOptProblem<T> * problem) : Base(problem)
    {
        this->defaultOptions();
    }

    bool callOptim( gsVector<T> & x,
                    gsOptProblem<T> & op,
                    optim::algo_settings_t & optimSettings)
    { return optim::de_prmm(x, gsOptimWrapper<T>(op), nullptr,optimSettings); }
};

template <typename T>
class gsOptimPSO : public gsOptim<T>
{
public:
    using Base = gsOptim<T>;

public:

    gsOptimPSO(gsOptProblem<T> * problem) : Base(problem)
    {
        this->defaultOptions();
    }

    bool callOptim( gsVector<T> & x,
                    gsOptProblem<T> & op,
                    optim::algo_settings_t & optimSettings)
    { return optim::pso(x, gsOptimWrapper<T>(op), nullptr,optimSettings); }

    void defaultOptions() override
    {
        Base::defaultOptions();
        m_options.addSwitch("CenterParticle","whether to add a particle that averages across the population in each generation.",true);
        m_options.addInt("nPop","size of population for each generation.",200);
        m_options.addInt("nGen","number of generations.",1000);
        m_options.addInt("OMPThreads","numbers of threads to use.",-1);
        m_options.addInt("InertiaMethod","1 for linear decreasing between w_min and w_max; 2 for dampening.",1);
        m_options.addInt("CheckFreq","how many generations to skip when evaluating whether the best candidate value has improved between generations (i.e., to check for potential convergence).",std::numeric_limits<index_t>::max()-1);
        // m_options.addReal("ParInitialW","Initial value for the inertia weight.",1.0);
        m_options.addReal("ParWDamp","initial value of the weight parameter w.",0.99);
        m_options.addReal("ParWMin","lower bound on the weight parameter w.",0.10);
        m_options.addReal("ParWMax","upper bound on the weight parameter w.",0.99);
        m_options.addInt("VelocityMethod","set velocity method (1 for fixed values or 2 for linear change from initial to final values).",1);
        m_options.addReal("ParCCog","value for c_c.",2.0);
        m_options.addReal("ParCSoc","value for c_s.",2.0);
        // m_options.addReal("ParInitialCCog","initial value for c_c.",2.5);
        m_options.addReal("ParFinalCCog","final value for c_c.",0.5);
        // m_options.addReal("ParInitialCSoc","initial value for c_s.",0.5);
        m_options.addReal("ParFinalCSoc","final value c_s.",2.5);
        m_options.addSwitch("ReturnPopulationMat","Whether to return the population matrix.",false);
    }

    /// Set the Upper and lower bounds of the uniform distributions used to generate the initial population
    void setBounds(const gsMatrix<T,Dynamic,2> & bounds)
    {
        m_optimSettings.pso_settings.initial_lb = bounds.col(0);
        m_optimSettings.pso_settings.initial_ub = bounds.col(1);
    }

    gsMatrix<T> getPopulationMatrix()
    {
        GISMO_ASSERT(m_optimSettings.pso_settings.return_position_mat,"The option ReturnPopulationMat was not set to true.");
        return m_optimSettings.pso_settings.position_mat;
    }

    void getOptions() override
    {
        Base::getOptions();
        m_optimSettings.pso_settings.center_particle = m_options.getSwitch("CenterParticle");
        m_optimSettings.pso_settings.n_pop = m_options.getInt("nPop");
        m_optimSettings.pso_settings.n_gen = m_options.getInt("nGen");

#   ifdef _OPENMP
        m_optimSettings.pso_settings.omp_n_threads = omp_get_num_threads(); // numbers of threads to use
#   else
        m_optimSettings.pso_settings.omp_n_threads = -1; // numbers of threads to use
#   endif

        m_optimSettings.pso_settings.inertia_method = m_options.getInt("InertiaMethod");
        m_optimSettings.pso_settings.check_freq = m_options.getInt("CheckFreq");
        // m_optimSettings.pso_settings.par_initial_w = m_options.getReal("ParInitialW");
        m_optimSettings.pso_settings.par_w_damp = m_options.getReal("ParWDamp");
        m_optimSettings.pso_settings.par_w_min = m_options.getReal("ParWMin");
        m_optimSettings.pso_settings.par_w_max = m_options.getReal("ParWMax");
        m_optimSettings.pso_settings.velocity_method = m_options.getInt("VelocityMethod");
        m_optimSettings.pso_settings.par_c_cog = m_options.getReal("ParCCog");
        m_optimSettings.pso_settings.par_c_soc = m_options.getReal("ParCSoc");
        // m_optimSettings.pso_settings.par_initial_c_cog = m_options.getReal("ParInitialCCog");
        m_optimSettings.pso_settings.par_final_c_cog   = m_options.getReal("ParFinalCCog");
        // m_optimSettings.pso_settings.par_initial_c_soc = m_options.getReal("ParInitialCSoc");
        m_optimSettings.pso_settings.par_final_c_soc   = m_options.getReal("ParFinalCSoc");
        m_optimSettings.pso_settings.return_position_mat = m_options.getSwitch("ReturnPopulationMat");
    }

protected:

    using Base::m_options;
    using Base::m_optimSettings;
};

template <typename T>
class gsOptimPSODV : public gsOptimPSO<T>
{
    using Base = gsOptimPSO<T>;
public:

    gsOptimPSODV(gsOptProblem<T> * problem) : Base(problem)
    {
        this->defaultOptions();
    }

    bool callOptim( gsVector<T> & x,
                    gsOptProblem<T> & op,
                    optim::algo_settings_t & optimSettings)
    { return optim::pso_dv(x, gsOptimWrapper<T>(op), nullptr,optimSettings); }
};

template <typename T>
class gsOptimSUMT : public gsOptim<T>
{
public:
    using Base = gsOptim<T>;

public:

    gsOptimSUMT(gsOptProblem<T> * problem) : Base(problem)
    {
        this->defaultOptions();
    }

    bool callOptim( gsVector<T> & x,
                    gsOptProblem<T> & op,
                    optim::algo_settings_t & optimSettings)
    { return optim::sumt(x, gsOptimWrapper<T>(op),nullptr,gsOptimWrapperConstraint<T>(op),nullptr,optimSettings); }

    void solve(const gsMatrix<T> & initialGuess) override
    {
        // Get the bounds
        this->getOptions();
        this->setConstraints();
        GISMO_ASSERT(initialGuess.cols()==1,"The initial guess should have vector format");
        gsVector<T> x = initialGuess.col(0);
        m_success = callOptim(x, *m_op, m_optimSettings);
        m_curDesign = x;
        m_numIterations = m_optimSettings.opt_iter;
        m_finalObjective = m_optimSettings.opt_fn_value;
    }

    void defaultOptions() override
    {
        Base::defaultOptions();
        m_options.addReal("ParEta","eta parameter.",10.0);
    }

    void getOptions() override
    {
        Base::getOptions();
        m_optimSettings.sumt_settings.par_eta = m_options.getReal("ParEta");
    }

protected:

    using Base::m_options;
    using Base::m_optimSettings;

    using Base::m_op;
    using Base::m_success;
    using Base::m_curDesign;
    using Base::m_numIterations;
    using Base::m_finalObjective;
};

} // end namespace gismo

// // note: statically compiled in header-only mode
// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsOptim.hpp)
// #endif
