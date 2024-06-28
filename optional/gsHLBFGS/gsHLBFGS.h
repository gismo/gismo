#include <functional>
#include <vector>
#include <cassert>

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsOptimizer/gsOptimizer.h>

#include "HLBFGS/HLBFGS.h"

// https://xueyuhanlang.github.io/software/HLBFGS/

#include<iostream>
#include <iomanip>// Header file needed to use setw

namespace gismo
{

struct gsHLBFGSObjective
{
    typedef real_t T;
    typedef gsEigen::Matrix<T, gsEigen::Dynamic, 1> Vector;
    // typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Matrix;

    gsHLBFGSObjective(gsOptProblem<T>* objective)
    :
    obj(objective)
    { }

    gsHLBFGSObjective()
    :
    obj(nullptr)
    { }

    void operator() (int N, real_t* x, real_t* prev_x, real_t* f, real_t* g) const
    {
        gsAsConstVector<real_t> u(x,N);
        *f = obj->evalObj(u);
        
        gsAsVector<real_t> Gvec(g,N);
        obj->gradObj_into(u,Gvec);
    }

    gsOptProblem<T> * obj;
};

template<typename T>
class gsHLBFGS : public gsOptimizer<T>
{
public:
    using Base = gsOptimizer<T>;

public:
    // gsHLBFGS(gsOptProblem<T> * problem)
    // :
    // Base(problem)
    // {
    //     this->defaultOptions();
    // }

    gsHLBFGS(gsOptProblem<T> * problem)
    :
    Base(problem)
    {
        this->defaultOptions();
    }


protected:
    void defaultOptions()
    {
        Base::defaultOptions();
        // See https://xueyuhanlang.github.io/software/HLBFGS/

        // DEPRECATED  
        m_options.addReal("MinGradientLength","Minimal gradient length",1e-9);
        m_options.addReal("MinStepLength","Minimal step length",1e-9);

        // Options from PARAMETERS
        m_options.addReal("tolF","Function tolerance (default 1e-4)",1e-4);
        m_options.addReal("tolL","Line search tolerance (default 1e-16)",1e-16);
        m_options.addReal("tolLG","Line search gradient tolerance (default 0.9)",0.9);
        m_options.addReal("minStepL","Minimal step size used in line-search (default 1e-20)",1e-20);
        m_options.addReal("maxStepL","Maximum step size used in line-search (default 1e+20)",1e20);
        m_options.addReal("tolRelG","Relative tolerance for the gradient (default 1e-9)",1e-9);
        m_options.addReal("tolG","Tolerance for the gradient (default 1e-10)",1e-10);

        // Options from INFO
        m_options.addInt("maxItL","the max number of evaluation in line-search (default 20)",20);
        m_options.addInt("strategy","The lbfgs strategy. 0: standard (default), 1: M1QN3 strategy [recommended].",0);
        m_options.addInt("hessUpdate","T: the update interval of Hessian [typical choices: 0-200] (default 10)",10);
        m_options.addInt("hessian","0: with hessian (default), 1: without hessian",0);
        m_options.addInt("ICFS","ICFS parameter (default 15)",15);
        m_options.addInt("linesearch","0: classical line-search (default), 1: modified line-search [not useful in practice]",0);
        m_options.addInt("preconditionedCG","0: disable preconditioned CG (default), 1: enable preconditioned CG",0);
        m_options.addInt("preconditionedCGpar","0 or 1 defines different methods for choosing beta in CG. (default 0)",1);
        // m_options.addInt("diagonalUpdate","internal usage. 0: only update the diag in USER_DEFINED_HLBFGS_UPDATE_H; 1: default. ",0);

        // Option M
        m_options.addInt("LBFGSUpdates","Number of LBFGS updates (typically 3-20, put to 0 for gradient descent)",20);
    }

    void getOptions()
    {
        Base::getOptions();

        T minGradientLength = m_options.askReal("MinGradientLength",-1);
        if (minGradientLength==-1)
            gsWarn<<"The option 'MinGradientLength' is deprecated, pease use 'tolRelG' instead\n";

        T minStepLength = m_options.askReal("MinStepLength",-1);
        if (minStepLength==-1)
            gsWarn<<"The option 'MinStepLength' is deprecated, pease use 'tolG' instead\n";

        // PARAMETERS
        m_hlbfgs_pars[0] = m_options.getReal("tolF");
        m_hlbfgs_pars[1] = m_options.getReal("tolL");
        m_hlbfgs_pars[2] = m_options.getReal("tolLG");
        m_hlbfgs_pars[3] = m_options.getReal("minStepL");
        m_hlbfgs_pars[4] = m_options.getReal("maxStepL");
        m_hlbfgs_pars[5] = (minGradientLength==-1) ? m_options.getReal("tolRelG") : minGradientLength;
        m_hlbfgs_pars[6] = (minStepLength==-1) ? m_options.getReal("tolG") : minStepLength;

        // INFO
        //      INFO[1,2] are output, INFO[4,5] come from Base
        m_hlbfgs_info[0]  = m_options.getInt("maxItL");
        m_hlbfgs_info[3]  = m_options.getInt("strategy");
        m_hlbfgs_info[4]  = static_cast<index_t>(m_maxIterations);
        m_hlbfgs_info[5]  = static_cast<index_t>(m_verbose);
        m_hlbfgs_info[6]  = m_options.getInt("hessUpdate");
        m_hlbfgs_info[7]  = m_options.getInt("hessian");
        m_hlbfgs_info[8]  = m_options.getInt("ICFS");
        m_hlbfgs_info[9]  = m_options.getInt("linesearch");
        m_hlbfgs_info[10] = m_options.getInt("preconditionedCG");
        m_hlbfgs_info[11] = m_options.getInt("preconditionedCGpar");
        // m_hlbfgs_info[12] = m_options.getInt("");

        // M
        m_M = m_options.getInt("LBFGSUpdates");

        // m_hlbfgs_info[3]:The lbfgs strategy. 0: standard, 1: M1QN3 strategy (recommended)
        // Gilbert, J. C., & Lemaréchal, C. (1989). Some numerical experiments with variable-storage
        // quasi-Newton algorithms. Mathematical programming, 45(1), 407-435.
    }

protected:

    static void static_func_grad(int N, T* x, T* prev_x, T* f, T* g)
    {
        (*local_func_grad)(N, x, prev_x, f, g);
    }

    static void static_newiter_callback(int iter, int call_iter, T *x, T* f, T *g, T* gnorm)
    {
        (*local_newiter_callback)(iter, call_iter, x, f, g, gnorm);
    }

public:
    // void run(const gsMatrix<T> & initialGuess)
    void run(std::vector<T> & sol)
    {
        // m_curDesign = initialGuess;

        INIT_HLBFGS(m_hlbfgs_pars, m_hlbfgs_info);
        this->getOptions();
        // std::function<void(index_t N, T* x, T* prev_x, T* f, T* g)>

        const std::function<void(int N, real_t* x, real_t* prev_x, real_t* f, real_t* g)> wrapfunc =
            [&](int N, real_t* x, real_t*, real_t* f, real_t* g) {
            std::vector<real_t> array_x(N), array_g(N);

            gsAsConstVector<real_t> u(x,N);
            *f = m_op->evalObj(u);
            
            gsAsVector<real_t> Gvec(g,N);
            m_op->gradObj_into(u,Gvec);
        };

        local_func_grad = &wrapfunc;

        const std::function<void(int iter, int call_iter, T *x, T* f, T *g, T* gnorm)> wrapcallback =
            [this](int iter, int call_iter, T *x, T* f, T *g, T* gnorm)
        {
            // TODO: give it a flag to decide if we need to print or not.
            if (m_verbose==2)
                gsInfo<<"# iter "<< iter << ": #func eval. " << call_iter << ", f = " << *f <<", ||g|| = " << *gnorm << std::endl;
        };

        local_newiter_callback = &wrapcallback;

        // WHAT ABOUT CONSTRAINTS????
        HLBFGS(
                sol.size(),
                m_M,
                sol.data(),
                static_func_grad,
//                obj,
                nullptr,
                HLBFGS_UPDATE_Hessian,
                static_newiter_callback,
                m_hlbfgs_pars,
                m_hlbfgs_info
              );
    }

    void solve(const gsMatrix<T> & initialGuess)
    {
        GISMO_ASSERT(initialGuess.cols()==1,"The initial guess should have vector format");
        std::vector<T> sol(initialGuess.rows());
        gsMatrix<T>::Map(sol.data(), initialGuess.rows(),1) = initialGuess;
        this->run(sol);

        m_curDesign = gsMatrix<T>::Map(sol.data(), sol.size(),1);
        m_numIterations = m_hlbfgs_info[2];
        m_finalObjective = m_op->evalObj(gsAsConstVector<T>(m_curDesign.data(),m_curDesign.rows()));

        if (m_verbose==1)
            gsInfo<<"HLBFGS finished in "<<m_numIterations<<" iterations, with final objective "<<m_finalObjective<<"\n";

    }

  /// gradientCheck subroutine, very useful to check the correctness of your
  // own analytic gradient. Author: Ye Ji (jiyess@outlook.com)
  void gradientCheck(const gsVector<T> &u) {
    // Get the analytic gradient
    std::vector<T> sol(u.size());
    gsAsVector<T> analyticGrad(sol);
    std::copy(u.begin(), u.end(), sol.begin());
    m_op->gradObj_into(sol, analyticGrad);

    // Finite difference calculation of gradient using central differences
    gsVector<T> numericalGrad(u.size());
    T h = sqrt(std::numeric_limits<T>::epsilon()) * u.cwiseAbs().maxCoeff();
    T forwardValue, backwardValue;

    std::vector<T> solForNumericalGrad(u.size());
    std::copy(u.begin(), u.end(), solForNumericalGrad.begin());

    // Iterate through each dimension
    for (int k = 0; k < u.size(); ++k) {
      // Compute function value at forward step
      solForNumericalGrad[k] += h;
      forwardValue = m_op->evalObj(solForNumericalGrad);

      // Compute function value at backward step
      solForNumericalGrad[k] -= 2.0 * h;
      backwardValue = m_op->evalObj(solForNumericalGrad);

      // Compute the numerical gradient using central difference formula
      numericalGrad(k) = (forwardValue - backwardValue) / (2.0 * h);

      // Reset the k-th component to its original value
      solForNumericalGrad[k] += h;
    }

    // Compare the analytic gradient and the numerical gradient
    gsInfo << "Analytical gradient:  Finite difference gradient: \n";

    int numElementsToPrint = std::min(analyticGrad.size(), 30);
    for (int i = 0; i < numElementsToPrint; ++i) {
      gsInfo << std::setw(5) << i << std::setw(20) << analyticGrad(i)
             << std::setw(20) << numericalGrad(i) << "\n";
    }

    if (u.size() > 30) {
      gsInfo << "(Displaying the first 30 components only)\n";
    }

    T relativeError =
        (analyticGrad - numericalGrad).norm() / analyticGrad.norm();
    gsInfo << "The relative error between the analytic gradient and the "
              "numerical gradient is: " << relativeError << "\n\n";
  }

// Members taken from Base
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

// HLBFGS options
protected:
    index_t m_hlbfgs_info[20] = {0};
    T  m_hlbfgs_pars[20] = {0};
    index_t m_M;

    static const std::function<void(int N, T* x, T* prev_x, T* f, T* g)> * local_func_grad;
    static const std::function<void(int iter, int call_iter, T *x, T* f, T *g, T* gnorm)> * local_newiter_callback;
};

template<typename T>
const std::function<void(int N, T* x, T* prev_x, T* f, T* g)> * gsHLBFGS<T>::local_func_grad = nullptr;
template<typename T>
const std::function<void(int iter, int call_iter, T *x, T* f, T *g, T* gnorm)> * gsHLBFGS<T>::local_newiter_callback = nullptr;


// using gsHLBFGS = gdc::GradientDescent<T, Objective, StepSize, Callback, FiniteDifferences>;

} //namespace gismo

