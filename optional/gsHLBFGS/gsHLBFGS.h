#include <functional>
#include <vector>
#include <cassert>

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsOptimizer/gsOptimizer.h>

#include "HLBFGS/HLBFGS.h"
/*
To do:
- Use Eigen
- change T, int
*/

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
        m_options.addReal("MinGradientLength","Minimal gradient length",1e-9);
        m_options.addReal("MinStepLength","Minimal step length",1e-9);
        m_options.addInt("LBFGSUpdates","Number of LBFGS updates (typically 3-20, put to 0 for gradient descent)",20);
    }

    void getOptions()
    {
        Base::getOptions();
        m_minGradientLength = m_options.getReal("MinGradientLength");
        m_minStepLength = m_options.getReal("MinStepLength");
        m_M = m_options.getInt("LBFGSUpdates");

        // m_hlbfgs_info[3]:The lbfgs strategy. 0: standard, 1: M1QN3 strategy (recommended)
        // Gilbert, J. C., & Lemaréchal, C. (1989). Some numerical experiments with variable-storage
        // quasi-Newton algorithms. Mathematical programming, 45(1), 407-435.
        m_hlbfgs_info[3] = 1;

        m_hlbfgs_info[4] = static_cast<index_t>(m_maxIterations);
        m_hlbfgs_info[5] = static_cast<index_t>(m_verbose);

        m_hlbfgs_pars[5] = m_minGradientLength;
        m_hlbfgs_pars[6] = m_minStepLength;
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
                m_M, // hardcoded??? -->>> change to an option of the class
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

// Options
protected:
    T m_minGradientLength;
    T m_minStepLength;

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

