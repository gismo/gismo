#include <functional>
#include <vector>
#include <cassert>

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsIO/gsOptionList.h>
#include <gsOptimizer/gsOptimizer.h>
#include <gsOptimizer/gsOptProblem.h>
#include "HLBFGS/HLBFGS.h"
// #include "gsHLBFGS/HLBFGS_wrapper.cpp"
//#include "stlbfgs/stlbfgs.h"
/*
To do:
- Use Eigen
- change T, int
- clean includes
- change .cpp to .hpp
*/

// https://xueyuhanlang.github.io/software/HLBFGS/

#include<iostream>
#include <iomanip>// Header file needed to use setw

namespace gismo
{

template<typename T>
static void static_newiter_callback(int iter, int call_iter, T *x, T* f, T *g, T* gnorm)
{
    // TODO: give it a flag to decide if we need to print or not.
    std::cout << "# iter "<< iter << ": #func eval. " << call_iter << ", f = " << *f <<
    ", ||g|| = " << *gnorm << std::endl;
    /*std::cout << "# iter "<< std::setw(4) << iter << ": #func eval. " << std::setw(4) << call_iter <<
        ", f = " << std::setiosflags(std::ios::fixed) << std::setiosflags(std::ios::right) << std::setprecision(6) << *f <<
        ", ||g|| = " << *gnorm << std::endl;*/
}


struct gsHLBFGSObjective
{
    typedef double T;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vector;
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


public:
    // const gsMatrix<T> & lambda() const { return m_lambda; }

    // void minimize(const Vector &initialGuess)
    // {
    //     m_result = Optimizer::minimize(initialGuess);
    //     Base::m_numIterations = m_result.iterations;
    //     Base::m_finalObjective = m_result.fval;
    // }

    // Result result() { return m_result; };

protected:
    void defaultOptions()
    {
        m_options.addInt("MaxIterations","Maximum iterations",1e3);
        m_options.addReal("MinGradientLength","Minimal gradient length",1e-9);
        m_options.addReal("MinStepLength","Minimal step length",1e-9);
        m_options.addSwitch("Verbose","Verbosity level", true);
        m_options.addInt("LBFGSUpdates","Number of LBFGS updates (typically 3-20, put to 0 for gradient descent)",20);
    }

    void getOptions()
    {
        m_maxIterations = m_options.getInt("MaxIterations");
        m_minGradientLength = m_options.getReal("MinGradientLength");
        m_minStepLength = m_options.getReal("MinStepLength");
        m_verbose = m_options.getSwitch("Verbose");
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

    static void static_func_grad(int N, T* x, T* prev_x, T* f, T* g) {
        (*local_func_grad)(N, x, prev_x, f, g);
    }

public:
    // void run(const gsMatrix<T> & initialGuess)
    void run(std::vector<T> & sol)
    {
        // m_curDesign = initialGuess;

        INIT_HLBFGS(m_hlbfgs_pars, m_hlbfgs_info);
        this->getOptions();
        // std::function<void(index_t N, T* x, T* prev_x, T* f, T* g)>

        const std::function<void(int N, double* x, double* prev_x, double* f, double* g)> wrapfunc =
            [&](int N, double* x, double*, double* f, double* g) {
            std::vector<double> array_x(N), array_g(N);

            gsAsConstVector<real_t> u(x,N);
            *f = m_op->evalObj(u);
            
            gsAsVector<real_t> Gvec(g,N);
            m_op->gradObj_into(u,Gvec);
        };

        local_func_grad = &wrapfunc;

        //gsHLBFGSObjective obj(m_op);

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

    }

// Members taken from Base
protected:
    using Base::m_op;
    using Base::m_numIterations;
    using Base::m_finalObjective;
    using Base::m_curDesign;
    using Base::m_options;

// Options
protected:
    index_t m_maxIterations;
    T m_minGradientLength;
    T m_minStepLength;
    bool m_verbose;

// HLBFGS options
protected:
    index_t m_hlbfgs_info[20] = {0};
    T  m_hlbfgs_pars[20] = {0};
    index_t m_M;

    static const std::function<void(int N, T* x, T* prev_x, T* f, T* g)> * local_func_grad;

};

template<typename T>
const std::function<void(int N, T* x, T* prev_x, T* f, T* g)> * gsHLBFGS<T>::local_func_grad = nullptr;


// using gsHLBFGS = gdc::GradientDescent<T, Objective, StepSize, Callback, FiniteDifferences>;

} //namespace gismo

