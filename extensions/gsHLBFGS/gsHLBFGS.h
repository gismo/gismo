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
- change real_t, int
- clean includes
- change .cpp to .hpp
*/

// https://xueyuhanlang.github.io/software/HLBFGS/

#include<iostream>
#include <iomanip>// Header file needed to use setw

namespace gismo
{

static void static_newiter_callback(int iter, int call_iter, real_t *x, real_t* f, real_t *g, real_t* gnorm)
{
    std::cout << "# iter "<< iter << ": #func eval. " << call_iter << ", f = " << *f <<
    ", ||g|| = " << *gnorm << std::endl;
    /*std::cout << "# iter "<< std::setw(4) << iter << ": #func eval. " << std::setw(4) << call_iter <<
        ", f = " << std::setiosflags(std::ios::fixed) << std::setiosflags(std::ios::right) << std::setprecision(6) << *f <<
        ", ||g|| = " << *gnorm << std::endl;*/
}

template<typename T>
class gsHLBFGS : public gsOptimizer<T>
{
public:
    using Base = gsOptimizer<T>;

    typedef std::function<void(const std::vector<real_t>& x, real_t &f, std::vector<real_t>& g)> func_grad_eval;

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

        func_grad =  [this](const std::vector<T>& X, T& F, std::vector<T>& G) {
            gsAsConstVector<real_t> u(X.data(),X.size());
            gsAsVector<real_t> Gvec(G.data(),G.size());
            F = m_op->evalObj(u);
            m_op->gradObj_into(u,Gvec);
            // objInterFreeFunc_into(X, F, G);
        };
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
        m_options.addInt("MaxIterations","Maximum iterations",0);
        m_options.addReal("MinGradientLength","Minimal gradient length",1e-9);
        m_options.addReal("MinStepLength","Minimal step length",1e-9);
        m_options.addInt("Verbose","Verbosity level",0);
    }

    void getOptions()
    {
        m_maxIterations = m_options.getInt("MaxIterations");
        m_minGradientLength = m_options.getReal("MinGradientLength");
        m_minStepLength = m_options.getReal("MinStepLength");
        m_verbose = m_options.getInt("Verbose");

        // m_hlbfgs_info[1] = static_cast<bool>(m_maxIterations);
        m_hlbfgs_info[4] = static_cast<bool>(m_maxIterations);
        m_hlbfgs_info[5] = static_cast<bool>(m_verbose);

        m_hlbfgs_pars[5] = gtol;
        m_hlbfgs_pars[6] = gtol;
    }

public:

protected:

    static void static_func_grad(int N, real_t* x, real_t* prev_x, real_t* f, real_t* g) {
        (*local_func_grad)(N, x, prev_x, f, g);
    }

public:
    // void run(const gsMatrix<T> & initialGuess)
    void run(std::vector<T> & sol)
    {
        // m_curDesign = initialGuess;

        INIT_HLBFGS(m_hlbfgs_pars, m_hlbfgs_info);

        // std::function<void(index_t N, real_t* x, real_t* prev_x, real_t* f, real_t* g)>

        std::function<void(int N, double* x, double* prev_x, double* f, double* g)> wrapfunc = [&](int N, double* x, double*, double* f, double* g) {
            std::vector<double> array_x(N), array_g(N);
            for (int i=0; i<N; i++)
                array_x[i] = x[i];

            func_grad(array_x, *f, array_g);

            for (int i=0; i<N; i++)
                g[i] = array_g[i];
        };

        local_func_grad = &wrapfunc;

        // WHAT ABOUT CONSTRAINTS????
        HLBFGS(
                sol.size(),
                5, // hardcoded??? -->>> change to an option of the class
                sol.data(),
                static_func_grad,
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
        this->getOptions();
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
    index_t m_verbose;

// HLBFGS options
protected:
    index_t m_hlbfgs_info[20] = {0};
    real_t  m_hlbfgs_pars[20] = {0};

    static const std::function<void(int N, real_t* x, real_t* prev_x, real_t* f, real_t* g)> * local_func_grad;

// to be removed
public:
    func_grad_eval func_grad;
    int maxiter = 10000; // Maximum number of quasi-Newton updates
    real_t gtol = 1e-10; // The iteration will stop when ||g||/max(1,||x||) <= gtol
    bool verbose = true;

};

template<typename T>
const std::function<void(int N, real_t* x, real_t* prev_x, real_t* f, real_t* g)> * gsHLBFGS<T>::local_func_grad = nullptr;


// using gsHLBFGS = gdc::GradientDescent<T, Objective, StepSize, Callback, FiniteDifferences>;

} //namespace gismo

