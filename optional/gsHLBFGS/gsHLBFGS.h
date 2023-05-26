#include <functional>
#include <vector>
#include <cassert>

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsIO/gsOptionList.h>
#include <gsOptimizer/gsOptimizer.h>
#include <gsOptimizer/gsOptProblem.h>
#include "HLBFGS/HLBFGS.h"
#include <gsUtils/gsStopwatch.h>
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


protected:
    void defaultOptions()
    {
        Base::defaultOptions();
        m_options.addReal("MinGradientLength","Minimal gradient length",1e-9);
        m_options.addReal("MinStepLength","Minimal step length",1e-10);
        m_options.addInt("LBFGSUpdates","Number of LBFGS updates (typically 3-20, put to 0 for gradient descent)",20);

        // see documentation in https://xueyuhanlang.github.io/software/HLBFGS/
        m_options.addReal("PARAMETER0","PARAMETERS0 upto PARAMETER6 are settable",1e-4);
        m_options.addInt("INFO0","INFO0 upto INFO12 are settable",20);
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

        for (int i = 0; i!=7; ++i)
            m_hlbfgs_pars[i] = m_options.askReal("PARAMETER"+util::to_string(i) , m_hlbfgs_pars[i]);

        for (int i = 0; i!=13; ++i)
            m_hlbfgs_info[i] = m_options.askInt("INFO"+util::to_string(i), m_hlbfgs_info[i]);
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

        const std::function<void(int N, double* x, double* prev_x, double* f, double* g)> wrapfunc =
            [&](int N, double* x, double*, double* f, double* g) {
            std::vector<double> array_x(N), array_g(N);

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

        // WHAT ABOUT CONSTRAINTS???? There are no.
        gsStopwatch time;
        time.restart();
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

        if (m_verbose==2)
            gsInfo<<"HLBFGS finished in "<<m_numIterations<<" iterations, with final objective "<<m_finalObjective<<"\n";

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
