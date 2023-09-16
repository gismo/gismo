#include <functional>
#include <vector>
#include <cassert>

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsIO/gsOptionList.h>
#include <gsOptimizer/gsOptimizer.h>
#include <gsOptimizer/gsOptProblem.h>
#include "HLBFGS/HLBFGS.h"
#include "HLBFGS/Lite_Sparse_Matrix.h"
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
    Base(problem), m_hm(0, 0, SYM_LOWER, CCS, FORTRAN_TYPE, true)
    {
        this->defaultOptions();
    }


protected:
    void defaultOptions()
    {
        Base::defaultOptions();

        // see documentation in https://xueyuhanlang.github.io/software/HLBFGS/
        // parameters 0 ... 6:
        m_options.addReal("FuncTol",    "function tolerance used in line-search", 1e-4);
        m_options.addReal("VarTol",     "variable tolerance used in line-search", 1e-16);
        m_options.addReal("GradTol",    "gradient tolerance used in line-search", 0.9);
        m_options.addReal("StpMin",     "stpmin used in line-search",             1e-20);
        m_options.addReal("StpMax",     "stpmax used in line-search",             1e+20);
        m_options.addReal("MinGradLen", "Minimal gradient length",                1e-9);
        m_options.addReal("MinStepLen", "Minimal step length",                    1e-10);

        // infos 0 ... 2, 6 ... 12:
        m_options.addInt("MaxEval",    "the max number of evaluation in line-search", 20);
        m_options.addInt("TotEval",    "the total number of evalfunc calls", 0);
        m_options.addInt("CurrInt",    "the current number of iterations", 0);
        //m_options.addInt("Strategy", "The lbfgs strategy. 0: standard, 1: M1QN3 strategy[8](recommended).", 0);
        // maxIterations: set by the parent class
        // verbose:       set by the parent class
        m_options.addInt("UpdateHess", "T: the update interval of Hessian. (typical choices: 0-200)", 10);
        m_options.addInt("SwitchHess", "0: without hessian, 1: with accurate hessian", 0);
        m_options.addInt("Icfs",       "icfs parameter", 15);
        m_options.addInt("LineSearch", "0: classical line-search; 1: modified line-search (it is not useful in practice)", 0);
        m_options.addInt("DissPre",    "0: Disable preconditioned CG; 1: Enable preconditioned CG", 0);
        m_options.addInt("BetaCG",     "0 or 1 defines different methods for choosing beta in CG.", 1);
        m_options.addInt("Diag",       "internal usage. 0: only update the diag in USER_DEFINED_HLBFGS_UPDATE_H; 1: default.", 1);

        m_options.addInt("LBFGSUpdates","Number of LBFGS updates (typically 3-20, put to 0 for gradient descent)",20);
    }

    void getOptions()
    {
        Base::getOptions();

        // parameters
        m_funcTol =    m_options.getReal("FuncTol");
        m_varTol =     m_options.getReal("VarTol");
        m_gradTol =    m_options.getReal("GradTol");
        m_stpMin =     m_options.getReal("StpMin");
        m_stpMax =     m_options.getReal("StpMax");
        m_minGradLen = m_options.getReal("MinGradLen");
        m_minStepLen = m_options.getReal("MinStepLen");

        // infos
        m_maxEval =    m_options.getInt("MaxEval");
        m_totEval =    m_options.getInt("TotEval");
        m_currInt =    m_options.getInt("CurrInt");
        m_updateHess = m_options.getInt("UpdateHess");
        m_switchHess = m_options.getInt("SwitchHess");
        m_icfs =       m_options.getInt("Icfs");
        m_lineSearch = m_options.getInt("LineSearch");
        m_dissPre =    m_options.getInt("DissPre");
        m_betaCG =     m_options.getInt("BetaCG");
        m_diag =       m_options.getInt("Diag");

        m_M = m_options.getInt("LBFGSUpdates");

        m_hlbfgs_info[0] = static_cast<index_t>(m_maxEval);
        m_hlbfgs_info[1] = static_cast<index_t>(m_totEval);
        m_hlbfgs_info[2] = static_cast<index_t>(m_currInt);
        // m_hlbfgs_info[3]:The lbfgs strategy. 0: standard, 1: M1QN3 strategy (recommended)
        // Gilbert, J. C., & Lemaréchal, C. (1989). Some numerical experiments with variable-storage
        // quasi-Newton algorithms. Mathematical programming, 45(1), 407-435.
        m_hlbfgs_info[3] = 1;
        m_hlbfgs_info[4] = static_cast<index_t>(m_maxIterations);
        m_hlbfgs_info[5] = static_cast<index_t>(m_verbose);
        m_hlbfgs_info[6] = static_cast<index_t>(m_updateHess);
        m_hlbfgs_info[7] = static_cast<index_t>(m_switchHess);
        m_hlbfgs_info[8] = static_cast<index_t>(m_icfs);
        m_hlbfgs_info[9] = static_cast<index_t>(m_lineSearch);
        m_hlbfgs_info[10] = static_cast<index_t>(m_dissPre);
        m_hlbfgs_info[11] = static_cast<index_t>(m_betaCG);
        m_hlbfgs_info[12] = static_cast<index_t>(m_diag);

        m_hlbfgs_pars[0] = m_funcTol;
        m_hlbfgs_pars[1] = m_varTol;
        m_hlbfgs_pars[2] = m_gradTol;
        m_hlbfgs_pars[3] = m_stpMin;
        m_hlbfgs_pars[4] = m_stpMax;
        m_hlbfgs_pars[5] = m_minGradLen;
        m_hlbfgs_pars[6] = m_minStepLen;

        for (int i = 0; i!=7; ++i)
        {
            m_hlbfgs_pars[i] = m_options.askReal("PARAMETER"+util::to_string(i) , m_hlbfgs_pars[i]);
            //gsInfo << "m_hlbfgs_pars[" << i << "]: " << m_hlbfgs_pars[i] << std::endl;
        }

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

    static void static_hessian(int N, double *x, double *prev_x, double *f,
                               double *g, HESSIAN_MATRIX& hessian)
    {
        (*local_hessian)(N, x, prev_x, f, g, hessian);
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
            [&](int N, double* x, double*, double* f, double* g)
            {
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

        const std::function<void(int N, double *x, double *prev_x, double *f, double *g, HESSIAN_MATRIX& hessian)> wraphess =
            [&](int N, double *x, double *prev_x, double *f, double *g, HESSIAN_MATRIX& hessian)
            {
                gsAsConstVector<real_t> u(x,N);

                static bool first = true;
                //if (first)
                {
                    *f = m_op->evalObj(u);                
                    gsAsVector<real_t> Gvec(g,N);
                    m_op->gradObj_into(u,Gvec);
                }
                first = false;

                std::vector<double> tmp(N*N);
                gsAsMatrix<real_t> Hm(tmp,N,N);
                m_op->hessObj_into(u,Hm);
                
                m_hm = Lite_Sparse_Matrix(N, N, SYM_LOWER, CCS, FORTRAN_TYPE, true);
                m_hm.begin_fill_entry();
                double * diag = m_hm.get_diag();
                gsAsVector<real_t> Hdiag(diag, N);
                Hdiag = Hm.diagonal();
                for (int i = 0; i < N; ++i)
                    for (int j = i+1; j < N; ++j)
                        m_hm.fill_entry(i, j, Hm(i,j));

                m_hm.end_fill_entry(); 
                hessian.set_diag  (m_hm.get_diag());
                hessian.set_values(m_hm.get_values());
                hessian.set_rowind(m_hm.get_rowind());
                hessian.set_colptr(m_hm.get_colptr());
                hessian.set_nonzeros(m_hm.get_nonzero());
            };

        local_hessian = &wraphess;

        // WHAT ABOUT CONSTRAINTS???? There are no.
        gsStopwatch time;
        time.restart();
        HLBFGS(
                sol.size(),
                m_M,
                sol.data(),
                static_func_grad,
//                obj,
                static_hessian,
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
    // parameters
    T m_funcTol;
    T m_varTol;
    T m_gradTol;
    T m_stpMin;
    T m_stpMax;
    T m_minGradLen;
    T m_minStepLen;

    // infos
    index_t m_maxEval;
    index_t m_totEval;
    index_t m_currInt;
    index_t m_updateHess;
    index_t m_switchHess;
    index_t m_icfs;
    index_t m_lineSearch;
    index_t m_dissPre;
    index_t m_betaCG;
    index_t m_diag;

// HLBFGS options
protected:
    index_t m_hlbfgs_info[20] = {0};
    T  m_hlbfgs_pars[20] = {0};
    index_t m_M;

    Lite_Sparse_Matrix m_hm;
    
    static const std::function<void(int N, T* x, T* prev_x, T* f, T* g)> * local_func_grad;
    static const std::function<void(int iter, int call_iter, T *x, T* f, T *g, T* gnorm)> * local_newiter_callback;
    static const std::function<void(int N, double *x, double *prev_x, double *f, double *g, HESSIAN_MATRIX& hessian)> * local_hessian;
};


template<typename T>
const std::function<void(int N, T* x, T* prev_x, T* f, T* g)> * gsHLBFGS<T>::local_func_grad = nullptr;
template<typename T>
const std::function<void(int iter, int call_iter, T *x, T* f, T *g, T* gnorm)> * gsHLBFGS<T>::local_newiter_callback = nullptr;
template<typename T>
const std::function<void(int N, double *x, double *prev_x, double *f, double *g, HESSIAN_MATRIX& hessian)> * gsHLBFGS<T>::local_hessian = nullptr;


} //namespace gismo
