/** @file gsLevenbergMarquardt.h
    @brief Provides declaration of the gradient descent method.
    This file is part of the G+Smo library.
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsIO/gsOptionList.h>
#include <gsOptimizer/gsOptimizer.h>
#include <gsOptimizer/gsOptProblem.h>
#define Eigen gsEigen
#include "lsqcpp.h"
#undef Eigen


namespace gismo
{

template<typename T>
struct gsLevenbergMarquardtObjective
{
    typedef gsEigen::Matrix<T, gsEigen::Dynamic, 1> Vector;
    typedef gsEigen::Matrix<T, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;

    gsLevenbergMarquardtObjective(gsOptProblem<T>* objective)
    :
    obj(objective)
    { }

    gsLevenbergMarquardtObjective()
    :
    obj(nullptr)
    { }

    T operator()(const Vector &xval, Vector &fval, Matrix &jacobian) const
    {
        gsAsConstVector<T> xvec(xval.data(),xval.size());
        const T val = obj->evalObj(xvec);
        fval.resize(1);
        fval(0,0) = val;

        jacobian.resize(1, obj->numDesignVars());
        gsAsVector<T> gvec(jacobian.data(), jacobian.size());
        obj->gradObj_into(xvec,gvec);
        return val;
    }

    gsOptProblem<T> * obj;
};

/**
 * @brief      This class describes the gradient descent method
 *
 * @tparam     T                  Real type
 * @tparam     StepSize           StepSize option, see external/gdcpp.h
 * @tparam     Callback           Callback option, see external/gdcpp.h
 * @tparam     FiniteDifferences  FiniteDifferences option, see external/gdcpp.h
 *
 * @ingroup    Optimizer
 */
template<typename T = real_t,
         typename StepSize=lsq::BarzilaiBorwein<T>,
         typename Callback=lsq::NoCallback<T>,
         typename FiniteDifferences=lsq::CentralDifferences<T> >
class gsLevenbergMarquardt : public gsOptimizer<T>
{
    using Base = gsOptimizer<T>;

    typedef typename lsq::LevenbergMarquardt<T, gsLevenbergMarquardtObjective<T>>::Result Result;

public:
    gsLevenbergMarquardt(gsOptProblem<T> * problem)
    :
    Base(problem),
    m_solver()
    {
        this->defaultOptions();
        gsLevenbergMarquardtObjective<T> obj(m_op);
        m_solver.setErrorFunction(obj);
    }


public:
    // const gsMatrix<T> & lambda() const { return m_lambda; }

    void minimize(const gsMatrix<T> &initialGuess)
    {
        m_result = m_solver.minimize(initialGuess);
        m_numIterations = m_result.iterations;
        m_finalObjective = m_result.fval.value(); //scalar assumed
        m_curDesign = m_result.xval;
    }

    Result result() { return m_result; };

protected:
    void defaultOptions()
    {
        Base::defaultOptions();
        m_options.addReal("MinGradLen","Minimal gradient length",1e-9);
        m_options.addReal("MinStepLen","Minimal step length",1e-9);
    }

    void getOptions()
    {
        Base::getOptions();
        m_solver.setMaxIterations(m_maxIterations);
        m_solver.setMinGradientLength( m_options.getReal("MinGradLen") ); // same notation as gsHLBFGS
        m_solver.setMinStepLength( m_options.getReal("MinStepLen") ); // same notation as gsHLBFGS
        m_solver.setVerbosity(m_verbose);

        // Set the minimum least squares error.
        m_solver.setMinError(0);

        // Set the parameters of the step method (Levenberg Marquardt).
        //m_solver.setMethodParameters({1.0, 2.0, 0.5, 100});
    }

public:

    void solve(const gsMatrix<T> &initialGuess)
    {
        this->getOptions();
        this->minimize(initialGuess);
    }

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

    Result m_result;

protected:
    T m_minGradientLength;

    lsq::LevenbergMarquardt<T, gsLevenbergMarquardtObjective<T> > m_solver;

};


// using gsLevenbergMarquardt = lsq::LevenbergMarquardt<T, Objective, StepSize, Callback, FiniteDifferences>;

} //namespace gismo
