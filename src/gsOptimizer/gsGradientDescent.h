/** @file gsGradientDescent.h

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
#include "gdcpp.h"
//#include "lsqcpp.h"

namespace gismo
{

template<typename T>
struct gsGradientDescentObjective
{
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> Vector;
    // typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Matrix;

    gsGradientDescentObjective(gsOptProblem<T>* objective)
    :
    obj(objective)
    { }

    gsGradientDescentObjective()
    :
    obj(nullptr)
    { }

    T operator()(const Vector & vx, Vector & vfgrad) const
    {
        vfgrad.resize(obj->numDesignVars());
        gsAsConstVector<T> xvec(vx.data(),vx.size());
        gsAsVector<T> gvec(vfgrad.data(),vfgrad.size());
        const T val = obj->evalObj(xvec);
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
         typename StepSize=gdc::BarzilaiBorwein<T>,
         typename Callback=gdc::NoCallback<T>,
         typename FiniteDifferences=gdc::CentralDifferences<T> >
class gsGradientDescent : public gsOptimizer<T>
{
    using Base = gsOptimizer<T>;

    typedef typename gdc::GradientDescent<T, gsGradientDescentObjective<T>, StepSize, Callback, FiniteDifferences>::Result Result;

public:
    gsGradientDescent(gsOptProblem<T> * problem)
    :
    Base(problem),
    m_solver()
    {
        this->defaultOptions();
        gsGradientDescentObjective<T> obj(m_op);
        m_solver.setObjective(obj);
    }


public:
    // const gsMatrix<T> & lambda() const { return m_lambda; }

    void minimize(const gsMatrix<T> &initialGuess)
    {
        m_result = m_solver.minimize(initialGuess);
        m_numIterations = m_result.iterations;
        m_finalObjective = m_result.fval;
        m_curDesign = m_result.xval;
    }

    Result result() { return m_result; };

protected:
    void defaultOptions()
    {
        Base::defaultOptions();
        m_options.addReal("MinGradientLength","Minimal gradient length",1e-9);
        m_options.addReal("MinStepLength","Minimal step length",1e-9);
    }

    void getOptions()
    {
        Base::getOptions();
        m_minGradientLength = m_options.getReal("MinGradientLength");
        m_minStepLength = m_options.getReal("MinStepLength");

        m_solver.setMaxIterations(m_maxIterations);
        m_solver.setMinGradientLength(m_minGradientLength);
        m_solver.setMinStepLength(m_minStepLength);
        m_solver.setVerbosity(m_verbose);
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
    T m_minStepLength;

    gdc::GradientDescent<T, gsGradientDescentObjective<T>, StepSize, Callback, FiniteDifferences> m_solver;

};


// using gsGradientDescent = gdc::GradientDescent<T, Objective, StepSize, Callback, FiniteDifferences>;

} //namespace gismo
