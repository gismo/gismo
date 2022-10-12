/** @file gsIpOpt.h

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
#include <gsIpopt/gsIpOpt.hpp>


#pragma once

namespace gismo
{

template <typename T> class gsIpOptTNLP;
class gsIpOptPrivate;
/**
   \brief Class defining an optimization problem
*/

template <typename T>
class gsIpOpt : public gsOptimizer<T>
{
public:
    using Base = gsOptimizer<T>;

    friend class gsIpOptTNLP<T>;

public:

    /** default constructor */
    gsIpOpt(gsOptProblem<T> * problem);

    /** default destructor */
    virtual ~gsIpOpt();

    void defaultOptions()
    {
        // m_options.addInt("MaxIterations","Maximum iterations",1e3);
        // m_options.addReal("MinGradientLength","Minimal gradient length",1e-9);
        // m_options.addReal("MinStepLength","Minimal step length",1e-9);
        // m_options.addInt("Verbose","Verbosity level 0: no output, 1: summary, 2: full output", 0);
        // m_options.addInt("LBFGSUpdates","Number of LBFGS updates (typically 3-20, put to 0 for gradient descent)",20);
    }

    void getOptions()
    {
        // m_maxIterations = m_options.getInt("MaxIterations");
        // m_minGradientLength = m_options.getReal("MinGradientLength");
        // m_minStepLength = m_options.getReal("MinStepLength");
        // m_verbose = m_options.getInt("Verbose");
        // m_M = m_options.getInt("LBFGSUpdates");

        // // m_hlbfgs_info[3]:The lbfgs strategy. 0: standard, 1: M1QN3 strategy (recommended)
        // // Gilbert, J. C., & Lemaréchal, C. (1989). Some numerical experiments with variable-storage
        // // quasi-Newton algorithms. Mathematical programming, 45(1), 407-435.
        // m_hlbfgs_info[3] = 1;

        // m_hlbfgs_info[4] = static_cast<index_t>(m_maxIterations);
        // m_hlbfgs_info[5] = static_cast<index_t>((m_verbose==2));

        // m_hlbfgs_pars[5] = m_minGradientLength;
        // m_hlbfgs_pars[6] = m_minStepLength;
    }

public:

    /// @brief Callback function is executed after every
    ///    iteration. Returning false causes premature termination of
    ///    the optimization
    virtual bool intermediateCallback() { return true;}

public:

    const gsMatrix<T> & lambda() const { return m_lambda; }

public:

    void solve (const gsMatrix<T> & initialGuess);

    std::ostream &print(std::ostream &os) const
    {
        os << "Design variables:" << m_op->numDesignVars()
           << "\nNumber of constraints: " << m_op->numConstraints()
            //<< "\ndesign lower:" << m_op->desLowerBounds().transpose()
            //<< "\ndesign upper:" << m_op->desUpperBounds().transpose()
            //<< "\nconstr lower:" << m_op->conLowerBounds().transpose()
            //<< "\nconstr upper:" << m_op->conUpperBounds().transpose()
           << "\nNumber of active sensitivities: " << m_op->numConJacNonZero()
           << "\nSparsity of sensitivities: " <<std::fixed<<std::setprecision(2)
           << (100.0*m_op->numConJacNonZero())/(m_op->numDesignVars()*m_op->numConstraints()) << " %"
            //<< "\nm_op->conJacRows():" << m_op->conJacRows().transpose()
            //<< "\nm_op->conJacCols():" << m_op->conJacCols().transpose()
            //<< "\nm_op->curDesign():" << m_op->curDesign().transpose()
            ;
        return os;
    }


protected:

    /// Lagrange multipliers (set in the finalize_solution method)
    gsMatrix<T> m_lambda;

protected:

// Members taken from Base
protected:
    using Base::m_op;
    using Base::m_numIterations;
    using Base::m_finalObjective;
    using Base::m_curDesign;
    using Base::m_options;

    // Statistics
    int numIterations;
    T   finalObjective;

    gsIpOptPrivate * m_data;

private:

    /**@name Methods to block default compiler methods.
     * The compiler automatically generates the following three methods.
     *  Since the default compiler implementation is generally not what
     *  you want (for all but the most simple classes), we usually
     *  put the declarations of these methods in the private section
     *  and never implement them. This prevents the compiler from
     *  implementing an incorrect "default" behavior without us
     *  knowing. (See e.g. Scott Meyers book, "Effective C++")
     *
     */
    //@{
    gsIpOpt(const gsIpOpt & );
    gsIpOpt& operator=(const gsIpOpt & );
    //@}
};


} // end namespace gismo

// note: statically compiled in header-only mode
#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsIpOpt.hpp)
#endif
