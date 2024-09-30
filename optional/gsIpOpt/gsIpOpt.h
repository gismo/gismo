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

#pragma once

namespace gismo
{

template <typename T> class gsIpOptTNLP;
template <typename T> class gsIpOptPrivate;
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

    // // Currently not used
    // void defaultOptions()
    // {
    //     using Base::defaultOptions;
    // }

    void getOptions()
    {
        gsWarn<<"gsIpOpt.h has its options stored in filedata/options/ipopt.opt\n";
    }

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
    using Base::m_verbose;
    using Base::m_maxIterations;

    using Base::defaultOptions;
    using Base::getOptions;

    // Statistics

    gsIpOptPrivate<T> * m_data;

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
