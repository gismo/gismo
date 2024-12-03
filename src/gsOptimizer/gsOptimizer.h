/** @file gsOptimizer.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsCore/gsLinearAlgebra.h>
#include <gsOptimizer/gsOptProblem.h>
#include <gsIO/gsOptionList.h>

#pragma once

namespace gismo
{

/**
   \brief Class defining an optimizer
*/
template <typename T>
class gsOptimizer
{
public:

    /** default constructor */
    gsOptimizer()
    :
    m_op(nullptr)
    {};

    gsOptimizer(gsOptProblem<T> * problem)
    :
    m_op(problem)
    {

    }

    /** default destructor */
    virtual ~gsOptimizer() { };

public:

    virtual void defaultOptions()
    {
        m_options.addInt("MaxIterations","Maximum iterations",100);
        m_options.addInt("Verbose","Verbosity level",0);
    }

    virtual void getOptions()
    {
        m_maxIterations = m_options.getInt("MaxIterations");
        m_verbose = m_options.getInt("Verbose");
    }

    /// @brief Callback function is executed after every
    ///    iteration. Returning false causes premature termination of
    ///    the optimization
    virtual bool intermediateCallback() { return true;}

public:

    const gsMatrix<T> & currentDesign() const { return m_curDesign; }
    gsMatrix<T> & currentDesign() { return m_curDesign; }

    T currentObjValue() const
    {
        gsAsConstVector<T> tmp(m_curDesign.data(), m_op->numDesignVars());
        return m_op->evalObj(tmp);
    }

    T objective()    const { return m_finalObjective; }

    int iterations() const { return m_numIterations; }

    gsOptionList & options() { return m_options; }

public:

    virtual void solve (const gsMatrix<T> & initialGuess) = 0;
    virtual void solve ()
    {
        m_curDesign.resize(m_op->numDesignVars(),1);
        m_curDesign.setZero();
        this->solve(m_curDesign);
    }

    std::ostream &print(std::ostream &os) const
    {
        os << "gsOptimizer";
        return os;
    }


protected:

    gsOptProblem<T> * m_op;

    /// Current design variables (and starting point )
    gsMatrix<T>     m_curDesign;

    /// Options
    gsOptionList m_options;
    index_t m_maxIterations, m_verbose;

protected:

    // Statistics
    index_t     m_numIterations;
    T           m_finalObjective;

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
    gsOptimizer(const gsOptimizer & );
    gsOptimizer& operator=(const gsOptimizer & );
    //@}
};

} // end namespace gismo

// note: must be statically compiled in header-only mode
//#ifndef GISMO_BUILD_LIB
//#include GISMO_HPP_HEADER(gsOptimizer.hpp)
//#endif
