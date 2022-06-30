/** @file gsOptimizer.h

    @brief Provides declaration of an optimization problem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsCore/gsLinearAlgebra.h>

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
    gsOptimizer();

    /** default destructor */
    virtual ~gsOptimizer();

public:

    /// @brief Callback function is executed after every
    ///    iteration. Returning false causes premature termination of
    ///    the optimization
    virtual bool intermediateCallback() { return true;}

public:

    const gsMatrix<T> & currentDesign() const { return m_curDesign; }

    T currentObjValue() const
    {
        gsAsConstVector<T> tmp(m_curDesign.data(), m_numDesignVars);
        return m_op->evalObj(tmp);
    }

    T objective()    const { return finalObjective; }

    int iterations() const { return numIterations; }

public:

    void solve () = 0;

    std::ostream &print(std::ostream &os) const
    {
        os << "gsOptimizer";
        return os;
    }


protected:

    gsOptProblem<T> * m_op;

    /// Current design variables (and starting point )
    gsMatrix<T> m_curDesign;

protected:

    // Statistics
    int numIterations;
    T   finalObjective;

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
#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsOptimizer.hpp)
#endif
