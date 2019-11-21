/** @file gsOptProblem.h

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

template <typename T> class gsIpOptTNLP;
class gsOptProblemPrivate;

/**
   \brief Class defining an optimization problem
*/

template <typename T>
class gsOptProblem
{

    friend class gsIpOptTNLP<T>;

public:

    /** default constructor */
    gsOptProblem();

    /** default destructor */
    virtual ~gsOptProblem();

public:

    /// \brief Returns the gradient value of the objective function at design
    /// value \a u
    virtual T evalObj( const gsAsConstVector<T> & u ) const = 0;

    /// \brief Returns the gradient of the objective function at design value
    /// \a u
    /// By default it uses finite differences, overriding it should provide exact gradient.
    virtual void gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const;

    /// \brief Returns values of the constraints at design value \a u
    virtual void evalCon_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const = 0;

    /// \brief Returns Jacobian of the constraints at design value \a u.
    /// Format of \a result is sparse, complying to \a m_conJacRows
    /// and \a m_conJacCols
    virtual void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const = 0;

    /// \brief Returns Hessian Lagrangian of the constraints at design value
    virtual void hessLagr_into( const gsAsConstVector<T> &, gsAsVector<T> &) const
    {GISMO_NO_IMPLEMENTATION }

    /// @brief Computes the sparsity pattern of the constraint Jacobian matrix.
    ///
    /// Number of constraints and design variables need to be set
    /// before calling this.  By default the constraint Jacobian is
    /// set to full. Override this function to provide sparsity.
    virtual void computeJacStructure()
    {
        // full Jacobian
        m_conJacRows.resize(m_numConstraints*m_numDesignVars);
        m_conJacCols.resize(m_numConstraints*m_numDesignVars);
        index_t c = 0;
        for ( index_t j = 0; j < m_numDesignVars; ++j)
            for ( index_t i = 0; i < m_numConstraints; ++i)
            {
                m_conJacRows[c  ] = i;
                m_conJacCols[c++] = j;
            }

        m_numConJacNonZero = m_conJacRows.size();
    }

    //void evalSingleCon_into(int k, gsMatrix<T> & result);
    //void hassObj_into( gsMatrix<T> & result);


    /// @brief Callback function is executed after every
    ///    iteration. Returning false causes premature termination of
    ///    the optimization
    virtual bool intermediateCallback() { return true;}

public:

    int numDesignVars () const { return m_curDesign.size(); }

    int numConstraints() const { return m_conLowerBounds.size(); }

    int numConDerivs  () const { return m_conJacRows.size(); }

    T lowerCon(int k) const { return m_conLowerBounds[k]; }

    T upperCon(int k) const { return m_conUpperBounds[k]; }

    T lowerDesignVar(int k) const { return m_desLowerBounds[k]; }

    T upperDesignVar(int k) const { return m_desUpperBounds[k]; }

    const std::vector<index_t> & conJacRows() const { return m_conJacRows; }

    const std::vector<index_t> & conJacCols() const { return m_conJacCols; }

    const gsMatrix<T> & currentDesign() const { return m_curDesign; }

    const gsMatrix<T> & lambda() const { return m_lambda; }

    T currentObjValue() const
    {
        gsAsConstVector<T> tmp(m_curDesign.data(), m_numDesignVars);
        return evalObj(tmp);
    }

    T objective()    const { return finalObjective; }

    int iterations() const { return numIterations; }

public:

    void solve ();

    std::ostream &print(std::ostream &os) const
    {
        os << "Design variables:" << m_numDesignVars
           << "\nNumber of constraints: " << m_numConstraints
            //<< "\ndesign lower:" << m_desLowerBounds.transpose()
            //<< "\ndesign upper:" << m_desUpperBounds.transpose()
            //<< "\nconstr lower:" << m_conLowerBounds.transpose()
            //<< "\nconstr upper:" << m_conUpperBounds.transpose()
           << "\nNumber of active sensitivities: " << m_numConJacNonZero
           << "\nSparsity of sensitivities: " <<std::fixed<<std::setprecision(2)
           << (100.0*m_numConJacNonZero)/(m_numDesignVars*m_numConstraints) << " %"
            //<< "\nm_conJacRows:" << m_conJacRows.transpose()
            //<< "\nm_conJacCols:" << m_conJacCols.transpose()
            //<< "\nm_curDesign:" << m_curDesign.transpose()
            ;
        return os;
    }


protected:

    /// Number of design variables
    int m_numDesignVars;

    /// Number of constraints
    int m_numConstraints;

    /// Number of nonzero entries in the Constraint Jacobian
    int m_numConJacNonZero;

    /// Lower bounds for the design variables
    gsVector<T> m_desLowerBounds;

    /// Upper bounds for the design variables
    gsVector<T> m_desUpperBounds;

    /// Lower bounds for the constraints
    gsVector<T> m_conLowerBounds;

    /// Upper bounds for the constraints
    gsVector<T> m_conUpperBounds;

    /// Constraint Jacobian non-zero entries rows
    std::vector<index_t> m_conJacRows;

    /// Constraint Jacobian non-zero entries columns
    std::vector<index_t> m_conJacCols;

    /// Current design variables (and starting point )
    gsMatrix<T> m_curDesign;

    /// Lagrange multipliers (set in the finalize_solution method)
    gsMatrix<T> m_lambda;

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
    gsOptProblem(const gsOptProblem & );
    gsOptProblem& operator=(const gsOptProblem & );
    //@}

    gsOptProblemPrivate * m_data;
};


} // end namespace gismo

// note: statically compiled in header-only mode
// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsOptProblem.hpp)
// #endif
