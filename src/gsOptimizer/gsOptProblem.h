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
/**
   \brief Class defining an optimization problem
*/
template <typename T>
class gsOptProblem
{

public:

    // /** default constructor */
    // gsOptProblem();

    // /** default destructor */
    // virtual ~gsOptProblem();

public:

    /// \brief Returns the gradient value of the objective function at design
    /// value \a u
    virtual T evalObj( const gsAsConstVector<T> & u ) const
    {GISMO_NO_IMPLEMENTATION }

    /// \brief Returns the gradient of the objective function at design value
    /// \a u
    /// By default it uses finite differences, overriding it should provide exact gradient.
    virtual void gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
    {
        const index_t n = u.rows();
        //GISMO_ASSERT((index_t)m_numDesignVars == n*m, "Wrong design.");

        gsMatrix<T> uu = u;//copy
        gsAsVector<T> tmp(uu.data(), n);
        gsAsConstVector<T> ctmp(uu.data(), n);
        index_t c = 0;

        // for all partial derivatives (column-wise)
        for ( index_t i = 0; i!=n; i++ )
        {
            // to do: add m_desLowerBounds m_desUpperBounds check
            tmp[i]  += T(0.00001);
            const T e1 = this->evalObj(ctmp);
            tmp[i]   = u[i] + T(0.00002);
            const T e3 = this->evalObj(ctmp);
            tmp[i]   = u[i] - T(0.00001);
            const T e2 = this->evalObj(ctmp);
            tmp[i]   = u[i] - T(0.00002);
            const T e4 = this->evalObj(ctmp);
            tmp[i]   = u[i];
            result[c++]= ( 8 * (e1 - e2) + e4 - e3 ) / T(0.00012);
        }
    }


    /// \brief Returns values of the constraints at design value \a u
    virtual void evalCon_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
    {GISMO_NO_IMPLEMENTATION }

    /// \brief Returns Jacobian of the constraints at design value \a u.
    /// Format of \a result is sparse, complying to \a m_conJacRows
    /// and \a m_conJacCols
    virtual void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result ) const
    {GISMO_NO_IMPLEMENTATION }

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
    virtual void hessObj_into(const gsAsConstVector<T> &, gsAsMatrix<T> &) const
    {GISMO_NO_IMPLEMENTATION }

    /// @brief Callback function is executed after every
    ///    iteration. Returning false causes premature termination of
    ///    the optimization
    // virtual bool intermediateCallback() { return true;}

public:

    int numDesignVars () const { return m_curDesign.size(); }

    int numConstraints() const { return m_conLowerBounds.size(); }

    int numConJacNonZero() const { return m_numConJacNonZero; }

    int numConDerivs  () const { return m_conJacRows.size(); }

    T lowerCon(int k) const { return m_conLowerBounds[k]; }

    T upperCon(int k) const { return m_conUpperBounds[k]; }

    T lowerDesignVar(int k) const { return m_desLowerBounds[k]; }

    T upperDesignVar(int k) const { return m_desUpperBounds[k]; }

    const gsVector<T> & desLowerBounds() const {return m_desLowerBounds;}
    const gsVector<T> & desUpperBounds() const {return m_desUpperBounds;}

    const gsVector<T> & conLowerBounds() const {return m_conLowerBounds;}
    const gsVector<T> & conUpperBounds() const {return m_conUpperBounds;}

    const std::vector<index_t> & conJacRows() const { return m_conJacRows; }

    const std::vector<index_t> & conJacCols() const { return m_conJacCols; }

public:

    // void solve ();

    // std::ostream &print(std::ostream &os) const
    // {
    //     os << "Design variables: " << m_numDesignVars
    //        << "\nNumber of constraints: " << m_numConstraints
    //         //<< "\ndesign lower:" << m_desLowerBounds.transpose()
    //         //<< "\ndesign upper:" << m_desUpperBounds.transpose()
    //         //<< "\nconstr lower:" << m_conLowerBounds.transpose()
    //         //<< "\nconstr upper:" << m_conUpperBounds.transpose()
    //        << "\nNumber of active sensitivities: " << m_numConJacNonZero
    //        << "\nSparsity of sensitivities: " <<std::fixed<<std::setprecision(2)
    //        << (100.0*m_numConJacNonZero)/(m_numDesignVars*m_numConstraints) << " %"
    //         //<< "\nm_conJacRows:" << m_conJacRows.transpose()
    //         //<< "\nm_conJacCols:" << m_conJacCols.transpose()
    //         //<< "\nm_curDesign:" << m_curDesign.transpose()
    //         ;
    //     return os;
    // }


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
};


} // end namespace gismo

// note: must be statically compiled in header-only mode
#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsOptProblem.hpp)
#endif
