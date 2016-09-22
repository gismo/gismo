/** @file gsIterativeSolver.h

    @brief Abstract class for iterative solvers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsCore/gsExport.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsSolver/gsMatrixOp.h>

namespace gismo
{
/// @brief Abstract class for iterative solvers.
///
/// \ingroup Solver
class gsIterativeSolver
{
public:
    typedef gsMatrix<real_t>    VectorType;

    /// Constructor for general linear operator
    gsIterativeSolver( const gsLinearOperator<>::Ptr& mat, index_t max_iters=1000, real_t tol=1e-10 )
        : m_mat(mat), m_max_iters(max_iters), m_tol(tol), m_num_iter(0), m_initial_error(0.), m_error(0.)
    {
        GISMO_ASSERT(m_mat->rows() == m_mat->cols(), "Matrix is not square, current implementation requires this!");
    }

    virtual ~gsIterativeSolver()    {}

    /// @brief Solves the linear system and stores the solution in \a x
    ///
    /// Solves the linear system of equations
    /// \param[in] rhs      the right hand side of the linear system
    /// \param[in,out] x    starting value; the solution is stored in here
    /// \param[in] precond  the preconditioner used (default: identity preconditioner)
    ///
    /// \ingroup Solver
    void solve( const VectorType& rhs, VectorType& x, const gsLinearOperator<> & precond )
    {
        m_num_iter = 0;
        
        if (initIteration(rhs, x, precond))
        {
            m_error = 0.;
            return;
        }

        while (m_num_iter < m_max_iters)
        {
            m_num_iter++;
            if (step(x, precond))
                break;
        }
        
        finalizeIteration( rhs, x );

    }
    
    /// Solve system without preconditioner
    void solve( const VectorType& rhs, VectorType& x )
    {
        gsIdentityOp<> preConId(m_mat->rows());
        solve(rhs, x, preConId);
    }
    
    virtual bool initIteration( const VectorType& rhs, VectorType& x, const gsLinearOperator<>& precond ) = 0;
    virtual bool step( VectorType& x, const gsLinearOperator<>& precond ) = 0;
    virtual void finalizeIteration( const VectorType& rhs, VectorType& x ) {}

    /// Returns the size of the linear system
    index_t size() const                     { return m_mat->rows(); }

    /// Set the maximum number of iterations (default: 1000)
    void setMaxIterations(index_t max_iters) { m_max_iters = max_iters; }

    /// Set the tolerance for the error criteria (default: 1e-10)
    void setTolerance(real_t tol)            { m_tol = tol; }

    /// The number of iterations needed to reach the error criteria
    int iterations() const                   { return m_num_iter; }

    /// The error of the iterative method
    real_t error() const                     { return m_error; }

    /// The tolerance used in the iterative method
    real_t tolerance() const                 { return m_tol; }


protected:
    const gsLinearOperator<>::Ptr m_mat;
    index_t                       m_max_iters;
    real_t                        m_tol;
    index_t                       m_num_iter;
    real_t                        m_initial_error;
    real_t                        m_error;

};

} // namespace gismo
