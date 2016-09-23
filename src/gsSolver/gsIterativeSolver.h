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
#include <gsIO/gsOptionList.h>

namespace gismo
{
/// @brief Abstract class for iterative solvers.
///
/// \ingroup Solver
class gsIterativeSolver
{
public:
    typedef gsMatrix<real_t>    VectorType;

    typedef typename gsLinearOperator<>::Ptr LinOpPtr;
    
    gsIterativeSolver( const LinOpPtr& mat,
                       const LinOpPtr& precond)
    : m_mat(mat),
      m_precond(precond),
      m_max_iters(1000),
      m_tol(1e-10),
      m_num_iter(0),
      m_initial_error(0.),
      m_error(0.)
    {
        GISMO_ASSERT(m_mat->rows() == m_mat->cols(), "Matrix is not square.");
        if (!m_precond) m_precond = gsIdentityOp<>::make(m_mat->rows());
    }

    /// @brief  Contructor using any dense or sparse matrix
    ///
    /// @note: This does not copy the matrix. So, make sure that the
    /// matrix is not deleted before the solver.
    template<typename Derived>
    gsIterativeSolver( const Eigen::EigenBase<Derived> & mat,
                       const LinOpPtr& precond)
    : m_mat(makeMatrixOp(mat.derived())),
      m_precond(precond),
      m_max_iters(1000),
      m_tol(1e-10),
      m_num_iter(0),
      m_initial_error(0.),
      m_error(0.)
    {
        GISMO_ASSERT(m_mat->rows() == m_mat->cols(), "Matrix is not square.");
        if (NULL==m_precond.get()) m_precond = gsIdentityOp<>::make(m_mat->rows());
    }

    virtual ~gsIterativeSolver()    {}

    /// @brief Returns a list of default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addInt("MaxIterations", "Maximum number of iterations", 1000 );
        opt.addReal("Tolerance"   , "Numerical tolerance"         , 1e-10);
        return opt;
    }

    virtual void setOptions(const gsOptionList & opt)
    {
        m_max_iters = opt.askInt ("MaxIterations", m_max_iters );
        m_tol       = opt.askReal("Tolerance"    , m_tol       );
    }
    
    /// @brief Solves the linear system and stores the solution in \a x
    ///
    /// Solves the linear system of equations
    /// \param[in]     rhs      the right hand side of the linear system
    /// \param[in,out] x        starting value; the solution is stored in here
    ///
    /// \ingroup Solver
    void solve( const VectorType& rhs, VectorType& x )
    {
        GISMO_ASSERT( rhs.cols() == 1,
                      "Iterative solvers only work for single column right hand side." );
     
        GISMO_ASSERT( m_precond->rows() == m_mat->rows(),
                      "The preconditionner does not match the matrix." );

        GISMO_ASSERT( m_precond->cols() == m_mat->cols(),
                      "The preconditionner does not match the matrix." );
        
        m_num_iter = 0;
        
        /* // todo: check the following and uncomment
        m_rhsNorm = rhs.norm(); //m_initial_error
        if (0 == m_rhsNorm) // special case of zero rhs
        {
            x.setZero(rhs.rows()); // for sure zero is a solution
            m_error = 0.;
            return;
        }
        
        if ( 0 == x.size() ) // if no initial solution, start with zeros
            x.setZero(rhs.rows(), rhs.cols() );
        else
        {
           GISMO_ENSURE(m_mat->cols() == x.rows(), "Invalid initial solution");
           GISMO_ENSURE(rhs->cols() == x.cols()  , "Initial solution does not match right-hand side");
        }
        */

        if (initIteration(rhs, x))
        {
            m_error = 0.;
            return;
        }

        while (m_num_iter < m_max_iters)
        {
            m_num_iter++;
            if (step(x))
                break;
        }
        
        finalizeIteration( rhs, x );

    }

    virtual bool initIteration( const VectorType& rhs, VectorType& x ) = 0;
    virtual bool step( VectorType& x ) = 0;
    virtual void finalizeIteration( const VectorType& rhs, VectorType& x ) {}

    /// Returns the size of the linear system
    index_t size() const                                       { return m_mat->rows(); }

    /// Set the preconditionner. No copy is done, therefore \a precond
    /// must be a valid object while this iterative solver is used
    void setPreconditioner(const gsLinearOperator<> & precond) { setPreconditioner( memory::make_shared_not_owned( &precond ) ); }
    
    /// Set the preconditionner
    void setPreconditioner(const LinOpPtr & precond)           { m_precond = precond; }

    /// Set the maximum number of iterations (default: 1000)
    void setMaxIterations(index_t max_iters)                   { m_max_iters = max_iters; }

    /// Set the tolerance for the error criteria (default: 1e-10)
    void setTolerance(real_t tol)                              { m_tol = tol; }

    /// The number of iterations needed to reach the error criteria
    int iterations() const                                     { return m_num_iter; }

    /// The error of the iterative method
    real_t error() const                                       { return m_error; }

    /// The tolerance used in the iterative method
    real_t tolerance() const                                   { return m_tol; }


protected:
    const LinOpPtr     m_mat;
    /*const*/ LinOpPtr m_precond;
    index_t            m_max_iters;
    real_t             m_tol;
    index_t            m_num_iter;
    real_t             m_initial_error;
    real_t             m_error;

};

} // namespace gismo
