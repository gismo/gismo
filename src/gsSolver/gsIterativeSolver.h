/** @file gsIterativeSolver.h

    @brief Abstract class for iterative solvers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer, C. Hofreither, A. Manzaflaris, J. Sogn, S. Takacs
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
template<class T>
class gsIterativeSolver
{
public:
    typedef gsMatrix<T>    VectorType;

    typedef typename gsLinearOperator<T>::Ptr LinOpPtr;
    
    /// @brief Contructor using a linear operator to be solved for and
    ///  a preconditioner
    gsIterativeSolver( const LinOpPtr& mat,
                       const LinOpPtr& precond)
    : m_mat(mat),
      m_precond(precond),
      m_max_iters(1000),
      m_tol(1e-10),
      m_num_iter(0),
      m_rhs_norm(0.),
      m_error(0.)
    {
        GISMO_ASSERT(m_mat->rows() == m_mat->cols(), "Matrix is not square.");
        if (!m_precond) m_precond = gsIdentityOp<T>::make(m_mat->rows());
    }

    /// @brief Contructor using any dense or sparse matrix and a
    ///  preconditioner
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
      m_rhs_norm(0.),
      m_error(0.)
    {
        GISMO_ASSERT(m_mat->rows() == m_mat->cols(), "Matrix is not square.");
        if (!m_precond) m_precond = gsIdentityOp<T>::make(m_mat->rows());
    }

    virtual ~gsIterativeSolver()    {}

    /// @brief Returns a list of default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addInt   ("MaxIterations"    , "Maximum number of iterations", 1000 );
        opt.addReal  ("Tolerance"        , "Numerical tolerance"         , 1e-10);
        return opt;
    }

    /// @brief Set the options based on a gsOptionList
    virtual void setOptions(const gsOptionList & opt)
    {
        m_max_iters        = opt.askInt   ("MaxIterations"    , m_max_iters        );
        m_tol              = opt.askReal  ("Tolerance"        , m_tol              );
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
        if (initIteration(rhs, x)) return;

        while (m_num_iter < m_max_iters)
        {
            m_num_iter++;
            if (step(x)) break;
        }

        finalizeIteration(x);

    }

    /// @brief Solves the linear system and stores the solution in \a x and records
    /// the error histroy.
    ///
    /// Solves the linear system of equations
    /// \param[in]     rhs              the right hand side of the linear system
    /// \param[in,out] x                starting value; the solution is stored in here
    /// \param[out]    error_history    the error history is stored here
    ///
    /// \ingroup Solver  
    void solveDetailed( const VectorType& rhs, VectorType& x, VectorType& error_history )
    {
        if (initIteration(rhs, x))
        {
            error_history.resize(1);
            error_history[0] = m_error;
            return;
        }

        std::vector<T> tmp_error_hist;
        
        tmp_error_hist.clear();
        tmp_error_hist.reserve(m_max_iters / 3 );
        tmp_error_hist.push_back(m_error);        // store initial error (as provided by initIteration)

        while (m_num_iter < m_max_iters)
        {
            m_num_iter++;

            if (step(x)) break;

            tmp_error_hist.push_back(m_error);   // store initial error (as provided by step)
        }

        finalizeIteration(x);

        // move the error history to output variable
        error_history.swap( gsAsVector<T>(tmp_error_hist) );       
    }

    /// Init the iteration
    virtual bool initIteration( const VectorType& rhs, VectorType& x )
    {
        GISMO_ASSERT( rhs.cols() == 1,
                      "Iterative solvers only work for single column right hand side." );
     
        GISMO_ASSERT( m_precond->rows() == m_mat->rows(),
                      "The preconditionner does not match the matrix." );

        GISMO_ASSERT( m_precond->cols() == m_mat->cols(),
                      "The preconditionner does not match the matrix." );
        
        m_num_iter = 0;
        
        m_rhs_norm = rhs.norm();

        if (0 == m_rhs_norm) // special case of zero rhs
        {
            x.setZero(rhs.rows(),rhs.cols()); // for sure zero is a solution
            m_error = 0.;
            return true; // iteration is finished
        }
        
        if ( 0 == x.size() ) // if no initial solution, start with zeros
            x.setZero(rhs.rows(), rhs.cols() );
        else
        {
           GISMO_ENSURE(m_mat->cols() == x.rows(), "Invalid initial solution");
           GISMO_ENSURE(rhs.cols() == x.cols()   , "Initial solution does not match right-hand side");
        }
        return false; // iteration is not finished
    }

    virtual bool step( VectorType& x ) = 0;                     ///< Perform one step, requires initIteration
    virtual void finalizeIteration( VectorType& x ) {}          ///< Some post-processing might be required

    /// Returns the size of the linear system
    index_t size() const                                       { return m_mat->rows(); }
   
    /// Set the preconditionner
    void setPreconditioner(const LinOpPtr & precond)           { m_precond = precond; }

    /// Set the maximum number of iterations (default: 1000)
    void setMaxIterations(index_t max_iters)                   { m_max_iters = max_iters; }

    /// Set the tolerance for the error criteria (default: 1e-10)
    void setTolerance(T tol)                                   { m_tol = tol; }

    /// The number of iterations needed to reach the error criteria
    int iterations() const                                     { return m_num_iter; }

    /// The error of the iterative method
    T error() const                                            { return m_error; }

    /// The tolerance used in the iterative method
    T tolerance() const                                        { return m_tol; }

protected:
    const LinOpPtr m_mat;             ///< The matrix/operator to be solved for
    LinOpPtr       m_precond;         ///< The preconditioner
    index_t        m_max_iters;       ///< The upper bound for the number of iterations
    T              m_tol;             ///< The tolerance for m_error to be reached
    index_t        m_num_iter;        ///< The number of iterations performed
    T              m_rhs_norm;        ///< The norm of the right-hand-side
    T              m_error;           ///< The relative error as absolute_error/m_rhs_norm
};

} // namespace gismo
