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
template<class T=real_t>
class gsIterativeSolver
{
public:
    typedef memory::shared_ptr<gsIterativeSolver> Ptr;
    typedef memory::unique_ptr<gsIterativeSolver> uPtr;
    typedef T                                     ScalarType;
    typedef gsMatrix<T>                           VectorType;
    typedef typename gsLinearOperator<T>::Ptr     LinOpPtr;

    /// @brief Contructor using a linear operator to be solved for and
    ///  a preconditioner
    ///
    /// @param mat     The operator to be solved for as shared pointer to a gsLinearOperator
    /// @param precond The preconditioneras shared pointer to a gsLinearOperator,
    ///                a null pointer is defaulted to the identity
    gsIterativeSolver( const LinOpPtr& mat,
                       const LinOpPtr& precond)
    : m_mat(mat),
      m_precond(precond),
      m_max_iters(1000),
      m_tol(1e-10),
      m_num_iter(-1),
      m_rhs_norm(-1),
      m_error(-1)
    {
        GISMO_ASSERT(m_mat->rows()     == m_mat->cols(),     "The matrix is not square."                     );

        if (!m_precond) m_precond = gsIdentityOp<T>::make(m_mat->rows());
        GISMO_ASSERT(m_precond->rows() == m_precond->cols(), "The preconditioner is not square."             );
        GISMO_ASSERT(m_precond->rows() == m_mat->rows(),     "The preconditioner does not match the matrix: "
                                                             <<m_precond->rows()<<"!="<<m_mat->rows()        );
    }

    /// @brief Contructor using any dense or sparse matrix and a
    ///  preconditioner
    ///
    /// @param mat     The operator to be solved for as a (reference) to a matrix.
    /// @param precond The preconditioneras shared pointer to a gsLinearOperator,
    ///                a null pointer is defaulted to the identity
    ///
    /// @note This does not copy the matrix in \a mat. So, make sure that the
    /// matrix is not deleted before the solver. If you have a shared pointer to
    /// the matrix, you might use \ref makeMatrixOp() to obtain a shared pointer
    /// to a gsLinearOperator which can be supplied alternatively.
    template<typename Derived>
    gsIterativeSolver( const Eigen::EigenBase<Derived> & mat,
                       const LinOpPtr& precond)
    : m_mat(makeMatrixOp(mat.derived())),
      m_precond(precond),
      m_max_iters(1000),
      m_tol(1e-10),
      m_num_iter(-1),
      m_rhs_norm(-1),
      m_error(-1)
    {
        GISMO_ASSERT(m_mat->rows()     == m_mat->cols(),     "The matrix is not square."                     );

        if (!m_precond) m_precond = gsIdentityOp<T>::make(m_mat->rows());
        GISMO_ASSERT(m_precond->rows() == m_precond->cols(), "The preconditioner is not square."             );
        GISMO_ASSERT(m_precond->rows() == m_mat->rows(),     "The preconditioner does not match the matrix: "
                                                             <<m_precond->rows()<<"!="<<m_mat->rows()        );
    }

    virtual ~gsIterativeSolver()    {}

    /// @brief Returns a list of default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addInt   ("MaxIterations"    , "Maximum number of iterations", 1000       );
        opt.addReal  ("Tolerance"        , "Tolerance for the error criteria on the "
                                           "relative residual error",      1e-10      );
        return opt;
    }

    /// @brief Set the options based on a gsOptionList
    virtual gsIterativeSolver& setOptions(const gsOptionList & opt)
    {
        m_max_iters        = opt.askInt   ("MaxIterations"    , m_max_iters        );
        m_tol              = opt.askReal  ("Tolerance"        , m_tol              );
        return *this;
    }

    /// @brief Solves the linear system and stores the solution in \a x
    ///
    /// Solves the linear system of equations
    /// @param[in]     rhs      the right hand side of the linear system
    /// @param[in,out] x        starting value; the solution is stored in here
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
    /// @param[in]     rhs              the right hand side of the linear system
    /// @param[in,out] x                starting value; the solution is stored in here
    /// @param[out]    error_history    the error history is stored here
    void solveDetailed( const VectorType& rhs, VectorType& x, VectorType& error_history )
    {
        if (initIteration(rhs, x))
        {
            error_history.resize(1,1); //VectorType is actually gsMatrix
            error_history(0,0) = m_error;
            //gsDebug<<"Solution reached at iteration 0, err="<<m_error<<"\n";
            return;
        }

        std::vector<T> tmp_error_hist;
        tmp_error_hist.clear();
        tmp_error_hist.reserve(m_max_iters / 3 );
        tmp_error_hist.push_back(m_error); // store initial error (as provided by initIteration)

        while (m_num_iter < m_max_iters)
        {
            m_num_iter++;
            //gsDebug<<"Iteration : "<<std::setw(5)<<std::left<< m_num_iter<<"\n";

            if (step(x))
            {
                tmp_error_hist.push_back(m_error);
                //gsDebug<<"            err = "<<m_error<<" --> Solution reached.\n";
                break;
            }

            tmp_error_hist.push_back(m_error);
            //gsDebug<<"            err = "<<m_error<<"\n";
        }

        //if (m_num_iter == m_max_iters) gsDebug<<"Maximum number of iterations reached.\n";

        finalizeIteration(x);
        error_history = gsAsVector<T>(tmp_error_hist);
    }

    /// Init the iteration
    virtual bool initIteration( const VectorType& rhs, VectorType& x )
    {
        GISMO_ASSERT( rhs.cols() == 1,
                      "Iterative solvers only work for single column right-hand side." );
        GISMO_ASSERT( rhs.rows() == m_mat->rows(),
                      "The right-hand side does not match the matrix: "
                      << rhs.rows() <<"!="<< m_mat->rows() );

        m_num_iter = 0;

        m_rhs_norm = rhs.norm();

        if (0 == m_rhs_norm) // special case of zero rhs
        {
            x.setZero(rhs.rows(),rhs.cols()); // for sure zero is a solution
            m_error = 0.;
            return true; // iteration is finished
        }

        if ( 0 == x.size() ) // if no initial solution, start with zeros
            x.setZero(rhs.rows(), rhs.cols());
        else
        {
            GISMO_ASSERT( x.cols() == 1,
                      "Iterative solvers only work for single right-hand side and solution." );
            GISMO_ASSERT( x.rows() == m_mat->cols(),
                      "The initial guess does not match the matrix: "
                      << x.rows() <<"!="<< m_mat->cols() );
        }
        return false; // iteration is not finished
    }

    virtual bool step( VectorType& x ) = 0;                   ///< Perform one step, requires initIteration
    virtual void finalizeIteration( VectorType& ) {}          ///< Some post-processing might be required

    /// Returns the size of the linear system
    index_t size() const                                       { return m_mat->rows(); }

    /// Set the preconditionner
    void setPreconditioner(const LinOpPtr & precond)           { m_precond = precond; }

    /// Get the preconditioner
    LinOpPtr preconditioner() const                            { return m_precond; }

    /// Get the underlying matrix/operator to be solved for
    LinOpPtr underlying() const                                { return m_mat; }

    /// Set the maximum number of iterations (default: 1000)
    void setMaxIterations(index_t max_iters)                   { m_max_iters = max_iters; }

    /// Set the tolerance for the error criteria on the relative residual error (default: 1e-10)
    void setTolerance(T tol)                                   { m_tol = tol; }

    /// The number of iterations needed to reach the error criteria
    index_t iterations() const                                 { return m_num_iter; }

    /// @brief The relative residual error of the current iterate
    ///
    /// This is the Euclidean norm of the residual, devided by the Euclidean
    /// norm of the right-hand side.
    T error() const                                            { return m_error; }

    /// The chosen tolerance for the error criteria on the relative residual error
    T tolerance() const                                        { return m_tol; }

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const = 0;

    /// Prints the object as a string with extended details.
    virtual std::string detail() const
    {
        std::ostringstream os;
        print(os);
        os << " Tolerance            : " << tolerance() << "\n";
        os << " Solver error         : " << error() << "\n";
        os << " Number of iterations : " << iterations() << " (max="<<m_max_iters<<")\n";
        return os.str();
    }

protected:
    const LinOpPtr m_mat;             ///< The matrix/operator to be solved for
    LinOpPtr       m_precond;         ///< The preconditioner
    index_t        m_max_iters;       ///< The upper bound for the number of iterations
    T              m_tol;             ///< The tolerance for m_error to be reached
    index_t        m_num_iter;        ///< The number of iterations performed
    T              m_rhs_norm;        ///< The norm of the right-hand-side
    T              m_error;           ///< The relative error as absolute_error/m_rhs_norm
};

/// \brief Print (as string) operator for iterative solvers
/// \relates gsIterativeSolver
template<class T>
std::ostream &operator<<(std::ostream &os, const gsIterativeSolver<T>& b)
{return b.print(os); }

/// \brief This wrapper class allows \a gsIterativeSolver to be used as \a gsLinearOperator
template <class SolverType>
class gsIterativeSolverOp GISMO_FINAL : public gsLinearOperator<typename SolverType::ScalarType>
{
    typedef typename SolverType::ScalarType T;
    typedef typename gsLinearOperator<T>::Ptr LinOpPtr;
public:
    /// Shared pointer for gsIterativeSolverOp
    typedef memory::shared_ptr<gsIterativeSolverOp> Ptr;

    /// Unique pointer for gsIterativeSolverOp
    typedef memory::unique_ptr<gsIterativeSolverOp> uPtr;

    /// Constructor taking the underlying matrix/operator and the preconditioner
    template<class MatrixType>
    gsIterativeSolverOp(const MatrixType& matrix, const LinOpPtr& preconditioner = LinOpPtr())
        : m_solver(matrix, preconditioner) {}

    /// Make function taking the underlying matrix/operator and the preconditioner
    template<class MatrixType>
    static uPtr make(const MatrixType& matrix, const LinOpPtr& preconditioner = LinOpPtr())
    { return uPtr( new gsIterativeSolverOp(matrix,preconditioner) ); }

    /// Make function taking a matrix OR a shared pointer
    static uPtr make(SolverType solver) { return memory::make_unique( new gsIterativeSolverOp(solver) ); }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & res) const
    {
        res.setZero(input.rows(), input.cols());
        m_solver.solve(input, res);
        gsDebug << ( (m_solver.error() < m_solver.tolerance())
                     ? "gsIterativeSolverOp reached desired error bound after "
                     : "msIterativeSolverOp did not reach desired error bound within " )
                << m_solver.iterations() << " iterations.\n";
        gsDebug << input.rows() << "; " << input.norm() << "; " << res.norm() << "\n";
    }

    index_t rows() const { return m_solver.underlying()->rows(); }

    index_t cols() const { return m_solver.underlying()->cols(); }

    /// Access the solver class
    SolverType& solver()                { return m_solver; }

    /// Access the solver class
    const SolverType& solver() const    { return m_solver; }

private:
    mutable SolverType m_solver;  // The iterative solvers change their state during solving
};


} // namespace gismo
