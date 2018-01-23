/** @file gsNewtonIterator.h

    @brief Performs Newton iterations to solve a nonlinear equation system.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Goyal, A. Mantzaflaris, O. Weeger
*/


#pragma once

#include <gsAssembler/gsAssembler.h>


namespace gismo
{

/** 
    @brief Performs Newton iterations to solve a nonlinear system of PDEs.
    
    \tparam T coefficient type
    
    \ingroup Pde
*/
template <class T> 
class gsNewtonIterator
{
public:

    /// Constructor giving access to the gsAssemblerBase object to
    /// create a linear system per iteration
    gsNewtonIterator(  gsAssembler<T> & assembler,
                       const gsMultiPatch<T> & initialSolution)
    : m_assembler(assembler),
      m_curSolution(initialSolution),
      m_numIterations(0),
      m_maxIterations(100),
      m_tolerance(1e-12),
      m_converged(false)
    { 

    }

    gsNewtonIterator(gsAssembler<T> & assembler)
    : m_assembler(assembler),
      m_numIterations(0),
      m_maxIterations(100),
      m_tolerance(1e-12),
      m_converged(false)
    { 

    }


public:

    /// \brief Applies Newton method and Performs Newton iterations
    /// until convergence or maximum iterations.
    void solve();

    /// \brief Solves linear system obtained using linear elasticity
    /// in first step and computes residual
    void firstIteration();

    /// \brief Solves linear system in each iteration based on last
    /// solution and computes residual
    void nextIteration();

public:

    /// \brief Returns the latest configuration
    const gsMultiPatch<T> & solution() const { return m_curSolution; }

    /// \brief Tells whether the Newton method converged
    bool converged() const {return m_converged;}

    /// \brief Returns the number of Newton iterations performed
    index_t numIterations() const { return m_numIterations; }

    /// \brief Returns the tolerance value used
    T tolerance() const {return m_tolerance;}

    /// \brief Returns the error after solving the nonlinear system
    T residue()   const {return m_residue;}

    /// \brief Set the maximum number of Newton iterations allowed
    void setMaxIterations(index_t nIter) {m_maxIterations = nIter;}

    /// \brief Set the tolerance for convergence
    void setTolerance(T tol) {m_tolerance = tol;}

protected:

    virtual void solveLinearProblem(gsMatrix<T> &updateVector);

    virtual void solveLinearProblem(const gsMultiPatch<T> & currentSol, gsMatrix<T> &updateVector);

    virtual T getResidue() {return m_assembler.rhs().norm();}
protected:

    /// \brief gsAssemblerBase object to generate the linear system
    /// for each iteration
    gsAssembler<T> & m_assembler;

    /// \brief The latest/current solution
    gsMultiPatch<T>     m_curSolution;

    /// \brief Solution of the linear system in each iteration
    gsMatrix<T>         m_updateVector;

    /// Linear solver employed
    //gsSparseSolver<>::LU  m_solver;
    //typename gsSparseSolver<T>::BiCGSTABDiagonal m_solver;
    //typename gsSparseSolver<>::CGDiagonal m_solver;
    gsSparseSolver<>::LU  m_solver;

protected:

    /// \brief Number of Newton iterations performed
    index_t m_numIterations;

    /// \brief Maximum number of Newton iterations allowed
    index_t m_maxIterations;

    /// \brief Tolerance value to decide convergence
    T       m_tolerance;

protected:

    /// \brief Convergence result
    bool m_converged;

    /// \brief Final error
    T m_residue;

    /// \brief Norm of the current Newton update vector
	T m_updnorm;

};


} // namespace gismo


namespace gismo
{

template <class T>
void gsNewtonIterator<T>::solveLinearProblem(gsMatrix<T>& updateVector)
{
    // Construct the linear system
    m_assembler.assemble();

    // gsDebugVar( m_assembler.matrix().toDense() );
    // gsDebugVar( m_assembler.rhs().transpose() );

    // Compute the newton update
    m_solver.compute( m_assembler.matrix() );
    updateVector = m_solver.solve( m_assembler.rhs() );
    
    // gsDebugVar(updateVector);
}

template <class T>
void gsNewtonIterator<T>::solveLinearProblem(const gsMultiPatch<T> & currentSol, gsMatrix<T>& updateVector)
{
    // Construct linear system for next iteration
    m_assembler.assemble(currentSol);

    // gsDebugVar( m_assembler.matrix().toDense() );
    // gsDebugVar( m_assembler.rhs().transpose() );
    
    // Compute the newton update
    m_solver.compute( m_assembler.matrix() );
    updateVector = m_solver.solve( m_assembler.rhs() );

    // gsDebugVar(updateVector);
}


template <class T> 
void gsNewtonIterator<T>::solve()
{
    firstIteration();

    const T initResidue = m_residue;
	const T initUpdate = m_updnorm;

    // ----- Iterations start -----
    for (m_numIterations = 1; m_numIterations < m_maxIterations; ++m_numIterations)
    {
        gsDebug << "Newton iteration " << m_numIterations 
            << " residue " << math::abs(m_residue)
            << ".\n";
        nextIteration();
        
        // termination criteria
        if ( math::abs(m_updnorm / initUpdate)  < m_tolerance ||
             math::abs(m_residue / initResidue) < m_tolerance )
        {
            m_converged = true;
            break;
        }
    }
}


template <class T> 
void gsNewtonIterator<T>::firstIteration()
{
    // ----- First iteration -----
    m_converged = false;

    // Solve 
    solveLinearProblem(m_updateVector);

    // Construct initial solution
    m_assembler.constructSolution(m_updateVector, m_curSolution);

    // Homogenize Dirichlet dofs (values are now copied in m_curSolution)
    m_assembler.homogenizeFixedDofs(-1);

    // Compute initial residue
    m_residue = getResidue();
    m_updnorm = m_updateVector   .norm();

	gsDebug<<"Iteration: "<< 0
               <<", residue: "<< m_residue
               <<", update norm: "<< m_updnorm
               <<"\n";
}

template <class T> 
void gsNewtonIterator<T>::nextIteration()
{
    // Solve the linaer system of the current iteration
    solveLinearProblem(m_curSolution, m_updateVector);

    // Update the deformed solution
    m_assembler.updateSolution(m_updateVector, m_curSolution);
    
    // Compute residue
    m_residue = getResidue();
    m_updnorm = m_updateVector.norm();
    
    gsDebug<<"Iteration: "<< m_numIterations
           <<", residue: "<< m_residue
           <<", update norm: "<< m_updnorm
           <<"\n";
}



} // namespace gismo

