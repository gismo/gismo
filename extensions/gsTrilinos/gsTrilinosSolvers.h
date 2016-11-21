/** @file gsTrilinosSolvers.h

    @brief Wrappers for Trilinos solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsTrilinos/Vector.h>

namespace gismo
{

namespace trilinos
{

/** @namespace gismo::trilinos::solver

    @brief This namespace contains wrappers for Trilinos linear
    system solvers and eigenvalue solvers
*/
namespace solver
{

struct AbstractSolverPrivate;
struct AbstractSolverBelosPrivate;

/** 
    
 */
class GISMO_EXPORT AbstractSolver
{
public:
    
     AbstractSolver();
     AbstractSolver(const SparseMatrix & A, const Vector & b);
     AbstractSolver(const SparseMatrix & A, const Vector & b, 
                    const std::string solver_teuchosUser);

    ~AbstractSolver();
    
    const Vector & solve();

    void getSolution(gsVector<real_t> & sol, const int rank = 0) const;
    
protected:
        virtual void solveProblem() = 0;
    
protected:
        AbstractSolverPrivate        * my = NULL;
        AbstractSolverBelosPrivate   * myBelos = NULL;
        std::string SolverTeuchosUser = "";
};

/** 
    
 */
class GISMO_EXPORT GMRES : public AbstractSolver
{
public:
    typedef AbstractSolver Base;
    
public:
    
    GMRES(const SparseMatrix & A, const Vector & b)
    : Base(A, b), m_tolerance(10e-6), m_maxIter(50)
    { }
    
private:

    double m_tolerance;
    int    m_maxIter;

private:
    
    void solveProblem();
};

/** 
    
 */
class GISMO_EXPORT KLU : public AbstractSolver
{
public:
    typedef AbstractSolver Base;
    
public:
    
    KLU(const SparseMatrix & A, const Vector & b)
    : Base(A, b)
    { }
    
private:
    
    void solveProblem();
};


class GISMO_EXPORT SuperLU : public AbstractSolver
{
public:
    typedef AbstractSolver Base;
    
public:
    
    SuperLU(const SparseMatrix & A, const Vector & b)
    : Base(A, b)
    { }
    
private:
    
    void solveProblem();
};

class GISMO_EXPORT Belos_solver : public AbstractSolver
{
public:
    typedef AbstractSolver Base;

    Belos_solver(const SparseMatrix & A, const Vector & b, 
                 const std::string solver_teuchosUser = "Belos")
                 : Base(A, b, solver_teuchosUser), blocksize(1), maxiters(500)
    { }

private:

    int blocksize;   // blocksize
    int maxiters;  // maximum number of iterations allowed per linear system

    void solveProblem();

};


};// namespace solver
};// namespace trilinos
};// namespace gismo


