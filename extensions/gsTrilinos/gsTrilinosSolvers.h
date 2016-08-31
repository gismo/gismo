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

/** 
    
 */
class GISMO_EXPORT AbstractSolver
{
public:
    
    AbstractSolver(const SparseMatrix & A );

    ~AbstractSolver();
    
    const Vector & solve( const Vector & b );

    void getSolution(gsVector<> & sol, const int rank = 0) const;
    
protected:
        virtual void solveProblem() = 0;
    
protected:
        AbstractSolverPrivate * my;
};


/** 
    
 */
class GISMO_EXPORT GMRES : public AbstractSolver
{
public:
    typedef AbstractSolver Base;
    
public:
    
    GMRES(const SparseMatrix & A )
    : Base(A), m_tolerance(10e-6), m_maxIter(50)
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
    
    KLU(const SparseMatrix & A )
    : Base(A)
    { }
    
private:
    
    void solveProblem();
};


class GISMO_EXPORT SuperLU : public AbstractSolver
{
public:
    typedef AbstractSolver Base;
    
public:
    
    SuperLU(const SparseMatrix & A )
    : Base(A)
    { }
    
private:
    
    void solveProblem();
};



};// namespace solver
};// namespace trilinos
};// namespace gismo


