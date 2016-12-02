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
#include <gsIO/gsOptionList.h>

namespace gismo
{

namespace trilinos
{

enum BelosSolverMode
{
    BiCGStab                =  1, ///< BiCGStab solver
    BlockCG                 =  2, ///< Block GC solver
/// BlockGCRODR             =  3, ///< Block Recycling GMRES solver
    BlockGmres              =  4, ///< Block GMRES solver
    FixedPoint              =  5, ///< Fixed-point solver
    GCRODR                  =  6, ///< Recycling GMRES solver
    GmresPoly               =  7, ///< Hybrid block GMRES solver
    LSQR                    =  8, ///< LSQR solver
    Minres                  =  9, ///< Minres solver
    PCPG                    = 10, ///< PCPG solver
    PseudoBlockCG           = 11, ///< Pseudo Block CG solver
    PseudoBlockGmres        = 12, ///< Pseudo Block GMRES solver
    PseudoBlockStochasticCG = 13, ///< Pseudo Block CG solver
    PseudoBlockTFQMR        = 14, ///< Pseudo Block TFQMR solver
    RCG                     = 15, ///< RCG solver
    TFQMR                   = 16  ///< TFQMR solver
    
};

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
    
     AbstractSolver();
     explicit AbstractSolver(const SparseMatrix & A);

    ~AbstractSolver();
    
    const Vector & solve(const Vector & b);

    void getSolution(gsVector<real_t> & sol, const int rank = 0) const;
    
protected:
        virtual void solveProblem() = 0;
    
protected:
        AbstractSolverPrivate        * my;
};

/** 
    
 */
class GISMO_EXPORT GMRES : public AbstractSolver
{
public:
    typedef AbstractSolver Base;
    
public:
    
    explicit GMRES(const SparseMatrix & A)
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
    
    explicit KLU(const SparseMatrix & A)
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
    
    explicit SuperLU(const SparseMatrix & A)
    : Base(A)
    { }
    
private:
    
    void solveProblem();
};

struct BelosSolverPrivate;

template<int mode>
class BelosSolver : public AbstractSolver
{
public:
    typedef AbstractSolver Base;

    explicit BelosSolver(const SparseMatrix & A);

    ~BelosSolver();

    /// Blocksize to be used by iterative solver
    void set(const std::string & name, const int & value);
    void set(const std::string & name, const double & value);
    void set(const std::string & name, const std::string & value);

    /// Hermitian problem type used to be by iterative solver
    /// : if matrix is symmetric, specifies it in the linear problem
    void setHermitian();

    std::string printValidParams() const;
    int numIterations() const;

private:

    void solveProblem();

private:

    BelosSolverPrivate * myBelos;

    int maxiters;  // maximum number of iterations allowed per linear system
};


};// namespace solver
};// namespace trilinos
};// namespace gismo


