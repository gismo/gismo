/** @file gsSmoother.cpp

    @brief Provides Multigrid smoothers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/
#include <gsSolver/gsSmoother.h>
#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{

// smoothers

void dampedRichardsonSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau = (real_t)(1))
{
    assert( A.rows() == x.rows() && x.rows() == f.rows() );
    assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);

    gsMatrix<real_t> temp = f - A * x;
    x += tau * temp;
}

void dampedRichardsonSweepBoundary(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau, int numBdNodes)
{
    assert( A.rows() == x.rows() && x.rows() == f.rows() );
    assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);

    gsMatrix<real_t> temp = f - A * x;
    for (int i = 0; i < A.outerSize(); ++i)
    {
        if (i < numBdNodes || i >= A.outerSize() - numBdNodes)
            x(i) += tau * temp(i);
    }
}

// assumes symmetric matrix
void kaczmarzSweepBoundary(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau, int numBdNodes)
{
    assert( A.rows() == x.rows() && x.rows() == f.rows() );
    assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);

    for (int eq = 0; eq < A.outerSize(); ++eq)
    {
        if (!(eq < numBdNodes || eq >= A.outerSize() - numBdNodes))
            continue;

        real_t Asq = 0;
        real_t Aeqx = 0;

        for (Eigen::SparseMatrix<real_t>::InnerIterator it(A,eq); it; ++it)
        {
            Asq += it.value() * it.value();         // compute squared norm of A_eq
            Aeqx += it.value() * x(it.index());     // compute A_eq . x
        }

        const real_t beta = (f(eq) - Aeqx) / Asq;

        for (Eigen::SparseMatrix<real_t>::InnerIterator it(A,eq); it; ++it)
        {
            x(it.index()) += beta * it.value();     // add beta * A_eq
        }
    }
}


void jacobiSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
{
    assert( A.rows() == x.rows() && x.rows() == f.rows() );
    assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);

    gsMatrix<real_t> temp = f - A * x;
    temp.array() /= A.diagonal().array();
    x += temp;
}

void dampedJacobiSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau = (real_t)(1.0/2.0))
{
    assert( A.rows() == x.rows() && x.rows() == f.rows() );
    assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);

    gsMatrix<real_t> temp = f - A * x;
    temp.array() /= A.diagonal().array();
    x += tau * temp;
}


void dampedPreRichardsonSweep(const Eigen::SparseMatrix<real_t>& A, const Eigen::SparseMatrix<real_t>& P, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau)
{
    assert( A.rows() == x.rows() && x.rows() == f.rows() );
    assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);

    gsSparseSolver<>::CGDiagonal solver;
    gsMatrix<real_t> corr;
    corr = f - A * x;
    corr = solver.compute( P ).setTolerance(1e-3).solve(corr).eval();
    //corr.array() *= P.diagonal().array();
    x += tau * corr;
}


void dampedPreJacobiSweep(const Eigen::SparseMatrix<real_t>& A, const Eigen::SparseMatrix<real_t>& P, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau = (real_t)(1.0/2.0))
{
    assert( A.rows() == x.rows() && x.rows() == f.rows() );
    assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);

    //gsMatrix<real_t> temp = f - A * x;
    //temp.array() /= A.diagonal().array();
    //x = tau * temp;

#if 0
    // original version
    gsMatrix<real_t> corr;
    corr = f - A * x;
    corr.array() /= A.diagonal().array();
    
    corr = gsSparseSolver<>::CGDiagonal( P ).setTolerance(1e-3).solve(corr).eval();
    corr.array() *= P.diagonal().array();
#else
    // symmetricized version
    gsMatrix<real_t> corr;
    corr = f - A * x;

    corr.array() /= A.diagonal().array().sqrt();
    corr.array() *= P.diagonal().array().sqrt();

    corr = gsSparseSolver<>::CGDiagonal( P ).setTolerance(1e-3).solve(corr).eval();

    corr.array() *= P.diagonal().array().sqrt();
    corr.array() /= A.diagonal().array().sqrt();
#endif

    x += tau * corr;
}


void gaussSeidelSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
{
    assert( A.rows() == x.rows() && x.rows() == f.rows() );
    assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);

    // A is supposed to be symmetric, so it doesn't matter if it's stored in row- or column-major order
    for (int i = 0; i < A.outerSize(); ++i)
    {
        real_t diag = real_t(0.0);
        real_t sum = real_t(0.0);

        for (Eigen::SparseMatrix<real_t>::InnerIterator it(A,i); it; ++it)
        {
            sum += it.value() * x( it.index() );        // compute A.x
            if (it.index() == i)
                diag = it.value();
        }

        x(i) += (f(i) - sum) / diag;
    }
}

void reverseGaussSeidelSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
{
    assert( A.rows() == x.rows() && x.rows() == f.rows() );
    assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);

    // A is supposed to be symmetric, so it doesn't matter if it's stored in row- or column-major order
    for (int i = A.outerSize() - 1; i >= 0; --i)
    {
        real_t diag = real_t(0.0);
        real_t sum = real_t(0.0);

        for (Eigen::SparseMatrix<real_t>::InnerIterator it(A,i); it; ++it)
        {
            sum += it.value() * x( it.index() );        // compute A.x
            if (it.index() == i)
                diag = it.value();
        }

        x(i) += (f(i) - sum) / diag;
    }
}

void gaussSeidelSingleBlock(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, gsVector<index_t>& DoFs)
{
    //Sorting from lowest to highest
    DoFs.sortByColumn(0);
    const index_t size = DoFs.rows();

    GISMO_ASSERT ( DoFs(size-1)< A.cols(), "The given DoF is higher than the size of the matrix");

    gsMatrix<real_t> Dblock(size, size);
    gsMatrix<real_t> residual(size, 1);
    for (int i = 0; i< size; i++)
    {
        residual(i,0) = f(DoFs(i),0) - (A.block(DoFs(i), 0, 1, A.cols())*x).value();
        for (int j = 0; j< size; j++)
            Dblock(i,j) = A.coeff(DoFs(i), DoFs(j));
    }
    //Multiply residual with the inverse of the diagonal block
    residual = Dblock.inverse()*residual;

    //Update
    for (int i = 0; i< size; i++)
       x(DoFs(i),0) += residual(i,0);
}

void preGaussSeidelSweep(const Eigen::SparseMatrix<real_t>& A, const Eigen::SparseMatrix<real_t>& P, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, real_t tau, bool reverse)
{
    assert( A.rows() == x.rows() && x.rows() == f.rows() );
    assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);

    gsMatrix<real_t> corr;
    corr.setZero( x.rows(), 1 );

    //const int start_i = reverse ? A.outerSize() - 1 : 0                 ;
    //const int end_i   = reverse ? 0                 : A.outerSize() - 1 ;
    //const int inc_i   = reverse ? -1                : 1                 ;

    // A is supposed to be symmetric, so it doesn't matter if it's stored in row- or column-major order
    for (int i = 0; i < A.outerSize(); ++i)
    {
        real_t diag = real_t(0.0);
        real_t sum = real_t(0.0);

        for (Eigen::SparseMatrix<real_t>::InnerIterator it(A,i); it; ++it)
        {
            sum += it.value() * (x( it.index() ) + corr( it.index() ));        // compute A.x
            if (it.index() == i)
                diag = it.value();
        }

        corr(i) = (f(i) - sum) / diag;
    }

    
    corr = gsSparseSolver<>::CGDiagonal( P )
            .setTolerance(1e-3).solve(corr).eval();
    x.array() += tau * /*P.diagonal().array() * */ corr.array();
}


////////////////////////////////////////////////////////////////////////////////
// smoother class implementations
////////////////////////////////////////////////////////////////////////////////


void gsRichardsonSmoother::apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
{
    dampedRichardsonSweep(A, x, f, m_damping);
}


void gsJacobiSmoother::apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
{
    dampedJacobiSweep(A, x, f, m_damping);
}


void gsGaussSeidelSmoother::apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
{
    gaussSeidelSweep(A, x, f);
}

void gsGaussSeidelSmoother::applyT(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
{
    reverseGaussSeidelSweep(A, x, f);
}

void gsGaussSeidelBlockSmoother::apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
{
    assert( A.rows() == x.rows() && x.rows() == f.rows() );
    assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);
    assert( m_blockInfo.size() > 0);

    const index_t numberOfBlocks = m_blockInfo.size();
    for (index_t k = 0; k < numberOfBlocks; k++)
        gaussSeidelSingleBlock(A, x, f, m_blockInfo[k]);
}

void gsGaussSeidelBlockSmoother::applyT(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
{
    assert( A.rows() == x.rows() && x.rows() == f.rows() );
    assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);
    assert( m_blockInfo.size() > 0);

    const index_t numberOfBlocks = m_blockInfo.size();
    for (index_t k = 0; k < numberOfBlocks; k++)
        gaussSeidelSingleBlock(A, x, f, m_blockInfo[numberOfBlocks - 1 - k]);
}


gsILUTSmoother::gsILUTSmoother(const Eigen::SparseMatrix<real_t>& K)
{
    m_ilu.setFillfactor(1);
    m_ilu.compute(K);
}

void gsILUTSmoother::apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
{
    x += m_ilu.solve( f - A * x ).eval();
}


} // namespace gismo

