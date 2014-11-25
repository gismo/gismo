//
// gsSmoother.cpp
//
// Multigrid smoothers.
//
// Clemens Hofreither
//
 
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

    Eigen::ConjugateGradient< gsSparseMatrix<real_t> > solver;
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
    
    corr = Eigen::ConjugateGradient< gsSparseMatrix<real_t> >( P ).setTolerance(1e-3).solve(corr).eval();
    corr.array() *= P.diagonal().array();
#else
    // symmetricized version
    gsMatrix<real_t> corr;
    corr = f - A * x;

    corr.array() /= A.diagonal().array().sqrt();
    corr.array() *= P.diagonal().array().sqrt();

    corr = Eigen::ConjugateGradient< gsSparseMatrix<real_t> >( P ).setTolerance(1e-3).solve(corr).eval();

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

    
    corr = Eigen::ConjugateGradient< gsSparseMatrix<real_t> >( P )
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


void gsDampedPrecRichardsonSmoother::apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
{
    assert( A.rows() == x.rows() && x.rows() == f.rows() );
    assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);

    gsMatrix<real_t> corr;
    corr = f - A * x;
    corr = m_solver.solve(corr).eval();
    ////corr.array() *= P.diagonal().array();
    x += m_damping * corr;
}



void gsGaussSeidelSmoother::apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
{
    gaussSeidelSweep(A, x, f);
}

void gsGaussSeidelSmoother::applyT(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
{
    reverseGaussSeidelSweep(A, x, f);
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

