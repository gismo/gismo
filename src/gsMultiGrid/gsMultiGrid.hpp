/** @file gsMultiGrid.hpp

    @brief Multigrid solver for isogeometric discretizations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#include <gsMultiGrid/gsMultiGrid.h>
#include <gsSolver/gsMatrixOp.h>

namespace gismo
{

template<class T>
gsMultiGridOp<T>::gsMultiGridOp(SpMatrix fineMatrix, std::vector< SpMatrixRowMajor > transferMatrices, OpPtr coarseSolver )
{
    const index_t sz = transferMatrices.size();
    std::vector<SpMatrixRowMajorPtr> transferMatrixPtrs(sz);
    for (index_t i=0; i<sz; ++i)
        transferMatrixPtrs[i] = transferMatrices[i].moveToPtr();

    init(fineMatrix.moveToPtr(),give(transferMatrixPtrs),coarseSolver);
}

template<class T>
gsMultiGridOp<T>::gsMultiGridOp(SpMatrixPtr fineMatrix, std::vector< SpMatrixRowMajorPtr > transferMatrices, OpPtr coarseSolver )
{
    init(give(fineMatrix),give(transferMatrices),give(coarseSolver));
}

template<class T>
gsMultiGridOp<T>::gsMultiGridOp( const std::vector<OpPtr>& ops, const std::vector<OpPtr>& prolongation,
                                          const std::vector<OpPtr>& restriction, OpPtr coarseSolver)
    : n_levels( ops.size() ), m_ops(ops), m_smoother(n_levels), m_prolong(prolongation), m_restrict(restriction),
      m_numPreSmooth(1), m_numPostSmooth(1), m_numCycles(1), m_damping(1)
{
    GISMO_ASSERT ( prolongation.size() == restriction.size(), "The number of prolongation and restriction operators differ." );
    GISMO_ASSERT ( ops.size() == prolongation.size()+1, "The number of prolongation and restriction operators do not fit to the number of operators." );

    if (coarseSolver)
        m_coarseSolver = coarseSolver;
    else
        initCoarseSolver();
}

template<class T>
void gsMultiGridOp<T>::init(SpMatrixPtr fineMatrix, std::vector<SpMatrixRowMajorPtr> transferMatrices, OpPtr coarseSolver)
{
    GISMO_ASSERT ( fineMatrix->rows() == fineMatrix->cols(), "gsMultiGridOp need quadratic matrices." );

    const index_t sz = transferMatrices.size();

    n_levels = sz+1;
    m_ops.resize(n_levels);
    m_smoother.resize(n_levels);
    m_prolong.resize(sz);
    m_restrict.resize(sz);
    m_numPreSmooth = 1;
    m_numPostSmooth = 1;
    m_numCycles = 1;
    m_damping = 1;

    SpMatrixPtr mat = fineMatrix;
    m_ops[n_levels-1] = makeMatrixOp(mat);

    for ( index_t i = n_levels - 2; i >= 0; --i )
    {
        SpMatrixPtr newMat = SpMatrixPtr(new SpMatrix(
            transferMatrices[i]->transpose() * *mat * *(transferMatrices[i])
        ));
        m_ops[i] = makeMatrixOp(newMat);
        mat = newMat; // copies just the smart pointers
    }

    for ( index_t i=0; i<sz; ++i )
    {
        m_prolong[i] = makeMatrixOp(transferMatrices[i]);
        m_restrict[i] = makeMatrixOp(transferMatrices[i]->transpose()); // note that this works as we know that
                                                           // m_prolong does not get destroyed before
                                                           // m_restrict, so the shared pointer stored in m_prolong
                                                           // will make sure that the matrix will not get destroyed
    }

    if (coarseSolver)
        m_coarseSolver = coarseSolver;
    else
        gsMultiGridOp<T>::initCoarseSolver();

}

template<class T>
void gsMultiGridOp<T>::initCoarseSolver()
{
    if (n_levels > 1)
    {
        const gsMatrixOp<SpMatrix>* matrOp = dynamic_cast< const gsMatrixOp<SpMatrix>* >( m_ops[0].get() );
        if (matrOp)
        {
            const SpMatrix & matr = matrOp->matrix();
            m_coarseSolver = makeSparseLUSolver(matr);
        }
        else
        {
            // Fallback for other types of operators (matrix free implementations, etc.)
            // Warn if we do this for big matrices...
            if (m_ops[0]->rows() > 50)
                gsWarn << "gsMultiGridOp::initCoarseSolver(): The coarse grid solver is constructed based on "
                    "gsLinearOperator::toMatrix(). This might be inefficient. Consider providing matrices of type "
                    "gsSparseMatrix<T> or an exact solver for the coarset grid level to gsMultiGridOp constructor.\n";
            gsMatrix<T> coarse_dense;
            m_ops[0]->toMatrix( coarse_dense );
            m_coarseSolver = makePartialPivLUSolver( coarse_dense );
        }
    }
}

template<class T>
void gsMultiGridOp<T>::smoothingStep(index_t level, const gsMatrix<T>& rhs, gsMatrix<T>& x) const
{
    GISMO_ASSERT (m_smoother[level], "Smoother is not defined. Define it using setSmoother." );

    // pre-smooth
    for (index_t i = 0; i < m_numPreSmooth; ++i)
    {
        m_smoother[level]->step( rhs, x );
    }

    // post-smooth
    for (index_t i = 0; i < m_numPostSmooth; ++i)
    {
        m_smoother[level]->stepT( rhs, x );
    }

}

template<class T>
void gsMultiGridOp<T>::multiGridStep(index_t level, const gsMatrix<T>& rhs, gsMatrix<T>& x) const
{
    GISMO_ASSERT ( 0 <= level && level < n_levels, "The given level is not feasible." );
    GISMO_ASSERT ( n_levels > 1, "Multigrid is only available if at least two grids are present. Use smoothingStep for running the smoother only." );

    if (level == 0)
    {
        solveCoarse(rhs, x);
    }
    else
    {
        const index_t lf = level;
        const index_t lc = lf - 1;

        GISMO_ASSERT (m_smoother[lf], "Smoother is not defined. Define it using setSmoother." );

        gsMatrix<T> fineRes, fineCorr, coarseRes, coarseCorr;

        // pre-smooth
        for (index_t i = 0; i < m_numPreSmooth; ++i)
        {
            m_smoother[lf]->step( rhs, x );
        }

        // compute fine residual
        m_ops[lf]->apply( x, fineRes );
        fineRes -= rhs;

        // restrict residual to coarse grid
        restrictVector( lf, fineRes, coarseRes );

        // obtain coarse-grid correction by recursing
        coarseCorr.setZero( nDofs(lc), coarseRes.cols() );
        for (index_t i = 0; i < ((lc==0 && m_numCycles>0) ? 1 : m_numCycles); ++i)      // coarse solve is never cycled
        {
            multiGridStep( lc, coarseRes, coarseCorr );
        }

        // prolong correction
        prolongVector( lc, coarseCorr, fineCorr );

        // apply correction
        x -= m_damping * fineCorr;

        // post-smooth
        for (index_t i = 0; i < m_numPostSmooth; ++i)
        {
            m_smoother[lf]->stepT( rhs, x );
        }
    }
}

template<class T>
void gsMultiGridOp<T>::fullMultiGrid(const std::vector< gsMatrix<T> >& rhs, const std::vector< gsMatrix<T> >& dirichletIntp, gsMatrix<T>& result) const
{

    GISMO_ASSERT (dirichletIntp.size() == (unsigned)n_levels, "The number of dirichletIntp does not correspond to the number of levels!");
    GISMO_ASSERT (rhs.size() == (unsigned)n_levels, "The number of rhs does not correspond to the number of levels!");

    std::vector< gsMatrix<T> > u(n_levels);

    // solve coarse problem
    solveCoarse( rhs[0], u[0] );

    for (index_t i = 1; i <= finestLevel(); ++i)
    {
        // transfer result from previous level
        prolongVector(i-1, u[i-1], u[i]);
        u[i] += dirichletIntp[i];

        // run one multigrid step
        multiGridStep(i, rhs[i], u[i]);
    }

    result = u[ finestLevel() ];
}

template<class T>
void gsMultiGridOp<T>::cascadicMultiGrid(const std::vector< gsMatrix<T> >& rhs, const std::vector< gsMatrix<T> >& dirichletIntp, gsMatrix<T>& result) const
{

    GISMO_ASSERT (dirichletIntp.size() == (unsigned)n_levels, "The number of dirichletIntp does not correspond to the number of levels!");
    GISMO_ASSERT (rhs.size() == (unsigned)n_levels, "The number of rhs does not correspond to the number of levels!");

    std::vector< gsMatrix<T> > u(n_levels);

    // solve coarse problem
    solveCoarse( rhs[0], u[0] );

    for (index_t i = 1; i <= finestLevel(); ++i)
    {
        // transfer result from previous level
        prolongVector(i-1, u[i-1], u[i]);
        u[i] += dirichletIntp[i];

        // smooth
        smoothingStep(i, rhs[i], u[i]);
    }

    result = u[ finestLevel() ];
}

template<class T>
void gsMultiGridOp<T>::restrictVector(index_t lf, const gsMatrix<T>& fine, gsMatrix<T>& coarse) const
{
    GISMO_ASSERT ( 0 < lf && lf < n_levels, "The given level is not feasible." );
    GISMO_ASSERT ( fine.rows() == nDofs(lf), "The dimensions do not fit." );

    m_restrict[lf-1]->apply( fine, coarse );
}

template<class T>
void gsMultiGridOp<T>::prolongVector(index_t lc, const gsMatrix<T>& coarse, gsMatrix<T>& fine) const
{
    GISMO_ASSERT ( 0 <= lc && lc < n_levels - 1, "The given level is not feasible." );
    GISMO_ASSERT ( coarse.rows() == nDofs(lc), "The dimensions do not fit." );

    m_prolong[lc]->apply( coarse, fine );
}

template<class T>
void gsMultiGridOp<T>::setSmoother(index_t lvl, const PrecondPtr& sm)
{
    GISMO_ASSERT ( 0 <= lvl && lvl < n_levels, "The given level is not feasible." );
    m_smoother[lvl] = sm;
}

template<class T>
const typename gsMultiGridOp<T>::SpMatrix& gsMultiGridOp<T>::matrix(index_t lvl) const
{
    const gsMatrixOp<SpMatrix>* matrOp = dynamic_cast< const gsMatrixOp<SpMatrix>* >( m_ops[lvl].get() );
    GISMO_ASSERT( matrOp, "Matrices are not available for matrix-free multigrid solvers." );
    //return matrOp->matrix(); does not work because we must not return a temporary
    return *(matrOp->matrixPtr());
}

template<class T>
gsOptionList gsMultiGridOp<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addInt   ("NumPreSmooth"                , "Number of pre-smoothing steps",                             1      );
    opt.addInt   ("NumPostSmooth"               , "Number of post-smoothing steps",                            1      );
    opt.addInt   ("NumCycles"                   , "Number of cycles (usually 1 for V-cycle or 2 for W-cycle)", 1      );
    opt.addReal  ("CorarseGridCorrectionDamping", "Damping of the coarse-grid correction (usually 1)",      (T)1      );
    return opt;
}

template<class T>
void gsMultiGridOp<T>::setOptions(const gsOptionList & opt)
{
    Base::setOptions(opt);
    m_numPreSmooth     = opt.askInt   ("NumPreSmooth"                , m_numPreSmooth    );
    m_numPostSmooth    = opt.askInt   ("NumPostSmooth"               , m_numPostSmooth   );
    m_numCycles        = opt.askInt   ("NumCycles"                   , m_numCycles       );
    m_damping          = opt.askReal  ("CorarseGridCorrectionDamping", m_damping         );
}

}
