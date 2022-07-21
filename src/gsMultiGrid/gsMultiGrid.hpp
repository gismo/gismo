/** @file gsMultiGrid.hpp

    @brief Multigrid preconditioner

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#include <gsSolver/gsMatrixOp.h>

namespace gismo
{

template<class T>
gsMultiGridOp<T>::gsMultiGridOp(
    SpMatrix fineMatrix,
    std::vector<SpMatrixRowMajor> transferMatrices,
    OpPtr coarseSolver
)
    : m_nLevels(transferMatrices.size()+1), m_ops(m_nLevels), m_smoother(m_nLevels),
      m_prolong(m_nLevels-1), m_restrict(m_nLevels-1), m_coarseSolver(coarseSolver),
      m_numPreSmooth(1), m_numPostSmooth(1), m_symmSmooth(true), m_damping(1)
{
    m_numCycles.setConstant(m_nLevels-1,1);
    std::vector<SpMatrixRowMajorPtr> transferMatrixPtrs(m_nLevels-1);
    for (index_t i=0; i<m_nLevels-1; ++i)
        transferMatrixPtrs[i] = transferMatrices[i].moveToPtr();

    init(fineMatrix.moveToPtr(),give(transferMatrixPtrs));
}

template<class T>
gsMultiGridOp<T>::gsMultiGridOp(
    SpMatrixPtr fineMatrix,
    std::vector<SpMatrixRowMajorPtr> transferMatrices,
    OpPtr coarseSolver
)
    : m_nLevels(transferMatrices.size()+1), m_ops(m_nLevels), m_smoother(m_nLevels),
      m_prolong(m_nLevels-1), m_restrict(m_nLevels-1), m_coarseSolver(coarseSolver),
      m_numPreSmooth(1), m_numPostSmooth(1), m_symmSmooth(true), m_damping(1)
{
    m_numCycles.setConstant(m_nLevels-1,1);
    init(give(fineMatrix),give(transferMatrices));
}

template<class T>
gsMultiGridOp<T>::gsMultiGridOp(
    const std::vector<OpPtr>& ops,
    const std::vector<OpPtr>& prolongation,
    const std::vector<OpPtr>& restriction,
    OpPtr coarseSolver
)
    : m_nLevels(ops.size()), m_ops(ops), m_smoother(m_nLevels),
      m_prolong(prolongation), m_restrict(restriction), m_coarseSolver(coarseSolver),
      m_numPreSmooth(1), m_numPostSmooth(1), m_symmSmooth(true), m_damping(1)
{
    m_numCycles.setConstant(m_nLevels-1,1);
    GISMO_ASSERT (prolongation.size() == restriction.size(),
        "gsMultiGridOp: The number of prolongation and restriction operators differ.");
    GISMO_ASSERT (ops.size() == prolongation.size()+1,
        "gsMultiGridOp: The number of prolongation and restriction operators do not fit to the number of operators.");
}

template<class T>
void gsMultiGridOp<T>::init(SpMatrixPtr fineMatrix, std::vector<SpMatrixRowMajorPtr> transferMatrices)
{
    GISMO_ASSERT (fineMatrix->rows() == fineMatrix->cols(), "gsMultiGridOp needs quadratic matrices.");

    SpMatrixPtr mat = fineMatrix;
    m_ops[m_nLevels-1] = makeMatrixOp(mat);

    for ( index_t i = m_nLevels - 2; i >= 0; --i )
    {
        SpMatrixPtr newMat = SpMatrixPtr(new SpMatrix(
            transferMatrices[i]->transpose() * *mat * *(transferMatrices[i])
        ));
        m_ops[i] = makeMatrixOp(newMat);
        mat = newMat; // copies just the smart pointers
    }

    for ( index_t i=0; i<m_nLevels-1; ++i )
    {
        m_prolong[i] = makeMatrixOp(transferMatrices[i]);
        // Note that the following works since we know that m_prolong does not
        // get destroyed beforem_restrict, so the shared pointer stored in
        // m_prolong will make sure that the underlying matrix will not get
        // destroyed. Thus, there will be no dangling pointer.
        m_restrict[i] = makeMatrixOp(transferMatrices[i]->transpose());
    }
}

// This function is const since it is called from "apply", which itself is const.
// Since this function realizes a late initialization for the coarse solver, it is
// semantically const. We want late initialization since this gives the caller a
// better chance to use an alternate solver.
template<class T>
void gsMultiGridOp<T>::initCoarseSolver() const
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
            gsWarn << "gsMultiGridOp::initCoarseSolverIfUninitialized(): The coarse grid solver is constructed "
                "based on gsLinearOperator::toMatrix(). This might be inefficient. Consider providing matrices of "
                "type gsSparseMatrix<T> or an exact solver for the coarset grid level to gsMultiGridOp constructor.\n";
        Matrix coarse_dense;
        m_ops[0]->toMatrix( coarse_dense );
        m_coarseSolver = makePartialPivLUSolver( coarse_dense );
    }
}

template<class T>
void gsMultiGridOp<T>::smoothingStep(index_t level, const Matrix& rhs, Matrix& x) const
{
    GISMO_ASSERT (m_smoother[level], "gsMultiGridOp::smoothingStep: "
        "Smoother is not defined for level "<<level<<". Define it using setSmoother.");

    // pre-smooth
    for (index_t i = 0; i < m_numPreSmooth; ++i)
    {
        m_smoother[level]->step( rhs, x );
    }

    // post-smooth
    for (index_t i = 0; i < m_numPostSmooth; ++i)
    {
        if (m_symmSmooth)
            m_smoother[level]->stepT( rhs, x );
        else
            m_smoother[level]->step( rhs, x );
    }

}

template<class T>
void gsMultiGridOp<T>::multiGridStep(index_t level, const Matrix& rhs, Matrix& x) const
{
    GISMO_ASSERT ( 0 <= level && level < m_nLevels, "TgsMultiGridOp: the given level is not feasible." );

    if (level == 0)
    {
        solveCoarse(rhs, x);
    }
    else
    {
        const index_t lf = level;
        const index_t lc = lf - 1;

        GISMO_ASSERT (m_smoother[lf], "gsMultiGridOp::multiGridStep: "
            "Smoother is not defined for level "<<lf<<". Define it using setSmoother.");

        Matrix fineRes, fineCorr, coarseRes, coarseCorr;

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
        for (index_t i = 0; i < m_numCycles[lc]; ++i)
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
            if (m_symmSmooth)
                m_smoother[lf]->stepT( rhs, x );
            else
                m_smoother[lf]->step( rhs, x );
        }
    }
}

template<class T>
void gsMultiGridOp<T>::fullMultiGrid(
    const std::vector<Matrix>& rhs,
    const std::vector<Matrix>& fixedValues,
    gsMatrix<T>& result
) const
{

    GISMO_ASSERT (fixedValues.size() == static_cast<size_t>(m_nLevels),
        "gsMultiGridOp::fullMultiGrid: The size of fixedValues does not correspond to the number of levels.");
    GISMO_ASSERT (rhs.size() == static_cast<size_t>(m_nLevels),
        "gsMultiGridOp::fullMultiGrid: The size of rhs does not correspond to the number of levels.");

    std::vector<Matrix> u(m_nLevels);

    // solve coarse problem
    solveCoarse( rhs[0], u[0] );

    for (index_t i = 1; i <= finestLevel(); ++i)
    {
        // transfer result from previous level
        prolongVector(i-1, u[i-1], u[i]);
        u[i] += fixedValues[i];

        // run one multigrid step
        multiGridStep(i, rhs[i], u[i]);
    }

    result = u[ finestLevel() ];
}

template<class T>
void gsMultiGridOp<T>::cascadicMultiGrid(
    const std::vector<Matrix>& rhs,
    const std::vector<Matrix>& fixedValues,
    Matrix& result
) const
{

    GISMO_ASSERT (fixedValues.size() == static_cast<size_t>(m_nLevels),
        "gsMultiGridOp::cascadicMultiGrid: The size of fixedValue does not correspond to the number of levels.");
    GISMO_ASSERT (rhs.size() == static_cast<size_t>(m_nLevels),
        "gsMultiGridOp::cascadicMultiGrid: The size of rhs does not correspond to the number of levels.");

    std::vector<Matrix> u(m_nLevels);

    // solve coarse problem
    solveCoarse( rhs[0], u[0] );

    for (index_t i = 1; i <= finestLevel(); ++i)
    {
        // transfer result from previous level
        prolongVector(i-1, u[i-1], u[i]);
        u[i] += fixedValues[i];

        // smooth
        smoothingStep(i, rhs[i], u[i]);
    }

    result = u[ finestLevel() ];
}

template<class T>
void gsMultiGridOp<T>::restrictVector(index_t lf, const Matrix& fine, Matrix& coarse) const
{
    GISMO_ASSERT ( 0 < lf && lf < m_nLevels, "gsMultiGrid: The given level is not feasible." );
    GISMO_ASSERT ( fine.rows() == nDofs(lf), "gsMultiGrid: The dimensions do not fit." );

    m_restrict[lf-1]->apply( fine, coarse );
}

template<class T>
void gsMultiGridOp<T>::prolongVector(index_t lc, const Matrix& coarse, Matrix& fine) const
{
    GISMO_ASSERT ( 0 <= lc && lc < m_nLevels - 1, "gsMultiGrid: The given level is not feasible." );
    GISMO_ASSERT ( coarse.rows() == nDofs(lc), "gsMultiGrid: The dimensions do not fit." );

    m_prolong[lc]->apply( coarse, fine );
}

template<class T>
void gsMultiGridOp<T>::setSmoother(index_t lvl, const PrecondPtr& sm)
{
    GISMO_ASSERT ( 0 <= lvl && lvl < m_nLevels, "gsMultiGrid: The given level is not feasible." );
    m_smoother[lvl] = sm;
}

template<class T>
const typename gsMultiGridOp<T>::SpMatrix& gsMultiGridOp<T>::matrix(index_t lvl) const
{
    const gsMatrixOp<SpMatrix>* matrOp = dynamic_cast< const gsMatrixOp<SpMatrix>* >( m_ops[lvl].get() );
    GISMO_ASSERT ( matrOp, "Matrices are not available for matrix-free multigrid solvers." );
    //return matrOp->matrix(); does not work because we must not return a temporary
    return *(matrOp->matrixPtr());
}

template<class T>
gsOptionList gsMultiGridOp<T>::defaultOptions()
{
    gsOptionList opt = Base::defaultOptions();
    opt.addInt   ("NumPreSmooth"                , "Number of pre-smoothing steps",                                         1 );
    opt.addInt   ("NumPostSmooth"               , "Number of post-smoothing steps",                                        1 );
    opt.addInt   ("NumCycles"                   , "Number of cycles (usually 1 for V-cycle or 2 for W-cycle)",             1 );
    opt.addReal  ("CoarseGridCorrectionDamping" , "Damping of the coarse-grid correction (usually 1)", (gsOptionList::Real)1 );
    opt.addSwitch("SymmSmooth"                  , "Iff true, stepT is called for post-smoothing",           true );

    return opt;
}

template<class T>
void gsMultiGridOp<T>::setOptions(const gsOptionList & opt)
{
    Base::setOptions(opt);
    m_numPreSmooth     = opt.askInt   ("NumPreSmooth"                , m_numPreSmooth    );
    m_numPostSmooth    = opt.askInt   ("NumPostSmooth"               , m_numPostSmooth   );
    m_damping          = opt.askReal  ("CoarseGridCorrectionDamping" , m_damping         );
    m_symmSmooth       = opt.askSwitch("SymmSmooth"                  , m_symmSmooth      );

    const index_t nc   = opt.askInt   ("NumCycles"                   , -1                );
    if (nc > -1)
    {
        m_numCycles.setConstant(m_nLevels-1, nc);
        // The direct solver on coarsest level is only invoked once
        m_numCycles[0] = 1;
    }
}

}
