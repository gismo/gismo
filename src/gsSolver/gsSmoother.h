    /** @file gsSmoother.h

    @brief Provides Multigrid smoothers.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither
*/
 
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsExport.h>
#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

// smoothers


/// @brief A struct containing an enumeration of smoothing algorithms and functions to apply them for use in gsMultiGrid. (OBSOLETE, see gsSmoother below)
/// 
/// \ingroup Solver
struct Smoother
{
    /// List of implemented smoothers
    enum SmootherType {
        Richardson,     ///< damped Richardson
        Jacobi,         ///< damped Jacobi
        GaussSeidel,    ///< Gauss-Seidel
        ILU             ///< incomplete LU
    };

    /// Return the name of the given smoother \a type.
    static std::string name(SmootherType type)
    {
        switch (type)
        {
            case Richardson:    return "Richardson";
            case Jacobi:        return "Jacobi";
            case GaussSeidel:   return "Gauss-Seidel";
            case ILU:           return "ILU";
            default:
                throw std::runtime_error("unknown smoother");
        }
    }
};



/// Abstract base class for multigrid smoothers
class GISMO_EXPORT gsSmoother
{
public:
    virtual ~gsSmoother() { }

    /// Apply the smoother for the equation Ax=f and update the current iterate x.
    virtual void apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f) = 0;

    /// Apply the transposed smoother for the equation Ax=f and update the current iterate x.
    virtual void applyT(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
    {
        apply(A, x, f);     // by default, assume symmetric smoother
    }
};



/// Damped Richardson smoother
class GISMO_EXPORT gsRichardsonSmoother : public gsSmoother
{
public:
    gsRichardsonSmoother(real_t damping = real_t(1))
        : m_damping(damping)
    { }

    void setDamping(real_t d)
    {
        m_damping = d;
    }

    virtual void apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

private:
    real_t m_damping;

};



/// Damped Jacobi smoother
class GISMO_EXPORT gsJacobiSmoother : public gsSmoother
{
public:
    gsJacobiSmoother(real_t damping = real_t(1))
        : m_damping(damping)
    { }

    void setDamping(real_t d)
    {
        m_damping = d;
    }

    virtual void apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

private:
    real_t m_damping;

};


/// Gauss-Seidel smoother
class GISMO_EXPORT gsGaussSeidelSmoother : public gsSmoother
{
public:
    virtual void apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);
    virtual void applyT(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);
};

/// Block Gauss-Seidel smoother
class GISMO_EXPORT gsGaussSeidelBlockSmoother : public gsSmoother
{
public:
    /// Each elemet in \a blockInfo contains a vector with the indices of the DoFs with are grouped
    /// together in one block
    gsGaussSeidelBlockSmoother(std::vector<gsVector<index_t> > &blockInfo)
        : m_blockInfo(blockInfo)
    { }

    virtual void apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);
    virtual void applyT(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

private:
    std::vector<gsVector<index_t> > m_blockInfo;
};


/// ILUT (incomplete LU with thresholding) smoother
class GISMO_EXPORT gsILUTSmoother : public gsSmoother
{
public:
    gsILUTSmoother(const Eigen::SparseMatrix<real_t>& K);

    virtual void apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

private:
    Eigen::IncompleteLUT<real_t> m_ilu;
};



/// Generic smoother which applies an arbitrary linear operator to the residual
class GISMO_EXPORT gsOperatorSmoother : public gsSmoother
{
public:
    /// Takes ownership of the operator
    gsOperatorSmoother(gsLinearOperator * op, real_t damping = 1.0)
        : m_op(op), m_damping(damping)
    {
    }

    virtual ~gsOperatorSmoother()
    {
        delete m_op;
    }

    virtual void apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f)
    {
        assert( A.rows() == x.rows() && x.rows() == f.rows() );
        assert( A.cols() == A.rows() && x.cols() == 1 && f.cols() == 1);

        m_residual.noalias() = m_damping * (f - A * x);
        m_op->apply(m_residual, m_temp);
        x += m_temp;
    }

private:
    gsLinearOperator * m_op;
    real_t m_damping;

    gsMatrix<> m_residual, m_temp;      // keep temporary storage to avoid allocations
};



/// Update \a x with a forward Gauss-Seidel sweep
GISMO_EXPORT void gaussSeidelSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

/// Update \a x with a backward Gauss-Seidel sweep
GISMO_EXPORT void reverseGaussSeidelSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

/// Preforms a block Gauss-Seidel on the degrees of freedom in DoFs.
GISMO_EXPORT void gaussSeidelSingleBlock(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, gsVector<index_t>& DoFs);


} // namespace gismo

