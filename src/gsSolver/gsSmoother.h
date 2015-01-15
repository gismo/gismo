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

namespace gismo
{

// smoothers


/// A struct containing an enumeration of smoothing algorithms and functions to apply them for use in gsMultiGrid. (OBSOLETE, see gsSmoother below)
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


class gsSmoother
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


class GISMO_EXPORT gsDampedPrecRichardsonSmoother : public gsSmoother
{
public:
    gsDampedPrecRichardsonSmoother(const Eigen::SparseMatrix<real_t>* P, real_t damping = real_t(1))
        : m_prec(P), m_damping(damping)
    {
        m_solver.compute(*m_prec);
        //m_solver.setTolerance(1e-9);
    }

    void setDamping(real_t d)
    {
        m_damping = d;
    }

    virtual void apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

private:
    const Eigen::SparseMatrix<real_t>* m_prec;
    real_t m_damping;
    Eigen::SparseLU< Eigen::SparseMatrix<real_t> > m_solver;

};



class GISMO_EXPORT gsGaussSeidelSmoother : public gsSmoother
{
public:
    virtual void apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);
    virtual void applyT(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);
};


class GISMO_EXPORT gsILUTSmoother : public gsSmoother
{
public:
    gsILUTSmoother(const Eigen::SparseMatrix<real_t>& K);

    virtual void apply(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

private:
    Eigen::IncompleteLUT<real_t> m_ilu;
};


GISMO_EXPORT void reverseGaussSeidelSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);
GISMO_EXPORT void gaussSeidelSweep(const Eigen::SparseMatrix<real_t>& A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsSmoother.cpp)
#endif
