/** @file gsGMRes.h

    @brief Preconditioned iterative solver using the generalized minimal residual method.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsSolver/gsIterativeSolver.h>

namespace gismo
{

/// @brief The generalized minimal residual (GMRES) method.
///
/// \ingroup Solver    
class GISMO_EXPORT gsGMRes: public gsIterativeSolver<real_t>
{
public:
    typedef gsIterativeSolver<real_t> Base;
    
    typedef gsMatrix<real_t>  VectorType;
    
    typedef Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<gsGMRes> Ptr;
    typedef memory::unique_ptr<gsGMRes> uPtr;

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    explicit gsGMRes( const OperatorType& mat, const LinOpPtr & precond = LinOpPtr() )
    : Base(mat, precond) {}

    bool initIteration( const VectorType& rhs, VectorType& x );
    bool step( VectorType& x );
    void finalizeIteration( VectorType& x );

private:

    /// Solves the Upper triangular system Ry = gg
    /// and stores the solution in the private member y.
    void solveUpperTriangular(const gsMatrix<real_t> & R, const gsMatrix<real_t> & gg)
    {
       y = R.triangularView<Eigen::Upper>().solve(gg);
    }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "gsGMRes\n";
        return os;
    }

private:
    using Base::m_mat;
    using Base::m_precond;
    using Base::m_max_iters;
    using Base::m_tol;
    using Base::m_num_iter;
    using Base::m_rhs_norm;
    using Base::m_error;


    gsMatrix<real_t> tmp, g, g_tmp, h_tmp, y, w;
    gsMatrix<real_t> residual;
    gsMatrix<real_t> H_prev, H, Omega, Omega_prev, Omega_tmp, Omega_prev_tmp;
    std::vector<gsMatrix<real_t> > v;
    real_t beta;
};

} // namespace gismo
