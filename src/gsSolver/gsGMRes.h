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
template<class T = real_t>
class gsGMRes : public gsIterativeSolver<T>
{
public:
    typedef gsIterativeSolver<T> Base;

    typedef gsMatrix<T>  VectorType;

    typedef typename Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<gsGMRes> Ptr;
    typedef memory::unique_ptr<gsGMRes> uPtr;

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    explicit gsGMRes( const OperatorType& mat, const LinOpPtr& precond = LinOpPtr() )
    : Base(mat, precond) {}

    /// @brief Make function using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    static uPtr make( const OperatorType& mat, const LinOpPtr& precond = LinOpPtr() )
    { return uPtr( new gsGMRes(mat, precond) ); }

    bool initIteration( const VectorType& rhs, VectorType& x );
    bool step( VectorType& x );
    void finalizeIteration( VectorType& x );

private:

    /// Solves the Upper triangular system Ry = gg
    /// and stores the solution in the private member y.
    void solveUpperTriangular(const VectorType& R, const VectorType& gg)
    {
       y = R.template triangularView<Eigen::Upper>().solve(gg);
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


    gsMatrix<T> tmp, g, g_tmp, h_tmp, y, w;
    gsMatrix<T> residual;
    gsMatrix<T> H_prev, H, Omega, Omega_prev, Omega_tmp, Omega_prev_tmp;
    std::vector< gsMatrix<T> > v;
    T beta;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsGMRes.hpp)
#endif
