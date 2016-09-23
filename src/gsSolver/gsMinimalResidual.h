/** @file gsMinimalResidual.h

    @brief Preconditioned iterative solver using the minimal residual method.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, S. Takacs
*/
#pragma once

#include <gsSolver/gsIterativeSolver.h>

namespace gismo
{

class GISMO_EXPORT gsMinimalResidual : public gsIterativeSolver
{

public:
    typedef gsIterativeSolver Base;
    
    typedef gsMatrix<real_t>  VectorType;

    typedef Base::LinOpPtr LinOpPtr;
        
    /// Constructor using a matrix (operator) and optionally a preconditionner
    template< typename OperatorType >
    explicit gsMinimalResidual( const OperatorType& mat,
                                const LinOpPtr& precond = LinOpPtr())
        : Base(mat, precond) { }
    
    bool initIteration( const VectorType& rhs, VectorType& x );

    bool step( VectorType& x );

private:
    using Base::m_mat;
    using Base::m_precond;
    using Base::m_max_iters;
    using Base::m_tol;
    using Base::m_num_iter;
    using Base::m_rhs_norm;
    using Base::m_error;

    gsMatrix<real_t> negResidual,
                     vPrev, v, vNew,
                     wPrev, w, wNew, AwPrev, Aw, AwNew,
                     zNew, z, Az;

    real_t eta,
           gammaPrev, gamma, gammaNew,
           sPrev, s, sNew,
           cPrev, c, cNew;
};

} // namespace gismo

