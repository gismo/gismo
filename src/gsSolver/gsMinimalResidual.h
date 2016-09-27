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

class GISMO_EXPORT gsMinimalResidual : public gsIterativeSolver<real_t>
{

public:
    typedef gsIterativeSolver<real_t> Base;
    
    typedef gsMatrix<real_t>  VectorType;

    typedef Base::LinOpPtr LinOpPtr;
        
    /// Constructor using a matrix (operator) and optionally a preconditionner
    template< typename OperatorType >
    explicit gsMinimalResidual( const OperatorType& mat,
                                const LinOpPtr& precond = LinOpPtr())
        : Base(mat, precond), m_inexact_residual(false) { }
    
    bool initIteration( const VectorType& rhs, VectorType& x );
    void finalizeIteration( VectorType& x );
    
    bool step( VectorType& x );

    /// @brief Returns a list of default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addSwitch( "InexactResidual",
                       "If true, the residual is estimated, not accurately computed.",
                       false );
        return opt;
    }

    void setOptions(const gsOptionList & opt)
    {
        Base::setOptions(opt);
        m_inexact_residual = opt.askSwitch("InexactResidual", m_inexact_residual);
    }

    void setSloppyResidual( bool flag )     { m_inexact_residual = flag ;}

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

    bool m_inexact_residual;
};

} // namespace gismo

