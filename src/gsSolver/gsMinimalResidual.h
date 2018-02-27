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

/** @brief The minimal residual (MinRes) method.
  *
  * \ingroup Solver
  */
template<class T=real_t>
class gsMinimalResidual : public gsIterativeSolver<T>
{

public:
    typedef gsIterativeSolver<T> Base;

    typedef gsMatrix<T>  VectorType;

    typedef typename Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<gsMinimalResidual> Ptr;
    typedef memory::unique_ptr<gsMinimalResidual> uPtr;

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    explicit gsMinimalResidual( const OperatorType& mat,
                                const LinOpPtr& precond = LinOpPtr())
        : Base(mat, precond), m_inexact_residual(false) {}

    /// @brief Make function using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    static uPtr make( const OperatorType& mat, const LinOpPtr& precond = LinOpPtr() )
    { return uPtr( new gsMinimalResidual(mat, precond) ); }

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

    /// @brief Set the options based on a gsOptionList
    gsMinimalResidual& setOptions(const gsOptionList & opt)
    {
        Base::setOptions(opt);
        m_inexact_residual = opt.askSwitch("InexactResidual", m_inexact_residual);
        return *this;
    }

    /// @brief If true, the residual is estimated, not accurately computed.
    void setInexactResidual( bool flag )     { m_inexact_residual = flag; }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "gsMinimalResidual\n";
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

    gsMatrix<T> negResidual,
                     vPrev, v, vNew,
                     wPrev, w, wNew, AwPrev, Aw, AwNew,
                     zNew, z, Az;

    T eta,
           gammaPrev, gamma, gammaNew,
           sPrev, s, sNew,
           cPrev, c, cNew;

    bool m_inexact_residual;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMinimalResidual.hpp)
#endif
