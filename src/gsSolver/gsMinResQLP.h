/** @file gsMinResQLP.h

    @brief Preconditioned iterative solver using the minimal residual QLP method.

    Implementation from the algorithm
    MINRES-QLP: A KRYLOV SUBSPACE METHOD FOR INDEFINITE OR SINGULAR SYMMETRIC SYSTEMS
    BY: SOU-CHENG T. CHOI , CHRISTOPHER C. PAIGE , AND MICHAEL A. SAUNDERS

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#pragma once

#include <gsSolver/gsIterativeSolver.h>

//TODOs: Create accessors for values of interest, like, Arnorm, xnorm etc.
//TODOs: Add the same flaging system the matlab code uses
//TODOs: Merge the minres and minresQLP code like petsc and matlab has
namespace gismo
{

/** @brief The minimal residual (MinRes-QLP) method.
  *
  * \ingroup Solver
  */
template<class T=real_t>
class gsMinResQLP : public gsIterativeSolver<T>
{

public:
    typedef gsIterativeSolver<T> Base;

    typedef gsMatrix<T>  VectorType;

    typedef typename Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<gsMinResQLP> Ptr;
    typedef memory::unique_ptr<gsMinResQLP> uPtr;

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    explicit gsMinResQLP( const OperatorType& mat,
                                const LinOpPtr& precond = LinOpPtr())
        : Base(mat, precond), m_inexact_residual(false) {}

    /// @brief Make function using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    static uPtr make( const OperatorType& mat, const LinOpPtr& precond = LinOpPtr() )
    { return uPtr( new gsMinResQLP(mat, precond) ); }

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
    gsMinResQLP& setOptions(const gsOptionList & opt)
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
        os << "gsMinResQLP\n";
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

    gsMatrix<T> m_rhs, r1, r2, r3, resvec, Aresvec,
                v, xl2, negResidual,
                wl, w, wl2;


    T beta1, betal, beta, tau, taul, taul2, phi, betan, gmin, cs, sn,cr1,sr1,
        cr2,sr2, deltan,gamma,gammal, gammal2, gammal3, eta,etal, etal2,
        vepln, veplnl, veplnl2, u, ul, ul2, ul3, rnorm, rnorml, xnorm, xl2norm,
        Axnorm, Anorm, Acond, relres, alpha, pnorm, dbar, delta,
        epsilon, epsilonn, gbar, maxxnorm, gammal_QLP, gamma_QLP, vepln_QLP,
        abs_gamma, gminl, gminl2, Arnorml, Arnorm, relAresl, relAres,
        rootl, ul_QLP, u_QLP;;


    index_t QLPiter;
    bool m_inexact_residual;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsMinResQLP.hpp)
#endif
