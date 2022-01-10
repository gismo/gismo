/** @file gsBiCgStab.h

    @brief Biconjugate gradient stabilized solver

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): R. Tielen, S. Takacs
*/

#pragma once

#include <gsSolver/gsIterativeSolver.h>

namespace gismo
{

/// @brief Biconjugate gradient stabilized solver
///
/// \ingroup Solver
template<class T = real_t>
class gsBiCgStab : public gsIterativeSolver<T>
{
public:
    typedef gsIterativeSolver<T> Base;

    typedef gsMatrix<T>  VectorType;

    typedef typename Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<gsBiCgStab> Ptr;
    typedef memory::unique_ptr<gsBiCgStab> uPtr;

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    explicit gsBiCgStab( const OperatorType& mat, const LinOpPtr& precond = LinOpPtr() )
    : Base(mat, precond), m_restartThereshold(1e-32) {}

    /// @brief Make function using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    static uPtr make( const OperatorType& mat, const LinOpPtr& precond = LinOpPtr() )
    { return uPtr( new gsBiCgStab(mat, precond) ); }

    /// @brief Returns a list of default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addReal("RestartThereshold", "Thereshold for restarting gsBiCgStab solver.", 1e-32 );
        return opt;
    }

    /// @brief Set the options based on a gsOptionList
    gsBiCgStab& setOptions(const gsOptionList& opt)
    {
        Base::setOptions(opt);
        m_restartThereshold = opt.askReal("RestartThereshold", m_restartThereshold);
        return *this;
    }

    bool initIteration( const VectorType& rhs, VectorType& x );
    bool step( VectorType& x );

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "gsBiCgStab\n";
        return os;
    }

private:
    using Base::m_mat;
    using Base::m_precond;
    using Base::m_max_iters;
    using Base::m_tol;
    using Base::m_num_iter;
    using Base::m_initial_error;
    using Base::m_current_error;

    VectorType m_res;
    VectorType m_r0;
    VectorType m_tmp;
    VectorType m_v;
    VectorType m_p;
    VectorType m_y;
    VectorType m_s;
    VectorType m_z;
    VectorType m_t;

    T m_alpha;
    T m_rho;
    T m_w;

    T m_restartThereshold;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBiCgStab.hpp)
#endif
