/** @file gsGradientMethod.h

    @brief Gradient iteration

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsIterativeSolver.h>

namespace gismo
{

/// @brief The gradient method.
///
/// Iterates x^{new} = x - damping * precond * ( mat * x - rhs ),
///
/// where by default damping = 1 and precond = identity.
///
/// \ingroup Solver
template<class T = real_t>
class GISMO_EXPORT gsGradientMethod : public gsIterativeSolver<T>
{
public:
    typedef gsIterativeSolver<T> Base;

    typedef gsMatrix<T>  VectorType;

    typedef typename Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<gsGradientMethod> Ptr;
    typedef memory::unique_ptr<gsGradientMethod> uPtr;

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    /// @param damping The damping parameter (step width)
    template< typename OperatorType >
    explicit gsGradientMethod( const OperatorType& mat,
                               const LinOpPtr& precond = LinOpPtr(),
                               T damping = (T)1 )
    : Base(mat, precond), m_damping(damping) {}

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    /// @param damping The damping parameter (step width)
    template< typename OperatorType >
    static uPtr make( const OperatorType& mat,
                      const LinOpPtr& precond = LinOpPtr(),
                      T damping = (T)1 )
    { return uPtr( new gsGradientMethod(mat, precond, damping) ); }

    /// @brief Returns the chosen damping parameter (step width)
    T getDamping() { return m_damping; }

    /// @brief Set the damping parameter (step width)
    void setDamping(T damping) { m_damping = damping; }

    /// @brief Returns a list of default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addReal("Damping", "Damping (step width)", 1 );
        return opt;
    }

    /// @brief Set the options based on a gsOptionList
    gsGradientMethod& setOptions(const gsOptionList & opt)
    {
        Base::setOptions(opt);
        m_damping = opt.askReal("Damping", m_damping);
        return *this;
    }

    bool initIteration( const VectorType& rhs, VectorType& x );
    bool step( VectorType& x );

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "gsGradientMethod\n";
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

    VectorType m_res;
    VectorType m_tmp;
    VectorType m_update;
    T m_damping;

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsGradientMethod.hpp)
#endif
