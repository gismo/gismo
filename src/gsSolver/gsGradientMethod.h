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
class gsGradientMethod : public gsIterativeSolver<T>
{
public:
    typedef gsIterativeSolver<T> Base;

    typedef gsMatrix<T>  VectorType;

    typedef typename Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<gsGradientMethod> Ptr;
    typedef memory::unique_ptr<gsGradientMethod> uPtr;

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat       The operator to be solved for, see gsIterativeSolver for details
    /// @param precond   The preconditioner, defaulted to the identity
    ///
    /// In each step, the step width is chosen such that the residual is minimized (adaptive
    /// step sizes)
    template< typename OperatorType >
    explicit gsGradientMethod( const OperatorType& mat,
                               const LinOpPtr& precond = LinOpPtr())
    : Base(mat, precond), m_adapt_step_size(true), m_step_size(0) {}

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat       The operator to be solved for, see gsIterativeSolver for details
    /// @param precond   The preconditioner, defaulted to the identity
    /// @param step_size The step size
    template< typename OperatorType >
    explicit gsGradientMethod( const OperatorType& mat,
                               const LinOpPtr& precond,
                               T step_size )
    : Base(mat, precond), m_adapt_step_size(false), m_step_size(step_size) {}

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat       The operator to be solved for, see gsIterativeSolver for details
    /// @param precond   The preconditioner, defaulted to the identity
    ///
    /// In each step, the step width is chosen such that the residual is minimized (adaptive
    /// step sizes)
    template< typename OperatorType >
    static uPtr make( const OperatorType& mat,
                      const LinOpPtr& precond = LinOpPtr())
    { return uPtr( new gsGradientMethod(mat, precond) ); }

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat       The operator to be solved for, see gsIterativeSolver for details
    /// @param precond   The preconditioner, defaulted to the identity
    /// @param step_size The step size
    template< typename OperatorType >
    static uPtr make( const OperatorType& mat,
                      const LinOpPtr& precond,
                      T step_size )
    { return uPtr( new gsGradientMethod(mat, precond, step_size) ); }

    /// @brief Returns true iff adaptive step sizes are activated
    bool adaptiveStepSize() const { return m_adapt_step_size;                           }

    /// @brief Returns the chosen step size
    T stepSize() const            { return m_step_size;                                 }

    /// @brief Activate adaptive step size. Then in each step, the step size
    /// is chosen such that the norm of the residual is minimized.
    void setAdaptiveStepSize()    { m_adapt_step_size = true;  m_step_size = 0;         }

    /// @brief Set the step size
    void setStepSize(T step_size) { m_adapt_step_size = false; m_step_size = step_size; }

    /// @brief Returns a list of default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addSwitch("AdaptiveStepSize", "Adaptive step sizes (to minimize residual)", true );
        opt.addReal  ("StepSize", "Step size (requires AdaptiveStepSize to be false)" , (T)0    );
        return opt;
    }

    /// @brief Set the options based on a gsOptionList
    gsGradientMethod& setOptions(const gsOptionList & opt)
    {
        Base::setOptions(opt);
        m_adapt_step_size = opt.askSwitch("AdaptiveStepSize", m_adapt_step_size);
        m_step_size       = opt.askReal  ("StepSize",         m_step_size      );
        if (m_adapt_step_size) m_step_size = 0; // for consistency of output, see above
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
    bool m_adapt_step_size;
    T m_step_size;

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsGradientMethod.hpp)
#endif
