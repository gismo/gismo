/** @file gsConjugateGradient.h

    @brief Conjugate gradient solver

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/

#pragma once

#include <gsSolver/gsIterativeSolver.h>

namespace gismo
{

/// @brief The conjugate gradient method.
///
/// The conjugate gradient implementation from Eigen, adapted to allow for more
/// general preconditioners and better iteration control. Also capable of using
/// a gsLinearOperator as matrix.
///
/// \ingroup Solver
template<class T = real_t>
class gsConjugateGradient : public gsIterativeSolver<T>
{
public:
    typedef gsIterativeSolver<T> Base;

    typedef gsMatrix<T>  VectorType;

    typedef typename Base::LinOpPtr LinOpPtr;

    typedef memory::shared_ptr<gsConjugateGradient> Ptr;
    typedef memory::unique_ptr<gsConjugateGradient> uPtr;

    /// @brief Constructor using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    explicit gsConjugateGradient( const OperatorType& mat,
                                  const LinOpPtr& precond = LinOpPtr() )
    : Base(mat, precond), m_calcEigenvals(false) {}

    /// @brief Make function using a matrix (operator) and optionally a preconditionner
    ///
    /// @param mat     The operator to be solved for, see gsIterativeSolver for details
    /// @param precond The preconditioner, defaulted to the identity
    template< typename OperatorType >
    static uPtr make( const OperatorType& mat, const LinOpPtr& precond = LinOpPtr() )
    { return uPtr( new gsConjugateGradient(mat, precond) ); }

    /// @brief Returns a list of default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addSwitch("CalcEigenvalues", "Additionally to solving the system,"
                      " CG computes the eigenvalues of the Lanczos matrix", false );
        return opt;
    }

    /// @brief Set the options based on a gsOptionList
    gsConjugateGradient& setOptions(const gsOptionList& opt)
    {
        Base::setOptions(opt);
        m_calcEigenvals = opt.askSwitch("CalcEigenvalues", m_calcEigenvals);
        return *this;
    }

    bool initIteration( const VectorType& rhs, VectorType& x );
    bool step( VectorType& x );

    /// @brief specify if you want to store data for eigenvalue estimation
    /// @param flag true stores the coefficients of the lancos matrix, false not.
    void setCalcEigenvalues( bool flag )     { m_calcEigenvals = flag ;}

    /// @brief returns the condition number of the (preconditioned) system matrix
    T getConditionNumber();

    /// @brief returns the eigenvalues of the Lanczos matrix
    void getEigenvalues( VectorType& eigs );

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "gsConjugateGradient\n";
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
    VectorType m_update;
    VectorType m_tmp;
    T m_abs_new;

    bool m_calcEigenvals;

    std::vector<T> m_delta, m_gamma;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsConjugateGradient.hpp)
#endif
