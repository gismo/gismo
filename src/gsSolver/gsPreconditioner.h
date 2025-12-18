/** @file gsPreconditioner.h

    @brief Simple abstract class for (discrete) linear preconditioners.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, A. Manzaflaris, C. Hofreither, S. Takacs, C. Hofer
*/
#pragma once

#include <gsIO/gsOptionList.h>
#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/// @brief Simple abstract class for perconditioners.
///
/// The class represents an iteration method in the form
///
/// \f$ x_{new} = x_{old} + P(f - A*x_{old}).\f$
///
/// If the number of steps is set to 1, the member function \a apply represents
/// the application of \f$ P \f$. If \a setNumOfSweeps is used to set the
/// number of sweeps to some \f$ \nu>0 \f$, the member function \a apply
/// realizes \f$ I-(I-PA)^{\nu} A^{-1} \f$. The same interpretation applies
/// if the object is used as a gsLinearOperator.
///
/// The member function \a step represents one step of the above iteration method.
///
/// Usually, the step operation can be performed in an optimized way.
///
/// The derived classes have to contain the functions: \a step, \a stepT,
/// \a rows, \a cols.
///
/// \ingroup Solver
template<class T>
class gsPreconditionerOp : public gsLinearOperator<T>
{
public:

    /// Shared pointer for gsLinearOperator
    typedef memory::shared_ptr<gsPreconditionerOp> Ptr;

    /// Unique pointer for gsLinearOperator
    typedef memory::unique_ptr<gsPreconditionerOp> uPtr;

    /// Base class
    typedef gsLinearOperator<T> Base;

    /// Base class
    typedef typename gsLinearOperator<T>::Ptr BasePtr;

    gsPreconditionerOp() : m_num_of_sweeps(1) {}

    /**
     * @brief Apply the method for given right hand side and current iterate
     * @param rhs Right hand side vector
     * @param x   Current iterate vector
     */
    virtual void step(const gsMatrix<T> & rhs, gsMatrix<T> & x) const = 0;

    /**
     * @brief Apply the transposed variant of the method for given right hand
     *        side and current iterate
     * @param rhs Right hand side vector
     * @param x   Current iterate vector
     *
     * @warning Derived classes *must* overwrite this virtual function if the
     * preconditioner is not symmetric.
     */
    virtual void stepT(const gsMatrix<T> & rhs, gsMatrix<T> & x) const
    { step(rhs, x); } // Assume symmetry.

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        x.setZero(this->rows(),input.cols()); // we assume quadratic matrices
        for ( index_t i = 0; i < m_num_of_sweeps; ++i )
            step(input,x);
    }

    /// Return the underlying operator \f$ A \f$.
    virtual BasePtr underlyingOp() const = 0;

    /// Set the number of sweeps to be applied in the member function apply
    void setNumOfSweeps( index_t n )
    {
        GISMO_ASSERT ( n > 0, "Number of sweeps needs to be positive." );
        m_num_of_sweeps = n;
    }

    /// Get the number of sweeps to be applied in the member function \a apply
    index_t numOfSweeps()
    {
        return m_num_of_sweeps;
    }

    /// Get the default options as gsOptionList object
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addInt( "NumOfSweeps", "Number of sweeps to be applied in the member function apply", 1 );
        return opt;
    }

    /// Set options based on a gsOptionList object
    virtual void setOptions(const gsOptionList & opt)
    {
        Base::setOptions(opt);
        m_num_of_sweeps = opt.askInt( "NumOfSweeps", m_num_of_sweeps );
    }

    /**
     *  @brief Estimates the largest eigenvalue of  \f$ PA \f$
     *  @param steps  Number of steps to be performed.
     */
    T estimateLargestEigenvalueOfPreconditionedSystem(index_t steps = 10) const
    {
        gsMatrix<T> x, tmp;
        x.setRandom(this->rows(),1);
        for (index_t i=0; i<steps; ++i)
        {
            x.array() /= math::sqrt( x.col(0).dot(x.col(0)) );
            underlyingOp()->apply(x,tmp);
            x.setZero(this->rows(),1);
            this->step(tmp,x);
        }
        return math::sqrt( x.col(0).dot(x.col(0)) );
    }

protected:
    index_t m_num_of_sweeps;

}; // gsPreconditionerOp

/// @brief Simple class allowing to construct a preconditioner from a
/// linear operator.
///
/// The class represents an iteration method in the form
///
/// \f$ x_{new} = x_{old} + \tau P (f - A*x_{old}).\f$
///
/// @warning The implemenation pf stepT assumes that P is symmetric.
///
/// \ingroup Solver
template<class T>
class gsPreconditionerFromOp GISMO_FINAL : public gsPreconditionerOp<T>
{
public:

    /// Shared pointer for gsLinearOperator
    typedef memory::shared_ptr<gsPreconditionerFromOp> Ptr;

    /// Unique pointer for gsLinearOperator
    typedef memory::unique_ptr<gsPreconditionerFromOp> uPtr;

    /// Base class
    typedef gsLinearOperator<T> Base;

    /// Base class
    typedef typename gsLinearOperator<T>::Ptr BasePtr;

    /**
     * @brief Constructor
     * @param underlying      The underlying operator \f$ A \f$.
     * @param preconditioner  The operator \f$ P \f$ to be used as preconditioner.
     * @param tau             A damping parameter, defaulted to 1.
     */
    gsPreconditionerFromOp( BasePtr underlying, BasePtr preconditioner, T tau = (T)1)
        : m_underlying(give(underlying)), m_preconditioner(give(preconditioner)), m_tau(tau)
    {
        GISMO_ASSERT( m_underlying->rows() == m_underlying->cols()
            && m_preconditioner->rows() == m_preconditioner->cols()
            && m_underlying->rows() == m_preconditioner->rows(),
            "The dimensions do not agree." );
    }

    /**
     * @brief Make function returning a smart pointer
     * @param underlying      The underlying operator \f$ A \f$.
     * @param preconditioner  The underlying preconditioner \f$ P \f$.
     * @param tau             A damping parameter, defaulted to 1.
     */
    static uPtr make( BasePtr underlying, BasePtr preconditioner, T tau = (T)1)
    { return uPtr( new gsPreconditionerFromOp(give(underlying), give(preconditioner), tau) ); }

    void step(const gsMatrix<T> & rhs, gsMatrix<T> & x) const
    {
        GISMO_ASSERT( m_underlying->rows() == x.rows() && x.rows() == rhs.rows() && x.cols() == rhs.cols(),
            "The dimensions do not agree." );

        m_underlying->apply(x, m_res);
        m_res -= rhs;
        m_preconditioner->apply(m_res, m_corr);
        x -= m_tau * m_corr;
    }

    // virtual void stepT(const gsMatrix<T> & rhs, gsMatrix<T> & x) const { step( rhs, x); } // Assume symmetry.

    void apply(const gsMatrix<T> & rhs, gsMatrix<T> & x) const
    {
        GISMO_ASSERT( m_underlying->rows() == rhs.rows(), "The dimensions do not agree." );

        // special treatment for the first step
        m_preconditioner->apply(rhs, x);
        if (m_tau != (T)1)
            x *= m_tau;

        for (index_t i=1; i<m_num_of_sweeps; ++i)
        {
            m_underlying->apply(x, m_res);
            m_res -= rhs;
            m_preconditioner->apply(m_res, m_corr);
            x -= m_tau * m_corr;
        }
    }

    /// Set damping parameter
    void setDamping(const T tau) { m_tau = tau;  }

    /// Get damping parameter
    T getDamping() const         { return m_tau; }

    /// Get the default options as gsOptionList object
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addReal( "Damping", "Damping parameter of the operator preconditioner", 1 );
        return opt;
    }

    /// Set options based on a gsOptionList object
    void setOptions(const gsOptionList & opt)
    {
        Base::setOptions(opt);
        m_tau = opt.askReal( "Damping", m_tau );
    }

    BasePtr underlyingOp() const { return m_underlying; }
    index_t rows() const { return m_preconditioner->rows(); }
    index_t cols() const { return m_preconditioner->cols(); }

protected:
    BasePtr m_underlying;
    BasePtr m_preconditioner;
    T m_tau;
    mutable gsMatrix<T> m_res;
    mutable gsMatrix<T> m_corr;
    using gsPreconditionerOp<T>::m_num_of_sweeps;
}; // gsPreconditionerFromOp

/// @brief Wrapper for classes inheriting from gsPreconditionerOp
/// This wrapper adapts any Gismo preconditioner operator for use with Eigen's iterative solvers
template <class T>
class gsPreconditionerWrapper
{
    typedef T Scalar;
    typedef gsEigen::Matrix<Scalar,gsEigen::Dynamic,1> Vector;
    typedef typename Vector::StorageIndex StorageIndex;

public:
    typedef gsEigen::SparseMatrix<Scalar,0,index_t> EigenSparseMatrixType;

    /// Default constructor
    gsPreconditionerWrapper()
    :
    m_preconditioner(nullptr)
    {}

    /// Constructor with preconditioner
    explicit gsPreconditionerWrapper(typename gsPreconditionerOp<Scalar>::uPtr preconditioner)
    :
    m_preconditioner(give(preconditioner))
    {}

    /// Move constructor
    gsPreconditionerWrapper(gsPreconditionerWrapper&& other) noexcept
    :
    m_preconditioner(give(other.m_preconditioner))/* ,
    m_stored_matrix(give(other.m_stored_matrix)) */
    {}

    /// Move assignment operator
    gsPreconditionerWrapper& operator=(gsPreconditionerWrapper&& other) noexcept
    {
        if (this != &other) {
            m_preconditioner = give(other.m_preconditioner);
            // m_stored_matrix = give(other.m_stored_matrix);
        }
        return *this;
    }

    /// Assignment from preconditioner
    gsPreconditionerWrapper& operator=(typename gsPreconditionerOp<Scalar>::uPtr preconditioner)
    {
        m_preconditioner = give(preconditioner);
        return *this;
    }

    /// Set preconditioner
    void set(typename gsPreconditionerOp<Scalar>::uPtr preconditioner)
    {
        m_preconditioner = give(preconditioner);
    }

    StorageIndex rows() const
    {
        return m_preconditioner ? gsEigen::internal::convert_index<StorageIndex>(m_preconditioner->rows())
                                : -1;
    }

    StorageIndex cols() const
    {
        return m_preconditioner ? gsEigen::internal::convert_index<StorageIndex>(m_preconditioner->cols())
                                : -1;
    }

    template<typename EigenMatType>
    gsPreconditionerWrapper& analyzePattern(const EigenMatType& )
    {
        return *this;
    }

    template<typename EigenMatType>
    gsPreconditionerWrapper& factorize(const EigenMatType& mat)
    {
        // Store a copy of the matrix for later use
        // m_stored_matrix = mat;
        return *this;
    }

    template<typename EigenMatType>
    gsPreconditionerWrapper& compute(const EigenMatType& mat)
    {
        return analyzePattern(mat).factorize(mat);
    }

    template<typename Rhs>
    inline const Vector solve(const gsEigen::MatrixBase<Rhs>& b) const
    {
        GISMO_ENSURE(m_preconditioner != nullptr, "gsPreconditionerWrapper: No preconditioner set.\n"
            "Please call solver.preconditioner() = <preconditioner> or solver.setPreconditioner(<preconditioner>)"
            " before calling solve().");

        // Convert Eigen vector to Gismo format
        const gsMatrix<Scalar> & bmat = b;
        gsMatrix<Scalar> result(b.rows(), b.cols());

        // Apply Gismo preconditioner
        m_preconditioner->apply(bmat, result);

        // Convert back to Eigen format
        return result;
    }

    gsEigen::ComputationInfo info() { return gsEigen::Success; }

    /// Access to the underlying preconditioner
    const gsPreconditionerOp<Scalar>* preconditioner() const { return m_preconditioner.get(); }
          gsPreconditionerOp<Scalar>* preconditioner()       { return m_preconditioner.get(); }

protected:
    typename gsPreconditionerOp<Scalar>::uPtr m_preconditioner;
    // NOTE: Intentionally kept commented-out. When extending gsPreconditionerWrapper
    // to support matrix-dependent Eigen preconditioners that need access to the
    // original system matrix, this member can be uncommented and used to store
    // a copy of that matrix during initialization.
    // EigenSparseMatrixType m_stored_matrix;
};

} // namespace gismo
