/** @file gsLinearOperator.h

    @brief Simple abstract class for (discrete) linear operators.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, A. Manzaflaris, C. Hofreither, S. Takacs, C. Hofer
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

/// @brief Simple abstract class for discrete operators.
///
/// Simple abstract class for discrete operators.
/// The derived classes have to contain the functions: apply(), cols(), and rows().
///
/// \ingroup Solver
template<class T>
class gsLinearOperator
{
public:

    /// Shared pointer for gsLinearOperator
    typedef typename memory::shared_ptr<gsLinearOperator> Ptr;

    /// Unique pointer for gsLinearOperator   
    typedef memory::unique_ptr<gsLinearOperator> uPtr;

    /// Identity operator
    static gsIdentityOp<T> Identity(const index_t dim) {return gsIdentityOp<T>(dim); }
    
    virtual ~gsLinearOperator() {}

    /**
     * @brief apply the operator on the input vector and store the result in x
     * @param input Input vector
     * @param x     result vector
     */
    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const = 0;

    /// Returns the number of rows of the operator
    virtual index_t rows() const = 0;

    /// Returns the number of columns of the operator
    virtual index_t cols() const = 0;

    // NOTE: this is rather inefficient and is only provided for debugging and testing purposes
    void toMatrix(gsMatrix<T>& result)
    {
        GISMO_ASSERT(rows() == cols(),
                     "gsLinearOperator::toMatrix is only implemented for square operators");

        gsMatrix<T> eye = gsMatrix<T>::Identity(cols(), cols());
        this->apply(eye, result);
    }
    
    /// Get the default options as gsOptionList object
    // This implementation provides an empty object
    static gsOptionList defaultOptions()               { return gsOptionList(); }

    /// Set options based on a gsOptionList object
    // This implementation does not read any input
    virtual void setOptions(const gsOptionList & opt)  {}
    
}; // gsLinearOperator

/// @brief Simple abstract class for Steppable discrete operators.
///
/// The class represents an iteration method in the form x_new = x_old + P^{-1}(f - A*x_old).
/// The member function \a apply represents the application of P^{-1}.
/// The member function \a step represents one step of the above iteration method.
///
/// If \a setNumOfSweeps is used to set the number of sweeps to some tau>0, the member
/// function \a apply represents I-(I-P^{-1}A)^{tau} A^{-1}.
///
/// Usually, the step operation can be performed in an optimized way.
///
/// The derived classes have to contain the functions: step(), cols(), and rows().
///
/// \ingroup Solver
template<class T>
class gsSteppableOperator : public gsLinearOperator<T>
{
public:

    /// Shared pointer for gsLinearOperator
    typedef typename memory::shared_ptr<gsSteppableOperator> Ptr;

    /// Unique pointer for gsLinearOperator
    typedef memory::unique_ptr<gsSteppableOperator> uPtr;

    /// Base class
    typedef gsLinearOperator<T> Base;
    
    gsSteppableOperator() : m_num_of_sweeps(1) {}
    
    virtual ~gsSteppableOperator() {}

    /**
     * @brief apply the method on given right hand side and current iterate
     * @param rhs Right hand side vector
     * @param x   Current iterate vector
     */
    virtual void step(const gsMatrix<T> & rhs, gsMatrix<T> & x) const = 0;

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        x.setZero(this->rows(),input.cols()); // we assume quadratic matrices
        for ( index_t i = 0; i < m_num_of_sweeps; ++i )
            step(input,x);
    }
    
    /// Set the number of sweeps to be applied in the member function \a apply
    void setNumOfSweeps( index_t n )
    {
        GISMO_ASSERT ( n > 0, "Number of sweeps needs to be positive." );
        m_num_of_sweeps = n;
    }
    
    /// Get the number of sweeps to be applied in the member function \a apply
    index_t getNumOfSweeps()
    {
        return m_num_of_sweeps;
    }

    /// Get the default options as gsOptionList object
    static gsOptionList defaultOptions()
    {
        gsOptionList opt = Base::defaultOptions();
        opt.addInt( "NumOfSweeps", "Number of sweeps to be applied in the member function \a apply", 1 );
        return opt;
    }

    /// Set options based on a gsOptionList object
    virtual void setOptions(const gsOptionList & opt)
    {
        Base::setOptions(opt);
        m_num_of_sweeps = opt.askInt( "NumOfSweeps", m_num_of_sweeps );
    }
    

protected:
    index_t m_num_of_sweeps;
    
}; // gsSteppableOperator

/// @brief Allows an operator to be multiplied with a scalar
///
/// \ingroup Solver
template<class T>
class gsScaledOp : public gsLinearOperator<T>
{
public:
    /// Shared pointer for gsScaledOp
    typedef typename memory::shared_ptr<gsScaledOp> Ptr;

    /// Unique pointer for gsScaledOp
    typedef memory::unique_ptr<gsScaledOp> uPtr;

    /// Shared pointer for gsLinearOperator
    typedef typename gsLinearOperator<T>::Ptr BasePtr;

    /// Constructor taking a shared pointer to a linear operator and a scalar
    gsScaledOp(const BasePtr & linOp, T scalar = 1) : m_linOp(linOp), m_scalar(scalar)    {}

    /// Make command returing a shared pointer
    static Ptr make(const BasePtr & linOp, T scalar = 1) 
    { return memory::make_shared( new gsScaledOp(linOp, scalar) ); }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        m_linOp->apply(input, x);
        x *= m_scalar;
    }

    ///Returns the number of rows in the preconditioner
    index_t rows() const {return m_linOp->rows();}

    ///Returns the number of columns in the preconditioner
    index_t cols() const {return m_linOp->cols();}

private:
    const BasePtr m_linOp;
    const T m_scalar;
}; // gsScaladOp


/// @brief Identity operator
///
/// \ingroup Solver
template<class T>
class gsIdentityOp : public gsLinearOperator<T>
{
public:

    /// Shared pointer for gsIdentityOp
    typedef typename memory::shared_ptr<gsIdentityOp> Ptr;

    /// Unique pointer for gsIdentityOp   
    typedef memory::unique_ptr<gsIdentityOp> uPtr;
    
    /// Constructor taking the dimension of the identity operator
    gsIdentityOp(index_t dim) : m_dim(dim) {}

    /// Make command returing a shared pointer
    static Ptr make(index_t dim) { return memory::make_shared( new gsIdentityOp(dim) ); }

    void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const
    {
        x = input;
    }

    index_t rows() const {return m_dim;}

    index_t cols() const {return m_dim;}

private:
    const index_t m_dim;
};

} // namespace gismo
