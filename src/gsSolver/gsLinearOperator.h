/** @file gsLinearOperator.h

    @brief Simple abstract class for (discrete) linear operators.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsCore/gsLinearAlgebra.h>

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
    typedef typename memory::shared<gsLinearOperator>::ptr Ptr;

    /// Unique pointer for gsLinearOperator   
    typedef typename memory::unique<gsLinearOperator>::ptr uPtr;

    /// Identity operator
    static gsIdentityOp<T> Identity(const index_t dim) {return gsIdentityOp<T>(dim); }
    
    virtual ~gsLinearOperator() {}

    /**
     * @brief apply the operator on the input vector and store the result in x
     * @param input Input vector
     * @param x     result vector
     */
    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const = 0;

    ///Returns the number of rows of the operator
    virtual index_t rows() const = 0;

    ///Returns the number of columns of the operator
    virtual index_t cols() const = 0;

    // NOTE: this is rather inefficient and is only provided for debugging and testing purposes
    void toMatrix(gsMatrix<T>& result)
    {
        GISMO_ASSERT(rows() == cols(),
                     "gsLinearOperator::toMatrix is only implemented for square operators");

        gsMatrix<T> eye = gsMatrix<T>::Identity(cols(), cols());
        this->apply(eye, result);
    }
}; // gsLinearOperator

/// @brief Allows an operator to be multiplied with a scalar
template<class T>
class gsScaledOp : public gsLinearOperator<T>
{
public:
    /// Shared pointer for gsScaledOp
    typedef typename memory::shared<gsScaledOp>::ptr Ptr;

    /// Unique pointer for gsScaledOp
    typedef typename memory::unique<gsScaledOp>::ptr uPtr;

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


/// @brief Identity operator, must be square!
template<class T>
class gsIdentityOp : public gsLinearOperator<T>
{
public:

    /// Shared pointer for gsIdentityOp
    typedef typename memory::shared<gsIdentityOp>::ptr Ptr;

    /// Unique pointer for gsIdentityOp   
    typedef typename memory::unique<gsIdentityOp>::ptr uPtr;
    
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
