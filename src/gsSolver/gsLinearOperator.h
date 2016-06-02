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
class gsLinearOperator
{
public:

    /// Shared pointer for gsLinearOperator
    typedef memory::shared_ptr< gsLinearOperator > Ptr;

    /// Unique pointer for gsLinearOperator   
    typedef memory::unique< gsLinearOperator >::ptr uPtr;
    
    virtual ~gsLinearOperator() {}

    /**
     * @brief apply the operator on the input vector and store the result in x
     * @param input Input vector
     * @param x     result vector
     */
    virtual void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const = 0;

    ///Returns the number of rows in the preconditioner
    virtual index_t rows() const = 0;

    ///Returns the number of columns in the preconditioner
    virtual index_t cols() const = 0;

    // NOTE: this is rather inefficient and is only provided for debugging and testing purposes
    void toMatrix(gsMatrix<real_t>& result)
    {
        GISMO_ASSERT(rows() == cols(), "gsLinearOperator::toMatrix is only implemented for square operators");

        gsMatrix<real_t> eye = gsMatrix<real_t>::Identity(cols(), cols());
        this->apply(eye, result);
    }
}; // gsLinearOperator


/// @brief Identity preconditioner ("do nothing"), must be square!
class gsIdentityOp : public gsLinearOperator
{
public:

    /// Shared pointer for gsIdentityOp
    typedef memory::shared_ptr< gsIdentityOp > Ptr;

    /// Unique pointer for gsIdentityOp   
    typedef memory::unique< gsIdentityOp >::ptr uPtr;
    
    
    gsIdentityOp(index_t dim) : m_dim(dim) {}

    static Ptr make(index_t dim) { return shared( new gsIdentityOp(dim) ); }

    void apply(const gsMatrix<real_t> & input, gsMatrix<real_t> & x) const
    {
        x = input;
    }

    index_t rows() const {return m_dim;}

    index_t cols() const {return m_dim;}

private:
    const index_t m_dim;
};

} // namespace gismo
