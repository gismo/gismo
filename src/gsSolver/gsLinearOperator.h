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

}; // gsLinearOperator


/// @brief Identity preconditioner ("do nothing"), must be square!
class gsIdentityPreconditioner : public gsLinearOperator
{
public:

    gsIdentityPreconditioner(index_t dim) : m_dim(dim) {}


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
