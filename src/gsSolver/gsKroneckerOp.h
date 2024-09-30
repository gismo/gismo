/** @file gsKroneckerOp.h

    @brief Provides a linear operator representing the Kronecker product of linear operators

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/
#pragma once

#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/// @brief Class for representing a Kronecker product of operators of type \a gsLinearOperator
///
/// \f$ op = op_0 \otimes op_1 \otimes \cdots \otimes op_{n-1} \f$
///
/// where \f$ A \otimes B = ( a_{11} B \  a_{12} B \ ... ;  a_{21} B \  a_{22} B \ ... ; ... ) \f$.
///
/// \ingroup Solver
template <class T>
class gsKroneckerOp GISMO_FINAL : public gsLinearOperator<T>
{
    typedef typename gsLinearOperator<T>::Ptr BasePtr;
public:

    /// Shared pointer for gsKroneckerOp
    typedef memory::shared_ptr<gsKroneckerOp> Ptr;

    /// Unique pointer for gsKroneckerOp
    typedef memory::unique_ptr<gsKroneckerOp> uPtr;

    /// Kronecker product of a given list of operators
    gsKroneckerOp(std::vector< BasePtr > ops)
        : m_ops(give(ops))
    { }

    /// Convenience constructor for Kronecker product of two linear operators
    gsKroneckerOp(BasePtr op0, BasePtr op1)
        : m_ops(2)
    { m_ops[0] = give(op0); m_ops[1] = give(op1); }

    /// Convenience constructor for Kronecker product of three linear operators
    gsKroneckerOp(BasePtr op0, BasePtr op1, BasePtr op2)
        : m_ops(3)
    { m_ops[0] = give(op0); m_ops[1] = give(op1); m_ops[2] = give(op2); }

    /// Make command returning a smart pointer
    static uPtr make(std::vector< BasePtr > ops)
    { return uPtr( new gsKroneckerOp(give(ops)) ); }

    /// Make command returning a smart pointer
    /// Convenience function for Kronecker product of two linear operators
    static uPtr make(BasePtr op0, BasePtr op1)
    { return uPtr( new gsKroneckerOp(give(op0),give(op1)) ); }

    /// Make command returning a smart pointer
    /// Convenience function for Kronecker product of three linear operators
    static uPtr make(BasePtr op0, BasePtr op1, BasePtr op2)
    { return uPtr( new gsKroneckerOp(give(op0),give(op1),give(op2)) ); }

    /// Add another operator at the end
    ///
    /// \f$ op_{new} = op_{old} \otimes op \f$
    void addOperator( BasePtr op )
    { m_ops.push_back(give(op)); }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & x) const;
    virtual index_t rows() const;
    virtual index_t cols() const;

    /// Return a vector of shared pointers to all operators
    const std::vector<BasePtr>& getOps() const { return m_ops; }

    /// Apply provided linear operators without the need of creating an object
    static void apply(const std::vector<BasePtr> & ops, const gsMatrix<T> & input, gsMatrix<T> & x);

private:
    std::vector<BasePtr> m_ops;
};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsKroneckerOp.hpp)
#endif
