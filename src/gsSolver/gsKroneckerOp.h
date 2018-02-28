/** @file gsKroneckerOp.h

    @brief Provides a linear operator representing the Kronecker product of linear operators

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/
#pragma once

#include <gsCore/gsExport.h>
#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/// @brief Class for representing a Kronecker product of operators of type \a gsLinearOperator
///
/// \ingroup Solver
template <class T>
class gsKroneckerOp : public gsLinearOperator<T>
{
    typedef typename gsLinearOperator<T>::Ptr BasePtr;
public:

    /// Shared pointer for gsKroneckerOp
    typedef memory::shared_ptr<gsKroneckerOp> Ptr;

    /// Unique pointer for gsKroneckerOp
    typedef memory::unique_ptr<gsKroneckerOp> uPtr;

    /// Kronecker product of a given list of operators
    gsKroneckerOp(const std::vector< BasePtr >& ops)
        : m_ops(ops)
    { }

    /// Convenience constructor for Kronecker product of two linear operators
    gsKroneckerOp(const BasePtr & op0, const BasePtr & op1)
        : m_ops(2)
    { m_ops[0] = op0; m_ops[1] = op1; }

    /// Convenience constructor for Kronecker product of three linear operators
    gsKroneckerOp(const BasePtr & op0, const BasePtr & op1, const BasePtr & op2)
        : m_ops(3)
    { m_ops[0] = op0; m_ops[1] = op1; m_ops[2] = op2; }

    /// Make command returning a smart pointer
    static uPtr make(const std::vector< BasePtr >& ops)
    { return uPtr( new gsKroneckerOp(ops) ); }

    /// Make command returning a smart pointer
    /// Convenience function for Kronecker product of two linear operators
    static uPtr make(const BasePtr & op0, const BasePtr & op1)
    { return uPtr( new gsKroneckerOp(op0,op1) ); }

    /// Make command returning a smart pointer
    /// Convenience function for Kronecker product of three linear operators
    static uPtr make(const BasePtr & op0, const BasePtr & op1, const BasePtr & op2)
    { return uPtr( new gsKroneckerOp(op0,op1,op2) ); }

    /// Add another operator at the end
    void addOperator( const BasePtr& op )
    { m_ops.push_back( op ); }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & result) const;
    virtual index_t rows() const;
    virtual index_t cols() const;

    /// Return a vector of shared pointers to all operators
    const std::vector<BasePtr>& getOps() const { return m_ops; }

    /// Apply provided linear operators without the need of creating an object
    static void applyKronecker(const std::vector<BasePtr> &, const gsMatrix<T> & input, gsMatrix<T> & result);

private:
    std::vector<BasePtr> m_ops;
};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsKronecker.hpp)
#endif
