/** @file gsConstantBasis.h

    @brief Provides declaration of a basis of constant functions,
    consisting of one constant function.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsConstantFunction.h>

namespace gismo
{


/** 
    @brief Class defining a dummy basis of constant functions. This is
    used for compatibility reasons.

    \tparam T value type

    \ingroup basis
    \ingroup Core
*/
template <class T>
class gsConstantBasis : public gsBasis<T>
{
public:
    /// Shared pointer for gsConstantBasis
    typedef memory::shared_ptr< gsConstantBasis > Ptr;

    /// Unique pointer for gsConstantBasis
    typedef memory::unique_ptr< gsConstantBasis > uPtr;

    gsConstantBasis(T x , short_t domainDim  = 1)
    : m_val(x), m_domainDim(domainDim)
    { }

    // compatibility constructor for gsRationalBasis
    gsConstantBasis(const gsBasis<T> * b, gsMatrix<T> weight)
    : m_val( weight.value() ), m_domainDim(b->dim())
    {
        GISMO_ASSERT(weight.size() == 1, "Something seems wrong.");
    }

    // compatibility constructor for gsTensorBasis
    gsConstantBasis(const std::vector<gsBasis<T>*> & rr)
    : m_val( 1.0 ), m_domainDim(1)
    { }

    GISMO_CLONE_FUNCTION(gsConstantBasis)

    static gsConstantBasis * New(std::vector<gsBasis<T>*> & bb )
    { 
        return new gsConstantBasis(bb);
    }

public:

    short_t domainDim() const   { return m_domainDim; }

    index_t size() const   { return 1; }

    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const
    {
        result.setZero(1,u.cols());
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                                              << ", expected "<< m_domainDim);
        result.setConstant(1, u.cols(), m_val);
    }

    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                                              << ", expected "<< m_domainDim);
        result = gsMatrix<T>::Zero(m_domainDim, u.cols() );
    }

    virtual void anchors_into(gsMatrix<T>& result) const
    {
        result.setZero(1,1);
    }
    
    std::ostream &print(std::ostream &os) const
    {
        os << m_val; 
        return os; 
    }

    memory::unique_ptr<gsGeometry<T> > makeGeometry( gsMatrix<T> coefs ) const
    {
        coefs *= m_val;
        return memory::unique_ptr<gsGeometry<T> >(new gsConstantFunction<T>(coefs.transpose(), m_domainDim));
    }

public:

    T value() const { return m_val;}

private:

    T m_val;

    short_t m_domainDim;
};


} // namespace gismo
