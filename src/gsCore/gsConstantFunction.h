/** @file gsConstantFunction.h

    @brief Provides declaration of ConstantFunction class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunction.h>

namespace gismo
{

/** 
    @brief Class defining a constant function

    \tparam T value type
*/

template <class T>
class gsConstantFunction : public gsFunction<T>
{
public:
    explicit gsConstantFunction(const gsVector<T>& val, int domainDim = 1)
        : m_val(val), m_domainDim(domainDim)
    { }

    gsConstantFunction(T x, int domainDim  = 1)
        : m_domainDim(domainDim)
    {
        m_val.resize(1);
        m_val(0) = x;
    }

    gsConstantFunction(T x, T y, int domainDim)
        : m_domainDim(domainDim)
    {
        m_val.resize(2);
        m_val(0) = x;
        m_val(1) = y;
    }

    gsConstantFunction(T x, T y, T z, int domainDim)
        : m_domainDim(domainDim)
    {
        m_val.resize(3);
        m_val(0) = x;
        m_val(1) = y;
        m_val(2) = z;
    }

    virtual gsConstantFunction * clone() const { return new gsConstantFunction(*this); }

    virtual int domainDim() const   { return m_domainDim ; }
    virtual int targetDim() const   { return m_val.rows(); }
    virtual int       dim() const   { return m_domainDim ; }
        
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result = m_val.rowwise().replicate( u.cols() );
    }

    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        result = gsMatrix<T>::Zero( this->targetDim(), this->domainDim() * u.cols() );
    }

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
        os << m_val; return os; 
    }
  
private:

    gsVector<T> m_val;

    int m_domainDim;
};

}
