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
    @brief Class defining a globally constant function

    \tparam T value type

    \ingroup function
    \ingroup Core
*/

template <class T>
class gsConstantFunction : public gsFunction<T>
{
public:
    gsConstantFunction() { }

    explicit gsConstantFunction(const gsVector<T>& val, int domainDim)
        : m_val(val), m_domainDim(domainDim)
    { }


    ///  Constructs a constant function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R \f$
    explicit gsConstantFunction(T x, int domainDim)
        : m_domainDim(domainDim)
    {
        m_val.resize(1);
        m_val(0) = x;
    }

    /// Constructs a constant function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^2 \f$
    gsConstantFunction(T x, T y, int domainDim)
        : m_domainDim(domainDim)
    {
        m_val.resize(2);
        m_val(0) = x;
        m_val(1) = y;
    }

    /// Constructs a constant Function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^3 \f$
    gsConstantFunction(T x, T y, T z, int domainDim)
        : m_domainDim(domainDim)
    {
        m_val.resize(3);
        m_val(0) = x;
        m_val(1) = y;
        m_val(2) = z;
    }

    /// Constructs a constant Function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^4 \f$
    gsConstantFunction(T x, T y, T z, T w,  int domainDim)
        : m_domainDim(domainDim)
    {
        m_val.resize(4);
        m_val(0) = x;
        m_val(1) = y;
        m_val(2) = z;
        m_val(3) = w;
    }

    /// Compatibility constructor
    gsConstantFunction(const gsConstantBasis<T> & cb, const gsMatrix<T> & coef)
    : m_val( cb.value()*coef ), m_domainDim(1)
    { }


    // Documentation in gsFunction class
    virtual gsConstantFunction * clone() const { return new gsConstantFunction(*this); }

    // Documentation in gsFunction class
    virtual int domainDim() const   { return m_domainDim ; }

    // Documentation in gsFunction class
    virtual int targetDim() const   { return m_val.rows(); }

    const gsVector<T> & value() const { return m_val;}

    T value(size_t i) const { return m_val[i];}

    void setValue(T val) { m_val.setConstant(val);}

    void setValue(const gsVector<T> & val) { m_val = val;}

    // Documentation in gsFunction class
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                                              << ", expected "<< m_domainDim);
        result = m_val.rowwise().replicate( u.cols() );
    }

    // Documentation in gsFunction class
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                                              << ", expected "<< m_domainDim);
        result = gsMatrix<T>::Zero( this->targetDim()*this->domainDim(), u.cols() );
    }

    // Documentation in gsFunction class
    virtual void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                                              << ", expected "<< m_domainDim);
        result = gsMatrix<T>::Zero( (this->domainDim()*(this->domainDim()+1))/2,
                                    this->targetDim()*u.cols() );
    }

    // Documentation in gsFunction class
    virtual std::ostream &print(std::ostream &os) const
    {
        os << m_val; 
        return os; 
    }

private:

    /// Global value of this function
    gsVector<T> m_val;

    /// Spatial dimension of the domain of definition of this function
    int m_domainDim;
};

}
