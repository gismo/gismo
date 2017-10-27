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
class gsConstantFunction : public gsGeometry<T>
{
public:
    typedef gsGeometry<T> Base;

    /// Shared pointer for gsConstantFunction
    typedef memory::shared_ptr< gsConstantFunction > Ptr;

    /// Unique pointer for gsConstantFunction
    typedef memory::unique_ptr< gsConstantFunction > uPtr;

    /// Returns a null function 
    static const gsConstantFunction Zero(int domDim, int tarDim)
    { return gsConstantFunction(gsVector<T>::Zero(tarDim),domDim); }

    gsConstantFunction() { }

    explicit gsConstantFunction(const gsVector<T>& val, int domainDim)
    :  m_domainDim(domainDim)
    {
        m_coefs = val.transpose();
    }


    ///  Constructs a constant function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R \f$
    explicit gsConstantFunction(T x, int domainDim)
        : m_domainDim(domainDim)
    {
        m_coefs.resize(1,1);
        m_coefs.at(0) = x;
    }

    /// Constructs a constant function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^2 \f$
    gsConstantFunction(T x, T y, int domainDim)
        : m_domainDim(domainDim)
    {
        m_coefs.resize(1,2);
        m_coefs.at(0) = x;
        m_coefs.at(1) = y;
    }

    /// Constructs a constant Function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^3 \f$
    gsConstantFunction(T x, T y, T z, int domainDim)
        : m_domainDim(domainDim)
    {
        m_coefs.resize(1,3);
        m_coefs.at(0) = x;
        m_coefs.at(1) = y;
        m_coefs.at(2) = z;
    }

    /// Constructs a constant Function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^4 \f$
    gsConstantFunction(T x, T y, T z, T w,  int domainDim)
        : m_domainDim(domainDim)
    {
        m_coefs.resize(1,4);
        m_coefs.at(0) = x;
        m_coefs.at(1) = y;
        m_coefs.at(2) = z;
        m_coefs.at(2) = w;
    }

    /// Compatibility constructor
    gsConstantFunction(const gsConstantBasis<T> & cb, const gsMatrix<T> & coef)
    : m_domainDim(1)
    {
        m_coefs = cb.value()*coef;
    }

    GISMO_CLONE_FUNCTION(gsConstantFunction, virtual)

    const gsConstantFunction<T> & piece(const index_t k) const
    {
        // same on all pieces
        GISMO_UNUSED(k);
        return *this; 
    }

    // Documentation in gsFunction class
    virtual int domainDim() const   { return m_domainDim ; }

    // Documentation in gsFunction class
    virtual int targetDim() const   { return m_coefs.cols(); }

    const gsVector<T> & value() const { return m_coefs.transpose();}

    T value(size_t i) const { return m_coefs.at(i);}

    void setValue(T val, int domainDim)
    { m_coefs.setConstant(1,1,val); m_domainDim = domainDim;}

    void setValue(const gsVector<T> & val, int domainDim)
    { m_coefs = val.transpose(); m_domainDim = domainDim;}

    // Documentation in gsFunction class
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                                              << ", expected "<< m_domainDim);
        result = m_coefs.transpose().rowwise().replicate( u.cols() );
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
        os << m_coefs.transpose(); 
        return os; 
    }

    virtual const gsBasis<T> & basis() const {GISMO_NO_IMPLEMENTATION}
    virtual gsBasis<T> & basis() {GISMO_NO_IMPLEMENTATION}

    void compute(const gsMatrix<T> & in, gsFuncData<T> & out) const
    { gsFunction<T>::compute(in, out); }
    
private:

    /// Global value of this function
    using Base::m_coefs;

    /// Spatial dimension of the domain of definition of this function
    int m_domainDim;
};

}
