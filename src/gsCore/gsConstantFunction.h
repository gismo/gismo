/** @file gsConstantFunction.h

    @brief Provides declaration of ConstantFunction class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/

#pragma once

#include <gsCore/gsGeometry.h>

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
    static const gsConstantFunction Zero(short_t domDim, short_t tarDim)
    { return gsConstantFunction(gsVector<T>::Zero(tarDim),domDim); }

    /// Returns a uPtr to a null function
    static uPtr makeZero(short_t domDim, short_t tarDim)
    { return uPtr(new gsConstantFunction(domDim, tarDim)); }

    gsConstantFunction() { }

    /// Constructs a constant function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^{\text{dim(val)}} \f$
    gsConstantFunction(const gsVector<T>& val, short_t domainDim);


    ///  Constructs a constant function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R \f$
    gsConstantFunction(T x, short_t domainDim)
        : m_domainDim(domainDim)
    {
        m_coefs.resize(1,1);
        m_coefs(0,0) = x;
    }

    /// Constructs a constant function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^2 \f$
    gsConstantFunction(T x, T y, short_t domainDim);

    /// Constructs a constant Function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^3 \f$
    gsConstantFunction(T x, T y, T z, short_t domainDim);

    /// Constructs a constant Function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^4 \f$
    gsConstantFunction(T x, T y, T z, T w,  short_t domainDim);

    /// Compatibility constructor
    gsConstantFunction(const gsConstantBasis<T> & cb, const gsMatrix<T> & coef);

    /// Copy constructor
    gsConstantFunction(const gsConstantFunction<T> & o);

    /// Move constructor
    gsConstantFunction(gsConstantFunction<T> && o);

    /// Assignment operator
    gsConstantFunction<T> & operator=(const gsConstantFunction<T> & o);

    /// Move assignment operator
    gsConstantFunction<T> & operator=(gsConstantFunction<T> && o);

    /// Constructs a constant function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^{\text{dim(val)}} \f$
    static uPtr make(const gsVector<T>& val, short_t domainDim)
    { return uPtr(new gsConstantFunction(val, domainDim)); }

    ///  Constructs a constant function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R \f$
    static uPtr make(T x, short_t domainDim)
    { return uPtr(new gsConstantFunction(x, domainDim)); }

    /// Constructs a constant function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^2 \f$
    static uPtr make(T x, T y, short_t domainDim)
    { return uPtr(new gsConstantFunction(x, y, domainDim)); }

    /// Constructs a constant Function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^3 \f$
    static uPtr make(T x, T y, T z, short_t domainDim)
    { return uPtr(new gsConstantFunction(x, y, z, domainDim)); }

    /// Constructs a constant Function \f$ \mathbb R^{\text{domainDim}} \to \mathbb R^4 \f$
    static uPtr make(T x, T y, T z, T w,  short_t domainDim)
    { return uPtr(new gsConstantFunction(x, y, z, w, domainDim)); }

    GISMO_CLONE_FUNCTION(gsConstantFunction)

    const gsConstantFunction<T> & piece(const index_t) const override
    {
        // same on all pieces
        return *this;
    }

    // Documentation in gsFunction class
    virtual short_t domainDim() const override { return m_domainDim ; }

    // Documentation in gsFunction class
    virtual short_t targetDim() const override
    { return static_cast<short_t>(m_coefs.cols()); }

    const gsVector<T> value() const { return m_coefs.transpose();}

    T value(size_t i) const { return m_coefs.at(i);}

    void setValue(T val, short_t domainDim)
    { m_coefs.setConstant(1,1,val); m_domainDim = domainDim;}

    void setValue(const gsVector<T> & val, short_t domainDim)
    { m_coefs = val.transpose(); m_domainDim = domainDim;}

    // Documentation in gsFunction class
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override;

    // Documentation in gsFunction class
    virtual void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override;

    // Documentation in gsFunction class
    virtual void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override;

    void evalAllDers_into(const gsMatrix<T> & u, int n,
                          std::vector<gsMatrix<T> > & result,
                          bool sameElement = false) const override;

    // Documentation in gsFunction class
    virtual std::ostream &print(std::ostream &os) const override
    {
        os << m_coefs.transpose();
        return os;
    }

    virtual const gsBasis<T> & basis() const override {GISMO_NO_IMPLEMENTATION}
    virtual gsBasis<T> & basis() override {GISMO_NO_IMPLEMENTATION}

    void compute(const gsMatrix<T> & in, gsFuncData<T> & out) const override
    { gsFunction<T>::compute(in, out); }

private:

    /// Global value of this function
    using Base::m_coefs;

    /// Spatial dimension of the domain of definition of this function
    short_t m_domainDim;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsConstantFunction.hpp)
#endif
