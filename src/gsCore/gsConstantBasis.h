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
//#include <gsCore/gsBasis.h>

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
    gsConstantBasis(T x , int domainDim  = 1)
    : m_val(x), m_domainDim(domainDim)
    { }

    // compatibility constructor for gsRationalBasis
    gsConstantBasis(const gsBasis<T> * b, const gsMovable<gsMatrix<T> >& weight)
    : m_val( 1.0 ), m_domainDim(b->dim())
    {
        GISMO_ASSERT(weight.ref().size() == 1, "Something seems wrong.");
    }

    // compatibility constructor for gsTensorBasis
    gsConstantBasis(const std::vector<gsBasis<T>*> & rr)
    : m_val( 1.0 ), m_domainDim(1)
    { }
      
    gsConstantBasis * clone() const { return new gsConstantBasis(*this); }

public:

    int dim() const   { return m_domainDim; }

    int size() const   { return 1; }

    void active_into(const gsMatrix<T> & u, gsMatrix<unsigned>& result) const
    {
        result.setZero(1,u.cols());
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                                              << ", expected "<< m_domainDim);
        result.setConstant(m_val);
    }

    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                                              << ", expected "<< m_domainDim);
        result = gsMatrix<T>::Zero(m_domainDim, u.cols() );
    }


    std::ostream &print(std::ostream &os) const
    {
        os << m_val; 
        return os; 
    }

    gsGeometry<T> * makeGeometry( const gsMatrix<T> & coefs )      const { return NULL; }
    gsGeometry<T> * makeGeometry( gsMovable< gsMatrix<T> > coefs ) const { return NULL; }

public:

    T value() const { return m_val;}

private:

    T m_val;

    int m_domainDim;
};


} // namespace gismo
