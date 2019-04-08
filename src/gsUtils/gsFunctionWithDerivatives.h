/** @file gsFunctionWithDerivatives.h

    @brief A function with explicitly set derivatives

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <gsCore/gsFunction.h>


#pragma once

namespace gismo {

template <typename T=real_t>
class gsFunctionWithDerivatives : public gsFunction<T>
{
protected:
    const gsFunction<T> *m_values;
    const gsFunction<T> *m_derivs;
    const gsFunction<T> *m_deriv2;
public:
    /// Shared pointer for gsFunctionWithDerivatives
    typedef memory::shared_ptr< gsFunctionWithDerivatives > Ptr;

    /// Unique pointer for gsFunctionWithDerivatives
    typedef memory::unique_ptr< gsFunctionWithDerivatives > uPtr;

    gsFunctionWithDerivatives()
    : m_values(NULL),
      m_derivs(NULL),
      m_deriv2(NULL)
    { }

    gsFunctionWithDerivatives(const gsFunction<T> &func, const gsFunction<T> &deriv)
    : m_values(&func),
      m_derivs(&deriv),
      m_deriv2(NULL)
    {
        GISMO_ASSERT(checkDimensions(), "Dimensions do not fit");
    }

    gsFunctionWithDerivatives(const gsFunction<T> &func,
                              const gsFunction<T> &deriv,
                              const gsFunction<T> &deriv2)
    : m_values(&func  ),
      m_derivs(&deriv ),
      m_deriv2(&deriv2)
    {
        GISMO_ASSERT(checkDimensions(), "Dimensions do not fit");
    }

    void eval_into(const gsMatrix<T> &u, gsMatrix<T> &result) const
    {
        GISMO_ASSERT(m_derivs,"Not initialized");
        m_values->eval_into(u,result);
    }
    
    void deriv_into(const gsMatrix<T> &u, gsMatrix<T> &result) const
    {
        GISMO_ASSERT(m_derivs,"No first derivative available");
        m_derivs->eval_into(u,result);
    }
    
    void deriv2_into(const gsMatrix<T> &u, gsMatrix<T> &result) const
    {
        GISMO_ASSERT(m_deriv2,"No second derivative available");
        m_deriv2->eval_into(u,result);
    }

    gsMatrix<T> laplacian(const gsMatrix<T> &u) const
    {
        gsMatrix<T>  secDer;
        deriv2_into(u,secDer);
        return secDer.topRows(m_values->domainDim()).colwise().sum();
    }

    const gsFunction<T> & getDeriv () const
    {
        return *m_derivs;
    }
    const gsFunction<T> & getDeriv2 () const
    {
        return *m_deriv2;
    }

    const gsFunctionWithDerivatives & piece(const index_t) const
    {
        // same on all pieces
        return *this; 
    }

    GISMO_UPTR_FUNCTION_NO_IMPLEMENTATION(gsFunction<T>, clone)
 
    short_t targetDim () const
    {
        return m_values->targetDim();
    }
    short_t domainDim () const
    {
        return m_values->domainDim();
    }


    /*
    bool check()
    {
        if (m_derivs)
        {

        }
        if (m_deriv2)
        {

        }
    }
    */
    
private:

    gsFunctionWithDerivatives(const gsFunctionWithDerivatives<T> &func);
        
    bool checkDimensions()
    {
        bool ok=true;

        const short_t parDim=m_values->domainDim();
        const short_t tarDim=m_values->targetDim();

        if (m_derivs)
        {
            ok = ok &&  m_derivs->domainDim() == parDim;
            ok = ok &&  m_derivs->targetDim() == parDim*tarDim;
        }
        if (m_deriv2)
        {
            ok = ok &&  m_deriv2->domainDim() == parDim;
            ok = ok &&  m_deriv2->targetDim() == (parDim+1)*parDim*tarDim/2;
        }
        return ok;
    }
};


}
