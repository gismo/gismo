/** @file gsPiecewiseFunction.h

    @brief Provides declaration of a gsPiecewiseFunction class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFunction.h>

namespace gismo
{

template <class T> class gsFunction;

/** @brief A function depending on an index \a i, typically referring
    to a patch/sub-domain.

    \tparam T argument and value type
    
    \ingroup Core
*/

template <class T>
class gsPiecewiseFunction
{
public:
    typedef std::size_t size_t;
    typedef typename std::vector<gsFunction<T>*>::iterator fiterator;
    typedef typename std::vector<gsFunction<T>*>::const_iterator const_fiterator;

public:

    gsPiecewiseFunction() { }
    
    explicit gsPiecewiseFunction(const gsFunction<T> & func, size_t numPieces = 1)
    { 
        gsFunction<T> * cf = func.clone();
        m_funcs.resize(numPieces, cf);
    }

    ~gsPiecewiseFunction()
    {
        fiterator it = std::unique(m_funcs.begin(), m_funcs.end() );
        freeAll(m_funcs.begin(), it);
    }    

public:

    /// Add a piece
    void addPiece(const gsFunction<T> & func)
    { 
        m_funcs.push_back( func.clone() );
    }
    
    /// Evaluate at \a i on evaluation points \a u, output at \a result.
    void eval_into(size_t i, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(i< m_funcs.size(), "Wrong piece index");
        m_funcs[i]->eval_into(u,result);
    }

    /// Derivative at \a i on evaluation points \a u, output at \a result.
    void deriv_into(size_t i, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(i< m_funcs.size(), "Wrong piece index");
        m_funcs[i]->deriv_into(u,result);
    }

    const gsFunction<T> & piece(size_t i) const 
    { 
        GISMO_ASSERT(i< m_funcs.size(), "Wrong piece index");
        return *m_funcs[i]; 
    }

    inline const gsFunction<T> & operator [] (size_t & i) const 
    {
        GISMO_ASSERT(i< m_funcs.size(), "Wrong piece index");
        return *m_funcs[i];
    }
    
    inline gsFunction<T> & operator [] (size_t i) 
    {
        GISMO_ASSERT(i<m_funcs.size(), "Wrong piece index");
        return *m_funcs[i];
    }

    std::ostream &print(std::ostream &os) const
    {
        os << "Piecewise Function with "<<m_funcs.size() <<" pieces.\n";
        return os; 
    }

    friend std::ostream & operator<<(std::ostream &os, const gsPiecewiseFunction &pwf)
    {
        return pwf.print(os);
    }

private:
    
    std::vector<gsFunction<T> * > m_funcs;
};

}

