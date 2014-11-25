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

/** @brief A function depending on an index \a i, typically refering
    to a patch/sub-domain.

    \tparam T argument and value type
*/

template <class T>
class gsPiecewiseFunction
{
public:
    typedef std::size_t size_t;
    typedef typename std::vector<const gsFunction<T>*>::iterator fiterator;
    typedef typename std::vector<const gsFunction<T>*>::const_iterator const_fiterator;

public:

    gsPiecewiseFunction()
    { 
        const gsConstantFunction<T> * cf = new gsConstantFunction<T>(0);
        m_funcs.resize(1, cf);
    }

    explicit gsPiecewiseFunction(size_t numPieces)
    { 
        const gsConstantFunction<T> * cf = new gsConstantFunction<T>(0);
        m_funcs.resize(numPieces, cf);
    }

    explicit gsPiecewiseFunction(const gsFunction<T> & func, size_t numPieces = 1)
    { 
        const gsFunction<T> * cf = func.clone();
        m_funcs.resize(numPieces, cf);
    }
    ~gsPiecewiseFunction()
    {
        fiterator it = std::unique(m_funcs.begin(), m_funcs.end() );
        freeAll(m_funcs.begin(), it);
    }
    

public:

    /// Add a piece
    void addPiece(const gsFunction<T> & func) const
    { 
        m_funcs.push_back( func.clone() );
    }

    
    /// Evaluate at \a i on evaluation points \a u, output at \a result.
    void eval_into(size_t i, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_funcs[i]->eval_into(u,result);
    }

    /// Derivative at \a i on evaluation points \a u, output at \a result.
    void deriv_into(size_t i, const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        m_funcs[i]->deriv_into(u,result);
    }

    const gsFunction<T> & piece(size_t i) const { return *m_funcs[i]; }

    gsFunction<T> & piece(size_t i)             { return *m_funcs[i]; }
    
    inline const gsFunction<T> & operator [] (size_t & i) const 
    {
        return *m_funcs[i];
    }
    
    inline gsFunction<T> & operator [] (size_t i) 
    {
        ///TODO: this is ugly and dangerous
        return *const_cast<gsFunction<T>*>(m_funcs[i]);
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
    
    std::vector<const gsFunction<T> * > m_funcs;
};

}

