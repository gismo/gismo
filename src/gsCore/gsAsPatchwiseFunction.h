/** @file gsAsPatchwiseFunction.h

    @brief Provides declaration of a gsAsPatchwiseFunction class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo
{

/** @brief A function depending on an index \a i, typically referring
    to a patch/sub-domain.

    \tparam T argument and value type
    
    \ingroup Core
*/

template <class T>
class gsAsPatchwiseFunction
{
private:
    gsAsPatchwiseFunction();

public:

    gsAsPatchwiseFunction(const gsFunction<T> & func, index_t numPieces)
    : m_numPieces(numPieces), m_func_ptr(func)
    { }

public:

    /// Clones the function object, making a deep copy.
    gsPatchwiseFunction * clone() const
    {
        gsPatchwiseFunction * c = new gsPatchwiseFunction<T>();
        for(size_t i=0; i<m_funcs.size();i++)
            c->addPiece(piece(i));
        return c;
    }

    index_t size() const {return m_numPieces;}
    
    const gsFunction<T> & piece(size_t i) const 
    { 
        GISMO_ASSERT(i< m_funcs.size(), "Wrong piece index");
        return m_func;
    }

    std::ostream &print(std::ostream &os) const
    {
        os << "AsPatchwise Function with "<<m_numPieces <<" pieces "<<m_func<<"\n";
        return os; 
    }

private:
    
    index_t m_numPieces;
    
    const gsFunction<T> & m_func;
};

} // namespace gismo

