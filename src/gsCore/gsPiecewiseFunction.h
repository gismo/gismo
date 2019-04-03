/** @file gsPiecewiseFunction.h

    @brief Provides declaration of a gsPiecewiseFunction class.

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
    to a patch/sub-domain. On each patch a different gsFunction object is used.

    \tparam T argument and value type
    
    \ingroup Core
*/

template <class T>
class gsPiecewiseFunction : public gsFunctionSet<T>
{
public:
    typedef gsFunctionSet<T>                           Base;
    typedef typename std::vector<gsFunction<T>*>       FunctionContainer;
    typedef typename FunctionContainer::iterator       fiterator;
    typedef typename FunctionContainer::const_iterator const_fiterator;

    /// Shared pointer for gsPiecewiseFunction
    typedef memory::shared_ptr< gsPiecewiseFunction > Ptr;

    /// Unique pointer for gsPiecewiseFunction
    typedef memory::unique_ptr< gsPiecewiseFunction > uPtr;

public:

    explicit gsPiecewiseFunction(index_t npieces = 0) : gsFunctionSet<T>()
    { m_funcs.reserve(2+npieces); }

    gsPiecewiseFunction(const gsFunction<T> & func) : gsFunctionSet<T>(func)
    {
        m_funcs.push_back(func.clone().release());
        //m_funcs.resize(n, func.clone());
    }

    gsPiecewiseFunction(const gsPiecewiseFunction & other) : gsFunctionSet<T>(other)
    {
        m_funcs.resize(other.m_funcs.size() );
        cloneAll( other.m_funcs.begin(), other.m_funcs.end(),
                  m_funcs.begin() );
    }

    gsPiecewiseFunction(FunctionContainer & funcs) : gsFunctionSet<T>()
    {
        m_funcs.swap(funcs); // funcs are consumed
    }

    ~gsPiecewiseFunction()
    {
        freeAll(m_funcs);
    }

    #if EIGEN_HAS_RVALUE_REFERENCES
    /// Move constructor
    gsPiecewiseFunction(gsPiecewiseFunction&& other)
    : m_funcs(give(other.m_funcs)) {}

    /// Assignment operator
    gsPiecewiseFunction& operator= ( const gsPiecewiseFunction& other )
    {
        freeAll(m_funcs);
        m_funcs.resize(other.m_funcs.size() );
        cloneAll( other.m_funcs.begin(), other.m_funcs.end(),
                  m_funcs.begin() );
        return *this;
    }

    /// Move assignment operator
    gsPiecewiseFunction& operator= ( gsPiecewiseFunction&& other )
    {
        freeAll(m_funcs);
        m_funcs = give(other.m_funcs);
        return *this;
    }
#else
    /// Assignment operator (uses copy-and-swap idiom)
    gsPiecewiseFunction & operator= ( gsPiecewiseFunction other )
    {
        this->swap( other );
        return *this;
    }
#endif

    /// \brief Swap with another gsPiecewiseFunction
    void swap(gsPiecewiseFunction & other)
    {
        m_funcs.swap( other.m_funcs );
    }

    GISMO_CLONE_FUNCTION(gsPiecewiseFunction)

    short_t domainDim () const {return m_funcs.front()->domainDim();};
    short_t targetDim () const {return m_funcs.front()->targetDim();};
    
    /// Add a piece
    void addPiece(const gsFunction<T> & func)
    { 
        m_funcs.push_back( func.clone().release() );
    }

    void addPiecePointer(gsFunction<T> * func)
    { 
        m_funcs.push_back( func );
    }

    void addPiecePointer(typename gsFunction<T>::uPtr func)
    {
        m_funcs.push_back(func.release());
    }

    const gsFunction<T> & piece(const index_t i) const 
    { 
        GISMO_ASSERT(static_cast<size_t>(i) < m_funcs.size(), "Wrong piece index");
        return *m_funcs[i]; 
    }

    index_t size() const {return m_funcs.size();}

    std::ostream &print(std::ostream &os) const
    {
        os << "Piecewise Function with "<<m_funcs.size() <<" pieces.\n";
        return os; 
    }

    friend std::ostream & operator<<(std::ostream & os, const gsPiecewiseFunction & pwf)
    {
        return pwf.print(os);
    }

    /// Clear (delete) all functions
    void clear()
    {
        freeAll(m_funcs);
        m_funcs.clear();
    }

protected:
    
    FunctionContainer m_funcs;
};

} // namespace gismo

