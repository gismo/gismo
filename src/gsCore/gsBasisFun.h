/** @file gsBasisFun.h

    @brief Provides definition of the BasisFun class.

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

template<class T> class gsBasis;

/** 
    @brief Represents an individual function in a function set, or a
           certain component of a vector-valued function

    It is constructed using a function set and the index of a single
    function in that basis.  Note that it does not own the underlying
    function set (eg. basis). The lifetime of the parent object should
    be at least the lifetime of gsBasisFun.

    \ingroup function
    \ingroup Core
*/  
template<class T>
class gsBasisFun : public gsFunction<T>
{
private:
    /// Default empty constructor
    gsBasisFun() { }
    
public:
    /// Shared pointer for gsBasisFun
    typedef memory::shared_ptr< gsBasisFun > Ptr;

    /// Unique pointer for gsBasisFun
    typedef memory::unique_ptr< gsBasisFun > uPtr;

    /// Construct a basis function by a pointer to a basis and an index i
    gsBasisFun(const gsBasis<T> & basis, unsigned const i);

    ~gsBasisFun() { } //destructor

    GISMO_CLONE_FUNCTION(gsBasisFun)

public:
  
    short_t domainDim () const {return m_basis.domainDim();}

    short_t targetDim () const {return m_basis.targetDim();}

    gsMatrix<T> support() const;

    void eval_into (const gsMatrix<T>& u, gsMatrix<T>& result ) const;

    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result ) const;

    /// The gsBasisFun points to the i-th basis function of m_basis
    /// after calling this setter.
    void setFunction( unsigned const & i );

    /// Point to the first basis function of the basis
    void first();

    /// Points to the first basis function of the basis
    /// or return false if this is the last basis function
    bool next();

    /// Return false if the the function iteration is invalidated
    bool valid();

    unsigned index() const { return m_index; }
    
    // temporary hack
    virtual const gsBasisFun & piece(const index_t) const
    {
        return *this; 
    }


// Data members
private:
    const gsBasis<T> & m_basis;
    unsigned m_index;

}; // class gsBasisFun


template<class T>
gsBasisFun<T>::gsBasisFun(const gsBasis<T> & basis, const unsigned i )
: m_basis(basis), m_index(i) 
{ 
    GISMO_ASSERT( i<unsigned(m_basis.size()),"Invalid basis function index" );
}

template<class T> void
gsBasisFun<T>::setFunction( unsigned const & i )
{
    GISMO_ASSERT( i<unsigned(m_basis.size()),"Invalid basis function index" );
    m_index = i;
}

template<class T> void
gsBasisFun<T>::first()
{
    m_index = 0;
}

template<class T> bool
gsBasisFun<T>::next()
{
    return ( ++m_index  < static_cast<unsigned>(m_basis.size()) );
}

template<class T> bool
gsBasisFun<T>::valid()
{
    return ( m_index  < static_cast<unsigned>(m_basis.size()) );
}

template<class T> gsMatrix<T>
gsBasisFun<T>::support()  const
{
    return m_basis.support(m_index);
}

template<class T> void
gsBasisFun<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result )  const
{
    m_basis.evalSingle_into(m_index, u, result);
}

template<class T> void
gsBasisFun<T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result )  const
{
    m_basis.derivSingle_into(m_index, u, result);
}


} // namespace gismo
