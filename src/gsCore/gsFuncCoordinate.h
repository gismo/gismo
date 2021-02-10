/** @file gsFuncCoordinate.h

    @brief Provides definition of the FuncCoordinate class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, A. Limkilde
*/

#pragma once

#include <gsCore/gsFunction.h>
#include <gsCore/gsFunctionSet.h>

namespace gismo
{

/**
    @brief Represents a certain component of a vector-valued function

    It is constructed using a vector valued function and a coordinate.
    Note that it does not own the underlying function.
    The lifetime of the parent object should
    be at least the lifetime of gsFuncCoordinate.

    \ingroup function
    \ingroup Core
*/
template<class T>
class gsFuncCoordinate : public gsFunction<T>
{
public:
    /// Default empty constructor
    gsFuncCoordinate() { }

public:
    /// Shared pointer for gsFuncCoordinate
    typedef memory::shared_ptr< gsFuncCoordinate > Ptr;

    /// Unique pointer for gsFuncCoordinate
    typedef memory::unique_ptr< gsFuncCoordinate > uPtr;

    /// Construct a basis function by a pointer to a basis and an index i
    gsFuncCoordinate(const gsFunctionSet<T> & function, unsigned const i);

    ~gsFuncCoordinate() { } //destructor

    GISMO_CLONE_FUNCTION(gsFuncCoordinate)

public:

    short_t domainDim () const {return m_function->domainDim();}

    short_t targetDim () const {return 1;} // It is a coordinate of a vector function

    gsMatrix<T> support() const;

    void eval_into (const gsMatrix<T>& u, gsMatrix<T>& result ) const;

    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result ) const;

    void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result ) const;

    /// The gsFuncCoordinate points to the i-th coordinate
    /// after calling this setter.
    void setCoordinate( unsigned const & i );

    /// Point to the first coordinate
    void first();

    /// Points to the first coordinate
    /// or return false if this is the last coordinate
    bool next();

    /// Return false if the the function iteration is invalidated
    bool valid();

    unsigned index() const { return m_index; }

    // // temporary hack
    virtual const gsFuncCoordinate & piece(const index_t k) const
    {
        m_pieces[k] = gsFuncCoordinate(m_function->piece(k),m_index);
        return m_pieces[k];
    }


// Data members
private:
    const gsFunctionSet<T> * m_function;
    unsigned m_index;

    mutable std::map<index_t,gsFuncCoordinate> m_pieces;
}; // class gsFuncCoordinate


template<class T>
gsFuncCoordinate<T>::gsFuncCoordinate(const gsFunctionSet<T> & function, const unsigned i )
: m_function(&function), m_index(i)
{
    GISMO_ASSERT( i<unsigned(m_function->targetDim()),"Invalid coordinate" );
}

template<class T> void
gsFuncCoordinate<T>::setCoordinate( unsigned const & i )
{
    GISMO_ASSERT( i<unsigned(m_function->targetDim()),"Invalid coordinate" );
    m_index = i;
}

template<class T> void
gsFuncCoordinate<T>::first()
{
    m_index = 0;
}

template<class T> bool
gsFuncCoordinate<T>::next()
{
    return ( ++m_index  < static_cast<unsigned>(m_function->targetDim()) );
}

template<class T> bool
gsFuncCoordinate<T>::valid()
{
    return ( m_index  < static_cast<unsigned>(m_function->targetDim()) );
}

template<class T> gsMatrix<T>
gsFuncCoordinate<T>::support()  const
{
    return m_function->support();
}

template<class T> void
gsFuncCoordinate<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result )  const
{
    gsMatrix<T> tmp;
    m_function->eval_into(u, tmp);
    result = tmp.row(m_index);
}

template<class T> void
gsFuncCoordinate<T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result )  const
{
    gsMatrix<T> tmp;
    m_function->deriv_into(u, tmp);
    const index_t stride = domainDim();
    result = tmp.middleRows(m_index*stride,stride);
}

template<class T> void
gsFuncCoordinate<T>::deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result )  const
{
    gsInfo << "deriv2 is called!\n";
    gsMatrix<T> tmp;
    m_function->deriv2_into(u, tmp);
    const index_t dim = domainDim();
    const index_t stride = dim*(dim+1)/2;
    result = tmp.middleRows(m_index*stride,stride);
}


} // namespace gismo
