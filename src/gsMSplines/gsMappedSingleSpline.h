/** @file gsMappedSingleSpline.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsCore/gsGeometry.h>
#include <gsMSplines/gsMappedSpline.h>
#include <gsMSplines/gsMappedSingleBasis.h>

namespace gismo
{

template<short_t d,class T> class gsMappedSpline;

/**
   @brief Class gsMappedSingleSpline represents an individual .....of a

   Note that it does not own the underlying spline: if you delete
   the object the spline is not deleted. The lifetime of the
   underlying spline should be at least the lifetime of gsMappedSingleSpline.

*/
template<short_t d,class T>
class gsMappedSingleSpline : public gsGeoTraits<d,T>::GeometryBase
{
private:
    typedef typename gsGeoTraits<d,T>::GeometryBase Base;

    typedef gsMappedSingleBasis<d,T> Basis;

public:
    /// Shared pointer for gsMappedSingleSpline
    typedef memory::shared_ptr< gsMappedSingleSpline > Ptr;

    /// Unique pointer for gsMappedSingleSpline
    typedef memory::unique_ptr< gsMappedSingleSpline > uPtr;

private:
    /// Default empty constructor
    gsMappedSingleSpline() : m_spline(nullptr) { }

public:

    /// Construct a spline function by a pointer to a spline and an index i
    gsMappedSingleSpline(gsMappedSpline<d,T> * spline, unsigned const & i = 0)
    : m_spline(spline),m_index(i)
    {
        GISMO_ASSERT( i<unsigned(m_spline->nPatches()),"Invalid spline function index" );
    }

    gsMappedSingleSpline( const gsMappedSingleSpline& other ) : Base( other )
    {
        m_spline = other.m_spline;
        m_index = other.m_index;
    }

    static uPtr make(   const gsMappedSingleSpline& other)
    { return uPtr( new gsMappedSingleSpline( other ) ); }

    ~gsMappedSingleSpline() { } //destructor

    GISMO_BASIS_ACCESSORS

public:

    short_t domainDim() const
    {
        return d;
    }

    /// Evaluates the non-zero spline functions at value u.
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const
    {
        // m_spline->evalGlobal_into(m_index,u,result);
        m_spline->eval_into(m_index,u,result);
    }

    /// Evaluates the (partial) derivatives of non-zero spline functions at (the columns of) u.
    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
    {
        m_spline->deriv_into(m_index,u,result);
    }

    /// Evaluates the (partial) derivatives of the nonzero spline functions at points \a u into \a result.
    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result ) const
    {
        m_spline->deriv2_into(m_index,u,result);
    }

    /// @brief Evaluate the nonzero spline functions and their derivatives up
    /// to order \a n at points \a u into \a result.
    void evalAllDers_into(const gsMatrix<T> & u, int n, std::vector<gsMatrix<T> >& result) const
    {
        m_spline->evalAllDers_into(m_index,u,n,result);
    }

    GISMO_CLONE_FUNCTION(gsMappedSingleSpline)

    std::ostream &print(std::ostream &os) const
    {
        GISMO_UNUSED(os);
        GISMO_NO_IMPLEMENTATION;
    }

    //////////////////////////////////////////////////
    // Virtual member that may be implemented or not by the derived class
    //////////////////////////////////////////////////

    /// Returns the polynomial degree.
    short_t degree(short_t i) const
    {
        return m_spline->degree(m_index,i);
    }

    /// The gsBasisFun points to the i-th spline function of m_spline
    /// after calling this setter.
    void setPiece( unsigned const & i ) const
    {
        GISMO_ASSERT( i<unsigned(m_spline->nPatches()),"Invalid spline index" );
        m_index = i;
    }
    
// Data members
private:
    gsMappedSpline<d,T> * m_spline;
    mutable index_t m_index;

}; // class gsMappedSingleSpline


} // namespace gismo
