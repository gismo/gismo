
#pragma once

#include <gsCore/gsGeometry.h>
#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

/** 
    \brief
    Abstract base class representing a surface.

    \ingroup geometry
*/

template<class T>
class gsBulk : public gsGeometry<T>
{

public: 
  /// Shared pointer for gsBulk
    typedef memory::shared_ptr< gsBulk > Ptr;
//  typedef memory::unique_ptr< gsBulk > LocalPtr;

    typedef T Scalar_t;
public:

    /// @name Constructors
    /// @{

    /// Default empty constructor
    gsBulk() : gsGeometry<T>() { }
    
    /// Constructor which copies the given coefficient matrix \a
    /// coefs.
    gsBulk( const gsMatrix<T> & coefs ) :
        gsGeometry<T>( coefs )
    { 
        GISMO_ASSERT( this->m_coefs.size() >= 1,
        "Coefficient matrix cannot be empty.\n");
        // GISMO_ASSERT( coefs.cols() >= 2, 
        // "Surface must be embedded in dimension at least two.\n");
    }

    /// Constructor which takes ownership of the given coefficient
    /// matrix \a coefs.
    gsBulk( gsMovable< gsMatrix<T> > coefs ) :
        gsGeometry<T>( coefs )
    { 
        GISMO_ASSERT( this->m_coefs.size() >= 1,
        "Coefficient matrix cannot be empty.\n");
        // GISMO_ASSERT( this->m_coefs.size() >= 2,
        // "Surface must be embedded in dimension at least two.\n");
    }

    /// @}

    virtual gsBulk * clone() const = 0;

    void toMesh(gsMesh<T> & msh, int npoints = 3375) const;

    virtual gsGeometryEvaluator<Scalar_t> * evaluator(unsigned flags) const;

}; // class gsBulk

}; // namespace gismo


#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsBulk.hpp)
#endif
