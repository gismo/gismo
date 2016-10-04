/** @file gsSurface.h

    @brief Provides declaration of Surface abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsGeometryEvaluator.h>

namespace gismo
{

/** 
    \brief
    Abstract base class representing a surface.

    \ingroup geometry
    \ingroup Core
*/

template<class T>
class gsSurface : public gsGeometry<T>
{

public: 
  /// Shared pointer for gsSurface
    typedef typename memory::shared< gsSurface >::ptr Ptr;
//  typedef memory::unique_ptr< gsSurface > LocalPtr;

    typedef T Scalar_t;
public:

    /// @name Constructors
    /// @{

    /// Default empty constructor
    gsSurface() : gsGeometry<T>() { }
    
    /// Constructor which copies the given coefficient matrix \a
    /// coefs.
    gsSurface(const gsBasis<T> & basis, const gsMatrix<T> & coefs ) :
    gsGeometry<T>(basis, coefs )
    { 
        GISMO_ASSERT( this->m_coefs.size() >= 1,
        "Coefficient matrix cannot be empty.\n");
        // GISMO_ASSERT( coefs.cols() >= 2, 
        // "Surface must be embedded in dimension at least two.\n");
    }

    /// Constructor which takes ownership of the given coefficient
    /// matrix \a coefs.
    gsSurface(const gsBasis<T> & basis, gsMovable< gsMatrix<T> > coefs ) :
    gsGeometry<T>(basis, coefs )
    { 
        GISMO_ASSERT( this->m_coefs.size() >= 1,
        "Coefficient matrix cannot be empty.\n");
        // GISMO_ASSERT( this->m_coefs.size() >= 2,
        // "Surface must be embedded in dimension at least two.\n");
    }

    /// @}

    virtual gsSurface * clone() const = 0;

    void toMesh(gsMesh<T> & msh, int npoints = 625) const;

    virtual gsGeometryEvaluator<Scalar_t> * evaluator(unsigned flags) const;

}; // class gsSurface

} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsSurface.hpp)
#endif

