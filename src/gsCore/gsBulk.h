/** @file gsBulk.h

    @brief Provides declaration of a 4D bulk.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo
{

/** 
    \brief
    Abstract base class representing a 4D bulk.

    \ingroup geometry
    \ingroup Core
*/

template<class T>
class gsBulk : public gsGeometry<T>
{

public:
    /// Shared pointer for gsBulk
    typedef memory::shared_ptr< gsBulk > Ptr;

    /// Unique pointer for gsBulk
    typedef memory::unique_ptr< gsBulk > uPtr;

    typedef T Scalar_t;
public:

    /// @name Constructors
    /// @{

    /// Default empty constructor
    gsBulk() : gsGeometry<T>() { }
    
    /// Constructor which copies the given coefficient matrix \a
    /// coefs.
    gsBulk(const gsBasis<T> & basis, gsMatrix<T> coefs ) :
    gsGeometry<T>(basis, give(coefs))
    { 
        GISMO_ASSERT( this->m_coefs.size() >= 1,
        "Coefficient matrix cannot be empty.\n");
        // GISMO_ASSERT( coefs.cols() >= 2, 
        // "Surface must be embedded in dimension at least two.\n");
    }

    /// @}

    GISMO_UPTR_FUNCTION_PURE(gsBulk, clone)

    short_t domainDim() const { return 4; }

}; // class gsBulk

} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBulk.hpp)
#endif
