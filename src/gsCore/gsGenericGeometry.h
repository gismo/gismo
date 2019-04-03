/** @file gsGenericGeometry.h

    @brief Provides declaration of gsGenericGeometry class

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

/**
   \brief The generic geometry class represents a function defined as
   coefficients times basis function defined in a basis.

   This is a generic implementation with minimal functionality. For
   common geometry types (B-splines, tensor-product B-splines,
   hierarchical splines) specific classes are implemented, which
   provide additional functionalities.

 */
#include <gsCore/gsGeometry.h>

namespace gismo
{

template<short_t d, class T>
class gsGenericGeometry : public gsGeoTraits<d,T>::GeometryBase
{
public:
    typedef gsBasis<T> Basis;

    typedef typename gsGeoTraits<d,T>::GeometryBase Base;

    /// Shared pointer for gsGenericGeometry
    typedef memory::shared_ptr< gsGenericGeometry > Ptr;

    /// Unique pointer for gsGenericGeometry
    typedef memory::unique_ptr< gsGenericGeometry > uPtr;

public:
    gsGenericGeometry(const gsBasis<T> & basis, 
                      gsMatrix<T> coefs)
    : Base (basis, give(coefs))
    { 
        GISMO_ASSERT( this->m_basis->dim() == static_cast<int>(d), 
                      "Incoherent basis dimension in gsGenericGeometry");
    }

    GISMO_CLONE_FUNCTION(gsGenericGeometry)

    GISMO_BASIS_ACCESSORS

};

} // namespace gismo

