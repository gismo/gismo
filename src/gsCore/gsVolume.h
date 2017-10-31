/** @file gsVolume.h

    @brief Provides declaration of Volume abstract interface.

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
    Provides declaration of Volume abstract interface.

    \ingroup geometry
    \ingroup Core
*/

template<class T>
class gsVolume : public gsGeometry<T>
{

public:
    /// Shared pointer for gsVolume
    typedef memory::shared_ptr< gsVolume > Ptr;

    /// Unique pointer for gsVolume
    typedef memory::unique_ptr< gsVolume > uPtr;

    typedef T Scalar_t;
public:

    /// @name Constructors
    /// @{

    /// Default empty constructor
    gsVolume() : gsGeometry<T>() { }
    
    /// Constructor which copies the given coefficient matrix \a
    /// coefs.
    gsVolume(const gsBasis<T> & basis, gsMatrix<T> coefs ) :
    gsGeometry<T>(basis, give(coefs) )
    { 
        GISMO_ASSERT( this->m_coefs.size() >= 1,
        "Coefficient matrix cannot be empty.\n");
        // GISMO_ASSERT( coefs.cols() >= 2, 
        // "Surface must be embedded in dimension at least two.\n");
    }

    /// @}

    //GISMO_CLONE_FUNCTION_FORWARD(gsVolume)
private: virtual gsVolume * doClone() const = 0; \
public: inline uPtr clone() const { return uPtr(doClone()); }

    int domainDim() const { return 3; }

    void toMesh(gsMesh<T> & msh, int npoints = 3375) const;

    virtual gsGeometryEvaluator<Scalar_t> * evaluator(unsigned flags) const;

}; // class gsVolume

} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsVolume.hpp)
#endif

