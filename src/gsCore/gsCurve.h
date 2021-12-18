/** @file gsCurve.h

    @brief Provides declaration of Curve abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsGeometry.h>

namespace gismo
{

/** 
    \brief
    Abstract base class representing a curve.

    \ingroup geometry
    \ingroup Core
*/

template<class T>
class gsCurve : public gsGeometry<T>
{
    
public:
    /// Shared pointer for gsCurve
    typedef memory::shared_ptr< gsCurve > Ptr;

    /// Unique pointer for gsCurve
    typedef memory::unique_ptr< gsCurve > uPtr;
    
    typedef T Scalar_t;
public:

    /// @name Constructors
    /// @{

    /// Default empty constructor
    gsCurve() : gsGeometry<T>() { }

    /// Constructor which copies the given coefficient matrix \a
    /// coefs.
    gsCurve(const gsBasis<T> & basis, gsMatrix<T> coefs ) :
    gsGeometry<T>(basis, give(coefs))
    {
        GISMO_ASSERT( this->m_coefs.size() >= 1,
        "Coefficient matrix cannot be empty.\n");
    }

    /// @}

//    GISMO_UPTR_FUNCTION_PURE(gsCurve, clone)
private: virtual gsCurve * clone_impl() const = 0;
public: inline uPtr clone() const { return uPtr(clone_impl()); }

    short_t domainDim() const { return 1; }
    
    short_t degree() const;

    void toMesh(gsMesh<T> & msh, int npoints = 100) const;

    /// Reverse the coefficients
    void reverse()
    {
        this->m_coefs = this->m_coefs.colwise().reverse().eval();
        this->basis().reverse();
    }

    virtual bool isOn(gsMatrix<T> const &, T = 1e-3) const
    { GISMO_NO_IMPLEMENTATION }

}; // class gsCurve

} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCurve.hpp)
#endif
