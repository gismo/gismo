/** @file gsGenericGeometry.h

    @brief Provides declaration of GenericGeometry abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsCurve.h>
#include <gsCore/gsSurface.h>
#include <gsCore/gsVolume.h>
#include <gsCore/gsBulk.h>

namespace gismo
{

// Generic traits for geometry with dimension known at compile time
template <unsigned d, typename T>
struct gsGeoTraits
{
    typedef gsGeometry<T> GeometryBase;
};

// Traits for curve
template <typename T>
struct gsGeoTraits<1,T>
{
    typedef gsCurve<T> GeometryBase;
};

// Traits for surface
template <typename T>
struct gsGeoTraits<2,T>
{
    typedef gsSurface<T> GeometryBase;
};

// Traits for volume
template <typename T>
struct gsGeoTraits<3,T>
{
    typedef gsVolume<T> GeometryBase;
};

// Traits for bulk
template <typename T>
struct gsGeoTraits<4,T>
{
    typedef gsBulk<T> GeometryBase;
};


/** \brief  Base class for geometries with statically known basis types.
 *
 * This is an internal implementation class which all concrete geometries
 * derive from. It allows the compile-time specification of the associated
 * basis type. Note that the basis dimension is also statically known.
 *
 * \tparam Basis_t  type of the basis for this geometry
 *
 *
 * \ingroup Core
 */
template<class Basis_t>
class gsGenericGeometry : 
        public gsGeoTraits<Basis_t::Dim, 
                           typename Basis_t::Scalar_t>::GeometryBase
{
public: 
    typedef typename Basis_t::Scalar_t Scalar_t;

    typedef typename gsGeoTraits<Basis_t::Dim,Scalar_t>::GeometryBase Base;

public:

    // Look at gsGeometry base constructor for a brief description
    gsGenericGeometry() : Base() { }

    // Look at gsGeometry base constructor for a brief description
    gsGenericGeometry(const Basis_t & basis, const gsMatrix<Scalar_t> & coefs )
    : Base (basis,coefs)
    {
        //todo: gsBasis as constructor argument, dynamic cast in assertion

        GISMO_ASSERT( basis.size() == coefs.rows(), 
                      "The coefficient matrix of the geometry (rows="<<coefs.rows()<<") does not match the number of basis functions in its basis("<< basis.size() <<").");
    }
    
    // Look at gsGeometry base constructor for a brief description
    gsGenericGeometry(const Basis_t & basis, gsMovable< gsMatrix<Scalar_t> > coefs )
    : Base (basis,coefs)
    { 
        GISMO_ASSERT( basis.size() == this->m_coefs.rows(), 
                      "The coefficient matrix of the geometry (rows="<<this->m_coefs.rows()<<") does not match the number of basis functions in its basis("<< basis.size() <<").");
    }
        
public:
    
    // for rational bases, we WILL store the coefficients in projective form later
    bool isProjective() const
    { return Basis_t::IsRational; }
    
    // Look at gsGeometry base constructor for a brief description
    Basis_t & basis()       { return static_cast<Basis_t&>(*this->m_basis); }
    
    // Look at gsGeometry base constructor for a brief description
    const Basis_t & basis() const { return static_cast<const Basis_t&>(*this->m_basis); }
    
}; // class gsGenericGeometry



} // namespace gismo
