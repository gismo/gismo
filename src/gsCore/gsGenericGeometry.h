/** @file gsGenericGeometry.h

    @brief Provides declaration of GenericGeometry abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo
{


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
//template<unsigned d, typename T>
template<typename Basis_t>
class gsGenericGeometry : 
        public gsGeoTraits<Basis_t::Dim, 
                           typename Basis_t::Scalar_t>::GeometryBase
//public gsGeoTraits<d,T>::GeometryBase
{
public: 
    typedef typename Basis_t::Scalar_t Scalar_t;

    typedef Basis_t Basis;

    typedef typename gsGeoTraits<Basis_t::Dim,Scalar_t>::GeometryBase Base;

public:

    // Look at gsGeometry base constructor for a brief description
    gsGenericGeometry() : Base() { }

    // Look at gsGeometry base constructor for a brief description
    gsGenericGeometry(const Basis_t & basis, const gsMatrix<Scalar_t> & coefs )
    : Base (basis,coefs)
    { }
    
    // Look at gsGeometry base constructor for a brief description
    gsGenericGeometry(const Basis_t & basis, gsMovable< gsMatrix<Scalar_t> > coefs )
    : Base (basis,coefs)
    { }
        
public:

    GISMO_BASIS_ACCESSORS
    
}; // class gsGenericGeometry



} // namespace gismo
