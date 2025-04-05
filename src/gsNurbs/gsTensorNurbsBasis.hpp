/** @file gsTensorNurbsBasis.hpp

    @brief Implementation of d-D NURBS basis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsNurbs/gsTensorNurbs.h>

namespace gismo
{

template<short_t d, class T>
typename gsTensorNurbsBasis<d,T>::gsGeoPtr
gsTensorNurbsBasis<d,T>::makeGeometry( gsMatrix<T> coefs ) const
{ return gsGeoPtr(new GeometryType(*this, give(coefs))); }

template<short_t d, class T>
memory::unique_ptr<typename gsTensorNurbsBasis<d,T>::BoundaryBasisType>
gsTensorNurbsBasis<d,T>::boundaryBasis(boxSide const & s)
{
    typename Src_t::BoundaryBasisType::uPtr bb = m_src->boundaryBasis(s);
    gsMatrix<index_t> ind = m_src->boundary(s);
    
    gsMatrix<T> ww( ind.size(),1);
    for ( index_t i=0; i<ind.size(); ++i)
        ww(i,0) = m_weights( (ind)(i,0), 0);

    return memory::unique_ptr<BoundaryBasisType>(
        new BoundaryBasisType(bb.release(), give(ww)) );
}

} // namespace gismo
