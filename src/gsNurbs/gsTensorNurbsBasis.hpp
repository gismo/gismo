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

} // namespace gismo
