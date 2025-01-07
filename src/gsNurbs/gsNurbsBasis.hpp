/** @file gsNurbsBasis.hpp

    @brief Implementation of 1D NURBS basis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsNurbs/gsNurbs.h>
#include <gsNurbs/gsTensorNurbs.h>

namespace gismo
{

template <class T>
typename gsNurbsBasis<T>::gsGeoPtr
gsNurbsBasis<T>::makeGeometry( gsMatrix<T> coefs ) const
{ return gsGeoPtr(new GeometryType(*this, give(coefs))); }

template <class T>
typename gsNurbsBasis<T>::gsBasisPtr
gsNurbsBasis<T>::create(std::vector<KnotVectorType> cKV, gsMatrix<T> weights)
{
    const index_t dd = cKV.size();
    switch (dd)
    {
    case 1:
        return gsBasisPtr(new gsNurbsBasis<T>(give(cKV.front()), give(weights)));
        break;
    case 2:
        return gsBasisPtr(new gsTensorNurbsBasis<2,T>(give(cKV), give(weights)));
        break;
    case 3:
        return gsBasisPtr(new gsTensorNurbsBasis<3,T>(give(cKV), give(weights)));
        break;
    case 4:
        return gsBasisPtr(new gsTensorNurbsBasis<4,T>(give(cKV), give(weights)));
        break;
    }
    GISMO_ERROR("Dimension should be between 1 and 4.");
}

} // namespace gismo
