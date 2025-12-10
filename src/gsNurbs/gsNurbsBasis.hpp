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
#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>

namespace gismo
{

template <class T>
gsTensorNurbsBasis<1,T>::~gsTensorNurbsBasis() { }

template <class T>
typename gsTensorNurbsBasis<1,T>::gsGeoPtr
gsTensorNurbsBasis<1,T>::makeGeometry( gsMatrix<T> coefs ) const
{ return gsGeoPtr(new GeometryType(*this, give(coefs))); }

template <class T>
typename gsTensorNurbsBasis<1,T>::gsBasisPtr
gsTensorNurbsBasis<1,T>::create(std::vector<KnotVectorType> cKV, gsMatrix<T> weights)
{
    const index_t dd = cKV.size();
    switch (dd)
    {
    case 1:
        return gsBasisPtr(new gsTensorNurbsBasis<1,T>(give(cKV.front()), give(weights)));
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

namespace internal {
/// Get a NurbsBasis from XML data
template<class T>
class gsXml< gsNurbsBasis<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsNurbsBasis<T>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return "NurbsBasis"; }

    static gsNurbsBasis<T> * get (gsXmlNode * node)
    {
        return getRationalBasisFromXml<gsNurbsBasis<T> >(node);
    }

    static gsXmlNode * put (const gsNurbsBasis<T> & obj,
                            gsXmlTree & data )
    {
        return putRationalBasisToXml(obj,data);
    }
};

}

} // namespace gismo
