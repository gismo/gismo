/** @file gsNurbsXml.hpp

    @brief XML serialization for the gsNurbs types (implementation
    header).

    Moved here from gsIO/gsXmlUtils.hpp: each module owns its gsXml
    specializations (modularization stream S3, step A6).

    ONLY include this file from gsXmlRegistration_.cpp (lib mode) or through
    gsNurbsXmlRegistration.h in header-only mode - a second inclusion in another
    translation unit of the library would poison the symbol visibility
    of the explicit instantiations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>
#include <gsNurbs/gsNurbs.h>
#include <gsNurbs/gsTensorNurbs.h>
#include <gsNurbs/gsTensorNurbsBasis.h>

namespace gismo {
namespace internal {

template<short_t d, class T>
class gsXml< gsTensorNurbsBasis<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTensorNurbsBasis<TMPLA2(d,T)>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return "TensorNurbsBasis"+to_string(d); }

    static gsTensorNurbsBasis<d,T> * get (gsXmlNode * node)
    {
        return getRationalBasisFromXml< gsTensorNurbsBasis<d,T> >(node);
    }

    static gsXmlNode * put (const gsTensorNurbsBasis<d,T> & obj,
                            gsXmlTree & data )
    {
        return putRationalBasisToXml< gsTensorNurbsBasis<d,T> >(obj,data);
    }
};

template<class T>
class gsXml< gsNurbs<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsNurbs<T>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return "Nurbs"; }

    static gsNurbs<T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsNurbs<T> >(node);
    }

    static gsXmlNode * put (const gsNurbs<T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml(obj,data);
    }
};

template<short_t d, class T>
class gsXml< gsTensorNurbs<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsTensorNurbs<TMPLA2(d,T)>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return "TensorNurbs"+to_string(d); }

    static gsTensorNurbs<d,T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsTensorNurbs<d,T> >( node );
    }

    static gsXmlNode * put (const gsTensorNurbs<d,T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml(obj,data);
    }
};

} // namespace internal
} // namespace gismo
