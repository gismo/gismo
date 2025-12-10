/** @file gsRationalTHBSpline.h

    @brief Represents a rational truncated hierarchical B-Spline patch

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, C. Karampatzakis
*/
 
#pragma once

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>

namespace gismo
{
namespace internal
{

/// Get a RationalTHBSPline from XML data
template<short_t d, class T>
class gsXml< gsRationalTHBSpline<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsRationalTHBSpline<TMPLA2(d,T)>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return "RationalTHBSpline"+to_string(d); }

    static gsRationalTHBSpline<d,T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsRationalTHBSpline<d,T> >( node );
    }

    static gsXmlNode * put (const gsRationalTHBSpline<d,T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml(obj,data);
    }
};


} // namespace internal
} // namespace gismo
