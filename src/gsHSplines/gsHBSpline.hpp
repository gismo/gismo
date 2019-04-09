/** @file gsHBSpline.hpp

    @brief Provides implementation of gsHBSpline class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>

namespace gismo
{

namespace internal
{


/// Get a HBSpline from XML data
template<class T, short_t d>
class gsXml< gsHBSpline<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsHBSpline<TMPLA2(d,T)>);
    static std::string tag () { return "Geometry"; }
    static std::string type () { return "HBSpline"+to_string(d); }

    static gsHBSpline<d,T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsHBSpline<d,T> >(node);
    }
    
    static gsXmlNode * put (const gsHBSpline<d,T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml< gsHBSpline<d,T> >(obj,data);
    }
};


} // end namespace internal

} // end namespace gismo
