/** @file gsRationalTHBSplineBasis.h

    @brief Provides declaration of RationalTHBSplineBasis abstract interface.

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
    
/// Get a RationalTHBSplineBasis from XML data
template<short_t d, class T>
class gsXml< gsRationalTHBSplineBasis<d,T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsRationalTHBSplineBasis<TMPLA2(d,T)>);
    static std::string tag () { return "Basis"; }
    static std::string type () { return "RationalTHBSplineBasis"+to_string(d); }

    static gsRationalTHBSplineBasis<d,T> * get (gsXmlNode * node)
    {
        return getRationalBasisFromXml< gsRationalTHBSplineBasis<d,T> >(node);
    }

    static gsXmlNode * put (const gsRationalTHBSplineBasis<d,T> & obj,
                            gsXmlTree & data )
    {
        return putRationalBasisToXml< gsRationalTHBSplineBasis<d,T> >(obj,data);
    }
};

} // namespace internal
} // namespace gismo
