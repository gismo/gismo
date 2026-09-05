/** @file gsPdeXml.hpp

    @brief XML serialization for the gsPde types (implementation
    header).

    Moved here from gsIO/gsXmlUtils.hpp: each module owns its gsXml
    specializations (modularization stream S3, step A6).

    ONLY include this file from gsXmlRegistration_.cpp (lib mode) or through
    gsPdeXmlRegistration.h in header-only mode - a second inclusion in another
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
#include <gsPde/gsPoissonPde.h>
#include <gsFunctionExpr/gsFunctionExpr.h>

namespace gismo {
namespace internal {

template<class T>
class gsXml< gsPoissonPde<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsPoissonPde<T>);
    static std::string tag () { return "Pde"; }
    static std::string type () { return "PoissonPde"; }

    static gsPoissonPde<T> * get (gsXmlNode * node)
    {
        assert( ( !strcmp( node->name(),"Pde") ) &&
                ( !(
                    strcmp( node->first_attribute("type")->value(),"PoissonPde")
                    && strcmp( node->first_attribute("type")->value(),"SurfacePoissonPde")
                     )) );

        // Read the dimension
        GISMO_ASSERT( node->first_attribute("dim"), "xml reader: No dim found" ) ;
        short_t d = atoi( node->first_attribute("dim")->value() );


        unsigned tDim = 0;
        gsXmlAttribute * targetDim = node->first_attribute("targetDim");

        if ( targetDim )
            tDim = atoi( targetDim->value() );

        if ( tDim >= 1 )
        {
            gsXmlNode * tmp = node->first_node("rhs");
            gsFunctionExpr<T>  rhs_fnct;
            internal::gsXml<gsFunctionExpr<T> >::get_into(tmp, rhs_fnct);

            tmp = node->first_node("solution");
            if ( tmp )
            {
                gsFunctionExpr<T> msol;
                internal::gsXml<gsFunctionExpr<T> >::get_into(tmp, msol);

                return new gsPoissonPde<T>(rhs_fnct, d, msol );
            }
            else
            {
                return new gsPoissonPde<T>( rhs_fnct, d );
            }
        }

        // Read right hand side function
        gsXmlNode * tmp = node->first_node("rhs");
        gsFunctionExpr<T> rhs(tmp->value(), d);

        // Read exact solution, if one exists in the file
        tmp = node->first_node("solution");
        if ( tmp )
        {
            gsFunctionExpr<T> sol(tmp->value(), d);
            return new gsPoissonPde<T>(rhs, d, sol );
        }
        else
        {
            return new gsPoissonPde<T>( rhs, d );
        }
    }

    static gsXmlNode * put (const gsPoissonPde<T> &,
                            gsXmlTree & )
    {
        return NULL;
    }

};





} // namespace internal
} // namespace gismo
