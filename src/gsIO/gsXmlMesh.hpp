/** @file gsXmlMesh.hpp

    @brief XML serialization for the gsMesh types (implementation
    header).

    Moved here from gsIO/gsXmlUtils.hpp: each module owns its gsXml
    specializations (modularization stream S3, step A6).

    ONLY include this file from gsMeshXml_.cpp (lib mode) or through
    gsMeshXml_.cpp in header-only mode - a second inclusion in another
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
#include <gsUtils/gsMesh/gsMesh.h>

namespace gismo {
namespace internal {

template<class T>
class gsXml< gsMesh<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsMesh<T>);
    static std::string tag () { return "Mesh"; }
    static std::string type () { return ""; } // no inheritance

    static gsMesh<T> * get (gsXmlNode * node)
    {
        gsMesh<T> * m = new gsMesh<T>;

        GISMO_ASSERT( ( !strcmp( node->name(),"Mesh") ), "Problem in name of the node." );
        if ( !strcmp(node->first_attribute("format")->value(),"off") )
        {
            std::istringstream str;
            str.str( node->value() );

            std::string line;
            getline(str, line);
            if ( line.compare(0,3,"OFF") != 0)
                return nullptr;

            std::istringstream lnstream;
            getline(str, line);
            lnstream.str(line);
            unsigned nverts, nfaces, nedges(0);
            lnstream >> std::ws >>  nverts >> std::ws >> nfaces >> std::ws >> nedges ;

            T x, y, z;
            for (unsigned i=0; i<nverts; ++i)
            {
                gsGetReal(str, x);
                gsGetReal(str, y);
                gsGetReal(str, z);
                m->addVertex(x,y,z);
            }

            unsigned c = 0;
            std::vector<int> face;
            for (unsigned i=0; i<nfaces; ++i)
            {
                gsGetInt(str, c);
                face.resize(c);
                for (unsigned j=0; j<c; ++j)
                    gsGetInt(str, face[j]);
                m->addFace(face);
            }
            m->cleanMesh();
            return m;
        }
        
        gsDebug<<"Reader not implemented.\n";
        gsWarn<<"Problem in reading "<<node->first_attribute("format")->value()<<" file to gsMesh.\n";
        return m;
    }

    static gsXmlNode * put (const gsMesh<T> &,
                            gsXmlTree & )
    {
        return nullptr;
    }
};
} // namespace internal
} // namespace gismo
