/** @file gsCoreXml.hpp

    @brief XML serialization for the gsCore container types (implementation
    header).

    Moved here from gsIO/gsXmlUtils.hpp: each module owns its gsXml
    specializations (modularization stream S3, step A6).

    ONLY include this file from gsCoreXml_.cpp (lib mode) or through
    gsCoreXmlRegistration.h in header-only mode - a second inclusion in another
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
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsMultiBasis.h>

namespace gismo {
namespace internal {

template<class T>
class gsXml< gsMultiPatch<T> >
{
private:
    gsXml() { }
    typedef gsMultiPatch<T> Object;

public:
    GSXML_COMMON_FUNCTIONS(Object);
    static std::string tag () { return "MultiPatch"; }
    static std::string type () { return ""; }

    GSXML_GET_POINTER(Object);

    static void get_into (gsXmlNode * node, Object & obj)
    {
        GISMO_ASSERT( !strcmp( node->name(),"MultiPatch"),
                      "Something went wrong. Expected Multipatch tag." );

        // the geometry patches should be siblings of node
        gsXmlNode * toplevel = node->parent();

        const int d = atoi( node->first_attribute("parDim")->value() );

        gsXmlNode * tmp = node->first_node("patches");
        std::istringstream str ;
        str.str( tmp->value() );

        std::vector< gsGeometry<T> *> patches;
        std::map<int,int> ids;
        if ( ! strcmp( tmp->first_attribute("type")->value(),"id_range") )
        {
            int first, last;
            gsGetInt(str, first);
            gsGetInt(str, last);
            for ( int i = first; i<=last; ++i )
            {
                GISMO_ASSERT( searchId(i, toplevel) != NULL,
                              "No Geometry with Id "<<i<<" found in the XML data.");
                patches.push_back( getById< gsGeometry<T> >( toplevel, i ) );
                patches.back()->setId(i);
                ids[i] = i - first;
            }
        }
        else if ( ! strcmp( tmp->first_attribute("type")->value(),"id_index") )
        {
            int c = 0;
            for (int pindex; gsGetInt(str, pindex);)
            {
                GISMO_ASSERT( searchId(pindex, toplevel) != NULL,
                              "No Geometry with Id "<<pindex<<" found in the XML data.");
                patches.push_back( getById< gsGeometry<T> >( toplevel, pindex ) );
                patches.back()->setId(pindex);
                ids[pindex] = c++;
            }
        }
        else
        {
            gsWarn<<"Unknown tag in XML multipatch object.\n";
        }


        //patches: 2 0 1
        // before offset range: 5 3 4

        // Boundaries and interfaces are also 3,4,5 so we need to translate them t0 0,1,2

        // Read boundary
        std::vector< patchSide > boundaries;
        for (gsXmlNode * child = node->first_node("boundary"); child;
                child = child->next_sibling("boundary"))
        {
            std::vector< patchSide > tmp_boundaries;
            if (child)
            {
                getBoundaries(child, ids, tmp_boundaries);
                boundaries.insert( boundaries.end(), tmp_boundaries.begin(), tmp_boundaries.end() );
            }
        }
        // Remove duplicates (keeps the first one)
        std::sort(boundaries.begin(), boundaries.end());
        boundaries.erase(std::unique(boundaries.begin(), boundaries.end()), boundaries.end());

        // Read interfaces
        std::vector< boundaryInterface > interfaces;
        for (gsXmlNode * child = node->first_node("interfaces"); child;
                child = child->next_sibling("interfaces"))
        {
            std::vector< boundaryInterface > tmp_interfaces;
            if (child)
            {
                getInterfaces(child, d, ids, tmp_interfaces);
                interfaces.insert( interfaces.end(), tmp_interfaces.begin(), tmp_interfaces.end() );
            }
        }
        // Remove duplicates (keeps the first one)
        std::sort(interfaces.begin(), interfaces.end());
        interfaces.erase(std::unique(interfaces.begin(), interfaces.end()), interfaces.end());


        obj = gsMultiPatch<T>(patches, boundaries, interfaces);
    }

    static gsXmlNode * put (const gsMultiPatch<T> & obj,
                            gsXmlTree & data)
    {
        // First insert all geometries
        int max_id = data.maxId();
        gsXmlNode * tmp;
        std::map<index_t, index_t> id_map;
        for ( typename gsMultiPatch<T>::const_iterator it = obj.begin();
              it != obj.end(); ++it )
        {
            tmp = gsXml<gsGeometry<T> >::put(**it,data);
            data.appendToRoot(tmp);
            id_map[obj.findPatchIndex(*it)] =  std::stoi(tmp->first_attribute("id")->value());
        }

        std::ostringstream str;
        str<< max_id+1 <<" "<< data.maxId();
        tmp = internal::makeNode("patches" , str.str(), data);
        tmp->append_attribute( internal::makeAttribute("type", "id_range", data) );
        str.clear(); str.str("");

        // Make MultiPatch node
        gsXmlNode * mp_node = internal::makeNode("MultiPatch" , data);
        mp_node->append_attribute( internal::makeAttribute("parDim", obj.parDim() , data) );
        mp_node->append_node(tmp);

        appendBoxTopology(obj, mp_node, id_map, data);

        if (obj.numBoxProperties()!=0)
            gsWarn<<"Multi-patch object has box properties that are not written to XML\n";

        return mp_node;
    }

};

template <class T>
class gsXml< gsMultiBasis<T> >
{
private:
    gsXml() { }
    typedef gsMultiBasis<T> Object;

public:
    GSXML_COMMON_FUNCTIONS(Object);
    GSXML_GET_POINTER(Object);
    static std::string tag() { return "MultiBasis"; }
    static std::string type() { return ""; }

    static void get_into(gsXmlNode* node, Object & result)
    {
        GISMO_ASSERT( !strcmp( node->name(), "MultiBasis" ),
                      "Something went wrong. Expected MultiBasis tag." );

        gsXmlNode* topLevel = node->parent();

        const int d = atoi( node->first_attribute("parDim")->value() );

        gsXmlNode* patchNode = node->first_node("patches");
        std::istringstream iss;
        iss.str( patchNode->value() );

        typename gsMultiBasis<T>::BasisContainer bases;
        std::map<int, int> ids;
        if ( !strcmp( patchNode->first_attribute("type")->value(), "id_range") )
        {
            int first, last;
            gsGetInt(iss, first);
            gsGetInt(iss, last);
            for (int i = first; i <= last; ++i)
            {
                bases.push_back( getById< gsBasis<T> >( topLevel, i ) );
                ids[i] = i - first;
            }
        }
        else if ( !strcmp( patchNode->first_attribute("type")->value(), "id_index") )
        {
            int c = 0;
            for ( int pindex; gsGetInt(iss, pindex); )
            {
                bases.push_back( getById< gsBasis<T> >( topLevel, pindex ) );
                ids[pindex] = c++;
            }
        }
        else
        {
            gsWarn << "unknown tag in XML multipatch object \n";
        }

        // Read boundary
        std::vector< patchSide > boundaries;
        gsXmlNode * tmp = node->first_node("boundary");
        if (tmp)
            getBoundaries(tmp, ids, boundaries);

        // Read interfaces
        std::vector< boundaryInterface > interfaces;
        tmp = node->first_node("interfaces");
        if (tmp)
            getInterfaces(tmp, d, ids, interfaces);

        gsBoxTopology topology( d, bases.size(), boundaries, interfaces);

        result = gsMultiBasis<T>(bases, topology);
        freeAll(bases);
    }

    static gsXmlNode* put(const gsMultiBasis<T>& obj,
                          gsXmlTree& data)
    {
        // Insert all the basis
        std::map<index_t, index_t> id_map;
        int max_id = data.maxId();
        for ( typename gsMultiBasis<T>::const_iterator it = obj.begin();
              it != obj.end(); ++it )
        {
            gsXmlNode* basisXml = gsXml< gsBasis<T> >::put(**it, data);
            data.appendToRoot( basisXml );
            id_map[obj.findBasisIndex(*it) ] =  std::stoi(basisXml->first_attribute("id")->value());
        }

        std::ostringstream oss;
        oss<<  max_id+1 <<" "<< data.maxId();
        gsXmlNode* node = internal::makeNode("patches", oss.str(), data);
        node->append_attribute( internal::makeAttribute("type", "id_range", data) );
        oss.clear();
        oss.str("");

        gsXmlNode* mbNode = internal::makeNode(tag(), data);
        mbNode->append_attribute( internal::makeAttribute("parDim", obj.dim(), data) );
        mbNode->append_node(node);

        appendBoxTopology(obj.topology(), mbNode, id_map, data);

        return mbNode;
    }

};

} // namespace internal
} // namespace gismo
