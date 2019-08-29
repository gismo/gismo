/** @file gsXml.cpp

    @brief Provides implementation of input/output XML utilities struct.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <fstream>
#include <iomanip>      // std::setprecision

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsBoxTopology.h>

#include <gsIO/gsXmlUtils.h>

#include <rapidxml/rapidxml.hpp> // External file
#include <rapidxml/rapidxml_print.hpp> // External file

namespace gismo {

namespace internal {


/* Helpers to allocate XML data  */
    
char * makeValue( const std::string & value, gsXmlTree & data)
{
    return data.allocate_string( value.c_str() );
}

gsXmlAttribute * makeAttribute( const std::string & name, const std::string & value, gsXmlTree & data)
{
    return data.allocate_attribute( 
	data.allocate_string(name.c_str() ), 
	data.allocate_string(value.c_str()) );
}   

gsXmlAttribute * makeAttribute( const std::string & name, const unsigned & value, gsXmlTree & data)
{
    char tmp[16];
    sprintf (tmp,"%d",value);
    return data.allocate_attribute( 
	   data.allocate_string(name.c_str() ), 
	   data.allocate_string(tmp) );
}
    

gsXmlNode * makeNode( const std::string & name, gsXmlTree & data)
{
    return data.allocate_node(rapidxml::node_element ,
	   data.allocate_string(name.c_str()) );
}
    

gsXmlNode * makeNode( const std::string & name, const std::string & value, gsXmlTree & data)
{
    return data.allocate_node(rapidxml::node_element ,
	   data.allocate_string(name.c_str() ), 
	   data.allocate_string(value.c_str()) );
}

gsXmlNode * makeComment(const std::string & comment, gsXmlTree & data)
{
    return data.allocate_node(rapidxml::node_comment, 0, 
                       data.allocate_string(comment.c_str()));
}
    
std::string to_string(const unsigned & i)
{
    char tmp[4];
    sprintf (tmp,"%d",i);
    return tmp;
}

int countByTag(const std::string & tag, 
                      gsXmlNode * root )
{
    if ( ! root )
    {
        gsWarn<< "Invalid root node.\n";
        return 0;
    }

    int c = 0;
    for (gsXmlNode * child = root->first_node( tag.c_str() ) ; 
         child; child = child->next_sibling( tag.c_str() ) )
        c++;
    return c;
}

int countByTagType(const std::string & tag, 
                   const std::string & type, 
                   gsXmlNode * root )
{
    if ( ! root )
    {
        gsWarn<< "Invalid root node.\n";
        return 0;
    }
    if ( type == "" )
        return countByTag( tag.c_str(), root );
    else
    {
        int c = 0;
        for (gsXmlNode * child = root->first_node( tag.c_str() ) ; 
             child; child = child->next_sibling( tag.c_str() ) )
            if ( !strcmp( child->first_attribute("type")->value(), type.c_str() ) )
                c++;
        return c;
    }
}

gsXmlNode *  firstByTag(const std::string & tag, 
                        gsXmlNode * root )
{
    if ( ! root )
    {
        gsWarn<< "Invalid root node.\n";
        return NULL;
    }
    else
    {
        return root->first_node( tag.c_str() );
    }
}

gsXmlNode *  firstByTagType(const std::string & tag, 
                            const std::string & type, 
                            gsXmlNode * root )
{
    if ( ! root )
    {
        gsWarn<< "Invalid root node.\n";
        return NULL;
    }
    if ( type == "" )
        return root->first_node( tag.c_str() );
    else
    {
        for (gsXmlNode * child = root->first_node( tag.c_str() ) ; 
             child;   child = child->next_sibling( tag.c_str() ) )
            if ( !strcmp( child->first_attribute("type")->value(), 
                          type.c_str() ) )
                return child;
        return NULL;
    }   
}

gsXmlNode *  anyByTag(const std::string & tag, 
                      gsXmlNode * root )
{
    if ( ! root )
    {
        gsWarn<< "Invalid root node.\n";
        return NULL;
    }
    // Searching up to third level of the XML tree
    for (gsXmlNode * child = root->first_node() ; 
         child; child = child->next_sibling() )
    {
        if (!strcmp( child->name(), tag.c_str() ) ) 
            return child;
        // Level 2
        for (gsXmlNode * child2 = child->first_node() ; 
             child2; child2 = child2->next_sibling() )
        {
            if ( !strcmp( child2->name(), tag.c_str() ) )
                return child2;
            // Level 3
            for (gsXmlNode * child3 = child2->first_node() ; 
                 child3; child3 = child3->next_sibling() )
                if ( !strcmp( child3->name(), tag.c_str() ) )
                    return child3;
        }
    }
    return NULL;
}

/// Appends a box topology into node, used for gsMultiPatch and gsMultiBasis.
void appendBoxTopology(const gsBoxTopology& topology,
                       gsXmlNode* node,
                       gsXmlTree& data)
{
    std::ostringstream oss;

    if ( topology.nInterfaces() != 0 )
    {
        for ( gsBoxTopology::const_iiterator it = topology.iBegin();
              it != topology.iEnd(); ++it )
        {
            oss << it->first().patch  << " " << int(it->first().side()) << " "
                << it->second().patch << " " << int(it->second().side()) << " "
                << it->dirMap().transpose() << " "
                << it->dirOrientation().transpose() << "\n";
        }
        node->append_node(internal::makeNode("interfaces", oss.str(), data));
        oss.clear();
        oss.str("");
    }

    if ( topology.nBoundary() != 0)
    {
        for ( gsBoxTopology::const_biterator it = topology.bBegin();
              it != topology.bEnd(); ++it )
        {
            oss << it->patch << " " << int(it->side()) << "\n";
        }
        node->append_node(internal::makeNode("boundary", oss.str(), data));
        oss.clear();
        oss.str("");
    }
}

void getInterfaces(gsXmlNode* node,
                   const int d,
                   std::map<int, int>& ids,
                   std::vector< boundaryInterface > & result)
{
    GISMO_ASSERT( node != NULL, "Node does not exist");
    result.clear();
    
    // temporaries for interface reading
    gsVector<index_t> dirMap(d);
    gsVector<bool>    dirOrient(d);

    std::istringstream iss;
    iss.str( node->value() );
    
    gsVector<short_t> p(4); // { patch, side, patch, side }
    
    // Read interface (groups or size 4 + 2*d)
    while ( iss >> std::ws >> p[0] ) // While there are more ints (groups or size 4+d-1)
    {
        for ( int i = 1; i != 4; ++i )
        {
            if ( !(iss >> std::ws >> p[i] >> std::ws) )
            {
                gsWarn << "Error reading interface.\n";
            }
        }
        
        // Get ids
        p[0] = ids[ p[0] ];
        p[2] = ids[ p[2] ];
        
        // Read the matching direction permutation
        for ( int i = 0; i != d; ++i )
        {
            if ( !(iss >> std::ws >> dirMap[i]) )
            {
                gsWarn << "Error reading interface direction map.\n";
            }
        }
        
        // Read the interface orientation
        for ( int i = 0; i != d; ++i )
        {
            if ( !(iss >> std::ws >> dirOrient[i]) )
            {
                gsWarn << "Error reading interface orientation.\n";
            }
        }
        
        result.push_back( boundaryInterface(p, dirMap, dirOrient) );
        
//            // OLD format: read in Orientation flags
//            gsVector<bool> orient(d-1);// orientation flags
//            int k;
//            for ( int i=0; i!=d-1; ++i)
//            {
//                if ( !(str >> std::ws >> k >> std::ws) )
//                    gsWarn<<"Error reading interface orientation.\n";
//                orient[i]= (k>0);
//            }
//            interfaces.push_back( boundaryInterface(p,orient) ) ;
    }
}

void getBoundaries(gsXmlNode * node, std::map<int, int> & ids,
        std::vector<patchSide> & result)
{
    GISMO_ASSERT(node != NULL, "Node does not exist");
    result.clear();
    
    std::istringstream iss;
    iss.str(node->value());
    int patch, side;
    
    while (iss >> std::ws >> patch)
    {
        if (ids.size() != 0)
        {
            patch = ids[patch];
        }
        iss >> std::ws >> side;
        result.push_back(patchSide(patch, side));
    }
}

}// end namespace internal

}// end namespace gismo
