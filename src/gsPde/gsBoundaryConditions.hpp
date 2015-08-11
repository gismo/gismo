/** @file gsBoundaryConditions.hpp

    @brief Implementation file for the gsBoundaryConditions class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsIO/gsXml.h>
#include <gsCore/gsFunctionExpr.h>

namespace gismo
{

namespace internal
{

/// @brief I/O for boundary conditions from file
template<class T>
class gsXml< gsBoundaryConditions<T> >
{
private:
    gsXml() { }
    typedef gsBoundaryConditions<T> Object;
public:
    GSXML_COMMON_FUNCTIONS(Object);
    static std::string tag  () { return "boundaryConditions"; }
    static std::string type () { return ""; }

    GSXML_GET_POINTER(Object);

    static void get_into (gsXmlNode * node, Object & result)
    {
        GISMO_ASSERT( !strcmp( node->name(), tag().c_str() ),
                      "Something went wrong. Expected tag "<< tag() );
        
        // There should exist a sibling of type MultiPatch
        gsXmlNode * toplevel = node->parent();
        gsXmlNode * mp       = NULL;

        const gsXmlAttribute * mp_at = node->first_attribute("multipatch");
        if ( mp_at )
        {
            const int d = atoi( mp_at->value() );
            mp = searchId(d, toplevel);
        }
        else
        {
            // multipatch not referenced, grab the first one in the
            // file
            mp = toplevel->first_node("MultiPatch");
        }

        if ( mp == NULL || strcmp( mp->name(), "MultiPatch" ) )
            gsWarn <<"Did not find a mulitpatch object.\n";

        gsXmlNode * tmp = mp->first_node("patches");
        GISMO_ASSERT(tmp, "No pathes tag");

        std::istringstream str;
        str.str( tmp->value() );

        // Resolve ID numbers
        std::map<int,int> ids;
        if ( ! strcmp( tmp->first_attribute("type")->value(), "id_range") )
        {
            int first, last;
            str >> std::ws >>  first >>  std::ws >> last >> std::ws ;
            for ( int i = first; i<=last; ++i )
            {
                GISMO_ASSERT( searchId(i, toplevel) != NULL, 
                              "Invalid reference to node Id");
                ids[i] = i - first;
            }
        }
        else if ( ! strcmp( mp->first_attribute("type")->value(),"id_index") )
        {
            int c = 0;
            for (int pindex; str >> pindex;)
            {
                GISMO_ASSERT( searchId(pindex, toplevel) != NULL, 
                              "Invalid reference to node Id");
                ids[pindex] = c++;
            }
        }
        else
        {
            gsWarn<<"Unknown tag in XML multipatch object.\n";
        }

        // Read boundary
        gsXmlNode * boundaryNode = mp->first_node("boundary");
        std::vector< patchSide > boundaries;
        if (boundaryNode)
        {
            getBoundaries(boundaryNode, ids, boundaries);
        }

        // Read function data
        int count = countByTag("Function", node);
        std::vector<typename gsFunctionExpr<T>::Ptr> func(count);// todo: gsFunction::Ptr
        for (gsXmlNode * child = node->first_node("Function"); 
             child; child = child->next_sibling("Function") )
        {
            const int i = atoi( child->first_attribute("index")->value() );
            func[i]     = shared(new gsFunctionExpr<T>);
            getFunctionFromXml(child, *func[i]);            
        }

        // Read boundary conditions
        for (gsXmlNode * child = node->first_node("bc"); 
             child; child = child->next_sibling("bc") )
        {
            const int uIndex = atoi( child->first_attribute("unknown")->value() );
            const int fIndex = atoi( child->first_attribute("function")->value() );
        
            str.clear(); //str.str("");
            str.str( child->value() );
            if ( !strcmp(child->first_attribute("type")->value(), "dirichlet") )
            {
                for (int bIndex; str >> bIndex;) 
                    result.addCondition( boundaries[bIndex], 
                                         condition_type::dirichlet, 
                                         func[fIndex], uIndex ); //,parametric
            }
            else if ( !strcmp(child->first_attribute("type")->value(), "neumann") )
            {		       
                for (int bIndex; str >> bIndex;) 
                    result.addCondition( boundaries[bIndex],
                                         condition_type::neumann, 
                                         func[fIndex], uIndex ); //,parametric
            }
        }
        
        T val;
        for (gsXmlNode * child = node->first_node("cv"); 
             child; child = child->next_sibling("cv") )
        {   
            str.clear();
            str.str( child->value() );       
            GISMO_ENSURE( str >>  val, "No value");
            const int uIndex = atoi( child->first_attribute("unknown")->value() );
            const int cIndex = atoi( child->first_attribute("corner")->value() );
            int pIndex       = atoi( child->first_attribute("patch")->value() );
            pIndex = ids[pIndex];
            
            result.addCornerValue(cIndex, val, pIndex, uIndex);                
        }
    }
    
    static gsXmlNode * put (const Object & obj, 
                            gsXmlTree & data )
    {
        // Check if the last node is a multipatch
        //gsXmlNode * mp = deta.getRoot()->last_node("MultiPatch");

        gsWarn<<"To do\n";

        gsXmlNode * BCs = internal::makeNode("boundaryConditions" , data);
        data.appendToRoot(BCs);

        // collect function pointers

        // for all bcs, append bc, cv

        return BCs;
    }
};

} // end namespace internal

} // end namespace gismo
