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
#include <gsUtils/gsSortedVector.h>

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

        gsXmlNode * tmp = node->first_node("patches");
        GISMO_ASSERT(tmp, "No pathes tag");

        std::istringstream str;
        str.str( tmp->value() );
        // Resolve ID numbers
        std::map<int,int> ids;
        if ( ! strcmp( tmp->first_attribute("type")->value(), "id_range") )
        {
            int first, last;
            gsGetInt(str, first);
            gsGetInt(str, last);
            for ( int i = first; i<=last; ++i )
                ids[i] = i - first;
        }
        else if ( ! strcmp( tmp->first_attribute("type")->value(),"id_index") )
        {
            int c = 0;
            for (int pindex; gsGetInt(str, pindex);)
                ids[pindex] = c++;
        }
        else
        {
            gsWarn<<"Incomplete tag \"patch\" in boundaryConditions.\n";
        }

        // Read function inventory
        int count = countByTag("Function", node);
        std::vector<typename gsFunctionExpr<T>::Ptr> func(count);// todo: gsFunction::Ptr
        for (gsXmlNode * child = node->first_node("Function"); 
             child; child = child->next_sibling("Function") )
        {
            const int i = atoi( child->first_attribute("index")->value() );
            func[i]     = memory::make_shared(new gsFunctionExpr<T>);
            getFunctionFromXml(child, *func[i]);
        }

        // Read boundary conditions
        std::vector< patchSide > boundaries;
        for (gsXmlNode * child = node->first_node("bc"); 
             child; child = child->next_sibling("bc") )
        {
            const int uIndex = atoi( child->first_attribute("unknown")->value() );
            const int fIndex = atoi( child->first_attribute("function")->value() );
            //const int cIndex = atoi( child->first_attribute("comp")->value() );

            getBoundaries(child, ids, boundaries);
            
            const gsXmlAttribute * bcat = child->first_attribute("type");
            GISMO_ASSERT( NULL != bcat, "No type provided");
            const char * bctype = bcat->value();
            for (std::vector<patchSide>::const_iterator it = boundaries.begin();
                 it != boundaries.end(); ++it )
                result.add(it->patch, it->side(),
                           bctype, func[fIndex],
                           uIndex, false );//parametric
        }
        
        T val(0);
        for (gsXmlNode * child = node->first_node("cv"); 
             child; child = child->next_sibling("cv") )
        {   
            str.clear();
            str.str( child->value() );
            GISMO_ENSURE( gsGetReal(str, val), "No value");
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
        //gsXmlNode * mp = data.getRoot()->last_node("MultiPatch");

        gsWarn<<"To do\n";
            
        gsXmlNode * BCs = internal::makeNode("boundaryConditions" , data);
        data.appendToRoot(BCs);

        // inventory of functions
        typedef typename Object::const_bciterator bctype_it;
        typedef typename Object::const_iterator   bc_it;
        typedef typename gsSortedVector<typename gsFunction<T>::Ptr>::const_iterator fun_it;
        gsSortedVector<typename gsFunction<T>::Ptr> fun;
        typedef typename std::map<int,std::vector<const boundary_condition<T>*> >
            ::const_iterator inv_it;
        typedef typename std::vector<const boundary_condition<T>*>::const_iterator bcptr_it;

        std::map<int,std::vector<const boundary_condition<T>*> > fi;
        for (bctype_it it = obj.beginAll(); it!=obj.endAll(); ++it)
            for (bc_it bc = it->second.begin(); bc != it->second.end(); ++bc)
            {
                fun.push_sorted_unique(bc->function());
                fi[fun.getIndex(bc->function()) ].push_back(&(*bc));
            }
        
        int c = 0;        
        for (fun_it fit = fun.begin(); fit!=fun.end(); ++fit)
        {
            //gsXmlNode * ff = putFunctionFromXml(BCs, **fit);
            //ff index = c
            ++c;
        }
        
        // for all bcs, append bc, cv
        std::ostringstream oss;
        c = 0;
        for (inv_it it = fi.begin(); it!=fi.end(); ++it)
        {
            bcptr_it bc = it->second.begin();
            const std::string & label = (*bc)->ctype();
            for (; bc!=it->second.end(); ++it)
            {
            
            //unknown=, type, function=c
            //ps
//            oss << it->patch << " " << int(it->side()) << "\n";
            }
            //node->append_node(internal::makeNode("boundary", oss.str(), data));
        }
                    
        return BCs;
    }
};

} // end namespace internal

} // end namespace gismo
