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
#include <gsCore/gsConstantFunction.h>
#include <gsUtils/gsSortedVector.h>
#include <gsIO/gsXmlGenericUtils.hpp>

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

    static void get_into(gsXmlNode * node, Object & result)
    {
        GISMO_ASSERT(!strcmp(node->name(), tag().c_str()),
                "Something went wrong. Expected tag "<< tag());

        std::istringstream str;
        std::map<int, int> ids;

        // Check if any of the BCs is defined on a boundary set name
        const int mp_index = atoi(node->first_attribute("multipatch")->value());
        gsXmlNode* toplevel = node->parent();
        std::vector< patchSide > allboundaries;
        for (gsXmlNode * child = node->first_node("bc"); child;
                child = child->next_sibling("bc"))
        {
            std::map<int, int> tmp_ids;;

            const gsXmlAttribute * att_name = child->first_attribute("name");
            if (NULL != att_name)
            {
                gsXmlNode* mp_node = searchId(mp_index, toplevel);
                GISMO_ASSERT( mp_node != NULL,
                              "No Multipatch with Id "<<mp_index<<" found in the XML data.");

                gsXmlNode * tmp = mp_node->first_node("patches");
                std::istringstream tmp_str ;
                tmp_str.str( tmp->value() );
                // Handle id_range or id_index for multipatches. This is needed to assign the right indices for the BCs
                if ( ! strcmp( tmp->first_attribute("type")->value(),"id_range") )
                {
                    int first, last;
                    gsGetInt(tmp_str, first);
                    gsGetInt(tmp_str, last);
                    for ( int i = first; i<=last; ++i )
                        tmp_ids[i] = i - first;
                }
                else if ( ! strcmp( tmp->first_attribute("type")->value(),"id_index") )
                {
                    int c = 0;
                    for (int pindex; gsGetInt(tmp_str, pindex);)
                        tmp_ids[pindex] = c++;
                }
                else
                {
                    gsWarn<<"Unknown tag in XML multipatch object.\n";
                }

                for (gsXmlNode * child = mp_node->first_node("boundary"); child;
                        child = child->next_sibling("boundary"))
                {
                    std::vector< patchSide > tmp_boundaries;
                    if (child)
                    {
                        getBoundaries(child, tmp_ids, tmp_boundaries);
                        allboundaries.insert( allboundaries.end(), tmp_boundaries.begin(), tmp_boundaries.end() );
                    }
                }
                break;
            }
        }

        //gsXmlNode * tmp = node->first_node("patches");
        //GISMO_ASSERT(tmp, "No pathes tag");
//
//        std::istringstream str;
//        str.str( tmp->value() );
//        // Resolve ID numbers
//        std::map<int,int> ids;
//        if ( ! strcmp( tmp->first_attribute("type")->value(), "id_range") )
//        {
//            int first, last;
//            gsGetInt(str, first);
//            gsGetInt(str, last);
//            for ( int i = first; i<=last; ++i )
//                ids[i] = i - first;
//        }
//        else if ( ! strcmp( tmp->first_attribute("type")->value(),"id_index") )
//        {
//            int c = 0;
//            for (int pindex; gsGetInt(str, pindex);)
//                ids[pindex] = c++;
//        }
//        else
//        {
//            gsWarn<<"Incomplete tag \"patch\" in boundaryConditions.\n";
//        }

        // Read function inventory
        int count = countByTag("Function", node);
        std::vector<typename gsFunctionExpr<T>::Ptr> func(count); // todo: gsFunction::Ptr
        for (gsXmlNode * child = node->first_node("Function"); child; child =
                child->next_sibling("Function"))
        {
            const int i = atoi(child->first_attribute("index")->value());
            func[i] = memory::make_shared(new gsFunctionExpr<T>);
            internal::gsXml<gsFunctionExpr<T> >::get_into(child, *func[i]);
        }

        // Read boundary conditions
        std::vector<patchSide> boundaries;
        for (gsXmlNode * child = node->first_node("bc"); child;
                child = child->next_sibling("bc"))
        {
            const int uIndex = atoi(child->first_attribute("unknown")->value());
            const int fIndex = atoi(
                    child->first_attribute("function")->value());

            const gsXmlAttribute * comp = child->first_attribute("component");
            int cIndex = -1;
            if (NULL != comp)
                cIndex = atoi( comp->value() );

            const gsXmlAttribute * att_ispar = child->first_attribute("parametric");
            bool ispar = false;
            if (NULL != att_ispar)
                ispar = atoi( att_ispar->value() );

            const gsXmlAttribute * att_name = child->first_attribute("name");
            if (NULL != att_name)
            {
                boundaries.clear();
                std::string name = att_name->value();
                for (typename std::vector<patchSide>::const_iterator it=allboundaries.begin(); it!=allboundaries.end(); it++)
                    if (it->label()==name)
                        boundaries.push_back(*it);
            }
            else
                getBoundaries(child, ids, boundaries);

            if (boundaries.size() == 0) {
              gsWarn << "Boundary condition without boundary to apply to. The"
                        " following bc will be unused\n" << *child
                     << std::endl;
            }

            const gsXmlAttribute * bcat = child->first_attribute("type");
            GISMO_ASSERT(NULL != bcat, "No type provided");
            const char * bctype = bcat->value();
            for (std::vector<patchSide>::const_iterator it = boundaries.begin();
                    it != boundaries.end(); ++it)
                result.add(it->patch, it->side(), bctype,
                           func[fIndex], uIndex,cIndex, ispar);
        }

        T val(0);
        for (gsXmlNode * child = node->first_node("cv"); child;
                child = child->next_sibling("cv"))
        {
            str.clear();
            str.str(child->value());
            GISMO_ENSURE(gsGetReal(str, val), "No value");

            // Unknown is optional, otherwise 0
            const gsXmlAttribute * unk = child->first_attribute("unknown");
            int uIndex = 0;
            if (NULL != unk)
                uIndex = atoi( unk->value() );

            // Component is optional, otherwise -1
            const gsXmlAttribute * comp = child->first_attribute("component");
            int cIndex = -1;
            if (NULL != comp)
                cIndex = atoi( comp->value() );

            const int cornIndex = atoi(child->first_attribute("corner")->value());
            int pIndex = atoi(child->first_attribute("patch")->value());

            result.addCornerValue(cornIndex, val, pIndex, uIndex, cIndex);
        }
    }

    static gsXmlNode * put(const Object & obj, gsXmlTree & data)
    {
        // Check if the last node is a multipatch
        //gsXmlNode * mp = data.getRoot()->last_node("MultiPatch");

        gsXmlNode * BCs = internal::makeNode("boundaryConditions", data);
        //data.appendToRoot(BCs);
        gsXmlAttribute * multi = internal::makeAttribute("multipatch", "0",
                data);

        BCs->append_attribute(multi);

        // inventory of functions
        typedef typename gsBoundaryConditions<T>::const_bciterator bctype_it;
        typedef typename Object::const_iterator bc_it;

        std::vector<typename gsFunctionSet<T>::Ptr> fun;
        //gsSortedVector<typename gsFunction<T>::Ptr> fun;
        typedef typename std::vector<const boundary_condition<T>*> bctype_vec;
        typedef typename std::map<int, bctype_vec> bctype_map;
        std::map<std::string, bctype_map> fi;

        for (bctype_it it = obj.beginAll(); it != obj.endAll(); ++it)
        {
            std::string label = it->first;
            bctype_map map;
            for (bc_it bc = it->second.begin(); bc != it->second.end(); ++bc)
            {
                typename gsFunctionSet<T>::Ptr ptr = bc->function();
                bool contains = std::find(fun.begin(), fun.end(), ptr)
                        != fun.end();
                if (!contains)
                {
                    fun.push_back(ptr);
                }
                int index = std::find(fun.begin(), fun.end(), ptr)
                        - fun.begin();
//                fun.push_sorted_unique(ptr);
//                int index = fun.getIndex(ptr);
                std::vector<const boundary_condition<T>*> vec = map[index];
                const boundary_condition<T>* b = &(*bc);
                vec.push_back(b);
                map[index] = vec;
            }
            std::pair<std::string, bctype_map> pair(label, map);
            fi.insert(pair);
        }

        int count = 0;
        typedef typename std::vector<typename gsFunctionSet<T>::Ptr>::const_iterator fun_it;
        for (fun_it fit = fun.begin(); fit != fun.end(); ++fit)
        {
            gsXmlNode * ff = putFunctionToXml<T>(*fit, data, count);
            BCs->append_node(ff);
            ++count;
        }

        // for all bcs, append bc, cv
        typedef typename std::map<std::string, bctype_map>::const_iterator bctype_map_it;
        typedef typename std::map<int, bctype_vec>::const_iterator bctype_iv_it;
        typedef typename bctype_vec::const_iterator bctype_vec_it;

        count = 0;
        for (bctype_map_it it = fi.begin(); it != fi.end(); ++it)
        {
            std::string label = it->first;
            //gsDebug << "Label='" << label << "'\n";
            bctype_map map = it->second;

            for (bctype_iv_it bcV = map.begin(); bcV != map.end(); ++bcV)
            {
                int index = bcV->first;
                bctype_vec vec = bcV->second;
                //gsDebug << "index='" << index << "'\n";
                //gsDebug << "vec='" << vec.size() << "'\n";
                gsXmlNode * bcNode = internal::makeNode("bc", data);
                gsXmlAttribute * typeNode = internal::makeAttribute("type",
                        label, data);
                gsXmlAttribute * indexNode = internal::makeAttribute("function",
                        index, data);
                bcNode->append_attribute(typeNode);
                bcNode->append_attribute(indexNode);
                bool first = true;
                std::ostringstream oss;
                for (bctype_vec_it bc = vec.begin(); bc != vec.end(); ++bc)
                {
                    const boundary_condition<T> b = (**bc);
                    //gsDebug << "iterate over boundary condition with '"
                    //        << b.m_label << "'\n";
                    if (first)
                    {
                        gsXmlAttribute * unknownNode = internal::makeAttribute(
                                "unknown", b.m_unknown, data);
                        gsXmlAttribute * componentNode = internal::makeAttribute("component",
                                b.unkComponent(), data);
                        bcNode->append_attribute(componentNode);
                        bcNode->append_attribute(unknownNode);
                        first = false;
                    }
                    else
                    {
                        oss << " ";
                    }
                    oss << b.ps.patch << " " << b.ps.m_index;
                }
                char * value = data.allocate_string(oss.str().c_str());
                bcNode->value(value);
                BCs->append_node(bcNode);
                ++count;
            }
        }
        typename gsBoundaryConditions<T>::const_citerator ci;
        for (ci = obj.cornerValues().begin(); ci != obj.cornerValues().end();
                ci++)
        {
            corner_value<T> c = *ci;
            gsXmlNode * cvNode = internal::makeNode("cv", data);
            gsXmlAttribute * unknownNode = internal::makeAttribute("unknown",
                    c.unknown, data);
            gsXmlAttribute * componentNode = internal::makeAttribute("component",
                    c.component, data);
            gsXmlAttribute * patchNode = internal::makeAttribute("patch",
                    c.patch, data);
            gsXmlAttribute * cornerNode = internal::makeAttribute("corner",
                    c.corner.m_index, data);
            cvNode->append_attribute(unknownNode);
            cvNode->append_attribute(componentNode);
            cvNode->append_attribute(patchNode);
            cvNode->append_attribute(cornerNode);
            std::ostringstream oss;
            oss << c.value;
            char * value = data.allocate_string(oss.str().c_str());
            cvNode->value(value);
            /// gsDebug << "Corner value='" << c.value << ", " << c.patch << ", " << c.unknown << "'\n";
            BCs->append_node(cvNode);
        }
        return BCs;
    }
};

} // end namespace internal

} // end namespace gismo
