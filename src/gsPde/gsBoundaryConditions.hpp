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

    static void get_into(gsXmlNode * node, Object & result)
    {
        GISMO_ASSERT(!strcmp(node->name(), tag().c_str()),
                "Something went wrong. Expected tag "<< tag());

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
        std::istringstream str;
        std::map<int, int> ids;

        // Read function inventory
        int count = countByTag("Function", node);
        std::vector<typename gsFunctionExpr<T>::Ptr> func(count); // todo: gsFunction::Ptr
        for (gsXmlNode * child = node->first_node("Function"); child; child =
                child->next_sibling("Function"))
        {
            const int i = atoi(child->first_attribute("index")->value());
            func[i] = memory::make_shared(new gsFunctionExpr<T>);
            getFunctionFromXml(child, *func[i]);
        }

        // Read boundary conditions
        std::vector<patchSide> boundaries;
        for (gsXmlNode * child = node->first_node("bc"); child;
                child = child->next_sibling("bc"))
        {
            const int uIndex = atoi(child->first_attribute("unknown")->value());
            const int fIndex = atoi(
                    child->first_attribute("function")->value());
            //const int cIndex = atoi( child->first_attribute("comp")->value() );

            getBoundaries(child, ids, boundaries);
            
            const gsXmlAttribute * bcat = child->first_attribute("type");
            GISMO_ASSERT(NULL != bcat, "No type provided");
            const char * bctype = bcat->value();
            for (std::vector<patchSide>::const_iterator it = boundaries.begin();
                    it != boundaries.end(); ++it)
                result.add(it->patch, it->side(), bctype, func[fIndex], uIndex,
                        false);        //parametric
        }
        
        T val(0);
        for (gsXmlNode * child = node->first_node("cv"); child;
                child = child->next_sibling("cv"))
        {
            str.clear();
            str.str(child->value());
            GISMO_ENSURE(gsGetReal(str, val), "No value");
            const int uIndex = atoi(child->first_attribute("unknown")->value());
            const int cIndex = atoi(child->first_attribute("corner")->value());
            int pIndex = atoi(child->first_attribute("patch")->value());
            pIndex = ids[pIndex];
            
            result.addCornerValue(cIndex, val, pIndex, uIndex);
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

        std::vector<typename gsFunction<T>::Ptr> fun;
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
                typename gsFunction<T>::Ptr ptr = bc->function();
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
        typedef typename std::vector<typename gsFunction<T>::Ptr>::const_iterator fun_it;
        for (fun_it fit = fun.begin(); fit != fun.end(); ++fit)
        {
            gsXmlNode * ff = putFunctionFromXml(*fit, data, count);
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
            gsXmlAttribute * patchNode = internal::makeAttribute("patch",
                    c.patch, data);
            gsXmlAttribute * cornerNode = internal::makeAttribute("corner",
                    c.corner.m_index, data);
            cvNode->append_attribute(unknownNode);
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

private:
    static gsXmlNode * putFunctionFromXml(
            const typename gsFunction<T>::Ptr & obj, gsXmlTree & data,
            int index)
    {
        gsXmlNode * result = internal::makeNode("Function", data);
        if (typeid(*obj) == typeid(gsFunctionExpr<T> ))
        {
            gsFunctionExpr<T> * ptr2 =
                    dynamic_cast<gsFunctionExpr<T> *>(obj.get());
            gsFunctionExpr<T> expr = *ptr2;
            result = putFunctionExprFromXml(expr, result, data);
        }
        gsXmlAttribute * indexNode = internal::makeAttribute("index", index,
                data);
        result->append_attribute(indexNode);
        return result;
    }

    static gsXmlNode * putFunctionExprFromXml(const gsFunctionExpr<T> & obj,
            gsXmlNode * result, gsXmlTree & data)
    {
        std::string typeStr = gsXml<gsFunctionExpr<T> >::type();
        gsXmlAttribute * type = internal::makeAttribute("type", typeStr, data);
        result->append_attribute(type);
        gsXmlAttribute * dim = internal::makeAttribute("dim", obj.domainDim(),
                data);
        result->append_attribute(dim);
        // set value
        std::ostringstream stream;
        bool first = true;
        for (short_t i = 0; i < obj.targetDim(); ++i)
        {
            if (!first)
            {
                stream << ", ";
            }
            stream << obj.expression(i);
            first = false;
        }
        //val = stream.str();
        char * value = data.allocate_string(stream.str().c_str());
        result->value(value);
        return result;
    }

};

} // end namespace internal

} // end namespace gismo
