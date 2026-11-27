/** @file gsWriteOpenNurbs.cpp

    @brief Writes the objects of a G+Smo XML tree to the Rhinoceros 3DM file format.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Angelos Mantzaflaris
*/

#include <gsOpennurbs/gsWriteOpenNurbs.h>
#include <gsOpennurbs/gsWriteOpenNurbs.hpp>   // writeON_* (... ONX_Model&, ...) templates

#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsUtils/gsMesh/gsMesh.h>

#include <onurbs/opennurbs.h>

#include <sstream>
#include <string>

namespace gismo
{

namespace extensions
{

namespace
{

/// Name of the 3dm object made from \a node: its label, else tag+id, else tag+count
std::string objectName(internal::gsXmlNode * node, const char * tag, index_t count)
{
    std::ostringstream nm;
    if (const rapidxml::xml_attribute<char> * lab = node->first_attribute("label"))
        nm << lab->value();
    else if (const rapidxml::xml_attribute<char> * id = node->first_attribute("id"))
        nm << tag << id->value();
    else
        nm << tag << count;
    return nm.str();
}

} // anonymous namespace

bool gsWriteOpenNurbs(const char * fn, const internal::gsXmlTree & data)
{
    internal::gsXmlNode * root = data.first_node("xml");
    if (!root)
    {
        gsWarn << "gsWriteOpenNurbs: no <xml> root node found.\n";
        return false;
    }

    ONX_Model model;
    writeON_Init(model);

    index_t count = 0;
    for (internal::gsXmlNode * child = root->first_node(); child; child = child->next_sibling())
    {
        const char * tag = child->name();
        const std::string tagStr(tag);
        const std::string nm = objectName(child, tag, count);
        ++count;

        if (tagStr == "Geometry")
        {
            memory::unique_ptr< gsGeometry<real_t> > g(internal::gsXml< gsGeometry<real_t> >::get(child));
            if (!g)
            {
                gsWarn << "gsWriteOpenNurbs: could not read <Geometry> " << nm << ".\n";
                continue;
            }

            const bool dimOk = (g->geoDim() == 2 || g->geoDim() == 3);
            if (const gsBSpline<real_t> * c = dynamic_cast<const gsBSpline<real_t>*>(g.get()))
            {
                if (dimOk)
                {
                    writeON_NurbsCurve(*c, model, nm);
                    continue;
                }
            }
            else if (const gsTensorBSpline<2,real_t> * s = dynamic_cast<const gsTensorBSpline<2,real_t>*>(g.get()))
            {
                if (dimOk)
                {
                    writeON_NurbsSurface(*s, model, nm);
                    continue;
                }
            }

            const rapidxml::xml_attribute<char> * type = child->first_attribute("type");
            gsWarn << "gsWriteOpenNurbs: skipping geometry of type "
                   << (type ? type->value() : "unknown") << ".\n";
        }
        else if (tagStr == "MultiPatch")
        {
            // The patches of a multipatch are separate root-level <Geometry> nodes
            // and are written by the branch above; topology has no 3dm counterpart.
        }
        else if (tagStr == "Mesh")
        {
            memory::unique_ptr< gsMesh<real_t> > m(internal::gsXml< gsMesh<real_t> >::get(child));
            if (!m)
            {
                gsWarn << "gsWriteOpenNurbs: could not read <Mesh> " << nm << ".\n";
                continue;
            }
            writeON_Mesh(*m, model, nm);
        }
        else
        {
            gsWarn << "gsWriteOpenNurbs: skipping <" << tag << ">\n";
        }
    }

    if (model.m_object_table.Count() == 0)
    {
        ON::End();
        gsWarn << "gsWriteOpenNurbs: nothing to write to " << fn << ".\n";
        return false;
    }

    return writeON_Write3dm(model, fn);
}

} // namespace extensions

} // namespace gismo
