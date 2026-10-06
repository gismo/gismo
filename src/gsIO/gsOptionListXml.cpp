#include <gsIO/gsOptionListXml.h>

namespace gismo
{

namespace internal
{

GISMO_EXPORT void gsXml<gsOptionList>::get_into(gsXmlNode * node, gsOptionList & result)
{
    // get all child-nodes
    gsXmlNode * tmp = node->first_node();
    while ( tmp )
    {
        const char* name = tmp->name();

        const std::string label = tmp->first_attribute("label")->value();
        const std::string desc = tmp->first_attribute("desc")->value();
        const std::string val = tmp->first_attribute("value")->value();

        if (strcmp("int", name) == 0)
        {
            std::istringstream str;
            str.str( val );
            index_t myVal;
            gsGetInt(str, myVal);
            result.addInt(label, desc, myVal);
        }
        else if (strcmp("real", name) == 0)
        {
            std::istringstream str;
            str.str( val );
            real_t myVal;
            gsGetReal(str, myVal);
            result.addReal(label, desc, myVal);
        }
        else if (strcmp("bool", name) == 0)
        {
            std::istringstream str;
            str.str( val );
            index_t myVal;
            gsGetInt(str, myVal);
            result.addSwitch(label, desc, (0 != myVal) );
        }
        else
        {
            result.addString(label, desc, val);
        }
        tmp =  tmp->next_sibling();
    }
}

GISMO_EXPORT gsXmlNode *
gsXml<gsOptionList>::put (const gsOptionList & obj, gsXmlTree & data)
{
    // Append data
    gsXmlNode * optionList = internal::makeNode("OptionList", data);

    // /*
    // iterate over all strings
    std::vector<gsOptionList::OptionListEntry> entries = obj.getAllEntries();
    std::vector<gsOptionList::OptionListEntry>::const_iterator it;
    for (it = entries.begin(); it != entries.end(); it++)
    {
        const gsOptionList::OptionListEntry & entry = *it;
        gsXmlNode * node_str = internal::makeNode(entry.type, data);
        gsXmlAttribute * attr_label = internal::makeAttribute("label", entry.label, data);
        gsXmlAttribute * attr_desc = internal::makeAttribute("desc", entry.desc, data);
        gsXmlAttribute * attr_val = internal::makeAttribute("value", entry.val, data);
        node_str->insert_attribute(0, attr_label);
        node_str->insert_attribute(0, attr_desc);
        node_str->insert_attribute(0, attr_val);
        optionList->insert_node(0, node_str);
    }
    // */

    /*
      gsXmlNode * tmp;
      gsXmlAttribute * atr;
      gsOptionList::StringTable::const_iterator it1;
      for ( it1 = obj.m_strings.begin(); it1 != obj.m_strings.end(); it1++ )
      {
      tmp = internal::makeNode("string", data);
      atr = internal::makeAttribute("label", it1->first, data);
      tmp->insert_attribute(0, atr);
      atr = internal::makeAttribute("desc", it1->second.second, data);
      tmp->insert_attribute(0, atr);
      atr = internal::makeAttribute("value", it1->second.first, data);
      tmp->insert_attribute(0, atr);
      optionList->insert_node(0, tmp);
      }
      gsOptionList::IntTable::const_iterator it2;
      for ( it2 = obj.m_ints.begin(); it2 != obj.m_ints.end(); it2++ )
      {
      tmp = internal::makeNode("int", data);
      atr = internal::makeAttribute("label", it2->first, data);
      tmp->insert_attribute(0, atr);
      atr = internal::makeAttribute("desc", it2->second.second, data);
      tmp->insert_attribute(0, atr);
      atr = internal::makeAttribute("value", it2->second.first, data);
      tmp->insert_attribute(0, atr);
      optionList->insert_node(0, tmp);
      }
      gsOptionList::RealTable::const_iterator it3;
      for ( it3 = obj.m_reals.begin(); it3 != obj.m_reals.end(); it3++ )
      {
      tmp = internal::makeNode("real", data);
      atr = internal::makeAttribute("label", it3->first, data);
      tmp->insert_attribute(0, atr);
      atr = internal::makeAttribute("desc", it3->second.second, data);
      tmp->insert_attribute(0, atr);
      atr = internal::makeAttribute("value", it3->second.first, data);
      tmp->insert_attribute(0, atr);
      optionList->insert_node(0, tmp);
      }
      gsOptionList::SwitchTable::const_iterator it4;
      for ( it4 = obj.m_switches.begin(); it4 != obj.m_switches.end(); it4++ )
      {
      tmp = internal::makeNode("switch", data);
      atr = internal::makeAttribute("label", it4->first, data);
      tmp->insert_attribute(0, atr);
      atr = internal::makeAttribute("desc", it4->second.second, data);
      tmp->insert_attribute(0, atr);
      atr = internal::makeAttribute("value", it4->second.first, data);
      tmp->insert_attribute(0, atr);
      optionList->insert_node(0, tmp);
      }
    */

    return optionList;
}

} // namespace internal

}
