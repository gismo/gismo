/** @file gsOptionList.cpp

    @brief Provides a list of labeled parameters/options that can be
    set and accessed easily

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H. Weiner, S. Takacs
*/


#include <cstring>
#include <string>
#include <sstream>

#include <gsIO/gsOptionList.h>

#include <gsIO/gsXml.h>

namespace gismo
{

std::string gsOptionList::getString(const std::string & label) const
{
    StringTable::const_iterator it = m_strings.find(label);
    GISMO_ENSURE(it!=m_strings.end(), "Invalid request (getString): "<<label<<" is not a string; it is "<<getInfo(label)<<".");
    return it->second.first;
}

index_t gsOptionList::getInt(const std::string & label) const
{
    IntTable::const_iterator it = m_ints.find(label);
    GISMO_ENSURE(it!=m_ints.end(), "Invalid request (getInt): "<<label<<" is not not an int; it is "<<getInfo(label)<<".");
    return it->second.first;
}

real_t gsOptionList::getReal(const std::string & label) const
{
    RealTable::const_iterator it = m_reals.find(label);
    GISMO_ENSURE(it!=m_reals.end(), "Invalid request (getReal): "<<label<<" is not a real; it is "<<getInfo(label)<<".");
    return it->second.first;
}

bool gsOptionList::getSwitch(const std::string & label) const
{
    SwitchTable::const_iterator it = m_switches.find(label);
    GISMO_ENSURE(it!=m_switches.end(), "Invalid request (getSwitch): "<<label<<" is not not a switch; it is "<<getInfo(label)<<".");
    return it->second.first;
}

std::vector<std::string> gsOptionList::getMultiString(const std::string & gn) const
{
    GISMO_ASSERT(hasGroup(gn), "Invalid request (getMultiString): The group " + gn + " does not exist.");

    std::vector<std::string> result;
    const std::string search = gn + ".";
    size_t sz = static_cast<size_t>(getInt(search + "Size"));
    result.reserve(sz);
    // add strings to vector
    for (size_t i = 0; i < sz; ++i)
        result.push_back(getString(search + util::to_string(i)));

    return result;
}

std::vector<index_t> gsOptionList::getMultiInt(const std::string & gn) const
{
    GISMO_ASSERT(hasGroup(gn), "Invalid request (getMultiInt): The group " + gn + " does not exist.");

    std::vector<index_t> result;

    const std::string search = gn + ".";

    // add integers to vector
    for (IntTable::const_iterator it = m_ints.begin(); it != m_ints.end(); it++)
        if (util::starts_with(it->first, search) && !util::ends_with(it->first, "Size"))
            result.push_back(it->second.first);

    return result;
}

std::vector<real_t> gsOptionList::getMultiReal(const std::string & gn) const
{
    GISMO_ASSERT(hasGroup(gn), "Invalid request (getMultiReal): The group " + gn + " does not exist.");

    std::vector<real_t> result;

    const std::string search = gn + ".";

    // add reals to vector
    for (RealTable::const_iterator it = m_reals.begin(); it != m_reals.end(); it++)
        if (util::starts_with(it->first, search) && !util::ends_with(it->first, "Size"))
            result.push_back(it->second.first);

    return result;
}

std::string gsOptionList::askString(const std::string & label,
                                    const std::string & value) const
{
    StringTable::const_iterator it = m_strings.find(label);
#if defined(GISMO_EXTRA_DEBUG)
    if ( it == m_strings.end() && exists(label) )
        gsWarn << "Invalid request (askString): "<<label<<" is given, but not a string; it is "<<getInfo(label)<<".\n";
#endif
    return ( it == m_strings.end() ? value : it->second.first);
}

index_t gsOptionList::askInt(const std::string & label,
                         const index_t & value) const
{
    IntTable::const_iterator it = m_ints.find(label);
#if defined(GISMO_EXTRA_DEBUG)
    if ( it == m_ints.end() && exists(label) )
        gsWarn << "Invalid request (askInt): "<<label<<" is given, but not an int; it is "<<getInfo(label)<<".\n";
#endif
    return ( it == m_ints.end() ? value : it->second.first);
}

bool gsOptionList::askSwitch(const std::string & label,
                             const bool & value) const
{
    SwitchTable::const_iterator it = m_switches.find(label);
#if defined(GISMO_EXTRA_DEBUG)
    if ( it == m_switches.end() && exists(label) )
        gsWarn << "Invalid request (askSwitch): "<<label<<" is given, but not a switch; it is "<<getInfo(label)<<".\n";
#endif
    return ( it == m_switches.end() ? value : it->second.first);
}

real_t gsOptionList::askReal(const std::string & label,
                             const real_t & value) const
{
    RealTable::const_iterator it = m_reals.find(label);
#if defined(GISMO_EXTRA_DEBUG)
    if ( it == m_reals.end() && exists(label) )
        gsWarn << "Invalid request (askReal): "<<label<<" is given, but not a real; it is "<<getInfo(label)<<".\n";
#endif
    return ( it == m_reals.end() ? value : it->second.first);
}

void gsOptionList::setString(const std::string & label,
                             const std::string & value)
{
    StringTable::iterator it = m_strings.find(label);
    GISMO_ENSURE(it!=m_strings.end(), "Invalid request (setString): "<<label<<" is not a string; it is "<<getInfo(label)<<".");
    it->second.first = value;
}

void gsOptionList::setInt(const std::string & label,
                          const index_t & value)
{
    IntTable::iterator it = m_ints.find(label);
    GISMO_ENSURE(it!=m_ints.end(), "Invalid request (setInt): "<<label<<" is not a int; it is "<<getInfo(label)<<".");
    it->second.first = value;
}

void gsOptionList::setReal(const std::string & label,
                           const real_t & value)
{
    RealTable::iterator it = m_reals.find(label);
    GISMO_ENSURE(it!=m_reals.end(), "Invalid request (setReal): "<<label<<" is not a real; it is "<<getInfo(label)<<".");
    it->second.first = value;
}

void gsOptionList::setSwitch(const std::string & label,
                             const bool & value)
{
    SwitchTable::iterator it = m_switches.find(label);
    GISMO_ENSURE(it!=m_switches.end(), "Invalid request (setSwitch): "<<label<<" is not a switch; it is "<<getInfo(label)<<".");
    it->second.first = value;
}

void gsOptionList::addString(const std::string & label,
                             const std::string & desc,
                             const std::string & value)
{
    GISMO_ENSURE( !( isInt(label) || isReal(label) || isSwitch(label) ),
        "Invalid request (addString): Option "<<label<<" already exists, but not as a string; it is "<<getInfo(label)<<"." );
    //GISMO_ASSERT( !exists(label), "Option "<<label<<" already exists." );
    m_strings[label] = std::make_pair(value,desc);
}


void gsOptionList::addInt(const std::string & label,
                          const std::string & desc,
                          const index_t & value)
{
    GISMO_ENSURE( !( isString(label) || isReal(label) || isSwitch(label) ),
        "Invalid request (addInt): Option "<<label<<" already exists, but not as an int; it is "<<getInfo(label)<<"." );
    //GISMO_ASSERT( !exists(label), "Option "<<label<<" already exists." );
    m_ints[label] = std::make_pair(value,desc);
}

void gsOptionList::addReal(const std::string & label,
                           const std::string & desc,
                           const real_t & value)
{
    GISMO_ENSURE( !( isString(label) || isInt(label) || isSwitch(label) ),
         "Invalid request (addReal): Option "<<label<<" already exists, but not as a real; it is "<<getInfo(label)<<"." );
    //GISMO_ASSERT( !exists(label), "Option "<<label<<" already exists." );
    m_reals[label] = std::make_pair(value,desc);
}

void gsOptionList::addMultiInt(const std::string & label,
                               const std::string & desc,
                               const std::vector<index_t> & values)
{
    GISMO_ENSURE( !( isString(label) || isReal(label) || isSwitch(label) ),
                  "Invalid request (addMultiInt): Option "<<label<<" already exists, but not as an multiint; it is "<<getInfo(label)<<"." );

    for (size_t i = 0; i < values.size(); ++i)
    {
        addInt(label + "." + util::to_string(i), desc, values[i]);
    }
    addInt(label + ".Size", desc, values.size());
}

void gsOptionList::addSwitch(const std::string & label,
                             const std::string & desc,
                             const bool & value)
{
    GISMO_ENSURE( !( isString(label) || isInt(label) || isReal(label) ),
         "Invalid request (addSwitch): Option "<<label<<" already exists, but not as a switch; it is "<<getInfo(label)<<"." );
    //GISMO_ASSERT( !exists(label), "Option "<<label<<" already exists." );
    m_switches[label] = std::make_pair(value,desc);
}

void gsOptionList::remove(const std::string & label)
{
    if (m_strings .erase(label)) return;
    if (m_ints    .erase(label)) return;
    if (m_reals   .erase(label)) return;
    m_switches.erase(label);
}

void gsOptionList::update(const gsOptionList & other, gsOptionList::updateType type)
{
    // add strings to list
    for ( StringTable::const_iterator it1 = other.m_strings.begin(); it1 != other.m_strings.end(); it1++ )
    {
        if (exists(it1->first))          setString(it1->first,it1->second.first);
        else if(type == addIfUnknown)    addString(it1->first,it1->second.second,it1->second.first);
    }

    // add integers to list
    for ( IntTable::const_iterator it2 = other.m_ints.begin(); it2 != other.m_ints.end(); it2++ )
    {
        if (exists(it2->first))         setInt(it2->first,it2->second.first);
        else if(type == addIfUnknown)   addInt(it2->first,it2->second.second,it2->second.first);
    }

    // add reals to list
    for ( RealTable::const_iterator it3 = other.m_reals.begin(); it3 != other.m_reals.end(); it3++ )
    {
        if (exists(it3->first))         setReal(it3->first,it3->second.first);
        else if(type == addIfUnknown)   addReal(it3->first,it3->second.second,it3->second.first);
    }

    // add bools to list
    for ( SwitchTable::const_iterator it4 = other.m_switches.begin(); it4 != other.m_switches.end(); it4++ )
    {
        if (exists(it4->first))         setSwitch(it4->first,it4->second.first);
        else if(type == addIfUnknown)   addSwitch(it4->first,it4->second.second,it4->second.first);
    }
}

gsOptionList gsOptionList::wrapIntoGroup(const std::string & gn) const
{
    GISMO_ASSERT( !hasGroup(gn), "Invalid request (wrapIntoGroup): A group cannot be wrapped into a group with the same name." );

    const std::string prepend = gn+".";
    gsOptionList result;

    // add strings to list
    for ( StringTable::const_iterator it1 = m_strings.begin(); it1 != m_strings.end(); it1++ )
        result.addString(prepend+it1->first,it1->second.second,it1->second.first);

    // add integers to list
    for ( IntTable::const_iterator it2 = m_ints.begin(); it2 != m_ints.end(); it2++ )
        result.addInt(prepend+it2->first,it2->second.second,it2->second.first);

    // add reals to list
    for ( RealTable::const_iterator it3 = m_reals.begin(); it3 != m_reals.end(); it3++ )
        result.addReal(prepend+it3->first,it3->second.second,it3->second.first);

    // add bools to list
    for ( SwitchTable::const_iterator it4 = m_switches.begin(); it4 != m_switches.end(); it4++ )
        result.addSwitch(prepend+it4->first,it4->second.second,it4->second.first);

    return result;
}

gsOptionList gsOptionList::getGroup(const std::string & gn) const
{
    GISMO_ASSERT( hasGroup(gn), "Invalid request (getGroup): The group "+gn+" does not exist." );

    gsOptionList result;

    const std::string search = gn+".";
    const index_t len = search.length();

    // add strings to list
    for ( StringTable::const_iterator it1 = m_strings.begin(); it1 != m_strings.end(); it1++ )
        if ( util::starts_with(it1->first,search) )
            result.addString(it1->first.substr(len),it1->second.second,it1->second.first);

    // add integers to list
    for ( IntTable::const_iterator it2 = m_ints.begin(); it2 != m_ints.end(); it2++ )
        if ( util::starts_with(it2->first,search) )
            result.addInt(it2->first.substr(len),it2->second.second,it2->second.first);

    // add reals to list
    for ( RealTable::const_iterator it3 = m_reals.begin(); it3 != m_reals.end(); it3++ )
        if ( util::starts_with(it3->first,search) )
            result.addReal(it3->first.substr(len),it3->second.second,it3->second.first);

    // add bools to list
    for ( SwitchTable::const_iterator it4 = m_switches.begin(); it4 != m_switches.end(); it4++ )
        if ( util::starts_with(it4->first,search) )
            result.addSwitch(it4->first.substr(len),it4->second.second,it4->second.first);

    return result;
}

bool gsOptionList::hasGlobals() const
{
    // check strings
    for ( StringTable::const_iterator it1 = m_strings.begin(); it1 != m_strings.end(); it1++ )
        if ( std::string::npos == it1->first.find('.') ) return true;

    // check integers
    for ( IntTable::const_iterator it2 = m_ints.begin(); it2 != m_ints.end(); it2++ )
        if ( std::string::npos == it2->first.find('.') ) return true;

    // check reals
    for ( RealTable::const_iterator it3 = m_reals.begin(); it3 != m_reals.end(); it3++ )
        if ( std::string::npos == it3->first.find('.') ) return true;

    // check bools
    for ( SwitchTable::const_iterator it4 = m_switches.begin(); it4 != m_switches.end(); it4++ )
        if ( std::string::npos == it4->first.find('.') ) return true;

    return false;
}

bool gsOptionList::hasGroup(const std::string & gn) const
{
    const std::string search = gn+".";

    // check strings
    for ( StringTable::const_iterator it1 = m_strings.begin(); it1 != m_strings.end(); it1++ )
        if ( util::starts_with(it1->first,search) ) return true;

    // check integers
    for ( IntTable::const_iterator it2 = m_ints.begin(); it2 != m_ints.end(); it2++ )
        if ( util::starts_with(it2->first,search) ) return true;

    // check reals
    for ( RealTable::const_iterator it3 = m_reals.begin(); it3 != m_reals.end(); it3++ )
        if ( util::starts_with(it3->first,search) ) return true;

    // check bools
    for ( SwitchTable::const_iterator it4 = m_switches.begin(); it4 != m_switches.end(); it4++ )
        if ( util::starts_with(it4->first,search) ) return true;

    return false;
}

namespace internal
{
template <typename It>
inline gsOptionList::OptionListEntry makeOptionListEntry(const char * type, It it)
{
    gsOptionList::OptionListEntry entry;
    entry.type = type;
    entry.label = it->first;
    entry.val = util::to_string(it->second.first);
    entry.desc = it->second.second;
    return entry;
}
} // namespace internal

std::vector<gsOptionList::OptionListEntry> gsOptionList::getAllEntries() const
{
    std::vector<OptionListEntry> result;
    result.reserve(size());
    const char * XML_STR = "string";
    const char * XML_INT = "int";
    const char * XML_REAL = "real";
    const char * XML_BOOL = "bool";

    // handle strings
    for ( StringTable::const_iterator it1 = m_strings.begin(); it1 != m_strings.end(); it1++ )
        result.push_back( internal::makeOptionListEntry(XML_STR, it1) );

    // handle integers
    for ( IntTable::const_iterator it2 = m_ints.begin(); it2 != m_ints.end(); it2++ )
        result.push_back( internal::makeOptionListEntry(XML_INT, it2) );

    // handle reals
    for ( RealTable::const_iterator it3 = m_reals.begin(); it3 != m_reals.end(); it3++ )
        result.push_back( internal::makeOptionListEntry(XML_REAL, it3) );

    // handle bools
    for ( SwitchTable::const_iterator it4 = m_switches.begin(); it4 != m_switches.end(); it4++ )
        result.push_back( internal::makeOptionListEntry(XML_BOOL, it4) );

    return result;
}

std::ostream & gsOptionList::OptionListEntry::print(std::ostream & os, index_t slot_label) const
{
    const index_t slot_val = 8;
    const index_t sz_label = label.size();
    const index_t sz_val = val.size();
    if (sz_label<=slot_label && sz_val<=slot_val)
        os <<"* "<<std::setw(slot_label)<<std::left<<label<<std::setw(8)<<std::right<<("("+type)<<") = "
            <<std::setw(slot_val)<<std::left<<val<<" "<<desc<<"\n";
    else
        os <<"* "<<std::setw(slot_label)<<std::left<<label<<std::setw(8)<<std::right<<("("+type)<<") = "
            <<val<<"\n"<<std::setw(slot_label+slot_val+15)<<" "<<desc<<"\n";
    return os;
}

std::ostream & gsOptionList::print(std::ostream & os) const
{
    typedef std::vector<OptionListEntry> DataTable;
    DataTable data = getAllEntries();
    os<<"Options ("<<data.size()<<"):\n";
    std::sort( data.begin(), data.end() );
    index_t slot_label = 15;
    for ( DataTable::const_iterator it = data.begin(); it != data.end(); it++ )
        slot_label = std::max( slot_label, (index_t)it->label.size() );
    slot_label = std::min( slot_label, (index_t)35 );
    for ( DataTable::const_iterator it = data.begin(); it != data.end(); it++ )
        it->print(os, slot_label);
    return os;
}

std::string gsOptionList::getInfo(const std::string& label) const
{
    // find in strings
    StringTable::const_iterator it1 = m_strings.find(label);
    if ( it1 != m_strings.end() )
        return "a string (value:\"" + it1->second.first + "\")";

    // find in integers
    IntTable::const_iterator it2 = m_ints.find(label);
    if ( it2 != m_ints.end() )
        return "an int (value:" + util::to_string(it2->second.first) + ")";

    // find in reals
    RealTable::const_iterator it3 = m_reals.find(label);
    if ( it3 != m_reals.end() )
        return "a real (value:" + util::to_string(it3->second.first) + ")";

    // find in bools
    SwitchTable::const_iterator it4 = m_switches.find(label);
    if ( it4 != m_switches.end() )
        return "a switch (value:" + util::to_string(it4->second.first) + ")";

    return "undefined";
}

bool gsOptionList::exists(const std::string & label) const
{
    if ( m_strings.find(label)  != m_strings.end()  )  return true;
    if ( m_ints.find(label)     != m_ints.end()     )  return true;
    if ( m_reals.find(label)    != m_reals.end()    )  return true;
    if ( m_switches.find(label) != m_switches.end() )  return true;
    return false;
}

bool gsOptionList::isString(const std::string & label) const
{
    return m_strings.find(label) != m_strings.end();
}

bool gsOptionList::isInt(const std::string & label) const
{
    return m_ints.find(label) != m_ints.end();
}

bool gsOptionList::isReal(const std::string & label) const
{
    return m_reals.find(label) != m_reals.end();
}

bool gsOptionList::isSwitch(const std::string & label) const
{
    return m_switches.find(label) != m_switches.end();
}



namespace internal
{

void gsXml<gsOptionList>::get_into(gsXmlNode * node, gsOptionList & result)
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

gsXmlNode *
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

} //namespace gismo
