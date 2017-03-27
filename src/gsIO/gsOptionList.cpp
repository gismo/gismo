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
    GISMO_ASSERT(it!=m_strings.end(), "Invalid request (getString): "<<label<<" is not given or not a string.");
    return it->second.first;
}

int gsOptionList::getInt(const std::string & label) const
{
    IntTable::const_iterator it = m_ints.find(label);
    GISMO_ASSERT(it!=m_ints.end(), "Invalid request (getInt): "<<label<<" is not given or not an int.");
    return it->second.first;
}

bool gsOptionList::getSwitch(const std::string & label) const
{
    SwitchTable::const_iterator it = m_switches.find(label);
    GISMO_ASSERT(it!=m_switches.end(), "Invalid request (getSwitch): "<<label<<" is not given or not a switch.");
    return it->second.first;
}

real_t gsOptionList::getReal(const std::string & label) const
{
    RealTable::const_iterator it = m_reals.find(label);
    GISMO_ASSERT(it!=m_reals.end(), "Invalid request (getReal): "<<label<<" is not given or not a real.");
    return it->second.first;
}

std::string gsOptionList::askString(const std::string & label,
                                    const std::string & val) const
{
    StringTable::const_iterator it = m_strings.find(label);
#if defined(GISMO_EXTRA_DEBUG)
    if ( it == m_strings.end() && exists(label) )
        gsWarn << "Invalid request (askString): "<<label<<" is given, but not a string.\n";
#endif
    return ( it == m_strings.end() ? val : it->second.first);
}

int gsOptionList::askInt(const std::string & label, const int & val) const
{
    IntTable::const_iterator it = m_ints.find(label);
#if defined(GISMO_EXTRA_DEBUG)
    if ( it == m_ints.end() && exists(label) )
        gsWarn << "Invalid request (askInt): "<<label<<" is given, but not an int.\n";
#endif
    return ( it == m_ints.end() ? val : it->second.first);
}

bool gsOptionList::askSwitch(const std::string & label, const bool & res) const
{
    SwitchTable::const_iterator it = m_switches.find(label);
#if defined(GISMO_EXTRA_DEBUG)
    if ( it == m_switches.end() && exists(label) )
        gsWarn << "Invalid request (askSwitch): "<<label<<" is given, but not a switch.\n";
#endif
    return ( it == m_switches.end() ? res : it->second.first);
}

real_t gsOptionList::askReal(const std::string & label, const real_t & val) const
{
    RealTable::const_iterator it = m_reals.find(label);
#if defined(GISMO_EXTRA_DEBUG)
    if ( it == m_reals.end() && exists(label) )
        gsWarn << "Invalid request (askReal): "<<label<<" is given, but not a real.\n";
#endif
    return ( it == m_reals.end() ? val : it->second.first);
}

void gsOptionList::setString(const std::string & label,
                             const std::string & value)
{
    StringTable::iterator it = m_strings.find(label);
    if ( it == m_strings.end() )
    {
        if ( exists(label) )
        {
            printInfo(label);
            GISMO_ERROR("Invalid request (setString): "<<label<<" is not a string.");
        }
        GISMO_ERROR("Invalid request (setString): "<<label<<" does not exist.");
    }
    it->second.first = value;
}

void gsOptionList::setInt(const std::string & label,
                          const int & value)
{
    IntTable::iterator it = m_ints.find(label);
    if ( it == m_ints.end() )
    {
        if ( exists(label) )
        {
            printInfo(label);
            GISMO_ERROR("Invalid request (setInt): "<<label<<" is not an int.");
        }
        GISMO_ERROR("Invalid request (setInt): "<<label<<" does not exist.");
    }
    it->second.first = value;
}

void gsOptionList::setReal(const std::string & label,
                           const real_t & value)
{
    RealTable::iterator it = m_reals.find(label);
    if ( it == m_reals.end() )
    {
        if ( exists(label) )
        {
            printInfo(label);
            GISMO_ERROR("Invalid request (setReal): "<<label<<" is not a real.");
        }
        GISMO_ERROR("Invalid request (setReal): "<<label<<" does not exist.");
    }
    it->second.first = value;
}

void gsOptionList::setSwitch(const std::string & label,
                             const bool & value)
{
    SwitchTable::iterator it = m_switches.find(label);
    if ( it == m_switches.end() )
    {
        if ( exists(label) )
        {
            printInfo(label);
            GISMO_ERROR("Invalid request (setSwitch): "<<label<<" is not a switch.");
        }
        GISMO_ERROR("Invalid request (setSwitch): "<<label<<" does not exist.");
    }
    it->second.first = value;
}

void gsOptionList::addString(const std::string & label,
                             const std::string & desc,
                             const std::string & value)
{
    GISMO_ASSERT( !( isInt(label) || isReal(label) || isSwitch(label) ),
        "Invalid request (addString): Option "<<label<<" already exists, but not as a string." );
    //GISMO_ASSERT( !exists(label), "Option "<<label<<" already exists." );
    m_strings[label] = std::make_pair(value,desc);
}


void gsOptionList::addInt(const std::string & label,
                          const std::string & desc,
                          const int & value)
{
    GISMO_ASSERT( !( isString(label) || isReal(label) || isSwitch(label) ),
        "Invalid request (addInt): Option "<<label<<" already exists, but not as an int." );
    //GISMO_ASSERT( !exists(label), "Option "<<label<<" already exists." );
    m_ints[label] = std::make_pair(value,desc);
}

void gsOptionList::addReal(const std::string & label,
                           const std::string & desc,
                           const real_t & value)
{
    GISMO_ASSERT( !( isString(label) || isInt(label) || isSwitch(label) ),
         "Invalid request (addReal): Option "<<label<<" already exists, but not as a real." );
    //GISMO_ASSERT( !exists(label), "Option "<<label<<" already exists." );
    m_reals[label] = std::make_pair(value,desc);
}

void gsOptionList::addSwitch(const std::string & label,
                             const std::string & desc,
                             const bool & value)
{
    GISMO_ASSERT( !( isString(label) || isInt(label) || isReal(label) ),
         "Invalid request (addSwitch): Option "<<label<<" already exists, but not as a switch." );
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
    gsOptionList::StringTable::const_iterator it1;
    for ( it1 = other.m_strings.begin(); it1 != other.m_strings.end(); it1++ )
    {
        if (exists(it1->first))                        setString(it1->first,it1->second.first);
        else if(type == gsOptionList::addIfUnknown)    addString(it1->first,it1->second.second,it1->second.first);
    }

    // add integers to list
    gsOptionList::IntTable::const_iterator it2;
    for ( it2 = other.m_ints.begin(); it2 != other.m_ints.end(); it2++ )
    {
        if (exists(it2->first))                       setInt(it2->first,it2->second.first);
        else if(type == gsOptionList::addIfUnknown)   addInt(it2->first,it2->second.second,it2->second.first);
    }

    // add reals to list
    gsOptionList::RealTable::const_iterator it3;
    for ( it3 = other.m_reals.begin(); it3 != other.m_reals.end(); it3++ )
    {
        if (exists(it3->first))                       setReal(it3->first,it3->second.first);
        else if(type == gsOptionList::addIfUnknown)   addReal(it3->first,it3->second.second,it3->second.first);
    }

    // add bools to list
    gsOptionList::SwitchTable::const_iterator it4;
    for ( it4 = other.m_switches.begin(); it4 != other.m_switches.end(); it4++ )
    {
        if (exists(it4->first))                       setSwitch(it4->first,it4->second.first);
        else if(type == gsOptionList::addIfUnknown)   addSwitch(it4->first,it4->second.second,it4->second.first);
    }
}

gsOptionList gsOptionList::wrapIntoGroup(const std::string & gn) const
{
    GISMO_ASSERT( !hasGroup(gn), "Invalid request (wrapIntoGroup): A group cannot be wrapped into a group with the same name." );

    const std::string prepend = gn+".";
    gsOptionList result;

    // add strings to list
    gsOptionList::StringTable::const_iterator it1;
    for ( it1 = m_strings.begin(); it1 != m_strings.end(); it1++ )
        result.addString(prepend+it1->first,it1->second.second,it1->second.first);

    // add integers to list
    gsOptionList::IntTable::const_iterator it2;
    for ( it2 = m_ints.begin(); it2 != m_ints.end(); it2++ )
        result.addInt(prepend+it2->first,it2->second.second,it2->second.first);

    // add reals to list
    gsOptionList::RealTable::const_iterator it3;
    for ( it3 = m_reals.begin(); it3 != m_reals.end(); it3++ )
        result.addReal(prepend+it3->first,it3->second.second,it3->second.first);

    // add bools to list
    gsOptionList::SwitchTable::const_iterator it4;
    for ( it4 = m_switches.begin(); it4 != m_switches.end(); it4++ )
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
    gsOptionList::StringTable::const_iterator it1;
    for ( it1 = m_strings.begin(); it1 != m_strings.end(); it1++ )
        if ( util::starts_with(it1->first,search) )
            result.addString(it1->first.substr(len),it1->second.second,it1->second.first);

    // add integers to list
    gsOptionList::IntTable::const_iterator it2;
    for ( it2 = m_ints.begin(); it2 != m_ints.end(); it2++ )
        if ( util::starts_with(it2->first,search) )
            result.addInt(it2->first.substr(len),it2->second.second,it2->second.first);

    // add reals to list
    gsOptionList::RealTable::const_iterator it3;
    for ( it3 = m_reals.begin(); it3 != m_reals.end(); it3++ )
        if ( util::starts_with(it3->first,search) )
            result.addReal(it3->first.substr(len),it3->second.second,it3->second.first);

    // add bools to list
    gsOptionList::SwitchTable::const_iterator it4;
    for ( it4 = m_switches.begin(); it4 != m_switches.end(); it4++ )
        if ( util::starts_with(it4->first,search) )
            result.addSwitch(it4->first.substr(len),it4->second.second,it4->second.first);

    return result;
}

bool gsOptionList::hasGlobals() const
{
    gsOptionList::StringTable::const_iterator it1;
    for ( it1 = m_strings.begin(); it1 != m_strings.end(); it1++ )
        if ( std::string::npos == it1->first.find('.') ) return true;

    gsOptionList::IntTable::const_iterator it2;
    for ( it2 = m_ints.begin(); it2 != m_ints.end(); it2++ )
        if ( std::string::npos == it2->first.find('.') ) return true;

    gsOptionList::RealTable::const_iterator it3;
    for ( it3 = m_reals.begin(); it3 != m_reals.end(); it3++ )
        if ( std::string::npos == it3->first.find('.') ) return true;

    gsOptionList::SwitchTable::const_iterator it4;
    for ( it4 = m_switches.begin(); it4 != m_switches.end(); it4++ )
        if ( std::string::npos == it4->first.find('.') ) return true;

    return false;
}

bool gsOptionList::hasGroup(const std::string & gn) const
{
    const std::string search = gn+".";

    gsOptionList::StringTable::const_iterator it1;
    for ( it1 = m_strings.begin(); it1 != m_strings.end(); it1++ )
        if ( util::starts_with(it1->first,search) ) return true;

    gsOptionList::IntTable::const_iterator it2;
    for ( it2 = m_ints.begin(); it2 != m_ints.end(); it2++ )
        if ( util::starts_with(it2->first,search) ) return true;

    gsOptionList::RealTable::const_iterator it3;
    for ( it3 = m_reals.begin(); it3 != m_reals.end(); it3++ )
        if ( util::starts_with(it3->first,search) ) return true;

    gsOptionList::SwitchTable::const_iterator it4;
    for ( it4 = m_switches.begin(); it4 != m_switches.end(); it4++ )
        if ( util::starts_with(it4->first,search) ) return true;

    return false;
}

// /*
// todo: can we implement the loops in a more intelligent way?
// e.g. re-factor this function by moving the loops into its own function
std::vector<gsOptionList::OptionListEntry> gsOptionList::getAllEntries() const
{
    std::vector<gsOptionList::OptionListEntry> result;
    result.reserve(4*size());
    const char * XML_STR = "string";
    const char * XML_INT = "int";
    const char * XML_REAL = "real";
    const char * XML_BOOL = "bool";
    // add strings to list
    gsOptionList::StringTable::const_iterator it1;
    for ( it1 = m_strings.begin(); it1 != m_strings.end(); it1++ )
    {
        gsOptionList::OptionListEntry entry;
        entry.type = XML_STR;
        entry.label = it1->first;
        std::stringstream str;
        str.str( it1->second.first );
        entry.val = str.str();
        entry.desc = it1->second.second;
        result.push_back(entry);
    }
    // add integers to list
    gsOptionList::IntTable::const_iterator it2;
    for ( it2 = m_ints.begin(); it2 != m_ints.end(); it2++ )
    {
        gsOptionList::OptionListEntry entry;
        entry.type = XML_INT;
        entry.label = it2->first;
        std::stringstream str;
        str << it2->second.first;
        entry.val = str.str();
        entry.desc = it2->second.second;
        result.push_back(entry);
    }
    // add reals to list
    gsOptionList::RealTable::const_iterator it3;
    for ( it3 = m_reals.begin(); it3 != m_reals.end(); it3++ )
    {
        gsOptionList::OptionListEntry entry;
        entry.type = XML_REAL;
        entry.label = it3->first;
        std::stringstream str;
        str << it3->second.first;
        entry.val = str.str();
        entry.desc = it3->second.second;
        result.push_back(entry);
    }
    // add bools to list
    gsOptionList::SwitchTable::const_iterator it4;
    for ( it4 = m_switches.begin(); it4 != m_switches.end(); it4++ )
    {
        gsOptionList::OptionListEntry entry;
        entry.type = XML_BOOL;
        entry.label = it4->first;
        std::stringstream str;
        str << it4->second.first;
        entry.val = str.str();
        entry.desc = it4->second.second;
        result.push_back(entry);
    }
    return result;
}
//*/


#define OL_PRINT_INFO(it,type)                                          \
    os<<"* "<<std::setw(17)<<std::left<<it->first <<std::setw(12)<<std::right<<" ("#type") = " \
    <<std::setw(7)<<std::left<<it->second.first<<" "<<it->second.second<<"\n"
//<<std::boolalpha

std::ostream & gsOptionList::print(std::ostream & os) const
{
    os<<"Options ("<<size()<<"):\n";
    for (StringTable::const_iterator it1 = m_strings.begin();it1!=m_strings.end();++it1)
        OL_PRINT_INFO(it1,string);
    for (IntTable::const_iterator it2 = m_ints.begin();it2!=m_ints.end();++it2)
        OL_PRINT_INFO(it2,int);
    for (RealTable::const_iterator it3 = m_reals.begin();it3!=m_reals.end();++it3)
        OL_PRINT_INFO(it3,real);
    for (SwitchTable::const_iterator it4 = m_switches.begin();it4!=m_switches.end();++it4)
        OL_PRINT_INFO(it4,switch);
    return os;
}

void gsOptionList::printInfo(const std::string& label) const
{
    std::ostream & os = gsInfo;
    StringTable::const_iterator it1 = m_strings.find(label);
    if ( it1 != m_strings.end() )
    {
        OL_PRINT_INFO(it1,string);
        return;
    }
    IntTable::const_iterator it2 = m_ints.find(label);
    if ( it2 != m_ints.end() )
    {
        OL_PRINT_INFO(it2,int);
        return;
    }
    RealTable::const_iterator it3 = m_reals.find(label);
    if ( it3 != m_reals.end() )
    {
        OL_PRINT_INFO(it3,real);
        return;
    }
    SwitchTable::const_iterator it4 = m_switches.find(label);
    if ( it4 != m_switches.end() )
    {
        OL_PRINT_INFO(it4,switch);
        os<<"* "<<it4->second.second<<" (switch) \n  "<<it4->first <<" = "<<(it4->second.first ? "ON" : "OFF")<<"\n";
        return;
    }
    gsInfo <<"Problem: "<< label <<" does not exist.\n";
}

#undef OL_PRINT_INFO

bool gsOptionList::exists(const std::string & label) const
{
    StringTable::const_iterator it1 = m_strings.find(label);
    if ( it1 != m_strings.end() )   return true;
    IntTable::const_iterator it2 = m_ints.find(label);
    if ( it2 != m_ints.end() )      return true;
    RealTable::const_iterator it3 = m_reals.find(label);
    if ( it3 != m_reals.end() )     return true;
    SwitchTable::const_iterator it4 = m_switches.find(label);
    if ( it4 != m_switches.end() )  return true;
    return false;
}

bool gsOptionList::isString(const std::string & label) const
{
    StringTable::const_iterator it = m_strings.find(label);
    return it != m_strings.end();
}

bool gsOptionList::isInt(const std::string & label) const
{
    IntTable::const_iterator it = m_ints.find(label);
    return it != m_ints.end();
}

bool gsOptionList::isReal(const std::string & label) const
{
    RealTable::const_iterator it = m_reals.find(label);
    return it != m_reals.end();
}

bool gsOptionList::isSwitch(const std::string & label) const
{
    SwitchTable::const_iterator it = m_switches.find(label);
    return it != m_switches.end();
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
            int myVal;
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
            int myVal;
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
