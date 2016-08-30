/** @file gsOptionList.cpp

    @brief Provides a list of labeled parameters/options that can be
    set and accessed easily

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#include <string>

#include <gsIO/gsOptionList.h>

#include <gsIO/gsXml.h>

namespace gismo
{

std::string gsOptionList::getString(const std::string & label,
                                    const std::string & val) const
{
    StringTable::const_iterator it = m_strings.find(label);
    return ( it == m_strings.end() ? val : it->second.first);
}

int gsOptionList::getInt(const std::string & label, const int & val) const
{
    IntTable::const_iterator it = m_ints.find(label);
    return ( it == m_ints.end() ? val : it->second.first);
}

bool gsOptionList::getSwitch(const std::string & label, const bool & res) const
{
    SwitchTable::const_iterator it = m_switches.find(label);
    return ( it == m_switches.end() ? res : it->second.first);
}

real_t gsOptionList::getReal(const std::string & label, const real_t & val) const
{
    RealTable::const_iterator it = m_reals.find(label);
    return ( it == m_reals.end() ? val : it->second.first);
}

void gsOptionList::setString(const std::string& label,
                             const std::string & value)
{
    StringTable::iterator it = m_strings.find(label);
    if ( it == m_strings.end() )
    {
        if ( exists(label) )
        {
            printInfo(label);
            GISMO_ERROR("Invalid request: "<<label <<" is not a string.");
        }
        // m_strings[label] = std::make_pair(value,""); return;
        GISMO_ERROR("Invalid request: "<<label <<" does not exist.");
    }
    it->second.first = value;
}

void gsOptionList::setInt(const std::string& label,
                          const int & value)
{
    IntTable::iterator it = m_ints.find(label);
    if ( it == m_ints.end() )
    {
        if ( exists(label) )
        {
            printInfo(label);
            GISMO_ERROR("Invalid request: "<<label <<" is not a string.");
        }
        // m_strings[label] = std::make_pair(value,""); return;
        GISMO_ERROR("Invalid request: "<<label <<" does not exist.");
    }
    it->second.first = value;
}

void gsOptionList::setReal(const std::string& label,
                           const real_t & value)
{
    StringTable::iterator it = m_strings.find(label);
    if ( it == m_strings.end() )
    {
        if ( exists(label) )
        {
            printInfo(label);
            GISMO_ERROR("Invalid request: "<<label <<" is not a string.");
        }
        // m_strings[label] = std::make_pair(value,""); return;
        GISMO_ERROR("Invalid request: "<<label <<" does not exist.");
    }
    it->second.first = value;
}

void gsOptionList::setSwitch(const std::string& label,
                             const bool & value)
{
    SwitchTable::iterator it = m_switches.find(label);
    if ( it == m_switches.end() )
    {
        if ( exists(label) )
        {
            printInfo(label);
            GISMO_ERROR("Invalid request: "<<label <<" is not a switch.");
        }
        // m_strings[label] = std::make_pair(value,""); return;
        GISMO_ERROR("Invalid request: "<<label <<" does not exist.");
    }
    it->second.first = value;
}

void gsOptionList::addString(const std::string& label,
                             const std::string& desc,
                             const std::string & value)
{
    if ( exists(label) )
    {
        printInfo(label);
        GISMO_ERROR("Option "<<label<<" already exists.");        
    }
    m_strings[label] = std::make_pair(value,desc);
}


void gsOptionList::addInt(const std::string& label,
                          const std::string& desc,
                          const int & value)
{
    if ( exists(label) )
    {
        printInfo(label);
        GISMO_ERROR("Option "<<label<<" already exists.");        
    }
    m_ints[label] = std::make_pair(value,desc);
}

void gsOptionList::addReal(const std::string& label,
                          const std::string& desc,
                          const real_t & value)
{
    if ( exists(label) )
    {
        printInfo(label);
        GISMO_ERROR("Option "<<label<<" already exists.");        
    }
    m_reals[label] = std::make_pair(value,desc);
}

void gsOptionList::addSwitch(const std::string& label,
                             const std::string& desc,
                             const bool & value)
{
    if ( exists(label) )
    {
        printInfo(label);
        GISMO_ERROR("Option "<<label<<" already exists.");        
    }
    m_switches[label] = std::make_pair(value,desc);
}

bool gsOptionList::exists(const std::string& label)
{
    StringTable::iterator it1 = m_strings.find(label);
    if ( it1 != m_strings.end() )   return true;
    IntTable::iterator it2 = m_ints.find(label);
    if ( it2 != m_ints.end() )      return true;
    RealTable::iterator it3 = m_reals.find(label);
    if ( it3 != m_reals.end() )     return true;
    SwitchTable::iterator it4 = m_switches.find(label);
    if ( it4 != m_switches.end() )  return true;
    return false;
}

void gsOptionList::print(std::ostream & os)
{
    os<<"Options ("<<size()<<"):\n";
    for (StringTable::iterator it1 = m_strings.begin();it1!=m_strings.end();++it1)
        os<<"* "<<it1->second.second<<" (string) \n  "<<it1->first <<" = "<<it1->second.first<<"\n";
    for (IntTable::iterator it2 = m_ints.begin();it2!=m_ints.end();++it2)
        os<<"* "<<it2->second.second<<" (int) \n  "<<it2->first <<" = "<<it2->second.first<<"\n";
    for (RealTable::iterator it3 = m_reals.begin();it3!=m_reals.end();++it3)
        os<<"* "<<it3->second.second<<" (real) \n  "<<it3->first <<" = "<<it3->second.first<<"\n";
    for (SwitchTable::iterator it4 = m_switches.begin();it4!=m_switches.end();++it4)
        os<<"* "<<it4->second.second<<" (switch) \n  "<<it4->first <<" = "<<(it4->second.first ? "ON" : "OFF")<<"\n";
}

void gsOptionList::printInfo(const std::string& label)
{
    std::ostream & os = gsInfo;
    StringTable::iterator it1 = m_strings.find(label);
    if ( it1 != m_strings.end() )
    {
        os<<"* "<<it1->second.second<<" (string) \n  "<<it1->first <<" = "<<it1->second.first<<"\n";
        return;
    }
    IntTable::iterator it2 = m_ints.find(label);
    if ( it2 != m_ints.end() )
    {
        os<<"* "<<it2->second.second<<" (int) \n  "<<it2->first <<" = "<<it2->second.first<<"\n";
        return;
    }
    RealTable::iterator it3 = m_reals.find(label);
    if ( it3 != m_reals.end() )
    {
        os<<"* "<<it3->second.second<<" (real) \n  "<<it3->first <<" = "<<it3->second.first<<"\n";
        return;
    }
    SwitchTable::iterator it4 = m_switches.find(label);
    if ( it4 != m_switches.end() )
    {
        os<<"* "<<it4->second.second<<" (switch) \n  "<<it4->first <<" = "<<(it4->second.first ? "ON" : "OFF")<<"\n";
        return;
    }
    gsInfo <<"Problem: "<< label <<" does not exist.\n";
}


namespace internal
{

/** \brief Read OptionList from XML data
    \ingroup Nurbs
*/
template<>
class gsXml< gsOptionList >
{
private:
    gsXml() { }

public:
    GSXML_COMMON_FUNCTIONS(gsOptionList)
    GSXML_GET_POINTER(gsOptionList)
    static std::string tag () { return "OptionList"; }
    static std::string type() { return ""; }

    static void get_into(gsXmlNode * node, gsOptionList & result)
    {
        // read in data
        // result.add(..)
    }

    static gsXmlNode * put (const gsOptionList & obj, gsXmlTree & data)
    {
        gsXmlNode * tmp = internal::makeNode("OptionList", data);

        // Append data

        return tmp;
    }
};

}// namespace internal


};//namespace gismo

