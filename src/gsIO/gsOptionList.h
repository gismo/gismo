/** @file gsOptionList.h

    @brief Provides a list of labeled parameters/options that can be
    set and accessed easily

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H. Weiner
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsIO/gsXml.h>

namespace gismo
{

/** 
    @brief Class which holds a list of parameters/options, and
    provides easy access to them.
    
    Every parameter has a unique label, and its value can be retrieved
    using that.
    
    \ingroup IO
*/
class GISMO_EXPORT gsOptionList
{
public:
    
    /// \brief Reads value for option \a label from options. If \a
    /// label is not found, the function throws
    std::string getString(const std::string & label) const;
    /// @copydoc gsOptionList::getString
    int         getInt   (const std::string & label) const;
    /// @copydoc gsOptionList::getString
    real_t      getReal  (const std::string & label) const;
    /// @copydoc gsOptionList::getString
    bool        getSwitch(const std::string & label) const;

    /// \brief Reads value for option \a label from options. If \a
    /// label is not found, it defaults to \a val (otherwise \a val is
    /// not used)
    std::string askString(const std::string & label, const std::string& val = "") const;
    /// @copydoc gsOptionList::getString
    int         askInt   (const std::string & label, const int &        val = 0 ) const;
    /// @copydoc gsOptionList::getString
    real_t      askReal  (const std::string & label, const real_t &     val = 0 ) const;
    /// @copydoc gsOptionList::getString
    bool        askSwitch(const std::string & label, const bool &   val = false ) const;

    /// \brief Sets an existing option \a label to be equal to \a
    /// value
    void setString(const std::string & label, const std::string & value);
    /// @copydoc gsOptionList::setString
    void setInt   (const std::string & label, const int & res          );
    /// @copydoc gsOptionList::setString
    void setReal  (const std::string & label, const real_t & res       );
    /// @copydoc gsOptionList::setString
    void setSwitch(const std::string & label, const bool & res         );

    /// \brief Adds a new option named \a label, with description \a
    /// desc and current value \a value
    void addString(const std::string& label, const std::string& desc, const std::string& value);
    /// @copydoc gsOptionList::addString
    void addInt   (const std::string& label, const std::string& desc, const int& value);
    /// @copydoc gsOptionList::addString
    void addReal  (const std::string& label, const std::string& desc, const real_t& value);
    /// @copydoc gsOptionList::addString
    void addSwitch(const std::string& label, const std::string& desc, const bool& value);

    /// \brief Prints this list of options to stream \a os
    std::ostream & print(std::ostream & os) const;

    /// \brief Returns the length of this list of options
    int size() const
    {return m_strings.size()+m_ints.size()+m_reals.size()+m_switches.size();}

    // /*
    typedef struct {
        std::string type;
        std::string label;
        std::string desc;
        std::string val;
    } OptionListEntry;
    
    std::vector<OptionListEntry> getAllEntries() const;
    //*/

    gsOptionList& operator=(const gsOptionList& other)
    { // Note: implcitly degerated operator was buggy on some platforms
        if (this != &other)
        {
            m_strings  = other.m_strings;
            m_ints     = other.m_ints;
            m_reals    = other.m_reals;
            m_switches = other.m_switches;            
        }
        return *this;
    }
    
    /// \brief Updates the object using the data from \a other. Fieleds unkown to \a other,
    /// are kept unchanged. Fields in \a other unknown to this are ignored.
    void update(const gsOptionList& other);

private:

    /// \brief Prints information regarding the option nnamed \a label
    void printInfo(const std::string & label) const;

    /// \brief Returns true iff an option named \a label exists
    bool exists(const std::string & label) const;
    
private:
    friend class internal::gsXml<gsOptionList>;
    
    // Format: std::pair<Value,Description>
    typedef std::pair<std::string,std::string> StringOpt;
    typedef std::pair<int        ,std::string> IntOpt;
    typedef std::pair<real_t     ,std::string> RealOpt;
    typedef std::pair<bool       ,std::string> SwitchOpt;

    // Format: std::map<Label, std::pair<Value,Description> >
    typedef std::map<std::string,StringOpt> StringTable;
    typedef std::map<std::string,IntOpt>    IntTable;
    typedef std::map<std::string,RealOpt>   RealTable;
    typedef std::map<std::string,SwitchOpt> SwitchTable;
    
    StringTable m_strings;  ///< String-valued options/parameters
    IntTable    m_ints;     ///< Integer-valued options/parameters
    RealTable   m_reals;    ///< Real-valued options/parameters
    SwitchTable m_switches; ///< Switches (ON/OFF) options/parameters

}; // class gsOptionList

/// Print (as string) operator to be used by all derived classes
inline std::ostream &operator<<(std::ostream &os, const gsOptionList& b)
{return b.print(os); }

namespace internal
{

/** \brief Read OptionList from XML data
    \ingroup IO
*/
template<>
class GISMO_EXPORT gsXml<gsOptionList>
{
private:
    gsXml();
public:
    GSXML_COMMON_FUNCTIONS(gsOptionList)
    GSXML_GET_POINTER(gsOptionList)
    static std::string tag () { return "OptionList"; }
    static std::string type() { return ""; }

    static void get_into(gsXmlNode * node, gsOptionList & result);
    static gsXmlNode * put (const gsOptionList & obj, gsXmlTree & data);
};

}

} // namespace gismo
