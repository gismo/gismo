/** @file gsOptionList.h

    @brief Provides a list of labeled parameters/options that can be
    set and accessed easily

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H. Weiner, S. Takacs
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

    /// \brief Reads value for option \a label from options.
    ///
    /// If \a label is not found, the function throws.
    std::string getString(const std::string & label) const;
    /// @copydoc gsOptionList::getString()
    index_t     getInt   (const std::string & label) const;
    /// @copydoc gsOptionList::getString()
    real_t      getReal  (const std::string & label) const;
    /// @copydoc gsOptionList::getString()
    bool        getSwitch(const std::string & label) const;

    /// \brief Reads values of an option-group \a gn into a std::vector.
    ///
    /// If \a gn is not found, the function throws.
    std::vector<std::string> getMultiString(const std::string & gn) const;
    /// @copydoc gsOptionList::getMultiString()
    std::vector<index_t>     getMultiInt   (const std::string & gn) const;
    /// @copydoc gsOptionList::getMultiString()
    std::vector<real_t>      getMultiReal  (const std::string & gn) const;

    /// \brief Reads value for option \a label from options.
    ///
    /// If \a label is not found, it defaults to \a value (otherwise \a value is not used).
    std::string askString(const std::string & label, const std::string & value = ""    ) const;
    /// @copydoc gsOptionList::askString()
    index_t     askInt   (const std::string & label, const index_t &     value = 0     ) const;
    /// @copydoc gsOptionList::askString()
    real_t      askReal  (const std::string & label, const real_t &      value = 0     ) const;
    /// @copydoc gsOptionList::askString()
    bool        askSwitch(const std::string & label, const bool &        value = false ) const;

    /*/// \brief Reads values of an option-group \a gn into a std::vector.
    ///
    /// If \a label is not found, it defaults to \a value (otherwise \a value is not used).
    //std::vector<std::string> askMultiString(const std::string & gn, const std::vector<std::string> & values = std::vector<std::string>()) const;
    /// @copydoc gsOptionList::askMultiString
    //std::vector<index_t>     askMultiInt   (const std::string & gn, const std::vector<index_t> &     values = std::vector<index_t>()    ) const;
    /// @copydoc gsOptionList::askMultiString
    //std::vector<real_t>      askMultiReal  (const std::string & gn, const std::vector<real_t> &      values = std::vector<real_t>()     ) const;*/

    /// \brief Sets an existing option \a label to be equal to \a value.
    ///
    /// If \a label is not found, the function throws.
    void setString(const std::string & label, const std::string & value);
    /// @copydoc gsOptionList::setString()
    void setInt   (const std::string & label, const index_t &     value);
    /// @copydoc gsOptionList::setString()
    void setReal  (const std::string & label, const real_t &      value);
    /// @copydoc gsOptionList::setString()
    void setSwitch(const std::string & label, const bool &        value);

    /*/// \brief Sets values of option-group \a gn from values of a std::vector.
    ///
    /// If \a gn is not found, the function throws.
    //void setMultiString(const std::string & gn, const std::vector<std::string> & values );
    /// @copydoc gsOptionList::setMultiString
    //void setMultiInt   (const std::string & gn, const std::vector<index_t> &     values );
    /// @copydoc gsOptionList::setMultiString
    //void setMultiReal  (const std::string & gn, const std::vector<real_t> &      values );*/

    /// \brief Adds a option named \a label, with description \a desc
    /// and value \a value.
    ///
    /// If an option with \a label already exists with the same type,
    /// the function overwrites it. If it has another type, the function
    /// throws.
    void addString(const std::string & label, const std::string & desc, const std::string & value );
    /// @copydoc gsOptionList::addString()
    void addInt   (const std::string & label, const std::string & desc, const index_t &     value );
    /// @copydoc gsOptionList::addString()
    void addReal  (const std::string & label, const std::string & desc, const real_t &      value );
    /// @copydoc gsOptionList::addString()
    void addSwitch(const std::string & label, const std::string & desc, const bool &        value );

    /*/// \brief Adds an option-group \a gn containing values of a std::vector.
    ///
    /// If \a gn is not found, the function throws.
    //void addMultiString(const std::string & label, const std::string & desc, const std::vector<std::string> & values);
    /// @copydoc gsOptionList::addMultiString()*/
    /// \brief Adds an option-group \a gn containing values of a std::vector.
    ///
    /// If \a gn is not found, the function throws.
    void addMultiInt(const std::string & label, const std::string & desc, const std::vector<index_t> &  values);
    /*/// @copydoc gsOptionList::addMultiString()
    //void addMultiReal  (const std::string & label, const std::string & desc, const std::vector<real_t> & values);*/

    /// \brief Removes the option named \a label (if it exists).
    void remove(const std::string& label);

    /// \brief Options for gsOptionList::update
    enum updateType {
        ignoreIfUnknwon = 0,
        addIfUnknown = 1
    };

    /// \brief Updates the object using the data from \a other.
    ///
    /// Options which do not exist in \a other, are kept unchanged.
    /// Options in \a other which do not exist in this, are kept unchanged if
    /// \a type is set to gsOptionList::ignoreIfUnknwon (default) or are added
    /// if \a type is set to gsOptionList::addIfUnknown.
    void update(const gsOptionList& other, updateType type = ignoreIfUnknwon);

    /// \brief Creates a new gsOptionList where all labels are wrapped into a groupname \a gn.
    ///
    /// Wrapping means that the label is prepended with the groupname and a dot. So, the label
    /// "Tolerance" wrapped into the group "IterativeSolver" is "InterativeSolver.Tolerance"
    gsOptionList wrapIntoGroup(const std::string & gn) const;

    /// \brief Creates a new gsOptionList, where only the options from the group \a gn are taken.
    /// In the result, the groupname and the corresponding dot are removed.
    ///
    /// If the groupname is "IterativeSolver", then a label "InterativeSolver.Tolerance" becomes
    /// "Tolerance" and a label "Basis.Degree" is ignored.
    gsOptionList getGroup(const std::string & gn) const;

    /// \brief Checks if there are labels that do not belong to a group.
    ///
    /// This is the case if there is a label which does not contain a dot.
    bool hasGlobals() const;

    /// \brief Checks if there are labels that belong to the group \a gn.
    ///
    /// This is the case if there is a label which starts with the groupname and a dot.
    bool hasGroup(const std::string & gn) const;

    /// \brief Prints this list of options to stream \a os
    std::ostream & print(std::ostream & os) const;

    /// \brief Returns the length of this list of options
    size_t size() const
    {return m_strings.size()+m_ints.size()+m_reals.size()+m_switches.size();}

    /// getAllEntries() returns a vector of those. Contains the name of its type
    /// (\a type), its label (\a label), its description (\a desc) and the string
    /// representation of its value (\a val).
    struct OptionListEntry {
        std::string type;                               ///< Type (as string)
        std::string label;                              ///< Label
        std::string desc;                               ///< Description
        std::string val;                                ///< Value (as string)
        std::ostream & print(std::ostream & os, index_t label_slot = 19) const;
    };

    /// Provides a list of all entries as vector of gsOptionList::OptionListEntry structs
    std::vector<OptionListEntry> getAllEntries() const;

    gsOptionList& operator=(const gsOptionList & other)
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

#if EIGEN_HAS_RVALUE_REFERENCES
    // Note: functions cannot be defaulted/deleted in VS2013 or older

    gsOptionList() {}

    gsOptionList(const gsOptionList & other) :
        m_strings(other.m_strings), 
        m_ints(other.m_ints), 
        m_reals(other.m_reals), 
        m_switches(other.m_switches) { }

    gsOptionList(gsOptionList && other) :
        m_strings(std::move(other.m_strings)),
        m_ints(std::move(other.m_ints)), 
        m_reals(std::move(other.m_reals)), 
        m_switches(std::move(other.m_switches)) { }

    gsOptionList& operator=(gsOptionList && other)
    {
        m_strings  = std::move(other.m_strings);
        m_ints     = std::move(other.m_ints);
        m_reals    = std::move(other.m_reals);
        m_switches = std::move(other.m_switches);
        return *this;
    }
#endif

    /// Swaps contents
    void swap(gsOptionList& other)
    {
        m_strings .swap(other.m_strings );
        m_ints    .swap(other.m_ints    );
        m_reals   .swap(other.m_reals   );
        m_switches.swap(other.m_switches);
    }

private:

    /// \brief Gives information regarding the option named \a label
    std::string getInfo(const std::string & label) const;

    /// \brief Returns true iff an option named \a label exists
    bool exists(const std::string & label) const;

    /// \brief Returns true iff a string named \a label exists
    bool isString(const std::string & label) const;
    /// \brief Returns true iff an int named \a label exists
    bool isInt(const std::string & label) const;
    /// \brief Returns true iff a real named \a label exists
    bool isReal(const std::string & label) const;
    /// \brief Returns true iff a switch named \a label exists
    bool isSwitch(const std::string & label) const;

private:
    friend class internal::gsXml<gsOptionList>;

    // Format: std::pair<Value,Description>
    typedef std::pair<std::string,std::string> StringOpt;
    typedef std::pair<index_t    ,std::string> IntOpt;
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

/// Objects of class gsOptionList can be printed using the standard io-streams
inline std::ostream &operator<<(std::ostream &os, const gsOptionList& b)
{ return b.print(os); }

/// Objects of class gsOptionList::OptionListEntry can be printed using the standard io-streams
inline std::ostream &operator<<(std::ostream &os, const gsOptionList::OptionListEntry& b)
{ return b.print(os); }

/// Objects of class gsOptionList::OptionListEntry can be ordered by label
inline bool operator< ( const gsOptionList::OptionListEntry& a, const gsOptionList::OptionListEntry& b )
{ return a.label < b.label; }


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
