/** @file gsOptionList.hpp

    @brief Provides a list of labeled parameters/options that can be
    set and accessed easily

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H. Weiner, S. Takacs
*/

namespace gismo
{

template<typename T>
T gsOptionList::getReal(const std::string & label) const
{
    RealTable::const_iterator it = m_reals.find(label);
    GISMO_ENSURE(it!=m_reals.end(), "Invalid request (getReal): "<<label<<" is not a real; it is "<<getInfo(label)<<".");
    return (T)(it->second.first);
}

template<typename T>
std::vector<T> gsOptionList::getMultiReal(const std::string & gn) const
{
    GISMO_ASSERT(hasGroup(gn), "Invalid request (getMultiReal): The group " + gn + " does not exist.");

    std::vector<T> result;

    const std::string search = gn + ".";

    // add reals to vector
    for (RealTable::const_iterator it = m_reals.begin(); it != m_reals.end(); it++)
        if (util::starts_with(it->first, search) && !util::ends_with(it->first, "Size"))
            result.push_back((T)it->second.first);

    return result;
}

template<typename T>
T gsOptionList::askReal(const std::string & label,
                             const T & value) const
{
    RealTable::const_iterator it = m_reals.find(label);
#if defined(GISMO_EXTRA_DEBUG)
    if ( it == m_reals.end() && exists(label) )
        gsWarn << "Invalid request (askReal): "<<label<<" is given, but not a real; it is "<<getInfo(label)<<".\n";
#endif
    return (T)( it == m_reals.end() ? value : it->second.first);
}

template<typename T>
void gsOptionList::setReal(const std::string & label,
                           const T & value)
{
    RealTable::iterator it = m_reals.find(label);
    GISMO_ENSURE(it!=m_reals.end(), "Invalid request (setReal): "<<label<<" is not a real; it is "<<getInfo(label)<<".");
    it->second.first = (real_t)value;
}

template<typename T>
void gsOptionList::addReal(const std::string & label,
                           const std::string & desc,
                           const T & value)
{
    GISMO_ENSURE( !( isString(label) || isInt(label) || isSwitch(label) ),
         "Invalid request (addReal): Option "<<label<<" already exists, but not as a real; it is "<<getInfo(label)<<"." );
    //GISMO_ASSERT( !exists(label), "Option "<<label<<" already exists." );
    m_reals[label] = std::make_pair( (real_t)(value),desc);
}



} //namespace gismo
