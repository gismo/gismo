/** @file gsRegistry.hpp

    @brief Implementation of gsRegistry.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsSystem/gsRegistry.h>
#include <gsSystem/gsModuleManager.h>

namespace gismo
{

template <class Base>
gsRegistry<Base>::gsRegistry()
{ }

template <class Base>
gsRegistry<Base> & gsRegistry<Base>::get()
{
    // Heap-allocated and never destroyed: registrations may fire from
    // static initializers of other libraries, so this object must
    // outlive every static-destruction order.
    // The manager is created first: its built-in modules may call get()
    // for this very registry, which must not happen inside the
    // initialiser of r.
    gsModuleManager::get();
    static gsRegistry * r = new gsRegistry;
    return *r;
}

template <class Base>
bool gsRegistry<Base>::add(const std::string & key, tFactory f)
{
    if ( NULL == f || m_entries.find(key) != m_entries.end() )
        return false;
    Entry e;
    e.factory = f;
    e.owner   = gsModuleManager::get().registering();
    m_entries.insert( std::make_pair(key, e) );
    return true;
}

template <class Base>
bool gsRegistry<Base>::has(const std::string & key) const
{ return m_entries.find(key) != m_entries.end(); }

template <class Base>
typename gsRegistry<Base>::BasePtr
gsRegistry<Base>::make(const std::string & key) const
{
    typename std::map<std::string, Entry>::const_iterator it = m_entries.find(key);
    if ( it == m_entries.end() )
        return BasePtr();
    return (it->second.factory)();
}

template <class Base>
std::vector<std::string> gsRegistry<Base>::keys() const
{
    std::vector<std::string> res;
    res.reserve(m_entries.size());
    for (typename std::map<std::string, Entry>::const_iterator it = m_entries.begin();
         it != m_entries.end(); ++it)
        res.push_back(it->first);
    return res;
}

template <class Base>
size_t gsRegistry<Base>::removeModule(const std::string & module)
{
    size_t count = 0;
    for (typename std::map<std::string, Entry>::iterator it = m_entries.begin();
         it != m_entries.end(); )
    {
        if ( it->second.owner == module )
        {
            m_entries.erase(it++);
            ++count;
        }
        else
            ++it;
    }
    return count;
}

template <class Base>
std::string gsRegistry<Base>::names() const
{
    std::string res;
    for (typename std::map<std::string, Entry>::const_iterator it = m_entries.begin();
         it != m_entries.end(); ++it)
        res += (res.empty() ? "" : ", ") + it->first;
    return res.empty() ? "(none)" : res;
}

} // namespace gismo
