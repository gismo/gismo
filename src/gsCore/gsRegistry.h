/** @file gsRegistry.h

    @brief Generic name-to-factory registry for runtime-discoverable
    implementations (optimizers, readers, solvers, ...).

    One registry instance exists per Base type and per process: the
    owning module of \a Base provides the explicit instantiation
    (CLASS_TEMPLATE_INST gsRegistry<Base>;) in one of its compiled
    translation units, so that the GISMO_EXPORT'ed get() resolves to a
    single authoritative object even when modules are loaded as
    separate shared libraries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsDebug.h>

#include <map>
#include <string>
#include <vector>

namespace gismo
{

/**
   @brief Name-to-factory registry for implementations of \a Base.

   Optional modules register a factory on load (statically or from
   their gismo_module_register() entry point); users query by name:

   \code
   if ( gsRegistry<gsOptimizer<real_t> >::get().has("gsOptim-BFGS") )
       auto opt = gsRegistry<gsOptimizer<real_t> >::get().make("gsOptim-BFGS");
   \endcode

   \ingroup Core
*/
template <class Base>
class gsRegistry
{
public:
    typedef memory::unique_ptr<Base> BasePtr;
    typedef BasePtr (*tFactory)();

    /// The unique registry instance for \a Base (per process)
    static gsRegistry & get();

    /// Registers \a factory under \a name (idempotent: the first
    /// registration of a name wins, repeats are ignored)
    void add(const std::string & name, tFactory factory)
    { m_factories.insert( std::make_pair(name, factory) ); }

    /// Returns true if \a name is registered
    bool has(const std::string & name) const
    { return m_factories.find(name) != m_factories.end(); }

    /// Creates the implementation registered under \a name
    BasePtr make(const std::string & name) const
    {
        typename FactoryMap::const_iterator it = m_factories.find(name);
        GISMO_ENSURE( it != m_factories.end(), "gsRegistry: \"" << name
                      << "\" is not registered. Available: " << names() );
        return (it->second)();
    }

    /// Returns the registered names
    std::vector<std::string> list() const
    {
        std::vector<std::string> res;
        res.reserve(m_factories.size());
        for (typename FactoryMap::const_iterator it = m_factories.begin();
             it != m_factories.end(); ++it)
            res.push_back(it->first);
        return res;
    }

    /// Returns the registered names as a single string (diagnostics)
    std::string names() const
    {
        std::string res;
        for (typename FactoryMap::const_iterator it = m_factories.begin();
             it != m_factories.end(); ++it)
            res += (res.empty() ? "" : ", ") + it->first;
        return res.empty() ? "(none)" : res;
    }

private:
    gsRegistry() { }
    gsRegistry(const gsRegistry&);
    gsRegistry& operator=(const gsRegistry&);

    typedef std::map<std::string, tFactory> FactoryMap;
    FactoryMap m_factories;
};

/// Convenience factory for registering \a Derived under a base type:
/// registry.add("name", gsRegistryFactory<Base,Derived>);
template <class Base, class Derived>
typename gsRegistry<Base>::BasePtr gsRegistryFactory()
{ return typename gsRegistry<Base>::BasePtr( new Derived() ); }

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRegistry.hpp)
#endif
