/** @file gsRegistry.h

    @brief Generic key-to-factory registry for runtime-discoverable
    implementations (file readers, optimizers, solvers, ...).

    One registry instance exists per Base type and per process: libgismo
    provides the explicit instantiation (CLASS_TEMPLATE_INST
    gsRegistry<Base>;) in one of its translation units, so that get()
    resolves to a single object even when modules are loaded as separate
    shared libraries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>

#include <map>
#include <string>
#include <vector>

namespace gismo
{

/**
   @brief Type-independent part of every gsRegistry.

   Enlists itself in gsModuleManager on construction, so that the manager
   can ask every registry to drop the entries of a module being unloaded.

   @ingroup System
*/
class GISMO_EXPORT gsRegistryBase
{
public:
    virtual ~gsRegistryBase();

    /// Removes every entry registered by @a module; returns their number
    virtual size_t removeModule(const std::string & module) = 0;

protected:
    /// Enlists this registry in gsModuleManager
    gsRegistryBase();
};

/**
   @brief Key-to-factory registry for implementations of @a Base.

   Built-in and dynamic modules register a factory from their register
   function; users query by key. Each entry records the module that
   registered it (gsModuleManager::registering()), so that the entries of a
   module can be dropped when it is unloaded.

   @code
   if ( gsRegistry<gsFileIo>::get().has("ext:xml") )
       auto io = gsRegistry<gsFileIo>::get().make("ext:xml");
   @endcode

   @tparam Base type of the objects made by the factories

   @ingroup System
*/
template <class Base>
class gsRegistry : public gsRegistryBase
{
public:
    typedef memory::unique_ptr<Base> BasePtr;
    typedef BasePtr (*tFactory)();

    /// The unique registry instance for @a Base (per process)
    static gsRegistry & get();

    /// Registers @a f under @a key, owned by the module whose register
    /// function is running (empty outside one). The first registration of
    /// a key wins. Returns false, and changes nothing, if @a f is null or
    /// @a key is already registered.
    bool add(const std::string & key, tFactory f);

    /// Returns true if @a key is registered
    bool has(const std::string & key) const;

    /// Creates the object registered under @a key; a null pointer if
    /// there is none
    BasePtr make(const std::string & key) const;

    /// Returns the registered keys in lexicographic order
    std::vector<std::string> keys() const;

    /// Removes every entry whose owner is @a module; returns their number
    size_t removeModule(const std::string & module);

private:
    gsRegistry();
    gsRegistry(const gsRegistry &);
    gsRegistry & operator=(const gsRegistry &);

    /// Returns the registered keys as a single string (diagnostics)
    std::string names() const;

    struct Entry
    {
        tFactory    factory;
        std::string owner;
    };

    std::map<std::string, Entry> m_entries;
};

/// Factory that makes a @a Derived as a @a Base, for use with
/// gsRegistry<Base>::add
template <class Base, class Derived>
typename gsRegistry<Base>::BasePtr gsRegistryFactory()
{ return typename gsRegistry<Base>::BasePtr( new Derived() ); }

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRegistry.hpp)
#endif
