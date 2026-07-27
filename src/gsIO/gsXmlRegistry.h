/** @file gsXmlRegistry.h

    @brief Runtime registry for XML readers/writers of abstract base
    types (gsGeometry, gsBasis, gsFunction, ...).

    Inverts the previous design where gsIO/gsXmlUtils.hpp enumerated
    every concrete type: each module registers the (base, type-string)
    pairs it can read and a priority-ordered chain of "try-put" thunks
    per base. gsIO stays type-blind; adding a serializable type in any
    module requires no gsIO edits.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsDebug.h>
#include <gsIO/gsXml.h>

#include <algorithm>
#include <map>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

namespace gismo
{
namespace internal
{

/**
   @brief Type-erased storage of the XML reader/writer registrations.

   Keys are strings only, so registrations from runtime-loaded module
   libraries are safe. One instance per process: get() is defined in
   gsXmlRegistry.cpp (compiled into libgismo) and GISMO_EXPORT'ed;
   header-only builds fall back to an inline definition.

   \ingroup IO
*/
class gsXmlRegistry
{
public:
    /// Type-erased function pointer (cast back by gsXmlDispatch)
    typedef void (*AnyFn)();

#ifdef GISMO_BUILD_LIB
    // One instance per process, defined in gsXmlRegistry.cpp inside
    // libgismo; module libraries resolve the exported symbol
    static GISMO_EXPORT gsXmlRegistry & get();
#else
    // Header-only mode: single binary, inline local static suffices
    static gsXmlRegistry & get()
    { static gsXmlRegistry * r = new gsXmlRegistry; return *r; }
#endif

    /// Registers a reader for objects of \a type under \a baseId
    /// (idempotent: the first registration wins)
    void addGet(const std::string & baseId, const std::string & type,
                AnyFn fn)
    { m_getters[baseId].insert( std::make_pair(type, fn) ); }

    /// Registers a try-put thunk for \a baseId at \a priority (lower
    /// runs first; equal priorities run in registration order)
    void addPut(const std::string & baseId, int priority, AnyFn fn)
    {
        std::vector<PutEntry> & chain = m_putters[baseId];
        for (size_t i = 0; i != chain.size(); ++i)
            if (chain[i].second == fn) return; // idempotent
        chain.push_back( std::make_pair(priority, fn) );
    }

    /// Returns the reader for (\a baseId, \a type), or NULL
    AnyFn findGet(const std::string & baseId,
                  const std::string & type) const
    {
        GetterMaps::const_iterator b = m_getters.find(baseId);
        if (b == m_getters.end()) return NULL;
        std::map<std::string,AnyFn>::const_iterator t =
            b->second.find(type);
        return t == b->second.end() ? NULL : t->second;
    }

    /// Fills \a out with the put chain of \a baseId in priority order
    void putChain(const std::string & baseId,
                  std::vector<AnyFn> & out) const
    {
        out.clear();
        PutterMaps::const_iterator b = m_putters.find(baseId);
        if (b == m_putters.end()) return;
        std::vector<PutEntry> chain = b->second; // stable sort a copy
        std::stable_sort(chain.begin(), chain.end(), lessPriority);
        for (size_t i = 0; i != chain.size(); ++i)
            out.push_back(chain[i].second);
    }

    /// Known type strings for \a baseId (diagnostics)
    std::string knownTypes(const std::string & baseId) const
    {
        std::string res;
        GetterMaps::const_iterator b = m_getters.find(baseId);
        if (b != m_getters.end())
            for (std::map<std::string,AnyFn>::const_iterator t =
                     b->second.begin(); t != b->second.end(); ++t)
                res += (res.empty() ? "" : ", ") + t->first;
        return res.empty() ? "(none)" : res;
    }

private:
    gsXmlRegistry() { }
    gsXmlRegistry(const gsXmlRegistry&);

    typedef std::pair<int,AnyFn> PutEntry;
    static bool lessPriority(const PutEntry & a, const PutEntry & b)
    { return a.first < b.first; }

    typedef std::map<std::string, std::map<std::string,AnyFn> > GetterMaps;
    typedef std::map<std::string, std::vector<PutEntry> >       PutterMaps;
    GetterMaps m_getters;
    PutterMaps m_putters;
};

/**
   @brief Typed dispatch used by the abstract-base gsXml
   specializations: reads resolve the node's "type" attribute against
   the registry; writes walk the registered dynamic_cast chain.
*/
template <class Base>
struct gsXmlDispatch
{
    typedef Base*      (*GetFn)(gsXmlNode*);
    typedef gsXmlNode* (*PutFn)(const Base&, gsXmlTree&);

    static std::string baseId() { return typeid(Base*).name(); }

    static Base * get(gsXmlNode * node)
    {
        gsXmlAttribute * at = node->first_attribute("type");
        const std::string type = at ? at->value() : "";
        gsXmlRegistry::AnyFn fn =
            gsXmlRegistry::get().findGet(baseId(), type);
        if (!fn)
        {
            gsWarn << "gsXml: no registered reader for type \"" << type
                   << "\". Known: "
                   << gsXmlRegistry::get().knownTypes(baseId()) << "\n";
            return NULL;
        }
        return reinterpret_cast<GetFn>(fn)(node);
    }

    static gsXmlNode * put(const Base & obj, gsXmlTree & data)
    {
        std::vector<gsXmlRegistry::AnyFn> chain;
        gsXmlRegistry::get().putChain(baseId(), chain);
        for (size_t i = 0; i != chain.size(); ++i)
            if (gsXmlNode * n = reinterpret_cast<PutFn>(chain[i])(obj, data))
                return n;
        gsWarn << "gsXml: no registered writer accepted an object of base "
               << typeid(Base).name() << "\n";
        return NULL;
    }
};

/**
   @brief Registration helper: enrolls gsXml<Derived> as reader/writer
   for base type \a Base.

   The get thunk performs the Derived* -> Base* upcast; the put thunk
   dynamic_casts and returns NULL on mismatch (chain semantics identical
   to the previous GSXML_PUT_DYNAMIC_CAST macros).
*/
template <class Base, class Derived>
struct gsXmlRegisterAs
{
    static Base * getThunk(gsXmlNode * node)
    { return gsXml<Derived>::get(node); }

    static gsXmlNode * tryPutThunk(const Base & b, gsXmlTree & data)
    {
        const Derived * p = dynamic_cast<const Derived*>(&b);
        return p ? gsXml<Derived>::put(*p, data) : NULL;
    }

    /// Registers the reader for \a Base under gsXml<Derived>::type()
    static void enrollGet()
    {
        gsXmlRegistry::get().addGet(
            gsXmlDispatch<Base>::baseId(), gsXml<Derived>::type(),
            reinterpret_cast<gsXmlRegistry::AnyFn>(
                (typename gsXmlDispatch<Base>::GetFn) &getThunk));
    }

    /// Registers the writer for \a Base at \a putPriority
    static void enrollPut(int putPriority)
    {
        gsXmlRegistry::get().addPut(
            gsXmlDispatch<Base>::baseId(), putPriority,
            reinterpret_cast<gsXmlRegistry::AnyFn>(
                (typename gsXmlDispatch<Base>::PutFn) &tryPutThunk));
    }

    /// Registers both reader and writer (the common case)
    static void enroll(int putPriority)
    { enrollGet(); enrollPut(putPriority); }
};

/*
  Registration macros for use inside a class's explicit-instantiation
  translation unit (the *_.cpp that also does CLASS_TEMPLATE_INST
  internal::gsXml<Derived>), inside namespace gismo. Template arguments
  containing commas must be wrapped with TMPLA2/TMPLA3 (see gsXml.h).
  This placement is deliberate: only the instantiation unit may include
  the class .hpp, so only there are the gsXml specializations visible.
*/
#ifndef GISMO_COMMA
#define GISMO_COMMA ,
#endif
#define GISMO_XML_REG_IMPL2(Regcall, id)                                     namespace { struct gsXmlReg_##id                                         { gsXmlReg_##id() { Regcall; } } const gsXmlRegObj_##id; }
#define GISMO_XML_REG_IMPL(Regcall, id) GISMO_XML_REG_IMPL2(Regcall, id)

/// Registers Derived under Base for reading and writing (put priority p)
#define GISMO_XML_REGISTER(Base, Derived, p) GISMO_XML_REG_IMPL(              (internal::gsXmlRegisterAs<Base GISMO_COMMA Derived>::enroll(p)), __COUNTER__)
/// Read-only registration
#define GISMO_XML_REGISTER_GET(Base, Derived) GISMO_XML_REG_IMPL(             (internal::gsXmlRegisterAs<Base GISMO_COMMA Derived>::enrollGet()), __COUNTER__)
/// Write-only registration (put priority p)
#define GISMO_XML_REGISTER_PUT(Base, Derived, p) GISMO_XML_REG_IMPL(          (internal::gsXmlRegisterAs<Base GISMO_COMMA Derived>::enrollPut(p)), __COUNTER__)

} // namespace internal
} // namespace gismo
