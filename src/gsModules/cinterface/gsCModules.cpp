/** @file gsCModules.cpp

    @brief C interface: which modules are available.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gsCore/gsConfig.h>
#include <gsModules/gsModuleManager.h>
#include <gsModules/cinterface/gsCModules.h>
#include <gsCInterface/gsCError.h>

#include <cstring>
#include <string>

#ifndef GISMO_MODULES
#define GISMO_MODULES ""
#endif

using namespace gismo;

namespace {

// exact match of \a name in the semicolon-separated list
bool inList(const char * list, const char * name)
{
    const size_t n = std::strlen(name);
    for (const char * p = list; *p; )
    {
        const char * e = std::strchr(p, ';');
        const size_t len = e ? static_cast<size_t>(e - p) : std::strlen(p);
        if (len == n && 0 == std::strncmp(p, name, n))
            return true;
        if (!e) break;
        p = e + 1;
    }
    return false;
}

} // namespace

extern "C"
{

GISMO_EXPORT int gsC_hasModule(const char * name)
{
    GISMO_CAPI_BEGIN
    if (!name) return 0;
    if (inList(GISMO_MODULES, name)) return 1;
    return gsModuleManager::get().isLoaded(name) ? 1 : 0;
    GISMO_CAPI_END(0)
}

GISMO_EXPORT const char * gsC_modules(void)
{
    return GISMO_MODULES;
}

GISMO_EXPORT int gsC_loadModules(void)
{
    GISMO_CAPI_BEGIN
    return static_cast<int>(gsModuleManager::get().loadAll());
    GISMO_CAPI_END(-1)
}

} // extern "C"
