/** @file gsModule.h

    @brief Entry point macros of optional modules and host-side
    registration of the built-in modules.

    A module is compiled in one of two modes:

    - dynamic: built by gismo_add_module_library(), which defines
      GISMO_MODULE_DYNAMIC. The module is a shared library exporting the C
      symbols gismo_module_info and gismo_module_register.
    - built-in: the module is part of libgismo and exports the C symbol
      gismo_module_register_<name>, called by gsRegisterBuiltinModule().

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsConfig.h>
#include <gsSystem/gsModuleInfo.h>
#include <gsSystem/gsModuleConfig.h>

/// Exports a symbol from a dynamic module, independently of the
/// export settings of libgismo
#if defined(_WIN32) || defined(__CYGWIN__)
#  define GISMO_MODULE_EXPORT __declspec(dllexport)
#else
#  define GISMO_MODULE_EXPORT __attribute__((visibility("default")))
#endif

/**
   @def GISMO_MODULE_ENTRY(name, version)
   @brief Defines the entry points of an optional module; to be followed
   by the body of the register function.

   @a name is a bare identifier equal to the CMake name of the module,
   @a version is a string literal. The macro is used once per module, at
   namespace scope (global or inside namespace gismo). All entry symbols
   have C linkage, so their names do not depend on the enclosing namespace.

   - With GISMO_MODULE_DYNAMIC (defined by gismo_add_module_library()) it
     defines the exported functions
     <tt>const gismo_module_info_t * gismo_module_info()</tt>, returning
     the ABI information of the module, and
     <tt>void gismo_module_register()</tt>, whose body follows the macro.
   - Otherwise it defines <tt>void gismo_module_register_<name>()</tt>,
     without export attribute since it stays inside libgismo; @a version
     is unused.

   @code
   GISMO_MODULE_ENTRY(gsMyModule, "1.0.0")
   {
       // register the factories of the module
   }
   @endcode
*/
#ifdef GISMO_MODULE_DYNAMIC
#define GISMO_MODULE_ENTRY(name, version)                                        \
    extern "C" GISMO_MODULE_EXPORT const ::gismo::gismo_module_info_t * gismo_module_info(); \
    extern "C" GISMO_MODULE_EXPORT void gismo_module_register();                \
    extern "C" GISMO_MODULE_EXPORT const ::gismo::gismo_module_info_t * gismo_module_info() \
    {                                                                            \
        static const ::gismo::gismo_module_info_t info =                         \
            { GISMO_MAJOR, GISMO_MINOR, #name, version,                          \
              GISMO_COEFF_TYPE_STR, GISMO_INDEX_TYPE_STR, {0} };                 \
        return &info;                                                            \
    }                                                                            \
    extern "C" GISMO_MODULE_EXPORT void gismo_module_register()
#else
#define GISMO_MODULE_ENTRY(name, version)                                        \
    extern "C" void gismo_module_register_##name();                              \
    extern "C" void gismo_module_register_##name()
#endif

namespace gismo
{

/// Names of the modules compiled into the library, terminated by NULL
const char * const * gsBuiltinModuleNames();

/// Calls the register function of the built-in module @a name; false if there is none
bool gsRegisterBuiltinModule(const char * name);

/// Calls gismo_module_register_<name>() for each built-in module, in list order
void gsRegisterBuiltinModules();

} // namespace gismo
