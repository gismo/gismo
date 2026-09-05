/** @file gsCModules.h

    @brief C interface: which modules are available.

    Feature detection for the language bindings: a module is available
    when it was compiled into the library (in-tree modules and the
    optional modules of GISMO_OPTIONAL) or loaded at runtime through
    gsModuleManager.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsCInterface/gsCTypes.h>

#ifdef __cplusplus
extern "C"
{
#endif

/// 1 if \a name is a module compiled into the library or loaded at
/// runtime, 0 otherwise (also on error, see gsCLastError())
GISMO_EXPORT int gsC_hasModule(const char * name);

/// Semicolon-separated names of the modules compiled into the library
GISMO_EXPORT const char * gsC_modules(void);

/// Loads the runtime modules of the module directory (GISMO_MODULE_DIR
/// or the install location); returns the number loaded, -1 on error
GISMO_EXPORT int gsC_loadModules(void);

#ifdef __cplusplus
}
#endif
