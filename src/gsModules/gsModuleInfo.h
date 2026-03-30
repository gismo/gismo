/** @file gsModuleInfo.h

    @brief Plain-old-data struct exchanged during the module version handshake

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <cstdint>

namespace gismo
{

/**
 * @brief POD struct exchanged between libgismo and a MODULE .so during loading.
 *
 * Both the host (gsModuleLoader) and every MODULE .so include this header so
 * that they share an identical struct layout.  Fields may only be added at the
 * end; the reserved[] padding keeps the total size fixed so that a newer module
 * can be loaded by an older host that ignores the extra fields.
 *
 * @ingroup gismo_modules
 */
struct gismo_module_info_t
{
    uint32_t    abi_major;        ///< GISMO_MAJOR at module compile time
    uint32_t    abi_minor;        ///< GISMO_MINOR at module compile time
    const char* module_name;      ///< null-terminated, e.g. "gsOptim"
    const char* module_version;   ///< null-terminated semver, e.g. "1.0.0"
    const char* coeff_type;       ///< e.g. "double"  (matches real_t)
    const char* index_type;       ///< e.g. "int"     (matches index_t)
    uint8_t     reserved[96];     ///< future expansion – must be zero-filled
};

} // namespace gismo
