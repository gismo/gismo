/** @file gsModuleInfo.h

    @brief POD struct exchanged during the module ABI handshake.

    Both the host (gsModuleManager) and every module library include this
    header, so they share an identical struct layout. Fields may only be
    added at the end; the reserved[] padding keeps the total size fixed.

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

/// @brief ABI handshake data returned by a module's C entry point
/// \c gismo_module_info(). \ingroup Modules
struct gismo_module_info_t
{
    uint32_t    abi_major;      ///< GISMO_MAJOR at module compile time
    uint32_t    abi_minor;      ///< GISMO_MINOR at module compile time
    const char* module_name;    ///< e.g. "gsOptim"
    const char* module_version; ///< semver, e.g. "1.0.0"
    const char* coeff_type;     ///< must match the host's real_t
    const char* index_type;     ///< must match the host's index_t
    uint8_t     reserved[96];   ///< future expansion - zero-filled
};

} // namespace gismo
