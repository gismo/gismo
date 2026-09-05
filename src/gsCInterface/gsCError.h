/** @file gsCError.h

    @brief Error reporting for the G+Smo C interface.

    C++ exceptions must never cross the extern "C" boundary. Every entry
    point of the C interface is wrapped: on an exception the function
    returns a sentinel value (NULL for pointers, -1 for int, NaN for
    double, nothing for void) and the message is stored per-thread;
    callers retrieve it with gsCLastError().

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsExport.h>

#ifdef __cplusplus
extern "C"
{
#endif

/** Returns the error message of the last failed C-interface call on
    this thread, or an empty string if the last call succeeded. The
    returned pointer stays valid until the next C-interface call on the
    same thread. */
GISMO_EXPORT const char * gsCLastError(void);

/** Clears the per-thread error state. */
GISMO_EXPORT void gsCClearError(void);

#ifdef __cplusplus
} // extern "C"

// ---- implementation side (C++ only) -----------------------------------

#include <exception>

namespace gismo { namespace cinterface {
GISMO_EXPORT void setLastError(const char * msg);
}}

/// Opens the exception firewall of a C entry point
#define GISMO_CAPI_BEGIN gsCClearError(); try {

/// Closes the firewall of a value-returning entry point;
/// \a retval is returned on error (NULL / -1 / NaN, see file docs)
#define GISMO_CAPI_END(retval)                                          \
    } catch (const std::exception & e)                                  \
    { gismo::cinterface::setLastError(e.what()); return retval; }       \
    catch (...)                                                         \
    { gismo::cinterface::setLastError("gsCInterface: unknown error");   \
      return retval; }

/// Closes the firewall of a void entry point
#define GISMO_CAPI_END_VOID                                             \
    } catch (const std::exception & e)                                  \
    { gismo::cinterface::setLastError(e.what()); }                      \
    catch (...)                                                         \
    { gismo::cinterface::setLastError("gsCInterface: unknown error"); }

#endif // __cplusplus
