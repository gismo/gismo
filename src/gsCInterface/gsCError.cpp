/** @file gsCError.cpp

    @brief Per-thread error storage for the G+Smo C interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gsCore/gsConfig.h>

#ifdef GISMO_BUILD_CINTERFACE

#include <gsCInterface/gsCError.h>

#include <string>

namespace {
thread_local std::string s_lastError;
}

namespace gismo { namespace cinterface {

GISMO_EXPORT void setLastError(const char * msg)
{ s_lastError = msg ? msg : ""; }

}}

extern "C"
{

GISMO_EXPORT const char * gsCLastError(void)
{ return s_lastError.c_str(); }

GISMO_EXPORT void gsCClearError(void)
{ s_lastError.clear(); }

}

#endif // GISMO_BUILD_CINTERFACE
