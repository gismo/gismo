/** @file gsModelingXmlRegistration.h

    @brief Registers gsModeling's XML readers/writers in gsXmlRegistry.

    The only file to touch when gsModeling gains a serializable type.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsIO/gsXmlRegistry.h>
#include <gsModeling/gsTrimSurface.h>

namespace gismo {
namespace internal {

template <class T>
void gsModelingXmlRegisterTypes()
{
    // historically writable through gsGeometry but not readable
    gsXmlRegisterAs<gsGeometry<T>, gsTrimSurface<T> >::enrollPut(210);
}

#ifndef GISMO_BUILD_LIB
namespace {
struct gsModelingXmlReg
{ gsModelingXmlReg() { gsModelingXmlRegisterTypes<real_t>(); } };
static const gsModelingXmlReg s_gsModelingXmlReg;
}
#endif

} // namespace internal
} // namespace gismo
