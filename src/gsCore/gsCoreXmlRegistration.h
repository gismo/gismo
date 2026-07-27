/** @file gsCoreXmlRegistration.h

    @brief Registers gsCore's XML readers/writers in gsXmlRegistry.

    The only file to touch when gsCore gains a serializable type.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsIO/gsXmlRegistry.h>
#include <gsCore/gsComposedGeometry.h>

namespace gismo {
namespace internal {

template <class T>
void gsCoreXmlRegisterTypes()
{
    // last in the historical put chain (after all concrete spline types)
    gsXmlRegisterAs<gsGeometry<T>, gsComposedGeometry<T> >::enroll(900);
}

#ifndef GISMO_BUILD_LIB
namespace {
struct gsCoreXmlReg
{ gsCoreXmlReg() { gsCoreXmlRegisterTypes<real_t>(); } };
static const gsCoreXmlReg s_gsCoreXmlReg;
}
#endif

} // namespace internal
} // namespace gismo
