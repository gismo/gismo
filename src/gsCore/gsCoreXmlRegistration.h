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
#include <gsCore/gsComposedBasis.h>
#include <gsCore/gsConstantFunction.h>
#include <gsCore/gsFunctionExpr.h>
#include <gsCore/gsComposedFunction.h>

#ifndef GISMO_BUILD_LIB
// header-only mode: the specializations must be visible to the registrars
#include <gsCore/gsCoreXml.hpp>
#endif

namespace gismo {
namespace internal {

template <class T>
void gsCoreXmlRegisterTypes()
{
    // last in the historical put chain (after all concrete spline types)
    gsXmlRegisterAs<gsGeometry<T>, gsComposedGeometry<T> >::enroll(900);
    gsXmlRegisterAs<gsBasis<T>, gsComposedBasis<T> >::enroll(900);

    gsXmlRegisterAs<gsFunction<T>, gsConstantFunction<T> >::enroll(100);
    gsXmlRegisterAs<gsFunction<T>, gsFunctionExpr<T> >    ::enroll(110);
    gsXmlRegisterAs<gsFunction<T>, gsComposedFunction<T> >::enroll(120);
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
