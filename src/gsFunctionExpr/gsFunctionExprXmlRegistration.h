/** @file gsFunctionExprXmlRegistration.h

    @brief Registers gsFunctionExpr's XML reader/writer in gsXmlRegistry
    (header-only-mode registrar; in library builds gsFunctionExpr_.cpp
    registers). Priority 110 = its historical position in the gsFunction
    dispatch chain (after gsConstantFunction, before gsComposedFunction).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsIO/gsXmlRegistry.h>
#include <gsFunctionExpr/gsFunctionExpr.h>

namespace gismo {
namespace internal {

template <class T>
void gsFunctionExprXmlRegisterTypes()
{
    gsXmlRegisterAs<gsFunction<T>, gsFunctionExpr<T> >::enroll(110);
}

#ifndef GISMO_BUILD_LIB
namespace {
struct gsFunctionExprXmlReg
{ gsFunctionExprXmlReg() { gsFunctionExprXmlRegisterTypes<real_t>(); } };
static const gsFunctionExprXmlReg s_gsFunctionExprXmlReg;
}
#endif

} // namespace internal
} // namespace gismo
