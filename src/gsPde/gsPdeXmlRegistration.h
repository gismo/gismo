/** @file gsPdeXmlRegistration.h

    @brief Registers gsPde's XML readers in gsXmlRegistry.

    The only file to touch when gsPde gains a serializable type.
    Reading only: writing PDEs was never supported.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsIO/gsXmlRegistry.h>
#include <gsPde/gsPoissonPde.h>
#include <gsPde/gsSurfacePoissonPde.h>

namespace gismo {
namespace internal {

template <class T>
void gsPdeXmlRegisterTypes()
{
    gsXmlRegisterAs<gsPde<T>, gsPoissonPde<T> >       ::enrollGet();
    gsXmlRegisterAs<gsPde<T>, gsSurfacePoissonPde<T> >::enrollGet();
}

#ifndef GISMO_BUILD_LIB
namespace {
struct gsPdeXmlReg
{ gsPdeXmlReg() { gsPdeXmlRegisterTypes<real_t>(); } };
static const gsPdeXmlReg s_gsPdeXmlReg;
}
#endif

} // namespace internal
} // namespace gismo
